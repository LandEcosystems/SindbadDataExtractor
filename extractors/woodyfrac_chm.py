#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Woody vegetation fraction within 1 km of the site from the Meta and WRI global 1 m
canopy height map v1 (Tolan et al. 2024).

woody is the share of canopy height grid cells of at least 2 m among the grid cells
whose centres lie within 1 km of the site, tree5 the share of at least 5 m and
mean_height their mean height in m. Grid cells are read remotely from the public
cloud optimised GeoTIFFs, whose tile index is downloaded once to tiles_path. The
cache holds one row per site, the same format as meta_chm_1km_sites.csv of
add_woody_chm.py.
"""
import json
import os
import sys
import urllib.request

import numpy as np

import extractors.footprint_utils as fu

ROOT = "https://dataforgood-fb-data.s3.amazonaws.com/forests/v1/alsgedi_global_v6_float"
WOODY_M, TREE_M = 2, 5
R_MERC = 6378137.0

GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif", "GDAL_CACHEMAX": 256,
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}
_tiles = {}


def load_tiles(path):
    from shapely.geometry import shape
    if path not in _tiles:
        with fu.file_lock(path):
            if not os.path.exists(path):
                urllib.request.urlretrieve(f"{ROOT}/tiles.geojson", path)
        feats = json.load(open(path))["features"]
        _tiles[path] = [(f["properties"]["tile"], shape(f["geometry"]))
                        for f in feats]
    return _tiles[path]


def merc(lat, lon):
    x = np.radians(lon) * R_MERC
    y = R_MERC * np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))
    return x, y


def fetch_site(site, lat, lon, info):
    import rasterio
    from rasterio.windows import Window, from_bounds
    from shapely.geometry import box
    half_lat, half_lon = fu.half_box(lat)
    bb = box(lon - half_lon, lat - half_lat, lon + half_lon, lat + half_lat)
    hits = [t for t, geom in load_tiles(info["tiles_path"]) if geom.intersects(bb)]
    x, y = merc(lat, lon)
    # Web Mercator stretches distances by 1 / cos(lat).
    h = fu.RADIUS_M / np.cos(np.radians(lat))
    heights = []
    with rasterio.Env(**GDAL_ENV):
        for t in hits:
            with rasterio.open(f"/vsicurl/{ROOT}/chm/{t}.tif") as src:
                win = from_bounds(x - h, y - h, x + h, y + h, src.transform)
                win = win.round_offsets().round_lengths()
                win = win.intersection(Window(0, 0, src.width, src.height))
                a = src.read(1, window=win)
                tr = src.window_transform(win)
            px = tr.c + (np.arange(a.shape[1]) + 0.5) * tr.a
            py = tr.f + (np.arange(a.shape[0]) + 0.5) * tr.e
            plon = np.degrees(px / R_MERC)
            plat = np.degrees(2 * np.arctan(np.exp(py / R_MERC)) - np.pi / 2)
            lon2d, lat2d = np.meshgrid(plon, plat)
            heights.append(a[fu.haversine(lat, lon, lat2d, lon2d) <= fu.RADIUS_M])
    rec = {"site": site, "n_tiles": len(hits)}
    if heights:
        v = np.concatenate(heights).astype(float)
        rec.update({"n_grids": v.size, "woody": float((v >= WOODY_M).mean()),
                    "tree5": float((v >= TREE_M).mean()),
                    "mean_height": float(v.mean())})
    return [rec]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    row = rows.iloc[-1]
    year = int(info.get("years", [2019])[0])
    return {k: {year: float(row.get(k, np.nan))}
            for k in ["woody", "tree5", "mean_height"]}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
