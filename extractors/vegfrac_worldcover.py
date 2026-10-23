#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESA WorldCover 2021 v200 cover fractions within 1 km of the site.

Each fraction is the share of valid 10 m pixels of a class among the pixels whose
centres lie within 1 km of the site. Windows are cut through the Planetary Computer
data API on the 10 m grid without resampling, since the blob storage behind it is not
always reachable. The cache holds one row per site with the share of every class, the
same format as worldcover_2021_v200_1km_sites.csv of add_vegfrac_worldcover.py.

Tree is tree cover and mangroves. Grass is grassland, herbaceous wetland and moss and
lichen. Vegetation is tree, shrub, grass and cropland.
"""
import io
import json
import sys

import numpy as np

import extractors.footprint_utils as fu

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
CROP = ("https://planetarycomputer.microsoft.com/api/data/v1/item/bbox/"
        "{minx},{miny},{maxx},{maxy}/{w}x{h}.tif?collection=esa-worldcover&"
        "item={item}&assets=map&resampling=nearest")
PX = 1 / 12000
CLASSES = {10: "tree", 20: "shrub", 30: "grass", 40: "crops", 50: "built",
           60: "bare", 70: "snow", 80: "water", 90: "wetland", 95: "mangrove",
           100: "moss"}
GROUPS = {
    "tree": ["tree", "mangrove"],
    "shrub": ["shrub"],
    "grass": ["grass", "wetland", "moss"],
    "crop": ["crops"],
    "bare": ["bare"],
    "veg": ["tree", "mangrove", "shrub", "grass", "wetland", "moss", "crops"],
}


def fetch_site(site, lat, lon, info):
    import rasterio
    year = str(info.get("years", [2021])[0])
    half_lat, half_lon = fu.half_box(lat)
    # Snap to the 10 m grid so the data API returns native pixels.
    minx = np.floor((lon - half_lon) / PX) * PX
    maxx = np.ceil((lon + half_lon) / PX) * PX
    miny = np.floor((lat - half_lat) / PX) * PX
    maxy = np.ceil((lat + half_lat) / PX) * PX
    query = json.dumps({"collections": ["esa-worldcover"],
                        "bbox": [minx, miny, maxx, maxy], "limit": 20}).encode()
    items = json.loads(fu.fetch_url(STAC, query,
                                    {"Content-Type": "application/json"}))
    items = [f for f in items["features"]
             if f["properties"]["start_datetime"].startswith(year)]
    codes, plats, plons = [], [], []
    for f in items:
        b = f["bbox"]
        x0, x1 = max(minx, b[0]), min(maxx, b[2])
        y0, y1 = max(miny, b[1]), min(maxy, b[3])
        w, h = int(round((x1 - x0) / PX)), int(round((y1 - y0) / PX))
        if w <= 0 or h <= 0:
            continue
        url = CROP.format(minx=f"{x0:.8f}", miny=f"{y0:.8f}", maxx=f"{x1:.8f}",
                          maxy=f"{y1:.8f}", w=w, h=h, item=f["id"])
        with rasterio.open(io.BytesIO(fu.fetch_url(url))) as src:
            a = src.read(1)
            valid = src.read(2) > 0
            t = src.transform
        lon2d, lat2d = np.meshgrid(t.c + (np.arange(w) + 0.5) * t.a,
                                   t.f + (np.arange(h) + 0.5) * t.e)
        codes.append(a[valid])
        plats.append(lat2d[valid])
        plons.append(lon2d[valid])
    rec = {"site": site, "n_items": len(items)}
    if codes:
        codes, plats, plons = map(np.concatenate, (codes, plats, plons))
        c = codes[fu.haversine(lat, lon, plats, plons) <= fu.RADIUS_M]
        rec["n_pixels"] = int(c.size)
        for code, name in CLASSES.items():
            rec[name] = float((c == code).mean()) if c.size else np.nan
    return [rec]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    row = rows.iloc[-1]
    year = int(info.get("years", [2021])[0])
    out = {}
    for g, cols in GROUPS.items():
        vals = row.reindex(cols).astype(float)
        out[g] = {year: float(vals.sum()) if vals.notna().any() else np.nan}
    return out


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
