#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vegetation fraction within 1 km of the site from Sentinel-2 L2A at 10 m.

For each site the year is split into six two-month bins and the clearest L2A scene of
each bin with at most max_cloud % cloud is used. Per scene, red and near infrared are
read for the pixels within 1 km of the site, the scene classification layer drops
cloud, shadow, snow and saturated pixels, and water counts as zero cover. NDVI is
scaled linearly to a cover fraction between NDVI_SOIL and NDVI_VEG.

veg is the mean over the footprint of the per-pixel peak cover over the scenes.
share_ndvi_gt_0.15 is the share of footprint pixels whose peak NDVI is above 0.15.
The per-pixel peak NDVI of each site is saved in peak_dir when it is set, so other
thresholds can be tried without reading the scenes again.

Scenes and bands come from the Planetary Computer cloud optimised GeoTIFFs, signed
with a SAS token that is renewed every 30 minutes. The cache holds one row per site,
the same format as s2_fvc_2019_1km_sites_v2.csv of add_vegfrac_s2.py.
"""
import json
import os
import sys
import time

import numpy as np

import extractors.footprint_utils as fu

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-2-l2a"
NDVI_SOIL, NDVI_VEG = 0.15, 0.90
NDVI_MIN = 0.15
# Scene classification: 4 vegetation, 5 bare, 6 water, 7 unclassified.
SCL_KEEP, SCL_WATER = (4, 5, 6, 7), 6

GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}
_token = {"value": None, "time": 0.0}


def fetch_json(url, data=None, headers=None):
    return json.loads(fu.fetch_url(url, data, headers, wait=5))


def sign(href):
    # The read token lasts about an hour, so renew it well before then.
    if time.time() - _token["time"] > 1800:
        _token["value"] = fetch_json(SAS)["token"]
        _token["time"] = time.time()
    return f"/vsicurl/{href}?{_token['value']}"


def year_of(info):
    return int(info.get("years", [2019])[0])


def scenes(lat, lon, info):
    """Clearest L2A scene of each two-month bin of the year at the site."""
    year = year_of(info)
    query = json.dumps({
        "collections": ["sentinel-2-l2a"],
        "intersects": {"type": "Point", "coordinates": [lon, lat]},
        "datetime": f"{year}-01-01/{year}-12-31",
        "query": {"eo:cloud_cover": {"lt": info.get("max_cloud", 30)}},
        "limit": 500,
    }).encode()
    items = fetch_json(STAC, query, {"Content-Type": "application/json"})["features"]
    best = {}
    for it in items:
        b = (int(it["properties"]["datetime"][5:7]) - 1) // 2
        cloud = it["properties"]["eo:cloud_cover"]
        if b not in best or cloud < best[b]["properties"]["eo:cloud_cover"]:
            best[b] = it
    return list(best.values())


def offset(item):
    """Processing baseline 04.00 and later adds 1000 to the reflectance digits."""
    return 1000 if float(item["properties"]["s2:processing_baseline"]) >= 4 else 0


def scene_cover(item, lat, lon):
    """Cover and NDVI of the pixels within the radius, on the scene 10 m grid."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import transform
    from rasterio.windows import bounds as win_bounds, from_bounds
    off = offset(item)
    a = item["assets"]
    r = fu.RADIUS_M
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(sign(a["B04"]["href"])) as src:
            crs = str(src.crs)
            x, y = transform("EPSG:4326", src.crs, [lon], [lat])
            x, y = x[0], y[0]
            win = from_bounds(x - r - 10, y - r - 10, x + r + 10, y + r + 10,
                              src.transform)
            win = win.round_offsets().round_lengths()
            red = src.read(1, window=win, boundless=True, fill_value=0)
            tr = src.window_transform(win)
            bounds = win_bounds(win, src.transform)
        with rasterio.open(sign(a["B08"]["href"])) as src:
            nir = src.read(1, window=win, boundless=True, fill_value=0)
        with rasterio.open(sign(a["SCL"]["href"])) as src:
            w20 = from_bounds(*bounds, src.transform).round_offsets().round_lengths()
            scl = src.read(1, window=w20, boundless=True, fill_value=0,
                           out_shape=red.shape, resampling=Resampling.nearest)
    red = (red.astype(float) - off) / 10000
    nir = (nir.astype(float) - off) / 10000
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = (nir - red) / (nir + red)
    fvc = np.clip((ndvi - NDVI_SOIL) / (NDVI_VEG - NDVI_SOIL), 0, 1)
    fvc[scl == SCL_WATER] = 0
    fvc[~np.isin(scl, SCL_KEEP) | (red <= 0) | (nir <= 0)] = np.nan
    px = tr.c + (np.arange(red.shape[1]) + 0.5) * tr.a
    py = tr.f + (np.arange(red.shape[0]) + 0.5) * tr.e
    xx, yy = np.meshgrid(px, py)
    inside = np.hypot(xx - x, yy - y) <= r
    return crs, inside, fvc, ndvi


def stack_scenes(site, lat, lon, info, scene_fn):
    """
    Run scene_fn on the clearest scenes, keeping those on the grid of the first one.
    Returns the cache record so far, the stacked outputs and the footprint mask.
    """
    rec = {"site": site}
    try:
        scs = scenes(lat, lon, info)
    except RuntimeError:
        return rec, [], None
    outs, used, crs0, inside0 = [], [], None, None
    for it in scs:
        try:
            crs, inside, *res = scene_fn(it, lat, lon)
        except Exception:
            continue
        # Scenes from a neighbouring UTM zone are on a different grid.
        if crs0 is None:
            crs0, inside0 = crs, inside
        if crs != crs0 or res[0].shape != inside0.shape:
            continue
        outs.append(res)
        used.append(it["properties"]["datetime"][:10].replace("-", ""))
    rec["n_scenes"] = len(used)
    rec["dates"] = " ".join(sorted(used))
    return rec, outs, inside0


def peak(arr):
    return np.nanmax(np.where(np.isfinite(arr), arr, -1), axis=0)


def fetch_site(site, lat, lon, info):
    rec, outs, inside = stack_scenes(site, lat, lon, info, scene_cover)
    if not outs:
        return [rec]
    arr = np.stack([o[0] for o in outs])
    ok = np.isfinite(arr).any(axis=0) & inside
    pcover = np.where(ok, peak(arr), np.nan)
    pndvi = np.nanmax(np.where(np.isfinite(arr), np.stack([o[1] for o in outs]), -1),
                      axis=0)
    if info.get("peak_dir"):
        os.makedirs(info["peak_dir"], exist_ok=True)
        np.save(os.path.join(info["peak_dir"], f"{site}.npy"),
                pndvi[ok].astype(np.float16))
    rec.update({"n_pixels": int(inside.sum()), "n_valid": int(ok.sum()),
                "veg": float(np.nanmean(pcover[ok])) if ok.any() else np.nan,
                f"share_ndvi_gt_{NDVI_MIN}": float((pndvi[ok] > NDVI_MIN).mean())
                if ok.any() else np.nan,
                "share_ndvi_gt_0.3": float((pndvi[ok] > 0.3).mean()) if ok.any()
                else np.nan})
    return [rec]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    row = rows.iloc[-1]
    year = year_of(info)
    return {"veg": {year: float(row.get("veg", np.nan))},
            "veg_ndvi_gt015": {year: float(row.get(f"share_ndvi_gt_{NDVI_MIN}",
                                                   np.nan))}}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
