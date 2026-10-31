#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NDVI of every clear Landsat 4 to 9 scene at the site, averaged over the 30 m grid
cells whose centres lie within 1 km of the site.

Scenes come from the Landsat Collection 2 Level-2 surface reflectance of the
Planetary Computer, signed with a SAS token that is renewed every 30 minutes. Only
scenes of the tiers in the json, by default T1, and with at most max_cloud % cloud
over the whole scene are read, from first to last. Per scene, red and near infrared
are read for the grid cells within the radius, and a grid cell counts as clear when
the QA_PIXEL bits of fill, dilated cloud, cirrus, cloud, cloud shadow and snow are
all zero and both reflectances lie between 0 and 1. The gaps of Landsat 7 after the
failure of its scan line corrector in 2003 are fill and so not clear. The NDVI of
Landsat 4, 5 and 7 and that of Landsat 8 and 9 come from different sensors and are
not harmonised.

The cache is a folder with one csv per site, with a row per scene holding the mean
NDVI of the clear grid cells, their number and the number of grid cells within the
radius, and a json per site that records the search, which marks the site as done.
Sites are independent, so prefetch.py --part i/n splits them over jobs of a cluster
array, and the scenes of a site are read in parallel threads.

In the cube a day takes the mean NDVI of its scenes with at least min_clear of the
grid cells clear, weighted by their clear grid cells, and is NaN on days without
such a scene. NDVI_QC is the clear fraction of the grid cells of those scenes.
"""
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS = "https://planetarycomputer.microsoft.com/api/sas/v1/token/landsat-c2-l2"
GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}
# QA_PIXEL bits 0 to 5: fill, dilated cloud, cirrus, cloud, cloud shadow and snow.
QA_NOT_CLEAR = 0b111111
SR_SCALE, SR_OFFSET = 2.75e-5, -0.2
_token = {"value": None, "time": 0.0}
_token_lock = threading.Lock()


def sign(href):
    # The read token lasts about an hour, so renew it well before then.
    with _token_lock:
        if time.time() - _token["time"] > 1800:
            _token["value"] = json.loads(fu.fetch_url(SAS))["token"]
            _token["time"] = time.time()
        return f"/vsicurl/{href}?{_token['value']}"


def scenes(lat, lon, info):
    """Scenes of the tiers in the json at the site, with their cloud cover filter."""
    body = {"collections": ["landsat-c2-l2"], "limit": 1000,
            "intersects": {"type": "Point", "coordinates": [lon, lat]},
            "datetime": f"{info.get('first', '1982-01-01')}/"
                        f"{info.get('last', '2030-12-31')}",
            "query": {"eo:cloud_cover": {"lte": info.get("max_cloud", 80)}}}
    tiers = set(info.get("tiers", ["T1"]))
    items = []
    while body is not None:
        page = json.loads(fu.fetch_url(STAC, json.dumps(body).encode(),
                                       {"Content-Type": "application/json"}))
        items += [it for it in page["features"]
                  if it["properties"].get("landsat:collection_category") in tiers]
        nxt = [ln for ln in page.get("links", []) if ln["rel"] == "next"]
        body = nxt[0]["body"] if nxt else None
    return sorted(items, key=lambda it: it["properties"]["datetime"])


def scene_row(item, lat, lon):
    """Mean NDVI of the clear grid cells within the radius of the site in a scene."""
    import rasterio
    from rasterio.warp import transform
    from rasterio.windows import from_bounds
    bands = {}
    with rasterio.Env(**GDAL_ENV):
        for band in ("red", "nir08", "qa_pixel"):
            with rasterio.open(sign(item["assets"][band]["href"])) as src:
                xs, ys = transform("EPSG:4326", src.crs, [lon], [lat])
                r = fu.RADIUS_M + 30
                win = from_bounds(xs[0] - r, ys[0] - r, xs[0] + r, ys[0] + r,
                                  src.transform).round_offsets().round_lengths()
                bands[band] = src.read(1, window=win, boundless=True, fill_value=0)
                tr = src.window_transform(win)
    rows, cols = np.indices(bands["red"].shape)
    x, y = tr * (cols + 0.5, rows + 0.5)
    inside = (x - xs[0]) ** 2 + (y - ys[0]) ** 2 <= fu.RADIUS_M ** 2
    red = bands["red"] * SR_SCALE + SR_OFFSET
    nir = bands["nir08"] * SR_SCALE + SR_OFFSET
    clear = inside & ((bands["qa_pixel"] & QA_NOT_CLEAR) == 0) & \
        (red > 0) & (red <= 1) & (nir > 0) & (nir <= 1)
    ndvi = (nir[clear] - red[clear]) / (nir[clear] + red[clear])
    p = item["properties"]
    return {"date": p["datetime"][:10], "scene": item["id"], "platform": p["platform"],
            "cloud_cover": p.get("eo:cloud_cover"),
            "ndvi": float(ndvi.mean()) if clear.any() else np.nan,
            "n_clear": int(clear.sum()), "n_total": int(inside.sum())}


def paths(info, site):
    base = os.path.join(info["cache_path"], site)
    return f"{base}.csv", f"{base}.json"


def is_done(info, site, lat, lon):
    """The search of the site was completed at these coordinates."""
    _, meta = paths(info, site)
    if not os.path.exists(meta):
        return False
    m = json.load(open(meta))
    return max(abs(m["latitude"] - lat), abs(m["longitude"] - lon)) <= fu.COORD_TOL_DEG


def fetch_site(site, lat, lon, info, workers=8):
    """Read every scene of a site and write its csv and json."""
    items = scenes(lat, lon, info)
    failed = []

    def one(item):
        try:
            return scene_row(item, lat, lon)
        except Exception as e:
            failed.append(item["id"])
            logger.warning(f"{site}: could not read {item['id']}: {e}")
            return None

    with ThreadPoolExecutor(workers) as ex:
        rows = [r for r in ex.map(one, items) if r is not None]
    csv, meta = paths(info, site)
    df = pd.DataFrame(rows, columns=["date", "scene", "platform", "cloud_cover", "ndvi",
                                     "n_clear", "n_total"])
    df.insert(0, "site", site)
    df["latitude"], df["longitude"] = lat, lon
    df["fetched_utc"] = fu.now_utc()
    df.round(5).to_csv(f"{csv}.tmp", index=False)
    os.replace(f"{csv}.tmp", csv)
    # A site with failed scenes is not marked as done, so a later run reads it again.
    if not failed:
        with open(meta, "w") as f:
            json.dump({"site": site, "latitude": lat, "longitude": lon,
                       "n_scenes": len(items), "first": info.get("first", ""),
                       "last": info.get("last", ""), "max_cloud": info.get("max_cloud", 80),
                       "fetched_utc": fu.now_utc()}, f)
    return len(items), len(failed)


def prefetch(sites, lats, lons, info, workers=None, part="", **kwargs):
    """Fetch every site of this part that is not done yet, one after the other."""
    from extractors.clms_dekad import select_part
    os.makedirs(info["cache_path"], exist_ok=True)
    todo = select_part([s for s in zip(sites, lats, lons)
                        if not is_done(info, *s)], part)
    n_workers = workers or info.get("n_workers", 8)
    print(f"{len(todo)} sites to fetch into {info['cache_path']} with {n_workers} "
          "threads per site", flush=True)
    for site, lat, lon in todo:
        t = time.time()
        try:
            n, bad = fetch_site(site, lat, lon, info, n_workers)
        except Exception as e:
            print(f"failed {site}: {e}", flush=True)
            continue
        print(f"done {site}: {n} scenes, {bad} failed, {time.time() - t:.0f} s",
              flush=True)


def report(sites, lats, lons, info):
    done = [s for s in zip(sites, lats, lons) if is_done(info, *s)]
    print(f"{len(done)} of {len(sites)} sites done in {info['cache_path']}")


def site_rows(info, site, lat, lon):
    csv, _ = paths(info, site)
    rows = fu.read_cache(csv)
    if len(rows) == 0:
        return rows
    near = (np.abs(rows["latitude"] - lat) <= fu.COORD_TOL_DEG) & \
           (np.abs(rows["longitude"] - lon) <= fu.COORD_TOL_DEG)
    return rows[near]


def daily_series(rows, days, min_clear):
    """Clear grid cell weighted NDVI and clear fraction of the days with scenes."""
    ndvi = pd.Series(np.nan, index=days)
    qc = pd.Series(np.nan, index=days)
    if len(rows):
        r = rows[(rows["n_total"] > 0) & (rows["n_clear"] / rows["n_total"] >= min_clear)
                 & rows["ndvi"].notna()].copy()
        r["w"] = r["ndvi"] * r["n_clear"]
        g = r.groupby(pd.to_datetime(r["date"]))[["w", "n_clear", "n_total"]].sum()
        g = g[g.index.isin(days)]
        ndvi[g.index] = g["w"] / g["n_clear"]
        qc[g.index] = g["n_clear"] / g["n_total"]
    return ndvi.values, qc.values


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    per_day = 24 if bxtr.temporal_resolution == "hourly" else 1
    days = pd.DatetimeIndex(date_[::per_day]).normalize()

    rows = site_rows(info, bxtr.site, bxtr.lat, bxtr.lon)
    if len(rows) == 0:
        logger.warning(f"::MISSING:: {dataset} has no Landsat scenes for {bxtr.site} "
                       f"in {info['cache_path']}. NaN values will be set.")
    ndvi, qc = daily_series(rows, days, info.get("min_clear", 0.5))
    series = {"NDVI": ndvi, "NDVI_QC": qc}

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        values = np.repeat(series[bxtr.vars[tar_name]["sourceVariableName"]], per_day)
        data = xr.Dataset({tar_name: xr.DataArray(data=values, dims=['time'],
                                                  coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    src_dataset = fu.with_fetch_attrs(src_dataset, rows, info)
    for v in src_dataset.data_vars:
        src_dataset[v].attrs["n_days_with_scenes"] = int(np.isfinite(ndvi).sum())
    return src_dataset


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
