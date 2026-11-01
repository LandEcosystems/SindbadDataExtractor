#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vegetation fraction within 1 km of the site from Copernicus FCOVER 300 m V1.

FCOVER is the fraction of ground covered by green vegetation every 10 days, so the
vegetation fraction is the per grid cell maximum over the 36 dekads of the year,
averaged over the 300 m grid cells whose centres lie within 1 km of the site. Windows
are read from the cloud optimised GeoTIFFs on the Copernicus Data Space S3 endpoint
with the keys of the cdse profile in ~/.aws/credentials. prefetch opens each dekad
file once and reads the windows of all missing sites from it. The dekad file list is
saved next to the cache. The cache holds one row per site, the same format as
fcover300_2015/site_means.csv of add_vegfrac_fcover.py.

The windows read from a dekad file are kept in a part file of that dekad,
<cache>_parts_<year>/<dekad>_<tag>.npz, with the sites, their coordinates and the
values of each window. A site gets its row once every dekad of the year has its
window at the coordinates of the site. A dekad that fails to read is logged and read
again by the next run, which reads only the dekads and sites that no part file holds.

FCOVER has no tree layer, so only vegetation is provided.
"""
import glob
import hashlib
import json
import os
import sys
import urllib.parse

import numpy as np

import extractors.footprint_utils as fu

ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
GDAL_ENV = {"AWS_PROFILE": "cdse", "AWS_S3_ENDPOINT": "eodata.dataspace.copernicus.eu",
            "AWS_HTTPS": "YES", "AWS_VIRTUAL_HOSTING": "FALSE",
            "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}


def year_of(info):
    return int(info.get("years", [2015])[0])


def dekad_paths(info):
    """S3 paths of the FCOVER layer for every dekad of the year, listed once."""
    year = year_of(info)
    listing = f"{os.path.splitext(info['cache_path'])[0]}_dekads_{year}.json"
    with fu.file_lock(listing):
        if os.path.exists(listing):
            return json.load(open(listing))
        flt = (f"contains(Name,'c_gls_FCOVER300_{year}') and endswith(Name,'_cog') "
               "and contains(Name,'V1.0.1')")
        url = (f"{ODATA}?$filter={urllib.parse.quote(flt)}&$top=100"
               "&$orderby=ContentDate/Start%20asc")
        out = []
        for p in json.loads(fu.fetch_url(url))["value"]:
            stamp = p["Name"].split("_")[3]
            tif = f"c_gls_FCOVER300-FCOVER_{stamp}_GLOBE_PROBAV_V1.0.1.tiff"
            out.append([stamp, f"/vsis3{p['S3Path']}/{tif}"])
        with open(listing, "w") as f:
            json.dump(out, f)
        return out


def dekad_values(path, sites, lats, lons):
    """FCOVER of the grid cells within the radius of every site for one dekad."""
    import rasterio
    from rasterio.windows import Window
    vals = []
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(path) as src:
            scale, offset, nodata = src.scales[0], src.offsets[0], src.nodata
            for lat, lon in zip(lats, lons):
                col0, row0, w, h, inside = fu.geo_window(src.transform, src.width,
                                                         src.height, lat, lon)
                raw = src.read(1, window=Window(col0, row0, w, h)).astype(float)
                bad = raw == nodata if nodata is not None else np.zeros_like(raw, bool)
                v = raw * scale + offset
                v[bad | (v < 0) | (v > 1)] = np.nan
                vals.append(v[inside])
    return vals


def site_record(site, lat, lon, arr):
    """Peak over the dekads per grid cell, then its mean over the footprint."""
    n_obs = np.isfinite(arr).sum(axis=0)
    peak = np.where(n_obs > 0, np.nanmax(np.where(np.isfinite(arr), arr, -1), axis=0),
                    np.nan)
    return {"site": site, "n_grids": arr.shape[1],
            "n_grids_valid": int((n_obs > 0).sum()),
            "median_dekads_valid": float(np.median(n_obs)),
            "veg": np.nanmean(peak) if np.isfinite(peak).any() else np.nan,
            "mean_fcover": np.nanmean(arr) if np.isfinite(arr).any() else np.nan,
            "latitude": lat, "longitude": lon}


def fetch_sites(sites, lats, lons, info, workers=None):
    from concurrent.futures import ThreadPoolExecutor
    paths = [p for _, p in dekad_paths(info)]
    with ThreadPoolExecutor(workers or info.get("n_workers", 4)) as ex:
        per_dekad = list(ex.map(lambda p: dekad_values(p, sites, lats, lons), paths))
    return [site_record(s, la, lo, np.vstack([d[i] for d in per_dekad]))
            for i, (s, la, lo) in enumerate(zip(sites, lats, lons))]


def fetch_site(site, lat, lon, info):
    return fetch_sites([site], [lat], [lon], info, workers=1)


def part_dir(info):
    return f"{os.path.splitext(info['cache_path'])[0]}_parts_{year_of(info)}"


def coord_key(site, lat, lon):
    return site.lower(), round(float(lat), 6), round(float(lon), 6)


def write_part(info, stamp, sites, lats, lons, vals):
    """The windows of one dekad file at a batch of sites, as one part file."""
    keys = "\n".join(f"{s}|{la}|{lo}" for s, la, lo in zip(sites, lats, lons))
    tag = hashlib.sha1(keys.encode()).hexdigest()[:10]
    path = os.path.join(part_dir(info), f"{stamp}_{tag}.npz")
    os.makedirs(part_dir(info), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "wb") as f:
        np.savez(f, site=np.array(sites, dtype=str), latitude=np.array(lats, float),
                 longitude=np.array(lons, float),
                 offsets=np.cumsum([0] + [len(v) for v in vals]),
                 values=np.concatenate(vals) if vals else np.zeros(0))
    os.replace(tmp, path)


def read_parts(info):
    """{dekad: {(site, lat, lon): window values}} of the part files."""
    out = {}
    for path in sorted(glob.glob(os.path.join(part_dir(info), "*.npz"))):
        stamp = os.path.basename(path).split("_")[0]
        with np.load(path) as p:
            off, val = p["offsets"], p["values"]
            have = out.setdefault(stamp, {})
            for i, key in enumerate(zip(p["site"], p["latitude"], p["longitude"])):
                have[coord_key(*key)] = val[off[i]:off[i + 1]]
    return out


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    from concurrent.futures import ThreadPoolExecutor
    sites, lats, lons = fu.missing_sites(info, sites, lats, lons)
    if not sites:
        print(f"all sites are in {info['cache_path']}")
        return
    dekads = dekad_paths(info)
    keys = [coord_key(*s) for s in zip(sites, lats, lons)]
    have = read_parts(info)
    jobs = []
    for stamp, path in dekads:
        todo = [i for i, k in enumerate(keys) if k not in have.get(stamp, {})]
        if todo:
            jobs.append((stamp, path, todo))
    print(f"{len(sites)} sites to fetch into {info['cache_path']}, {len(jobs)} of "
          f"{len(dekads)} dekads to read", flush=True)

    def one(job):
        stamp, path, todo = job
        s, la, lo = ([x[i] for i in todo] for x in (sites, lats, lons))
        try:
            vals = dekad_values(path, s, la, lo)
        except Exception as e:
            print(f"failed {stamp}: {type(e).__name__}: {e}", flush=True)
            return
        write_part(info, stamp, s, la, lo, vals)
        print(f"done {stamp}, {len(todo)} sites", flush=True)

    with ThreadPoolExecutor(workers or info.get("n_workers", 4)) as ex:
        list(ex.map(one, jobs))
    have = read_parts(info)
    rows, waiting = [], 0
    for (s, la, lo), k in zip(zip(sites, lats, lons), keys):
        per = [have.get(stamp, {}).get(k) for stamp, _ in dekads]
        if any(v is None for v in per):
            waiting += 1
            continue
        rows.append(site_record(s, la, lo, np.vstack(per)))
    fu.append_cache(info["cache_path"], rows)
    print(f"{len(rows)} sites added to {info['cache_path']}, {waiting} wait for dekads "
          "that failed", flush=True)


def components(rows, info):
    return {"veg": {year_of(info): float(rows.iloc[-1]["veg"])}}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
