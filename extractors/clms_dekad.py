#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ten-day series of a Copernicus Land Monitoring Service layer, such as NDVI or LAI at
300 m or 1 km, averaged over the grid cells whose centres lie within 1 km of the
site.

The json lists the series of the product, each with the OData filter of its
products in the Copernicus Data Space catalogue and the first and last day it is
used for, so that a product made of several sensors or versions, such as LAI 300 m
from PROBA-V and then Sentinel-3 OLCI, is one series in time. Products released in
stages, as RT0 to RT6 of LAI 300 m OLCI, keep the latest stage of every dekad. The
layer is read from the cloud optimised GeoTIFF of each product on the Copernicus
Data Space S3 endpoint, with the keys of the cdse profile in ~/.aws/credentials.

prefetch lists the products once, in <cache_path>/listing.json, and then opens the
file of every dekad once and reads the windows of all missing sites from it. The
cache is a folder with one csv per year, with a row per dekad and site holding the
mean of the valid grid cells, their number and the number of grid cells within
the radius. Years are independent, so prefetch.py --part i/n splits them over jobs
of a cluster array. A grid cell is valid when it is not nodata and its value lies
within valid_range of the json, which drops the flag values above the range.

In the cube every day of a dekad, or every hour of those days, takes the value of
the dekad. <layer>_QC is the fraction of the grid cells within the radius that are
valid.
"""
import calendar
import glob
import json
import os
import re
import urllib.parse
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
GDAL_ENV = {"AWS_PROFILE": "cdse", "AWS_S3_ENDPOINT": "eodata.dataspace.copernicus.eu",
            "AWS_HTTPS": "YES", "AWS_VIRTUAL_HOSTING": "FALSE",
            "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}
STAMP = re.compile(r"_(\d{12})_")
STAGE = re.compile(r"-RT(\d)_")
PAGE = 1000


def fetch_json(url):
    return json.loads(fu.fetch_url(url))


def catalogue(flt):
    """Name and S3 path of every product matching an OData filter, by date."""
    out, skip = [], 0
    while True:
        url = (f"{ODATA}?$filter={urllib.parse.quote(flt)}&$top={PAGE}&$skip={skip}"
               "&$orderby=ContentDate/Start%20asc")
        page = fetch_json(url)["value"]
        out += [(p["Name"], p["S3Path"], p["Id"]) for p in page]
        if len(page) < PAGE:
            return out
        skip += PAGE


def layer_file(pid, layer):
    """Name of the GeoTIFF of a layer inside a product, from its node listing."""
    # Staged products name the stage after the layer, as -LAI-RT6_.
    pattern = re.compile(rf"-{layer}(-RT\d)?_")
    for node in fetch_json(f"{ODATA}({pid})/Nodes")["result"]:
        for sub in fetch_json(node["Nodes"]["uri"])["result"]:
            if pattern.search(sub["Name"]) and sub["Name"].endswith(".tiff"):
                return sub["Name"]
    return None


def dekad_bounds(stamp):
    """First and last day of the dekad holding the date of a product stamp."""
    d = pd.Timestamp(stamp[:8])
    first = 1 if d.day <= 10 else 11 if d.day <= 20 else 21
    last = 10 if first == 1 else 20 if first == 11 else \
        calendar.monthrange(d.year, d.month)[1]
    return d.replace(day=first), d.replace(day=last)


def listing(info, workers=8):
    """
    Dekad, first day and S3 path of the layer file of every product of the series,
    made once and kept in <cache_path>/listing.json.
    """
    path = os.path.join(info["cache_path"], "listing.json")
    os.makedirs(info["cache_path"], exist_ok=True)
    with fu.file_lock(path):
        if os.path.exists(path):
            return json.load(open(path))
        chosen = {}
        for s in info["series"]:
            first, last = pd.Timestamp(s["first"]), pd.Timestamp(s["last"])
            for name, s3, pid in catalogue(s["odata"]):
                m = STAMP.search(name)
                if not m:
                    continue
                start = dekad_bounds(m.group(1))[0]
                if not first <= start <= last:
                    continue
                key = start.strftime("%Y%m%d")
                stage = int(STAGE.search(name).group(1)) if STAGE.search(name) else 0
                if key not in chosen or stage >= chosen[key][0]:
                    chosen[key] = (stage, name, s3, pid)
        keys = sorted(chosen)
        with ThreadPoolExecutor(workers) as ex:
            files = list(ex.map(lambda k: layer_file(chosen[k][3], info["layer"]), keys))
        out = [[k, chosen[k][1], f"/vsis3{chosen[k][2]}/{f}"]
               for k, f in zip(keys, files) if f]
        missing = [chosen[k][1] for k, f in zip(keys, files) if not f]
        if missing:
            print(f"no {info['layer']} file in {len(missing)} products: {missing[:5]}")
        with open(path, "w") as f:
            json.dump(out, f)
        return out


def dekad_rows(dekad, product, path, sites, lats, lons, info):
    """Cache rows of one dekad for the given sites, read from one open file."""
    import rasterio
    from rasterio.windows import Window
    lo, hi = info.get("valid_range", [-np.inf, np.inf])
    rows = []
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(path) as src:
            scale, offset, nodata = src.scales[0], src.offsets[0], src.nodata
            for site, lat, lon in zip(sites, lats, lons):
                col0, row0, w, h, inside = fu.geo_window(src.transform, src.width,
                                                         src.height, lat, lon)
                raw = src.read(1, window=Window(col0, row0, w, h))[inside]
                v = raw.astype(float) * scale + offset
                ok = (v >= lo) & (v <= hi)
                if nodata is not None:
                    ok &= raw != nodata
                rows.append({"dekad": dekad, "site": site, "latitude": lat,
                             "longitude": lon,
                             "value": float(v[ok].mean()) if ok.any() else np.nan,
                             "n_valid": int(ok.sum()), "n_total": int(inside.sum()),
                             "product": product})
    return rows


def year_path(info, year):
    return os.path.join(info["cache_path"], f"{year}.csv")


def todo_of_year(info, dekads, sites, lats, lons):
    """Dekads of one year with the sites, latitudes and longitudes still missing."""
    cache = fu.read_cache(year_path(info, dekads[0][0][:4]))
    done = set()
    if len(cache):
        done = set(zip(cache["dekad"].astype(str), cache["site"].astype(str).str.lower()))
    out = []
    for dekad, product, path in dekads:
        keep = [i for i, s in enumerate(sites) if (dekad, s.lower()) not in done]
        if keep:
            out.append((dekad, product, path, [sites[i] for i in keep],
                        [lats[i] for i in keep], [lons[i] for i in keep]))
    return out


def select_part(items, part):
    """Items of part i of n, as i/n with i from 1 to n."""
    if not part:
        return items
    i, n = (int(x) for x in part.split("/"))
    if not 1 <= i <= n:
        raise ValueError(f"part {part} is not i/n with 1 <= i <= n")
    return items[i - 1::n]


def prefetch(sites, lats, lons, info, workers=None, part="", **kwargs):
    """
    Read every missing site and dekad of the years of this part, one thread per
    open dekad file, writing each dekad to the csv of its year as it finishes.
    """
    n_workers = workers or info.get("n_workers", 8)
    by_year = {}
    for d in listing(info, n_workers):
        by_year.setdefault(d[0][:4], []).append(d)
    years = select_part(sorted(by_year), part)
    print(f"{info['layer']}: {len(years)} of {len(by_year)} years in this part, "
          f"{len(sites)} sites, {n_workers} threads", flush=True)

    def run(job):
        dekad, product, path, s, la, lo = job
        try:
            rows = dekad_rows(dekad, product, path, s, la, lo, info)
        except Exception as e:
            print(f"failed {dekad} {path}: {e}", flush=True)
            return
        fu.append_cache(year_path(info, dekad[:4]), rows)
        print(f"done {dekad}, {len(rows)} sites", flush=True)

    for year in years:
        jobs = todo_of_year(info, by_year[year], sites, lats, lons)
        print(f"{year}: {len(jobs)} of {len(by_year[year])} dekads to read", flush=True)
        with ThreadPoolExecutor(n_workers) as ex:
            list(ex.map(run, jobs))


def report(sites, lats, lons, info):
    """Print how many dekads of the listing are cached for all the sites."""
    path = os.path.join(info["cache_path"], "listing.json")
    if not os.path.exists(path):
        print(f"no listing in {info['cache_path']}, nothing cached")
        return
    by_year = {}
    for d in json.load(open(path)):
        by_year.setdefault(d[0][:4], []).append(d)
    open_years = []
    for year, dekads in sorted(by_year.items()):
        if todo_of_year(info, dekads, sites, lats, lons):
            open_years.append(year)
    print(f"{info['layer']}: {len(by_year) - len(open_years)} of {len(by_year)} years "
          f"complete for {len(sites)} sites in {info['cache_path']}")
    if open_years:
        print("    incomplete years: " + " ".join(open_years))


def site_rows(info, site, lat, lon):
    """Cached rows of a site in all years, at its coordinates."""
    parts = []
    for path in sorted(glob.glob(os.path.join(info["cache_path"], "[0-9]" * 4 + ".csv"))):
        rows = fu.site_rows(fu.read_cache(path), site)
        if len(rows):
            near = (np.abs(rows["latitude"] - lat) <= fu.COORD_TOL_DEG) & \
                   (np.abs(rows["longitude"] - lon) <= fu.COORD_TOL_DEG)
            parts.append(rows[near])
    if not parts:
        return pd.DataFrame(columns=["dekad", "value", "n_valid", "n_total"])
    return pd.concat(parts).drop_duplicates("dekad", keep="last")


def daily_series(rows, days):
    """Value and valid fraction of every day, from the dekad that holds it."""
    value = pd.Series(np.nan, index=days)
    qc = pd.Series(np.nan, index=days)
    for r in rows.itertuples():
        first, last = dekad_bounds(str(r.dekad))
        value[first:last] = r.value
        qc[first:last] = r.n_valid / r.n_total if r.n_total else np.nan
    return value.values, qc.values


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
        logger.warning(f"::MISSING:: {dataset} has no {info['layer']} for {bxtr.site} "
                       f"in {info['cache_path']}. NaN values will be set.")
    value, qc = daily_series(rows, days)
    series = {info["layer"]: value, f"{info['layer']}_QC": qc}

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        src = bxtr.vars[tar_name]["sourceVariableName"]
        values = np.repeat(series[src], per_day)
        data = xr.Dataset({tar_name: xr.DataArray(data=values, dims=['time'],
                                                  coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    src_dataset = fu.with_fetch_attrs(src_dataset, rows, info)
    for v in src_dataset.data_vars:
        src_dataset[v].attrs["n_dekads"] = int(rows["value"].notna().sum())
    return src_dataset


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
