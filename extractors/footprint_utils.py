#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared tools for the footprint extractors, which summarise a gridded product over the
pixels whose centres lie within RADIUS_M great-circle distance of a site.

Remote products are slow or rate limited, so every footprint extractor keeps a per-
site cache in a csv file given by cache_path in its json. The cache is filled once,
ideally for all sites with prefetch.py, and extract only reads it. When fetch_missing
is true in the json, a site missing from the cache is fetched inside extract under a
lock, so parallel workers never fetch the same site twice. An empty cache_path means
the source is local and is read on every call.

A provider module defines fetch_site(site, lat, lon, info) returning a list of cache
rows, components(rows, info) returning {component: {year: value}}, and optionally
prefetch(sites, lats, lons, info, **kwargs) for a faster batch download.
"""
import fcntl
import os
import time
import urllib.request
from contextlib import contextmanager

import numpy as np
import pandas as pd
import xarray as xr

import logging
logger = logging.getLogger(__name__)

RADIUS_M = 1000.0
EARTH_RADIUS_M = 6371008.8
# Cached rows whose coordinates differ more than this from the site are reported.
COORD_TOL_DEG = 1e-3


def haversine(lat0, lon0, lat, lon):
    p0, p = np.radians(lat0), np.radians(lat)
    a = (np.sin((p - p0) / 2) ** 2
         + np.cos(p0) * np.cos(p) * np.sin(np.radians(lon - lon0) / 2) ** 2)
    return 2 * EARTH_RADIUS_M * np.arcsin(np.sqrt(a))


def half_box(lat, radius=RADIUS_M):
    """Half sizes in degrees of latitude and longitude of a box around the radius."""
    half_lat = np.degrees(radius / EARTH_RADIUS_M)
    return half_lat, half_lat / np.cos(np.radians(lat))


def geo_window(t, width, height, lat, lon, radius=RADIUS_M):
    """
    Pixel window and footprint mask on a north-up geographic grid with affine t.
    Returns col0, row0, the window width and height, and the boolean mask of pixels
    whose centres lie within radius of the site.
    """
    half_lat, half_lon = half_box(lat, radius)
    col0 = max(int(np.floor((lon - half_lon - t.c) / t.a)), 0)
    col1 = min(int(np.ceil((lon + half_lon - t.c) / t.a)), width)
    row0 = max(int(np.floor((lat + half_lat - t.f) / t.e)), 0)
    row1 = min(int(np.ceil((lat - half_lat - t.f) / t.e)), height)
    plons = t.c + (np.arange(col0, col1) + 0.5) * t.a
    plats = t.f + (np.arange(row0, row1) + 0.5) * t.e
    lon2d, lat2d = np.meshgrid(plons, plats)
    inside = haversine(lat, lon, lat2d, lon2d) <= radius
    return col0, row0, col1 - col0, row1 - row0, inside


def fetch_url(url, data=None, headers=None, retries=6, wait=3, timeout=180):
    """Bytes of url, retrying with a growing wait on network errors."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data, headers or {})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except OSError:
            time.sleep(wait * (attempt + 1))
    raise RuntimeError(f"failed {url}")


@contextmanager
def file_lock(path):
    """
    Exclusive lock held on path + '.lock'. Every call opens its own file handle, so
    flock also keeps threads of one process apart.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path + ".lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def read_cache(path):
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=["site"])
    return pd.read_csv(path)


def append_cache(path, rows):
    """Add rows to the cache csv atomically, keeping the columns of both."""
    if not rows:
        return
    new = pd.DataFrame(rows)
    with file_lock(path):
        old = read_cache(path)
        out = pd.concat([old, new], ignore_index=True) if len(old) else new
        tmp = f"{path}.tmp{os.getpid()}"
        out.round(5).to_csv(tmp, index=False)
        os.replace(tmp, path)


def cached_sites(path):
    return set(read_cache(path)["site"].astype(str))


def site_rows(cache, site):
    """Cached rows of one site, matched without case, keeping the last fetch."""
    if len(cache) == 0:
        return cache
    return cache[cache["site"].astype(str).str.lower() == site.lower()]


def rows_for_site(provider, info, site, lat, lon):
    """
    Cache rows for a site. A missing site is fetched only when fetch_missing is true,
    and the cache is checked again under the lock in case another worker fetched it.
    """
    path = info.get("cache_path", "")
    if not path:
        return pd.DataFrame(provider.fetch_site(site, lat, lon, info))
    rows = site_rows(read_cache(path), site)
    if len(rows) == 0 and info.get("fetch_missing", False):
        with file_lock(f"{path}.{site}"):
            rows = site_rows(read_cache(path), site)
            if len(rows) == 0:
                logger.info(f"fetching {site} for {path}")
                try:
                    fetched = provider.fetch_site(site, lat, lon, info)
                except Exception as e:
                    logger.warning(f"could not fetch {site} for {path}: {e}")
                    fetched = []
                for r in fetched:
                    r.setdefault("latitude", lat)
                    r.setdefault("longitude", lon)
                append_cache(path, fetched)
                rows = site_rows(read_cache(path), site)
    if len(rows) and {"latitude", "longitude"} <= set(rows.columns):
        dlat = np.nanmax(np.abs(rows["latitude"].astype(float) - lat))
        dlon = np.nanmax(np.abs(rows["longitude"].astype(float) - lon))
        if max(dlat, dlon) > COORD_TOL_DEG:
            logger.warning(f"cached rows of {site} in {path} were fetched at a "
                           f"location {max(dlat, dlon):.4f} deg away from the site")
    return rows


def missing_sites(info, sites, lats, lons):
    """Sites, latitudes and longitudes not yet in the cache of info."""
    done = {s.lower() for s in cached_sites(info["cache_path"])}
    keep = [i for i, s in enumerate(sites) if s.lower() not in done]
    return [sites[i] for i in keep], [lats[i] for i in keep], [lons[i] for i in keep]


def prefetch_sites(provider, info, sites, lats, lons, workers=None):
    """Fetch every site missing from the cache with a pool of threads."""
    from concurrent.futures import ThreadPoolExecutor
    todo = list(zip(*missing_sites(info, sites, lats, lons)))
    n_workers = workers or info.get("n_workers", 4)
    print(f"{len(todo)} of {len(sites)} sites to fetch into {info['cache_path']} "
          f"with {n_workers} workers", flush=True)

    def one(args):
        site, lat, lon = args
        try:
            rows = provider.fetch_site(site, lat, lon, info)
        except Exception as e:
            print(f"failed {site}: {e}", flush=True)
            return
        for r in rows:
            r.setdefault("latitude", lat)
            r.setdefault("longitude", lon)
        append_cache(info["cache_path"], rows)
        print(f"done {site}", flush=True)

    with ThreadPoolExecutor(n_workers) as ex:
        list(ex.map(one, todo))


def site_pft(pft_path, site):
    """PFT of the site from a csv with Site_ID and PFT columns, or None."""
    if not pft_path or not os.path.exists(pft_path):
        return None
    df = pd.read_csv(pft_path)
    sel = df.loc[df["Site_ID"].str.lower() == site.lower(), "PFT"]
    return str(sel.values[0]).strip() if len(sel) else None


def insitu_fractions(veg_path, site, pft, tree_pfts):
    """
    In-situ vegetation and tree fraction of the site from FLUXNET_veg_frac.csv, or
    NaN. Vegetation is 1 - frac_soil, or overstory + understory when soil is missing.
    Tree is frac_overstory at forest and savanna PFTs and 0 elsewhere. Several rows
    of a site are averaged.
    """
    if not veg_path or not os.path.exists(veg_path):
        return np.nan, np.nan
    c = pd.read_csv(veg_path, encoding="utf-8-sig").dropna(subset=["frac_overstory"])
    c = c[c["site_id"].str.lower() == site.lower()]
    if len(c) == 0:
        return np.nan, np.nan
    over = c["frac_overstory"] / 100
    under = c["frac_understory"].fillna(0) / 100
    veg = np.where(c["frac_soil"].notna(), 1 - c["frac_soil"] / 100,
                   np.minimum(1.0, over + under))
    veg, tree = float(np.mean(veg)), float(over.mean())
    if pft not in tree_pfts:
        tree = 0.0
    return veg, tree


def collapse(series, placement):
    """For constant placement, the mean over the years under the key None."""
    if placement == "constant":
        vals = [v for v in series.values() if np.isfinite(v)]
        return {None: float(np.mean(vals)) if vals else np.nan}
    return series


def merge_insitu(comps, ins_veg, ins_tree, placement):
    """
    In-situ merged tree and vegetation per year: in situ where the site has it,
    otherwise the product. Tree is capped at vegetation, and both are clipped to 0-1.
    """
    veg_p = collapse(comps.get("veg", {}), placement)
    tree_p = collapse(comps.get("tree", {}), placement)
    has = np.isfinite(ins_veg)
    keys = [None] if has and placement == "constant" else sorted(
        set(veg_p) | set(tree_p), key=lambda k: -1 if k is None else k)
    veg_m, tree_m = {}, {}
    for k in keys:
        veg = ins_veg if has else veg_p.get(k, np.nan)
        tree = ins_tree if has else tree_p.get(k, np.nan)
        veg_m[k] = float(np.clip(veg, 0, 1))
        tree_m[k] = float(np.clip(np.minimum(tree, veg), 0, 1))
    return {"veg": veg_m, "tree": tree_m}


def place(date_, series, placement):
    """
    Values over the time steps of date_. constant: one value broadcast over all time
    steps. annual: each year's value on every time step of that year. year_end: each
    year's value on the last time step of that year. Time steps without a value are
    NaN.
    """
    out = np.full(len(date_), np.nan)
    if placement == "constant":
        if None in series:
            out[:] = series[None]
        return out
    years = date_.astype("M8[Y]").astype(int) + 1970
    for yr, val in series.items():
        in_year = np.flatnonzero(years == int(yr))
        if len(in_year) == 0:
            continue
        if placement == "annual":
            out[in_year] = val
        else:
            out[in_year[-1]] = val
    return out


def extract_footprint(provider, dataset, site_info, config):
    """
    Generic extract for footprint products. Each variable in the json names the
    component it takes and whether it merges the product with the in-situ fractions.
    """
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    placement = info.get("time_placement", "constant")
    date_ = bxtr.get_date_vec()

    rows = rows_for_site(provider, info, bxtr.site, bxtr.lat, bxtr.lon)
    comps = provider.components(rows, info) if len(rows) else {}
    if len(rows) == 0:
        logger.warning(f"::MISSING:: {dataset} has no cached data for {bxtr.site} in "
                       f"{info.get('cache_path', '')}. NaN values will be set. Run "
                       "prefetch.py or set fetch_missing to true.")
    merged = None
    if any(bxtr.vars[v]["merge_insitu"] for v in bxtr.vars_list):
        pft = site_pft(info.get("pft_path", ""), bxtr.site)
        ins_veg, ins_tree = insitu_fractions(info.get("insitu_path", ""), bxtr.site,
                                             pft, info.get("tree_pfts", []))
        merged = merge_insitu(comps, ins_veg, ins_tree, placement)

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        comp = bxtr.vars[tar_name]["component"]
        if bxtr.vars[tar_name]["merge_insitu"]:
            series = merged.get(comp, {})
        else:
            series = collapse(comps.get(comp, {}), placement)
        data = xr.Dataset({tar_name: xr.DataArray(
            data=place(date_, series, placement), dims=['time'],
            coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    return bxtr.merge_and_format(src_data)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
