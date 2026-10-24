#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared tools for the footprint extractors, which summarise a gridded product over the
grid cells whose centres lie within RADIUS_M great-circle distance of a site.

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
# IGBP classes by number, whose abbreviations are the BADM PFTs.
IGBP = {1: "ENF", 2: "EBF", 3: "DNF", 4: "DBF", 5: "MF", 6: "CSH", 7: "OSH", 8: "WSA",
        9: "SAV", 10: "GRA", 11: "WET", 12: "CRO", 13: "URB", 14: "CVM", 15: "SNO",
        16: "BSV", 17: "WAT"}
IGBP_CODE = {v: k for k, v in IGBP.items()}
IGBP_NAMES = {1: "evergreen needleleaf forests", 2: "evergreen broadleaf forests",
              3: "deciduous needleleaf forests", 4: "deciduous broadleaf forests",
              5: "mixed forests", 6: "closed shrublands", 7: "open shrublands",
              8: "woody savannas", 9: "savannas", 10: "grasslands",
              11: "permanent wetlands", 12: "croplands", 13: "urban and built-up lands",
              14: "cropland and natural vegetation mosaics", 15: "permanent snow and ice",
              16: "barren", 17: "water bodies"}
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
    Grid cell window and footprint mask on a north-up geographic grid with affine t.
    Returns col0, row0, the window width and height, and the boolean mask of grid
    cells whose centres lie within radius of the site.
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


def nc_footprint(path, name, lat, lon, radius=RADIUS_M):
    """
    Raw values of variable name in a global north-up latitude longitude NetCDF over
    the box around the radius, with the footprint mask and the window row and column
    of the cell holding the site. The latitude and longitude axes are cut to the box
    and other axes are read whole. Scale factors and fill values are not applied,
    since some BGI files carry a scale factor that should divide.
    """
    import netCDF4
    from rasterio.transform import from_origin
    with netCDF4.Dataset(path) as nc:
        v = nc[name]
        v.set_auto_maskandscale(False)
        lat_dim = [d for d in v.dimensions if d in ("lat", "latitude")][0]
        lon_dim = [d for d in v.dimensions if d in ("lon", "longitude")][0]
        lats = np.asarray(nc[lat_dim][:], dtype=float)
        lons = np.asarray(nc[lon_dim][:], dtype=float)
        dlat = (lats[0] - lats[-1]) / (len(lats) - 1)
        dlon = (lons[-1] - lons[0]) / (len(lons) - 1)
        t = from_origin(lons[0] - dlon / 2, lats[0] + dlat / 2, dlon, dlat)
        col0, row0, w, h, inside = geo_window(t, len(lons), len(lats), lat, lon, radius)
        cut = {lat_dim: slice(row0, row0 + h), lon_dim: slice(col0, col0 + w)}
        data = v[tuple(cut.get(d, slice(None)) for d in v.dimensions)]
    row = int(np.floor((lat - t.f) / t.e)) - row0
    col = int(np.floor((lon - t.c) / t.a)) - col0
    return data, inside, (row, col)


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


def _fetch_one(args):
    """Rows of one site for prefetch_sites, or None when the fetch fails."""
    import importlib
    module, info, site, lat, lon = args
    try:
        rows = importlib.import_module(module).fetch_site(site, lat, lon, info)
    except Exception as e:
        print(f"failed {site}: {e}", flush=True)
        return site, None
    for r in rows:
        r.setdefault("latitude", lat)
        r.setdefault("longitude", lon)
    return site, rows


def prefetch_sites(provider, info, sites, lats, lons, workers=None, processes=False):
    """
    Fetch every site missing from the cache with a pool of threads, or of processes
    for readers such as HDF5 that are not thread safe. Rows are written to the cache
    as each site finishes.
    """
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
    todo = [(provider.__name__, info) + s
            for s in zip(*missing_sites(info, sites, lats, lons))]
    n_workers = workers or info.get("n_workers", 4)
    print(f"{len(todo)} of {len(sites)} sites to fetch into {info['cache_path']} "
          f"with {n_workers} {'processes' if processes else 'threads'}", flush=True)
    pool = ProcessPoolExecutor if processes else ThreadPoolExecutor
    with pool(n_workers) as ex:
        for site, rows in ex.map(_fetch_one, todo):
            if rows is not None:
                append_cache(info["cache_path"], rows)
                print(f"done {site}", flush=True)


def site_pft(pft_path, site):
    """PFT of the site from a csv with Site_ID and PFT columns, or None."""
    if not pft_path or not os.path.exists(pft_path):
        return None
    df = pd.read_csv(pft_path)
    sel = df.loc[df["Site_ID"].str.lower() == site.lower(), "PFT"]
    return str(sel.values[0]).strip() if len(sel) else None


def class_legend(names, pfts):
    """
    CF flag attributes of a class variable from {code: name} and {code: BADM PFT}.
    flag_meanings joins the names with underscores for spaces, and flag_pft gives
    the BADM PFT of each code, or filled for a class without one, whose PFT comes
    from the datasets in rs_fill and then pft_path through resolve_pft.
    """
    codes = sorted(names)
    return {"flag_values": codes,
            "flag_meanings": " ".join(names[c].replace(" ", "_") for c in codes),
            "flag_pft": " ".join(pfts.get(c, "filled") for c in codes)}


IGBP_LEGEND = class_legend(IGBP_NAMES, IGBP)


def badm_by_year(classes_by_year, to_badm):
    """BADM PFT of each year from {year: class}, without classes that have none."""
    out = {}
    for y, c in classes_by_year.items():
        pft = to_badm(c)
        if pft:
            out[int(y)] = pft
    return out


def fill_sources(config, keys, site, lat, lon):
    """
    BADM PFT by year of the site from the PFT datasets keys of config, as a list of
    (key, {year: PFT}) in the order of keys. A dataset that is not selected in the
    run is read from its json in config_dir_path or ancillary_config_dir_path. Each
    provider gives its PFTs through badm_series(rows).
    """
    import importlib
    import json
    out = []
    for key in keys:
        info = config["dataset"].get(key)
        if isinstance(info, str):
            from utils.shared_utils import dataset_json_path
            path = dataset_json_path(config.get("config_dir_path", ""), info,
                                     config.get("ancillary_config_dir_path", ""))
            with open(path) as f:
                info = json.load(f)
        if not info:
            logger.warning(f"PFT fill dataset {key} is not in the main json")
            continue
        prov = importlib.import_module(f"extractors.{info['extractor']}")
        rows = rows_for_site(prov, info, site, lat, lon)
        out.append((key, prov.badm_series(rows) if len(rows) else {}))
    return out


def resolve_pft(year, sources):
    """
    BADM PFT of year and the label of its source, from sources as a list of (label,
    {year: PFT}) in priority order. The first source with any PFT gives it, at the
    year itself or else its nearest year. Returns None and None without any.
    """
    for label, src in sources:
        if src:
            keys = np.array(sorted(src))
            return src[int(keys[np.argmin(np.abs(keys - year))])], label
    return None, None


def pft_mode(years, sources, pft_path, site):
    """
    Most frequent PFT over years through resolve_pft and the label of its source.
    Without any PFT in sources, the PFT of the site in the csv at pft_path, and
    otherwise undefined.
    """
    got = [resolve_pft(y, sources) for y in years]
    got = [g for g in got if g[0]]
    if got:
        top = pd.Series([p for p, _ in got]).mode().iloc[0]
        return top, next(label for p, label in got if p == top)
    pft = site_pft(pft_path, site)
    if pft in IGBP_CODE:
        return pft, "pft_path"
    return "undefined", "none"


def merge_pft(date_, pft_path, site, years, sources, placement):
    """
    IGBP number of the PFT of the site in the csv at pft_path on every time step.
    Where the csv has no row for the site or a PFT that is not an IGBP class, such as
    TBD, each of years takes its PFT from sources through resolve_pft, placed with
    placement.
    """
    pft = site_pft(pft_path, site)
    if pft in IGBP_CODE:
        return np.full(len(date_), float(IGBP_CODE[pft]))
    rs = {}
    for y in years:
        p, _ = resolve_pft(y, sources)
        if p in IGBP_CODE:
            rs[y] = float(IGBP_CODE[p])
    return place(date_, rs, placement)


def insitu_fractions(veg_path, site, pft, tree_pfts, woody_pfts=None):
    """
    In-situ vegetation, tree and woody fraction of the site from
    FLUXNET_veg_frac.csv, or NaN. Vegetation is 1 - frac_soil, or overstory +
    understory when soil is missing. Tree is frac_overstory at the PFTs in tree_pfts
    and 0 elsewhere, and woody is frac_overstory at the PFTs in woody_pfts, which
    defaults to tree_pfts, and 0 elsewhere. Several rows of a site are averaged.
    """
    out = {"veg": np.nan, "tree": np.nan, "woody": np.nan}
    if not veg_path or not os.path.exists(veg_path):
        return out
    c = pd.read_csv(veg_path, encoding="utf-8-sig").dropna(subset=["frac_overstory"])
    c = c[c["site_id"].str.lower() == site.lower()]
    if len(c) == 0:
        return out
    over = c["frac_overstory"] / 100
    under = c["frac_understory"].fillna(0) / 100
    veg = np.where(c["frac_soil"].notna(), 1 - c["frac_soil"] / 100,
                   np.minimum(1.0, over + under))
    over = float(over.mean())
    woody_pfts = tree_pfts if woody_pfts is None else woody_pfts
    out["veg"] = float(np.mean(veg))
    out["tree"] = over if pft in tree_pfts else 0.0
    out["woody"] = over if pft in woody_pfts else 0.0
    return out


def collapse(series, placement):
    """For constant placement, the mean over the years under the key None."""
    if placement == "constant":
        vals = [v for v in series.values() if np.isfinite(v)]
        return {None: float(np.mean(vals)) if vals else np.nan}
    return series


def merge_insitu(comps, ins, targets, placement, full=1.0):
    """
    In-situ merged series per year for each variable in targets, which maps a
    variable to its product component and the in-situ fraction, veg, tree or woody,
    that replaces it. The in-situ fractions of ins are used where the site has an
    in-situ vegetation fraction, otherwise the product. Tree and woody are capped at
    the merged vegetation fraction where it is known. full is the value of a
    complete cover in the units of the product, such as 100 for percent, and all
    values are clipped to 0 to full.
    """
    has = np.isfinite(ins["veg"])
    comps = {c: collapse(comps.get(c, {}), placement)
             for c in {"veg"} | {c for c, _ in targets.values()}}
    keys = [None] if has and placement == "constant" else sorted(
        set().union(*comps.values()), key=lambda k: -1 if k is None else k)
    veg = {k: ins["veg"] * full if has else comps["veg"].get(k, np.nan)
           for k in keys}
    out = {}
    for var, (comp, kind) in targets.items():
        out[var] = {}
        for k in keys:
            val = ins[kind] * full if has else comps[comp].get(k, np.nan)
            if kind != "veg":
                val = np.fmin(val, veg[k])
            out[var][k] = float(np.clip(val, 0, full))
    return out


def place(date_, series, placement):
    """
    Values over the time steps of date_. constant: one value broadcast over all time
    steps. annual: each year's value on every time step of that year. annual_extend:
    as annual, with every other year taking the value of the nearest year that has
    one. year_start and year_end: each year's value on the first or last time step of
    that year. Time steps without a value are NaN.
    """
    out = np.full(len(date_), np.nan)
    if placement == "constant":
        if None in series:
            out[:] = series[None]
        return out
    years = date_.astype("M8[Y]").astype(int) + 1970
    if placement == "annual_extend":
        known = {int(y): v for y, v in series.items() if np.isfinite(v)}
        if known:
            keys = np.array(sorted(known))
            for yr in np.unique(years):
                near = keys[np.argmin(np.abs(keys - yr))]
                out[years == yr] = known[int(near)]
        return out
    for yr, val in series.items():
        in_year = np.flatnonzero(years == int(yr))
        if len(in_year) == 0:
            continue
        if placement == "annual":
            out[in_year] = val
        elif placement == "year_start":
            out[in_year[0]] = val
        else:
            out[in_year[-1]] = val
    return out


def extract_footprint(provider, dataset, site_info, config):
    """
    Generic extract for footprint products. Each variable in the json names the
    component it takes and whether it merges the product with the in-situ fractions.
    A provider with site_attrs(rows, info) adds the attributes it returns to the
    output.
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
    # merge_insitu of a variable is true to replace its component by the in-situ
    # fraction of the same name, or names the in-situ fraction, veg, tree or woody.
    targets = {}
    for v in bxtr.vars_list:
        kind = bxtr.vars[v]["merge_insitu"]
        if kind:
            comp = bxtr.vars[v]["component"]
            targets[v] = (comp, comp if kind is True else kind)
    merged = {}
    if targets:
        pft = site_pft(info.get("pft_path", ""), bxtr.site)
        ins = insitu_fractions(info.get("insitu_path", ""), bxtr.site, pft,
                               info.get("tree_pfts", []), info.get("woody_pfts"))
        full = 1 / bxtr.vars[next(iter(targets))].get("source2sindbadUnit", 1)
        merged = merge_insitu(comps, ins, targets, placement, full)

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        comp = bxtr.vars[tar_name]["component"]
        if tar_name in merged:
            series = merged[tar_name]
        else:
            series = collapse(comps.get(comp, {}), placement)
        data = xr.Dataset({tar_name: xr.DataArray(
            data=place(date_, series, placement), dims=['time'],
            coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    if len(rows) and hasattr(provider, "site_attrs"):
        src_dataset = src_dataset.assign_attrs(**provider.site_attrs(rows, info))
    return src_dataset


def classes(rows, kind):
    """
    Class of each map year in column kind of the cache rows, without missing years.
    """
    if kind not in rows:
        return {}
    rows = rows.drop_duplicates("year", keep="last")
    cls = dict(zip(rows["year"].astype(int), rows[kind].astype(float)))
    return {y: c for y, c in cls.items() if np.isfinite(c)}


def class_badm_series(rows, to_badm):
    """
    BADM PFT by year of the grid class, or of the majority class in years where the
    grid class has none, for the fill of other PFT datasets.
    """
    grid = badm_by_year(classes(rows, "grid"), to_badm)
    return {**badm_by_year(classes(rows, "majority"), to_badm), **grid}


def extract_classes(provider, dataset, site_info, config):
    """
    Generic extract for class maps whose cache rows hold the majority class within
    the radius and the class of the grid cell holding the site, as grid, per map
    year. The component of a variable is majority or grid for the series over time,
    placed with time_placement, or majority_<year> or grid_<year> for the class of
    one map year without a time dimension.

    PFTs come from to_badm of the provider. Where a class has none, such as forest
    of unknown leaf type, the other of majority and grid is used, then the datasets
    listed in rs_fill of the json in order, and then the csv at pft_path, each
    through resolve_pft. A component ending in _merged, such as grid_merged, puts the
    csv first instead, as IGBP numbers through merge_pft. The output attributes
    PFT_grid_<suffix> and PFT_majority_<suffix> are the most frequent PFT over the
    map years of each chain, PFT is PFT_grid_<suffix>, and PFT_source names the
    dataset or csv that gave it.
    """
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    placement = info.get("time_placement", "annual_extend")
    pft_path = info.get("pft_path", "")

    rows = rows_for_site(provider, info, bxtr.site, bxtr.lat, bxtr.lon)
    cls = {k: classes(rows, k) if len(rows) else {} for k in ("majority", "grid")}
    if not cls["grid"]:
        logger.warning(f"::MISSING:: {dataset} has no data for {bxtr.site}. NaN "
                       "values will be set.")
    years = sorted(set(cls["grid"]) | set(cls["majority"]))
    years = years or [int(y) for y in info.get("years", [])]
    own = {k: badm_by_year(cls[k], provider.to_badm) for k in cls}
    # Other datasets are read only when the product has no PFT in any year.
    fills = []
    if not (own["grid"] or own["majority"]):
        fills = fill_sources(config, info.get("rs_fill", []), bxtr.site, bxtr.lat,
                             bxtr.lon)
    chain = {"grid": [(f"{dataset} grid", own["grid"]),
                      (f"{dataset} majority", own["majority"])] + fills,
             "majority": [(f"{dataset} majority", own["majority"]),
                          (f"{dataset} grid", own["grid"])] + fills}

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        kind, _, year = bxtr.vars[tar_name]["component"].partition("_")
        if year == "merged":
            series = merge_pft(date_, pft_path, bxtr.site, years, chain[kind],
                               placement)
            data = xr.DataArray(data=series, dims=['time'], coords={'time': date_})
        elif not year:
            series = place(date_, cls[kind], placement)
            data = xr.DataArray(data=series, dims=['time'], coords={'time': date_})
        else:
            val = cls[kind].get(int(year), np.nan)
            data = xr.DataArray(np.full((1, 1), val), dims=['latitude', 'longitude'],
                                coords={'latitude': [bxtr.lat],
                                        'longitude': [bxtr.lon]})
        data = data.to_dataset(name=tar_name)
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)

    sfx = info.get("var_suffix", dataset)
    grid_pft, grid_src = pft_mode(years, chain["grid"], pft_path, bxtr.site)
    maj_pft, _ = pft_mode(years, chain["majority"], pft_path, bxtr.site)
    return src_dataset.assign_attrs(**{f"PFT_grid_{sfx}": grid_pft,
                                       f"PFT_majority_{sfx}": maj_pft,
                                       "PFT": grid_pft, "PFT_source": grid_src})


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
