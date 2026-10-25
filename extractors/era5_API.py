#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hourly ERA5 and ERA5-Land series at the sites, fetched once by prefetch.py into one
NetCDF file per site and source, for the climate forcing of SINDBAD.

The json of a product lists its sources, each with the variables it gives. A source
is one of
    earthmover: the Icechunk repo of Earthmover, read with arraylake after
        arraylake auth login. One site is read at a time from the group laid out for
        time series.
    cds: a time-series dataset of the Copernicus Climate Data Store, named by id,
        such as reanalysis-era5-land-timeseries, with the key in ~/.cdsapirc. One request per
        site, and the CDS runs max_jobs requests of a user at a time.
    edh: a zarr store of Earth Data Hub, with the key in ~/.netrc. The sites are read
        together for each spatial chunk of the store, so that every chunk is read
        once, since Earth Data Hub allows a fixed number of requests per month. With
        accumulated true, the values accumulate from 00 UTC over the day as in the
        original ERA5-Land and are differenced to hourly sums.

Each source takes the grid cell nearest the site. For a store that is NaN over the
sea, such as ERA5-Land, the nearest cell with data is taken, and its distance from
the site is kept as grid_km. The full record of each source is read, from its first
to its last time step, in UTC.

The files are written to the folder of cache_path as <site>_<source>.nc, with the
variables of the source in their original names and units and the global attributes
site, source, dataset, grid_latitude, grid_longitude, fetched_utc and citation.
cache_path is an index csv with one row per site and source, which prefetch.py reads
to find the sites still missing. A site counts as cached when every source of the
product has a row.
"""
import os
import shutil
import tempfile
import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

# Encoding of the cached variables.
ENCODING = {"dtype": "float32", "zlib": True, "complevel": 4}
EDH_OPTIONS = {"client_kwargs": {"trust_env": True}}
# Polling interval of the CDS jobs in seconds.
CDS_POLL_S = 10
# Attempts and first wait in seconds for a remote read that fails, doubled each time.
TRIES = 5
WAIT_S = 30


def retry(fn, what):
    """fn() retried TRIES times with waits of WAIT_S seconds, doubled each time."""
    for attempt in range(TRIES):
        try:
            return fn()
        except Exception as e:
            if attempt == TRIES - 1:
                raise
            wait = WAIT_S * 2 ** attempt
            print(f"{what}: attempt {attempt + 1} failed, retrying in {wait} s: {e}",
                  flush=True)
            time.sleep(wait)


def cache_folder(info):
    return os.path.dirname(info["cache_path"])


def site_file(info, site, source):
    return os.path.join(cache_folder(info), f"{site}_{source}.nc")


def complete_sites(cache, info):
    """Sites with a row for every source of the product in the index."""
    if len(cache) == 0:
        return set()
    need = set(info["sources"])
    got = cache.groupby(cache["site"].astype(str).str.lower())["source"].agg(set)
    return {s for s, srcs in got.items() if need <= srcs}


def missing(info, source, sites, lats, lons):
    """Sites without a row for source in the index."""
    cache = fu.read_cache(info["cache_path"])
    done = set()
    if len(cache):
        done = set(cache.loc[cache["source"] == source, "site"].astype(str).str.lower())
    keep = [i for i, s in enumerate(sites) if s.lower() not in done]
    return [sites[i] for i in keep], [lats[i] for i in keep], [lons[i] for i in keep]


def write_site(info, source, site, lat, lon, ds, grid_lat, grid_lon, grid_km):
    """Write the hourly series of one site and source and add its index row."""
    cfg = info["sources"][source]
    ds = ds.drop_vars([c for c in ("latitude", "longitude") if c in ds.coords])
    if "valid_time" in ds.dims:
        ds = ds.rename({"valid_time": "time"})
    fetched = fu.now_utc()
    ds.attrs = {"site": site, "site_latitude": lat, "site_longitude": lon,
                "source": source, "dataset": cfg["dataset"],
                "grid_latitude": grid_lat, "grid_longitude": grid_lon,
                "grid_km": round(grid_km, 3), "fetched_utc": fetched,
                "citation": info.get("citation", ""),
                "time_zone": "UTC"}
    path = site_file(info, site, source)
    tmp = path + f".tmp{os.getpid()}{threading.get_ident()}"
    ds.to_netcdf(tmp, encoding={v: ENCODING for v in ds.data_vars})
    os.replace(tmp, path)
    nan = {v: int(ds[v].isnull().sum()) for v in ds.data_vars}
    fu.append_cache(info["cache_path"], [{
        "site": site, "source": source, "file": os.path.basename(path),
        "latitude": lat, "longitude": lon, "grid_latitude": grid_lat,
        "grid_longitude": grid_lon, "grid_km": round(grid_km, 3),
        "start": str(ds.time.values[0])[:13], "end": str(ds.time.values[-1])[:13],
        "n_hours": ds.sizes["time"], "n_nan": ";".join(f"{k}:{v}" for k, v in nan.items()),
        "fetched_utc": fetched}])


def nearest(lat_grid, lon_grid, lat, lon, valid=None):
    """Indices of the grid cell nearest the site, among the valid cells when given."""
    lon360 = lon % 360
    if valid is None:
        dlon = np.abs((lon_grid - lon360 + 180) % 360 - 180)
        return int(np.abs(lat_grid - lat).argmin()), int(dlon.argmin())
    ii, jj = np.nonzero(valid)
    d = fu.haversine(lat, lon360, lat_grid[ii], lon_grid[jj])
    k = int(d.argmin())
    return int(ii[k]), int(jj[k])


# Earthmover

def fetch_earthmover(info, sites, lats, lons, workers):
    from arraylake import Client
    cfg = info["sources"]["earthmover"]
    repo = retry(lambda: Client().get_repo(cfg["repo"]), "earthmover open")
    ds = retry(lambda: xr.open_zarr(repo.readonly_session(cfg.get("branch", "main")).store,
                                    group=cfg["group"], consolidated=False, chunks=None),
               "earthmover open")
    lat_grid, lon_grid = ds.latitude.values, ds.longitude.values

    def one(args):
        site, lat, lon = args
        try:
            i, j = nearest(lat_grid, lon_grid, lat, lon)
            p = retry(lambda: ds[cfg["variables"]].isel(latitude=i, longitude=j).load(),
                      f"earthmover {site}")
            write_site(info, "earthmover", site, lat, lon, p, float(lat_grid[i]),
                       float(lon_grid[j]),
                       float(fu.haversine(lat, lon % 360, lat_grid[i], lon_grid[j]))
                       / 1000)
            print(f"earthmover done {site}", flush=True)
        except Exception as e:
            print(f"earthmover failed {site}: {e}", flush=True)

    todo = list(zip(*missing(info, "earthmover", sites, lats, lons)))
    print(f"earthmover: {len(todo)} sites to fetch", flush=True)
    with ThreadPoolExecutor(workers or cfg.get("n_workers", 4)) as ex:
        list(ex.map(one, todo))


# Copernicus Climate Data Store

def cds_client():
    from ecmwf.datastores import Client
    cfg = dict(line.split(":", 1) for line in open(os.path.expanduser("~/.cdsapirc"))
               if ":" in line)
    return Client(url=cfg["url"].strip(), key=cfg["key"].strip())


def cds_end_date(client, dataset):
    """Last day of the dataset from its catalogue entry."""
    col = client.get_collection(dataset)
    end = col.end_datetime
    return str(end)[:10] if end else None


def read_cds_result(zip_path):
    """Merge the NetCDF files of a CDS time-series zip into one dataset."""
    folder = tempfile.mkdtemp(dir=os.path.dirname(zip_path))
    try:
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(folder)
        files = sorted(os.path.join(folder, f) for f in os.listdir(folder)
                       if f.endswith(".nc"))
        ds = xr.merge([xr.open_dataset(f).load() for f in files], compat="override")
    finally:
        shutil.rmtree(folder, ignore_errors=True)
    return ds


def fetch_cds(info, sites, lats, lons, workers):
    cfg = info["sources"]["cds"]
    client = cds_client()
    end = cfg.get("end_date") or cds_end_date(client, cfg["id"])
    date = f"{cfg['start_date']}/{end}"
    todo = list(zip(*missing(info, "cds", sites, lats, lons)))
    print(f"cds: {len(todo)} sites to fetch, {date}, {cfg.get('max_jobs', 3)} at a "
          "time", flush=True)
    running = {}
    while todo or running:
        while todo and len(running) < cfg.get("max_jobs", 3):
            site, lat, lon = todo.pop(0)
            req = {"variable": cfg["request_variables"],
                   "location": {"latitude": float(lat),
                                "longitude": float((lon + 180) % 360 - 180)},
                   "date": [date], "data_format": "netcdf"}
            try:
                running[site] = (client.submit(cfg["id"], req), lat, lon)
            except Exception as e:
                print(f"cds failed to submit {site}: {e}", flush=True)
        time.sleep(CDS_POLL_S)
        for site in list(running):
            job, lat, lon = running[site]
            try:
                status = job.status
            except Exception as e:
                print(f"cds status error {site}: {e}", flush=True)
                continue
            if status == "successful":
                del running[site]
                try:
                    zip_path = os.path.join(cache_folder(info), f"{site}_cds.zip")
                    job.download(zip_path)
                    ds = read_cds_result(zip_path)
                    os.remove(zip_path)
                    glat = float(ds["latitude"].values)
                    glon = float(ds["longitude"].values)
                    ds = ds[[v for v in cfg["variables"] if v in ds]]
                    write_site(info, "cds", site, lat, lon, ds, glat, glon,
                               float(fu.haversine(lat, lon, glat, glon)) / 1000)
                    print(f"cds done {site}", flush=True)
                except Exception as e:
                    print(f"cds failed {site}: {e}", flush=True)
            elif status in ("failed", "rejected", "dismissed"):
                del running[site]
                print(f"cds {status} {site}", flush=True)


# Earth Data Hub

def deaccumulate(values, hours):
    """
    Hourly sums from values accumulated from 00 UTC, where the value at 00 UTC is
    the total of the previous day. The first time step stays NaN.
    """
    out = np.full_like(values, np.nan)
    out[1:] = values[1:] - values[:-1]
    restart = hours == 1
    out[restart] = values[restart]
    return out


def fetch_edh(info, sites, lats, lons, workers):
    cfg = info["sources"]["edh"]
    ds = retry(lambda: xr.open_dataset(cfg["url"], engine="zarr", chunks=None,
                                       storage_options=EDH_OPTIONS), "edh open")
    lat_grid, lon_grid = ds.latitude.values, ds.longitude.values
    clat, clon = ds[cfg["variables"][0]].encoding["chunks"][1:]
    todo = list(zip(*missing(info, "edh", sites, lats, lons)))
    # Sites grouped by the spatial chunk of the nearest grid cell.
    tiles = {}
    for site, lat, lon in todo:
        i, j = nearest(lat_grid, lon_grid, lat, lon)
        tiles.setdefault((i // clat, j // clon), []).append((site, lat, lon))
    print(f"edh: {len(todo)} sites in {len(tiles)} chunks, {cfg['variables']}",
          flush=True)
    hours = ds.valid_time.dt.hour.values

    def one(key):
        members = tiles[key]
        try:
            cells = []
            for site, lat, lon in members:
                if cfg.get("land_only", False):
                    # Nearest valid cell, searched in the chunk and its neighbours.
                    i, j = nearest_valid(ds, cfg, lat_grid, lon_grid, lat, lon)
                else:
                    i, j = nearest(lat_grid, lon_grid, lat, lon)
                cells.append((i, j))
            ii = xr.DataArray([c[0] for c in cells], dims="site")
            jj = xr.DataArray([c[1] for c in cells], dims="site")
            data = retry(lambda: ds[cfg["variables"]].isel(latitude=ii,
                                                           longitude=jj).load(),
                         f"edh chunk {key}")
            for k, (site, lat, lon) in enumerate(members):
                p = data.isel(site=k)
                if cfg.get("accumulated", False):
                    for v in cfg["variables"]:
                        p[v] = (p[v].dims, deaccumulate(p[v].values, hours),
                                p[v].attrs)
                i, j = cells[k]
                write_site(info, "edh", site, lat, lon, p.drop_vars("site", errors="ignore"),
                           float(lat_grid[i]), float(lon_grid[j]),
                           float(fu.haversine(lat, lon % 360, lat_grid[i], lon_grid[j]))
                           / 1000)
            print(f"edh done chunk {key}: {', '.join(m[0] for m in members)}", flush=True)
        except Exception as e:
            print(f"edh failed chunk {key}: {e}", flush=True)

    with ThreadPoolExecutor(workers or cfg.get("n_workers", 2)) as ex:
        list(ex.map(one, list(tiles)))


def nearest_valid(ds, cfg, lat_grid, lon_grid, lat, lon, half=0.5):
    """
    Nearest grid cell with data within half degrees of the site, from the last time
    step of the first variable, for stores that are NaN over the sea.
    """
    lat_sel = np.nonzero(np.abs(lat_grid - lat) <= half)[0]
    dlon = np.abs((lon_grid - lon % 360 + 180) % 360 - 180)
    lon_sel = np.nonzero(dlon <= half / max(np.cos(np.radians(lat)), 0.1))[0]
    if len(lat_sel) == 0 or len(lon_sel) == 0:
        return nearest(lat_grid, lon_grid, lat, lon)
    vals = ds[cfg["variables"][0]].isel(valid_time=-1, latitude=lat_sel,
                                        longitude=lon_sel).load().values
    valid = np.isfinite(vals)
    if not valid.any():
        return nearest(lat_grid, lon_grid, lat, lon)
    i, j = nearest(lat_grid[lat_sel], lon_grid[lon_sel], lat, lon, valid)
    return int(lat_sel[i]), int(lon_sel[j])


FETCHERS = {"earthmover": fetch_earthmover, "cds": fetch_cds, "edh": fetch_edh}


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    """Fetch every source of the product at the same time, each in its own thread."""
    os.makedirs(cache_folder(info), exist_ok=True)
    # The threads would import zarr at the same time, which can fail with a
    # circular import, so it is imported here first.
    import zarr  # noqa: F401
    if "earthmover" in info["sources"]:
        import arraylake  # noqa: F401
    threads = [threading.Thread(target=FETCHERS[src], args=(info, sites, lats, lons,
                                                            workers), name=src)
               for src in info["sources"]]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    done = complete_sites(fu.read_cache(info["cache_path"]), info)
    print(f"{len(done)} of {len(sites)} sites have every source in "
          f"{info['cache_path']}", flush=True)


def extract(dataset, site_info, config):
    logger.warning(f"{dataset}: the extractor of the ERA5 forcing is not written yet, "
                   "only prefetch.py uses this module.")
    return None


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
