#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Forcing of the CLIFF gap filling: a gridded meteorology at the nearest grid cells of
every site, filled once by prefetch.py as the pixel files that pycliff-downscale
reads, one per variable and site:
    <cache_path>/<var>/<site>/<var>_Z_BIG_<start>_<end>_pix_ALL_of_<n>.nc
with var(time, nclosest), lat(nclosest) and lon(nclosest), time in hours since
1900-01-01 UTC and the values in the units of ERA5. The variables are the ERA5 names
that PyCLIFF uses, t2m, q, tp, ssrd, ssr, strd, ws10, rH and sp, and q is taken at
the pressure level of the json. A forcing other than ERA5 is written under the same
names and units, so PyCLIFF needs no change to use it. db of the json names the
forcing in the PyCLIFF products, as ERA5.e1.

source of the json says where the values come from:
    pycliff_extract: pycliff-extract, run with the python of the PyCLIFF environment
        and the checkout pycliff_dir, reads the monthly files below pbase once per
        variable for all missing sites, so the forcing is ERA5 of the cluster by
        default. Another gridded product with the same file layout only needs
        another pbase.
    <module>:<function>: any other forcing, read by a function of this repository
        called as function(var, site, lat, lon, info). It returns the UTC times as
        datetime64, the values as an array of grid cells by times, and the latitude
        and longitude of each grid cell, or None when it has no data for the site.
A file is missing when the site has none for the variable, or when its last year is
before end_year of the json, in which case the earlier file is replaced. Variables
are independent, so prefetch.py --part i/n splits them over jobs of a cluster array.
"""
import glob
import importlib
import os
import re
import subprocess
import time

import numpy as np

import logging
logger = logging.getLogger(__name__)

PLEV = {"q": 1000}
VARIABLES = ["t2m", "q", "tp", "ssrd", "ssr", "strd", "ws10", "rH", "sp"]
FILE = "{var}/{site}/{var}_Z_BIG_{start}_{end}_pix_ALL_of_{n}.nc"
END_YEAR = re.compile(r"_Z_BIG_\d+_(\d+)_pix_ALL_of_")


def variables(info):
    return info.get("variables", VARIABLES)


def pixel_path(info, var, site):
    return os.path.join(info["cache_path"], FILE.format(
        var=var, site=site, start=info["start_year"], end=info["end_year"],
        n=info.get("n_pixels", 3)))


def existing(info, var, site):
    return sorted(glob.glob(os.path.join(info["cache_path"], var, site,
                                         f"{var}_Z_BIG_*_pix_ALL_of_*.nc")))


def is_missing(info, var, site):
    """No file of the variable at the site, or one that ends before end_year."""
    files = existing(info, var, site)
    if not files:
        return True
    m = END_YEAR.search(os.path.basename(files[-1]))
    return not m or int(m.group(1)) < int(info["end_year"])


def drop_older(info, var, site):
    """Remove the files of a site that the newest extraction replaced."""
    keep = pixel_path(info, var, site)
    if os.path.exists(keep):
        for f in existing(info, var, site):
            if f != keep:
                os.remove(f)


def write_pixels(path, var, times, values, lats, lons):
    """A pixel file in the layout of pycliff-extract."""
    import netCDF4
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("nclosest", values.shape[0])
        ds.createDimension("time", values.shape[1])
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = "hours since 1900-01-01 00:00:00"
        tv.calendar = "proleptic_gregorian"
        base = np.datetime64("1900-01-01T00:00:00", "s")
        tv[:] = (times.astype("datetime64[s]") - base) / np.timedelta64(1, "h")
        ds.createVariable("nclosest", "i8", ("nclosest",))[:] = \
            np.arange(1, values.shape[0] + 1)
        ds.createVariable("lon", "f8", ("nclosest",))[:] = lons
        ds.createVariable("lat", "f8", ("nclosest",))[:] = lats
        dv = ds.createVariable(var, "f8", ("time", "nclosest"), zlib=True, complevel=1)
        dv[:] = np.asarray(values, dtype=float).T
        ds.source = "SindbadDataExtractor cliff_forcing"
    os.replace(tmp, path)


def pycliff_env(info):
    """
    Environment of a PyCLIFF call: its output root, and the checkout pycliff_dir of
    the json first on PYTHONPATH, so that this checkout runs whatever the python of
    the PyCLIFF environment has installed.
    """
    env = dict(os.environ, PYCLIFF_OUT=info.get("pycliff_out") or info["cache_path"])
    if info.get("pycliff_dir"):
        env["PYTHONPATH"] = os.pathsep.join(
            [info["pycliff_dir"]] + [p for p in [os.environ.get("PYTHONPATH")] if p])
    return env


def extract_with_pycliff(info, var, todo, workers):
    """Run pycliff-extract for one variable and the missing sites."""
    req_dir = os.path.join(info["cache_path"], "_requests")
    os.makedirs(req_dir, exist_ok=True)
    batch = os.path.join(req_dir, f"{var}.txt")
    plev = PLEV.get(var)
    with open(batch, "w") as f:
        for site, lat, lon in todo:
            line = (f" -x {lon} -y {lat} -n {info.get('n_pixels', 3)} "
                    f"--year1={info['start_year']} --year2={info['end_year']}")
            if plev is not None:
                line += f" --plev={plev}"
            f.write(f"{line} {var} {pixel_path(info, var, site)}\n")
    cmd = [info["pycliff_python"], "-m", "pycliff.extract.cli", "-b", batch,
           "--pbase", info["pbase"], "--workers", str(workers), var]
    # Sites with an earlier extraction only read the months after it.
    if any(existing(info, var, s) for s, _, _ in todo):
        cmd += ["--append-from", info["cache_path"]]
    env = pycliff_env(info)
    t = time.time()
    run = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(f"pycliff-extract {var} failed:\n{run.stderr[-3000:]}")
    for site, _, _ in todo:
        drop_older(info, var, site)
    print(f"{var}: {len(todo)} sites extracted in {time.time() - t:.0f} s", flush=True)


def extract_with_function(info, var, todo):
    """Write the pixel files of one variable from a function of this repository."""
    module, func = info["source"].split(":")
    fetch = getattr(importlib.import_module(module), func)
    for site, lat, lon in todo:
        try:
            got = fetch(var, site, lat, lon, info)
        except Exception as e:
            print(f"{var}: failed {site}: {e}", flush=True)
            continue
        if got is None:
            print(f"{var}: no data for {site}", flush=True)
            continue
        times, values, plats, plons = got
        write_pixels(pixel_path(info, var, site), var, np.asarray(times),
                     np.atleast_2d(values), plats, plons)
        drop_older(info, var, site)
    print(f"{var}: {len(todo)} sites written", flush=True)


def prefetch(sites, lats, lons, info, workers=None, part="", **kwargs):
    """Extract the missing sites of every variable of this part."""
    from extractors.clms_dekad import select_part
    n_workers = workers or info.get("n_workers", 8)
    for var in select_part(variables(info), part):
        todo = [s for s in zip(sites, lats, lons) if is_missing(info, var, s[0])]
        print(f"{var}: {len(todo)} of {len(sites)} sites to extract into "
              f"{info['cache_path']}", flush=True)
        if not todo:
            continue
        if info.get("source", "pycliff_extract") == "pycliff_extract":
            extract_with_pycliff(info, var, todo, n_workers)
        else:
            extract_with_function(info, var, todo)


def report(sites, lats, lons, info):
    """Print the sites that miss a pixel file of any variable."""
    missing = {s for s in sites for v in variables(info) if is_missing(info, v, s)}
    print(f"{info.get('db', '')}: {len(sites) - len(missing)} of {len(sites)} sites "
          f"have every variable up to {info['end_year']} in {info['cache_path']}")
    if missing:
        print("    missing: " + " ".join(sorted(missing)))


def complete(info, site):
    """The site has the pixel files of every variable."""
    return all(existing(info, v, site) for v in variables(info))


if __name__ == '__main__':
    print(__doc__)
