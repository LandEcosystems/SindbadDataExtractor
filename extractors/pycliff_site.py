#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tower meteorology of a site gap filled by CLIFF, run for each site of an eddy
covariance source while the source is extracted, with PyCLIFF left unchanged.

For a site with tower meteorology, extract runs pycliff-downscale for the site
with the python of the PyCLIFF environment and the PyCLIFF checkout pycliff_dir first
on PYTHONPATH, and then reads its gap-filled product with the reader of
pycliff_gapfilled.py:
    <VAR>_F       the tower value where its QC is at least qc_threshold of the
                  json, 0.85 by default, else the meteorology of the forcing
                  downscaled by the selected method
    <VAR>_F_QC    hourly: 0 tower, 1 selected method, 2 nearest grid cell;
                  daily: the fraction of the day from the tower
The tower data are the hourly CLIFF input file of the site, found with site_input
of the json, a pattern with {site}. extract_product.py sets site_input and fn for
each source from the pycliff block of the catalog, so that every source has its
own product per forcing, pycliff_era5.<source> or pycliff_crujra.<source>. A site
whose file is missing, or whose TA or SW_IN has fewer than min_tower_days days with
data, gets no product.

The forcing is a cliff_forcing dataset, whose cache_path is pix_dir of the json and
whose db names it, so that CLIFF runs with ERA5 or any other forcing written in the
same layout. The hourly product runs PyCLIFF at hourly and the daily product at
daily resolution, as fn.HR and fn.DD, with the fixed seed, units, time mode and
methods of the json. The products are written below pycliff_out, which PyCLIFF takes
as PYCLIFF_OUT.

PyCLIFF skips a site that is done with the same settings, but its check does not
see changed input files. A marker per site, run and forcing in
<pycliff_out>/sindbad_inputs/<fn>/<db>/<site>.txt therefore keeps the size and time
of the tower file and the pixel files, and PyCLIFF runs with --force when they
changed. A lock per marker keeps parallel extract jobs from running the same fit
twice. forcing_variables of the json lists the pixel files of the forcing, the
ERA5 variables by default.

skip_forcing_variables of the json lists forcing variables that PyCLIFF leaves out,
with --skip-grid-vars, while their pixel files are not ready, such as ws10 of ERA5.
Their products are not written and read as NaN. They are also left out of the
marker, so taking a variable off the list changes the marker and fits every site
again with it, while the other variables come out the same, each with its own seed.

PyCLIFF refuses a resolution finer than the step of the forcing, so a forcing on
6 hour steps such as CRUJRA only has a daily product. PyCLIFF fits the variables of
a site in PYCLIFF_JOBS processes, one by default, as when many sites share the
cpus of a job; run_site.py sets it when a few sites have more cpus.

The log of pycliff-downscale goes line by line into the log of the site, as it
runs, so a long fit can be followed there.
"""
import glob
import json
import os

import numpy as np

import extractors.footprint_utils as fu
import extractors.pycliff_gapfilled as pg
from extractors.cliff_forcing import VARIABLES, existing, pycliff_env, run_logged

import logging
logger = logging.getLogger(__name__)

TOWER_CHECK = ("TA", "SW_IN")


def tower_file(info, site):
    files = sorted(glob.glob(info["site_input"].format(site=site)))
    return files[-1] if files else None


def tower_days(path):
    """Days with any TA and with any SW_IN in the tower file, the smaller of both."""
    import netCDF4
    with netCDF4.Dataset(path) as ds:
        n = []
        for v in TOWER_CHECK:
            if v not in ds.variables:
                return 0
            x = np.ma.filled(np.ma.asarray(ds.variables[v][:], dtype=float), np.nan)
            x = np.ravel(x)
            per_day = 24 if len(x) % 24 == 0 else 1
            n.append(int(np.isfinite(x.reshape(-1, per_day)).any(axis=1).sum()))
    return min(n)


def run_name(info, freq):
    return f"{info['fn']}.{'HR' if freq == 'hourly' else 'DD'}"


def inputs(info, site, tower):
    """Tower file and pixel files of a site, with their size and time."""
    skip = set(info.get("skip_forcing_variables", []))
    paths = [tower] + [f for v in info.get("forcing_variables", VARIABLES)
                       if v not in skip
                       for f in existing({"cache_path": info["pix_dir"]}, v, site)]
    return [(p, os.path.getsize(p), os.path.getmtime(p)) for p in paths]


def source_files(info, site, freq):
    """Inputs of the fit, whose changes extract_product.py includes in its hash."""
    tower = tower_file(info, site)
    return [p for p, _, _ in inputs(info, site, tower)] if tower else []


def run_pycliff(info, site, tower, freq):
    """Fit and write the PyCLIFF product of a site unless it is current."""
    fn = run_name(info, freq)
    marker = os.path.join(info["pycliff_out"], "sindbad_inputs", fn, info["db"],
                          f"{site}.txt")
    os.makedirs(os.path.dirname(marker), exist_ok=True)
    stamp = json.dumps(inputs(info, site, tower))
    # markers written before the forcing was part of the path; its pixel files
    # are in the stamp, so another forcing never matches it
    legacy = os.path.join(info["pycliff_out"], "sindbad_inputs", fn, f"{site}.txt")
    with fu.file_lock(marker):
        seen = marker if os.path.exists(marker) else legacy
        changed = not os.path.exists(seen) or open(seen).read() != stamp
        cmd = [info["pycliff_python"], "-m", "pycliff.downscale.cli", site,
               "--site-file", tower, "--pix-dir", info["pix_dir"], "--fn-version", fn,
               "--db-id", info["db"], "--resolution", freq,
               "--seed", str(info.get("seed", 42)), "--units", info.get("units", "si"),
               "--time-mode", info.get("time_mode", "lst"),
               "--sites-csv", info["sites_csv"],
               "--jobs", os.environ.get("PYCLIFF_JOBS", "1"),
               "--no-figures", "--no-analysis", "--force" if changed else "--skip-done"]
        if info.get("methods"):
            cmd += ["--methods", info["methods"]]
        if info.get("skip_forcing_variables"):
            cmd += ["--skip-grid-vars", ",".join(info["skip_forcing_variables"])]
        if info.get("qc_threshold") is not None:
            cmd += ["--qc-threshold", str(info["qc_threshold"])]
        env = dict(pycliff_env(info), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1")
        logger.info("pycliff-downscale: " + " ".join(cmd))
        code, last = run_logged(cmd, env,
                                lambda line: logger.info(f"pycliff {fn}: {line}"))
        if code:
            raise RuntimeError(f"pycliff-downscale {site} {fn} failed:\n{last}")
        with open(marker, "w") as f:
            f.write(stamp)


def extract(dataset, site_info, config):
    info = config["dataset"][dataset]
    freq = config["temporal_resolution"]
    site = site_info["site_ID"]
    tower = tower_file(info, site)
    if tower is None:
        logger.warning(f"::MISSING:: {dataset} has no tower input of {info['fn']} for "
                       f"{site} at {info['site_input']}. No product is made.")
        return None
    days = tower_days(tower)
    if days < info.get("min_tower_days", 365):
        logger.warning(f"::MISSING:: {site} has {days} days of tower TA and SW_IN in "
                       f"{tower}, fewer than {info.get('min_tower_days', 365)}. No "
                       f"product is made.")
        return None
    run_pycliff(info, site, tower, freq)
    reader = dict(info, source_dir=os.path.join(info["pycliff_out"], "output",
                                                "gapfilled"))
    config = dict(config, dataset={**config["dataset"], dataset: reader})
    ds = pg.extract(dataset, site_info, config)
    for v in ds.data_vars:
        ds[v].attrs.update({"pycliff_tower_file": tower, "pycliff_db": info["db"],
                            "pycliff_forcing": info["pix_dir"],
                            "pycliff_qc_threshold": info.get("qc_threshold", 0.85)})
    return ds


if __name__ == '__main__':
    print(__doc__)
