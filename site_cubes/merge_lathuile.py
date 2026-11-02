#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge the two sets of LaThuile2007 site files into one read-only copy in
ancillary_data/legacy, so that the meteorology and the fluxes of a site are in one
file.

usage:
    python merge_lathuile.py <legacy folder> [--ltl07 folder] [--v007 folder]
                             [--workers n] [--sites A,B]

The sources are fluxnetBGI2021.LTL07.{DD,HR} of --ltl07, the files of the CLIFF
input with the meteorology but no NEE, GPP or RECO, and fluxnetBGI.LaThuile.{DD,HR}
of --v007, the fluxnetBGI v0.07 extraction of the FLUXCOM v0.04 site cube with the
fluxes, whose daily files have no meteorology. Each merged file is the v0.07
file of the site with every variable that only the LTL07 file has, set on the time
axis of the v0.07 file. A variable of both keeps the v0.07 values and takes the LTL07
value where the v0.07 one is missing. The coordinates are those of the LTL07 file,
which the site registry was built from.

The merged files are written to LaThuile2007/LaThuile2007.{DD,HR}/data/{daily,hourly}
below the legacy folder, as <site>.<start>.<end>.<daily|hourly>.nc, with
merge_report.csv, one row per site, resolution and variable of both files with the
number of time steps with values in both, the largest difference there, and the
number of time steps filled from the LTL07 file. coords_report.csv lists the sites
whose coordinates differ between the copies. The script stops for a site whose LTL07
times are not all within the v0.07 times. It then writes SOURCE.md and MD5SUMS and
makes the folder read only, so it is run once, as the legacy copies are made.
"""
import argparse
import glob
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import xarray as xr

LTL07 = "/Net/Groups/BGI/work_3/cliff/input/fluxnetBGI2021_siteCube-v0.03_iter-5"
V007 = "/Net/Groups/BGI/scratch/skoirala/fluxnetBGI_v0.07"
OLD, NEW, OUT = "fluxnetBGI2021.LTL07", "fluxnetBGI.LaThuile", "LaThuile2007"
RES = {"DD": "daily", "HR": "hourly"}
COORDS = ("latitude", "longitude")


def folder(root, prefix, res):
    return os.path.join(root, f"{prefix}.{res}", "data", RES[res])


def site_files(path):
    return {os.path.basename(f).split(".")[0]: f
            for f in sorted(glob.glob(os.path.join(path, "*.nc")))}


def coords_of(ds):
    return tuple(float(np.ravel(ds[c].values)[0]) for c in COORDS)


def on_new_grid(da, new):
    """A variable of the LTL07 file on the time axis and coordinates of the v0.07
    file, its dimensions in the order of the LTL07 file."""
    order = da.dims
    da = da.squeeze([c for c in COORDS if c in da.dims], drop=True)
    if "time" in da.dims:
        da = da.reindex(time=new["time"])
    extra = {c: new[c].values for c in COORDS if c in order}
    if extra:
        da = da.expand_dims(extra)
    return da.transpose(*order)


def merge_site(job):
    """Write the merged file of one site and resolution, and return its report rows."""
    site, res, old_path, new_path, out_dir = job
    old = xr.open_dataset(old_path).load()
    new = xr.open_dataset(new_path).load()
    outside = np.setdiff1d(old["time"].values, new["time"].values)
    if len(outside):
        raise ValueError(f"{site} {res}: {len(outside)} LTL07 times are not in "
                         f"{new_path}, from {outside[0]}")
    out = new.copy()
    rows = []
    for v in old.data_vars:
        o = on_new_grid(old[v], new)
        if v not in new:
            out[v] = o
            out[v].attrs = old[v].attrs
            continue
        n_vals, o_vals = new[v].values, o.values
        both = np.isfinite(n_vals) & np.isfinite(o_vals) if n_vals.dtype.kind == "f" \
            else np.zeros(n_vals.shape, bool)
        fill = ~np.isfinite(n_vals) & np.isfinite(o_vals) if n_vals.dtype.kind == "f" \
            else np.zeros(n_vals.shape, bool)
        if fill.any():
            out[v] = new[v].where(~fill, o)
            out[v].attrs = new[v].attrs
        rows.append({"site": site, "res": res, "variable": v, "n_both": int(both.sum()),
                     "max_abs_diff": float(np.abs(n_vals[both] - o_vals[both]).max())
                     if both.any() else np.nan, "n_filled_from_ltl07": int(fill.sum())})
    lat, lon = coords_of(old)
    moved = coords_of(new) != (lat, lon)
    out = out.assign_coords(latitude=[lat], longitude=[lon])
    out.attrs.update({
        "latitude": str(lat), "longitude": str(lon),
        "merged_from": f"{new_path} with the variables of {old_path} it lacks",
        "merged_by": "SindbadDataExtractor site_cubes/merge_lathuile.py"})
    years = pd.DatetimeIndex(out["time"].values).year
    path = os.path.join(out_dir, f"{site}.{years.min()}.{years.max()}.{RES[res]}.nc")
    enc = {v: {"zlib": True, "complevel": 4} for v in out.data_vars}
    out.to_netcdf(path + ".tmp", encoding=enc)
    os.replace(path + ".tmp", path)
    coord_row = {"site": site, "res": res, "ltl07_latitude": lat, "ltl07_longitude": lon,
                 "v007_latitude": coords_of(new)[0], "v007_longitude": coords_of(new)[1]}
    return rows, coord_row if moved else None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("legacy")
    p.add_argument("--ltl07", default=LTL07, help="folder of the LTL07 files")
    p.add_argument("--v007", default=V007, help="folder of the v0.07 LaThuile files")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--sites", default="", help="comma separated sites, for a test")
    a = p.parse_args()
    out_root = os.path.join(a.legacy, OUT)
    if os.path.exists(os.path.join(out_root, "SOURCE.md")):
        raise SystemExit(f"{out_root} holds a merged copy already, see its SOURCE.md")
    jobs = []
    for res in RES:
        old = site_files(folder(a.ltl07, OLD, res))
        new = site_files(folder(a.v007, NEW, res))
        missing = sorted(set(old) ^ set(new))
        if missing:
            print(f"{res}: sites in only one copy, left out: {missing}")
        out_dir = folder(out_root, OUT, res)
        os.makedirs(out_dir, exist_ok=True)
        sites = sorted(set(old) & set(new))
        if a.sites:
            sites = [s for s in sites if s in a.sites.split(",")]
        jobs += [(s, res, old[s], new[s], out_dir) for s in sites]
    print(f"{len(jobs)} files to merge with {a.workers} workers", flush=True)
    rows, coords = [], []
    with ProcessPoolExecutor(a.workers) as ex:
        for job, (r, c) in zip(jobs, ex.map(merge_site, jobs)):
            rows += r
            if c:
                coords.append(c)
            print(f"done {job[0]} {job[1]}", flush=True)
    report = pd.DataFrame(rows)
    report.to_csv(os.path.join(out_root, "merge_report.csv"), index=False)
    pd.DataFrame(coords).to_csv(os.path.join(out_root, "coords_report.csv"), index=False)
    differ = report[report["max_abs_diff"] > 0]
    print(f"{len(differ)} site variables differ where both copies have values, "
          f"{int(report['n_filled_from_ltl07'].sum())} values filled from LTL07, "
          f"{len(coords)} files with other coordinates in the v0.07 copy")
    if a.sites:
        print("a test of some sites: no SOURCE.md, MD5SUMS or write protection")
        return
    md5 = subprocess.run("find . -type f -name '*.nc' | sort | xargs md5sum", shell=True,
                         cwd=out_root, capture_output=True, text=True, check=True).stdout
    open(os.path.join(out_root, "MD5SUMS"), "w").write(md5)
    n = {r: len(glob.glob(os.path.join(folder(out_root, OUT, r), "*.nc"))) for r in RES}
    stamp = pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ")
    open(os.path.join(out_root, "SOURCE.md"), "w").write(
        "# LaThuile2007 site files\n\n"
        "The LaThuile2007 source of the site cubes, with the meteorology and the fluxes "
        "of a site in one\nfile, merged by site_cubes/merge_lathuile.py of "
        "SindbadDataExtractor. Never edit these files.\n\n"
        f"- Fluxes and time axis: {a.v007}/{NEW}.{{DD,HR}}\n"
        f"- Added variables, gaps and coordinates: {a.ltl07}/{OLD}.{{DD,HR}}\n"
        f"- Merged: {stamp} by {os.environ.get('USER', '')}\n"
        f"- Files: {n['DD']} daily, {n['HR']} hourly\n"
        "- Reports: merge_report.csv, coords_report.csv\n"
        "- Checksums: MD5SUMS\n")
    subprocess.run(["chmod", "-R", "a-w", out_root], check=True)
    print(f"merged copy written to {out_root}")


if __name__ == '__main__':
    main()
