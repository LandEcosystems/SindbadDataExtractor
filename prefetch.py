#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fill the site caches of the footprint extractors once, before run_extraction.py.

Remote products are fetched for all sites missing from the cache_path of each
selected dataset, in batches that suit each source. Sites already in the cache are
skipped, so an interrupted run resumes where it stopped and later runs only fetch new
sites.

usage:
    python prefetch.py <main json> [dataset ...] [--sites-csv csv] [--version v]
                       [--workers n] [--step submit|status|process]

The datasets are keys of the dataset block of the main json, by default all datasets
whose json sets cache_path. Their jsons are read from the folder of the main json,
or from the ancillary folder next to it.
Sites and coordinates come from --sites-csv with site, latitude and longitude columns
when given, otherwise from the FLUXCOM site cube of --version, as in
run_extraction.py. --step is only used by vf_modis, whose AppEEARS request is
asynchronous.

--check fetches nothing. It lists for each dataset the sites missing from its cache
and the cached sites whose coordinates differ from those of the site list by more
than footprint_utils.COORD_TOL_DEG, which the extraction would use with a warning.
Rows cached before the coordinates were stored with them cannot be compared and are
counted as without stored coordinates.

Providers whose cache needs several sources per site, such as era5_API, define
complete_sites(cache, info), which --check uses to count the cached sites.

--cache-dir moves the folder that holds the caches of the datasets to another folder,
for example to prefetch on a local computer and copy the caches to the cluster
afterwards. The folder replaced is the common folder of the cache_path of all
datasets with a cache, and the paths in PATH_KEYS below it keep their place relative
to it, so a local copy needs the same layout.
"""
import argparse
import importlib
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from utils.shared_utils import dataset_json_path  # noqa: E402
import extractors.footprint_utils as fu  # noqa: E402

# Paths of a dataset json that point into the cache folder.
PATH_KEYS = ["cache_path", "peak_dir", "work_dir", "tiles_path", "code_dir",
             "shuttle_cache_path"]


def site_table(args, main):
    if args.sites_csv:
        df = pd.read_csv(args.sites_csv)
    else:
        import fluxcom.providers as prov
        from utils.workflow_extraction import get_site_info
        version = args.version or main["FLUXNET_version"]
        cube = args.cubepath
        sites = sorted(prov.EddyProvider(version=version, cubepath=cube).sites)
        df = pd.DataFrame([get_site_info(s) for s in sites]).rename(
            columns={"site_ID": "site"})
    if main.get("sites_toRun", ["all"]) != ["all"] and not args.all_sites:
        df = df[df["site"].isin(main["sites_toRun"])]
    return (df["site"].astype(str).tolist(), df["latitude"].astype(float).tolist(),
            df["longitude"].astype(float).tolist())


def move_cache_dir(infos, cache_dir):
    """Replace the common folder of the cache paths by cache_dir in all PATH_KEYS."""
    caches = [v["cache_path"] for v in infos.values() if v.get("cache_path")]
    if not caches:
        return
    root = os.path.commonpath([os.path.dirname(c) for c in caches])
    print(f"cache folder {root} -> {cache_dir}")
    for info in infos.values():
        for key in PATH_KEYS:
            path = info.get(key, "")
            if path and os.path.commonpath([root, path]) == root:
                info[key] = os.path.join(cache_dir, os.path.relpath(path, root))


def check(name, info, sites, lats, lons, provider=None):
    """
    Print the sites of the list that the cache of a dataset lacks or misplaces. A
    provider with complete_sites(cache, info) decides which sites count as cached.
    """
    path = info["cache_path"]
    cache = fu.read_cache(path)
    if not os.path.exists(path):
        print(f"{name}: no cache at {path}, all {len(sites)} sites missing")
        return
    if hasattr(provider, "complete_sites"):
        have = provider.complete_sites(cache, info)
    else:
        have = {s.lower() for s in cache["site"].astype(str)}
    missing = [s for s in sites if s.lower() not in have]
    moved, unknown = [], 0
    for site, lat, lon in zip(sites, lats, lons):
        rows = fu.site_rows(cache, site)
        if len(rows) == 0:
            continue
        if not {"latitude", "longitude"} <= set(rows.columns) or \
                rows[["latitude", "longitude"]].isna().all().any():
            unknown += 1
            continue
        d = max(np.nanmax(np.abs(rows["latitude"].astype(float) - lat)),
                np.nanmax(np.abs(rows["longitude"].astype(float) - lon)))
        if d > fu.COORD_TOL_DEG:
            moved.append(f"{site} ({d:.4f} deg)")
    print(f"{name}: {len(sites) - len(missing)} of {len(sites)} sites cached, "
          f"{len(missing)} missing, {len(moved)} at other coordinates, {unknown} "
          f"without stored coordinates, in {path}")
    if missing:
        print("    missing: " + " ".join(missing))
    if moved:
        print("    other coordinates: " + ", ".join(moved))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("main_json")
    p.add_argument("datasets", nargs="*")
    p.add_argument("--sites-csv", default="")
    p.add_argument("--version", default="")
    p.add_argument("--cubepath",
                   default="/Net/Groups/BGI/work_4/scratch/fluxcom/sitecube_proc/"
                           "v0.07/site_cube")
    p.add_argument("--all-sites", action="store_true",
                   help="ignore sites_toRun of the main json")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--step", default="process", choices=["submit", "status", "process"])
    p.add_argument("--check", action="store_true",
                   help="only report the sites missing from each cache")
    p.add_argument("--cache-dir", default="",
                   help="folder that replaces the cache folder of the jsons")
    args = p.parse_args()

    main = json.load(open(args.main_json))
    conf_dir = os.path.dirname(os.path.abspath(args.main_json))
    paths = {k: dataset_json_path(conf_dir, f) for k, f in main["dataset"].items()}
    infos = {k: json.load(open(p)) for k, p in paths.items() if os.path.exists(p)}
    if args.cache_dir:
        move_cache_dir(infos, os.path.abspath(os.path.expanduser(args.cache_dir)))
    datasets = args.datasets or [k for k, v in infos.items() if v.get("cache_path")]
    sites, lats, lons = site_table(args, main)
    print(f"{len(sites)} sites, datasets: {datasets}")

    for name in datasets:
        info = infos[name]
        if not info.get("cache_path"):
            print(f"{name} reads a local source and has nothing to prefetch")
            continue
        try:
            provider = importlib.import_module("extractors." + info["extractor"])
        except ImportError as e:
            if not args.check:
                raise
            print(f"{name}: extractor {info['extractor']} cannot be imported here, {e}")
            provider = None
        if args.check:
            check(name, info, sites, lats, lons, provider)
            continue
        if not hasattr(provider, "prefetch"):
            print(f"{name}: extractor {info['extractor']} has no prefetch")
            continue
        print(f"--- {name} ({info['extractor']}) -> {info['cache_path']}", flush=True)
        provider.prefetch(sites, lats, lons, info, workers=args.workers,
                          step=args.step)


if __name__ == '__main__':
    main()
