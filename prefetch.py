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
whose json sets cache_path. Their jsons are read from the folder of the main json.
Sites and coordinates come from --sites-csv with site, latitude and longitude columns
when given, otherwise from the FLUXCOM site cube of --version, as in
run_extraction.py. --step is only used by vf_modis, whose AppEEARS request is
asynchronous.
"""
import argparse
import importlib
import json
import os
import sys

import pandas as pd

sys.path.append(os.getcwd())


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
    args = p.parse_args()

    main = json.load(open(args.main_json))
    conf_dir = os.path.dirname(os.path.abspath(args.main_json))
    infos = {k: json.load(open(os.path.join(conf_dir, f)))
             for k, f in main["dataset"].items()
             if os.path.exists(os.path.join(conf_dir, f))}
    datasets = args.datasets or [k for k, v in infos.items() if v.get("cache_path")]
    sites, lats, lons = site_table(args, main)
    print(f"{len(sites)} sites, datasets: {datasets}")

    for name in datasets:
        info = infos[name]
        if not info.get("cache_path"):
            print(f"{name} reads a local source and has nothing to prefetch")
            continue
        provider = importlib.import_module("extractors." + info["extractor"])
        if not hasattr(provider, "prefetch"):
            print(f"{name}: extractor {info['extractor']} has no prefetch")
            continue
        print(f"--- {name} ({info['extractor']}) -> {info['cache_path']}", flush=True)
        provider.prefetch(sites, lats, lons, info, workers=args.workers,
                          step=args.step)


if __name__ == '__main__':
    main()
