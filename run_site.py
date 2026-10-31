#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The whole workflow for one site, from what the cluster already holds: the site
registry, the prefetch caches and the sources. It runs the steps of the Snakefile
for that site only, in their order:
    1. the CLIFF forcing of the site, prefetched when its pixel files are missing
    2. the cliffin cubes, the hourly tower input of PyCLIFF, of the collections
       whose PyCLIFF products read a cube
    3. every store product of every cube of the site, with stage_product.py --sites
    4. the cubes of the site, with run_extraction.py for that site only

usage:
    python run_site.py <site> [--collections Shuttle,LaThuile2007] [--freqs daily]
                       [--force] [--dry-run] [--catalog json]

By default every cube whose registry flag the site has, at every frequency. The
cliffin cubes needed by the PyCLIFF products always come along. Each command is
printed before it runs, and --dry-run only prints them.

Nothing here writes the markers of the Snakefile, STAGED and assembled.txt. A
later Snakemake run therefore still stages and assembles every product and cube,
and skips this site in the store, whose stage hash is unchanged. The registry is
read, never rebuilt; run Snakemake once, or its registry rule, for a new site.
The zarr stores are left out, since they combine all sites of a cube.
"""
import argparse
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils.pipeline as pl  # noqa: E402
import utils.sites as sites_reg  # noqa: E402
from utils.setup_config import fill_version, get_inp_config  # noqa: E402

REPO = pl.REPO


def cube_flag(cube, collection, freq):
    """Registry flag of the sites of a cube."""
    rec = get_inp_config(pl.recipe_path(cube, freq))
    return fill_version(rec.get("registry", {}), collection).get("flag", "all")


def has_flag(site, flag, reg):
    """Whether the site has the registry flag; a flag the registry lacks is no."""
    try:
        return site in sites_reg.sites_with(flag, reg)
    except KeyError as e:
        print(f"skipped: {e}", flush=True)
        return False


def cube_jobs(cat, reg, site, collections, freqs):
    """
    Cube jobs of the site, the cliffin cubes first, as (run, collection, freq,
    cube). A cliffin cube, hourly only, comes along whenever a cube of its
    collection is run, since the PyCLIFF products of that cube read it.
    """
    jobs = [(r, c, f, cube) for r, c, f, cube in pl.cube_jobs(cat)
            if (not collections or c in collections)
            and has_flag(site, cube_flag(cube, c, f), reg)]
    main = [j for j in jobs if j[0] != "cliffin" and (not freqs or j[2] in freqs)]
    with_main = {j[1] for j in main}
    cliff = [j for j in jobs if j[0] == "cliffin"
             and (j[1] in with_main or not freqs or j[2] in freqs)]
    return cliff + main


def products_of(cat, reg, site, job):
    """Store products a cube reads that cover the site, with their entries."""
    run, coll, freq, cube = job
    out = []
    for name in pl.recipe_products(cube, coll, freq):
        if name not in cat["products"]:
            continue
        entry = pl.product(cat, name)
        if freq not in entry["freqs"]:
            continue
        if has_flag(site, entry.get("sites", "all"), reg):
            out.append((name, entry))
    return out


def forcing_dataset(cat, entry, freq):
    """
    Prefetch dataset of the CLIFF forcing of a PyCLIFF product, the cliff_forcing
    json whose cache_path is the pix_dir of the product, with that json's main json.
    """
    with open(pl.dataset_json(cat, entry, freq)) as f:
        pix_dir = pl.map_paths(json.load(f), cat["path_map"]).get("pix_dir")
    base = pl.base_config(cat, entry, freq)
    with open(base["main"]) as f:
        names = json.load(f)["dataset"]
    from utils.shared_utils import dataset_json_path
    for key, fname in names.items():
        path = dataset_json_path(base["config_dir"], fname, base["ancillary_dir"])
        if not os.path.exists(path):
            continue
        with open(path) as f:
            info = pl.map_paths(json.load(f), cat["path_map"])
        if info.get("extractor") == "cliff_forcing" and info.get("cache_path") == pix_dir:
            return key, info, base["main"]
    return None


def run(cmd, dry):
    print("+ " + " ".join(cmd), flush=True)
    if dry:
        return 0
    return subprocess.run(cmd, cwd=REPO).returncode


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("site")
    p.add_argument("--collections", default="",
                   help="comma separated collections, such as Shuttle; default all")
    p.add_argument("--freqs", default="", help="daily, hourly or both; default both")
    p.add_argument("--catalog", default=pl.CATALOG)
    p.add_argument("--force", action="store_true", help="stage the site again")
    p.add_argument("--dry-run", action="store_true", help="print the commands only")
    a = p.parse_args()

    cat = pl.load_catalog(a.catalog)
    reg = sites_reg.load_registry(cat["registry"])
    if a.site not in reg.index:
        sys.exit(f"{a.site} is not in the registry {cat['registry']}")
    collections = [c for c in a.collections.split(",") if c]
    freqs = [f for f in a.freqs.split(",") if f]
    jobs = cube_jobs(cat, reg, a.site, collections, freqs)
    if not jobs:
        sys.exit(f"{a.site} is in no cube of {collections or 'any collection'}")
    print(f"{a.site}: cubes " + ", ".join(f"{r}.{c}.{f}" for r, c, f, _ in jobs), flush=True)
    py = cat["envs"]
    failed = []

    # 1. CLIFF forcing of the site
    site_csv = os.path.join(cat["store"], "_runs", "run_site", f"{a.site}.csv")
    forcings = {}
    for job in jobs:
        for name, entry in products_of(cat, reg, a.site, job):
            if entry.get("mode") == "pycliff":
                found = forcing_dataset(cat, entry, job[2])
                if found:
                    forcings[found[0]] = found
    for key, info, main_json in forcings.values():
        from extractors.cliff_forcing import is_missing, variables
        if not any(is_missing(info, v, a.site) for v in variables(info)):
            print(f"{key}: pixel files of {a.site} are there", flush=True)
            continue
        if not a.dry_run:
            os.makedirs(os.path.dirname(site_csv), exist_ok=True)
            row = reg.loc[a.site]
            with open(site_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["site", "latitude", "longitude"])
                w.writerow([a.site, row["latitude"], row["longitude"]])
        if run([py["prefetch"], "prefetch.py", main_json, key, "--sites-csv", site_csv],
               a.dry_run):
            failed.append(f"prefetch {key}")

    # 2-4. store products and cubes, the cliffin cubes first
    staged = set()
    for job in jobs:
        r, coll, freq, cube = job
        for name, entry in products_of(cat, reg, a.site, job):
            if (name, freq) in staged:
                continue
            staged.add((name, freq))
            cmd = [py[entry.get("env", "main")], "stage_product.py", name, freq,
                   "--catalog", a.catalog, "--sites", a.site, "--workers", "1"]
            if run(cmd + (["--force"] if a.force else []), a.dry_run):
                failed.append(f"stage {name} {freq}")
        cmd = [py["main"], "run_extraction.py", pl.recipe_path(cube, freq), coll, freq,
               cat["cubes_dir"], cat["fluxcom_cube"], r, a.site]
        if run(cmd, a.dry_run):
            failed.append(f"cube {r}.{coll}.{freq}")

    if failed:
        sys.exit(f"{a.site}: failed {failed}; see the logs in "
                 f"{os.path.join(cat['store'], '_runs')}")
    print(f"{a.site}: done", flush=True)


if __name__ == '__main__':
    main()
