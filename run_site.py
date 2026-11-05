#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The whole workflow for some sites, from what the cluster already holds: the site
registry, the prefetch caches and the sources. It runs the steps of the Snakefile
for those sites only, in their order:
    1. the CLIFF forcing of the sites, prefetched where pixel files are missing
    2. the cliffin cubes, the hourly tower input of PyCLIFF, of the collections
       whose PyCLIFF products read a cube
    3. every store product of every cube of the sites, with extract_product.py
    4. the cubes of the sites, with run_extraction.py for those sites only

usage:
    python run_site.py <site[,site...] | @file | all> [--collections Shuttle]
                       [--freqs daily] [--workers n] [--force] [--no-prefetch]
                       [--list] [--dry-run] [--catalog json]

A file given as @file lists sites, one or more per line, and all stands for every
site of the registry in a cube of the collections. --list prints the sites that
have a cube, one per line, and runs nothing; slurm_scripts/ submits a job per site
from it. By default every cube whose registry flag a site has, at every frequency.
The cliffin cubes needed by the PyCLIFF products always come along. --workers runs
that many sites at once within each product and cube, and the cpus the sites leave
free fit the variables of a site in parallel in PyCLIFF. --no-prefetch leaves
missing forcing missing, for jobs of single sites when the forcing was prefetched
for all sites at once. Each command is printed before it runs, and --dry-run only
prints them.

Nothing here writes the markers of the Snakefile, EXTRACTED and assembled.txt. A
later Snakemake run therefore still extracts and assembles every product and cube,
and skips these sites in the store, whose stage hash is unchanged. extract_product.py
merges its manifest under a lock, so jobs of different sites may run at once. The
registry is read, never rebuilt; run Snakemake once, or its registry rule, for a
new site. The zarr stores are left out, since they combine all sites of a cube.
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


def read_sites(arg):
    """Sites of a comma separated list, or of a file given as @file."""
    if arg.startswith("@"):
        text = open(arg[1:]).read().replace(",", " ")
        names = [s for line in text.splitlines() for s in line.split("#")[0].split()]
    else:
        names = [s.strip() for s in arg.split(",")]
    return list(dict.fromkeys(s for s in names if s))


def cube_flag(cube, collection, freq):
    """Registry flag of the sites of a cube."""
    rec = get_inp_config(pl.recipe_path(cube, freq))
    return fill_version(rec.get("registry", {}), collection).get("flag", "all")


def with_flag(sites, flag, reg):
    """The sites that have the registry flag; a flag the registry lacks has none."""
    try:
        have = set(sites_reg.sites_with(flag, reg))
    except KeyError as e:
        print(f"skipped: {e}", file=sys.stderr, flush=True)
        return []
    return [s for s in sites if s in have]


def cube_jobs(cat, reg, sites, collections, freqs):
    """
    Cube jobs of the sites, the cliffin cubes first, as (run, collection, freq,
    cube, sites). A cliffin cube, hourly only, comes along whenever a cube of its
    collection is run, since the PyCLIFF products of that cube read it.
    """
    jobs = []
    for r, c, f, cube in pl.cube_jobs(cat):
        if collections and c not in collections:
            continue
        here = with_flag(sites, cube_flag(cube, c, f), reg)
        if here:
            jobs.append((r, c, f, cube, here))
    main = [j for j in jobs if j[0] != "cliffin" and (not freqs or j[2] in freqs)]
    with_main = {j[1] for j in main}
    cliff = [j for j in jobs if j[0] == "cliffin"
             and (j[1] in with_main or not freqs or j[2] in freqs)]
    return cliff + main


def products_of(cat, reg, job):
    """Store products a cube reads, with their entries and the sites they cover."""
    run, coll, freq, cube, sites = job
    out = []
    for name in pl.recipe_products(cube, coll, freq):
        if name not in cat["products"]:
            continue
        entry = pl.product(cat, name)
        if freq not in entry["freqs"]:
            continue
        here = with_flag(sites, entry.get("sites", "all"), reg)
        if here:
            out.append((name, entry, here))
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


def run(cmd, dry, env=None):
    print("+ " + " ".join(cmd), flush=True)
    if dry:
        return 0
    return subprocess.run(cmd, cwd=REPO, env=env).returncode


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sites", help="comma separated sites, @file with one per line, or all")
    p.add_argument("--collections", default="",
                   help="comma separated collections, such as Shuttle; default all")
    p.add_argument("--freqs", default="", help="daily, hourly or both; default both")
    p.add_argument("--workers", type=int, default=1,
                   help="sites at once within each product and cube")
    p.add_argument("--catalog", default=pl.CATALOG)
    p.add_argument("--force", action="store_true", help="extract the sites again")
    p.add_argument("--no-prefetch", action="store_true",
                   help="do not prefetch missing CLIFF forcing")
    p.add_argument("--list", action="store_true",
                   help="print the sites that have a cube and run nothing")
    p.add_argument("--dry-run", action="store_true", help="print the commands only")
    a = p.parse_args()

    cat = pl.load_catalog(a.catalog)
    reg = sites_reg.load_registry(cat["registry"])
    sites = list(reg.index) if a.sites == "all" else read_sites(a.sites)
    unknown = [s for s in sites if s not in reg.index]
    if unknown:
        sys.exit(f"{unknown} are not in the registry {cat['registry']}")
    collections = [c for c in a.collections.split(",") if c]
    freqs = [f for f in a.freqs.split(",") if f]
    jobs = cube_jobs(cat, reg, sites, collections, freqs)
    if a.list:
        have = {s for job in jobs for s in job[4]}
        print("\n".join(s for s in sites if s in have))
        return
    if not jobs:
        sys.exit(f"{sites} are in no cube of {collections or 'any collection'}")
    print(f"{len(sites)} sites; cubes " + ", ".join(
        f"{r}.{c}.{f} ({len(s)} sites)" for r, c, f, _, s in jobs), flush=True)
    py = cat["envs"]
    workers = str(max(1, a.workers))
    failed = []

    # 1. CLIFF forcing of the sites
    forcings = {}
    for job in jobs:
        for name, entry, _ in products_of(cat, reg, job):
            if entry.get("mode") == "pycliff":
                found = forcing_dataset(cat, entry, job[2])
                if found:
                    forcings[found[0]] = found
    for key, info, main_json in forcings.values():
        from extractors.cliff_forcing import is_missing, variables
        todo = [s for s in sites if any(
            is_missing(info, v, s, reg.loc[s, "latitude"], reg.loc[s, "longitude"])
            for v in variables(info))]
        if not todo:
            print(f"{key}: pixel files of every site are there", flush=True)
            continue
        if a.no_prefetch:
            print(f"{key}: no pixel files of {todo}, their PyCLIFF products will be "
                  f"missing; prefetch them with prefetch_array.sh", flush=True)
            continue
        site_csv = os.path.join(cat["store"], "_runs", "run_site",
                                f"{key}_{os.getpid()}.csv")
        if not a.dry_run:
            os.makedirs(os.path.dirname(site_csv), exist_ok=True)
            with open(site_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["site", "latitude", "longitude"])
                for s in todo:
                    w.writerow([s, reg.loc[s, "latitude"], reg.loc[s, "longitude"]])
        if run([py["prefetch"], "prefetch.py", main_json, key, "--sites-csv", site_csv,
                "--workers", workers], a.dry_run):
            failed.append(f"prefetch {key}")

    # 2-4. store products and cubes, the cliffin cubes first
    extracted = {}
    for job in jobs:
        r, coll, freq, cube, job_sites = job
        for name, entry, here in products_of(cat, reg, job):
            done = extracted.setdefault((name, freq), set())
            todo = [s for s in here if s not in done]
            if not todo:
                continue
            done.update(todo)
            cmd = [py[entry.get("env", "main")], "extract_product.py", name, freq,
                   "--catalog", a.catalog, "--sites", ",".join(todo),
                   "--workers", workers]
            # the cpus that the sites leave free fit the variables of a site at once
            env = dict(os.environ, PYCLIFF_JOBS=str(max(1, int(workers) // len(todo))))
            if run(cmd + (["--force"] if a.force else []), a.dry_run, env):
                failed.append(f"extract {name} {freq}")
        cmd = [py["main"], "run_extraction.py", pl.recipe_path(cube, freq), coll, freq,
               cat["cubes_dir"], cat["fluxcom_cube"], r, ",".join(job_sites)]
        if run(cmd, a.dry_run):
            failed.append(f"cube {r}.{coll}.{freq}")

    if failed:
        sys.exit(f"failed {failed}; see the logs in "
                 f"{os.path.join(cat['store'], '_runs')}")
    print(f"{len(sites)} sites: done", flush=True)


if __name__ == '__main__':
    main()
