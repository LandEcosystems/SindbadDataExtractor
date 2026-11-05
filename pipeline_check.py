#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
What each step of a run of the site cubes has done, read from the files that the
steps write. Nothing is fetched or extracted.

usage:
    python pipeline_check.py [step ...] [--catalog json] [--products p ...]
                             [--hash] [--errors n]

The steps, in the order of README.md:
    setup     step 0: the files of setup_site_cubes.sh check, the last commit of
              each clone in tools and whether the environment imports its packages
    registry  step 1: the site lists, the registry with the number of sites of each
              flag, its marker, and the files it reads that are newer than it
    prefetch  step 2: prefetch_check.sh, after the number of sites that each
              prefetch product covers
    store     step 4a to 4c: for every product and resolution its EXTRACTED marker
              and the status of its sites in MANIFEST.csv, against the sites of its
              registry flag
    cubes     step 4d and 4e: for every cube its assembled.txt marker, its site
              files against the sites of its flag, and its zarr store
Without a step, all of them run. Each step ends with what to run next when
something is missing or out of date.

--products limits store to some products. --errors sets how many failed sites of a
product are listed, 5 by default.

--hash also computes the stage hash of every site of the products in store, as
extract_product.py does, and counts the sites that the next extraction would
extract again, such as those whose Shuttle zip, C2F flags or Landsat csv changed
after a prefetch. It reads the cache files of every site, so it takes minutes, and
it makes the run folders of the products in the store as extract_product.py does.
Snakemake reruns an extract job only when its inputs change, so such sites need
the --forcerun command that it prints.
"""
import argparse
import datetime
import glob
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils.pipeline as pl  # noqa: E402
import utils.sites as sites_reg  # noqa: E402

STEPS = ["setup", "registry", "prefetch", "store", "cubes"]
TAG = {"daily": "DD", "hourly": "HR"}


def mtime(path):
    return os.path.getmtime(path) if os.path.exists(path) else None


def stamp(t):
    if t is None:
        return "missing"
    return datetime.datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")


def header(text):
    print(f"\n== {text}")


def site_root(cat):
    """The site_cubes folder, two levels above the store."""
    return os.path.dirname(os.path.dirname(cat["store"].rstrip("/")))


def site_lists(cat):
    """Names of the site lists of the registry, as in the Snakefile."""
    made = cat["releases"] + ["LaThuile2007"]
    extra = sorted(os.path.basename(p)[:-4]
                   for p in glob.glob(os.path.join(cat["site_lists"], "*.csv"))
                   if os.path.basename(p)[:-4] not in made)
    return made + extra


def load_registry(cat):
    if not os.path.exists(cat["registry"]):
        return None
    return sites_reg.load_registry(cat["registry"])


def flag_sites(flag, reg):
    """Sites of a registry flag, or None when the registry lacks the flag."""
    try:
        return sites_reg.sites_with(flag, reg)
    except KeyError:
        return None


def check_setup(cat, root):
    header(f"step 0, setup of {root}")
    script = os.path.join(pl.REPO, "site_cubes", "setup_site_cubes.sh")
    subprocess.run(["bash", script, "check"], env={**os.environ, "ROOT": root})
    print("\nclones in tools:")
    for repo in sorted(glob.glob(os.path.join(root, "tools", "*", ".git"))):
        d = os.path.dirname(repo)
        last = subprocess.run(["git", "-C", d, "log", "-1", "--date=short",
                               "--format=%h %ad %s"], capture_output=True, text=True)
        state = subprocess.run(["git", "-C", d, "status", "-sb"],
                               capture_output=True, text=True).stdout.splitlines()
        dirty = len(state) - 1
        print(f"    {os.path.basename(d)}: {last.stdout.strip()}")
        print(f"        {state[0] if state else ''}"
              + (f", {dirty} files changed" if dirty > 0 else ""))
    print(f"this check runs from {pl.REPO}")
    py = cat["envs"]["main"]
    if not os.path.exists(py):
        print(f"\nenvironment: {py} is missing, run "
              "tools/SindbadDataExtractor/sinproc_installer/install_sinproc.sh")
        return
    test = subprocess.run(
        [py, "-c", "import fluxcom.providers, fluxnet_shuttle.shuttle, pycliff, "
                   "rasterio, snakemake, xarray, zarr"], capture_output=True, text=True)
    if test.returncode == 0:
        print(f"\nenvironment: {py} imports fluxcom, fluxnet_shuttle, pycliff, "
              "rasterio, snakemake, xarray and zarr")
    else:
        err = (test.stderr.strip().splitlines() or ["unknown error"])[-1]
        print(f"\nenvironment: {py} fails to import: {err}")


def check_registry(cat):
    header("step 1, site registry")
    print(f"site lists in {cat['site_lists']}:")
    lists = []
    for name in site_lists(cat):
        path = os.path.join(cat["site_lists"], f"{name}.csv")
        lists.append(path)
        n = len(pd.read_csv(path)) if os.path.exists(path) else "-"
        print(f"    {name + '.csv':22s} {n:>5} sites  {stamp(mtime(path))}")
    reg_path = cat["registry"]
    marker = os.path.join(cat["store"], "_markers", "registry.done")
    t_reg = mtime(reg_path)
    print(f"registry: {reg_path}, {stamp(t_reg)}")
    print(f"marker:   {marker}, {stamp(mtime(marker))}")
    if t_reg is None:
        print("next: ./run_workflow.sh registry")
        return
    reg = load_registry(cat)
    print(f"{len(reg)} sites, by flag:")
    for c in [c for c in reg.columns if c.startswith("in_")]:
        print(f"    {c:24s} {int(reg[c].fillna(0).astype(int).sum()):>5}")
    con = cat["consolidate"]
    inputs = {"site input": cat["site_input"]}
    inputs.update({f"site list {os.path.basename(p)}": p for p in lists})
    inputs.update({k: con[k] for k in ("fallback_pft", "pft_overrides",
                                       "coord_overrides") if con.get(k)})
    snaps = sorted(glob.glob(os.path.join(con["shuttle_dir"],
                                          "fluxnet_shuttle_snapshot_*.csv")))
    newer = [(k, p) for k, p in inputs.items() if (mtime(p) or 0) > t_reg]
    if snaps:
        print(f"newest Shuttle listing: {os.path.basename(snaps[-1])}, "
              f"{stamp(mtime(snaps[-1]))}")
    snap_newer = bool(snaps) and mtime(snaps[-1]) > t_reg
    for k, p in newer:
        print(f"newer than the registry: {k}, {p}")
    if snap_newer:
        print("newer than the registry: the Shuttle listing, which Snakemake does not "
              "track")
    if newer:
        print("next: ./run_workflow.sh registry")
    if snap_newer:
        print("next: ./run_workflow.sh registry --forcerun registry, if the new "
              "listing should add sites or move them")
    print(f"log: {os.path.join(cat['store'], '_logs', 'registry.log')}")


def check_prefetch(cat):
    header("step 2, prefetch caches")
    reg = load_registry(cat)
    if reg is None:
        print("no registry yet, run step 1 first")
        return
    by_flag = {}
    for name, e in cat["products"].items():
        if e.get("mode") == "prefetch":
            by_flag.setdefault(e.get("sites", "all"), []).append(name)
    print("sites that the prefetch products need, by registry flag. "
          "prefetch_check.sh counts against all registry sites:")
    for flag, names in by_flag.items():
        sites = flag_sites(flag, reg)
        n = len(sites) if sites is not None else "flag missing in the registry,"
        print(f"    {flag}: {n} sites, for {', '.join(names)}")
    sys.stdout.flush()
    subprocess.run(["bash", os.path.join(pl.REPO, "prefetch_check.sh")],
                   env={**os.environ, "PYTHON": sys.executable,
                        "REGISTRY": cat["registry"]})
    print("next: ./prefetch_all.sh <dataset ...> for the datasets with missing sites")


def changed_sites(cat, name, entry, freq, sites, manifest, reg):
    """Sites whose stage hash differs from that of their manifest row."""
    import contextlib
    import io
    import extract_product as ep
    old = manifest.set_index("site")["stage_hash"].astype(str) \
        if len(manifest) else pd.Series(dtype=str)
    out = []
    # the configuration prints a line per variable, which is not wanted here
    with contextlib.redirect_stdout(io.StringIO()):
        config = ep.build_config(cat, name, entry, freq)
        for site in sites:
            h = ep.stage_hash(entry, config, sites_reg.site_info(site, reg), freq)
            if site in old.index and old[site] != h:
                out.append(site)
    return out


def last_line(message):
    """The last line of a message, such as the error of a traceback."""
    lines = [x.strip() for x in str(message).splitlines() if x.strip()]
    return lines[-1] if lines else ""


def check_store(cat, products=None, hashes=False, n_errors=5):
    header("step 4a to 4c, product store")
    reg = load_registry(cat)
    if reg is None:
        print("no registry yet, run step 1 first")
        return
    no_marker, rerun, failed = [], [], []
    names = sorted(products or cat["products"])
    w = max(len(n) for n in names) + 1
    print(f"{'product':{w}s} {'freq':6s} {'sites':>5} {'ok':>5} {'none':>5} "
          f"{'error':>5} {'stale':>5} {'todo':>5}"
          + (f" {'changed':>7}" if hashes else "") + "  EXTRACTED")
    for name in names:
        entry = pl.product(cat, name)
        for freq in entry["freqs"]:
            d = pl.product_dir(cat, name, freq)
            marker = os.path.join(d, "EXTRACTED")
            sites = flag_sites(entry.get("sites", "all"), reg)
            if sites is None:
                print(f"{name:{w}s} {freq:6s} registry has no flag {entry['sites']}")
                continue
            man_path = os.path.join(d, "MANIFEST.csv")
            man = pd.read_csv(man_path, dtype={"stage_hash": str}) \
                if os.path.exists(man_path) else pd.DataFrame(
                    columns=["site", "status", "stage_hash", "message"])
            man = man[man["site"].isin(sites)]
            counts = man["status"].value_counts()
            todo = len(set(sites) - set(man["site"]))
            row = (f"{name:{w}s} {freq:6s} {len(sites):>5} "
                   + " ".join(f"{int(counts.get(s, 0)):>5}"
                              for s in ("ok", "none", "error", "stale"))
                   + f" {todo:>5}")
            if hashes:
                try:
                    changed = changed_sites(cat, name, entry, freq, sites, man, reg)
                    row += f" {len(changed):>7}"
                    if changed and os.path.exists(marker):
                        rerun.append(marker)
                except Exception as e:
                    row += f" {'failed':>7}"
                    print(f"    hash of {name} {freq} failed: {type(e).__name__}: {e}")
            print(f"{row}  {stamp(mtime(marker))}", flush=True)
            if not os.path.exists(marker):
                no_marker.append(f"{name} {freq}")
            bad = man[man["status"].isin(["error", "stale"])]
            if len(bad):
                failed.append((name, freq, bad))
    print("columns: ok extracted, none the extractor gave no data, error failed "
          "without an earlier file, stale failed and kept the earlier file, todo "
          "not extracted yet" + (", changed extracted again by the next run"
                                 if hashes else ""))
    for name, freq, bad in failed:
        print(f"failed sites of {name} {freq}, log "
              f"{os.path.join(cat['store'], '_logs', f'extract_{name}_{freq}.log')}:")
        for r in bad.head(n_errors).itertuples():
            print(f"    {r.site}: {r.status}, {last_line(r.message)}")
        if len(bad) > n_errors:
            print(f"    and {len(bad) - n_errors} more in MANIFEST.csv")
    if no_marker:
        print(f"{len(no_marker)} products not extracted yet, next: "
              "./run_workflow.sh store, or ./run_workflow.sh for the cubes too")
    if rerun:
        print("next, to extract the changed sites and rebuild their cubes:\n"
              "    ./run_workflow.sh --forcerun " + " ".join(rerun))


def check_cubes(cat):
    header("step 4d and 4e, cubes and zarr stores")
    reg = load_registry(cat)
    if reg is None:
        print("no registry yet, run step 1 first")
        return
    from utils.setup_config import fill_version, get_inp_config
    cubes = cat["cubes_dir"]
    missing, stale = [], []
    print(f"{'cube':34s} {'sites':>5} {'files':>5}  assembled.txt     zarr")
    for run, coll, freq, cube in pl.cube_jobs(cat):
        job = f"{run}.{coll}.{TAG[freq]}"
        rec = get_inp_config(pl.recipe_path(cube, freq))
        flag = fill_version(rec.get("registry", {}), coll).get("flag", "all")
        sites = flag_sites(flag, reg)
        if sites is None:
            print(f"{job:34s} registry has no flag {flag}")
            continue
        files = glob.glob(os.path.join(cubes, job, "data", freq, "*.nc"))
        have = {os.path.basename(f).split(".")[0] for f in files}
        marker = os.path.join(cubes, job, "assembled.txt")
        t_cube = mtime(marker)
        products = [os.path.join(pl.product_dir(cat, p, freq), "EXTRACTED")
                    for p in pl.recipe_products(cube, coll, freq)]
        newer = [p for p in products if (mtime(p) or 0) > (t_cube or 0)]
        zarr = ""
        if cube.get("zarr"):
            z = os.path.join(cubes, f"{coll}_{freq}.zarr")
            zarr = f"{stamp(mtime(z))}"
            if os.path.exists(z):
                try:
                    import xarray as xr
                    zarr += f", {xr.open_zarr(z).sizes.get('site', 0)} sites"
                except Exception as e:
                    zarr += f", unreadable: {type(e).__name__}"
                if t_cube and mtime(z) < t_cube:
                    zarr += ", older than the cube"
        print(f"{job:34s} {len(sites):>5} {len(have):>5}  {stamp(t_cube):16s}  {zarr}")
        if t_cube is None:
            missing.append(job)
        elif newer:
            stale.append((job, newer))
        lacking = sorted(set(sites) - have)
        if t_cube and lacking:
            print(f"    no file for {len(lacking)} sites: {', '.join(lacking[:10])}"
                  + (" ..." if len(lacking) > 10 else ""))
    for job, newer in stale:
        print(f"{job} is older than {len(newer)} of its products, the next run "
              "assembles it again")
    if missing or stale:
        print("next: ./run_workflow.sh")
    print(f"logs: {os.path.join(cubes, '_logs')}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("steps", nargs="*", help=" ".join(STEPS))
    p.add_argument("--catalog", default=pl.CATALOG)
    p.add_argument("--products", nargs="*", default=None)
    p.add_argument("--hash", action="store_true")
    p.add_argument("--errors", type=int, default=5)
    args = p.parse_args()
    unknown = [s for s in args.steps if s not in STEPS]
    if unknown:
        p.error(f"unknown steps {unknown}, choose from {STEPS}")
    cat = pl.load_catalog(args.catalog)
    steps = args.steps or STEPS
    if "setup" in steps:
        check_setup(cat, site_root(cat))
    if "registry" in steps:
        check_registry(cat)
    if "prefetch" in steps:
        check_prefetch(cat)
    if "store" in steps:
        check_store(cat, args.products, args.hash, args.errors)
    if "cubes" in steps:
        check_cubes(cat)
    print("\nthe jobs Snakemake would run: ./run_workflow.sh -n")


if __name__ == '__main__':
    main()
