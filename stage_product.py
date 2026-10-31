#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stage one product of the pipeline catalog into the store: run the extractor of its
dataset for the registry sites that the product covers, and write one NetCDF per site
to <store>/<product>/<freq>/<site>.nc with a MANIFEST.csv of the sites.

usage:
    python stage_product.py <product> <daily|hourly> [--catalog json] [--sites A,B]
                            [--workers n] [--force] [--path-map OLD=NEW ...]

The configuration of the dataset comes from the main json of the base config of the
product, daily, hourly or cliff, with only that dataset selected. The products of
the fluxcom site cube take the cube and FLUXNET_version of their catalog entry and
need the fluxcom package; all others run without it.

A site is staged again only when its stage hash changed. For legacy products, the
fluxcom releases and the LaThuile2007 copy, which are frozen once,
the hash covers the catalog entry, the dates and the coordinates of the site. For the
other products it also covers the expanded variable information of the dataset and
the source of its extractor and of utils.pipeline.SHARED_SOURCES, so that a changed
json or extractor stages the product again. --force stages every site.

Each file is written with zlib compression under a temporary name and renamed when
complete. The attributes of the extractor are kept. Every variable gets
stage_product, and data_path where the dataset json or the fluxcom cube gives one,
which the cube then keeps as the source path. The file gets stage_product,
stage_hash, staged_utc and stage_git, and source_latitude and source_longitude when
the source has coordinates of its own. The manifest records the status of every site,
its time range and source_km, the distance of the source coordinates from those of
the registry, and flags sites more than SOURCE_WARN_KM apart. A site whose staging
fails keeps its earlier file, marked stale in the manifest.
"""
import argparse
import datetime
import hashlib
import importlib
import json
import logging
import os
import subprocess
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import utils.pipeline as pl  # noqa: E402
import utils.setup_config as set_conf  # noqa: E402
import utils.sites as sites_reg  # noqa: E402

SOURCE_WARN_KM = 1.0
MANIFEST_COLUMNS = ["site", "status", "n_vars", "first", "last", "source_km",
                    "stage_hash", "staged_utc", "message"]


def build_config(cat, name, entry, freq):
    """Run configuration with only the dataset of the product selected."""
    base = pl.base_config(cat, entry, freq)
    dataset = entry["dataset"]
    overrides = {"sel_datasets": [dataset], "sel_gapfills": "none",
                 "resample_output": [], "diagnostic_plots": False,
                 "config_dir_path": base["config_dir"],
                 "ancillary_config_dir_path": base["ancillary_dir"]}
    run_dir = os.path.join(cat["store"], "_runs", name)
    config = set_conf.get_exp_configuration(
        base["main"], fn_version=entry.get("FLUXNET_version") or "store",
        temporal_resolution=freq, output_dir_path=run_dir, run_name="stage",
        overrides=overrides)
    config["fluxcom_cube_path"] = entry.get("cube", "")
    config["dataset"] = pl.map_paths(config["dataset"], cat["path_map"])
    config["dataset"][dataset].update(entry.get("dataset_overrides", {}))
    return config


def stage_hash(entry, config, info, freq):
    """
    Hash of everything that changes the staged file of a site. An extractor with
    source_files(info, site, freq), such as the reader of the PyCLIFF products, adds
    the size and time of the files it reads, so that a new run of the source stages
    the site again.
    """
    parts = {"entry": entry, "freq": freq, "site": info,
             "dates": [config["start_date"], config["end_date"]]}
    if entry.get("mode") != "legacy":
        dataset = entry["dataset"]
        info_ds = config["dataset"][dataset]
        parts["dataset"] = info_ds
        parts["code"] = [open(os.path.join(pl.REPO, f)).read()
                         for f in pl.extractor_files(info_ds["extractor"])]
        mod = importlib.import_module("extractors." + info_ds["extractor"])
        if hasattr(mod, "source_files"):
            parts["files"] = [(p, os.path.getsize(p), os.path.getmtime(p))
                              for p in mod.source_files(info_ds, info["site_ID"], freq)]
    text = json.dumps(parts, sort_keys=True, default=str)
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def git_version():
    try:
        out = subprocess.run(["git", "-C", pl.REPO, "describe", "--always", "--dirty"],
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def source_coords(data, entry, site):
    """
    Coordinates of the site in its source, from the source_latitude and
    source_longitude attributes of the data or of its variables, or from the fluxcom
    site cube, or None when the source has none.
    """
    for attrs in [data.attrs] + [data[v].attrs for v in data.data_vars]:
        if "source_latitude" in attrs and "source_longitude" in attrs:
            return float(attrs["source_latitude"]), float(attrs["source_longitude"])
    if entry.get("fluxcom"):
        from fluxcom.site_cube.site_cube_utils import get_site_lat_lon_utc
        lat_lon_utc = get_site_lat_lon_utc(site=site)
        return float(lat_lon_utc[0]), float(lat_lon_utc[1])
    return None


def time_range(data):
    """First and last time step with a value in any variable with time."""
    valid = None
    for v in data.data_vars:
        da = data[v]
        if "time" not in da.dims:
            continue
        has = da.notnull().any(dim=[d for d in da.dims if d != "time"])
        valid = has if valid is None else (valid | has)
    if valid is None or not bool(valid.any()):
        return "", ""
    times = data["time"].values[valid.values]
    return str(times[0])[:19], str(times[-1])[:19]


def write_site(data, path):
    """Write data to path with zlib compression, through a temporary file."""
    for v in data.variables:
        data[v].encoding = {}
    enc = {v: {"zlib": True, "complevel": 4} for v in data.data_vars
           if data[v].dtype.kind in "fiub"}
    tmp = path + ".tmp"
    data.to_netcdf(tmp, encoding=enc)
    os.replace(tmp, path)


def start_log(config, site):
    path = os.path.join(config["output_dir_path"]["log"], f"{site}.log")
    handler = logging.FileHandler(path, "w")
    handler.setFormatter(logging.Formatter(
        '%(levelname)s::%(filename)s::\n\t%(funcName)s::%(lineno)d::\n\t%(message)s'))
    log = logging.getLogger()
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    return handler


def stage_site(job):
    """Stage one site and return its manifest row."""
    site, info, h, config, name, entry, out_dir, git = job
    dataset = entry["dataset"]
    path = os.path.join(out_dir, f"{site}.nc")
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    row = {"site": site, "status": "", "n_vars": 0, "first": "", "last": "",
           "source_km": np.nan, "stage_hash": h, "staged_utc": now, "message": ""}
    handler = start_log(config, site)
    try:
        mod = importlib.import_module(
            "extractors." + config["dataset"][dataset]["extractor"])
        data = mod.extract(dataset=dataset, site_info=info, config=config)
        if data is None:
            if os.path.exists(path):
                os.remove(path)
            row.update(status="none", message="the extractor gave no data")
            return row
        var_info = config["dataset"][dataset]["variables"]
        for v in data.data_vars:
            data[v].attrs["stage_product"] = name
            src_path = var_info.get(v, {}).get("data_path", "") or \
                config["fluxcom_cube_path"]
            if src_path:
                data[v].attrs["data_path"] = src_path
        coords = source_coords(data, entry, site)
        if coords is not None:
            data.attrs["source_latitude"], data.attrs["source_longitude"] = coords
            km = float(haversine_km(info["latitude"], info["longitude"], *coords))
            row["source_km"] = round(km, 3)
            if km > SOURCE_WARN_KM:
                row["message"] = f"source coordinates {km:.1f} km from the registry"
        data.attrs.update({"stage_product": name, "stage_hash": h, "staged_utc": now,
                           "stage_git": git})
        write_site(data, path)
        first, last = time_range(data)
        row.update(status="ok", n_vars=len(data.data_vars), first=first, last=last)
    except Exception as e:
        logging.exception(f"staging {name} failed for {site}")
        status = "stale" if os.path.exists(path) else "error"
        row.update(status=status, message=f"{type(e).__name__}: {e}")
    finally:
        logging.getLogger().removeHandler(handler)
        handler.close()
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("product")
    p.add_argument("freq", choices=["daily", "hourly"])
    p.add_argument("--catalog", default=pl.CATALOG)
    p.add_argument("--sites", default="",
                   help="comma separated sites, staged instead of all of the product")
    p.add_argument("--workers", type=int, default=0,
                   help="parallel sites, by default njobs_parallel of the main json")
    p.add_argument("--force", action="store_true")
    p.add_argument("--path-map", action="append", default=[],
                   help="OLD=NEW, moves paths of the catalog and the dataset jsons")
    args = p.parse_args()

    path_map = dict(a.split("=", 1) for a in args.path_map)
    cat = pl.load_catalog(args.catalog, path_map)
    entry = pl.product(cat, args.product)
    if args.freq not in entry["freqs"]:
        sys.exit(f"{args.product} is staged at {entry['freqs']}, not {args.freq}")
    config = build_config(cat, args.product, entry, args.freq)
    reg = sites_reg.load_registry(cat["registry"])
    sites = sites_reg.sites_with(entry.get("sites", "all"), reg)
    todo = sites
    if args.sites:
        wanted = [s.strip() for s in args.sites.split(",") if s.strip()]
        unknown = [s for s in wanted if s not in sites]
        if unknown:
            sys.exit(f"{unknown} are not sites of {args.product}")
        todo = wanted

    out_dir = pl.product_dir(cat, args.product, args.freq)
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, "MANIFEST.csv")
    old = pd.DataFrame(columns=MANIFEST_COLUMNS)
    if os.path.exists(manifest_path):
        old = pd.read_csv(manifest_path, dtype={"stage_hash": str})
    old = old.set_index("site", drop=False)

    git = git_version()
    jobs, rows = [], {}
    for site in todo:
        info = sites_reg.site_info(site, reg)
        h = stage_hash(entry, config, info, args.freq)
        have = os.path.exists(os.path.join(out_dir, f"{site}.nc"))
        same = site in old.index and old.loc[site, "stage_hash"] == h and \
            old.loc[site, "status"] in ("ok", "none") and \
            (have or old.loc[site, "status"] == "none")
        if same and not args.force:
            rows[site] = old.loc[site].to_dict()
            continue
        jobs.append((site, info, h, config, args.product, entry, out_dir, git))
    print(f"{args.product} {args.freq}: {len(todo)} sites, {len(rows)} unchanged, "
          f"{len(jobs)} to stage", flush=True)

    workers = args.workers or int(config.get("njobs_parallel", 1))
    if workers > 1 and len(jobs) > 1:
        with Pool(min(workers, len(jobs)), maxtasksperchild=1) as pool:
            done = pool.map(stage_site, jobs, chunksize=1)
    else:
        done = [stage_site(j) for j in jobs]
    for row in done:
        if row["status"] == "stale":
            # The earlier file stays, so its row stays too, with the failure noted.
            prev = old.loc[row["site"]].to_dict() if row["site"] in old.index else {}
            row = {**row, **{k: prev.get(k, row[k]) for k in
                             ("n_vars", "first", "last", "source_km", "stage_hash",
                              "staged_utc")}}
        rows[row["site"]] = row

    # Rows of sites of the product that this run did not touch stay as they were.
    for site in sites:
        if site not in rows and site in old.index:
            rows[site] = old.loc[site].to_dict()
    manifest = pd.DataFrame([rows[s] for s in sites if s in rows],
                            columns=MANIFEST_COLUMNS)
    manifest.to_csv(manifest_path, index=False)
    counts = manifest["status"].value_counts().to_dict()
    print(f"written {manifest_path}: {counts}")
    far = manifest[manifest["source_km"] > SOURCE_WARN_KM]
    if len(far):
        print(f"{len(far)} sites whose source coordinates are more than "
              f"{SOURCE_WARN_KM} km from the registry: "
              + ", ".join(f"{r.site} {r.source_km} km" for r in far.itertuples()))
    bad = manifest[manifest["status"].isin(["error", "stale"])]
    for r in bad.itertuples():
        print(f"    {r.site}: {r.status}, {r.message}")


if __name__ == '__main__':
    main()
