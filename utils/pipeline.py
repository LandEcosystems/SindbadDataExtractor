#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The pipeline catalog, configuration/pipeline.json. It says where the registry, the
store and the cubes live, which python runs each environment, which products the
store holds and which cubes are assembled from them. stage_product.py and the
Snakefile read it through load_catalog.

A product is one dataset of the main jsons staged per site, with its mode, the
registry flag of its sites, its temporal resolutions and its environment. The
products of the FLUXCOM v0.07 releases are made from releases and release_datasets
as <prefix>.<release>, for example bgi_FLUXCOM_obs.FLUXNET2015.

path_map replaces the start of every path of the catalog and of the dataset jsons,
for example to run on a local copy of the cluster folders:
    "path_map": {"/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data":
                 "~/research/RnD/data/sindbad_ancillary_data"}
"""
import copy
import json
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CATALOG = os.path.join(REPO, "configuration", "pipeline.json")
# Modules that every extractor uses, part of the stage hash of the products that are
# not legacy.
SHARED_SOURCES = ["extractors/BasexTractor.py", "extractors/footprint_utils.py",
                  "utils/shared_utils.py"]


def map_paths(obj, path_map):
    """
    Copy of obj with every string that starts with a key of path_map moved, by the
    longest key that matches.
    """
    if isinstance(obj, dict):
        return {k: map_paths(v, path_map) for k, v in obj.items()}
    if isinstance(obj, list):
        return [map_paths(v, path_map) for v in obj]
    if isinstance(obj, str):
        for old in sorted(path_map, key=len, reverse=True):
            if obj.startswith(old):
                return os.path.expanduser(path_map[old]) + obj[len(old):]
    return obj


def release_products(cat):
    """Products of the FLUXCOM v0.07 releases, frozen once from the site cube."""
    out = {}
    for release in cat.get("releases", []):
        for prefix, spec in cat.get("release_datasets", {}).items():
            out[f"{prefix}.{release}"] = {
                "dataset": spec.get("dataset", prefix), "mode": "legacy",
                "fluxcom": True, "FLUXNET_version": release,
                "cube": cat["fluxcom_cube"], "sites": f"in_{release}",
                "freqs": spec["freqs"], "base": spec.get("base", ""), "env": "main"}
    return out


def pycliff_sources(cat):
    """
    Eddy covariance sources that CLIFF gap fills, with their registry flag and the
    pattern of their hourly CLIFF input files. A FLUXCOM v0.07 release reads its
    cliffin cube; other sources are given in the pycliff block of the catalog.
    """
    sources = {r: {"flag": f"in_{r}", "site_input": os.path.join(
        cat["cubes_dir"], f"cliffin.{r}.HR", "data", "hourly", "{site}.*.hourly.nc")}
        for r in cat.get("releases", [])}
    sources.update(cat.get("pycliff", {}).get("sources", {}))
    return sources


def pycliff_products(cat):
    """
    Products <dataset>.<source> that run CLIFF per site for every dataset of the
    pycliff block, one per forcing, and every source.
    """
    pc = cat.get("pycliff")
    if not pc:
        return {}
    return {f"{d}.{s}": {
        "dataset": d, "mode": "pycliff", "sites": src["flag"],
        "freqs": ["daily", "hourly"], "env": "main",
        "dataset_overrides": {"fn": s, "site_input": src["site_input"],
                              "sites_csv": cat["registry"]}}
        for d in pc.get("datasets", []) for s, src in pycliff_sources(cat).items()}


def load_catalog(path=CATALOG, path_map=None):
    """
    The catalog with the products of the releases and of PyCLIFF added and its paths
    moved by path_map.
    """
    with open(path) as f:
        cat = json.load(f)
    cat["path_map"] = {**cat.get("path_map", {}), **(path_map or {})}
    cat["products"] = {**release_products(cat), **pycliff_products(cat),
                       **cat.get("products", {})}
    return map_paths(cat, cat["path_map"])


def product(cat, name):
    """Catalog entry of a product, with its dataset and base config filled in."""
    if name not in cat["products"]:
        raise KeyError(f"{name} is not a product of the catalog")
    entry = copy.deepcopy(cat["products"][name])
    entry.setdefault("dataset", name)
    return entry


def base_config(cat, entry, freq):
    """Main json, config folder and ancillary folder of a product at freq."""
    base = cat["base_configs"][entry.get("base") or freq]
    return {k: os.path.join(REPO, v) for k, v in base.items()}


def product_dir(cat, name, freq):
    return os.path.join(cat["store"], name, freq)


def dataset_json(cat, entry, freq):
    """
    Path of the dataset json of a product at freq, as get_exp_configuration finds
    it.
    """
    from utils.shared_utils import dataset_json_path
    base = base_config(cat, entry, freq)
    with open(base["main"]) as f:
        names = json.load(f)["dataset"]
    return dataset_json_path(base["config_dir"], names[entry["dataset"]],
                             base["ancillary_dir"])


def extractor_files(extractor):
    """Repository paths of an extractor and of the modules every extractor uses."""
    return [f"extractors/{extractor}.py"] + SHARED_SOURCES


def code_files(cat, entry, freq):
    """Extractor modules a product runs, the code part of its stage hash."""
    with open(dataset_json(cat, entry, freq)) as f:
        extractor = json.load(f)["extractor"]
    return [os.path.join(REPO, p) for p in extractor_files(extractor)]


def recipe_path(cube, freq):
    return os.path.join(REPO, cube["recipe"].format(freq=freq))


def recipe_products(cube, collection, freq):
    """Store products that the recipe of a cube reads for one collection."""
    from utils.setup_config import fill_version, get_inp_config, get_selected_list
    rec = get_inp_config(recipe_path(cube, freq))
    store = fill_version(rec.get("store", {}), collection)
    sel = copy.deepcopy(rec["sel_datasets"])
    datasets = get_selected_list(sel, list(rec["dataset"]))
    return [store.get("products", {}).get(d, d) for d in datasets]


def cube_jobs(cat):
    """
    Run, collection, freq and cube spec of every cube to assemble. collections
    "releases" stands for all FLUXCOM v0.07 releases.
    """
    jobs = []
    for cube in cat.get("cubes", []):
        collections = cube["collections"]
        if collections == "releases":
            collections = cat["releases"]
        for collection in collections:
            for freq in cube["freqs"]:
                jobs.append((cube["run"], collection, freq, cube))
    return jobs
