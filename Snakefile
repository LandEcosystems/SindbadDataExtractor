# Snakefile of SindbadDataExtractor: the site registry, the product store and the
# cubes, as described in configuration/pipeline.json and utils/pipeline.py.
#
#   snakemake -n                          list what is missing or stale
#   snakemake --cores 20                  run on this node
#   snakemake --profile profiles/slurm    run through SLURM
#   snakemake store                       stage every product, no cubes
#   snakemake --config catalog=<json>     use another catalog, such as a copy with a
#                                         path_map for a local run
#
# Done once by hand, outside the workflow, and only read here:
#   prefetch.py, which fills the caches of the prefetch products, the ERA5 forcing
#     of CLIFF among them. A stage job needs the cache of its product to exist, and
#     adding registry sites needs prefetch.py again for the new sites.
#   the read-only copy of the LaThuile2007 files, with its SOURCE.md.
#
# CLIFF runs inside the workflow: the product <dataset>.<source> of the pycliff block
# of the catalog runs PyCLIFF per site on the hourly CLIFF input of the source, so
# its stage job waits for the cliffin cube of the source and for the forcing.
#
# Each stage job writes STAGED and each assemble job assembled.txt when it ends.
# These markers, not MANIFEST.csv or the site files, are the outputs of the rules,
# since Snakemake removes the outputs of a job before running it again and the
# manifest holds the stage hashes that let a job skip unchanged sites.

import json
import os
import sys

sys.path.insert(0, workflow.basedir)
import utils.pipeline as pl

CATALOG = config.get("catalog", pl.CATALOG)
CAT = pl.load_catalog(CATALOG)
PY = CAT["envs"]
STORE, CUBES = CAT["store"], CAT["cubes_dir"]
LISTS = CAT["site_lists"]
TAG = {"daily": "DD", "hourly": "HR"}
REGISTRY_DONE = os.path.join(STORE, "_markers", "registry.done")
LTL07_SOURCE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    CAT["ltl07"]["daily"]))), "SOURCE.md")
# Cube jobs by name, which is also the folder run_extraction.py writes.
JOBS = {f"{run}.{coll}.{TAG[freq]}": (run, coll, freq, cube)
        for run, coll, freq, cube in pl.cube_jobs(CAT)}
ZARRS = {f"{coll}_{freq}": job for job, (run, coll, freq, cube) in JOBS.items()
         if cube.get("zarr")}


def staged(product, freq):
    return os.path.join(STORE, product, freq, "STAGED")


def caches(entry, freq):
    """
    Prefetch caches of a product from its dataset json after the path_map: its own,
    and for CLIFF the forcing in pix_dir.
    """
    with open(pl.dataset_json(CAT, entry, freq)) as f:
        info = pl.map_paths(json.load(f), CAT["path_map"])
    return [info[k] for k in ("cache_path", "pix_dir") if info.get(k)]


def cliff_input(entry):
    """
    What the CLIFF input of a source needs: the assembled cliffin cube when the
    input is one, else the folder of the input files.
    """
    folder = os.path.dirname(entry["dataset_overrides"]["site_input"])
    job = os.path.dirname(os.path.dirname(folder))
    if os.path.dirname(job) == CUBES.rstrip("/"):
        return os.path.join(job, "assembled.txt")
    return ancient(folder)


def stage_inputs(wc):
    entry = pl.product(CAT, wc.product)
    ins = [REGISTRY_DONE]
    if entry.get("mode") != "legacy":
        ins.append(pl.dataset_json(CAT, entry, wc.freq))
        ins += pl.code_files(CAT, entry, wc.freq)
    ins += [ancient(c) for c in caches(entry, wc.freq)]
    if entry.get("mode") == "pycliff":
        ins.append(cliff_input(entry))
    if entry["dataset"].startswith("ltl07"):
        ins.append(ancient(LTL07_SOURCE))
    return ins


def assemble_inputs(wc):
    run, coll, freq, cube = JOBS[wc.job]
    return [REGISTRY_DONE, pl.recipe_path(cube, freq)] + \
        [staged(p, freq) for p in pl.recipe_products(cube, coll, freq)]


def assemble_args(wc):
    run, coll, freq, cube = JOBS[wc.job]
    return f"{pl.recipe_path(cube, freq)} {coll} {freq} {CUBES} {CAT['fluxcom_cube']} {run}"


def assemble_threads(wc):
    run, coll, freq, cube = JOBS[wc.job]
    from utils.setup_config import get_inp_config
    return int(get_inp_config(pl.recipe_path(cube, freq)).get("njobs_parallel", 1))


wildcard_constraints:
    product="[^/]+",
    freq="daily|hourly",
    job="[^/]+",
    name="[^/]+",


rule all:
    input:
        [os.path.join(CUBES, f"{name}.zarr") for name in ZARRS],
        [os.path.join(CUBES, job, "assembled.txt") for job in JOBS],


rule store:
    input:
        [staged(p, f) for p, e in CAT["products"].items() for f in e["freqs"]],


rule fluxcom_inventory:
    output:
        [os.path.join(LISTS, f"{r}.csv") for r in CAT["releases"]],
        os.path.join(STORE, "_inventory", "v0.07_releases.csv"),
    log:
        os.path.join(STORE, "_logs", "fluxcom_inventory.log"),
    params:
        py=PY["main"],
    shell:
        "{params.py} {workflow.basedir}/fluxcom_inventory.py --catalog {CATALOG} "
        "> {log} 2>&1"


rule ltl07_sites:
    input:
        ancient(LTL07_SOURCE),
    output:
        os.path.join(LISTS, "LaThuile2007.csv"),
    params:
        py=PY["main"],
        src=CAT["ltl07"]["daily"],
    shell:
        "{params.py} {workflow.basedir}/extractors/sitefile_nc.py "
        "--site-list {params.src} {output}"


rule registry:
    input:
        site_input=CAT["site_input"],
        lists=[os.path.join(LISTS, f"{n}.csv")
               for n in CAT["releases"] + ["LaThuile2007"]],
    output:
        REGISTRY_DONE,
    log:
        os.path.join(STORE, "_logs", "registry.log"),
    params:
        py=PY["prefetch"],
        out=CAT["registry"],
        lists=lambda wc, input: " ".join(
            f"--site-list {os.path.basename(p)[:-4]}={p}" for p in input.lists),
        options=" ".join(
            [f"--shuttle --shuttle-dir {CAT['consolidate']['shuttle_dir']}"]
            + [f"--{k.replace('_', '-')} {v}" for k, v in CAT["consolidate"].items()
               if k in ("pft_csv", "fallback_pft", "pft_overrides", "countries")]
            + [f"--coord-overrides {CAT['consolidate']['coord_overrides']}"
               if os.path.exists(CAT["consolidate"].get("coord_overrides", ""))
               else ""]),
    shell:
        # The registry is kept outside the outputs, so that a failed run leaves the
        # earlier one, which the moved-sites report of the log compares with.
        "{params.py} {workflow.basedir}/consolidate_sites.py {input.site_input} "
        "{params.out} {params.options} {params.lists} > {log} 2>&1 && date > {output}"


rule stage:
    input:
        stage_inputs,
    output:
        os.path.join(STORE, "{product}", "{freq}", "STAGED"),
    log:
        os.path.join(STORE, "_logs", "stage_{product}_{freq}.log"),
    threads: 20
    params:
        py=lambda wc: PY[pl.product(CAT, wc.product).get("env", "main")],
    shell:
        "{params.py} {workflow.basedir}/stage_product.py {wildcards.product} "
        "{wildcards.freq} --catalog {CATALOG} --workers {threads} > {log} 2>&1 "
        "&& date > {output}"


rule assemble:
    input:
        assemble_inputs,
    output:
        os.path.join(CUBES, "{job}", "assembled.txt"),
    log:
        os.path.join(CUBES, "_logs", "assemble_{job}.log"),
    threads: assemble_threads
    params:
        py=PY["main"],
        args=assemble_args,
    shell:
        "{params.py} {workflow.basedir}/run_extraction.py {params.args} > {log} 2>&1 "
        "&& date > {output}"


rule zarr:
    input:
        lambda wc: os.path.join(CUBES, ZARRS[wc.name], "assembled.txt"),
    output:
        directory(os.path.join(CUBES, "{name}.zarr")),
    log:
        os.path.join(CUBES, "_logs", "zarr_{name}.log"),
    params:
        py=PY["main"],
        nc_dir=lambda wc: os.path.join(CUBES, ZARRS[wc.name], "data",
                                       JOBS[ZARRS[wc.name]][2]),
        registry=CAT["registry"],
    shell:
        "{params.py} {workflow.basedir}/combine_sites_zarr.py {output} "
        "{params.nc_dir} --registry {params.registry} --overwrite > {log} 2>&1"
