# Processing tool for collecting site level data from various sources

A repository for the preprocessor of data to be used in SINDBAD MDI framework, and other research at the FLUXNET site scales


# Components:

The tool put together data from 3 main sources
- FLUXCOM data extractor that provides 
    - the fluxnet datasets LaThuile2007 and FLUXNET2015 releases
    - remote sensing cutouts from fluxnet sites

- Site level extractions of biomass and forest ages from Simon. These include the FLUXNET BADM as well as the GlobBiomass satellite retrievals.
- Site level extraction of land and soil characteristics from FIDC that includes soil properties among others
- Land cover properties from HILDA+ and Hansen
- CLIFF downscaling of several global dataset to gapfill the tower measurements

# Pipeline:

Data go through three layers, driven by Snakemake (`Snakefile`) and described in
`configuration/pipeline.json`:

1. **Site registry**: `consolidate_sites.py` writes `sites_consolidated.csv`, one row per
   site with its coordinates, UTC offset, PFT and `in_<network>` flags. The flags of the
   eddy covariance sources come from site lists made from the sources themselves; the
   FLUXNET Shuttle comes first for coordinates and PFT.
2. **Product store**: `stage_product.py <product> <freq>` runs the extractor of one dataset
   for its registry sites and writes `<store>/<product>/<freq>/<site>.nc` with a
   `MANIFEST.csv`. Unchanged sites are skipped. Only this layer reads sources.
3. **Cubes**: `run_extraction.py` with a recipe of `configuration/cubes/` reads the store
   (`extractors/store_reader.py`) and writes one cube per eddy covariance source, which
   `combine_sites_zarr.py` turns into a site x Ti zarr store.

The data live below `/Net/Groups/BGI/work_3/sindbad/site_cubes`, described in
`site_cubes/README.md`:
- `tools/`: git clones of this repository and of PyCLIFF, which the workflow runs.
- `ancillary_data/`: the sources, with `insitu/` (site registry, site lists, in-situ
  csvs and overrides), `reference/`, `prefetch/` (the caches of the remote products and
  the CLIFF forcing) and `legacy/` (the LaThuile2007 copy).
- `input/`: what the cubes are made of, `store/` with one file per product and site,
  and `pycliff/` with the PyCLIFF runs of each site.
- `output/`: the cubes, as site files and zarr stores.

Done once by hand, before Snakemake:
- `prefetch.py` fills the caches of the remote products, as before. The NDVI and LAI
  products of Copernicus (`ndvi_clms300`, `ndvi_clms1km`, `lai_clms300`, `lai_clms1km`)
  and of Landsat (`ndvi_landsat`) are large, so `prefetch_array.sh` splits them over a
  SLURM array, for example `sbatch --array=1-40 prefetch_array.sh ndvi_landsat`.
- The fluxcom package is legacy. `fluxcom_inventory.py` lists the releases of the
  v0.07 site cube and writes their site lists. Snakemake then stages each release
  once in the fluxcom environment.
- LaThuile2007 is a read-only copy of the fluxnetBGI2021 site cube iter-5 files, read by
  `extractors/sitefile_nc.py`. Write a `SOURCE.md` next to the copy.
- The forcings of the CLIFF gap filling, `cliff_forcing_era5` and
  `cliff_forcing_crujra`, are prefetched as well: `prefetch_array.sh` runs
  `pycliff-extract` once per variable for all registry sites. ERA5 comes from the
  hourly files of the cluster, CRUJRA v4.10 from the 6 hourly `trendy_met.zarr` of
  SINDBAD-TRENDY. Another forcing only needs another `cliff_forcing` json and a
  source json in PyCLIFF.

CLIFF then runs inside Snakemake: the products `pycliff_era5.<source>` and
`pycliff_crujra.<source>` run PyCLIFF per site on the hourly tower meteorology of each
eddy covariance source, and the cubes get the gap-filled meteorology as
`<VAR>_ERA5_gfld` and `<VAR>_CRUJRA_gfld`. No product is finer than its forcing:
CRUJRA has 6 hour steps, so it fills the daily cubes only, and it has no air
pressure, so it does not fill PA. `qc_threshold` of a pycliff json sets the QC from
which a tower value, and a day as its fraction of good hours, counts as tower data:
0.85 by default, 0.75 for CRUJRA.

```
./run_workflow.sh -n                  # what is missing or stale
./run_workflow.sh                     # run through SLURM on 75% of the idle cpus
./run_workflow.sh 50                  # on 50% of them
TOTAL_CPUS=800 ./run_workflow.sh 50   # on 50% of a fixed pool of 800 cpus
```

`run_workflow.sh` takes the percent, 75 by default, of the cpus that sinfo reports
idle on the partition, `big` or `PARTITION`, when it starts, or of `TOTAL_CPUS` when
that is set. The share stays fixed for the run. It loads the proxy module for
internet access and runs Snakemake with the SLURM profile and `--resources cpus=N`;
stage and assemble jobs ask for their threads as cpus, and no job asks for more than
N. A stage job runs up to 20 sites at once, so 400 cpus run 20 stage jobs side by
side. Snakemake runs on the node it is started from until the workflow ends, so
start it in screen or tmux.

`run_site.py` runs the same steps for some sites, to check them before a full run:
their CLIFF forcing when missing, the cliffin cubes, every store product and every
cube of the sites. It writes no Snakemake markers, so a later full run still runs
everything and skips the sites where nothing changed. `run_site.sh` submits it with
the proxy loaded and one worker per cpu. Jobs of different sites can run at once,
since the store manifests are merged under a lock.

```
python run_site.py DE-Hai --dry-run                  # the commands, in order
sbatch run_site.sh DE-Hai --collections Shuttle      # only the Shuttle cubes
sbatch --cpus-per-task=32 run_site.sh @sites_01.txt  # a list of sites, 32 at once
```

`slurm_scripts/` submits one `run_site.py` job per site within a share of the free
cpus and memory, 75% by default, refilled by a watcher in tmux, as described in
`slurm_scripts/README_SLURM.md`.

```
./slurm_scripts/run_sites.sh --dry-run               # what would be submitted
./slurm_scripts/run_sites.sh                         # the sites of slurm_scripts/env.sh
```

# Contact:

Sujan Koirala <skoirala@bgc-jena.mpg.de> (on the tool and data), Nuno Carvalhais <ncarvalhais@bgc-jena.mpg.de> (on the data)
