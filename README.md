# SindbadDataExtractor

Builds site-level data cubes at eddy covariance sites for SINDBAD and other site
studies. For every site it collects the tower fluxes and meteorology, the
meteorology gap-filled with CLIFF, and site properties from in-situ data, global
maps and remote sensing. It writes one NetCDF file per site and source and one zarr
store per source and resolution.

Everything runs on the BGI cluster. All data live below

```
SC=/Net/Groups/BGI/work_3/sindbad/site_cubes
```

`$SC/README.md`, a copy of [site_cubes/README.md](site_cubes/README.md), describes every
file and folder there. This README describes the steps that fill it.

## What comes out

One cube per eddy covariance source and resolution:

| Source | Sites, flag in the registry | Daily | Hourly |
|---|---|---|---|
| FLUXNET Shuttle, the maintained source | `in_OneFluxShuttle` | `sindbad.Shuttle.DD` | `sindbad.Shuttle.HR` |
| LaThuile2007, a frozen legacy copy | `in_LaThuile2007` | `sindbad.LaThuile2007.DD` | `sindbad.LaThuile2007.HR` |
| The 7 releases of the FLUXCOM v0.07 site cube: FLUXNET2015, AmeriFlux, ICOS, OzFLUX, ICOSDrought2018, WarmWinter2020, CCBY | `in_<release>` | `sindbad.<release>.DD` | `sindbad.<release>.HR` |

Each cube is a folder `$SC/output/<cube>/data/<daily|hourly>/` with one file
`<site>.<first year>.<last year>.<daily|hourly>.nc` per site. Each cube is also
combined into a zarr store `$SC/output/<source>_<daily|hourly>.zarr`, with the dimensions
site and Ti.

A cube holds the tower data of its source, the gap-filled meteorology from CLIFF as
`<VAR>_ERA5_gfld` and `<VAR>_CRUJRA_gfld`, and the site properties. The table of
products below lists them.

## The steps

There are five steps, in this order. Steps 1 and 2 are needed once, and again only
when sites or sources change. Run every command on the cluster from
`$SC/tools/SindbadDataExtractor`, the clone of this repository that the workflow
runs.

[pipeline_check.sh](pipeline_check.sh) shows what each step has done. It only reads
files and changes nothing, so you can run it at any time, also while jobs run:

```bash
./pipeline_check.sh                      # every step, in order
./pipeline_check.sh setup                # step 0
./pipeline_check.sh registry             # step 1
./pipeline_check.sh prefetch             # step 2
./pipeline_check.sh store                # step 4a to 4c
./pipeline_check.sh cubes                # step 4d and 4e
./pipeline_check.sh store --hash         # also count the sites the next run extracts again
./pipeline_check.sh store --products shuttle_obs ndvi_landsat
```

Each part ends with a `next:` line with the command to run when something is missing
or out of date. `./run_workflow.sh -n` lists the jobs Snakemake would run.

### Step 0: set up the site_cubes folder (once)

```bash
cd $SC
./setup_site_cubes.sh dirs        # make the folders
./setup_site_cubes.sh tools       # clone SindbadDataExtractor, PyCLIFF and fluxcom into tools/
tools/SindbadDataExtractor/sinproc_installer/install_sinproc.sh   # make tools/env
./rsync_ancillary_data.sh         # run on the local machine: copy ancillary_data, about 58 GB
./setup_site_cubes.sh legacy      # make the merged LaThuile2007 copy
./setup_site_cubes.sh check       # list what is there and what is missing
```

- **Environment.** `tools/env` is the only conda environment. Snakemake, the
  extraction, prefetch, fluxcom and PyCLIFF all run with `tools/env/bin/python`.
  fluxcom and PyCLIFF are installed from their clones in development mode, so a
  `git pull` in a clone takes effect without reinstalling.
- **Ancillary data.** This copies the site list that the registry starts from, the
  PFT and coordinate overrides, the in-situ csvs, the Shuttle zips and the static
  caches of the remote products.
- **LaThuile2007 copy.** `legacy` runs
  [site_cubes/merge_lathuile.py](site_cubes/merge_lathuile.py). It merges two sets of
  LaThuile2007 site files into `ancillary_data/legacy/LaThuile2007/`:
  - the LTL07 files of the CLIFF input, which have the meteorology but no NEE, GPP or
    RECO;
  - the fluxnetBGI v0.07 files, which have the fluxes.

  It then writes `SOURCE.md`, `merge_report.csv` and `MD5SUMS` and makes the folder
  read only. It refuses to run again once `SOURCE.md` exists. Nothing writes
  `SOURCE.md` by hand.

**Check:** `./pipeline_check.sh setup`. It prints three things:
- every file of `setup_site_cubes.sh check`, marked `ok` or `missing`. The kind
  column says how the file is made: `copy` by the rsync, `clone` for a repository,
  `made` by a later step.
- the last commit of each clone in `tools/`, and whether it is behind its remote or
  has local changes.
- whether `tools/env` imports fluxcom, fluxnet_shuttle, pycliff, rasterio,
  snakemake, xarray and zarr.

Files of kind `made` are missing until their step has run.

### Step 1: build the site registry

```bash
./run_workflow.sh registry
```

This runs three Snakemake rules:

1. `fluxcom_inventory` runs [fluxcom_inventory.py](fluxcom_inventory.py).
   - It reads the FLUXCOM v0.07 site cube with the fluxcom package.
   - For each release it writes `ancillary_data/insitu/site_lists/<release>.csv`,
     with the columns `Site_ID, Lat, Lon, UTC`.
   - It also writes `input/store/_inventory/v0.07_releases.csv`, which lists the
     variables each fluxcom provider offers per release.
   - The cube is frozen, so this rule runs once and its outputs are kept.
2. `ltl07_sites` runs `extractors/sitefile_nc.py --site-list`. It reads the daily
   LaThuile2007 copy and writes `site_lists/LaThuile2007.csv`.
3. `registry` runs [consolidate_sites.py](consolidate_sites.py). It merges the
   following into `ancillary_data/insitu/sites_consolidated.csv`, one row per site:
   - `ancillary_data/all_sites_multi_with_ru_flux_and_oneflux.csv`;
   - the site lists of steps 1 and 2, and any other csv in `site_lists/`;
   - the newest Shuttle listing in `prefetch/fluxnet_shuttle/`, which it downloads
     when that folder has none;
   - the PFT and coordinate overrides of `ancillary_data/insitu/`.

   For each site, the registry row holds:
   - the coordinates, which come from the Shuttle first, then the FLUXCOM releases,
     then LaThuile2007, then the input csv;
   - the UTC offset and the PFT;
   - an `in_<source>` flag for each source that holds the site.

   The flags decide which sites each product and cube covers.

**Check:** `./pipeline_check.sh registry`. It prints:
- each site list with its number of sites and date;
- the registry and its marker `input/store/_markers/registry.done`;
- the number of sites of each `in_` flag;
- every file the registry reads that is newer than the registry. This includes the
  Shuttle listing, which Snakemake does not track.

The log `input/store/_logs/registry.log` reports the merged sites, the sites whose
sources disagree on coordinates, and the sites that moved since the previous
registry.

Why the site lists are separate rules and not part of `consolidate_sites.py`:
- Reading the fluxcom cube needs the fluxcom package and is slow. As its own rule it
  runs once, while the registry is rebuilt whenever its inputs change.
- `consolidate_sites.py` only merges csvs of the form `Site_ID, Lat, Lon, UTC`.
  Another source is one more `--site-list NAME=csv`, with no change to the code.

### Step 2: prefetch the remote products

Some products come from remote services or from very large files, so they are
fetched into caches in `$SC/ancillary_data/prefetch/` before the workflow. The
workflow only reads these caches and never fills them.

```bash
./prefetch_check.sh          # what each cache holds and which sites are missing
./prefetch_all.sh            # submit a job for every dataset with missing sites
DRY=1 ./prefetch_all.sh      # only print the sbatch commands
./prefetch_all.sh ndvi_landsat          # only some datasets
```

[prefetch_all.sh](prefetch_all.sh) submits [prefetch_array.sh](prefetch_array.sh)
once per dataset. That script runs [prefetch.py](prefetch.py) for every site of the
registry. Sites already in a cache are skipped, so after a failure you can simply
run it again. A site counts as cached when its id matches and its coordinates are
within 0.01 deg. A site that moved further is therefore fetched again.

| Dataset | What | Job | Cache |
|---|---|---|---|
| `cliff_forcing_era5` | ERA5 hourly at the grid cells around each site, the forcing of PyCLIFF | array of 76, one year each, 1950 to 2025 | `prefetch/cliff_forcing/ERA5.e1/` |
| `cliff_forcing_crujra` | CRUJRA v4.10, 6 hourly, from SINDBAD-TRENDY | array of 10, one variable each | `prefetch/cliff_forcing/CRUJRA.v4.10/` |
| `ndvi_clms300`, `ndvi_clms1km` | Copernicus NDVI, 300 m and 1 km | array, one year each | `prefetch/clms_ndvi_300m/`, `clms_ndvi_1km/` |
| `lai_clms300`, `lai_clms1km` | Copernicus LAI, 300 m and 1 km | array, one year each | `prefetch/clms_lai_300m/`, `clms_lai_1km/` |
| `ndvi_landsat` | Landsat Collection 2 NDVI from Planetary Computer | array of 40, about 33 sites each | `prefetch/landsat_ndvi/`, a csv and a json per site |
| `shuttle_obs` | FLUXNET Shuttle zips | one job | `prefetch/fluxnet_shuttle/` |
| `jung_C2F` | quality flags of Jung et al. 2024, made from the Shuttle zips | one job, after `shuttle_obs` | `prefetch/jung_C2F/` |
| `fa_gami`, `agb_esacci`, `pft_cgls`, `pft_worldcover`, `vf_esacci`, `vf_cgls`, `vf_worldcover`, `vf_fcover`, `wf_chm`, `vf_s2`, `vf_s2sl2p` | footprint means of static maps | one job each | `prefetch/*_sites.csv` and folders |

Datasets that `prefetch_all.sh` leaves out:
- `vf_modis` uses an asynchronous AppEEARS request, so it runs by hand in three
  steps:
  ```bash
  P="$SC/tools/env/bin/python prefetch.py configuration/sindbad/sindbad_extraction_main.json vf_modis --sites-csv $SC/ancillary_data/insitu/sites_consolidated.csv --all-sites"
  $P --step submit     # send the request
  $P --step status     # repeat until every task is done
  $P --step process    # read the result into the cache
  ```
- `shuttle_clim` shares the cache of `shuttle_obs`.
- `era5_API` is paused and no cube reads it. The ERA5 of the cubes comes from the
  cluster files through `cliff_forcing_era5`.

Each task writes its log to `$SC/output/slurm_logs/prefetch_<dataset>_<job>_<task>.log`.

**Check:** `./pipeline_check.sh prefetch`. It first prints how many sites each
prefetch product needs, by registry flag, then the output of `prefetch_check.sh`:
the prefetch jobs still in SLURM and, for each dataset, the sites cached and
missing. How to read it:
- `prefetch_check.sh` counts every dataset against all registry sites. Shuttle and
  Jung data exist only for Shuttle sites, so `shuttle_obs`, `shuttle_clim` and
  `jung_C2F` always show the non-Shuttle sites as missing. Compare their cached count
  with the `in_OneFluxShuttle` count printed first.
- "Only at other coordinates" means the registry moved a site by more than 0.01 deg
  since it was cached. Check the registry coordinates of such sites before
  prefetching them again.

### Step 3 (optional): try a few sites

```bash
python run_site.py DE-Hai --dry-run                  # print the commands, in order
sbatch run_site.sh DE-Hai --collections Shuttle      # only the Shuttle cubes of a site
sbatch --cpus-per-task=32 run_site.sh @sites_01.txt  # a list of sites, 32 at once
```

[run_site.py](run_site.py) runs steps 4a to 4d below for the given sites only. It
first prefetches their CLIFF forcing where it is missing, unless you pass
`--no-prefetch`. It reads the registry and never rebuilds it. It writes no Snakemake
markers, so the full run of step 4 still runs every job. That run then skips these
sites in the store, since nothing about them changed.

**Check:** `./pipeline_check.sh store` and `./pipeline_check.sh cubes`. The sites
appear in the manifests and as cube files, while the markers stay missing until the
full run.

To run every site as its own SLURM job instead of through Snakemake, see
[slurm_scripts/README_SLURM.md](slurm_scripts/README_SLURM.md).

### Step 4: run the workflow

Start it in `screen` or `tmux`, since Snakemake keeps running on the login node until
the end:

```bash
./run_workflow.sh -n                  # dry run: list the jobs that would run
./run_workflow.sh                     # run on 75% of the idle cpus of the partition
./run_workflow.sh 50                  # on 50% of them
./run_workflow.sh 75 store            # only the store, no cubes
TOTAL_CPUS=800 ./run_workflow.sh 50   # on 50% of a fixed pool of 800 cpus
```

[run_workflow.sh](run_workflow.sh) works out the number of cpus to use:
- by default, the given percent of the cpus that `sinfo` reports idle on the
  partition, `big` or `PARTITION`, when the run starts;
- with `TOTAL_CPUS` set, the percent of that pool instead.

It then loads the proxy module and runs Snakemake with the SLURM profile in
[profiles/slurm](profiles/slurm/config.yaml) and that cpu budget. Arguments after the
percent go to Snakemake. `LOCAL=1` runs the jobs on this machine instead; see
"Running without SLURM" below.

The [Snakefile](Snakefile) and the catalog,
[configuration/pipeline.json](configuration/pipeline.json), define these jobs:

**4a. Extract the products into the store.** `extract` runs
[extract_product.py](extract_product.py) once per product and resolution, with up to
20 sites in parallel.
- For each site the product covers, it writes `$SC/input/store/<product>/<freq>/<site>.nc`.
  The product's `sites` flag in the catalog decides which sites it covers.
- It also writes `MANIFEST.csv`, with one row per site and its status: `ok`, `none`
  for no data, `error`, or `stale` for a failure that kept the earlier file.
- When the job ends it writes the marker `EXTRACTED`.
- A site whose stage hash did not change is skipped. The hash covers the site's
  registry row, the dates, the dataset json and the extractor code.

**4b. Assemble the CLIFF input cubes.** `assemble` runs `run_extraction.py` with
the recipes `configuration/cubes/cliffin_*_hourly.json`. It writes
`$SC/output/cliffin.<source>.HR/`, the hourly tower meteorology that PyCLIFF fills
gaps in. LaThuile2007 needs no cliffin cube, since PyCLIFF reads its legacy files
directly.

**4c. Run PyCLIFF.** The products `pycliff_era5.<source>` and
`pycliff_crujra.<source>` are extract jobs as well.
- Each job waits for the cliffin cube of its source and for the forcing cache, then
  runs PyCLIFF per site into `$SC/input/pycliff/`.
- A product is never finer than its forcing. CRUJRA has 6 hour steps, so it only
  fills the daily cubes.
- CRUJRA has no air pressure, so it does not fill PA.
- `qc_threshold` in a pycliff json sets the QC from which a tower value counts as
  tower data. For a day, the QC is the fraction of good hours. The default is 0.85,
  and 0.75 for CRUJRA.

**4d. Assemble the cubes.** `assemble` runs `run_extraction.py` with the recipes
`configuration/cubes/{Shuttle,LaThuile2007,v07}_{daily,hourly}.json`.
- A recipe lists the store products that make up the cube. The job reads them for
  every site of the source and writes `$SC/output/sindbad.<source>.<DD|HR>/`.
- It always rebuilds every site of the cube, then writes the marker `assembled.txt`.

**4e. Combine into zarr.** `zarr` runs [combine_sites_zarr.py](combine_sites_zarr.py)
on each cube with the registry and writes `$SC/output/<source>_<freq>.zarr`. The
`in_` flags of the registry become coordinates on the site dimension.

For one source, the jobs depend on each other in this order:

```
registry ─> extract cliffin_shuttle ─> assemble cliffin.Shuttle.HR ─┐
        │                                                           ├─> extract pycliff_{era5,crujra}.Shuttle ─┐
        │                       prefetch caches (ERA5, CRUJRA) ─────┘                                          │
        └─> extract shuttle_obs, shuttle_clim, jung_C2F, soilgrids, ndvi_landsat, ... ─────────────────────────┴─> assemble sindbad.Shuttle.{DD,HR} ─> Shuttle_{daily,hourly}.zarr
```

**Check:** `./pipeline_check.sh store cubes`, during the run or after it.
- `store` prints one row per product and resolution, with these columns:
  - the sites the product covers;
  - how many sites are `ok`, `none`, `error` and `stale`;
  - `todo`, the sites not extracted yet;
  - the time of its `EXTRACTED` marker.

  It then lists the failed sites with their messages.
- `cubes` prints one row per cube, with these columns:
  - the sites of its flag and its site files;
  - the time of its `assembled.txt` marker;
  - its zarr store, with the number of sites.

  It flags cubes older than their products, zarr stores older than their cube, and
  sites without a file.

Where to look when something fails:

| What | Where |
|---|---|
| Log of an extract job | `$SC/input/store/_logs/extract_<product>_<freq>.log` |
| Status of every site of a product | `$SC/input/store/<product>/<freq>/MANIFEST.csv` |
| Log of an assemble or zarr job | `$SC/output/_logs/assemble_<cube>.log`, `zarr_<name>.log` |
| SLURM log of any job | `$SC/output/slurm_logs/smk-<rule>-<jobid>.log` |

## What makes a step run again

Snakemake reruns a job when one of its inputs is newer than its marker, or when its
rule changed. Inside an extract job, `extract_product.py` then extracts only the
sites whose stage hash changed. The hash covers:
- the site's registry row and the dates;
- the dataset json and the extractor code;
- for the Shuttle, C2F, Landsat and PyCLIFF products, the size and time of the
  files the extractor reads for the site.

| Change | What reruns by itself | What to do |
|---|---|---|
| Extractor code or dataset json | the extract jobs of the product, only for the sites whose hash changed, then the cubes that use it | nothing |
| Cube recipe | that cube and its zarr store | nothing |
| `all_sites_multi_with_ru_flux_and_oneflux.csv`, a site list, a new csv in `site_lists/`, the PFT fallback or overrides, the coordinate overrides | the registry, then every extract job, only for the sites whose registry row changed, then every cube | nothing |
| A new Shuttle listing | nothing, since Snakemake does not track the listing | `./run_workflow.sh --forcerun registry` |
| New data in the Shuttle, C2F or Landsat cache | nothing, since the caches are inputs marked `ancient` | `./pipeline_check.sh store --hash`, then the `--forcerun` command it prints |
| New data in any other prefetch cache, for sites already extracted | nothing, and their hash does not change either | `extract_product.py <product> <freq> --force --sites A,B`, then `--forcerun` its `EXTRACTED` |
| A new version of PyCLIFF | nothing, since only the inputs of a fit are tracked | see "A new version of the code" below |

`--forcerun <marker>` runs that job again and then every job after it, so the cubes
are rebuilt as well. A site extracted while its cache was still empty gets NaN, so
always prefetch new sites before step 4.

## Running without SLURM

Every step also runs as plain processes on one machine:

| Step | Command |
|---|---|
| Registry | `LOCAL=1 ./run_workflow.sh registry` |
| Prefetch, all datasets | `LOCAL=1 ./prefetch_all.sh`: each dataset in one process here, one after the other, with `WORKERS` threads, 8 by default. Logs go to `$LOGS`, by default `$SC/output/slurm_logs/prefetch_<dataset>_local.log` |
| Prefetch, one dataset | `bash prefetch_array.sh ndvi_landsat`, or one part of it with `bash prefetch_array.sh ndvi_landsat --part 3/40` |
| Whole workflow | `LOCAL=1 ./run_workflow.sh 50`: Snakemake runs the jobs here on 50% of the cpus of this machine, or of `TOTAL_CPUS` |
| One product | `python extract_product.py shuttle_obs daily --workers 8`, with `--sites A,B` for some sites |
| Some sites, all their products and cubes | `python run_site.py DE-Hai,US-Ha1 --workers 2` |

On the cluster, take a node first instead of loading the login node, for example:

```bash
srun --partition=big --cpus-per-task=64 --mem=200G --time=24:00:00 --pty bash
cd $SC/tools/SindbadDataExtractor
LOCAL=1 ./run_workflow.sh 100
```

The scripts load the proxy module where the `module` command exists, which the
Shuttle, Copernicus and Landsat prefetch need for internet access.

On a computer outside the cluster, the global maps, the ERA5 files and the fluxcom
cube on `/Net` cannot be read. The steps that read them fail, while the others run
on local copies:
- Give a copy of the catalog a `path_map` that moves the cluster paths, for example
  `"path_map": {"/Net/Groups/BGI/work_3/sindbad/site_cubes": "~/site_cubes"}`.
- Pass the copy with `snakemake --config catalog=<json>`,
  `extract_product.py --catalog <json>` or `pipeline_check.py --catalog <json>`.
- `prefetch.py --cache-dir <folder>` fills the caches in another folder with the
  same layout, for example to copy them to the cluster afterwards.

## Adding sites, sources and updates

Every change below ends with the same two commands: `./run_workflow.sh -n`, to see
what will run, and then `./run_workflow.sh`. A full run after a small change is
cheap, since unchanged sites are skipped in the store. Only the cubes of the
changed sources are rebuilt in full.

### New sites

1. Add the sites in one of three ways:
   - as rows of `all_sites_multi_with_ru_flux_and_oneflux.csv`, with `Site_ID, Lat,
     Lon, UTC, PFT` and a 1 in the `in_<network>` columns that hold them;
   - in a csv of `site_lists/`, see "A new eddy covariance source" below;
   - they come with a new Shuttle listing, see the next section.
2. Rebuild the registry with `./run_workflow.sh registry`. Check with
   `./pipeline_check.sh registry`, and see `registry.log` for the new and moved
   sites.
3. Prefetch the new sites with `./prefetch_all.sh`, which fetches only the missing
   sites, and run the three steps of `vf_modis`. Check with
   `./pipeline_check.sh prefetch`.
4. Run the workflow. Only the new sites are extracted.

### New data of the Shuttle

1. Fetch the update with `./prefetch_all.sh shuttle_obs jung_C2F`.
   - This takes a new listing when the newest one is older than 24 hours.
   - It downloads every site whose product file changed, and the new sites.
   - The C2F flags are then made again for every site whose zip changed.
2. If the new listing adds or moves sites:
   - rebuild the registry with `./run_workflow.sh registry --forcerun registry`;
   - prefetch the other datasets for the new sites with `./prefetch_all.sh`.
   Every extract job then runs again, and `./run_workflow.sh` is enough.
3. Otherwise run `./pipeline_check.sh store --hash --products shuttle_obs
   shuttle_clim jung_C2F cliffin_shuttle`, then the `--forcerun` command it prints.
   Snakemake then runs, in order:
   - the extraction of the changed sites;
   - the cliffin cube;
   - PyCLIFF, for the sites whose tower input changed;
   - the Shuttle cubes and zarr stores.

More Landsat scenes work the same way. Prefetch them, then use
`./pipeline_check.sh store --hash --products ndvi_landsat`.

### A new product for all sites

1. Write the reader `extractors/<name>.py` with
   `extract(dataset, site_info, config)`, which returns the site's data as an
   xarray Dataset. For data that must be fetched first, also give it
   `prefetch(sites, lats, lons, info, ...)`, which fills the cache at the json's
   `cache_path`, and `source_files(info, site, freq)`, so that new cache data
   changes the hash.
2. Write its dataset json in `configuration/ancillary/`, with the extractor, the
   source path or `cache_path`, and the variables. Then add the name of the dataset
   to the `dataset` block of `configuration/sindbad/sindbad_extraction_main.json` and
   `configuration/hourly/hourly_extraction_main.json`.
3. Add the product to `products` of `configuration/pipeline.json`, for example
   `"name": {"mode": "prefetch", "sites": "all", "freqs": ["daily", "hourly"],
   "env": "main"}`.
4. Add it to `sel_datasets` of the cube recipes in `configuration/cubes/` that should
   hold it.
5. For a prefetched product, add a line to `prefetch_all.sh` and prefetch it.
6. Run the workflow. It extracts the product for all sites and rebuilds the cubes
   of the recipes that changed.

### A new eddy covariance source

1. Put its site list in `$SC/ancillary_data/insitu/site_lists/<NAME>.csv`, with the
   columns `Site_ID, Lat, Lon, UTC` and optionally `PFT`.
   - The next run adds the flag `in_<NAME>` to the registry.
   - For coordinates and PFT, the new source ranks after all others. Add
     `in_<NAME>` to `SOURCE_PRIORITY` in `consolidate_sites.py` to change that.
2. Write its readers for the tower fluxes and meteorology as in "A new product",
   step 1. Give them dataset jsons in `configuration/sindbad/` and
   `configuration/hourly/`, plus one in `configuration/cliff/` for the hourly CLIFF
   input.
3. Add its products to `products` of `pipeline.json` with
   `"sites": "in_<NAME>"`.
4. Copy `configuration/cubes/Shuttle_{daily,hourly}.json` to
   `<NAME>_{daily,hourly}.json`. In the copy, set:
   - `FLUXNET_version` to `<NAME>`;
   - the `registry` flag to `in_<NAME>`;
   - the tower products in `sel_datasets`.

   Add an entry for the cube to `cubes` of `pipeline.json`.
5. For CLIFF gap filling, make its hourly CLIFF input in one of two ways:
   - a cliffin recipe as `cliffin_Shuttle_hourly.json`, with an entry with
     `"run": "cliffin"` in `cubes`;
   - files of its own.

   Then add the source to `pycliff.sources` of `pipeline.json` with its flag and
   `site_input` pattern. The products `pycliff_era5.<NAME>` and
   `pycliff_crujra.<NAME>` then exist.
6. Check one site first with `python run_site.py <site> --collections <NAME>
   --dry-run`, then run the workflow.

### A new forcing for CLIFF

1. Add a `cliff_forcing` dataset json that writes the forcing in the layout of
   `prefetch/cliff_forcing/ERA5.e1/`. Add a `pycliff_<forcing>.json` like
   `pycliff_era5.json` in `configuration/sindbad/` and `configuration/hourly/`, and
   a source json in PyCLIFF.
2. Add the dataset to `pycliff.datasets` of `pipeline.json`, and to
   `pycliff.freqs` when it is coarser than hourly.
3. Add it to `sel_datasets` and `store.products` of the cube recipes, as
   `pycliff_crujra` is there.
4. Add a line to `prefetch_all.sh`, prefetch it, and run the workflow.

### A new version of the code

1. Run `git -C $SC/tools/<repo> pull`. The environment runs the new code without
   reinstalling.
2. `./run_workflow.sh -n` shows what the change reruns. Changed extractors, dataset
   jsons and recipes rerun by themselves.
3. A new PyCLIFF is not tracked. To fit the sites again, remove their markers below
   `$SC/input/pycliff/sindbad_inputs/<source>.<DD|HR>/<forcing>/<site>.txt`, run
   `extract_product.py pycliff_era5.<source> <freq> --force`, and then
   `./run_workflow.sh --forcerun` with the `EXTRACTED` marker of that product.

## Products

[configuration/pipeline.json](configuration/pipeline.json) lists every product, with
the sites it covers and its resolutions. The dataset jsons in `configuration/sindbad`,
`configuration/hourly`, `configuration/cliff` and `configuration/ancillary` give the
source paths and variables. The extractors are in [extractors/](extractors/).

| Group | Products | Sites |
|---|---|---|
| Tower data, Shuttle | `shuttle_obs`, `shuttle_clim`, `jung_C2F` (QC of Jung et al. 2024), `cliffin_shuttle` | `in_OneFluxShuttle` |
| Tower data, LaThuile2007 | `ltl07_obs`, `ltl07_clim` | `in_LaThuile2007` |
| Tower data, FLUXCOM v0.07 releases | `bgi_FLUXCOM_obs`, `bgi_FLUXCOM_clim`, `bgi_FLUXCOM_surface`, `jung_qc`, `merged_qc`, `nelson_tea`, `cliffin`, each as `<product>.<release>` | `in_<release>` |
| Gap-filled meteorology | `pycliff_era5.<source>`, `pycliff_crujra.<source>` | flag of the source |
| Global maps on the cluster | `soilgrids`, `sd_soilgrids`, `pft`, `vf_pft`, `vf_hilda`, `tf_hansen`, `fa_gami`, `atmCO2` | all |
| Remote products, prefetched | `pft_cgls`, `pft_worldcover`, `vf_esacci`, `vf_cgls`, `vf_worldcover`, `vf_modis`, `vf_fcover`, `wf_chm`, `vf_s2`, `vf_s2sl2p`, `ndvi_clms300`, `ndvi_clms1km`, `lai_clms300`, `lai_clms1km`, `ndvi_landsat`, `agb_esacci` | all |
| In-situ csvs | `sd_insitu`, `vf_insitu`, `tf_insitu`, `agb` | all |
| Computed | `swinpot_calc`, the potential shortwave radiation | all |

The FLUXCOM releases are frozen. Their products are extracted once with the fluxcom
package and never change, since the v0.07 site cube does not change.

## Other parts of the repository

| Path | What it is |
|---|---|
| [validation/shuttle_vs_fluxnet2015](validation/shuttle_vs_fluxnet2015) | Comparison of the Shuttle data with FLUXNET2015, with a dashboard |
| [slurm_scripts/](slurm_scripts/) | One `run_site.py` job per site, the alternative to Snakemake |
| [juicer/](juicer/) | Inventory of the global data on the cluster |
| `run_extraction.sh`, `run_extraction_cliff.sh`, `run_extraction_hourly.sh`, `check_sites_v006.py` | The runs from before the workflow, with fixed paths to old environments. Not used by the workflow |

# Contact

Sujan Koirala <skoirala@bgc-jena.mpg.de> (on the tool and data), Nuno Carvalhais
<ncarvalhais@bgc-jena.mpg.de> (on the data)
