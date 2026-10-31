# SINDBAD site cubes

Site-level data for SINDBAD and other work at FLUXNET sites, made by
[SindbadDataExtractor](https://github.com/LandEcosystems/SindbadDataExtractor) with its
Snakemake workflow. Contact: Sujan Koirala, skoirala@bgc-jena.mpg.de.

```
site_cubes/
├── README.md, setup_site_cubes.sh, rsync_ancillary_data.sh
├── tools/                the code: git clones of SindbadDataExtractor and PyCLIFF
├── ancillary_data/       the sources: copied by hand or fetched once
│   ├── all_sites_multi_with_ru_flux_and_oneflux.csv
│   ├── insitu/           site lists, registry, in-situ csvs, overrides
│   │   └── site_lists/
│   ├── reference/        maps used to check site coordinates
│   ├── prefetch/         caches of the remote products
│   └── legacy/           read-only copy of the LaThuile2007 site files
├── input/                what the cubes are made of, written by the workflow
│   ├── store/            one file per product and site
│   └── pycliff/          the PyCLIFF runs of each site
└── output/               the cubes, as site files and zarr stores
```

`tools` holds the code, `ancillary_data` everything the workflow reads, and `input`
and `output` everything it writes. Nothing in `input` or `output` is edited by hand,
and both can be rebuilt from `tools`, `ancillary_data` and the global data on /Net.

## Setting it up

```bash
cd /Net/Groups/BGI/work_3/sindbad/site_cubes
./setup_site_cubes.sh dirs      # make the folders
./setup_site_cubes.sh tools     # clone SindbadDataExtractor and PyCLIFF
# copy the files of the next section
./setup_site_cubes.sh legacy    # copy the LaThuile2007 files on the cluster
./setup_site_cubes.sh check     # list what is there and what is missing
```

`check` lists each file the workflow reads, marked `copy` when it is copied by hand,
`clone` for the repositories and `made` when a tool writes it.

## What goes where

### tools/

| Path | Content |
|---|---|
| `SindbadDataExtractor/` | This workflow. Snakemake, `prefetch.py` and `prefetch_array.sh` run from here |
| `PyCLIFF/` | The CLIFF gap filling, run by the workflow as `pycliff_dir` of its configs |

Update a tool with `git -C tools/<repo> pull`, then `snakemake -n` lists what the new
code changes. PyCLIFF runs with the python of its environment,
`/Net/Groups/BGI/scratch/skoirala/envs/pycliff`, which only provides the packages it
needs. The workflow puts `tools/PyCLIFF` first on the python path, so the code of
this clone is what runs. PyCLIFF writes its runs to `input/pycliff`, never into its
clone.

### ancillary_data/ (copied by hand)

Everything except `legacy` comes from the local `sindbad_ancillary_data` folder,
which has the same layout:

```bash
SRC=~/research/RnD/data/sindbad_ancillary_data
DST=<user>@<cluster>:/Net/Groups/BGI/work_3/sindbad/site_cubes/ancillary_data
rsync -av --exclude '*.lock' --exclude 'caches_before_*' \
    --exclude era5_API --exclude era5land_API "$SRC/" "$DST/"
```

That is about 58 GB, most of it the Shuttle zips.

| Path | Content |
|---|---|
| `all_sites_multi_with_ru_flux_and_oneflux.csv` | Multi-network site list, the input of the registry |
| `insitu/PFT_consolidated.csv`, `PFT.csv`, `PFT_overrides.csv` | Plant functional types of the sites, with reviewed fixes |
| `insitu/site_coords_overrides.csv` | Reviewed coordinates of sites whose sources disagree |
| `insitu/*_consolidated.csv`, `FLUXNET_veg_frac.csv`, `forest_age.csv`, `globbiomass_agb_100m.csv` | In-situ biomass, soil depth, tree, vegetation fraction and forest age |
| `reference/ne_10m_admin_0_map_units.zip` | Natural Earth countries, to check that coordinates lie in the country of the site ID |
| `prefetch/fluxnet_shuttle/` | FLUXNET Shuttle zips and their index, the maintained eddy covariance source |
| `prefetch/jung_C2F/` | Quality flags of Jung et al. 2024, one file per site, with the code that made them |
| `prefetch/*_sites.csv`, `fcover300_2015/`, `modis_vcf_2010/`, `s2_*/`, `meta_chm_tiles.geojson` | Footprint means of the remote land cover, vegetation fraction, biomass and canopy height products |

The ERA5 caches of the Climate Data Store API are left out. No cube uses them, since
the ERA5 of the cubes comes through PyCLIFF from the ERA5 files of the cluster.

### ancillary_data/legacy/ (copied on the cluster)

`./setup_site_cubes.sh legacy` copies the LaThuile2007 site files, daily and hourly,
from `/Net/Groups/BGI/work_3/cliff/input/fluxnetBGI2021_siteCube-v0.03_iter-5` to
`legacy/fluxnetBGI2021_siteCube-v0.03_iter-5/`. It writes `SOURCE.md` with the origin
and date and `MD5SUMS`, then makes the copy read only. LaThuile2007 is not extracted
again, so this copy is its only source.

### ancillary_data/ (written by tools)

| Path | Written by |
|---|---|
| `insitu/site_lists/<source>.csv` | the workflow: from the v0.07 fluxcom site cube for its releases, from `legacy/` for LaThuile2007 |
| `insitu/sites_consolidated.csv` | the workflow, through `consolidate_sites.py`: the site registry, one row per site |
| `prefetch/clms_ndvi_300m/`, `clms_ndvi_1km/`, `clms_lai_300m/`, `clms_lai_1km/` | `prefetch_array.sh`: Copernicus NDVI and LAI, one csv per year |
| `prefetch/landsat_ndvi/` | `prefetch_array.sh`: Landsat NDVI, one csv per site |
| `prefetch/cliff_forcing/ERA5.e1/` | `prefetch_array.sh`: ERA5 at the nearest grid cells of each site, the forcing of PyCLIFF |
| `prefetch/gami_v2.1_1km_sites.csv` | the first extraction of the forest age map |

### input/ (written by the workflow)

| Path | Content |
|---|---|
| `store/<product>/<daily\|hourly>/<site>.nc` | One file per product and site, with `MANIFEST.csv` and the `STAGED` marker |
| `store/_logs/`, `store/_markers/` | Logs and markers of the workflow jobs |
| `pycliff/` | The PyCLIFF output root: gap-filled meteorology and fitted models per site |

### output/ (written by the workflow)

| Path | Content |
|---|---|
| `<run>.<source>.<DD\|HR>/data/<freq>/<site>.<y0>.<y1>.<freq>.nc` | Cube site files, one folder per source and resolution |
| `cliffin.<source>.HR/` | Hourly tower meteorology of each source, the input of PyCLIFF |
| `<source>_<freq>.zarr` | The site files of a cube combined into one site by time zarr store |

The sources are the FLUXNET Shuttle, LaThuile2007 and the releases of the v0.07
fluxcom site cube: FLUXNET2015, AmeriFlux, ICOS, OzFLUX, ICOSDrought2018,
WarmWinter2020 and CCBY.

## Data read from elsewhere

The workflow also reads, without copying them here, the global maps and series of
`/Net/Groups/BGI/work_3/sindbad/data`, DataStructureMDI and BIOMASCAT, the ERA5 files
of `/Net/Groups/data_BGC/era5/e1/0d25_hourly`, and the v0.07 fluxcom site cube.
