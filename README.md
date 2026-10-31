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

The data live below `/Net/Groups/BGI/work_3/sindbad/site_cubes`:
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
- The forcing of the CLIFF gap filling, `cliff_forcing_era5`, is prefetched as well:
  `prefetch_array.sh cliff_forcing_era5` runs `pycliff-extract` once per ERA5 variable
  for all registry sites. Another forcing only needs another `cliff_forcing` json.

CLIFF then runs inside Snakemake: the product `pycliff_era5.<source>` runs PyCLIFF per
site on the hourly tower meteorology of each eddy covariance source, and the cubes get
the gap-filled meteorology as `<VAR>_ERA5_gfld`.

```
snakemake -n                          # what is missing or stale
snakemake --profile profiles/slurm    # run through SLURM
```

# Contact:

Sujan Koirala <skoirala@bgc-jena.mpg.de> (on the tool and data), Nuno Carvalhais <ncarvalhais@bgc-jena.mpg.de> (on the data)
