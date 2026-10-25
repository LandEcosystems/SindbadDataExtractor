# FLUXNET Shuttle against FLUXNET2015

Compares what the `fluxnet_shuttle` extractor gives for each site with the
FLUXNET2015 SINDBAD site files, variable by variable, on the years and time steps
both products cover.

Report: https://claude.ai/artifact/VGJ887jE89dzd37x3V5mWZ

## Inputs

Under `~/research/RnD/data`, or the folder in `SINDBAD_VALIDATION_DATA`:

- `fn/<site>.1979.2017.daily.nc`: daily FLUXNET2015 site files
- `hourly/<site>.1979.2025.hourly.nc`: hourly FLUXNET2015 site files
- `sindbad_ancillary_data/prefetch/fluxnet_shuttle/`: shuttle zips and their index,
  from `prefetch.py ... shuttle_obs`

## Running

```
./run_all.sh <work folder> [python] [parallel jobs]
```

The default python is the `sinproc-prefetch` env. Without the `fluxcom` package, the
stub in `fluxcom_stub/` is used; the shuttle extractor never reads the FLUXCOM cube.
A run over 161 sites with 4 jobs takes about an hour.

## Steps

| Script | Output in the work folder |
|---|---|
| `compare_site.py <site> <out>` | `out/<site>_{daily,hourly}.csv`, metrics per variable; `series/daily/<site>.csv`, paired daily values; `series/diurnal/<site>.csv`, mean paired cycle over the day per calendar month |
| `aggregate.py <work>` | `all_sites.csv`, `summary_by_variable.csv`, `summary_by_site.csv`, `status.csv`, `hourly_lags.csv` |
| `era_check.py <site> <out>` | `era/<site>.csv`, ERA5 P, VPD, TA and SW_IN of the zip against measured |
| `make_dashboard_data.py <work>` | `dashboard/*.csv`, the files the report reads |

## Method

- The shuttle is extracted with the real extractor and the
  `configuration/{sindbad,hourly}/fluxnet_shuttle_{obs,clim}.json` configs, so
  units, daily statistics and hourly aggregation match the FLUXNET2015 files.
- The period is the overlap of the FLUXNET2015 record, the years with flux or met
  data, and the shuttle years. Only time steps where both have a value are scored.
- A QC variable is scored only where its flux exists in both products, since the
  FLUXNET2015 files set QC to 0 where the flux is missing.
- Per site and variable: paired count, means, bias, relative bias, MAE, RMSE,
  correlation, slope, Nash-Sutcliffe efficiency with FLUXNET2015 as reference,
  share of values equal within 1e-4 of the range, steps in only
  one product. Hourly fluxes, TA and SW_IN also get the shift of -3 to +3 hours with
  the highest correlation.
- Across sites: percentiles of the within-site correlation and the correlation of
  the site means of the two products.

## Dashboard

`dashboard/index.html` is the page of the report, a dashboard whose datasets are the
csv files of `make_dashboard_data.py` attached under the same names without `.csv`.
To refresh the report, upload the new files and point each dataset at its new file.
