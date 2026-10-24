#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Forest age and the year of the last stand-replacing disturbance within 1 km of the
site from GAMI v2.1, Besnard et al. 2024, DOI:10.5880/GFZ.1.4.2023.006, at 100 m for
the map years 2010 and 2020.

All 20 ensemble members are read for the grid cells whose centres lie within 1 km of
the site. A grid cell is forest in a map year when every member gives it an age. The
age of a map year is the median over the forest grid cells of the ensemble median
age. Reading the compressed members takes about half a minute per site, so the cache
holds one row per site and map year, filled by prefetch.py or on demand with
fetch_missing.

The component of a variable names the map year after an underscore, and is one of
    age_<year>: forest age in that year, without a time dimension.
    disturbance_<year>: year of the disturbance, as the map year minus the rounded
        age, without a time dimension.
    intensity_<year>: 1 on the first time step of the disturbance year and 0 on every
        other time step, as the flags of disturbage_besnard2018.py.

An age at the GAMI cap of AGE_CAP years is old growth with no known disturbance, so
its disturbance year is NaN and its intensity is 0 throughout. Sites without any
GAMI forest in a map year get NaN in all variables of that year.
"""
import sys

import numpy as np
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

SOURCE_VARIABLE = "forest_age"
MAP_YEARS = [2010, 2020]
AGE_CAP = 300


def fetch_site(site, lat, lon, info):
    # Members, latitude, longitude and time, with fill values below 1.
    data, inside, _ = fu.nc_footprint(info["source_path"], SOURCE_VARIABLE, lat, lon)
    data = data.astype(float)
    data[data < 1] = np.nan
    vals = data[:, inside, :]
    rows = []
    for t, year in enumerate(MAP_YEARS):
        v = vals[..., t]
        forest = np.isfinite(v).all(axis=0)
        v = v[:, forest]
        rows.append({"site": site, "year": year, "n_grids": int(inside.sum()),
                     "n_valid": int(forest.sum()),
                     "age": float(np.median(np.median(v, axis=0))) if v.size else np.nan,
                     "age_std": float(v.std(axis=0).mean()) if v.size else np.nan})
    return rows


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers,
                      processes=True)


def ages(rows):
    """Median forest age of each map year, NaN where the site has no forest."""
    if "age" not in rows:
        return {}
    rows = rows.drop_duplicates("year", keep="last")
    return dict(zip(rows["year"].astype(int), rows["age"].astype(float)))


def disturbance_year(year, age):
    if not np.isfinite(age) or age >= AGE_CAP:
        return np.nan
    return float(year - np.round(age))


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()

    rows = fu.rows_for_site(sys.modules[__name__], info, bxtr.site, bxtr.lat, bxtr.lon)
    age_of = ages(rows) if len(rows) else {}
    if not np.isfinite(list(age_of.values())).any():
        logger.warning(f"::MISSING:: {dataset} has no data for {bxtr.site} in "
                       f"{info.get('cache_path', '')}. NaN values will be set. Run "
                       "prefetch.py or set fetch_missing to true.")

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        kind, year = bxtr.vars[tar_name]["component"].rsplit("_", 1)
        year = int(year)
        age = age_of.get(year, np.nan)
        dist = disturbance_year(year, age)
        if kind == "intensity":
            series = np.full(len(date_), np.nan if np.isnan(age) else 0.0)
            if np.isfinite(dist):
                series = fu.place(date_, {int(dist): 1.0}, "year_start")
                series = np.nan_to_num(series, nan=0.0)
            data = xr.DataArray(data=series, dims=['time'], coords={'time': date_})
        else:
            val = age if kind == "age" else dist
            data = xr.DataArray(np.full((1, 1), val), dims=['latitude', 'longitude'],
                                coords={'latitude': [bxtr.lat],
                                        'longitude': [bxtr.lon]})
        data = data.to_dataset(name=tar_name)
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    return fu.with_fetch_attrs(bxtr.merge_and_format(src_data), rows, info)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
