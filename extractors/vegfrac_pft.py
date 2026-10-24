#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Land cover fractions within 1 km of the site from the MODIS MCD12Q1 v061 IGBP
classification on the 1 km BGI grid, for every year from 2001 to 2024.

source_path is the file pattern with {year}. The share of each IGBP class is taken
over the classified cells whose centres lie within 1 km of the site. groups in the
json maps each component to its IGBP class numbers, and a component is the summed
share of its classes. With time_placement annual_extend, years outside the maps take
the nearest mapped year.

The PFT attribute of the output is the most frequent class over the years at the
cell holding the site, as the FLUXNET abbreviation, and fills sitePFT in the run
metadata. The source is local, so cache_path may be empty.
"""
import os
import sys

import numpy as np

import extractors.footprint_utils as fu

IGBP = fu.IGBP
YEARS = list(range(2001, 2025))
SOURCE_VARIABLE = "MCD12Q1_061_IGBP"


def years_of(info):
    return [int(y) for y in info.get("years", YEARS)]


def fetch_site(site, lat, lon, info):
    rows = []
    for year in years_of(info):
        path = info["source_path"].format(year=year)
        if not os.path.exists(path):
            continue
        data, inside, (r, c) = fu.nc_footprint(path, SOURCE_VARIABLE, lat, lon)
        classes = data[inside]
        valid = np.isin(classes, list(IGBP))
        row = {"site": site, "year": year, "site_class": int(data[r, c]),
               "n_valid": int(valid.sum()), "n_grids": int(inside.sum())}
        for k in IGBP:
            row[f"c{k}"] = float(np.mean(classes[valid] == k)) if valid.any() else np.nan
        rows.append(row)
    return rows


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers,
                      processes=True)


def components(rows, info):
    rows = rows.drop_duplicates("year", keep="last")
    out = {}
    for comp, classes in info["groups"].items():
        cols = [f"c{k}" for k in classes]
        out[comp] = dict(zip(rows["year"].astype(int),
                             rows[cols].sum(axis=1, min_count=len(cols)).astype(float)))
    return out


def site_attrs(rows, info):
    cls = rows["site_class"].astype(int)
    cls = cls[cls.isin(list(IGBP))]
    return {"PFT": IGBP[int(cls.mode().iloc[0])] if len(cls) else "undefined"}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
