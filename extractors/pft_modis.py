#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IGBP class of the site from the MODIS MCD12Q1 v061 classification on the 1 km BGI
grid, for every year from 2001 to 2024.

The class is the IGBP number from 1 to 17 of the cell holding the site. Its legend is
1 ENF, 2 EBF, 3 DNF, 4 DBF, 5 MF, 6 CSH, 7 OSH, 8 WSA, 9 SAV, 10 GRA, 11 WET, 12 CRO,
13 URB, 14 CVM, 15 SNO, 16 BSV and 17 WAT. The components are
    class: the class of each year, placed with time_placement, so that with
        annual_extend the years outside the maps take the nearest mapped year.
    class_<year>: the class of that one year, without a time dimension. The year is
        any of the mapped years, for example class_2010.
    class_merged: the PFT of the site in the csv at pft_path on every time step, or
        the class of each year where the csv has no IGBP class for the site.

The cells are read and cached as in vegfrac_pft.py, so both datasets may share one
cache_path. The PFT attribute of the output is the most frequent class over the
years as the BADM PFT, and fills sitePFT in the run metadata when no dataset before
it has one. A site without MODIS data takes its PFT from the datasets in rs_fill and
then from pft_path, through footprint_utils.resolve_pft, and PFT_source names the
one that gave it.
"""
import sys

import numpy as np
import xarray as xr

import extractors.footprint_utils as fu
import extractors.vegfrac_pft as vf

import logging
logger = logging.getLogger(__name__)

fetch_site = vf.fetch_site


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers,
                      processes=True)


def classes(rows):
    """IGBP class of the site cell in each year, without unclassified years."""
    rows = rows.drop_duplicates("year", keep="last")
    cls = dict(zip(rows["year"].astype(int), rows["site_class"].astype(int)))
    return {y: float(c) for y, c in cls.items() if c in vf.IGBP}


def badm_series(rows):
    return fu.badm_by_year(classes(rows), vf.IGBP.get)


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    pft_path = info.get("pft_path", "")

    rows = fu.rows_for_site(sys.modules[__name__], info, bxtr.site, bxtr.lat, bxtr.lon)
    cls = classes(rows) if len(rows) else {}
    if not cls:
        logger.warning(f"::MISSING:: {dataset} has no data for {bxtr.site} in "
                       f"{info.get('cache_path', '') or info['source_path']}. NaN "
                       "values will be set.")
    years = sorted(cls) or [int(y) for y in info.get("years", [])]
    chain = [(dataset, fu.badm_by_year(cls, vf.IGBP.get))]
    if not chain[0][1]:
        chain += fu.fill_sources(config, info.get("rs_fill", []), bxtr.site,
                                 bxtr.lat, bxtr.lon)

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        comp = bxtr.vars[tar_name]["component"]
        if comp == "class_merged":
            series = fu.merge_pft(date_, pft_path, bxtr.site, years, chain,
                                  info.get("time_placement", "annual_extend"))
            data = xr.DataArray(data=series, dims=['time'], coords={'time': date_})
        elif comp == "class":
            series = fu.place(date_, cls, info.get("time_placement", "annual_extend"))
            data = xr.DataArray(data=series, dims=['time'], coords={'time': date_})
        else:
            val = cls.get(int(comp.split("_")[1]), np.nan)
            data = xr.DataArray(np.full((1, 1), val), dims=['latitude', 'longitude'],
                                coords={'latitude': [bxtr.lat],
                                        'longitude': [bxtr.lon]})
        data = data.to_dataset(name=tar_name)
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    pft, source = fu.pft_mode(years, chain, pft_path, bxtr.site)
    return src_dataset.assign_attrs(PFT=pft, PFT_source=source)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
