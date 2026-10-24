#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESA WorldCover class of the site at 10 m, from the 2020 v100 and 2021 v200 maps.

The maps are read through the Planetary Computer as in vegfrac_worldcover.py, so the
provider needs no cluster access. For each map year two classes are kept:
    majority: the most frequent class among the valid grid cells whose centres lie
        within 1 km of the site, with its share of those grid cells as
        majority_share.
    grid: the class of the grid cell whose centre is nearest to the site, which is
        the cell holding the site.

The classes are the WorldCover codes 10 tree cover, 20 shrubland, 30 grassland, 40
cropland, 50 built-up, 60 bare and sparse vegetation, 70 snow and ice, 80 permanent
water, 90 herbaceous wetland, 95 mangroves and 100 moss and lichen. The component of
a variable is majority or grid for the series over time, placed with
time_placement, or majority_<year> or grid_<year> for the class of one map year
without a time dimension. A _merged component fills the PFT of pft_path with the
class of WorldCover through BADM_PFT. Tree cover has no BADM PFT, so its PFT comes
from the datasets in rs_fill and then pft_path, as described in
footprint_utils.extract_classes.
"""
import sys

import numpy as np

import extractors.footprint_utils as fu
import extractors.vegfrac_worldcover as wc

YEARS = [2020, 2021]
WC_NAMES = {10: "tree cover", 20: "shrubland", 30: "grassland", 40: "cropland",
            50: "built-up", 60: "bare and sparse vegetation", 70: "snow and ice",
            80: "permanent water bodies", 90: "herbaceous wetland", 95: "mangroves",
            100: "moss and lichen"}
# BADM PFT of each WorldCover class. Tree cover has no leaf type, so no PFT.
BADM_PFT = {95: "EBF", 20: "OSH", 30: "GRA", 100: "GRA", 40: "CRO", 50: "URB",
            60: "BSV", 70: "SNO", 80: "WAT", 90: "WET"}


def years_of(info):
    return [int(y) for y in info.get("years", YEARS)]


def fetch_site(site, lat, lon, info):
    rows = []
    for year in years_of(info):
        codes, plats, plons, n_items = wc.read_grids(lat, lon, year)
        row = {"site": site, "year": year, "n_items": n_items, "n_grids": 0,
               "majority": np.nan, "majority_share": np.nan, "grid": np.nan}
        if codes is not None:
            dist = fu.haversine(lat, lon, plats, plons)
            c = codes[dist <= fu.RADIUS_M]
            row["grid"] = int(codes[np.argmin(dist)])
            row["n_grids"] = int(c.size)
            if c.size:
                vals, counts = np.unique(c, return_counts=True)
                row["majority"] = int(vals[np.argmax(counts)])
                row["majority_share"] = float(counts.max() / c.size)
        rows.append(row)
    return rows


LEGEND = fu.class_legend(WC_NAMES, BADM_PFT)


def to_badm(code):
    """BADM PFT of a WorldCover class, or None."""
    return BADM_PFT.get(int(code))


def badm_series(rows):
    return fu.class_badm_series(rows, to_badm)


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def extract(dataset, site_info, config):
    return fu.extract_classes(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
