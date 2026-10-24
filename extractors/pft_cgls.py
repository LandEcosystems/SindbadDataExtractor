#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copernicus Global Land Service LC100 collection 3 class of the site at 100 m, from
the discrete classification maps of 2015 to 2019.

The maps are read remotely from the Zenodo GeoTIFFs as in vegfrac_cgls.py, so the
provider needs no cluster access. For each map year two classes are kept:
    majority: the most frequent class among the classified grid cells whose centres
        lie within 1 km of the site, with its share of those grid cells as
        majority_share.
    grid: the class of the grid cell holding the site.

The classes are the LC100 codes 20 shrubs, 30 herbaceous vegetation, 40 cropland, 50
urban, 60 bare and sparse vegetation, 70 snow and ice, 80 permanent water, 90
herbaceous wetland, 100 moss and lichen, 200 open sea, and the forests 111 to 116
closed and 121 to 126 open, each as evergreen needleleaf, evergreen broadleaf,
deciduous needleleaf, deciduous broadleaf, mixed and unknown. The variables are set
up as in footprint_utils.extract_classes, including the _merged components that
fill the PFT of pft_path with the class of LC100.

The PFTs of the classes are given by BADM_PFT. Forest of unknown leaf type has none
and is filled as described in footprint_utils.extract_classes, from rs_fill and
pft_path. The output attributes are set there too. PFT fills sitePFT in the run
metadata when pft_cgls comes first in the main json.
"""
import sys

import numpy as np

import extractors.footprint_utils as fu
from extractors.vegfrac_cgls import EPOCHS, GDAL_ENV

URL = ("/vsicurl/https://zenodo.org/records/{record}/files/"
       "PROBAV_LC100_global_v3.0.1_{epoch}_Discrete-Classification-map_EPSG-4326.tif")
# Codes 0 and 255 are unclassified and no data.
UNCLASSIFIED = (0, 255)
LC100_NAMES = {20: "shrubs", 30: "herbaceous vegetation", 40: "cropland",
               50: "urban and built up", 60: "bare and sparse vegetation",
               70: "snow and ice", 80: "permanent water bodies",
               90: "herbaceous wetland", 100: "moss and lichen",
               111: "closed forest evergreen needleleaf",
               112: "closed forest evergreen broadleaf",
               113: "closed forest deciduous needleleaf",
               114: "closed forest deciduous broadleaf", 115: "closed forest mixed",
               116: "closed forest unknown", 121: "open forest evergreen needleleaf",
               122: "open forest evergreen broadleaf",
               123: "open forest deciduous needleleaf",
               124: "open forest deciduous broadleaf", 125: "open forest mixed",
               126: "open forest unknown", 200: "open sea"}
# BADM PFT of each LC100 class. Closed and open forest share the leaf type, and LC100
# does not separate open from closed shrubs. Forest of unknown type has no PFT.
BADM_PFT = {111: "ENF", 121: "ENF", 112: "EBF", 122: "EBF", 113: "DNF", 123: "DNF",
            114: "DBF", 124: "DBF", 115: "MF", 125: "MF", 20: "OSH", 30: "GRA",
            40: "CRO", 50: "URB", 60: "BSV", 70: "SNO", 80: "WAT", 90: "WET",
            100: "GRA", 200: "WAT"}


def years_of(info):
    return [int(y) for y in info.get("years", list(EPOCHS))]


def fetch_site(site, lat, lon, info):
    import rasterio
    from rasterio.windows import Window
    rows = []
    with rasterio.Env(**GDAL_ENV):
        for year in years_of(info):
            record, epoch = EPOCHS[year]
            with rasterio.open(URL.format(record=record, epoch=epoch)) as src:
                col0, row0, w, h, inside = fu.geo_window(src.transform, src.width,
                                                         src.height, lat, lon)
                data = src.read(1, window=Window(col0, row0, w, h))
                r, c = src.index(lon, lat)
            row = {"site": site, "year": year, "n_grids": int(inside.sum()),
                   "majority": np.nan, "majority_share": np.nan,
                   "grid": float(data[r - row0, c - col0])}
            if row["grid"] in UNCLASSIFIED:
                row["grid"] = np.nan
            vals = data[inside]
            vals = vals[~np.isin(vals, UNCLASSIFIED)]
            if vals.size:
                codes, counts = np.unique(vals, return_counts=True)
                row["majority"] = int(codes[np.argmax(counts)])
                row["majority_share"] = float(counts.max() / vals.size)
            rows.append(row)
    return rows


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    # Zenodo throttles clients, so few workers.
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons,
                      workers or info.get("n_workers", 2))


LEGEND = fu.class_legend(LC100_NAMES, BADM_PFT)


def to_badm(code):
    """BADM PFT of an LC100 class, or None."""
    return BADM_PFT.get(int(code))


def badm_series(rows):
    return fu.class_badm_series(rows, to_badm)


def extract(dataset, site_info, config):
    return fu.extract_classes(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
