#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Depth to bedrock within 1 km of the site from SoilGrids250m BDTICM_M_250m_ll.tif.

The value is the median of the 250 m pixels whose centres lie within 1 km of the
site, in cm as in the source file. The file is local, so cache_path is empty and the
window is read on every call. The SoilGrids map has no time dimension and is placed
under the year 2017 of its release.
"""
import sys

import numpy as np

import extractors.footprint_utils as fu

RELEASE_YEAR = 2017


def fetch_site(site, lat, lon, info):
    import rasterio
    from rasterio.windows import Window
    with rasterio.open(info["source_path"]) as src:
        col0, row0, w, h, inside = fu.geo_window(src.transform, src.width, src.height,
                                                 lat, lon)
        data = src.read(1, window=Window(col0, row0, w, h)).astype(float)
        data[data == src.nodata] = np.nan
    vals = data[inside]
    n_valid = int(np.isfinite(vals).sum())
    return [{"site": site, "median": np.nanmedian(vals) if n_valid else np.nan,
             "n_valid": n_valid, "n_pixels": int(inside.sum())}]


def components(rows, info):
    return {"depth": {RELEASE_YEAR: float(rows.iloc[-1]["median"])}}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
