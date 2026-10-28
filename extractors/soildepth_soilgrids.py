#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Depth to bedrock within 1 km of the site from SoilGrids250m BDTICM_M_250m_ll.tif.

The value is the median of the 250 m grid cells whose centres lie within 1 km of the
site, in cm as in the source file. The file is local, so cache_path is empty and the
window is read on every call. The SoilGrids map has no time dimension and is placed
under the year 2017 of its release.

A variable with merge_insitu takes the depth in m of the site from the csv at
insitu_path, column insitu_column, and otherwise the SoilGrids depth clipped to
clip_m, so that depths at or below 0 become 0.01 m and depths at or above 50 m
become 49.99 m. The in-situ depth is not clipped.
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
             "n_valid": n_valid, "n_grids": int(inside.sum())}]


def components(rows, info):
    return {"depth": {RELEASE_YEAR: float(rows.iloc[-1]["median"])}}


def merge_insitu(comps, info, site, targets, placement, full):
    """
    In-situ depth where the csv has one, else the clipped SoilGrids depth, in the cm
    of the source. full is the number of cm in a m.
    """
    ins = fu.site_value(info.get("insitu_path", ""), site,
                        info.get("insitu_column", "soilDepth"))
    lo, hi = (v * full for v in info.get("clip_m", [0.01, 49.99]))
    out = {}
    for var, (comp, _) in targets.items():
        series = fu.collapse(comps.get(comp, {}), placement)
        if np.isfinite(ins):
            out[var] = {k: ins * full for k in series}
        else:
            out[var] = {k: float(np.clip(v, lo, hi)) if np.isfinite(v) else np.nan
                        for k, v in series.items()}
    return out


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
