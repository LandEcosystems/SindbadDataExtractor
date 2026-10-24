#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Soil properties within 1 km of the site from SoilGrids v2.0 on the 250 m BGI grid.

source_path is the file pattern with {property} and {depth}, for example
clay_0-5cm_mean.172800.86400.nc. Each property and depth layer is summarised over
the grid cells whose centres lie within 1 km of the site, as the mean, or the median
when statistic is median in the json.

The BGI files keep the ISRIC integers with a scale_factor attribute that should
divide rather than multiply. The raw integers are therefore divided by the ISRIC
conversion factor, which gives bdod in kg dm-3, cec in cmol(c) kg-1, cfvo in vol %,
clay, sand and silt in %, nitrogen and soc in g kg-1, phh2o in pH, ocd in kg m-3 and
ocs in kg m-2.

The maps have no time dimension. A variable with several depths in the json has the
depth_soilGrids dimension, numbered from 1 in the order of its depths. A variable
with one depth, such as the 0-30 cm ocs stock, has none. The source is local, so
cache_path may be empty.
"""
import sys

import numpy as np
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

# ISRIC conversion factors from the stored integers to conventional units.
FACTOR = {"bdod": 100, "cec": 10, "cfvo": 10, "clay": 10, "nitrogen": 100,
          "phh2o": 10, "sand": 10, "silt": 10, "soc": 10, "ocd": 10, "ocs": 10}
MISSING = -32768


def layers(info):
    """
    Property and depth pairs of all variables, from the raw json as read by
    prefetch.py or from the merged variable info of run_extraction.py.
    """
    base = info.get("base_info", {}).get("depths", [])
    out = []
    for v in info["variables"].values():
        prop = v["sourceVariableName"]
        out += [(prop, d) for d in v.get("depths", base) if (prop, d) not in out]
    return out


def fetch_site(site, lat, lon, info):
    rows = []
    for prop, depth in layers(info):
        path = info["source_path"].format(property=prop, depth=depth)
        data, inside, _ = fu.nc_footprint(path, prop, lat, lon)
        vals = data[inside].astype(float)
        vals[vals == MISSING] = np.nan
        vals /= FACTOR[prop]
        n_valid = int(np.isfinite(vals).sum())
        rows.append({"site": site, "property": prop, "depth": depth,
                     "mean": np.nanmean(vals) if n_valid else np.nan,
                     "median": np.nanmedian(vals) if n_valid else np.nan,
                     "n_valid": n_valid, "n_grids": int(inside.sum())})
    return rows


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers,
                      processes=True)


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    info = config["dataset"][dataset]
    stat = info.get("statistic", "mean")
    rows = fu.rows_for_site(sys.modules[__name__], info, bxtr.site, bxtr.lat, bxtr.lon)

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        prop = bxtr.vars[tar_name]["sourceVariableName"]
        depths = bxtr.vars[tar_name]["depths"]
        vals = np.full(len(depths), np.nan)
        if len(rows):
            sel = rows[rows["property"] == prop].drop_duplicates("depth", keep="last")
            vals = sel.set_index("depth")[stat].reindex(depths).astype(float).values
        if not np.isfinite(vals).any():
            logger.warning(f"::MISSING:: variable {tar_name} has no data in "
                           f"{info['source_path']} for {bxtr.site}. NaN values will be "
                           "set.")
        coords = {'latitude': [bxtr.lat], 'longitude': [bxtr.lon]}
        if len(depths) > 1:
            data = xr.DataArray(vals.reshape(-1, 1, 1),
                                dims=['depth_soilGrids', 'latitude', 'longitude'],
                                coords=dict(depth_soilGrids=np.arange(len(depths)) + 1,
                                            **coords))
        else:
            data = xr.DataArray(vals.reshape(1, 1), dims=['latitude', 'longitude'],
                                coords=coords)
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
