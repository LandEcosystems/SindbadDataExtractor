#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
import xarray as xr
import numpy as np

import logging
logger = logging.getLogger(__name__)
from extractors.BasexTractor import BasexTractor


def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    src_data = []
    for tar_name in bxtr.vars_list:
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        bxtr.log_var_start(tar_name)

        data = xr.open_dataset(bxtr.vars[tar_name]['data_path'])

        # Handle coordinate naming variations (lat/lon vs latitude/longitude)
        # For point data or malformed coordinates, use integer indexing (isel) instead of labeled indexing (sel)
        lat_coord = 'latitude' if 'latitude' in data.dims else ('lat' if 'lat' in data.dims else None)
        lon_coord = 'longitude' if 'longitude' in data.dims else ('lon' if 'lon' in data.dims else None)
        
        # Only use positional indexing if dimensions exist (data with spatial grid)
        if lat_coord and lon_coord:
            # Use isel to avoid coordinate index issues
            data = data.isel({lat_coord: 0, lon_coord: 0})
        
        #%% interpolate missing data
        data = data.where(np.isfinite(data)).interpolate_na(
            dim='time', method='nearest').compute()

        data = data.rename(**{src_name: tar_name})

        data = bxtr.convert_units(data, tar_name)

        src_data.append(data)
        bxtr.log_var_end(data[tar_name], tar_name, None)

    src_dataset = bxtr.merge_and_format(src_data)
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Provider at: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
