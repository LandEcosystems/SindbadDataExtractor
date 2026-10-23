#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
import xarray as xr
import numpy as np

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)

def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []

    for tar_name in bxtr.vars_list:

        hilda_data_all = xr.open_dataset(bxtr.vars[tar_name]['data_path']).sel(time = slice(int(bxtr.start_date[0:4]),int(bxtr.end_date[0:4])))
        hilda_data = hilda_data_all.sel(latitude = bxtr.lat, longitude = bxtr.lon, method='nearest').LULC_states

        bxtr.log_var_start(tar_name)

        data = xr.Dataset({tar_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})

        if len(hilda_data) > 0:
            classes = np.atleast_1d(bxtr.vars[tar_name]['class'])
            # Pixels without a HILDA+ state stay NaN rather than counting as absent.
            vegfrac = hilda_data.isin(classes).astype(float).where(hilda_data.notnull())
            annual_vegfrac = []

            for year_ in hilda_data.time:

                date_y = np.arange(str(int(year_.values)) + '-01' + '-01', str(int(year_.values) +1) + '-01' + '-01'  ,dtype=date_.dtype)

                annual_vegfrac.append(xr.Dataset({tar_name: xr.DataArray(data = float(vegfrac.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_y})}))

            data = xr.merge(annual_vegfrac)
        else:
            logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")

        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

    src_dataset =  bxtr.merge_and_format(src_data)   
                  
    return src_dataset
    

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)