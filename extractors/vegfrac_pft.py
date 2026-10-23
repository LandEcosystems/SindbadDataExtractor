#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
import pandas as pd
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
        
        PFT_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        PFT_dat = PFT_df.loc[PFT_df['Site_ID'] == bxtr.site]
                
        bxtr.log_var_start(tar_name)

        data = xr.Dataset({tar_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})

        if len(PFT_dat) == 0:
            PFT_ = 'undefined'
            logger.warning(f"::MISSING:: variable {tar_name} has no PFT data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. Setting PFT as undefined")

        else:
            PFT_ = PFT_dat['PFT'].values[0]
            if PFT_ in bxtr.vars[tar_name]['PFT_types']:
                data[tar_name].values = np.ones_like(data[tar_name].values)
            else:
                data[tar_name].values = np.zeros_like(data[tar_name].values)

        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

    src_dataset =  bxtr.merge_and_format(src_data)   
        
    src_dataset = src_dataset.assign_attrs(PFT= PFT_)
    return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)