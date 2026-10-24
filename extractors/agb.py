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
    years = date_.astype('M8[Y]').astype(int) + 1970

    src_data = []

    for tar_name in bxtr.vars_list:

        bxtr.log_var_start(tar_name)

        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        # Some source files spell site IDs in other cases, e.g. Ca-Gro or IT-Sro.
        src_site = src_df.loc[src_df['Site_ID'].str.lower() == bxtr.site.lower()]
        src_site = src_site.drop(['Site_ID', 'lon', 'lat'], axis=1, errors='ignore')
        sel_data = np.array(src_site, dtype=float)

        data = xr.Dataset({
            tar_name:
            xr.DataArray(data=np.nan, dims=['time'], coords={'time': date_})
        })
        dval = data[tar_name].values
        if len(sel_data) > 0:
            # Each row is a year and its biomass, placed on the last time step of
            # that year.
            for yr, yr_bio in sel_data:
                in_year = np.flatnonzero(years == yr)
                if len(in_year) > 0 and yr_bio > 0:
                    dval[in_year[-1]] = yr_bio
            data[tar_name].values = dval
        else:
            logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")

        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

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