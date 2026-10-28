#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
import pandas as pd
import xarray as xr
import numpy as np
# import ipdb

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []

    src_data_dict = {}

    for tar_name in bxtr.vars_list:
        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        if '_PFT' in tar_name:
            PFT_dat = src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
                ['Site_ID'], axis=1)
            if len(PFT_dat) == 0:
                sel_data = 'undefined'
            else:
                sel_data = PFT_dat['PFT'].values[0]
            print(sel_data)
        elif '_globBiomass' in tar_name:
            sel_data = np.array(src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
                ['Site_ID'], axis=1),
                                dtype=float)
            data = xr.Dataset({
                tar_name:
                xr.DataArray(data=np.nan, dims=['time'], coords={'time': date_})
            })

            sel_data = sel_data[2:] # ignore the latitude and longitude co-ordinates
            if len(sel_data) > 0:
                sel_data[sel_data <= 0] = np.nan
                if '_std' in tar_name:
                    sel_data = np.nanstd(sel_data)
                else:
                    sel_data = np.nanmean(sel_data)
                dval = data[tar_name].values
                for dI, _date in enumerate(date_):
                    if _date.astype(object).year == 2011:
                        dval[dI-1] = sel_data
                        break

                data[tar_name].values = dval
            else:
                logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")            
            sel_data = bxtr.convert_units(data, tar_name)
        else:
            sel_data = np.array(src_df.loc[src_df['Site_ID'] == bxtr.site].drop(
                ['Site_ID'], axis=1),
                                dtype=float)
            data = xr.Dataset({
                tar_name:
                xr.DataArray(data=np.nan, dims=['time'], coords={'time': date_})
            })

            if len(sel_data) > 0:
                sel_data[sel_data <= 0 ] = np.nan
                dval = data[tar_name].values
                for _sel_data in sel_data:
                    yr = _sel_data[0]
                    yr_bio = _sel_data[1]
                    for dI, _date in enumerate(date_):
                        if _date.astype(object).year == yr + 1:
                            dval[dI-1] = yr_bio
                            break
                # sel_data = np.nanmedian(sel_data)
                data[tar_name].values = dval
                # data[tar_name].values = np.ones_like(
                #     data[tar_name].values) * sel_data
            else:
                logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")
            sel_data = bxtr.convert_units(data, tar_name)
        src_data_dict[tar_name] = sel_data
        # ipdb.set_trace()


    # merge the data manually (rough method but works.. be careful of hardcoded variable names. things may break if there are changes in extractor's json settings)

    data = xr.where(np.isfinite(src_data_dict['agb_merged']['agb_merged']), src_data_dict['agb_merged']['agb_merged'], src_data_dict['agb_globBiomass']['agb_globBiomass'])
    tar_name = 'agb_merged'
    # needs renaming and resetting the attributes because xarray likes deleting things
    bxtr.log_var_start(tar_name)
    data = data.rename(tar_name)
    data = bxtr.convert_units(data, tar_name)

    agb_insitu = src_data_dict['agb_insitu']['agb_insitu']
    agb_insitu_values = agb_insitu.values
    has_insitu_data = False
    if sum(np.isnan(agb_insitu_values)) != len(agb_insitu_values):
        has_insitu_data = True
        logger.warning(f"insitu AGB data value/s available for {bxtr.site}. Using only the insitu data and ignoring glob-biomass data.")
        data.values = agb_insitu_values
    # tar_name = 'agb_merged'
    # # needs renaming and resetting the attributes because xarray likes deleting things
    # bxtr.log_var_start(tar_name)
    # data = data.rename(tar_name)
    # data = bxtr.convert_units(data, tar_name)
    # src_data_dict[tar_name][tar_name] = data.copy()

    src_data_dict[tar_name][tar_name] = data.copy()
    # ipdb.set_trace()


    tar_name = 'agb_merged_PFT'
    # needs renaming and resetting the attributes because xarray likes deleting things
    bxtr.log_var_start(tar_name)
    data = data.rename(tar_name)

    pft = src_data_dict[tar_name]
    if pft in bxtr.vars[tar_name]['PFT_types']:
        logger.warning(f"PFT ({pft}) in the selected PFTs ({bxtr.vars[tar_name]['PFT_types']}) for fixed biomass for {bxtr.site}. Using only fixed AGB of {bxtr.vars[tar_name]['AGB_scalar']} for agb_merged_PFT.")
        if has_insitu_data:
            dval = (agb_insitu_values / agb_insitu_values) * bxtr.vars[tar_name]['AGB_scalar']
        else:
            dval = data.values
            for dI, _date in enumerate(date_):
                if _date.astype(object).year == 2011:
                    dval[dI-1] = bxtr.vars[tar_name]['AGB_scalar']
                    break
        data.values=dval 
        # data.values = data.values * bxtr.vars[tar_name]['AGB_scalar']
        # if np.sum(np.isnan(data.values.flatten())) == len(data.values.flatten()):
            # data.values=np.ones_like(data.values) * bxtr.vars[tar_name]['AGB_scalar']
        # data.values = np.ones_like(data.values) * bxtr.vars[tar_name]['AGB_scalar']

    src_data_dict[tar_name] = data

    # only include the merged variables in this extractor to avoid confusion with agb.py. Otherwise, one could just do src_data_dict.values() and include everything
    src_dataset = bxtr.merge_and_format([src_data_dict['agb_merged_PFT'], src_data_dict['agb_merged']])
    for v in ('agb_merged_PFT', 'agb_merged'):
        bxtr.warn_merged_nan(src_data_dict[v], v)
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