#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 22 2021

@author: @dr-ko
"""
#%% Load libray
import xarray as xr
import numpy as np
import json
import logging
import os
logger = logging.getLogger(__name__)


#%% Process data
def dataset_json_path(config_dir, name, ancillary_dir=""):
    """
    Path of the dataset json name, from config_dir when it is there, otherwise from
    ancillary_dir. The jsons that are the same for every temporal resolution live in
    ancillary_dir, which defaults to the ancillary folder next to config_dir.
    """
    if not ancillary_dir:
        ancillary_dir = os.path.join(os.path.dirname(os.path.normpath(config_dir)),
                                     "ancillary")
    for folder in (config_dir, ancillary_dir):
        path = os.path.join(folder, name)
        if os.path.exists(path):
            return path
    return os.path.join(config_dir, name)


def data_structure_temporal_NoDepth(data, lat, lon):

    var_ = list(data.keys())
    if 'time' in data.dims and 'depth_FLUXNET' not in data.dims:
        ds_ = xr.DataArray(data.to_array().values.reshape(len(var_), -1, 1, 1),
                           dims=['variable', 'time', 'latitude', 'longitude'],
                           coords={
                               'variable': var_,
                               'time': data.time.values,
                               'latitude': [lat],
                               'longitude': [lon]
                           }).to_dataset(dim="variable")
    else:
        ds_ = data
    return ds_

def flatten_xr_hour_to_time(x_data):
    if 'hour' in x_data.dims:
        import pandas as pd
        time1 = x_data['time'].values[0]
        x_data = x_data.rename({
            'time': 'date'
        }).stack(time=('date', 'hour')).reset_index('time')
        times = pd.date_range(time1, periods=len(x_data['time']), freq='H')
        x_data['time'] = times
    else:
        x_data = x_data
    return x_data

def flatten_hour_to_time(_data):
    if isinstance(_data, list):
        _data = [flatten_xr_hour_to_time(__data) for __data in _data]
    else:
        _data = flatten_xr_hour_to_time(_data)
    return _data


def log_and_print(task,
                  site,
                  prod_name,
                  var_name,
                  src_name,
                  temporal_resolution,
                  config,
                  tar_label='target',
                  src_label='source'):
    logger.info(
        f"{config['FLUXNET_version']}:::{task}::{temporal_resolution}::{site} -> {tar_label}: {var_name.ljust(30,' ')}, {src_label}: {src_name.ljust(20,' ')} | {prod_name}"
    )
    print(
        f"{config['FLUXNET_version']}:::{task}::{temporal_resolution}::{site} -> {tar_label}: {var_name.ljust(30,' ')}, {src_label}: {src_name.ljust(20,' ')} | {prod_name}"
    )
    return


def log_and_print_sep(separator='-' * 160):
    logger.info(separator)
    print(separator)
    return


def log_datavar_info(dataset,
                     site,
                     t_reso,
                     src_var_name,
                     tar_var_name,
                     json_units_tar,
                     json_units_src,
                     bounds_in,
                     data,
                     transform,
                     partitioning='None'):
    if isinstance(data, xr.core.dataset.Dataset):
        data_ar = data[tar_var_name]
    else:
        data_ar = data

    data_units = data_ar.units
    tname = get_transform_name(transform)

    ## get data information
    data_a = data_ar.values
    data_min = np.nanmin(data_a)
    data_max = np.nanmax(data_a)
    ndata = data_a.size
    nan_c = np.isnan(data_a).sum()
    nan_c_p = round(nan_c * 100 / (ndata), 2)
    if nan_c_p > 90:
        logger.warning(f'nanwarn: {nan_c_p} % of data points are nan.')
    data_tmp = data_a[~np.isnan(data_a)]
    n_data_tmp = data_tmp.size
    n_out_bounds = ((bounds_in[0] > data_tmp) &
                    (data_tmp > bounds_in[1])).sum()
    n_out_bounds_perc = round(n_out_bounds * 100 / (n_data_tmp), 2)

    logger.info(f'\n Dataset:: {dataset}\
                    \n site: {site}\
                    \n reso:: {t_reso}\
                    \n var::\
                        \n \t src: {src_var_name}\
                        \n \t target: {tar_var_name}\
                        \n \t partitioning: {partitioning}\
                    \n units::\
                        \n \t data: {data_units}\
                        \n \t json-target:: {json_units_tar}\
                        \n \t json-source: {json_units_src}\
                    \n data (after processing)::\
                        \n \t shape: {data_a.shape}\
                        \n \t size: {data_a.size}\
                        \n \t min: {data_min}\
                        \n \t max: {data_max}\
                    \n json_bnds::\
                        \n \t min: {bounds_in[0]}\
                        \n \t max: {bounds_in[1]}\
                    \n counts::\
                        \n \t out_bounds_none_NaN: {n_out_bounds} ({n_out_bounds_perc}%)\
                        \n \t nans: {nan_c} ({nan_c_p}%)\
                    \n transform:: {tname}'
                )
    return


def save_json(json_path, data_to_save):
    with open(json_path, 'w') as fp:
        json.dump(data_to_save, fp, indent=2, sort_keys=True)
    return


def do_unit_conversion(data, var_name, src_units, tar_units, units_scalar):
    data_out = data.copy()
    if units_scalar == 1.:
        if src_units.lower() != tar_units.lower():
            logger.warning(
                f'Unit Conversion::\
                    \n {var_name}:\
                    \n units\
                        \n \t src: [{src_units}]\
                        \n \t target: [{tar_units}]\
                        \n \t scaled by: {units_scalar}\
                    \n Note: potentially inconsistent units. Make sure that transforms below take care of unit conversion.'
            )
        else:
            logger.warning(
                f'Unit Conversion::\
                    \n {var_name}:\
                    \n units\
                        \n \t src: [{src_units}]\
                        \n \t target: [{tar_units}]\
                        \n \t scaled by: {units_scalar}\
                    \n Note :: No unit conversion needed.'
            )

            data_out.attrs["units"] = src_units
    else:
        data_out.values = data_out.values * units_scalar
        logger.info(f'Unit Conversion::\
                    \n {var_name}:\
                    \n units\
                        \n \t src: [{src_units}]\
                        \n \t target: [{tar_units}]\
                        \n \t scaled by: {units_scalar}\
                    \n mean:\
                        \n \t before: {np.nanmean(data.values)}\
                        \n\t after: {np.nanmean(data_out.values)}'
        )
        data_out.attrs["units"] = tar_units
    return data_out


def get_transform_name(_transform):
    if _transform is None or _transform is False or len(_transform) == 0:
        tname = "None"
    else:
        tname = ",".join([__tr.class_name for __tr in _transform])
    return tname


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Utilities: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('shared functions: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
