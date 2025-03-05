#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
#%% Load library
import xarray as xr
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, date
import importlib
import logging
sys.path.append(os.path.join(os.getcwd(), '../'))
import copy
from extractors.cliff_gapfill import cliff_gapfill
import collections
import utils.shared_utils as shut


def get_site_info(site):
    import fluxcom as flx
    lat_lon_utc = flx.site_cube.site_cube_utils.get_site_lat_lon_utc(site=site)
    lat=lat_lon_utc[0]
    lon=lat_lon_utc[1]
    site_info = {"site_ID": site, "latitude": lat, "longitude": lon}
    return site_info


def load_extractor(mod_name):
    mod = importlib.import_module("extractors." + mod_name)
    return mod


def add_date_variables(dat):
    julian_base_date = date(1592, 10, 15)
    julian_base_date_str = julian_base_date.strftime("%Y-%m-%d")
    time = dat.time
    dates = pd.DatetimeIndex(time.values)
    year = dates.year.to_numpy()
    month = dates.month.to_numpy()
    day = dates.day.to_numpy()
    hour = dates.hour.to_numpy()

    julian_day = np.zeros_like(year)
    for _ind in range(len(year)):
        data_date = date(year[_ind], month[_ind], day[_ind])
        delta = data_date - julian_base_date
        julian_day[_ind] = delta.days

    x_year = xr.Variable(["time"],
                         year,
                         attrs={
                             'long_name': 'year',
                             '_FillValue': -9999
                         })
    dat = dat.assign(year=x_year)
    x_month = xr.Variable(["time"],
                          month,
                          attrs={
                              'long_name': 'month of year',
                              '_FillValue': -9999
                          })
    dat = dat.assign(month=x_month)
    x_day = xr.Variable(["time"],
                        day,
                        attrs={
                            'long_name': 'day of year',
                            '_FillValue': -9999
                        })
    dat = dat.assign(day=x_day)
    x_hour = xr.Variable(["time"],
                         hour,
                         attrs={
                             'long_name': 'hour of day',
                             '_FillValue': -9999
                         })
    dat = dat.assign(hour=x_hour)
    x_julian_day = xr.Variable(["time"],
                               julian_day,
                               attrs={
                                   'long_name':
                                   f'days since {julian_base_date_str}',
                                   'FillValue': np.nan
                               })
    dat = dat.assign(julian_day=x_julian_day)
    return dat


def get_bounds_str(data, bounds_in):
    bounds_str = f'bounds: data: [{round(np.nanmin(data.values),3)}, {round(np.nanmax(data.values),3)}], json: [{bounds_in[0]}, {bounds_in[1]}]'
    return bounds_str


def plot_dianostic_figures(ds, site_info, exp_settings, site, resampled=None):
    import matplotlib.pyplot as plt
    if resampled is None:
        data_freq = exp_settings["temporal_resolution"]
    else:
        data_freq = resampled
    site_fig_dir = os.path.join(exp_settings['output_dir_path']['figs'], data_freq)

    os.makedirs(site_fig_dir, exist_ok=True)
    skipvars = 'month day year hour julian_day'.split()
    plotvars = list(set(list(ds.keys())) - set(skipvars))
    fig_size = (12, 2.5)
    for data_key in exp_settings["sel_datasets"]:
        data_set = exp_settings["dataset"][data_key]
        for var_, var_info in data_set['variables'].items():
            if var_ in plotvars:
                dat_ = ds[var_]
                shut.log_and_print('plot',
                                   site,
                                   dat_.sourceDataProductName,
                                   var_,
                                   dat_.units,
                                   data_freq,
                                   exp_settings,
                                   tar_label='variable',
                                   src_label='units')

                if ('depth_FLUXNET' in dat_.dims) & ('time' in dat_.dims):
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=fig_size,
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    for layer_ in dat_.depth_FLUXNET.values:
                        ax.plot(
                            dat_.sel(
                                depth_FLUXNET=layer_).time.values.reshape(-1),
                            dat_.sel(depth_FLUXNET=layer_).values.reshape(-1),
                            lw=0.6,
                            label='Layer ' + str(layer_))
                    plt.legend(loc='best', fontsize=5)
                    x_l = ax.set_xlabel('time', size=8)
                elif ('time' in dat_.dims) & (
                        'depth_FLUXNET' not in dat_.dims) & ('depth_soilGrids'
                                                             not in dat_.dims):
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=fig_size,
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    if not var_.endswith("_gfld"):
                        ax.plot(dat_.time.values.reshape(-1),
                                dat_.values.reshape(-1),
                                lw=0.6,
                                ls=":")
                    else:

                        ax.plot(dat_.time.values.reshape(-1),
                                dat_.values.reshape(-1),
                                lw=0.6,
                                ls=":",
                                label=dat_.sourceDataProductName)
                        dat_ori_ = ds[var_info['with_gaps']]

                        ax.plot(dat_ori_.time.values.reshape(-1),
                                dat_ori_.values.reshape(-1),
                                lw=0.6,
                                ls="--",
                                label=dat_ori_.sourceDataProductName)
                        plt.legend(loc='best', fontsize=5)

                    x_l = ax.set_xlabel('time', size=8)
                elif 'depth_soilGrids' in dat_.dims:
                    if 'time' in dat_.dims:
                        ## time axis is automatically generated in the monthly data that are created using resampling
                        dat_ = dat_.mean(dim='time', keep_attrs=True)
                    fig, ax = plt.subplots(1,
                                           1,
                                           figsize=(4, 4),
                                           gridspec_kw={
                                               'wspace': 0.35,
                                               'hspace': 0.5
                                           })
                    ax.scatter(dat_.depth_soilGrids.values, dat_.values, s=10)
                    x_l = ax.set_xlabel('depth_soilGrids', size=8)
                ax.tick_params(labelsize=6)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                y_l = ax.set_ylabel(data_freq + " " + var_ + '\n[' +
                                    dat_.units + ']',
                                    size=8)
                t_l = plt.title(
                    f'{exp_settings["FLUXNET_version"]}:::{site_info["site_ID"]}:: {data_key} [exp-json]- src: {dat_.sourceDataProductName} \n{get_bounds_str(dat_,var_info["bounds"])}',
                    multialignment='center',
                    va='bottom',
                    size=8,
                    color='#339900',
                    y=1.06)
                plt.savefig(
                    site_fig_dir + "/" + site_info["site_ID"] + '.' + var_.lower() +
                    "." + exp_settings['start_date'].split('-')[0] + "-" +
                    str(int(exp_settings['end_date'].split('-')[0])) + '.' +
                    data_freq + '.png',
                    bbox_inches='tight',
                    # box_extra_artists=[x_l, y_l, t_l],
                    dpi=300)
                plt.close()
            else:
                logging.warning(
                    f'Variable {var_} from {data_key} dataset was not found in processed data.'
                )
    return


def add_gapfill_settings(exp_settings,
                         fill_src=None,
                         fill_tar=None,
                         dataset_key=None,
                         qc_thres=None):
    _exp_settings = copy.deepcopy(exp_settings)

    _exp_settings['sel_datasets'] = np.append(_exp_settings['sel_datasets'],
                                              dataset_key)
    _exp_settings['dataset'][dataset_key] = copy.deepcopy(
        _exp_settings['dataset'][fill_tar])
    _exp_settings['dataset'][dataset_key]['gfld_src'] = fill_src
    _exp_settings['dataset'][dataset_key]['gfld_tar'] = fill_tar

    src_suffix = _exp_settings['dataset'][fill_src]['var_suffix']
    if len(src_suffix.strip()) > 0:
        src_suffix = '_' + src_suffix

    tar_suffix = _exp_settings['dataset'][fill_tar]['var_suffix']
    if len(tar_suffix.strip()) > 0:
        tar_suffix = '_' + tar_suffix
    tar_info = copy.deepcopy(_exp_settings['dataset'][fill_tar]['variables'])
    src_info = copy.deepcopy(_exp_settings['dataset'][fill_src]['variables'])
    _v_tars = list(_exp_settings['dataset'][fill_tar]['variables'].keys())

    for _v_tar in _v_tars:
        _v_tar_gfld = _v_tar + src_suffix + '_gfld'
        _v_src = _v_tar_gfld.replace(tar_suffix, '').replace('_gfld', '')
        gfl_var_info = copy.deepcopy(tar_info[_v_tar])
        if _v_src in src_info.keys():
            gfl_var_info[
                'sourceDataProductName'] = f"{tar_info[_v_tar]['sourceDataProductName']} gapfld-w {src_info[_v_src]['sourceDataProductName']}"
            gfl_var_info[
                'publication'] = f"{tar_info[_v_tar]['publication']}; {src_info[_v_src]['publication']}"
            gfl_var_info[
                'description'] = f"{tar_info[_v_tar]['sourceDataProductName']}:: {tar_info[_v_tar]['description']};{src_info[_v_src]['sourceDataProductName']}:: {src_info[_v_src]['description']}"
            gfl_var_info[
                'sourceVariableName'] = f"{tar_info[_v_tar]['sourceDataProductName']}: {tar_info[_v_tar]['sourceVariableName']}, QC: {tar_info[_v_tar]['QC']}({qc_thres});{src_info[_v_src]['sourceDataProductName']}: {src_info[_v_src]['sourceVariableName']}"
            gfl_var_info[
                'data_path'] = f"{tar_info[_v_tar]['data_path']}; {src_info[_v_src]['data_path']}"
            gfl_var_info['with_gaps'] = _v_tar

        _exp_settings['dataset'][dataset_key]['variables'][
            _v_tar_gfld] = gfl_var_info
        _exp_settings['dataset'][dataset_key]['variables'].pop(_v_tar)
    return _exp_settings


def write_netcdf(ds, exp_settings, site, resampled=None):
    if resampled is None:
        data_freq = exp_settings["temporal_resolution"]
        ds = add_date_variables(ds)
    else:
        data_freq = resampled
    ncdir = os.path.join(exp_settings['output_dir_path']['nc_file'], data_freq)
    os.makedirs(ncdir, exist_ok=True)
    ncfile_name = os.path.join(
        ncdir, site + "." + exp_settings['start_date'].split('-')[0] + "." +
        str(int(exp_settings['end_date'].split('-')[0])) + '.' + data_freq +
        '.nc')
    ds.to_netcdf(ncfile_name, mode='w')
    shut.log_and_print_sep(f'ncFile: {data_freq} : {ncfile_name}'.center(
        150, ' '))
    return


def resample_data(_data, _res):

    if _res == 'monthly':
        _data_out = _data.resample(time='1M').mean()
    elif _res == 'annual':
        _data_out = _data.resample(time='1Y').mean()
    else:
        sys.exit(f'{_res} resampling is not implemented')
    return _data_out


def harmonize_data_attrs(site_info, data, _config, resampled=None):
    site = site_info['site_ID']

    #%Sort variable alphabeticatly
    var_selected_data = sorted(list(data.keys()))
    data = data[var_selected_data]
    if resampled is None:
        data_freq = _config["temporal_resolution"]
    else:
        data_freq = resampled

    #%% Add metadata
    for source_ in _config["sel_datasets"]:
        var_info = _config["dataset"][source_]['variables']
        var_selected_src = list(
            _config["dataset"][source_]['variables'].keys())
        var_selected = list(
            set(var_selected_src).intersection(set(var_selected_data)))

        for var_ in var_selected:
            units = var_info[var_]['variableUnit']
            t_units = ['mm'.lower(), 'gC m-2'.lower(), 'MJ m-2'.lower()]
            if units.lower() in t_units:
                if var_info[var_]['isCarbon'] or var_info[var_][
                        'isEnergy'] or var_info[var_]['isWater']:
                    if _config['temporal_resolution'] == 'daily':
                        units = units + ' d-1'
                    elif _config['temporal_resolution'] == 'hourly':
                        units = units + ' h-1'
                    else:
                        units = units

            sourceDataProductName = var_info[var_]['sourceDataProductName']
            if sourceDataProductName == 'FLUXNET':
                sourceDataProductName = f'{sourceDataProductName}-{_config["FLUXNET_version"]}'
            if 'time' in data[var_].dims and 'depth_soilGrids' in data[
                    var_].dims:
                nameLong = f"{data_freq} {var_info[var_]['nameLong']} from {sourceDataProductName}"
            else:
                nameLong = f"{var_info[var_]['nameLong']} from {sourceDataProductName}"

            src_path = _config["fluxcom_cube_path"]
            if 'data_path' in var_info[var_].keys():
                src_path = var_info[var_]['data_path']
            attrs_dict = dict(
                bounds=str([np.nanmin(data[var_]),
                            np.nanmax(data[var_])]),
                units=units,
                short_name=var_info[var_]['nameShort'],
                long_name=nameLong,
                sourceDataProductName=sourceDataProductName,
                publication=var_info[var_]['publication'],
                partitioning=var_info[var_]['partitioning'],
                sourceVariableName=var_info[var_]['sourceVariableName'],
                description=var_info[var_]['description'],
                data_path=src_path)

            data[var_] = data[var_].assign_attrs(
                collections.OrderedDict(sorted(attrs_dict.items())))
            shut.log_and_print('attrs',
                               site,
                               data[var_].sourceDataProductName,
                               var_,
                               data[var_].units,
                               data_freq,
                               _config,
                               tar_label='variable',
                               src_label='units')

    #%% Transform time dim to days since time dimension
    if resampled is None:
        sdate = _config['start_date'].split('-')
        if _config["temporal_resolution"] == 'daily':
            date_base = datetime(int(sdate[0]), int(sdate[1]), int(sdate[2]),
                                 00, 00, 00)
            date_arr = np.array([
                date_base + timedelta(days=i) for i in range(len(data['time']))
            ])
            data['time'] = date_arr.astype(np.datetime64)
        if _config["temporal_resolution"] == 'hourly':
            date_base = datetime(int(sdate[0]), int(sdate[1]), int(sdate[2]),
                                 00, 00, 00)
            date_arr = np.array([
                date_base + timedelta(hours=i)
                for i in range(len(data['time']))
            ])
            data['time'] = date_arr.astype(np.datetime64)

    #%% Add global attributes
    CurrentScript = os.path.basename(__file__)
    data = data.assign_attrs(
        title=f'site dataset for {_config["exp_name"]}',
        FLUXNET_version=_config["FLUXNET_version"],
        temporal_resolution=data_freq,
        last_disturbance_on=_config['last_disturbance_on'],
        PFT=_config['sitePFT'],
        PFT_extractor=_config['PFTextractor'],
        created_by='Sujan Koirala',
        contact='skoirala@bgc-jena.mpg.de',
        SITE_ID=site_info["site_ID"],
        latitude=site_info["latitude"],
        longitude=site_info["longitude"],
        start_year=_config['start_date'].split('-')[0],
        end_year=_config['end_date'].split('-')[0],
        creation_date=datetime.now().strftime("%d-%m-%Y %H:%M"),
        creation_script=CurrentScript)
    return data


def close_logger():
    log = logging.getLogger()  # root logger - Good to get it only once.
    for hdlr in log.handlers[:]:  # remove the existing file handlers
        if isinstance(hdlr, logging.FileHandler):  #fixed two typos here
            log.removeHandler(hdlr)
            hdlr.close()
    return


def update_out_dir(_exp_settings, data_variant):
    __exp_settings = copy.deepcopy(_exp_settings)
    if __exp_settings['diagnostic_plots']:
        __exp_settings['output_dir_path']['figs'] = os.path.join(
            __exp_settings['output_dir_path']['figs'], data_variant)
    __exp_settings['output_dir_path']['nc_file'] = os.path.join(
        __exp_settings['output_dir_path']['nc_file'], data_variant)
    return __exp_settings


def finalize_and_save(data_dict,
                      site,
                      site_info,
                      exp_settings,
                      gap_filled_data=None,
                      gap_filler_key=None):
    shut.log_and_print_sep(
        f'harmonizing, plotting, and saving for gapfiller: {gap_filler_key}'.
        center(150, '-'))
    _exp_settings = copy.deepcopy(exp_settings)
    _data_dict = copy.deepcopy(data_dict)
    try:
        data_values = _data_dict.values()
        if len(data_values) == 0:
            exit_msg = f"None of the selected datasets provide any data at {_exp_settings['temporal_resolution']} resolution. Selected datasets were [{','.join(_exp_settings['sel_datasets'])}]"
            logging.warning(exit_msg)
            raise Exception(exit_msg)
    except Exception as e:
        logging.info(
            f"Cannot compile data for: {site_info['site_ID']} for  {site_info['site_ID']}. Error: "
        )
        print(e)
        logging.info(e)
        data_values=None

    if gap_filled_data is None:
        if _exp_settings['do_gap_fill'] == False:
            _exp_settings = _exp_settings
        else:
            _exp_settings = update_out_dir(_exp_settings, _exp_settings["FLUXNET_version"])
            for _gf in _exp_settings['sel_gapfills']:
                _data_dict.pop(_exp_settings['gap_fill'][_gf]['source'], None)
    else:
        gap_filler = _exp_settings['gap_fill'][gap_filler_key]['source']
        gap_filled = _exp_settings['gap_fill'][gap_filler_key]['target']
        _data_dict['data_gfld'] = gap_filled_data
        _exp_settings = update_out_dir(
            _exp_settings, _exp_settings['dataset'][gap_filler]["src_cliff"])
        _exp_settings["start_date"] = _exp_settings['dataset'][gap_filler][
            "sel_src_start_date"]
        _exp_settings["end_date"] = _exp_settings['dataset'][gap_filler][
            "sel_src_end_date"]
        sel_data_gf_only = [gap_filler, gap_filled, f'{gap_filler_key}_gfld']
        for _gf in _exp_settings['sel_gapfills']:
            if _gf != gap_filler_key:
                _data_dict.pop(_exp_settings['gap_fill'][_gf]['source'], None)

    #%% Merge all data
    logging.info(f'Merging all data: {site_info["site_ID"]}')
    ds = xr.merge(_data_dict.values())

    # remove gap_fill_data from dictionary
    if 'data_gfld' in _data_dict.keys():
        _data_dict.pop('data_gfld')

    # slice the time series
    ds = ds.sel(
        time=slice(_exp_settings["start_date"], _exp_settings["end_date"]))

    # get only a subset to plot

    #%% resample, harmonize metadata, and save
    res_harmon = [None]
    if len(_exp_settings["resample_output"]) > 0:
        res_harmon = _exp_settings["resample_output"]
        res_harmon.insert(0, None)
    for res in res_harmon:
        if res is not None:
            ds_res = resample_data(ds, res)
        else:
            ds_res = ds
        ds_res_h = harmonize_data_attrs(site_info,
                                        ds_res,
                                        _exp_settings,
                                        resampled=res)
        write_netcdf(ds_res_h, _exp_settings, site, resampled=res)
        if gap_filler_key is not None:
            _exp_settings['sel_datasets'] = sel_data_gf_only
        if _exp_settings['diagnostic_plots']:
            plot_dianostic_figures(ds_res_h,
                                   site_info,
                                   _exp_settings,
                                   site,
                                   resampled=res)
        if res is None:
            res = _exp_settings['temporal_resolution']
        shut.log_and_print_sep(
            f'harmonizing, plotting, and saving complete {res}'.center(
                150, '-'))

    shut.log_and_print_sep()

    site_dict_d = ds_res_h.attrs

    # close datasets
    ds_res.close()
    ds.close()

    return site_dict_d


def compile_site_data(site, exp_settings=None):
    log_path = os.path.join(
        exp_settings['output_dir_path']['log'],
        f"{site}_{exp_settings['output_dir_path']['expName']}.log")
    filehandler = logging.FileHandler(log_path, 'w')
    formatter = logging.Formatter(
        '%(levelname)s::%(filename)s::\n\t%(funcName)s::%(lineno)d::\n\t%(message)s'
    )
    filehandler.setFormatter(formatter)
    log = logging.getLogger()
    log.addHandler(filehandler)  # set the new handler
    log.setLevel(logging.INFO)
    # logging.basicConfig(filename=log_path, filemode='w', level=logging.INFO)
    shut.log_and_print_sep()
    shut.log_and_print_sep(
        f'Initialization of the data processing: {site}'.center(150, ' '))
    shut.log_and_print_sep()
    logging.info(f'Processing data for site {site}')
    #%% Get site information
    site_info = get_site_info(site)

    #%% get all data
    data_all = {}
    sitePFT = ''
    PFTextractor = 'None'
    last_disturbance_on = ''
    for extractor in exp_settings["sel_datasets"]:
        extractor_name = exp_settings["dataset"][extractor]["extractor"]
        logging.info(
            f'Getting dataset| configuration: {extractor}, extractor: {extractor_name}, site: {site_info["site_ID"]}'
        )
        extractor_mod = load_extractor(extractor_name)
        extractor_function = getattr(extractor_mod, 'extract')
        if exp_settings["allow_extraction_errors"]:
            try:
                provided_data = extractor_function(site_info=site_info,
                                                config=exp_settings,
                                                dataset=extractor)
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                print(e)     
                logging.info(
                    f'Failed getting dataset| configuration: {extractor}, extractor: {extractor_name}, site: {site_info["site_ID"]} due to error {e}. Allowing to continue because allow_extraction_errors in exp[json] is set to {exp_settings["allow_extraction_errors"]}'
                )
                provided_data = None
        else:
                provided_data = extractor_function(site_info=site_info,
                                                config=exp_settings,
                                                dataset=extractor)

        if provided_data is not None:
            data_all[extractor] = provided_data
            if (len(sitePFT) == 0) & ('PFT' in provided_data.attrs):
                sitePFT = f"{provided_data.attrs['PFT']}"
                PFTextractor = f"{extractor} using {extractor_name}"
            if (len(last_disturbance_on) == 0) & ('last_disturbance_on'
                                                  in provided_data.attrs):
                last_disturbance_on = f"{provided_data.attrs['last_disturbance_on']}"

    exp_settings['last_disturbance_on'] = last_disturbance_on
    exp_settings['sitePFT'] = sitePFT
    exp_settings['PFTextractor'] = PFTextractor

    site_dict_d = finalize_and_save(data_all, site, site_info, exp_settings)
    if exp_settings['do_gap_fill'] == True:
        for gfik in exp_settings['sel_gapfills']:
            gfiv = exp_settings["gap_fill"][gfik]
            if gfiv['source'] in data_all.keys() and gfiv['target'] in data_all.keys():
                data_gapfilled = cliff_gapfill(fill_src=gfiv['source'],
                                            fill_tar=gfiv['target'],
                                            site_info=site_info,
                                            data_dict=data_all,
                                            qc_thres=gfiv['qc_thres'],
                                            config=exp_settings).gapfill()
                gf_exp_settings = add_gapfill_settings(exp_settings,
                                                    fill_src=gfiv['source'],
                                                    fill_tar=gfiv['target'],
                                                    dataset_key=f'{gfik}_gfld',
                                                    qc_thres=gfiv['qc_thres'])
                site_dict_d = finalize_and_save(data_all,
                                                site,
                                                site_info,
                                                gf_exp_settings,
                                                gap_filled_data=data_gapfilled,
                                                gap_filler_key=gfik)
    #%% Export site data
    site_dict = {}
    site_dict[site] = {
            'PFT': 'NA',
            'last_disturbance_on': 'NA',
            'latitude': site_info['latitude'],
            'longitude': site_info['longitude'],
            'status': 'success'
        }
    # site_dict[site] = {
    #         'start_year': exp_settings['start_date'].split('-')[0],
    #         'end_year': exp_settings['end_date'].split('-')[0],
    #         'latitude': site_info['latitude'],
    #         'longitude': site_info['longitude'],
    #         'status': 'success'
    #     }
        # add pft to the settings
    if site_dict_d is not None:
        if 'PFT' in list(site_dict_d.keys()):
            site_dict[site]['PFT'] = site_dict_d['PFT']
        if 'last_disturbance_on' in list(site_dict_d.keys()):
            site_dict[site]['last_disturbance_on'] = site_dict_d['last_disturbance_on']
    else:
        site_dict[site]['status'] = 'failed'


    shut.log_and_print_sep(
        f'Successful data processing for site: {site}'.center(150, ' '))

    close_logger()

    return site_dict


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Utilities: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('workflow: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)

# %%
