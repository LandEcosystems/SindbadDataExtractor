#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tower meteorology gap filled by PyCLIFF with ERA5, read from its per-site products
<source_dir>/<fn>.<HR|DD>/<db>/<site>.<start>.<end>.nc. fn names the eddy covariance
source whose tower data PyCLIFF filled, such as FLUXNET2015 or Shuttle, and HR or DD
the hourly or the daily PyCLIFF run. pycliff_site.py runs PyCLIFF per site and reads
its product with this module.

Each variable of the json reads the PyCLIFF variable of its sourceVariableName on
the time dimension of the file:
    <VAR>_F       the tower value where its QC is at least 0.85, else the fill
    <VAR>_F_QC    hourly: 0 tower, 1 selected method, 2 nearest ERA5 pixel;
                  daily: the fraction of the day from the tower
    <VAR>_SD      the spread of the series across the PyCLIFF methods
    PET           Priestley-Taylor potential evapotranspiration
The time of the file is an integer YYYYMMDDhhmm in local standard time, as the tower
data, and missing values are NaN. The values are taken as they are: with the default
settings units si and radiation MJ, PyCLIFF writes temperature in deg C, VPD and air
pressure in Pa, despite their hPa and kPa labels, radiation in MJ m-2 and
precipitation in mm per time step, the units of the tower data in the cubes. Every
variable gets the selected method of its PyCLIFF variable as pycliff_method, from
the methods attribute of the file, so that a fallback to the nearest pixel, RAW,
shows in the cube.
"""
import glob
import json
import os

import numpy as np
import pandas as pd
import xarray as xr

import logging
logger = logging.getLogger(__name__)

# Suffixes of the PyCLIFF variables after the name of the site variable.
SUFFIXES = ("_F_QC", "_F", "_DSC", "_SD")
FLAG_MEANINGS = "tower_value_QC_ge_0.85 filled_by_the_selected_method " \
    "filled_by_the_nearest_ERA5_pixel_RAW"


def site_file(info, site, freq):
    """Path of the PyCLIFF product of a site, or None."""
    run = "DD" if freq == "daily" else "HR"
    folder = os.path.join(info["source_dir"], f"{info['fn']}.{run}", info["db"])
    files = sorted(glob.glob(os.path.join(folder, f"{site}.*.nc")))
    return files[-1] if files else None


def source_files(info, site, freq):
    """Files read for a site, whose changes stage_product.py includes in its hash."""
    path = site_file(info, site, freq)
    return [path] if path else []


def base_variable(name):
    """Site variable of a PyCLIFF variable, TA_DayTime of TA_DayTime_F_QC."""
    for sfx in SUFFIXES:
        if name.endswith(sfx):
            return name[:-len(sfx)]
    return name


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    info = config["dataset"][dataset]
    freq = bxtr.temporal_resolution
    date_ = bxtr.get_date_vec()
    path = site_file(info, bxtr.site, freq)
    src, methods, times = None, {}, None
    if path is None:
        logger.warning(f"::MISSING:: {dataset} has no PyCLIFF product of {info['fn']} "
                       f"for {bxtr.site} in {info['source_dir']}. NaN values will be "
                       f"set.")
    else:
        with xr.open_dataset(path, decode_times=False) as f:
            src = f.load()
        methods = json.loads(src.attrs.get("methods", "{}"))
        times = pd.to_datetime(src["time"].values.astype(np.int64).astype(str),
                               format="%Y%m%d%H%M")
        if freq == "daily":
            times = times.floor("D")

    src_data, var_attrs = [], {}
    for tar_name in bxtr.vars_list:
        var = bxtr.vars[tar_name]
        bxtr.log_var_start(tar_name)
        name = var["sourceVariableName"]
        values = np.full(len(date_), np.nan)
        attrs = {}
        if src is not None and name in src and src[name].dims == ("time",):
            values = pd.Series(src[name].values.astype(float), index=times).groupby(
                level=0).first().reindex(pd.DatetimeIndex(date_)).values
            attrs = {"pycliff_file": path,
                     "pycliff_method": methods.get(base_variable(name), "")}
            if name.endswith("_F_QC") and freq == "hourly":
                attrs["flag_values"] = np.array([0., 1., 2.])
                attrs["flag_meanings"] = FLAG_MEANINGS
        elif src is not None:
            logger.warning(f"::MISSING:: {bxtr.site} has no {name} in {path} for "
                           f"{tar_name}. NaN values will be set.")
        var_attrs[tar_name] = attrs
        data = xr.Dataset({tar_name: xr.DataArray(
            data=values, dims=["time"], coords={"time": date_},
            attrs={"units": var["variableUnit"]})})
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    # merge_and_format makes new arrays, so the attributes are set afterwards.
    for v, attrs in var_attrs.items():
        src_dataset[v].attrs.update(attrs)
    return src_dataset


if __name__ == '__main__':
    print(__doc__)
