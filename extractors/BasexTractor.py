#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
from fluxcom import transformers
import fluxcom.providers as flpr
from fluxcom.core.variables import Variable
import xarray as xr
import utils.shared_utils as shut

import logging

logger = logging.getLogger(__name__)
import inspect
import numpy as np


class BasexTractor:
    def __init__(self, dataset, site_info, config):
        """
        instantiate a class with a set of attributes that are common to all data extractors
        """
        self.site = site_info["site_ID"]
        self.lat = site_info["latitude"]
        self.lon = site_info["longitude"]
        self.dataset = dataset
        self.flx_cubepath = config["fluxcom_cube_path"]
        self.version = config["FLUXNET_version"]
        self.vars = config["dataset"][dataset]["variables"]
        self.vars_list = list(config["dataset"][dataset]["variables"].keys())
        self.temporal_resolution = config["temporal_resolution"]
        self.vars = config["dataset"][dataset]["variables"]
        self.extractor = config["dataset"][dataset]["extractor"]
        self.config = config

        if ("src_start_date" in config["dataset"][dataset]) & (
                "src_end_date" in config["dataset"][dataset]):
            self.start_date = config["dataset"][dataset]["src_start_date"]
            self.end_date = config["dataset"][dataset]["src_end_date"]
        else:
            self.start_date = config["start_date"]
            self.end_date = config["end_date"]
        if "src_cliff" in config["dataset"][dataset]:
            self.source = config["dataset"][dataset]["src_cliff"]

    def is_resolution_supported(self, provided_reso=None):
        reso_heirarchy = ['hourly', 'daily', 'monthly', 'annual']
        if provided_reso is None:
            return True
        else:
            t_reso = self.temporal_resolution
            t_reso_ind = reso_heirarchy.index(t_reso)
            provided_reso_ind = reso_heirarchy.index(provided_reso)
            if provided_reso_ind > t_reso_ind:
                import sys
                logger.warning(
                    f"The selected dataset {self.dataset}, using {self.extractor} extractor, only provides {provided_reso} data but the resolution is {self.temporal_resolution} (coded in {self.extractor}.py in extractors). The following variables will not be included in {self.temporal_resolution} data for {self.site}: \n {list(self.vars.keys())}"
                )
                return False
            else:
                return True

    def convert_units(self, data, tar_name):
        """
        - sets the units attribute in data array when it does not exists.
        - converts the units.
        - takes  dataset or dataarray and returns the same
        - xarray classes: xarray.core.dataarray.DataArray for datarray and  xarray.core.dataset.Dataset for Dataset
        """
        if isinstance(data, xr.core.dataset.Dataset):
            data_tmp = data[tar_name]
            is_dataset = True
        else:
            is_dataset = False
            data_tmp = data
        if "units" not in data_tmp.attrs:
            data_tmp.attrs["units"] = self.vars[tar_name]["sourceVariableUnit"]
            logger.warning(
                f"Data array without units for {tar_name}. Setting the units to sourceVariableUnit [{self.vars[tar_name]['sourceVariableUnit']}] from json."
            )
        data_tmp = shut.do_unit_conversion(
            data_tmp,
            self.vars[tar_name]["sourceVariableName"],
            self.vars[tar_name]["sourceVariableUnit"],
            self.vars[tar_name]["variableUnit"],
            self.vars[tar_name]["source2sindbadUnit"],
        )
        if is_dataset:
            data[tar_name] = data_tmp
        else:
            data = data_tmp
        return data

    def get_date_vec(self):
        if self.temporal_resolution == "daily":
            dtype_metadata = "M8[D]"
        else:
            dtype_metadata = "M8[h]"
        date_ = np.arange(
            np.datetime64(self.start_date),
            np.datetime64(self.end_date) + np.timedelta64(1, "D"),
            dtype=dtype_metadata,
        )
        return date_

    def get_transform(self, tar_name, src_prov=None):
        if src_prov is None:
            return None
        else:
            src_units_l = self.vars[tar_name]["sourceVariableUnit"].lower()
            tar_units_l = self.vars[tar_name]["variableUnit"].lower()
            transform = []
            if self.temporal_resolution == "daily":
                if (self.vars[tar_name]["isEnergy"]
                        and src_units_l.startswith("w")
                        and tar_units_l.startswith("mj")):
                    transform.append(eDaySum())
                elif (self.vars[tar_name]["isCarbon"]
                      and src_units_l.startswith("umol")
                      and tar_units_l.startswith("gc")):
                    transform.append(cDaySum())
                elif self.vars[tar_name]["isWater"]:
                    transform.append(DaySum())
                elif "DayTime" in tar_name:
                    transform.append(
                        DaytimeMean(src_prov.get_data(Variable('SW_IN_POT'))))
                elif "DayMin" in tar_name:
                    transform.append(DayMin())
                elif "DayMax" in tar_name:
                    transform.append(DayMax())
                elif "DaySum" in tar_name:
                    transform.append(DaySum())
                elif "wDayMean" in tar_name:
                    transform.append(
                        wDayMean(
                            src_prov.get_data(
                                Variable(self.vars[tar_name]["weightVar"]))))
                else:
                    transform = transformers.hourly_to_daily
            elif self.temporal_resolution == "hourly":
                if (self.vars[tar_name]["isEnergy"]
                        and src_units_l.startswith("w")
                        and tar_units_l.startswith("mj")):
                    transform.append(wm_2_to_mjh_1())
                elif (self.vars[tar_name]["isCarbon"]
                      and src_units_l.startswith("umol")
                      and tar_units_l.startswith("gc")):
                    transform.append(umolm_2s_1_to_gCh_1())
                else:
                    transform = []
            return transform

    def log_var_end(self, data, tar_name, transforms):
        shut.log_datavar_info(
            self.dataset,
            self.site,
            self.temporal_resolution,
            self.vars[tar_name]["sourceVariableName"],
            tar_name,
            self.vars[tar_name]["variableUnit"],
            self.vars[tar_name]["sourceVariableUnit"],
            self.vars[tar_name]["bounds"],
            data,
            transforms,
            partitioning=self.vars[tar_name]["partitioning"],
        )
        if "_merged" in tar_name:
            self.warn_merged_nan(data, tar_name)
        return

    def warn_merged_nan(self, data, tar_name):
        """
        Loud warning when a merged variable has any NaN at the site, with the
        fraction of time steps that are NaN, since merged variables are meant to
        have a value everywhere.
        """
        da = data[tar_name] if tar_name in data else data
        vals = np.asarray(da.to_array() if hasattr(da, "to_array") else da, float)
        n_nan = int(np.isnan(vals).sum())
        if n_nan == 0:
            return
        frac = n_nan / vals.size
        bar = "!" * 85
        msg = (f"\n{bar}\n!!! NaN IN MERGED VARIABLE: {self.dataset} {tar_name} at "
               f"{self.site}\n!!! NaN fraction {frac:.4f} ({n_nan} of {vals.size} "
               f"values)\n!!! Check the in-situ csv and the product cache for this "
               f"site.\n{bar}")
        logger.warning(msg)
        print(f"{msg}\n", flush=True)

    def log_var_start(self, tar_name):
        shut.log_and_print(
            "xTrct",
            self.site,
            self.vars[tar_name]["sourceDataProductName"],
            tar_name,
            self.vars[tar_name]["sourceVariableName"],
            self.temporal_resolution,
            self.config
        )

    def merge_and_format(self, data):
        if self.temporal_resolution == "hourly":
            data = shut.flatten_hour_to_time(data)
        dataset = xr.merge(data)
        dataset = shut.data_structure_temporal_NoDepth(dataset, self.lat,
                                                       self.lon)
        return dataset


## class for transforming fluxcom-cube based data


class DayMin(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        return data.min(dim="hour", skipna=False)


class DayMax(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        return data.max(dim="hour", skipna=False)


class wm_2_to_mjh_1(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        data.values = data.values * 3600 * 1e-6
        return data


# class k_to_degC(transformers.BaseTransformer):
#     transform_dict = {"temporal_res": ("daily", "daily")}
#     @staticmethod
#     def apply_transform(data):
#         return data - 273.15


class k_to_degC(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        data_out = data.copy()
        data_out.values = data.values - 273.15
        # print("this function does not work or is called when added as a transform")
        data_out.attrs["units"] = "deg C"
        return data_out


class umolm_2s_1_to_gCh_1(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        data.values = data.values * ((3600 * 12.001) * 1e-6)
        return data


class DaySum(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        dataout = data.sum(dim="hour", skipna=False)
        return dataout


class eDaySum(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        dataout = data.sum(dim="hour", skipna=False)
        return dataout * 3600 * 1e-6


class cDaySum(transformers.temporal.BaseDayAgg):
    @staticmethod
    def day_agg(data):
        dataout = data.sum(dim="hour", skipna=False)
        return dataout * ((3600 * 12.001) * 1e-6)


class DaytimeMean(transformers.BaseTransformer):
    transform_dict = {"temporal_res": ("hourly", "daily")}

    def __init__(self, sw_in_pot):
        out = sw_in_pot.copy(deep=True)
        thres = sw_in_pot / 0.1
        tv1 = thres.values
        tv1[tv1 > 1] = 1
        thres.values = tv1
        thres = (thres / thres.mean(dim=["time", "hour"])) * (
            sw_in_pot / sw_in_pot.mean(dim=["time", "hour"]))
        tv2 = thres.values
        tv2[tv2 > 1] = 1
        out.values = tv2
        self.weight = out
        transform_dict = {"temporal_res": ("hourly", "daily")}

    def apply_transform(self, data):
        data_tmp = data.sum(dim="hour", skipna=False)
        weightVar = data.values * self.weight.values
        data.values = weightVar
        dataout = data.sum(dim="hour") / self.weight.sum(dim="hour")
        dataout = data_tmp / data_tmp * dataout
        return dataout  # .to_dataset(name=data.name)


class wDayMean(transformers.BaseTransformer):
    transform_dict = {"temporal_res": ("hourly", "daily")}

    def __init__(self, weight):
        self.weight = weight
        transform_dict = {"temporal_res": ("hourly", "daily")}

    def apply_transform(self, data):
        weightVar = data.values * self.weight.values
        data.values = weightVar
        dataout = data.sum(dim="hour", skipna=False) / self.weight.sum(
            dim="hour", skipna=False)
        return dataout  # .to_dataset(name=data.name)


def add_attrs(_xtract, _extra_fields):
    for exf, exv in _extra_fields.items():
        setattr(_xtract, exf, exv)
    return _xtract


if __name__ == "__main__":
    import inspect, os

    print("---------------------------------------------------")
    print(
        "Provider at: ",
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))),
    )
    print("Provider: ", inspect.getfile(inspect.currentframe()))
    print("---------------------------------------------------")
    print(__doc__)
