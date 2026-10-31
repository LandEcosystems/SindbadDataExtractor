#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transforms of fluxcom site cube data for BasexTractor.get_transform: the daily
aggregations and the unit scaling of hourly data. Only the extractors that read the
fluxcom site cube use them, so this is the only module of the extractors that needs
the fluxcom package, and BasexTractor imports it only when a transform is asked for.
"""
from fluxcom import transformers
from fluxcom.core.variables import Variable  # noqa: F401, used by get_transform


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
