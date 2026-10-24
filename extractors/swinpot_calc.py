#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Potential incoming shortwave radiation SW_IN_POT of the site, computed with the
algorithm of ONEFlux, get_rpot_with_solar_noon in oneflux_steps/common/common.c,
instead of read from the ONEFlux output as in oneflux_swinpot_csvload.py.

The radiation at the top of the atmosphere on a horizontal surface is computed for
every minute of each day in local standard time, with the solar constant of 1376 W
m-2 and the Fourier series of Spencer for declination and Earth-Sun distance, as a
curve centred on 12:00. The curve is shifted so that its centre falls on the solar
noon of the NOAA algorithm for the longitude and UTC offset of the site, and averaged
over each half hour. The values are written as MJ m-2 per time step, the sum of the
half hours of each hour or day, like the ONEFlux csv extractor.

The UTC offset of local standard time comes from site_info when the run gives one,
and otherwise from the UTC column of the csv at sites_path, such as
sites_consolidated.csv of consolidate_sites.py. A site without an offset gets NaN.
"""
import numpy as np
import pandas as pd
import xarray as xr

import logging
logger = logging.getLogger(__name__)

SOLAR_CONST = 1376.0
# The value of pi in the ONEFlux code.
PI = 3.141592654
MINUTES = 1440


def julian_day(year, month, day):
    year = np.where(month <= 2, year - 1, year)
    month = np.where(month <= 2, month + 12, month)
    a = np.floor(year / 100)
    b = 2 - a + np.floor(a / 4)
    return (np.floor(365.25 * (year + 4716)) + np.floor(30.6001 * (month + 1)) + day
            + b - 1524.5)


def equation_of_time(t):
    """Equation of time in minutes for t in Julian centuries since J2000."""
    sec = 21.448 - t * (46.8150 + t * (0.00059 - t * 0.001813))
    eps = (23.0 + (26.0 + sec / 60.0) / 60.0
           + 0.00256 * np.cos(np.radians(125.04 - 1934.136 * t)))
    l0 = np.radians(np.mod(280.46646 + t * (36000.76983 + 0.0003032 * t), 360.0))
    e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
    m = np.radians(357.52911 + t * (35999.05029 - 0.0001537 * t))
    y = np.tan(np.radians(eps) / 2.0) ** 2
    et = (y * np.sin(2 * l0) - 2 * e * np.sin(m) + 4 * e * y * np.sin(m) * np.cos(2 * l0)
          - 0.5 * y * y * np.sin(4 * l0) - 1.25 * e * e * np.sin(2 * m))
    return np.degrees(et) * 4.0


def solar_noon(days, lon, utc_offset):
    """Solar noon of each day in whole minutes of local standard time."""
    year = days.astype("datetime64[Y]").astype(int) + 1970
    month = days.astype("datetime64[M]").astype(int) % 12 + 1
    day = (days - days.astype("datetime64[M]")).astype(int) + 1
    jd = julian_day(year, month, day)
    century = lambda j: (j - 2451545.0) / 36525.0
    lon_west = -lon
    noon = 720 + lon_west * 4 - equation_of_time(century(jd + lon_west / 360.0))
    noon = 720 + lon_west * 4 - equation_of_time(century(jd - 0.5 + noon / 1440.0))
    return np.floor(noon + 60 * utc_offset).astype(int)


def minute_curves(days, lat):
    """Potential radiation in W m-2 per minute of each day, centred on 12:00."""
    doy = (days - days.astype("datetime64[Y]")).astype(int) + 1
    tthet = 2.0 * PI * (doy[:, None] - 1.0) / 365.0
    omega = -15.0 * np.abs(12.0 - np.arange(MINUTES) / 60.0)[None, :]
    decl = (0.006918 - 0.399912 * np.cos(tthet) + 0.070257 * np.sin(tthet)
            - 0.006758 * np.cos(2 * tthet) + 0.000907 * np.sin(2 * tthet)
            - 0.002697 * np.cos(3 * tthet) + 0.00148 * np.sin(3 * tthet))
    lat_r = lat * PI / 180.0
    cos_zen = (np.sin(decl) * np.sin(lat_r)
               + np.cos(decl) * np.cos(lat_r) * np.cos(omega * PI / 180.0))
    rpot = SOLAR_CONST * (1.00011 + 0.034221 * np.cos(tthet) + 0.00128 * np.sin(tthet)
                          + 0.000719 * np.cos(2 * tthet) + 0.000077 * np.sin(2 * tthet))
    return np.maximum(rpot * cos_zen, 0.0)


def sw_in_pot(start, end, lat, lon, utc_offset, step_minutes=30):
    """
    SW_IN_POT in W m-2, the mean over each step of step_minutes from the first day of
    start to the last day of end in local standard time, as days by steps.
    """
    days = np.arange(np.datetime64(start, "D"), np.datetime64(end, "D") + 1)
    curves = minute_curves(days, lat)
    # The shift fills the minutes it moves in with zero, as in ONEFlux.
    shift = solar_noon(days, lon, utc_offset) - 720
    out = np.zeros_like(curves)
    for k, s in enumerate(shift):
        if s >= 0:
            out[k, s:] = curves[k, :MINUTES - s]
        else:
            out[k, :MINUTES + s] = curves[k, -s:]
    return out.reshape(len(days), MINUTES // step_minutes, step_minutes).mean(axis=2)


def site_utc(site_info, sites_path):
    """UTC offset of the site from site_info, else from the csv at sites_path."""
    utc = site_info.get("UTC", np.nan)
    if utc is not None and np.isfinite(utc):
        return float(utc)
    if sites_path:
        sites = pd.read_csv(sites_path, keep_default_na=False, na_values=[""])
        sel = sites.loc[sites["site"].str.lower() == site_info["site_ID"].lower(),
                        "UTC"]
        if len(sel) and np.isfinite(sel.iloc[0]):
            return float(sel.iloc[0])
    return np.nan


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    per_day = 24 if bxtr.temporal_resolution == "hourly" else 1

    utc = site_utc(site_info, info.get("sites_path", ""))
    if np.isnan(utc):
        logger.warning(f"::MISSING:: {dataset} has no UTC offset for {bxtr.site} in "
                       f"site_info or {info.get('sites_path', '')}. NaN values will "
                       "be set.")
        values = np.full(len(date_), np.nan)
    else:
        half_hours = sw_in_pot(bxtr.start_date, bxtr.end_date, bxtr.lat, bxtr.lon, utc)
        # Mean W m-2 over each time step times its length in seconds, in MJ m-2.
        steps = half_hours.reshape(len(half_hours), per_day, -1).mean(axis=2)
        values = steps.ravel() * 86400 / per_day * 1e-6

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        data = xr.Dataset({tar_name: xr.DataArray(data=values.copy(), dims=['time'],
                                                  coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    for v in src_dataset.data_vars:
        src_dataset[v].attrs["UTC_offset"] = utc
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
