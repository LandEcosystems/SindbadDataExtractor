#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads per-site NetCDF files that this tool wrote in an earlier version, such as the
LaThuile 2007 extractions of the fluxnetBGI2021 site cube, which are kept as a
read-only copy instead of extracting them again from the old fluxcom package.

source_dir of the dataset json is the folder of the files of one temporal
resolution, holding <site>.<start>.<end>.<resolution>.nc. Each variable of the json
is read from the variable of the same name, or from the name that rename maps it
to. The values are taken as they are, since the files hold SINDBAD units already; a
units attribute that differs from variableUnit is reported. The time series is set
on the dates of the run, so years outside the file are NaN, and the latitude and
longitude are those of the site registry. Every variable keeps the path of its file
as legacy_file, the data_path of the file when it has one, and the coordinates of
the file as source_latitude and source_longitude, which stage_product.py compares
with the registry. A site without a file, or a variable missing from it, gives NaN
with a ::MISSING:: warning.

python extractors/sitefile_nc.py --site-list <source_dir> <csv> writes the site list
of the files, with the columns Site_ID, Lat, Lon and UTC, for consolidate_sites.py
--site-list.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import xarray as xr

import logging
logger = logging.getLogger(__name__)

# Unit suffixes that harmonize_data_attrs adds to fluxes of a temporal resolution.
RATE_SUFFIXES = (" d-1", " h-1")


def site_file(source_dir, site, freq):
    """Path of the file of a site, or None."""
    files = sorted(glob.glob(os.path.join(source_dir, f"{site}.*.{freq}.nc")))
    return files[-1] if files else None


def file_coords(src):
    """Latitude and longitude of a site file, from its coordinates or attributes."""
    out = []
    for name in ("latitude", "longitude"):
        if name in src.coords:
            out.append(float(np.ravel(src[name].values)[0]))
        elif name in src.attrs:
            out.append(float(src.attrs[name]))
        else:
            return None
    return tuple(out)


def same_units(file_units, json_units):
    a, b = str(file_units).strip().lower(), str(json_units).strip().lower()
    for sfx in RATE_SUFFIXES:
        a, b = a.removesuffix(sfx), b.removesuffix(sfx)
    return a == b


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    info = config["dataset"][dataset]
    freq = bxtr.temporal_resolution
    rename = info.get("rename", {})
    date_ = bxtr.get_date_vec()
    path = site_file(info["source_dir"], bxtr.site, freq)
    src, coords = None, None
    if path is None:
        logger.warning(f"::MISSING:: {dataset} has no {freq} file for {bxtr.site} in "
                       f"{info['source_dir']}. NaN values will be set.")
    else:
        with xr.open_dataset(path) as f:
            src = f.load()
        coords = file_coords(src)
        src = src.squeeze([d for d in ("latitude", "longitude") if d in src.dims],
                          drop=True)

    src_data, var_attrs = [], {}
    for tar_name in bxtr.vars_list:
        var = bxtr.vars[tar_name]
        bxtr.log_var_start(tar_name)
        name = rename.get(tar_name, tar_name)
        values = np.full(len(date_), np.nan)
        attrs = {}
        if src is not None and name in src and "time" in src[name].dims:
            da = src[name]
            times = da["time"].values.astype(date_.dtype)
            values = pd.Series(da.values.astype(float), index=times).groupby(
                level=0).first().reindex(date_).values
            if "units" in da.attrs and not same_units(da.attrs["units"],
                                                      var["variableUnit"]):
                logger.warning(f"{bxtr.site} {tar_name}: units {da.attrs['units']} in "
                               f"{path}, {var['variableUnit']} in the json.")
            attrs["legacy_file"] = path
            if da.attrs.get("data_path"):
                attrs["data_path"] = da.attrs["data_path"]
        elif src is not None:
            logger.warning(f"::MISSING:: {bxtr.site} has no {name} in {path} for "
                           f"{tar_name}. NaN values will be set.")
        if coords is not None:
            attrs["source_latitude"], attrs["source_longitude"] = coords
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


def site_list(source_dir, out_csv):
    """Write the sites of the files of source_dir with their coordinates."""
    rows = {}
    for f in sorted(glob.glob(os.path.join(source_dir, "*.nc"))):
        site = os.path.basename(f).split(".")[0]
        with xr.open_dataset(f) as src:
            coords = file_coords(src)
            utc = src.attrs.get("UTC", np.nan)
        lat, lon = coords if coords is not None else (np.nan, np.nan)
        rows[site] = {"Site_ID": site, "Lat": lat, "Lon": lon, "UTC": utc}
    pd.DataFrame(list(rows.values()), columns=["Site_ID", "Lat", "Lon", "UTC"]).to_csv(
        out_csv, index=False)
    print(f"{len(rows)} sites of {source_dir} written to {out_csv}")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--site-list", nargs=2, metavar=("SOURCE_DIR", "CSV"),
                   required=True)
    args = p.parse_args()
    site_list(*args.site_list)
