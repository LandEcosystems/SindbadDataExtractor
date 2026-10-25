#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eddy covariance fluxes and meteorology of the FLUXNET Shuttle, the ONEFlux products
that AmeriFlux, ICOS and TERN publish per site, fetched once by prefetch.py.

prefetch lists all data hubs with the fluxnet-shuttle package and keeps the listing
in the folder of cache_path as fluxnet_shuttle_snapshot_<time>.csv, reusing one
younger than snapshot_max_age_h of the json. It then
downloads the zip of every site of the site list that the listing has, unless the
index already holds the same product file, so a new release of a site replaces the
old one. cache_path is an index csv with one row per downloaded zip, holding the
data hub, the product file, its DOI and citation, the ONEFlux version, the years
and fetched_utc. The package is pip installed from github.com/fluxnet/shuttle and is
only needed for prefetch.

extract reads the half-hourly FLUXMET file of the zip of the site. The source
variable of the json, with its partitioning for GPP and RECO, picks a column of
COLUMNS, where NEE and its partitioned fluxes use the ustar variant nee_variant of
the json, USTAR50 by default as in the FLUXCOM site cube. The json may map a source
variable to another column with columns. Half hours are averaged to hours, and
summed for water fluxes. A QC column becomes the fraction of half hours with QC 0
or 1, the measured and good quality gap filled data. Random uncertainties of two
half hours are combined as for a mean of independent errors. Daily values follow
BasexTractor.get_transform: sums for energy, carbon and water fluxes, and the
daytime, minimum, maximum and weighted means for the variables that ask for them.
"""
import glob
import io
import os
import re
import time
import zipfile

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

# Source variable, with the partitioning for GPP and RECO, to the ONEFlux column.
# {nee} is replaced by nee_variant of the json.
COLUMNS = {
    "TA": "TA_F", "TA_QC": "TA_F_QC",
    "SW_IN": "SW_IN_F", "SW_IN_QC": "SW_IN_F_QC",
    "LW_IN": "LW_IN_F", "LW_IN_QC": "LW_IN_F_QC",
    "VPD": "VPD_F", "VPD_QC": "VPD_F_QC",
    "PA": "PA_F", "PA_QC": "PA_F_QC",
    "P": "P_F", "P_QC": "P_F_QC",
    "WS": "WS_F", "WS_QC": "WS_F_QC",
    "SW_IN_POT": "SW_IN_POT", "PPFD_IN": "PPFD_IN", "RH": "RH", "NETRAD": "NETRAD",
    "CO2": "CO2_F_MDS", "CO2_QC": "CO2_F_MDS_QC",
    "H": "H_F_MDS", "H_QC": "H_F_MDS_QC", "H_RANDUNC": "H_RANDUNC",
    "LE": "LE_F_MDS", "LE_QC": "LE_F_MDS_QC", "LE_RANDUNC": "LE_RANDUNC",
    "NEE": "NEE_VUT_{nee}", "NEE_QC": "NEE_VUT_{nee}_QC",
    "NEE_RANDUNC": "NEE_VUT_{nee}_RANDUNC",
    "GPP_NT": "GPP_NT_VUT_{nee}", "GPP_DT": "GPP_DT_VUT_{nee}",
    "RECO_NT": "RECO_NT_VUT_{nee}", "RECO_DT": "RECO_DT_VUT_{nee}",
    "GPP_NT_QC": "NEE_VUT_{nee}_QC", "GPP_DT_QC": "NEE_VUT_{nee}_QC",
    "RECO_NT_QC": "NEE_VUT_{nee}_QC", "RECO_DT_QC": "NEE_VUT_{nee}_QC",
}
# Carbon flux from umol m-2 s-1 to gC m-2 per hour, as in BasexTractor.
UMOL_TO_GC_H = 3600 * 12.001 * 1e-6
W_TO_MJ_H = 3600 * 1e-6


def cache_folder(info):
    return os.path.dirname(info["cache_path"])


def latest_snapshot(info):
    files = sorted(glob.glob(os.path.join(cache_folder(info),
                                          "fluxnet_shuttle_snapshot_*.csv")))
    return files[-1] if files else ""


def snapshot(info, listall):
    """
    Listing of the data hubs, made again when the newest one in the cache folder is
    older than snapshot_max_age_h of the json, 24 hours by default.
    """
    path = latest_snapshot(info)
    max_age = info.get("snapshot_max_age_h", 24) * 3600
    if path and os.path.getmtime(path) > time.time() - max_age:
        return path
    return listall(output_dir=cache_folder(info))


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    from fluxnet_shuttle.shuttle import download, listall
    folder = cache_folder(info)
    os.makedirs(folder, exist_ok=True)
    snap_path = snapshot(info, listall)
    snap = pd.read_csv(snap_path, keep_default_na=False)
    print(f"snapshot {snap_path}: " + ", ".join(
        f"{h} {n}" for h, n in snap["data_hub"].value_counts().items()), flush=True)
    by_id = {s.lower(): r for s, r in zip(snap["site_id"], snap.to_dict("records"))}
    have = fu.read_cache(info["cache_path"])
    have_files = set(have["file"].astype(str)) if len(have) else set()
    todo, absent = [], []
    for site, lat, lon in zip(sites, lats, lons):
        row = by_id.get(site.lower())
        if row is None:
            absent.append(site)
        elif row["fluxnet_product_name"] not in have_files:
            todo.append((site, lat, lon, row))
    print(f"{len(todo)} of {len(sites)} sites to download, {len(absent)} not on the "
          "shuttle", flush=True)
    batch = workers or info.get("n_workers", 10)
    for i in range(0, len(todo), batch):
        part = todo[i:i + batch]
        try:
            download(site_ids=[r["site_id"] for *_, r in part],
                     snapshot_file=snap_path, output_dir=folder,
                     user_info={"ameriflux": {}})
        except Exception as e:
            print(f"batch from {part[0][0]} failed: {e}", flush=True)
        for site, lat, lon, row in part:
            index_site(info, site, lat, lon, row, folder, snap_path)


def index_site(info, site, lat, lon, row, folder, snap_path):
    """Add the index row of a downloaded zip, or report the site as failed."""
    path = os.path.join(folder, row["fluxnet_product_name"])
    if not os.path.exists(path) or not zipfile.is_zipfile(path):
        print(f"failed {site}: {row['fluxnet_product_name']} not downloaded",
              flush=True)
        return
    fu.append_cache(info["cache_path"], [{
        "site": site, "site_id": row["site_id"], "data_hub": row["data_hub"],
        "file": row["fluxnet_product_name"], "product_id": row["product_id"],
        "citation": row["product_citation"],
        "oneflux_version": row["oneflux_code_version"],
        "first_year": row["first_year"], "last_year": row["last_year"],
        "igbp": row["igbp"], "latitude": lat, "longitude": lon,
        "shuttle_latitude": row["location_lat"],
        "shuttle_longitude": row["location_long"],
        "snapshot": os.path.basename(snap_path)}])
    print(f"done {site}", flush=True)


def site_record(info, site):
    """Index row of the newest product of a site, or None."""
    rows = fu.site_rows(fu.read_cache(info["cache_path"]), site)
    if len(rows) == 0:
        return None
    return rows.sort_values(["last_year", "fetched_utc"]).iloc[-1]


def zip_member(zf, kind):
    """Name of the member of kind, such as FLUXMET_HH, ERA5_HH or BIF, in a zip."""
    pat = re.compile(rf"_FLUXNET_{kind}_\d{{4}}-\d{{4}}_")
    names = [n for n in zf.namelist() if pat.search(n)]
    if not names:
        raise KeyError(f"no {kind} file in {zf.filename}")
    return names[0]


def read_member(zip_path, kind):
    """
    Table of a FLUXNET csv in the zip, indexed by TIMESTAMP_START when it has one.
    """
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(zip_member(zf, kind)) as f:
            df = pd.read_csv(io.TextIOWrapper(f), na_values=[-9999, "-9999"])
    if "TIMESTAMP_START" in df:
        df.index = pd.to_datetime(df.pop("TIMESTAMP_START").astype(int).astype(str),
                                  format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_END"], errors="ignore")
    return df


def read_subdaily(zip_path):
    """Half hourly FLUXMET table of a site, or the hourly one of sites such as US-Ha1
    and AU-Tum that report hourly."""
    try:
        return read_member(zip_path, "FLUXMET_HH")
    except KeyError:
        return read_member(zip_path, "FLUXMET_HR")


def bif_value(zip_path, variable):
    """First DATAVALUE of a BIF variable, such as UTC_OFFSET or IGBP, or None."""
    bif = read_member(zip_path, "BIF")
    sel = bif.loc[bif["VARIABLE"] == variable, "DATAVALUE"]
    return sel.iloc[0] if len(sel) else None


def column_for(info, src_name, partitioning):
    key = src_name
    if src_name in ("GPP", "RECO") and partitioning:
        key = f"{src_name}_{partitioning}"
    columns = {**COLUMNS, **info.get("columns", {})}
    col = columns.get(key, "")
    return col.format(nee=info.get("nee_variant", "USTAR50"))


def to_hourly(df, col, kind):
    """
    Hourly series of a column, kind being value, sum, qc or randunc. Steps of an hour
    are kept, half hours are combined in pairs starting on the full hour.
    """
    x = df[col].astype(float)
    if kind == "qc":
        x = (x <= 1).astype(float).where(x.notna())
    step = (df.index[1] - df.index[0]) if len(df) > 1 else pd.Timedelta("30min")
    if step >= pd.Timedelta("60min"):
        return x
    hours = x.index.floor("h")
    if kind == "sum":
        return x.groupby(hours).sum(min_count=2)
    if kind == "randunc":
        return np.sqrt((x ** 2).groupby(hours).sum(min_count=2)) / 2
    if kind == "qc":
        return x.groupby(hours).mean()
    return x.groupby(hours).mean().where(x.groupby(hours).count() == 2)


def daytime_weight(sw_in_pot):
    """Weights of DaytimeMean in BasexTractor for days by hours of SW_IN_POT."""
    thres = np.minimum(sw_in_pot / 0.1, 1)
    w = (thres / np.nanmean(thres)) * (sw_in_pot / np.nanmean(sw_in_pot))
    return np.minimum(w, 1)


def to_daily(v, tar_name, var, weights):
    """Daily values from days by hours v, following BasexTractor.get_transform."""
    src_u = var["sourceVariableUnit"].lower()
    tar_u = var["variableUnit"].lower()
    if var["isEnergy"] and src_u.startswith("w") and tar_u.startswith("mj"):
        return v.sum(axis=1) * W_TO_MJ_H
    if var["isCarbon"] and src_u.startswith("umol") and tar_u.startswith("gc"):
        return v.sum(axis=1) * UMOL_TO_GC_H
    if var["isWater"] or "DaySum" in tar_name:
        return v.sum(axis=1)
    if "DayTime" in tar_name:
        w = weights["SW_IN_POT"]
        out = np.nansum(v * w, axis=1) / np.nansum(w, axis=1)
        return np.where(np.isnan(v).any(axis=1), np.nan, out)
    if "DayMin" in tar_name:
        return v.min(axis=1)
    if "DayMax" in tar_name:
        return v.max(axis=1)
    if "wDayMean" in tar_name:
        w = weights[var["weightVar"]]
        return (v * w).sum(axis=1) / w.sum(axis=1)
    return v.mean(axis=1)


def hourly_scale(tar_name, var):
    src_u = var["sourceVariableUnit"].lower()
    tar_u = var["variableUnit"].lower()
    if var["isEnergy"] and src_u.startswith("w") and tar_u.startswith("mj"):
        return W_TO_MJ_H
    if var["isCarbon"] and src_u.startswith("umol") and tar_u.startswith("gc"):
        return UMOL_TO_GC_H
    return 1.0


def kind_of(src_name, var):
    if src_name.endswith("_QC"):
        return "qc"
    if src_name.endswith("_RANDUNC"):
        return "randunc"
    return "sum" if var["isWater"] else "value"


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    hours = np.arange(np.datetime64(bxtr.start_date, "h"),
                      np.datetime64(bxtr.end_date, "D") + np.timedelta64(1, "D"),
                      np.timedelta64(1, "h"))
    daily = bxtr.temporal_resolution == "daily"

    rec = site_record(info, bxtr.site)
    zip_path = os.path.join(cache_folder(info), str(rec["file"])) if rec is not None \
        else ""
    if not zip_path or not os.path.exists(zip_path):
        logger.warning(f"::MISSING:: {dataset} has no shuttle product for {bxtr.site} "
                       f"in {info['cache_path']}. NaN values will be set.")
        df = None
    else:
        df = read_subdaily(zip_path)

    def hourly_of(src_name, partitioning, var):
        col = column_for(info, src_name, partitioning)
        if df is None or not col or col not in df:
            return None, col
        x = to_hourly(df, col, kind_of(src_name, var))
        return x.reindex(pd.DatetimeIndex(hours)).values, col

    weights = {}
    if daily:
        sw, _ = hourly_of("SW_IN_POT", "", {"isWater": False})
        if sw is not None:
            weights["SW_IN_POT"] = daytime_weight(sw.reshape(-1, 24))

    src_data, columns = [], {}
    for tar_name in bxtr.vars_list:
        var = bxtr.vars[tar_name]
        bxtr.log_var_start(tar_name)
        h, col = hourly_of(var["sourceVariableName"], var["partitioning"], var)
        if h is None:
            if df is not None:
                logger.warning(f"::MISSING:: {bxtr.site} has no column {col or '?'} "
                               f"for {tar_name}. NaN values will be set.")
            h = np.full(len(hours), np.nan)
        if daily:
            v = h.reshape(-1, 24)
            if var.get("weightVar") and "wDayMean" in tar_name:
                wv = bxtr.vars.get(var["weightVar"], {})
                w, _ = hourly_of(wv.get("sourceVariableName", var["weightVar"]),
                                 wv.get("partitioning", ""), {"isWater": False})
                weights[var["weightVar"]] = w.reshape(-1, 24) if w is not None \
                    else np.full(v.shape, np.nan)
            values = to_daily(v, tar_name, var, weights)
        else:
            values = h * hourly_scale(tar_name, var)
        data = xr.Dataset({tar_name: xr.DataArray(
            data=values, dims=['time'], coords={'time': date_},
            attrs={"units": var["sourceVariableUnit"]})})
        data = bxtr.convert_units(data, tar_name)
        columns[tar_name] = col
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    for v, col in columns.items():
        src_dataset[v].attrs["source_column"] = col
    if rec is not None:
        attrs = {"product": rec["file"], "product_id": rec["product_id"],
                 "citation": rec["citation"], "data_hub": rec["data_hub"],
                 "oneflux_version": rec["oneflux_version"],
                 "nee_variant": info.get("nee_variant", "USTAR50")}
        rows = fu.site_rows(fu.read_cache(info["cache_path"]), bxtr.site)
        for v in src_dataset.data_vars:
            src_dataset[v].attrs.update(attrs)
        src_dataset = fu.with_fetch_attrs(src_dataset, rows, info)
    return src_dataset


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
