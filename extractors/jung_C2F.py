#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Daily C2F consistency flags of Jung et al. 2024, Biogeosciences 21, 1827-1846,
DOI:10.5194/bg-21-1827-2024, computed by prefetch.py from the FLUXNET Shuttle data of
each site with the published code c2f_fluxnet v1.0, DOI:10.5281/zenodo.10593332.

Prefetch runs in the conda env of sinproc_installer/sinproc.yml, with the unzipped
code at code_dir of the json. load_c2f patches the code and missingpy at import, so
that they run with numpy 2 and scikit-learn 1.8, which give slightly different
random forests and hence a few different flags than the old versions the code was
written for.

For every site of the site list with a zip in the index at
shuttle_cache_path, the half hours of the zip are turned into the daily input that
fluxcom_load of c2f_fluxnet builds from the FLUXCOM site cube: QC as the fraction of
half hours with QC 0 or 1, precipitation as a daily sum, night means where
SW_IN_POT is zero, the spread of the NEE ustar percentiles 16 to 84, and the
downscaled ERA5 air temperature and VPD of the ONEFlux file in place of the ERA5
cube. Each site is flagged on its own with the strictness of the json, loose by
default, after seeding the random number generator of numpy with seed, so that the
random forests give the same flags on every run. The across site step, which sets
the site and break flags, is not run.

The flags are written to the folder of cache_path as
<site>_jung_C2F_<strictness>.nc, with one <VAR>_C2F_Flag per variable that is 1 on
inconsistent days. cache_path is an index csv with one row per site and product
file, so a new shuttle release of a site is flagged again.

extract writes the flags as their own QC variables and never masks the fluxes, so
that they are applied when the data are used. It returns 1 minus the flag, so that 1
is consistent and 0 is inconsistent as in jung_QC, and repeats the daily value for
every hour at hourly resolution, as the paper advises to drop all sub daily values
of a flagged day. Days outside the flagged record are NaN.
"""
import os
import sys
import types

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu
import extractors.fluxnet_shuttle as fs

import logging
logger = logging.getLogger(__name__)

# Daily means of these ONEFlux columns, under the names of fluxcom_load.
MEAN_COLUMNS = {
    "TA": "TA_F_MDS", "H": "H_F_MDS", "VPD": "VPD_F_MDS", "LE": "LE_F_MDS",
    "GPP_NT": "GPP_NT_VUT_{nee}", "GPP_DT": "GPP_DT_VUT_{nee}",
    "RECO_NT": "RECO_NT_VUT_{nee}", "RECO_DT": "RECO_DT_VUT_{nee}",
    "NEE": "NEE_VUT_{nee}", "SW_IN": "SW_IN_F_MDS", "LW_OUT": "LW_OUT",
    "SW_IN_POT": "SW_IN_POT", "NETRAD": "NETRAD", "SWC_1": "SWC_F_MDS_1",
    "SWC_2": "SWC_F_MDS_2", "USTAR": "USTAR", "WS": "WS_F", "PPFD_IN": "PPFD_IN",
    "ERA__5__t2m": "TA_ERA", "ERA__5__vpd": "VPD_ERA",
}
QC_COLUMNS = {
    "TA_QC": "TA_F_MDS_QC", "H_QC": "H_F_MDS_QC", "VPD_QC": "VPD_F_MDS_QC",
    "LE_QC": "LE_F_MDS_QC", "NEE_QC": "NEE_VUT_{nee}_QC", "SW_IN_QC": "SW_IN_F_MDS_QC",
    "SWC_1_QC": "SWC_F_MDS_1_QC", "SWC_2_QC": "SWC_F_MDS_2_QC",
}
NIGHT_VARS = ["GPP_NT", "GPP_DT", "RECO_NT", "RECO_DT", "NEE"]
USTAR_PERCENTILES = ["16", "25", "50", "75", "84"]


def column(df, name, nee):
    """Half hourly series of a column, all NaN when the site lacks it."""
    col = name.format(nee=nee)
    if col in df:
        return df[col].astype(float)
    return pd.Series(np.nan, index=df.index)


def daily_mean(x, days):
    """Daily mean, NaN on days with any missing half hour."""
    g = x.groupby(days)
    return g.mean().where(g.count() == g.size())


def c2f_input(zip_path, nee="USTAR50"):
    """Daily input of the C2F QualityChecker from the zip of a shuttle site."""
    df = fs.read_subdaily(zip_path)
    days = df.index.floor("D")
    ds = pd.DataFrame(index=pd.DatetimeIndex(sorted(set(days)), name="time"))
    for name, col in MEAN_COLUMNS.items():
        ds[name] = daily_mean(column(df, col, nee), days)
    for name, col in QC_COLUMNS.items():
        qc = column(df, col, nee)
        ds[name] = daily_mean((qc <= 1).astype(float).where(qc.notna()), days)
    ds["P"] = column(df, "P_F", nee).groupby(days).sum(min_count=1)
    quant = pd.concat([column(df, f"NEE_VUT_{q}", nee).groupby(days).mean()
                       for q in USTAR_PERCENTILES], axis=1)
    ds["NEE_UstarUnc"] = quant.max(axis=1, skipna=False) - quant.min(axis=1,
                                                                      skipna=False)
    night = column(df, "SW_IN_POT", nee) == 0
    for name in NIGHT_VARS:
        x = column(df, MEAN_COLUMNS[name], nee).where(night)
        ds[name + "_NIGHT"] = x.groupby(days).mean()
    qc = column(df, QC_COLUMNS["NEE_QC"], nee)
    ds["NEE_QC_NIGHT"] = (qc <= 1).astype(float).where(qc.notna() & night).groupby(
        days).mean()

    out = xr.Dataset.from_dataframe(ds)
    out["noRain"] = np.logical_not(out["P"] > 0.5)
    out["posNETRAD"] = out["NETRAD"] > 0
    out["WUE_DT"] = out["GPP_DT"] * np.sqrt(out["VPD"])
    out["WUE_QC"] = np.fmin(out["NEE_QC"], out["VPD_QC"])
    out["WUE_NT"] = out["GPP_NT"] * np.sqrt(out["VPD"])
    out["LE_H"] = out["LE"] + out["H"]
    out["LE_H_QC"] = np.fmin(out["LE_QC"], out["H_QC"])
    return out


def flag_file(info, site):
    return os.path.join(os.path.dirname(info["cache_path"]),
                        f"{site}_jung_C2F_{info.get('strictness', 'loose')}.nc")


def flag_site(site, zip_path, info, qc_module):
    """C2F daily flags of one site, with one <VAR>_C2F_Flag per variable."""
    strict = info.get("strictness", "loose")
    params = getattr(qc_module.StrictnessParameters, strict)()
    ds = c2f_input(zip_path, info.get("nee_variant", "USTAR50")).assign_coords(
        site=site)
    np.random.seed(info.get("seed", 42))
    qc = qc_module.QualityChecker(strictness_parameters=params,
                                  njobsRF=info.get("n_jobs_rf", 8),
                                  log_file=flag_file(info, site) + ".log")
    qc.flag_multiple_sites([ds])
    flags = qc.get_qc_flags().astype(np.int8)
    if "site" in flags.dims:
        flags = flags.isel(site=0, drop=True)
    return flags.drop_vars([c for c in flags.coords if c != "time"])


class _NumpyWithIn1d(types.ModuleType):
    """numpy with in1d, which numpy 2.4 removed in favour of isin."""
    in1d = staticmethod(np.isin)

    def __getattr__(self, name):
        return getattr(np, name)


def patch_missingpy():
    """Let missingpy run with scikit-learn 1.8, which renamed force_all_finite of
    check_array to ensure_all_finite, and with numpy 2.4 or newer."""
    import missingpy.missforest as mf
    from sklearn.utils.validation import check_array

    def check(X, force_all_finite=True, **kwargs):
        return check_array(X, ensure_all_finite=force_all_finite, **kwargs)
    mf.check_array = check
    mf.np = _NumpyWithIn1d("numpy")


def load_c2f(code_dir):
    """The fluxnetqc package of c2f_fluxnet at code_dir, patched for numpy 2, which
    no longer sets an array element from an array of size one. calcWBvars gets P
    and ET as columns, so they are flattened first."""
    patch_missingpy()
    sys.path.insert(0, code_dir)
    import fluxnetqc
    import fluxnetqc.core as core
    calc = core.calcWBvars

    def calc_wb_vars(P, ET, *args, **kwargs):
        return calc(np.ravel(P), np.ravel(ET), *args, **kwargs)
    core.calcWBvars = calc_wb_vars
    return fluxnetqc


_QC_MODULE = None


def _flag_in_worker(site, zip_path, info):
    """flag_site in a worker process, which loads the C2F code once."""
    global _QC_MODULE
    if _QC_MODULE is None:
        _QC_MODULE = load_c2f(info["code_dir"])
    return flag_site(site, zip_path, info, _QC_MODULE)


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    """
    Flag the sites of the shuttle index that have no flags yet for their product
    file. Sites run in parallel in workers processes, n_workers of the json by
    default, each fitting its random forests with n_jobs_rf jobs. A site has only a
    few thousand days, so one job per forest and many sites at once is fastest. The
    seed is set per site and scikit-learn draws the seed of every tree before
    fitting, so the flags do not depend on the number of workers or jobs: DE-Hai,
    AR-Bal, BE-Vie and CH-Dav gave the same flags with 8, 2 and 1 jobs.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    shuttle = {"cache_path": info["shuttle_cache_path"]}
    done = fu.read_cache(info["cache_path"])
    done_files = set(done["file"].astype(str)) if len(done) else set()
    todo = []
    for site, lat, lon in zip(sites, lats, lons):
        rec = fs.site_record(shuttle, site)
        if rec is not None and str(rec["file"]) not in done_files:
            todo.append((site, lat, lon, rec))
    n = max(1, workers or info.get("n_workers", 1))
    print(f"{len(todo)} of {len(sites)} sites to flag with {n} workers", flush=True)
    folder = os.path.dirname(info["shuttle_cache_path"])
    with ProcessPoolExecutor(max_workers=n) as pool:
        jobs = {pool.submit(_flag_in_worker, site,
                            os.path.join(folder, str(rec["file"])), info):
                (site, lat, lon, rec) for site, lat, lon, rec in todo}
        for job in as_completed(jobs):
            site, lat, lon, rec = jobs[job]
            try:
                flags = job.result()
            except Exception as e:
                print(f"failed {site}: {type(e).__name__}: {e}", flush=True)
                continue
            write_flags(info, site, lat, lon, rec, flags)


def write_flags(info, site, lat, lon, rec, flags):
    """Write the flags of a site and add its row to the index at cache_path."""
    fetched = fu.now_utc()
    flags.attrs = {
        "site": site, "product": str(rec["file"]),
        "product_id": str(rec["product_id"]),
        "strictness": info.get("strictness", "loose"),
        "seed": info.get("seed", 42), "nee_variant": info.get("nee_variant",
                                                              "USTAR50"),
        "code": info.get("code_version", ""), "citation": info.get("citation", ""),
        "flag_meaning": "1 inconsistent, 0 consistent", "fetched_utc": fetched}
    path = flag_file(info, site)
    flags.to_netcdf(path + ".tmp")
    os.replace(path + ".tmp", path)
    counts = {f"n_{v.replace('_C2F_Flag', '')}": int(flags[v].sum())
              for v in flags.data_vars}
    fu.append_cache(info["cache_path"], [{
        "site": site, "file": str(rec["file"]),
        "flags": os.path.basename(path), "latitude": lat, "longitude": lon,
        "start": str(flags.time.values[0])[:10],
        "end": str(flags.time.values[-1])[:10], "n_days": flags.sizes["time"],
        **counts, "fetched_utc": fetched}])
    print(f"done {site}", flush=True)


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()
    days = date_.astype("datetime64[D]")
    path = flag_file(info, bxtr.site)
    flags = xr.open_dataset(path) if os.path.exists(path) else None
    if flags is None:
        logger.warning(f"::MISSING:: {dataset} has no C2F flags for {bxtr.site} at "
                       f"{path}. NaN values will be set.")

    src_data = []
    for tar_name in bxtr.vars_list:
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        bxtr.log_var_start(tar_name)
        values = np.full(len(date_), np.nan)
        if flags is not None and src_name in flags:
            f = flags[src_name].to_series()
            f.index = f.index.values.astype("datetime64[D]")
            values = 1.0 - f.reindex(days).values.astype(float)
        data = xr.Dataset({tar_name: xr.DataArray(data=values, dims=['time'],
                                                  coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = bxtr.merge_and_format(src_data)
    if flags is not None:
        attrs = {k: flags.attrs[k] for k in ("product", "product_id", "strictness",
                                             "seed", "code") if k in flags.attrs}
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
