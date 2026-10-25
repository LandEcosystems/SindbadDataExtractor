"""
Write the csv files the dashboard in dashboard/index.html reads, from a work folder
filled by run_all.sh: all_sites.csv, summary_by_variable.csv and status.csv of
aggregate.py, era/ of era_check.py and series/ of compare_site.py. The files go to
<work>/dashboard and are uploaded as the datasets named after them:

  variables.csv      agreement per resolution and variable
  site_variable.csv  agreement per resolution, variable and site
  site_status.csv    whether each site and resolution was compared
  era5_precip.csv    ERA5 precipitation at sites without a gauge
  hourly_lags.csv    hourly pairs that match best at a shift of whole hours
  series_index.csv   which monthly, daily and diurnal file holds each site
  monthly_<i>.csv    monthly means of the paired daily values
  daily_<i>.csv      paired daily values of the flux and meteo variables
  diurnal_<i>.csv    mean paired cycle over the day of each calendar month
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

from common import CACHE

FLUX = ["GPP_NT", "GPP_DT", "RECO_NT", "RECO_DT", "NEE", "LE", "H"]


def sig(df, cols, digits=4):
    """Round columns to significant digits to keep the files small."""
    fmt = f"%.{digits}g"
    for c in cols:
        df[c] = df[c].map(lambda v: float(fmt % v) if pd.notna(v) else v)
    return df


def kind_of(variable):
    if "_QC" in variable:
        return "QC"
    if "RANDUNC" in variable:
        return "Uncertainty"
    return "Flux" if variable in FLUX else "Meteo"


def chunked(out, prefix, frames, limit_mb):
    """Write site tables into files of at most about limit_mb, in site order."""
    chunks, part, size = {}, [], 0
    def flush():
        if part:
            name = f"{prefix}_{len(set(chunks.values())) + 1}"
            pd.concat(part, ignore_index=True).to_csv(f"{out}/{name}.csv",
                                                      index=False, float_format="%.4g")
            for d in part:
                chunks[d.site.iloc[0]] = name
    for d in frames:
        n = len(d.to_csv(index=False, float_format="%.4g"))
        if part and size + n > limit_mb * 2 ** 20:
            flush()
            part, size = [], 0
        part.append(d)
        size += n
    flush()
    return chunks


def write_series(work, out, ok, limit_mb=7):
    """
    Monthly, daily and diurnal series of every site from series/ of compare_site.py.
    Monthly means use months with at least 10 paired days and keep every variable.
    The daily files keep the flux and meteo variables only, to stay small enough
    to load on a click. Each kind is split into files of at most limit_mb.
    """
    keep = set(ok.loc[ok.kind.isin(["Flux", "Meteo"]), "variable"])
    monthly, daily, cyc, rows = [], [], [], {}
    for f in sorted(glob.glob(f"{work}/series/daily/*.csv")):
        site = os.path.basename(f)[:-4]
        d = pd.read_csv(f, index_col="date", parse_dates=True)
        m = d.groupby(d.index.to_period("M"))
        mm = m.mean().where(m.count() >= 10).dropna(how="all")
        mm.index = mm.index.astype(str)
        mm = mm.rename_axis("month").reset_index()
        mm.insert(0, "site", site)
        mm.insert(0, "id", site + "|" + mm.month)
        monthly.append(mm)
        cols = [c for c in d.columns if c.rsplit("__", 1)[0] in keep]
        dd = d[cols].dropna(how="all")
        dd.index = dd.index.strftime("%Y-%m-%d")
        dd = dd.rename_axis("date").reset_index()
        dd.insert(0, "site", site)
        dd.insert(0, "id", site + "|" + dd.date)
        daily.append(dd)
        rows[site] = {"site": site, "first_month": mm.month.min(),
                      "last_month": mm.month.max(), "n_months": len(mm),
                      "n_days": len(dd)}
    for f in sorted(glob.glob(f"{work}/series/diurnal/*.csv")):
        site = os.path.basename(f)[:-4]
        c = pd.read_csv(f)
        c.insert(0, "site", site)
        c.insert(0, "id", site + "|" + c.month.astype(str) + "|" + c.hour.astype(str))
        cyc.append(c)
        rows.setdefault(site, {"site": site})
    for kind, frames in [("monthly", monthly), ("daily", daily), ("diurnal", cyc)]:
        for site, name in chunked(out, kind, frames, limit_mb).items():
            rows[site][f"{kind}_file"] = name
    pd.DataFrame(list(rows.values())).sort_values("site").to_csv(
        f"{out}/series_index.csv", index=False)


def main(work):
    out = os.path.join(work, "dashboard")
    os.makedirs(out, exist_ok=True)
    a = pd.read_csv(f"{work}/all_sites.csv")
    ok = a[(a.status == "ok") & (a.n >= 30)].copy()
    ok["kind"] = ok.variable.map(kind_of)

    ok["id"] = ok.reso + "|" + ok.variable + "|" + ok.site
    cols = ["id", "site", "reso", "variable", "kind", "years", "n", "n_ref_only",
            "n_shuttle_only", "mean_ref", "mean_shuttle", "bias", "rel_bias", "rmse", "r",
            "slope", "nse", "frac_equal", "best_lag_h", "units", "source_column"]
    sv = sig(ok[cols].copy(), ["mean_ref", "mean_shuttle", "bias", "rmse"])
    for c in ["rel_bias", "r", "slope", "nse", "frac_equal"]:
        sv[c] = sv[c].round(4)
    sv.to_csv(f"{out}/site_variable.csv", index=False)

    s = pd.read_csv(f"{work}/summary_by_variable.csv")
    s = s.merge(ok.groupby(["reso", "variable"]).kind.first().reset_index(),
                on=["reso", "variable"])
    s["id"] = s.reso + "|" + s.variable
    s = s.rename(columns={"sites_r_lt_0.9": "sites_r_below_0_9"})
    keep = ["id", "reso", "variable", "kind", "n_sites", "n_steps", "within_r_median",
            "within_r_p10", "within_r_min", "nse_median", "nse_p10", "slope_median", "rel_bias_median",
            "abs_rel_bias_p90", "frac_equal_median", "sites_r_below_0_9", "across_r",
            "across_bias", "across_rmse", "across_mean_ref", "across_mean_shuttle",
            "worst_sites"]
    s = s[keep]
    sig(s, [c for c in keep if s[c].dtype == float]).to_csv(f"{out}/variables.csv",
                                                         index=False)

    st = pd.read_csv(f"{work}/status.csv")
    idx = pd.read_csv(CACHE).drop_duplicates("site", keep="last")
    years = a.groupby(["site", "reso"]).years.first().reset_index()
    st = st.merge(years, on=["site", "reso"], how="left").merge(
        idx[["site", "data_hub", "first_year", "last_year", "igbp"]], on="site",
        how="left")
    st["shuttle_years"] = st.first_year.astype(str) + "-" + st.last_year.astype(str)
    st["id"] = st.reso + "|" + st.site
    st = st.rename(columns={"years": "common_years"})
    st[["id", "site", "reso", "status", "common_years", "shuttle_years", "data_hub",
        "igbp"]].to_csv(f"{out}/site_status.csv", index=False)

    e = pd.concat([pd.read_csv(f) for f in glob.glob(f"{work}/era/*.csv")])
    p = a[(a.reso == "daily") & (a.variable == "P") & (a.status == "ok")][
        ["site", "mean_ref", "mean_shuttle"]]
    e = e[e.P_mm_per_year.isna()].merge(p, on="site", how="left")
    e["fluxnet2015_P_mm_per_year"] = (e.mean_ref * 365.25).round(0)
    e["shuttle_P_mm_per_year_common"] = (e.mean_shuttle * 365.25).round(0)
    e["ratio"] = (e.mean_shuttle / e.mean_ref).round(2)
    e["P_ERA_mm_per_year"] = e.P_ERA_mm_per_year.round(0)
    e[["site", "step_h", "P_ERA_mm_per_year", "fluxnet2015_P_mm_per_year",
       "shuttle_P_mm_per_year_common", "ratio"]].sort_values(
        "ratio", ascending=False).to_csv(f"{out}/era5_precip.csv", index=False)

    lag = ok[(ok.reso == "hourly") & ok.best_lag_h.notna() &
             (ok.best_lag_h != 0)].copy()
    lag["id"] = lag.site + "|" + lag.variable
    lag[["id", "site", "variable", "best_lag_h", "r_lag0", "r"]].round(3).to_csv(
        f"{out}/hourly_lags.csv", index=False)

    write_series(work, out, ok)
    for f in sorted(glob.glob(f"{out}/*.csv")):
        print(f"{os.path.basename(f)}: {os.path.getsize(f) // 1024} KB")


if __name__ == "__main__":
    main(sys.argv[1])
