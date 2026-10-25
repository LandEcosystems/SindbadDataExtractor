"""Compare the shuttle extraction of one site with the FLUXNET2015 daily and hourly
SINDBAD files, per variable, on the paired time steps of the common period."""
import json, sys, os, logging, warnings
import numpy as np, pandas as pd, xarray as xr

from common import CACHE, REF, REPO
logging.disable(logging.WARNING)
warnings.filterwarnings("ignore")
from utils.setup_config import get_variable_info
import extractors.fluxnet_shuttle as fs

LAG_VARS = ["SW_IN", "TA", "NEE", "LE", "H", "GPP_NT"]


def shuttle(reso, site, lat, lon, start, end):
    folder = {"daily": "sindbad", "hourly": "hourly"}[reso]
    cfg = {"temporal_resolution": reso, "start_date": start, "end_date": end,
           "fluxcom_cube_path": "", "FLUXNET_version": "shuttle", "dataset": {},
           "sel_datasets": ["shuttle_obs", "shuttle_clim"]}
    for k in ["obs", "clim"]:
        d = json.load(open(f"{REPO}/configuration/{folder}/fluxnet_shuttle_{k}.json"))
        d["cache_path"] = CACHE
        cfg["dataset"][f"shuttle_{k}"] = d
    _, cfg = get_variable_info(cfg)
    out = [fs.extract(k, {"site_ID": site, "latitude": lat, "longitude": lon}, cfg)
           for k in ["shuttle_obs", "shuttle_clim"]]
    return xr.merge([o for o in out if o is not None]).squeeze(drop=True)


def stats(a, b):
    """Metrics of shuttle b against reference a on paired finite values."""
    ok = np.isfinite(a) & np.isfinite(b)
    n = int(ok.sum())
    if n < 3:
        return {"n": n}
    a, b = a[ok], b[ok]
    d = b - a
    sa, sb = a.std(), b.std()
    r = np.corrcoef(a, b)[0, 1] if sa > 0 and sb > 0 else np.nan
    slope = np.cov(a, b)[0, 1] / a.var() if sa > 0 else np.nan
    tol = 1e-4 * max(np.abs(a).max(), 1e-9)
    # Nash-Sutcliffe efficiency of the shuttle with FLUXNET2015 as reference
    nse = 1 - (d ** 2).sum() / ((a - a.mean()) ** 2).sum() if sa > 0 else np.nan
    return {"n": n, "mean_ref": a.mean(), "mean_shuttle": b.mean(), "bias": d.mean(),
            "rel_bias": d.mean() / abs(a.mean()) if a.mean() != 0 else np.nan,
            "mae": np.abs(d).mean(), "rmse": np.sqrt((d ** 2).mean()),
            "std_ref": sa, "std_shuttle": sb, "r": r, "slope": slope, "nse": nse,
            "frac_equal": float((np.abs(d) <= tol).mean()),
            "max_abs_diff": np.abs(d).max()}


def qc_parent(name):
    """TA_QC_DayMin -> TA_DayMin, NEE_QC -> NEE, PPFD_IN_QC -> PPFD_IN."""
    if "_QC" not in name:
        return None
    return name.replace("_QC", "", 1)


def best_lag(a, b, lags=range(-3, 4)):
    rs = {}
    for k in lags:
        bb = np.roll(b, k)
        ok = np.isfinite(a) & np.isfinite(bb)
        if ok.sum() > 100:
            rs[k] = np.corrcoef(a[ok], bb[ok])[0, 1]
    return max(rs, key=rs.get) if rs else np.nan, rs.get(0, np.nan)


def compare(site, reso, rec, outdir):
    path = REF[reso].format(site=site)
    if not os.path.exists(path):
        return []
    ref = xr.open_dataset(path).squeeze(drop=True)
    lat, lon = float(ref.attrs["latitude"]), float(ref.attrs["longitude"])
    # reference period with any flux or met data
    core = [v for v in ["NEE", "GPP_NT", "LE", "H", "TA", "SW_IN"] if v in ref]
    valid = np.zeros(ref.sizes["time"], bool)
    for v in core:
        valid |= np.isfinite(ref[v].values)
    if not valid.any():
        return [{"site": site, "reso": reso, "status": "no reference data"}]
    t = ref.time.values[valid]
    y0 = max(int(str(t[0])[:4]), int(rec["first_year"]))
    y1 = min(int(str(t[-1])[:4]), int(rec["last_year"]))
    if y0 > y1:
        return [{"site": site, "reso": reso, "status": "no overlap",
                 "ref_years": f"{str(t[0])[:4]}-{str(t[-1])[:4]}",
                 "shuttle_years": f"{rec['first_year']}-{rec['last_year']}"}]
    sh = shuttle(reso, site, lat, lon, f"{y0}-01-01", f"{y1}-12-31")
    ref = ref.sel(time=slice(sh.time.values[0], sh.time.values[-1]))
    sh = sh.sel(time=ref.time)
    rows, paired = [], {}
    for v in sh.data_vars:
        if not v.endswith("_shuttle"):
            continue
        base = v[:-len("_shuttle")]
        if base not in ref or ref[base].dims != ("time",):
            continue
        a = ref[base].values.astype(float)
        b = sh[v].values.astype(float)
        n_ref_only = int((np.isfinite(a) & ~np.isfinite(b)).sum())
        n_shuttle_only = int((~np.isfinite(a) & np.isfinite(b)).sum())
        parent = qc_parent(base)
        if parent and parent in ref and f"{parent}_shuttle" in sh:
            # QC only where the variable it belongs to exists in both
            has = (np.isfinite(ref[parent].values) &
                   np.isfinite(sh[f"{parent}_shuttle"].values))
            a = np.where(has, a, np.nan)
            b = np.where(has, b, np.nan)
        ok = np.isfinite(a) & np.isfinite(b)
        row = {"site": site, "reso": reso, "variable": base, "status": "ok",
               "years": f"{y0}-{y1}",
               "n_ref": int(np.isfinite(a).sum()), "n_shuttle": int(np.isfinite(b).sum()),
               "first_paired": str(ref.time.values[ok][0])[:16] if ok.any() else "",
               "last_paired": str(ref.time.values[ok][-1])[:16] if ok.any() else "",
               "n_ref_only": n_ref_only, "n_shuttle_only": n_shuttle_only,
               "qc_masked_by": parent or "",
               "units": ref[base].attrs.get("units", ""),
               "source_column": sh[v].attrs.get("source_column", ""),
               **stats(a, b)}
        if reso == "hourly" and base in LAG_VARS:
            row["best_lag_h"], row["r_lag0"] = best_lag(a, b)
        rows.append(row)
        if ok.sum() >= 30:
            paired[f"{base}__ref"] = np.where(ok, a, np.nan)
            paired[f"{base}__shuttle"] = np.where(ok, b, np.nan)
    pd.DataFrame(rows).to_csv(f"{outdir}/{site}_{reso}.csv", index=False)
    write_series(site, reso, pd.DataFrame(paired, index=pd.DatetimeIndex(
        ref.time.values)), outdir)
    return rows


def write_series(site, reso, df, outdir):
    """
    Save the paired values of a site next to outdir: the daily series as they are in
    series/daily, and from the hourly ones the mean cycle over the day of each
    calendar month in series/diurnal, for hours with at least 10 paired values.
    """
    if df.empty:
        return
    base = os.path.join(os.path.dirname(os.path.abspath(outdir)), "series")
    if reso == "daily":
        os.makedirs(f"{base}/daily", exist_ok=True)
        df = df.dropna(how="all")
        df.index = df.index.strftime("%Y-%m-%d")
        df.index.name = "date"
        df.to_csv(f"{base}/daily/{site}.csv", float_format="%.4g")
        return
    os.makedirs(f"{base}/diurnal", exist_ok=True)
    g = df.groupby([df.index.month.rename("month"), df.index.hour.rename("hour")])
    cyc = g.mean().where(g.count() >= 10)
    cyc.dropna(how="all").to_csv(f"{base}/diurnal/{site}.csv", float_format="%.4g")


def main(site, outdir):
    os.makedirs(outdir, exist_ok=True)
    idx = pd.read_csv(CACHE)
    rec = idx[idx.site == site].iloc[-1]
    for reso in ["daily", "hourly"]:
        if os.path.exists(f"{outdir}/{site}_{reso}.csv"):
            continue
        try:
            rows = compare(site, reso, rec, outdir)
        except Exception as e:
            rows = [{"site": site, "reso": reso, "status": f"error {type(e).__name__}: {e}"}]
        if rows and rows[0].get("status") != "ok":
            pd.DataFrame(rows).to_csv(f"{outdir}/{site}_{reso}.csv", index=False)
    print("done", site, flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
