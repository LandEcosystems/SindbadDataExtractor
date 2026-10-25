"""Pool the per site comparisons into per variable summaries and an across site
comparison of the site means of the paired time steps. The argument is the work
folder whose out/ holds the per site csv files of compare_site.py."""
import glob, sys
import numpy as np, pandas as pd

S = sys.argv[1]
d = pd.concat([pd.read_csv(f) for f in glob.glob(f"{S}/out/*.csv")], ignore_index=True)
d.to_csv(f"{S}/all_sites.csv", index=False)
status = d.groupby(["site", "reso"]).status.first().reset_index()
status.to_csv(f"{S}/status.csv", index=False)
ok = d[(d.status == "ok") & (d.n >= 30)].copy()


def across(g):
    a, b = g.mean_ref.values, g.mean_shuttle.values
    r = np.corrcoef(a, b)[0, 1] if len(g) > 2 and a.std() > 0 else np.nan
    return pd.Series({
        "n_sites": len(g), "n_steps": int(g.n.sum()),
        "within_r_median": g.r.median(), "within_r_p10": g.r.quantile(0.1),
        "within_r_min": g.r.min(), "nse_median": g.nse.median(),
        "nse_p10": g.nse.quantile(0.1), "slope_median": g.slope.median(),
        "rel_bias_median": g.rel_bias.median(),
        "abs_rel_bias_p90": g.rel_bias.abs().quantile(0.9),
        "frac_equal_median": g.frac_equal.median(),
        "sites_r_lt_0.9": int((g.r < 0.9).sum()),
        "across_mean_ref": a.mean(), "across_mean_shuttle": b.mean(),
        "across_r": r, "across_bias": (b - a).mean(),
        "across_rmse": np.sqrt(((b - a) ** 2).mean()),
        "worst_sites": ", ".join(f"{s} ({x:.2f})" for s, x in
                                 g.nsmallest(5, "r")[["site", "r"]].values),
    })


summ = ok.groupby(["reso", "variable"]).apply(across, include_groups=False).reset_index()
summ.to_csv(f"{S}/summary_by_variable.csv", index=False)

site = ok.groupby(["reso", "site"]).agg(
    years=("years", "first"), n_vars=("variable", "size"),
    r_median=("r", "median"), r_min=("r", "min"), nse_median=("nse", "median"),
    worst_var=("r", lambda x: ok.loc[x.idxmin(), "variable"])).reset_index()
site.to_csv(f"{S}/summary_by_site.csv", index=False)

if "best_lag_h" in ok:
    lag = ok[ok.reso == "hourly"].dropna(subset=["best_lag_h"])
    lag.groupby("variable").best_lag_h.value_counts().unstack(fill_value=0).to_csv(
        f"{S}/hourly_lags.csv")
print(status.status.value_counts().to_string())
print(len(ok.site.unique()), "sites compared")
