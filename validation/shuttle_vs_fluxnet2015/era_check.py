"""ERA5 checks of the shuttle zip of one site: annual ERA5 precipitation against the
measured one where the site has a gauge, and the range of ERA5 VPD, TA and SW_IN
against the measured series."""
import os, sys, glob
import numpy as np, pandas as pd

from common import CACHE
import extractors.fluxnet_shuttle as fs

Z = os.path.dirname(CACHE) + "/"
site, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
f = glob.glob(Z + f"*_{site}_FLUXNET_*.zip")[0]
df = fs.read_subdaily(f)
step_h = (df.index[1] - df.index[0]) / pd.Timedelta("1h")
both = df.P.notna() & df.P_ERA.notna() if "P" in df else pd.Series(False, df.index)
years = len(df) * step_h / 8766
row = {"site": site, "file": f.split("/")[-1], "step_h": step_h,
       "P_ERA_mm_per_year": df.P_ERA.sum() / years,
       "P_mm_per_year": df.P[both].sum() / (both.sum() * step_h / 8766) if both.any() else np.nan,
       "P_ERA_over_P": df.P_ERA[both].sum() / df.P[both].sum() if both.any() and df.P[both].sum() > 0 else np.nan,
       "frac_P_F_from_ERA": float((df.P_F_QC == 2).mean()) if "P_F_QC" in df else np.nan,
       "VPD_ERA_max_hPa": df.VPD_ERA.max(), "VPD_ERA_mean_hPa": df.VPD_ERA.mean(),
       "VPD_mean_hPa": df.VPD_F_MDS.mean() if "VPD_F_MDS" in df else np.nan,
       "TA_ERA_mean": df.TA_ERA.mean(), "TA_mean": df.TA_F_MDS.mean() if "TA_F_MDS" in df else np.nan,
       "SW_IN_ERA_mean": df.SW_IN_ERA.mean(), "SW_IN_mean": df.SW_IN_F_MDS.mean() if "SW_IN_F_MDS" in df else np.nan}
pd.DataFrame([row]).to_csv(f"{out}/{site}.csv", index=False)
