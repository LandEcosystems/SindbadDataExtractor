#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MODIS Vegetation Continuous Fields MOD44B v061 fractions within 1 km of the site.

Each fraction is the mean of the 250 m MODIS pixels whose centres lie within 1 km of
the site. The pixel centres are computed on the MODIS sinusoidal grid and sent to
NASA AppEEARS as a point request, which returns the value of exactly those pixels.
The Earthdata login comes from ~/.netrc. AppEEARS requests take hours, so the data
cannot be fetched inside extract, and prefetch runs in three steps:
    python prefetch.py <main json> vf_modis --step submit
    python prefetch.py <main json> vf_modis --step status
    python prefetch.py <main json> vf_modis --step process
The request and the downloaded results are kept in work_dir. The cache holds one row
per site, the same format as modis_vcf_2010/site_means.csv of add_vegfrac_modis.py.

Water pixels (200) count as non-vegetated and fill pixels (253) are dropped. North of
60 N, MOD44B v061 reports 0 % tree, 0 % non-tree vegetation and 100 % non-vegetated
over whole tiles of boreal forest and tundra, so sites with such pixels get NaN.
Vegetation is tree + non-tree vegetation.
"""
import json
import netrc
import os
import sys
import urllib.request
from base64 import b64encode

import numpy as np
import pandas as pd

import extractors.footprint_utils as fu

API = "https://appeears.earthdatacloud.nasa.gov/api"
LAYERS = ["Percent_Tree_Cover", "Percent_NonTree_Vegetation", "Percent_NonVegetated"]
POINTS_PER_TASK = 1000

# MODIS sinusoidal grid for the 250 m products.
R_SIN = 6371007.181
PX = 231.65635826395834
X0, Y0 = -20015109.354, 10007554.677
WATER, FILL = 200, 253
BAD_NORTH_LAT = 60.0


def year_of(info):
    return int(info.get("years", [2010])[0])


def pixel_centres(site, lat, lon):
    x = R_SIN * np.radians(lon) * np.cos(np.radians(lat))
    y = R_SIN * np.radians(lat)
    ci, ri = int((x - X0) // PX), int((Y0 - y) // PX)
    n = int(np.ceil(fu.RADIUS_M / PX)) + 2
    cols, rows = np.meshgrid(np.arange(ci - n, ci + n + 1),
                             np.arange(ri - n, ri + n + 1))
    px = X0 + (cols + 0.5) * PX
    py = Y0 - (rows + 0.5) * PX
    plat = np.degrees(py / R_SIN)
    plon = np.degrees(px / (R_SIN * np.cos(np.radians(plat))))
    inside = fu.haversine(lat, lon, plat, plon) <= fu.RADIUS_M
    return pd.DataFrame({"site": site, "row": rows[inside], "col": cols[inside],
                         "latitude": plat[inside], "longitude": plon[inside],
                         "site_latitude": lat, "site_longitude": lon})


def request(path, token=None, data=None, method=None, auth=None):
    headers = {"Content-Type": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    if auth:
        headers["Authorization"] = f"Basic {auth}"
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(f"{API}/{path}", body, headers, method=method)
    with urllib.request.urlopen(req, timeout=300) as r:
        raw = r.read()
    return json.loads(raw) if raw[:1] in (b"{", b"[") else raw


def login():
    user, _, password = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
    auth = b64encode(f"{user}:{password}".encode()).decode()
    return request("login", auth=auth, method="POST", data={})["token"]


def paths(info):
    work = info["work_dir"]
    return os.path.join(work, "tasks.json"), os.path.join(work, "points.csv")


def submit(sites, lats, lons, info):
    tasks_path, points_path = paths(info)
    if os.path.exists(tasks_path):
        print(f"a request is already in {tasks_path}. Run the process step, or remove "
              "the work_dir files to submit a new one.")
        return
    sites, lats, lons = fu.missing_sites(info, sites, lats, lons)
    if not sites:
        print(f"all sites are in {info['cache_path']}")
        return
    year = year_of(info)
    pts = pd.concat([pixel_centres(s, la, lo) for s, la, lo in zip(sites, lats, lons)],
                    ignore_index=True)
    pts["id"] = [f"p{i}" for i in range(len(pts))]
    os.makedirs(info["work_dir"], exist_ok=True)
    pts.to_csv(points_path, index=False)
    token = login()
    tasks = []
    for k in range(0, len(pts), POINTS_PER_TASK):
        chunk = pts.iloc[k:k + POINTS_PER_TASK]
        task = {
            "task_type": "point",
            "task_name": f"fluxnet_vcf_{year}_{k // POINTS_PER_TASK:02d}",
            "params": {
                "dates": [{"startDate": f"01-01-{year}", "endDate": f"12-31-{year}"}],
                "layers": [{"product": "MOD44B.061", "layer": la} for la in LAYERS],
                "coordinates": [{"id": r.id, "category": r.site,
                                 "latitude": round(r.latitude, 7),
                                 "longitude": round(r.longitude, 7)}
                                for r in chunk.itertuples()],
            },
        }
        tasks.append(request("task", token=token, data=task, method="POST")["task_id"])
    with open(tasks_path, "w") as f:
        json.dump(tasks, f)
    print(f"{len(pts)} pixel centres for {len(sites)} sites submitted as "
          f"{len(tasks)} tasks")


def status(info):
    tasks_path, _ = paths(info)
    token = login()
    out = {t: request(f"task/{t}", token=token)["status"]
           for t in json.load(open(tasks_path))}
    print(pd.Series(out).value_counts().to_dict())
    return out


def download(info, token):
    tasks_path, _ = paths(info)
    frames = []
    for t in json.load(open(tasks_path)):
        dest = os.path.join(info["work_dir"], f"{t}.csv")
        if not os.path.exists(dest):
            files = request(f"bundle/{t}", token=token)["files"]
            f = next(f for f in files if f["file_name"].endswith("results.csv"))
            with open(dest, "wb") as fh:
                fh.write(request(f"bundle/{t}/{f['file_id']}", token=token))
        frames.append(pd.read_csv(dest))
    return pd.concat(frames, ignore_index=True)


def site_means(res, pts):
    """Per-site means of the AppEEARS pixel values, with the masking rules above."""
    cols = {la: next(c for c in res.columns if c.endswith(la)) for la in LAYERS}
    v = res[["ID"] + list(cols.values())].rename(columns={c: la for la, c in
                                                          cols.items()})
    v = v.groupby("ID").first().join(pts[["site", "latitude"]])
    for la in LAYERS:
        v[la] = v[la].astype(float)
        v.loc[v[la] == FILL, la] = np.nan
    all_bare = ((v["Percent_Tree_Cover"] == 0)
                & (v["Percent_NonTree_Vegetation"] == 0)
                & (v["Percent_NonVegetated"] == 100))
    bad_north = all_bare & (v["latitude"] > BAD_NORTH_LAT)
    # The artefact covers whole tiles, so the water pixels left at affected sites
    # would otherwise make them look unvegetated.
    bad_sites = v.loc[bad_north, "site"].unique()
    v.loc[v["site"].isin(bad_sites), LAYERS] = np.nan
    water = v["Percent_Tree_Cover"] == WATER
    v.loc[water, "Percent_Tree_Cover"] = 0
    v.loc[water, "Percent_NonTree_Vegetation"] = 0
    v.loc[water, "Percent_NonVegetated"] = 100
    v["veg"] = v["Percent_Tree_Cover"] + v["Percent_NonTree_Vegetation"]
    site = v.groupby("site")[LAYERS + ["veg"]].mean() / 100
    site["n_pixels"] = v.groupby("site").size()
    site["n_water"] = water.groupby(v["site"]).sum()
    site["n_bad_north"] = bad_north.groupby(v["site"]).sum()
    loc = pts.groupby("site")[["site_latitude", "site_longitude"]].first()
    site["latitude"] = loc["site_latitude"]
    site["longitude"] = loc["site_longitude"]
    return site.reset_index()


def process(info):
    if set(status(info).values()) != {"done"}:
        print("not all tasks are done yet")
        return
    res = download(info, login())
    _, points_path = paths(info)
    pts = pd.read_csv(points_path).set_index("id")
    site = site_means(res, pts)
    done = {s.lower() for s in fu.cached_sites(info["cache_path"])}
    site = site[~site["site"].str.lower().isin(done)]
    fu.append_cache(info["cache_path"], site.to_dict("records"))
    print(f"{len(site)} sites added to {info['cache_path']}")


def fetch_site(site, lat, lon, info):
    raise RuntimeError("MODIS VCF comes from asynchronous AppEEARS requests. Run "
                       "prefetch.py with --step submit, status and process.")


def prefetch(sites, lats, lons, info, workers=None, step="process", **kwargs):
    if step == "submit":
        submit(sites, lats, lons, info)
    elif step == "status":
        status(info)
    else:
        process(info)


def components(rows, info):
    row = rows.iloc[-1]
    year = year_of(info)
    return {"tree": {year: float(row["Percent_Tree_Cover"])},
            "nontree": {year: float(row["Percent_NonTree_Vegetation"])},
            "bare": {year: float(row["Percent_NonVegetated"])},
            "veg": {year: float(row["veg"])}}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
