#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESA CCI Land Cover PFT v2.0.81 cover fractions within 1 km of the site.

Each class is the mean of the 300 m grid cells whose centres lie within 1 km of the
site, read through the CEDA OPeNDAP server for every year in years. The cache holds
one row per site and year with the mean percentage of every class, the same format as
esacci_lc_pft_v2.0.81_1km_sites.csv of add_vegfrac_esacci.py.

Tree is the four TREES classes, shrub the four SHRUBS classes and grass the natural
and managed grass. Vegetation is tree + shrub + grass.
"""
import re
import sys

import numpy as np

import extractors.footprint_utils as fu

URL = ("https://dap.ceda.ac.uk/thredds/dodsC/neodc/esacci/land_cover/data/pft/v2.0.81/"
       "ESACCI-LC-L4-PFT-Map-300m-P1Y-{year}-v2.0.81.nc.ascii?")
CLASSES = ["TREES-BE", "TREES-BD", "TREES-NE", "TREES-ND", "SHRUBS-BE", "SHRUBS-BD",
           "SHRUBS-NE", "SHRUBS-ND", "GRASS-NAT", "GRASS-MAN", "BARE", "BUILT",
           "WATER", "SNOWICE"]
GROUPS = {
    "tree": [c for c in CLASSES if c.startswith("TREES")],
    "shrub": [c for c in CLASSES if c.startswith("SHRUBS")],
    "grass": ["GRASS-NAT", "GRASS-MAN"],
    "bare": ["BARE"],
}
GROUPS["veg"] = GROUPS["tree"] + GROUPS["shrub"] + GROUPS["grass"]
PX_PER_DEG = 360


def parse_ascii(text, shape):
    """Return {variable: array} from a DAP ascii response of 3D grids."""
    out, var, rows = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^([A-Z-]+)\.\1\[", line)
        if m:
            if var is not None:
                out[var] = np.array(rows, dtype=float).reshape(shape)
            var, rows = m.group(1), []
        elif var is not None and line.startswith("[0]["):
            rows.extend(float(v) for v in line.split(",")[1:])
    if var is not None:
        out[var] = np.array(rows, dtype=float).reshape(shape)
    return out


def site_year(site, lat, lon, year):
    half_lat, half_lon = fu.half_box(lat)
    row0 = int(np.floor((90 - (lat + half_lat)) * PX_PER_DEG))
    row1 = int(np.ceil((90 - (lat - half_lat)) * PX_PER_DEG))
    col0 = int(np.floor((lon - half_lon + 180) * PX_PER_DEG))
    col1 = int(np.ceil((lon + half_lon + 180) * PX_PER_DEG))
    sub = f"[0][{row0}:1:{row1 - 1}][{col0}:1:{col1 - 1}]"
    url = URL.format(year=year) + ",".join(c + sub for c in CLASSES)
    text = fu.fetch_url(url, retries=5, wait=2, timeout=120).decode()
    data = parse_ascii(text, (row1 - row0, col1 - col0))
    lats = 90 - (np.arange(row0, row1) + 0.5) / PX_PER_DEG
    lons = -180 + (np.arange(col0, col1) + 0.5) / PX_PER_DEG
    lon2d, lat2d = np.meshgrid(lons, lats)
    inside = fu.haversine(lat, lon, lat2d, lon2d) <= fu.RADIUS_M
    rec = {"site": site, "year": year, "n_grids": int(inside.sum())}
    for c in CLASSES:
        v = data[c][inside]
        v[(v < 0) | (v > 100)] = np.nan
        rec[c] = np.nanmean(v) if np.isfinite(v).any() else np.nan
    return rec


def years_of(info):
    y0, y1 = info.get("years", [1992, 2017])
    return list(range(int(y0), int(y1) + 1))


def fetch_site(site, lat, lon, info):
    return [site_year(site, lat, lon, y) for y in years_of(info)]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    rows = rows.drop_duplicates("year", keep="last")
    rows = rows[rows["year"].isin(years_of(info))]
    out = {}
    for g, cols in GROUPS.items():
        vals = rows[cols].sum(axis=1, min_count=1) / 100
        out[g] = dict(zip(rows["year"].astype(int), vals.astype(float)))
    return out


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
