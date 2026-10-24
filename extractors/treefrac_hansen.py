#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tree cover fraction within 1 km of the site from Hansen Global Forest Change v1.11
at 30 m, for every year from 2000 to 2023.

The cover of a year is the treecover2000 percentage of the grid cells whose centres
lie within 1 km of the site, set to 0 where lossyear marks a loss in or before that
year, and averaged over the grid cells that the datamask layer maps as land or water.
Gains are not added, since the gain layer only flags a gain at some time between
2000 and 2012. With time_placement annual_extend, years before 2000 take the 2000
cover and years after 2023 the 2023 cover.

The 10 degree tiles are read from source_dir. Sites outside the tile coverage from
80N to 60S get NaN. The source is local, so cache_path may be empty.
"""
import os
import sys

import numpy as np

import extractors.footprint_utils as fu

TILE = "{dir}/{layer}/Hansen_GFC-2023-v1.11_{layer}_{tile}.tif"
YEARS = list(range(2000, 2024))
PX_PER_DEG = 4000
TILE_PX = 40000


def years_of(info):
    return [int(y) for y in info.get("years", YEARS)]


def tile_name(top, left):
    ns = f"{top:02d}N" if top >= 0 else f"{-top:02d}S"
    ew = f"{left:03d}E" if left >= 0 else f"{-left:03d}W"
    return f"{ns}_{ew}"


def read_window(src_dir, layer, row0, row1, col0, col1):
    """
    Read global grid cell rows and cols, which may span several 10 degree tiles.
    """
    import rasterio
    from rasterio.windows import Window
    out = np.full((row1 - row0, col1 - col0), np.nan)
    for tr in range(row0 // TILE_PX, (row1 - 1) // TILE_PX + 1):
        for tc in range(col0 // TILE_PX, (col1 - 1) // TILE_PX + 1):
            r0, r1 = max(row0, tr * TILE_PX), min(row1, (tr + 1) * TILE_PX)
            c0, c1 = max(col0, tc * TILE_PX), min(col1, (tc + 1) * TILE_PX)
            path = TILE.format(dir=src_dir, layer=layer,
                               tile=tile_name(90 - 10 * tr, -180 + 10 * tc))
            if not os.path.exists(path):
                continue
            with rasterio.open(path) as src:
                win = Window(c0 - tc * TILE_PX, r0 - tr * TILE_PX, c1 - c0, r1 - r0)
                out[r0 - row0:r1 - row0, c0 - col0:c1 - col0] = src.read(1, window=win)
    return out


def fetch_site(site, lat, lon, info):
    half_lat, half_lon = fu.half_box(lat)
    row0 = int(np.floor((90 - (lat + half_lat)) * PX_PER_DEG))
    row1 = int(np.ceil((90 - (lat - half_lat)) * PX_PER_DEG))
    col0 = int(np.floor((lon - half_lon + 180) * PX_PER_DEG))
    col1 = int(np.ceil((lon + half_lon + 180) * PX_PER_DEG))
    cover, loss, mask = (read_window(info["source_dir"], layer, row0, row1, col0, col1)
                         for layer in ("treecover2000", "lossyear", "datamask"))
    lats = 90 - (np.arange(row0, row1) + 0.5) / PX_PER_DEG
    lons = -180 + (np.arange(col0, col1) + 0.5) / PX_PER_DEG
    lon2d, lat2d = np.meshgrid(lons, lats)
    inside = fu.haversine(lat, lon, lat2d, lon2d) <= fu.RADIUS_M
    valid = inside & np.isin(mask, (1, 2))
    rows = []
    for year in years_of(info):
        lost = (loss >= 1) & (loss <= year - 2000)
        vals = np.where(lost, 0, cover)[valid]
        rows.append({"site": site, "year": year,
                     "mean": float(vals.mean()) if valid.any() else np.nan,
                     "n_valid": int(valid.sum()), "n_grids": int(inside.sum())})
    return rows


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers,
                      processes=True)


def components(rows, info):
    rows = rows.drop_duplicates("year", keep="last")
    return {"tree": dict(zip(rows["year"].astype(int), rows["mean"].astype(float)))}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
