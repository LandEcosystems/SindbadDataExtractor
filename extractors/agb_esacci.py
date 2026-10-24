#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Above-ground biomass within 1 km of the site from ESA CCI Biomass v7.0, merged with
the in-situ AGB of FLUXNET BADM.

The CCI value is the mean of all 100 m AGB grid cells whose centres lie within 1 km
of the site, including zero grid cells, read remotely from the CEDA GeoTIFF tiles. It
is converted from Mg ha-1 dry matter to gC m-2 with a factor of 47, assuming a carbon
fraction of 0.47. The cache holds one row per site and map year with the mean in Mg
ha-1.

Each value goes on the last time step of its year. The components are
    cci: the CCI value.
    insitu: all rows of insitu_AGB.csv for the site, site IDs matched without case.
    merged: insitu at sites with any in-situ value in the time range, else cci. Sites
        are not mixed between sources.
    merged_PFT: merged with every value set to agb_fixed at sites whose PFT is in
        pft_fixed, keeping the time stamps of merged.
"""
import sys

import numpy as np
import pandas as pd
import xarray as xr

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

URL = ("/vsicurl/https://dap.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v7.0/"
       "geotiff/{year}/{tile}_ESACCI-BIOMASS-L4-AGB-MERGED-100m-{year}-fv7.0.tif")
YEARS = [2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2015, 2016, 2017]
PX_PER_DEG = 1125
TILE_PX = 11250
# Mg ha-1 dry matter to gC m-2 with a carbon fraction of 0.47, as for globBiomass.
MGHA_TO_GCM2 = 47.0

GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
            "GDAL_HTTP_MAX_RETRY": "4", "GDAL_HTTP_RETRY_DELAY": "2"}


def years_of(info):
    return [int(y) for y in info.get("years", YEARS)]


def tile_name(top, left):
    ns = f"N{top:02d}" if top >= 0 else f"S{-top:02d}"
    ew = f"E{left:03d}" if left >= 0 else f"W{-left:03d}"
    return ns + ew


def read_window(year, row0, row1, col0, col1):
    """
    Read global grid cell rows and cols, which may span several 10 degree tiles.
    """
    import rasterio
    from rasterio.errors import RasterioIOError
    from rasterio.windows import Window
    out = np.full((row1 - row0, col1 - col0), np.nan)
    for tr in range(row0 // TILE_PX, (row1 - 1) // TILE_PX + 1):
        for tc in range(col0 // TILE_PX, (col1 - 1) // TILE_PX + 1):
            r0, r1 = max(row0, tr * TILE_PX), min(row1, (tr + 1) * TILE_PX)
            c0, c1 = max(col0, tc * TILE_PX), min(col1, (tc + 1) * TILE_PX)
            tile = tile_name(90 - 10 * tr, -180 + 10 * tc)
            try:
                with rasterio.open(URL.format(year=year, tile=tile)) as src:
                    win = Window(c0 - tc * TILE_PX, r0 - tr * TILE_PX, c1 - c0, r1 - r0)
                    data = src.read(1, window=win)
            except RasterioIOError:
                # Tiles over open ocean are not published.
                continue
            out[r0 - row0:r1 - row0, c0 - col0:c1 - col0] = data
    return out


def site_year(site, lat, lon, year):
    import rasterio
    half_lat, half_lon = fu.half_box(lat)
    row0 = int(np.floor((90 - (lat + half_lat)) * PX_PER_DEG))
    row1 = int(np.ceil((90 - (lat - half_lat)) * PX_PER_DEG))
    col0 = int(np.floor((lon - half_lon + 180) * PX_PER_DEG))
    col1 = int(np.ceil((lon + half_lon + 180) * PX_PER_DEG))
    with rasterio.Env(**GDAL_ENV):
        data = read_window(year, row0, row1, col0, col1)
    lats = 90 - (np.arange(row0, row1) + 0.5) / PX_PER_DEG
    lons = -180 + (np.arange(col0, col1) + 0.5) / PX_PER_DEG
    lon2d, lat2d = np.meshgrid(lons, lats)
    vals = data[fu.haversine(lat, lon, lat2d, lon2d) <= fu.RADIUS_M]
    n_valid = int(np.isfinite(vals).sum())
    return {"site": site, "year": year,
            "mean": np.nanmean(vals) if n_valid else np.nan,
            "n_valid": n_valid, "n_grids": len(vals)}


def fetch_site(site, lat, lon, info):
    return [site_year(site, lat, lon, y) for y in years_of(info)]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    rows = rows.drop_duplicates("year", keep="last")
    rows = rows[rows["year"].isin(years_of(info))]
    return {"cci": dict(zip(rows["year"].astype(int),
                            rows["mean"].astype(float) * MGHA_TO_GCM2))}


def insitu_series(path, site):
    """In-situ AGB per year of the site from insitu_AGB.csv, in gC m-2."""
    df = pd.read_csv(path)
    df = df[df["Site_ID"].str.lower() == site.lower()]
    df = df[df["insitu_AGB"] > 0]
    return dict(zip(df["year"].astype(int), df["insitu_AGB"].astype(float)))


def extract(dataset, site_info, config):
    from extractors.BasexTractor import BasexTractor
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)
    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None
    info = config["dataset"][dataset]
    date_ = bxtr.get_date_vec()

    rows = fu.rows_for_site(sys.modules[__name__], info, bxtr.site, bxtr.lat, bxtr.lon)
    if len(rows) == 0:
        logger.warning(f"::MISSING:: {dataset} has no cached data for {bxtr.site} in "
                       f"{info.get('cache_path', '')}. NaN values will be set. Run "
                       "prefetch.py or set fetch_missing to true.")
    cci = fu.place(date_, components(rows, info)["cci"] if len(rows) else {},
                   "year_end")
    insitu = fu.place(date_, insitu_series(info["insitu_path"], bxtr.site), "year_end")
    merged = insitu if np.isfinite(insitu).any() else cci
    pft = fu.site_pft(info.get("pft_path", ""), bxtr.site)
    merged_pft = merged.copy()
    if pft in info.get("pft_fixed", []):
        merged_pft[np.isfinite(merged_pft)] = info["agb_fixed"]
    series = {"cci": cci, "insitu": insitu, "merged": merged,
              "merged_PFT": merged_pft}

    src_data = []
    for tar_name in bxtr.vars_list:
        bxtr.log_var_start(tar_name)
        data = xr.Dataset({tar_name: xr.DataArray(
            data=series[bxtr.vars[tar_name]["component"]].copy(), dims=['time'],
            coords={'time': date_})})
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)
    src_dataset = fu.with_fetch_attrs(bxtr.merge_and_format(src_data), rows, info)
    # merged_PFT also reads pft_path, the other in-situ components only insitu_path.
    for tar_name in bxtr.vars_list:
        comp = bxtr.vars[tar_name]["component"]
        if comp in ("insitu", "merged", "merged_PFT"):
            keys = ["insitu_path"] + (["pft_path"] if comp == "merged_PFT" else [])
            src_dataset = fu.with_insitu_paths(src_dataset, [tar_name], info, keys)
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
