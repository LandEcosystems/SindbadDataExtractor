#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copernicus Global Land Service LC100 collection 3 cover fractions within 1 km of the
site.

Each layer is the mean of the 100 m grid cells whose centres lie within 1 km of the
site, read remotely from the Zenodo GeoTIFFs. Zenodo allows about 133 requests a
minute, so prefetch opens each layer file once and reads the windows of all missing
sites from it, with few workers. The cache holds one row per site, year and layer,
the same format as cgls_lc100_v3_1km_sites.csv of add_vegfrac_cgls.py.

Grass is Grass + MossLichen. Vegetation is Tree + Shrub + Grass + Crops + MossLichen.
"""
import sys
import time

import numpy as np

import extractors.footprint_utils as fu

URL = ("/vsicurl/https://zenodo.org/records/{record}/files/"
       "PROBAV_LC100_global_v3.0.1_{epoch}_{layer}-CoverFraction-layer_EPSG-4326.tif")
EPOCHS = {2015: (3939038, "2015-base"), 2016: (3518026, "2016-conso"),
          2017: (3518036, "2017-conso"), 2018: (3518038, "2018-conso"),
          2019: (3939050, "2019-nrt")}
VEG_LAYERS = ["Tree", "Shrub", "Grass", "Crops", "MossLichen"]
LAYERS = VEG_LAYERS + ["Bare"]
GROUPS = {"tree": ["Tree"], "shrub": ["Shrub"], "grass": ["Grass", "MossLichen"],
          "crop": ["Crops"], "bare": ["Bare"], "veg": VEG_LAYERS}

GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
            "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
            "GDAL_HTTP_MAX_RETRY": "6", "GDAL_HTTP_RETRY_DELAY": "5"}


def years_of(info):
    return [int(y) for y in info.get("years", [2015])]


def layer_means(year, layer, sites, lats, lons):
    """Mean of one layer and year within the radius of every site."""
    import rasterio
    from rasterio.windows import Window
    record, epoch = EPOCHS[year]
    out = []
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(URL.format(record=record, epoch=epoch, layer=layer)) as src:
            for s, lat, lon in zip(sites, lats, lons):
                col0, row0, w, h, inside = fu.geo_window(src.transform, src.width,
                                                         src.height, lat, lon)
                for attempt in range(5):
                    try:
                        data = src.read(1, window=Window(col0, row0, w, h))
                        break
                    except rasterio.errors.RasterioIOError:
                        # Zenodo asks clients to wait about a minute when throttled.
                        if attempt == 4:
                            raise
                        time.sleep(60)
                data = data.astype(float)
                data[data == src.nodata] = np.nan
                v = data[inside]
                out.append({"site": s, "year": year, "layer": layer,
                            "n_grids": v.size, "n_valid": int(np.isfinite(v).sum()),
                            "mean": np.nanmean(v) if np.isfinite(v).any() else np.nan,
                            "latitude": lat, "longitude": lon})
    return out


def fetch_site(site, lat, lon, info):
    return [r for y in years_of(info) for la in LAYERS
            for r in layer_means(y, la, [site], [lat], [lon])]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    from concurrent.futures import ThreadPoolExecutor
    sites, lats, lons = fu.missing_sites(info, sites, lats, lons)
    if not sites:
        print(f"all sites are in {info['cache_path']}")
        return
    jobs = [(y, la) for y in years_of(info) for la in LAYERS]
    print(f"{len(sites)} sites, {len(jobs)} layer files", flush=True)

    def one(job):
        return layer_means(job[0], job[1], sites, lats, lons)

    # Rows are written once all layers are read, so a site is cached complete or not
    # at all.
    rows = []
    with ThreadPoolExecutor(workers or info.get("n_workers", 2)) as ex:
        for recs in ex.map(one, jobs):
            rows.extend(recs)
            print(f"done {recs[0]['year']} {recs[0]['layer']}", flush=True)
    fu.append_cache(info["cache_path"], rows)


def components(rows, info):
    rows = rows[rows["year"].isin(years_of(info)) & rows["layer"].isin(LAYERS)]
    wide = rows.drop_duplicates(["year", "layer"], keep="last").pivot(
        index="year", columns="layer", values="mean").reindex(columns=LAYERS)
    out = {}
    for g, cols in GROUPS.items():
        # Vegetation sums whatever layers are valid, the other groups need all.
        min_count = 1 if g == "veg" else len(cols)
        vals = wide[cols].sum(axis=1, min_count=min_count) / 100
        out[g] = dict(zip(wide.index.astype(int), vals.astype(float)))
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
