#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Vegetation fraction within 1 km of the site from Sentinel-2 L2A with the S2ToolBox
SL2P FCOVER neural network.

The network of Weiss and Baret, ATBD v1.1, takes the B03, B04, B05, B06, B07, B8A,
B11 and B12 surface reflectances and the cosines of the view zenith, sun zenith and
relative azimuth, and returns the green cover fraction. It runs on the 20 m grid, so
B03 and B04 are averaged from 10 m. The weights and normalisation bounds are those of
the Sentinel Hub custom script for FCOVER. Angles come from the 5 km angle grids of
the granule metadata at the node nearest the site.

Scenes and masking follow vegfrac_s2.py. veg is the mean over the 20 m footprint
grid cells of the per grid cell peak over the scenes. The cache holds one row per
site, the same format as s2_sl2p_fcover_2019_1km_sites.csv of add_vegfrac_s2_sl2p.py.
"""
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np

import extractors.footprint_utils as fu
from extractors.vegfrac_s2 import (GDAL_ENV, SCL_KEEP, SCL_WATER, offset, peak, sign,
                                   stack_scenes, year_of)

BANDS = ["B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12"]
# Band index of B8A in the granule metadata angle grids.
VIEW_BAND_ID = "8"

# Normalisation bounds of the inputs, in BANDS order, then cos view zenith and cos
# sun zenith. The cosine of the relative azimuth enters unscaled.
IN_MIN = np.array([0, 0, 0, 0.00663797254225, 0.0139727270189, 0.0266901380821,
                   0.0163880741923, 0, 0.918595400582, 0.342022871159])
IN_MAX = np.array([0.253061520472, 0.290393577911, 0.305398915249, 0.608900395798,
                   0.753827384323, 0.782011770669, 0.493761397883, 0.49302598446,
                   0.999999999991, 0.936206429175])
OUT_MIN, OUT_MAX = 0.000181230723879, 0.999638214715
# Hidden layer: bias, then weights for the 11 normalised inputs.
W1 = np.array([
    [-1.45261652206, -0.156854264841, 0.124234528462, 0.235625516229, -1.8323910258,
     -0.217188969888, 5.06933958064, -0.887578008155, -1.0808468167, -0.0323167041864,
     -0.224476137359, -0.195523962947],
    [-1.70417477557, -0.220824927842, 1.28595395487, 0.703139486363, -1.34481216665,
     -1.96881267559, -1.45444681639, 1.02737560043, -0.12494641532, 0.0802762437265,
     -0.198705918577, 0.108527100527],
    [1.02168965849, -0.409688743281, 1.08858884766, 0.36284522554, 0.0369390509705,
     -0.348012590003, -2.0035261881, 0.0410357601757, 1.22373853174, -0.0124082778287,
     -0.282223364524, 0.0994993117557],
    [-0.498002810205, -0.188970957866, -0.0358621840833, 0.00551248528107,
     1.35391570802, -0.739689896116, -2.21719530107, 0.313216124198, 1.5020168915,
     1.21530490195, -0.421938358618, 1.48852484547],
    [-3.88922154789, 2.49293993709, -4.40511331388, -1.91062012624, -0.703174115575,
     -0.215104721138, -0.972151494818, -0.930752241278, 1.2143441876, -0.521665460192,
     -0.445755955598, 0.344111873777]])
W2 = np.array([-0.0967998147811, 0.23080586765, -0.333655484884, -0.499418292325,
               0.0472484396749, -0.0798516540739])


def fcover(refl, view_zen, sun_zen, rel_az):
    """
    SL2P FCOVER for reflectances of shape (8, rows, cols) and scene angles in deg.
    """
    geo = np.array([np.cos(np.radians(view_zen)), np.cos(np.radians(sun_zen))])
    x = np.concatenate([refl, np.broadcast_to(geo[:, None, None],
                                              (2,) + refl.shape[1:])])
    x = 2 * (x - IN_MIN[:, None, None]) / (IN_MAX - IN_MIN)[:, None, None] - 1
    x = np.concatenate([x, np.full((1,) + refl.shape[1:], np.cos(np.radians(rel_az)))])
    hidden = np.tanh(W1[:, :1, None] + np.einsum("ni,irc->nrc", W1[:, 1:], x))
    out = W2[0] + np.einsum("n,nrc->rc", W2[1:], hidden)
    return np.clip(0.5 * (out + 1) * (OUT_MAX - OUT_MIN) + OUT_MIN, 0, 1)


def grid_value(node, path, r, c):
    rows = [list(map(float, v.text.split()))
            for v in node.find(path).iter("VALUES")]
    return np.array(rows)[r, c]


def angles(item, x, y):
    """Sun and B8A view angles at the 5 km angle grid node nearest to x, y."""
    href = item["assets"]["granule-metadata"]["href"]
    with urllib.request.urlopen(sign(href)[len("/vsicurl/"):], timeout=180) as f:
        root = ET.parse(f).getroot()
    geo = root.find(".//Tile_Geocoding")
    pos = geo.find("Geoposition[@resolution='10']")
    ulx, uly = float(pos.find("ULX").text), float(pos.find("ULY").text)
    ta = root.find(".//Tile_Angles")
    sun = ta.find("Sun_Angles_Grid")
    step = float(sun.find("Zenith/COL_STEP").text)
    r, c = int(round((uly - y) / step)), int(round((x - ulx) / step))
    sz = grid_value(sun, "Zenith/Values_List", r, c)
    sa = grid_value(sun, "Azimuth/Values_List", r, c)
    # Detectors view from alternate sides, so take the first one that covers the node
    # rather than averaging azimuths across detectors.
    for g in ta.findall(f"Viewing_Incidence_Angles_Grids[@bandId='{VIEW_BAND_ID}']"):
        vz = grid_value(g, "Zenith/Values_List", r, c)
        if np.isfinite(vz):
            return vz, grid_value(g, "Azimuth/Values_List", r, c), sz, sa
    raise ValueError("no detector covers the site")


def scene_cover(item, lat, lon):
    """SL2P FCOVER of the grid cells within the radius, on the scene 20 m grid."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.warp import transform
    from rasterio.windows import bounds as win_bounds, from_bounds
    off = offset(item)
    a = item["assets"]
    r = fu.RADIUS_M
    refl = []
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(sign(a["B05"]["href"])) as src:
            crs = str(src.crs)
            x, y = transform("EPSG:4326", src.crs, [lon], [lat])
            x, y = x[0], y[0]
            win = from_bounds(x - r - 20, y - r - 20, x + r + 20, y + r + 20,
                              src.transform)
            win = win.round_offsets().round_lengths()
            tr = src.window_transform(win)
            bounds = win_bounds(win, src.transform)
        shape = (int(win.height), int(win.width))
        for b in BANDS:
            with rasterio.open(sign(a[b]["href"])) as src:
                w = from_bounds(*bounds, src.transform).round_offsets().round_lengths()
                refl.append(src.read(1, window=w, boundless=True, fill_value=0,
                                     out_shape=shape, resampling=Resampling.average))
        with rasterio.open(sign(a["SCL"]["href"])) as src:
            w = from_bounds(*bounds, src.transform).round_offsets().round_lengths()
            scl = src.read(1, window=w, boundless=True, fill_value=0, out_shape=shape,
                           resampling=Resampling.nearest)
    vz, va, sz, sa = angles(item, x, y)
    refl = (np.stack(refl).astype(float) - off) / 10000
    fvc = fcover(refl, vz, sz, sa - va)
    fvc[scl == SCL_WATER] = 0
    fvc[~np.isin(scl, SCL_KEEP) | (refl <= 0).any(axis=0)] = np.nan
    px = tr.c + (np.arange(shape[1]) + 0.5) * tr.a
    py = tr.f + (np.arange(shape[0]) + 0.5) * tr.e
    xx, yy = np.meshgrid(px, py)
    inside = np.hypot(xx - x, yy - y) <= r
    return crs, inside, fvc


def fetch_site(site, lat, lon, info):
    rec, outs, inside = stack_scenes(site, lat, lon, info, scene_cover)
    if not outs:
        return [rec]
    arr = np.stack([o[0] for o in outs])
    ok = np.isfinite(arr).any(axis=0) & inside
    pcover = np.where(ok, peak(arr), np.nan)
    if info.get("peak_dir"):
        os.makedirs(info["peak_dir"], exist_ok=True)
        np.save(os.path.join(info["peak_dir"], f"{site}.npy"),
                pcover[ok].astype(np.float16))
    with np.errstate(all="ignore"):
        mean = np.nanmean(arr, axis=0)
    rec.update({"n_grids": int(inside.sum()), "n_valid": int(ok.sum()),
                "veg": float(np.nanmean(pcover[ok])) if ok.any() else np.nan,
                "veg_scene_mean": float(np.nanmean(mean[ok])) if ok.any()
                else np.nan})
    return [rec]


def prefetch(sites, lats, lons, info, workers=None, **kwargs):
    fu.prefetch_sites(sys.modules[__name__], info, sites, lats, lons, workers)


def components(rows, info):
    return {"veg": {year_of(info): float(rows.iloc[-1].get("veg", np.nan))}}


def extract(dataset, site_info, config):
    return fu.extract_footprint(sys.modules[__name__], dataset, site_info, config)


if __name__ == '__main__':
    import inspect
    print('---------------------------------------------------')
    print('Provider at: ', os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
