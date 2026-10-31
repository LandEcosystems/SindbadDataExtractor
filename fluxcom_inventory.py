#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inventory of the releases of the FLUXCOM v0.07 site cube and their site lists, the
first step of the one-off freeze of the fluxcom datasets. It needs the fluxcom
package, so it runs in the main environment of the cluster, sinproc008.

usage:
    python fluxcom_inventory.py [--catalog json] [--releases R ...] [--probe n]
                                [--path-map OLD=NEW ...]

For every release of the catalog, or of --releases, the site list
<site_lists>/<release>.csv gets the columns Site_ID, Lat, Lon and UTC, from
EddyProvider(version=release).sites and get_site_lat_lon_utc, for
consolidate_sites.py --site-list. A release without sites gets a list with only the
header, so that it still counts as given.

The inventory, <store>/_inventory/v0.07_releases.csv, has one row per release with
the number of sites and the variables that EddyProvider, the surface variables of
the eddy covariance provider, JungQCProvider and TEAProvider offer, or the error of
the provider. For the first --probe sites it reads NEE and gives the years with
data. The entries of the folder of the site cube are printed as well, to show
releases that the catalog does not name. A provider without variables for a
release means that the products of that dataset are left out of its freeze.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import utils.pipeline as pl  # noqa: E402


def variable_names(make):
    """Names of the variables of the provider made by make, or the error."""
    try:
        return sorted({v.name for v in make().variables}), ""
    except Exception as e:
        return [], f"{type(e).__name__}: {e}"


def nee_years(cube, release, site):
    """First and last year with NEE at a site, or the error."""
    from fluxcom.providers import EddyProvider
    from fluxcom.core.variables import Variable
    try:
        prov = EddyProvider(cubepath=cube, version=release, site=site,
                            NEE_partitioning_method=None)
        data = prov.get_data(Variable("NEE"))
        valid = data.notnull().any(dim=[d for d in data.dims if d != "time"])
        years = pd.DatetimeIndex(data["time"].values[valid.values]).year
        return (f"{years.min()}-{years.max()}" if len(years) else "no data"), ""
    except Exception as e:
        return "", f"{type(e).__name__}: {e}"


def release_inventory(cube, release, probe, site_list_path):
    from fluxcom.providers import EddyProvider, JungQCProvider, TEAProvider
    from fluxcom.providers import eddy_covariance as ec
    from fluxcom.site_cube.site_cube_utils import get_site_lat_lon_utc
    row = {"release": release, "n_sites": 0}
    try:
        sites = sorted(EddyProvider(version=release, cubepath=cube).sites)
    except Exception as e:
        sites = []
        row["error"] = f"{type(e).__name__}: {e}"
    rows = []
    for site in sites:
        lat_lon_utc = get_site_lat_lon_utc(site=site)
        rows.append({"Site_ID": site, "Lat": lat_lon_utc[0], "Lon": lat_lon_utc[1],
                     "UTC": lat_lon_utc[2] if len(lat_lon_utc) > 2 else np.nan})
    pd.DataFrame(rows, columns=["Site_ID", "Lat", "Lon", "UTC"]).to_csv(
        site_list_path, index=False)
    print(f"{release}: {len(sites)} sites written to {site_list_path}", flush=True)
    row["n_sites"] = len(sites)
    if not sites:
        return row
    first = sites[0]
    makers = {
        "eddy": lambda: EddyProvider(cubepath=cube, version=release, site=first,
                                     NEE_partitioning_method=None),
        "surface": lambda: ec.eddy_covariance.EddyProvider(
            cubepath=cube, version=release, site=first, NEE_partitioning_method=None),
        "jung_qc": lambda: JungQCProvider(cubepath=cube, version=release, site=first),
        "tea": lambda: TEAProvider(cubepath=cube, version=release),
    }
    for name, make in makers.items():
        names, err = variable_names(make)
        if name == "surface":
            names = [n for n in names if n.startswith(("SWC", "TS"))]
        row[f"{name}_vars"] = " ".join(names)
        row[f"{name}_error"] = err
    for site in sites[:probe]:
        years, err = nee_years(cube, release, site)
        row[f"NEE_years_{site}"] = years or err
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--catalog", default=pl.CATALOG)
    p.add_argument("--releases", nargs="*", default=None)
    p.add_argument("--probe", type=int, default=3)
    p.add_argument("--path-map", action="append", default=[],
                   help="OLD=NEW, moves paths of the catalog")
    args = p.parse_args()

    cat = pl.load_catalog(args.catalog, dict(a.split("=", 1) for a in args.path_map))
    cube = cat["fluxcom_cube"]
    print(f"entries of {cube}: {sorted(os.listdir(cube))[:50]}")
    os.makedirs(cat["site_lists"], exist_ok=True)
    rows = []
    for release in args.releases or cat["releases"]:
        path = os.path.join(cat["site_lists"], f"{release}.csv")
        rows.append(release_inventory(cube, release, args.probe, path))
    out_dir = os.path.join(cat["store"], "_inventory")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "v0.07_releases.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"written {out}")


if __name__ == '__main__':
    main()
