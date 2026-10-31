#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The site registry, the consolidated site list that consolidate_sites.py writes. It
gives the coordinates and UTC offset of every site and the sites of each network
flag, in place of the fluxcom site cube.

A main json names it in its registry block:
    "registry": {"path": "<sites_consolidated.csv>", "flag": "in_FLUXNET2015"}
flag selects the sites of a run. It may be a list of flags, of which any one is
enough, or "all".
"""
import os

import pandas as pd

_registries = {}


def load_registry(path):
    """The registry at path indexed by site, read once per process."""
    path = os.path.abspath(os.path.expanduser(path))
    if path not in _registries:
        reg = pd.read_csv(path)
        reg["site"] = reg["site"].astype(str).str.strip()
        _registries[path] = reg.set_index("site", drop=False)
    return _registries[path]


def site_info(site, reg):
    """
    site_ID, latitude, longitude and UTC of a site, the keys that get_site_info of
    the workflow gives from the fluxcom site cube.
    """
    if site not in reg.index:
        raise KeyError(f"{site} is not in the site registry")
    row = reg.loc[site]
    info = {"site_ID": site, "latitude": float(row["latitude"]),
            "longitude": float(row["longitude"])}
    if "UTC" in row and pd.notna(row["UTC"]):
        info["UTC"] = float(row["UTC"])
    return info


def sites_with(flag, reg):
    """Sorted sites with any of the flags set, or all sites for "all"."""
    flags = [flag] if isinstance(flag, str) else list(flag or [])
    if not flags or "all" in flags:
        return sorted(reg.index)
    missing = [f for f in flags if f not in reg.columns]
    if missing:
        raise KeyError(f"{missing} are not columns of the site registry")
    keep = (reg[flags].fillna(0).astype(int) == 1).any(axis=1)
    return sorted(reg.index[keep])
