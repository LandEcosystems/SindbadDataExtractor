#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads a product of the store, <store>/<product>/<freq>/<site>.nc, written by
stage_product.py. compile_site_data uses it for every dataset of a main json with a
store block, so that a cube is only a join of staged files:
    "store": {"path": "<store>", "products": {"bgi_FLUXCOM_obs":
              "bgi_FLUXCOM_obs.{FLUXNET_version}"}}
The product of a dataset is the dataset key unless products maps it to another.

Only the variables of the expanded dataset information are returned, so sel_variables
of the dataset json still applies. The coordinates of the file must be those of the
site registry, otherwise the product was staged for coordinates that have changed
since and the site stops with an error. The global attributes of the staging, those
starting with stage_ or source_, stay in the store; the attributes of the extractor,
such as PFT, are returned with the data.
"""
import os

import xarray as xr

import logging
logger = logging.getLogger(__name__)

COORD_TOL_DEG = 1e-6


def product_path(config, dataset, site):
    store = config["store"]
    product = store.get("products", {}).get(dataset, dataset)
    return os.path.join(store["path"], product, config["temporal_resolution"],
                        f"{site}.nc")


def extract(dataset, site_info, config):
    site = site_info["site_ID"]
    path = product_path(config, dataset, site)
    if not os.path.exists(path):
        msg = f"::MISSING:: {dataset} has no store file {path} for {site}."
        logger.warning(msg)
        print(msg)
        return None
    with xr.open_dataset(path) as ds:
        ds = ds.load()
    for c in ("latitude", "longitude"):
        if c in ds.coords and abs(float(ds[c].values.ravel()[0]) - site_info[c]) > \
                COORD_TOL_DEG:
            raise ValueError(f"{dataset} at {site}: {c} {float(ds[c].values[0])} of "
                             f"{path} differs from {site_info[c]} of the registry. "
                             f"Stage the product again.")
    wanted = list(config["dataset"][dataset]["variables"])
    keep = [v for v in wanted if v in ds.data_vars]
    missing = [v for v in wanted if v not in ds.data_vars]
    if missing:
        logger.warning(f"::MISSING:: {dataset} at {site} has no {missing} in {path}.")
    out = ds[keep]
    out.attrs = {k: v for k, v in ds.attrs.items()
                 if not k.startswith(("stage_", "staged_", "source_"))}
    return out if len(keep) else None


if __name__ == '__main__':
    print(__doc__)
