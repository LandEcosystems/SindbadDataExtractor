#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combine the per-site NetCDF files of run_extraction.py into one zarr store with a
site dimension. Only the dimensions and coordinates follow the layout of
FLUXNET_v2023_12_1D_addVar.zarr. Variable names, values and attributes are those of
the NetCDF files.

usage:
    python combine_sites_zarr.py <output zarr> <nc dir> [nc dir ...] [--sites s ...]
                                 [--overwrite]

Each nc dir is a data/<resolution> folder of a run, holding <site>.<start>.<end>.
<resolution>.nc files. When a site is in several folders, the file of the first
folder is used, so the folders are given in order of priority, for example one per
FLUXNET version.

In the store, time is renamed Ti and stored as days since the first day, the size one
latitude and longitude dimensions are dropped, and every variable gets the leading
site dimension. Other dimensions, such as depth_soilGrids, are kept. A site without
a variable gets NaN for it, and time steps outside the record of a site are NaN.
The date variables year, month, day, hour and julian_day are the same for all sites,
so they keep only the Ti dimension.

Attributes that are equal for all sites stay single values. Attributes that differ,
such as SITE_ID, PFT and latitude in the global attributes or bounds of a variable,
become lists in the order of the site dimension, with an empty string for sites
that lack them. Each variable is written as one chunk, compressed with blosc lz4,
and the metadata is consolidated at the end.
"""
import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

# Date variables of the site files, equal for all sites on the same time steps.
DATE_VARS = {"year", "month", "day", "hour", "julian_day"}
SITE_DIMS = ("latitude", "longitude")


def site_files(dirs, sites=None):
    """Site to file, taking the first folder that has the site."""
    out = {}
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "*.nc"))):
            site = os.path.basename(f).split(".")[0]
            if sites and site not in sites:
                continue
            if site in out:
                print(f"{site}: using {out[site]}, skipping {f}")
                continue
            out[site] = f
    return dict(sorted(out.items()))


def plain(v):
    """Attribute value as a json friendly python value."""
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    return v


def merge_attrs(per_site):
    """
    One attribute dict from a list of per-site dicts. Values equal for all sites stay
    single, the others become lists over sites with an empty string where missing.
    """
    keys = []
    for a in per_site:
        keys += [k for k in a if k not in keys and k != "_FillValue"]
    out = {}
    for k in keys:
        vals = [plain(a.get(k, "")) for a in per_site]
        same = all(k in a for a in per_site) and all(v == vals[0] for v in vals)
        out[k] = vals[0] if same else vals
    return out


def date_variable(per_site, template):
    """
    A date variable over Ti from the sites that have it, each time step taken from
    the first site covering it, in the type and attributes of the site files.
    """
    vals = np.full(template.shape, np.nan)
    for da in per_site:
        if da is not None:
            v = da.values.astype(float)
            vals = np.where(np.isnan(vals), v, vals)
    if not np.isnan(vals).any():
        vals = vals.astype(template.dtype)
    return xr.DataArray(vals, dims=template.dims, coords=template.coords,
                        attrs=dict(template.attrs))


def squeeze_site(da):
    """Drop the size one latitude and longitude dimensions and rename time to Ti."""
    for dim in SITE_DIMS:
        if dim in da.dims:
            da = da.isel({dim: 0}, drop=True)
    da = da.drop_vars([c for c in SITE_DIMS if c in da.coords])
    return da.rename({"time": "Ti"}) if "time" in da.dims else da


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("output")
    p.add_argument("nc_dirs", nargs="+")
    p.add_argument("--sites", nargs="*", default=None)
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()

    import numcodecs
    import zarr
    if os.path.exists(args.output) and not args.overwrite:
        sys.exit(f"{args.output} exists, use --overwrite to replace it")

    files = site_files(args.nc_dirs, args.sites)
    if not files:
        sys.exit(f"no site files in {args.nc_dirs}")
    sites = list(files)
    dss = [xr.open_dataset(f) for f in files.values()]
    print(f"{len(sites)} sites")

    # Union of the time steps and of the other coordinates over all sites.
    times = np.unique(np.concatenate([ds["time"].values for ds in dss]))
    coords = {"site": sites, "Ti": times}
    names = []
    for ds in dss:
        for c in ds.coords:
            if c not in ("time",) + SITE_DIMS and c not in coords:
                coords[c] = ds[c].values
        names += [v for v in ds.data_vars if v not in names]

    compressor = numcodecs.Blosc(cname="lz4", clevel=5,
                                 shuffle=numcodecs.Blosc.SHUFFLE)
    base = xr.Dataset(coords=coords)
    base["Ti"].encoding = {"units": f"days since {str(times[0])[:10]}",
                           "calendar": "standard", "dtype": "float64"}
    attrs = merge_attrs([ds.attrs for ds in dss])
    base.to_zarr(args.output, mode="w", consolidated=False)

    for n, name in enumerate(names):
        per_site = []
        template = None
        for ds in dss:
            if name in ds:
                da = squeeze_site(ds[name])
                da = da.reindex({d: coords[d] for d in da.dims})
                per_site.append(da)
                template = template if template is not None else da
            else:
                per_site.append(None)
        if name in DATE_VARS:
            out = date_variable(per_site, template)
        else:
            shape = (len(sites),) + template.shape
            data = np.full(shape, np.nan,
                           dtype=np.result_type(template.dtype, np.float32))
            for i, da in enumerate(per_site):
                if da is not None:
                    data[i] = da.values
            dims = ("site",) + template.dims
            out = xr.DataArray(data, dims=dims, coords={d: coords[d] for d in dims})
            out.attrs = merge_attrs([da.attrs if da is not None else {}
                                     for da in per_site])
        shape = out.shape
        out.encoding = {"compressor": compressor, "chunks": shape}
        # Keep the fill value of the site files for the integer date variables.
        if name in DATE_VARS and "_FillValue" in template.encoding:
            out.encoding["_FillValue"] = template.encoding["_FillValue"]
        # The coordinates are in the store already, so only the variable is added.
        out = out.to_dataset(name=name)
        out.drop_vars(list(out.coords)).to_zarr(args.output, mode="a",
                                                consolidated=False)
        print(f"{n + 1}/{len(names)} {name} {dict(zip(out[name].dims, shape))}",
              flush=True)

    # Appending variables replaces the group attributes, so they are set last.
    zarr.open_group(args.output, mode="r+").attrs.update(attrs)
    zarr.consolidate_metadata(args.output)
    for ds in dss:
        ds.close()
    print(f"written {args.output}")


if __name__ == '__main__':
    main()
