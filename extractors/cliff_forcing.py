#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Forcing of the CLIFF gap filling: a gridded meteorology at the nearest grid cells of
every site, filled once by prefetch.py as the pixel files that pycliff-downscale
reads, one per variable and site:
    <cache_path>/<var>/<site>/<var>_Z_BIG_<start>_<end>_pix_ALL_of_<n>.nc
with var(time, nclosest), lat(nclosest) and lon(nclosest), time in hours since
1900-01-01 UTC and the values in the units of ERA5. The variables are the ERA5 names
that PyCLIFF uses, t2m, q, tp, ssrd, ssr, strd, ws10, rH and sp, and q is taken at
the pressure level of the json. A forcing other than ERA5 is written under the same
names and units, so PyCLIFF needs no change to use it. db of the json names the
forcing in the PyCLIFF products, as ERA5.e1.

source of the json says where the values come from:
    pycliff_extract: the monthly files below pbase, ERA5 of the cluster by default.
        Another gridded product with the same file layout only needs another pbase.
        The cache is built per variable and year, see below. With pycliff_source,
        pycliff-extract reads that source of PyCLIFF instead, from the store of the
        json, such as crujra from trendy_met.zarr. Such a source keeps its own time
        step, 6 hours for CRUJRA, and marks its files with it, so PyCLIFF never runs
        it finer. It is run with the python of the PyCLIFF environment and the
        checkout pycliff_dir, one call per variable.
    <module>:<function>: any other forcing, read by a function of this repository
        called as function(var, site, lat, lon, info). It returns the UTC times as
        datetime64, the values as an array of grid cells by times, and the latitude
        and longitude of each grid cell, or None when it has no data for the site.
A file is missing when the site has none for the variable, or when its last year is
before end_year of the json, in which case the earlier file is replaced.

The monthly files are read per variable and year into year files:
    <cache_path>/<var>/_years/<var>_<year>_<tag>.nc
each with the values of one year at a set of grid cells, pixel(pixel) the flat index
lat_index * n_lon + lon_index of each cell and tag a hash of that set. A year whose
files fail to read, or whose time axis is not the hours of that year, writes nothing
and is read again by the next run, while the other years are kept. The cells are
chosen with the functions of PyCLIFF, from the checkout pycliff_dir, so a site gets
the cells pycliff-extract would give it. Only the cells of the sites with a missing
pixel file are read, and a year reads only the cells that no year file of that year
and no pixel file covering that year holds yet. New sites therefore read only their
own cells, and a later end_year only the new years.

read_strategy of the json is the way PyCLIFF reads the cells from a monthly file,
auto by default, which reads chunked files tile by tile through all hours. Files
chunked in small tiles, such as the ws10 files of 2020 and 2021 on the cluster, then
take tens of thousands of scattered reads each, so the ERA5 json sets slab, which
reads whole fields 24 hours at a time in the order they lie in the file.

A monthly file whose time does not decode to hours of its year, such as the ws10
files of 2024 and 2025 on the cluster, whose time is in seconds since 1970 but
labelled hours since 1900, is read with its raw time values and the first unit of
TIME_UNITS that puts them all on whole hours of that year. The log names the file,
its label and the unit taken. A file that no unit fits still fails the year.

Once every cell and year of a site is in year files or older pixel files, its pixel
file is written from them. Every prefetch run does this for every variable at its
end, under a lock, so the task of a cluster array that writes the last year file of
a variable also writes its pixel files. prefetch.py --part i/n splits the variable
and year pairs over the tasks, part i taking the pairs i, i + n and so on, so an
array of as many tasks as years gives each task one year of every variable.

The log of pycliff-extract is printed line by line as it runs, each line after the
variable, so it lands in the log of the prefetch job.
"""
import collections
import glob
import hashlib
import importlib
import os
import re
import subprocess
import sys
import time

import numpy as np

import extractors.footprint_utils as fu

import logging
logger = logging.getLogger(__name__)

PLEV = {"q": 1000}
VARIABLES = ["t2m", "q", "tp", "ssrd", "ssr", "strd", "ws10", "rH", "sp"]
FILE = "{var}/{site}/{var}_Z_BIG_{start}_{end}_pix_ALL_of_{n}.nc"
SPAN = re.compile(r"_Z_BIG_(\d+)_(\d+)_pix_ALL_of_")
YEARS = "_years"
# Units tried for the raw time of a monthly file whose own unit does not decode it
# to hours of its year.
TIME_UNITS = ["seconds since 1970-01-01", "hours since 1900-01-01",
              "hours since 1970-01-01", "minutes since 1970-01-01",
              "days since 1970-01-01"]
STEP_SECONDS = {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}


def variables(info):
    return info.get("variables", VARIABLES)


def years_of(info):
    return list(range(int(info["start_year"]), int(info["end_year"]) + 1))


def pixel_path(info, var, site):
    return os.path.join(info["cache_path"], FILE.format(
        var=var, site=site, start=info["start_year"], end=info["end_year"],
        n=info.get("n_pixels", 3)))


def existing(info, var, site):
    return sorted(glob.glob(os.path.join(info["cache_path"], var, site,
                                         f"{var}_Z_BIG_*_pix_ALL_of_*.nc")))


def is_missing(info, var, site):
    """No file of the variable at the site, or one that ends before end_year."""
    files = existing(info, var, site)
    if not files:
        return True
    m = SPAN.search(os.path.basename(files[-1]))
    return not m or int(m.group(2)) < int(info["end_year"])


def drop_older(info, var, site):
    """Remove the files of a site that the newest extraction replaced."""
    keep = pixel_path(info, var, site)
    if os.path.exists(keep):
        for f in existing(info, var, site):
            if f != keep:
                os.remove(f)


def write_pixels(path, var, times, values, lats, lons):
    """A pixel file in the layout of pycliff-extract."""
    import netCDF4
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("nclosest", values.shape[0])
        ds.createDimension("time", values.shape[1])
        tv = ds.createVariable("time", "f8", ("time",))
        tv.units = "hours since 1900-01-01 00:00:00"
        tv.calendar = "proleptic_gregorian"
        base = np.datetime64("1900-01-01T00:00:00", "s")
        tv[:] = (times.astype("datetime64[s]") - base) / np.timedelta64(1, "h")
        ds.createVariable("nclosest", "i8", ("nclosest",))[:] = \
            np.arange(1, values.shape[0] + 1)
        ds.createVariable("lon", "f8", ("nclosest",))[:] = lons
        ds.createVariable("lat", "f8", ("nclosest",))[:] = lats
        dv = ds.createVariable(var, "f8", ("time", "nclosest"), zlib=True, complevel=1)
        dv[:] = np.asarray(values, dtype=float).T
        ds.source = "SindbadDataExtractor cliff_forcing"
    os.replace(tmp, path)


def pycliff_env(info):
    """
    Environment of a PyCLIFF call: its output root, and the checkout pycliff_dir of
    the json first on PYTHONPATH, so that this checkout runs whatever the python of
    the PyCLIFF environment has installed.
    """
    env = dict(os.environ, PYCLIFF_OUT=info.get("pycliff_out") or info["cache_path"])
    if info.get("pycliff_dir"):
        env["PYTHONPATH"] = os.pathsep.join(
            [info["pycliff_dir"]] + [p for p in [os.environ.get("PYTHONPATH")] if p])
    return env


def run_logged(cmd, env, emit, tail=60):
    """
    Run a command and pass each line of its output, stderr included, to emit as
    it comes. Returns the exit code and the last lines, for an error message.
    """
    last = collections.deque(maxlen=tail)
    env = dict(env, PYTHONUNBUFFERED="1")
    with subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, bufsize=1) as proc:
        for line in proc.stdout:
            line = line.rstrip()
            last.append(line)
            emit(line)
    return proc.returncode, "\n".join(last)


def extract_with_pycliff(info, var, todo):
    """Run pycliff-extract for one variable of a pycliff_source and missing sites."""
    req_dir = os.path.join(info["cache_path"], "_requests")
    os.makedirs(req_dir, exist_ok=True)
    # one request file per process, so jobs of other sites never overwrite it
    batch = os.path.join(req_dir, f"{var}.{os.getpid()}.txt")
    with open(batch, "w") as f:
        for site, lat, lon in todo:
            f.write(f" -x {lon} -y {lat} -n {info.get('n_pixels', 3)} "
                    f"--year1={info['start_year']} --year2={info['end_year']} "
                    f"{var} {pixel_path(info, var, site)}\n")
    cmd = [info["pycliff_python"], "-m", "pycliff.extract.cli", "-b", batch, var,
           "--source", info["pycliff_source"], "--source-path", info["store"]]
    print("+ " + " ".join(cmd), flush=True)
    t = time.time()
    code, last = run_logged(cmd, pycliff_env(info),
                            lambda line: print(f"    {var}: {line}", flush=True))
    if code:
        raise RuntimeError(f"pycliff-extract {var} failed:\n{last}")
    for site, _, _ in todo:
        drop_older(info, var, site)
    os.remove(batch)
    print(f"{var}: {len(todo)} sites extracted in {time.time() - t:.0f} s", flush=True)


def extract_with_function(info, var, todo):
    """Write the pixel files of one variable from a function of this repository."""
    module, func = info["source"].split(":")
    fetch = getattr(importlib.import_module(module), func)
    for site, lat, lon in todo:
        try:
            got = fetch(var, site, lat, lon, info)
        except Exception as e:
            print(f"{var}: failed {site}: {e}", flush=True)
            continue
        if got is None:
            print(f"{var}: no data for {site}", flush=True)
            continue
        times, values, plats, plons = got
        write_pixels(pixel_path(info, var, site), var, np.asarray(times),
                     np.atleast_2d(values), plats, plons)
        drop_older(info, var, site)
    print(f"{var}: {len(todo)} sites written", flush=True)


# ---------------------------------------------------------------------------
# The monthly files of pbase, cached per variable and year


def pycliff_extract(info):
    """pycliff.extract.extract of the checkout pycliff_dir, imported here."""
    if info.get("pycliff_dir") and info["pycliff_dir"] not in sys.path:
        sys.path.insert(0, info["pycliff_dir"])
    from pycliff.extract import extract
    return extract


class Variable:
    """
    One variable of the monthly files: its name in the files, its pressure level,
    the grid and the cells PyCLIFF picks for a list of sites.
    """

    def __init__(self, px, info, var):
        self.px, self.info, self.var = px, info, var
        self.fvar = px.VARCOR.get(var, var)
        self.plev = str(PLEV.get(var, "."))
        first = self.files(int(info["start_year"]))[0]
        self.lons, self.lats = px.read_grid(first)

    def files(self, year):
        return self.px.file_list(self.fvar, self.plev, year, year, self.info["pbase"])

    def picks(self, sites, lats, lons):
        """Per site the (lon, lat) and flat index of its cells, nearest first."""
        n = int(self.info.get("n_pixels", 3))
        reqs = [self.px.Request(self.var, lo, la, n, 0, 0, self.plev, "")
                for la, lo in zip(lats, lons)]
        out = []
        for vals, idx in self.px.select_pixels(reqs, self.lons, self.lats):
            out.append((vals, [j * len(self.lons) + i for i, j in idx]))
        return out

    def flat_of(self, lon, lat):
        """Flat index of the cell at a grid lon and lat, or None off the grid."""
        lon = lon % 360.0
        i = int(round((lon - self.lons[0]) / (self.lons[1] - self.lons[0])))
        j = int(round((lat - self.lats[0]) / (self.lats[1] - self.lats[0])))
        i %= len(self.lons)
        if 0 <= j < len(self.lats) and abs(self.lons[i] - lon) < 1e-6 \
                and abs(self.lats[j] - lat) < 1e-6:
            return j * len(self.lons) + i
        return None

    def year_dir(self):
        return os.path.join(self.info["cache_path"], self.var, YEARS)

    def year_files(self):
        """{year: [(path, flat indices)]} of the year files."""
        import netCDF4
        out = {}
        pat = re.compile(rf"^{re.escape(self.var)}_(\d{{4}})_\w+\.nc$")
        for path in sorted(glob.glob(os.path.join(self.year_dir(), "*.nc"))):
            m = pat.match(os.path.basename(path))
            if not m:
                continue
            with netCDF4.Dataset(path) as ds:
                flat = np.asarray(ds["pixel"][:], dtype=np.int64)
            out.setdefault(int(m.group(1)), []).append((path, flat))
        return out

    def site_files(self):
        """
        Newest pixel file of each site as (path, first year, last year, flat indices
        of its cells), the files of the sites done before the year files.
        """
        import netCDF4
        out = []
        root = os.path.join(self.info["cache_path"], self.var)
        if not os.path.isdir(root):
            return out
        for site in sorted(os.listdir(root)):
            if site.startswith("_"):
                continue
            files = existing(self.info, self.var, site)
            m = SPAN.search(os.path.basename(files[-1])) if files else None
            if not m:
                continue
            try:
                with netCDF4.Dataset(files[-1]) as ds:
                    cells = [self.flat_of(lo, la) for lo, la in
                             zip(ds["lon"][:], ds["lat"][:])]
            except OSError as e:
                print(f"{self.var}: cannot read {files[-1]}: {e}", flush=True)
                continue
            if None not in cells:
                out.append((files[-1], int(m.group(1)), int(m.group(2)),
                            np.array(cells, dtype=np.int64)))
        return out

    def covered(self, year, year_files, site_files):
        """Flat indices of the cells of a year in year files or pixel files."""
        have = [f for _, f in year_files.get(year, [])]
        have += [f for _, y0, y1, f in site_files if y0 <= year <= y1]
        return np.unique(np.concatenate(have)) if have else np.zeros(0, np.int64)


def write_year(path, fvar, year, flat, lons, lats, values, pbase):
    """A year file: values(pixel, time) of one year at the cells flat."""
    import netCDF4
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    nlon = len(lons)
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
        ds.createDimension("pixel", len(flat))
        ds.createDimension("time", values.shape[1])
        ds.createVariable("pixel", "i8", ("pixel",))[:] = flat
        ds.createVariable("lon", "f8", ("pixel",))[:] = lons[flat % nlon]
        ds.createVariable("lat", "f8", ("pixel",))[:] = lats[flat // nlon]
        dv = ds.createVariable(fvar, "f8", ("pixel", "time"), zlib=True, complevel=1,
                               shuffle=True)
        dv[:] = values
        ds.year = year
        ds.pbase = pbase
        ds.source = "SindbadDataExtractor cliff_forcing"
    os.replace(tmp, path)


def decode_time(raw, unit):
    """datetime64[s] of raw time values in a unit such as hours since 1900-01-01."""
    step, _, start = unit.partition(" since ")
    epoch = np.datetime64(start.strip().replace(" ", "T"), "s")
    return epoch + np.round(raw * STEP_SECONDS[step.strip()]).astype(
        "int64").astype("timedelta64[s]")


def on_hours_of(t, year):
    """All times lie on whole hours of year."""
    t = np.asarray(t, dtype="datetime64[s]")
    return bool(len(t)) and bool(
        (t.astype("datetime64[Y]").astype(int) + 1970 == year).all()
        and (t.astype("int64") % 3600 == 0).all())


def read_month(job, year):
    """
    Values and times of the cells of a monthly file, read by PyCLIFF, whose time
    reader is wrapped here so that a file whose time does not decode to hours of
    year takes the first unit of TIME_UNITS that fits its raw values.
    """
    px = sys.modules["pycliff.extract.extract"]
    decode = getattr(px, "_decode_times", px.read_times)
    px._decode_times = decode

    def read_times(ds, time_names=None):
        tv = ds.variables[px._first_present(ds, time_names or px.TIME_NAMES)]
        try:
            t = decode(ds, time_names)
            if on_hours_of(t, year):
                return t
            why = f"decodes to {t[0]}"
        except (ValueError, OverflowError) as e:
            why = f"does not decode: {e}"
        raw = np.asarray(np.ma.filled(tv[:], np.nan), dtype=float)
        for unit in TIME_UNITS:
            t = decode_time(raw, unit)
            if on_hours_of(t, year):
                print(f"{job.path}: time labelled {getattr(tv, 'units', '')} {why}; "
                      f"read as {unit}", flush=True)
                return t
        raise ValueError(f"{job.path}: time labelled {getattr(tv, 'units', '')} {why}, "
                         f"and none of {TIME_UNITS} puts it on hours of {year}")

    px.read_times = read_times
    return px._read_pixels(job)


def extract_years(info, sites, lats, lons, part, n_workers):
    """
    Read the year files of the variable and year pairs of this part that the sites
    with a missing pixel file need. The months of all pairs are read by one pool of
    processes, and a year file is written as soon as its twelve months are in.
    """
    from concurrent.futures import ProcessPoolExecutor
    from extractors.clms_dekad import select_part
    px = pycliff_extract(info)
    pairs = select_part([(v, y) for v in variables(info) for y in years_of(info)], part)
    units = []
    for var in dict.fromkeys(v for v, _ in pairs):
        todo = [i for i, s in enumerate(sites) if is_missing(info, var, s)]
        print(f"{var}: {len(todo)} of {len(sites)} sites miss a pixel file", flush=True)
        if not todo:
            continue
        try:
            vv = Variable(px, info, var)
        except Exception as e:
            print(f"failed {var}: {type(e).__name__}: {e}", flush=True)
            continue
        need = np.unique(np.concatenate([f for _, f in vv.picks(
            [sites[i] for i in todo], [lats[i] for i in todo], [lons[i] for i in todo])]))
        year_files, site_files = vv.year_files(), vv.site_files()
        for v, year in pairs:
            if v != var:
                continue
            cells = np.setdiff1d(need, vv.covered(year, year_files, site_files))
            if len(cells):
                units.append((vv, year, cells))
    print(f"{len(units)} variable years to read in this part", flush=True)
    if not units:
        return

    def job(vv, path, cells):
        nlon = len(vv.lons)
        return px._Job(path, vv.fvar, cells // nlon, cells % nlon,
                       info.get("read_strategy", "auto"), 24, px._names())

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=max(1, n_workers)) as ex:
        futures = []
        for vv, year, cells in units:
            try:
                files = vv.files(year)
            except Exception as e:
                print(f"failed {vv.var} {year}: {e}", flush=True)
                futures.append(None)
                continue
            futures.append([ex.submit(read_month, job(vv, f, cells), year)
                            for f in files])
        for (vv, year, cells), months in zip(units, futures):
            if months is None:
                continue
            try:
                got = [m.result() for m in months]
                times = np.concatenate([g[1] for g in got])
                if not np.array_equal(times, px.expected_hours(year, year)):
                    raise ValueError("the time axis is not the hours of the year")
                values = np.concatenate([g[2] for g in got], axis=1)
            except Exception as e:
                print(f"failed {vv.var} {year}: {type(e).__name__}: {e}", flush=True)
                continue
            tag = hashlib.sha1(cells.tobytes()).hexdigest()[:10]
            write_year(os.path.join(vv.year_dir(), f"{vv.var}_{year}_{tag}.nc"),
                       vv.fvar, year, cells, vv.lons, vv.lats, values, info["pbase"])
            print(f"done {vv.var} {year}: {len(cells)} cells, "
                  f"{time.time() - t0:.0f} s", flush=True)


def combine(info, var, sites, lats, lons):
    """
    Write the pixel files of the sites of a variable whose cells are in year files or
    older pixel files for every year, under a lock so one task writes them.
    """
    import netCDF4
    px = pycliff_extract(info)
    root = os.path.join(info["cache_path"], var, YEARS)
    with fu.file_lock(os.path.join(root, "combine")):
        todo = [i for i, s in enumerate(sites) if is_missing(info, var, s)]
        if not todo:
            return
        vv = Variable(px, info, var)
        picks = vv.picks([sites[i] for i in todo], [lats[i] for i in todo],
                         [lons[i] for i in todo])
        year_files, site_files = vv.year_files(), vv.site_files()
        years = years_of(info)
        covered = {y: vv.covered(y, year_files, site_files) for y in years}
        build, lacking = [], collections.Counter()
        for k, (_, cells) in zip(todo, picks):
            gaps = [y for y in years if not np.isin(cells, covered[y]).all()]
            if gaps:
                lacking.update(gaps)
            else:
                build.append(k)
        if lacking:
            print(f"{var}: {len(todo) - len(build)} sites wait for the years "
                  f"{sorted(lacking)}", flush=True)
        if not build:
            return
        t0 = time.time()
        pick = dict(zip(todo, picks))
        cells = np.unique(np.concatenate([pick[k][1] for k in build]))
        row = {c: r for r, c in enumerate(cells)}
        hours = px.expected_hours(years[0], years[-1])
        starts = {y: int(np.searchsorted(hours, np.datetime64(f"{y}-01-01T00:00:00")))
                  for y in years + [years[-1] + 1]}
        starts[years[-1] + 1] = len(hours)
        store = np.lib.format.open_memmap(
            os.path.join(root, f"combine.{os.getpid()}.npy"), mode="w+",
            dtype=np.float64, shape=(len(cells), len(hours)))
        filled = np.zeros((len(cells), len(years)), bool)
        try:
            for iy, y in enumerate(years):
                h0, h1 = starts[y], starts[y + 1]
                for path, flat in year_files.get(y, []):
                    take = [(row[c], i) for i, c in enumerate(flat)
                            if c in row and not filled[row[c], iy]]
                    if not take:
                        continue
                    with netCDF4.Dataset(path) as ds:
                        ds.set_auto_mask(False)
                        vals = ds[vv.fvar][:]
                    r, i = map(list, zip(*take))
                    store[r, h0:h1] = vals[i]
                    filled[r, iy] = True
            for path, y0, y1, flat in site_files:
                spans = [(iy, y) for iy, y in enumerate(years) if y0 <= y <= y1]
                take = [(row[c], i) for i, c in enumerate(flat) if c in row
                        and not filled[row[c], [iy for iy, _ in spans]].all()]
                if not take:
                    continue
                src_hours = px.expected_hours(y0, y1)
                with netCDF4.Dataset(path) as ds:
                    ds.set_auto_mask(False)
                    name = next(n for n in (vv.fvar, var) if n in ds.variables)
                    vals = ds[name][:].T
                if vals.shape[1] != len(src_hours):
                    print(f"{var}: {path} has {vals.shape[1]} hours, not "
                          f"{len(src_hours)}, not used", flush=True)
                    continue
                for r, i in take:
                    for iy, y in spans:
                        if not filled[r, iy]:
                            s0 = int(np.searchsorted(
                                src_hours, np.datetime64(f"{y}-01-01T00:00:00")))
                            n = starts[y + 1] - starts[y]
                            store[r, starts[y]:starts[y + 1]] = vals[i, s0:s0 + n]
                            filled[r, iy] = True
            if not filled.all():
                raise RuntimeError(f"{var}: {int((~filled).sum())} cell years were "
                                   "not found in the files that cover them")
            for k in build:
                vals, flat = pick[k]
                write_pixels(pixel_path(info, var, sites[k]), vv.fvar, hours,
                             store[[row[c] for c in flat]],
                             [v[1] for v in vals], [v[0] for v in vals])
                drop_older(info, var, sites[k])
        finally:
            path = store.filename
            del store
            os.remove(path)
        print(f"{var}: {len(build)} pixel files written in {time.time() - t0:.0f} s",
              flush=True)


def prefetch(sites, lats, lons, info, workers=None, part="", **kwargs):
    """Fill the missing pixel files of the variables or variable years of this part."""
    from extractors.clms_dekad import select_part
    n_workers = workers or info.get("n_workers", 8)
    source = info.get("source", "pycliff_extract")
    if source == "pycliff_extract" and not info.get("pycliff_source"):
        extract_years(info, sites, lats, lons, part, n_workers)
        for var in variables(info):
            try:
                combine(info, var, sites, lats, lons)
            except Exception as e:
                print(f"failed to write the pixel files of {var}: "
                      f"{type(e).__name__}: {e}", flush=True)
        return
    for var in select_part(variables(info), part):
        todo = [s for s in zip(sites, lats, lons) if is_missing(info, var, s[0])]
        print(f"{var}: {len(todo)} of {len(sites)} sites to extract into "
              f"{info['cache_path']}", flush=True)
        if not todo:
            continue
        if source == "pycliff_extract":
            extract_with_pycliff(info, var, todo)
        else:
            extract_with_function(info, var, todo)


def report(sites, lats, lons, info):
    """Print the sites that miss a pixel file of any variable, and the year files."""
    missing = {s for s in sites for v in variables(info) if is_missing(info, v, s)}
    print(f"{info.get('db', '')}: {len(sites) - len(missing)} of {len(sites)} sites "
          f"have every variable up to {info['end_year']} in {info['cache_path']}")
    fu.print_sites(info.get("db", ""), "missing", sorted(missing), sites)
    for var in variables(info):
        files = glob.glob(os.path.join(info["cache_path"], var, YEARS, f"{var}_*.nc"))
        if files:
            years = {int(os.path.basename(f).split("_")[-2]) for f in files}
            print(f"    {var}: year files of {len(years)} years")


def complete(info, site):
    """The site has the pixel files of every variable."""
    return all(existing(info, v, site) for v in variables(info))


if __name__ == '__main__':
    print(__doc__)
