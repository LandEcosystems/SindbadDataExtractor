#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidate a site list collected from several networks into one row per site, with
FLUXNET style identifiers, coordinates, PFT and the networks that hold the site.

usage:
    python consolidate_sites.py <input csv> <output csv> [--pft-csv csv]
                                [--countries zip] [--shuttle]
                                [--shuttle-dir dir] [--shuttle-snapshot csv]
                                [--site-list NAME=csv ...] [--coord-overrides csv]
                                [--previous csv]

The input has one row per site and network, with the columns Site_ID, Lat, Lon, UTC,
PFT and one in_<network> column of 0 or 1 per network. The identifiers are of three
styles.
    fluxnet: two letter country code, a dash and three characters, as DE-Hai. Kept.
    sapfluxnet: three letter country code and parts joined by underscores, as
        ARG_MAZ. Kept, so that they still match the SAPFLUXNET data.
    name: any other name, such as Barro Colorado Island or Idamalayar, Kerala. The
        name is written in ASCII with every run of other characters replaced by an
        underscore, and gets the two letter code of the country holding the site and
        a dash in front, as PA-Barro_Colorado_Island and IN-Idamalayar_Kerala.

The country comes from the Natural Earth 10 m admin 0 map units at --countries, which
keep overseas regions such as French Guiana apart from their country. A site outside
every polygon, such as on a coast, takes the nearest polygon, and country_km gives
its distance to it.

Rows of the same site are merged in the order of SOURCE_PRIORITY. Coordinates come
from the first network in that order that has the site. When the coordinates of
the rows differ, the first row whose coordinates lie in the country of the identifier
is preferred, which drops rows with a wrong sign of latitude or longitude. PFT comes
from the first network in the order of PFT_PRIORITY that gives one, so that the
BADM PFT of the OneFlux shuttle wins. PFT_all lists every PFT of the site in that
order when the networks disagree, coord_spread_km the largest distance of any row
from the chosen coordinates, and the in_ columns are 1 when any row has the site in
that network.
source_ID lists the identifiers of the input that were merged into the site.
country is the country of the coordinates and id_country the country the identifier
names, which differ for the UK sites of FLUXNET and for RU-FLUX sites, whose prefix
is a Russian region such as KO for Komi.

UTC is the offset of local standard time from UTC in hours, which FLUXNET timestamps
use. time_zone is the IANA time zone of the coordinates from timezonefinder, and its
standard offset is the smaller of its offsets in January and July of the current
year, which leaves out summer time also where the zone defines winter time as the
exception, as in Ireland. UTC_input keeps the offset of the first network whose row
lies within SPREAD_WARN_KM of the chosen coordinates, since a row at another place
may be in another time zone. UTC takes the
standard offset of the zone where UTC_input is missing or one hour ahead of it, which
is summer time, and keeps UTC_input otherwise, so that offsets of zones that changed,
such as Santarem before 2008, still match the period of the data. UTC_source says
which of the two UTC holds, and the report lists the sites where they differ.

--pft-csv also writes the PFT of every site in the layout of PFT.csv, with the
columns Site_ID and PFT, for the pft_path of the extractors. --fallback-pft names an
older csv in that layout, whose PFT is added for the sites that have none, both those
missing from the input and those whose networks give no PFT. Values that are not a
PFT, such as TBD, are skipped. --pft-overrides names a csv with the columns Site_ID,
PFT and reason, whose PFT replaces that of the networks and of the fallback csv, for
sites that were checked by hand. PFT_source says whether the PFT of a site comes from
the networks, the fallback csv or the overrides.

--shuttle adds a row of the network in_OneFluxShuttle for every site the FLUXNET
Shuttle lists in AmeriFlux, ICOS and TERN, with its coordinates and IGBP class as
PFT. The listing is the newest snapshot in --shuttle-dir when it is less than a day
old, and is made again with the fluxnet-shuttle package otherwise, which keeps it in
that folder, the cache folder of the shuttle_obs prefetch. --shuttle-snapshot names a
listing to use instead. A listed site that is already in the input merges with its
rows as any other network does, so the shuttle only adds sites and the network flag.
The shuttle gives no UTC offset, which then comes from the time zone.

--site-list NAME=csv adds the sites of one source of eddy covariance data as rows of
the network in_NAME, such as a release of the FLUXCOM v0.07 site cube written by
fluxcom_inventory.py or LaThuile2007 from the copied legacy site files. A list has
the columns Site_ID, Lat, Lon and UTC, and may have PFT. The source lists and the
shuttle listing replace input networks that were collected by hand and are wrong or
out of date, as given in REPLACED: in_FLUXCOM by the FLUXCOM releases,
in_SINDBAD_BRK5 by in_FLUXNET2015, in_SINDBAD_LTL07 by in_LaThuile2007, and the
in_OneFluxShuttle of the input by the shuttle listing. An input network is dropped
only when all of its replacements are given, so a run without the lists, or with
only some of the FLUXCOM releases, keeps the old flags. An empty list counts as
given. The rows of a dropped network still give coordinates, UTC and PFT after all
others, and a site left without any network is removed and reported.

The FLUXNET Shuttle is the maintained source, so its coordinates come first, then
those of the FLUXCOM releases and of LaThuile2007, then those of the input networks.
coord_source names the network the coordinates come from. --coord-overrides names a
csv with the columns Site_ID, latitude, longitude and reason, whose coordinates
replace those of the networks for sites that were checked by hand. --previous names
an earlier output, by default the output file when it exists, and the report lists
the sites whose coordinates moved more than MOVE_WARN_KM from it. Their footprint
caches then need to be fetched again, which prefetch.py --check shows.
"""
import argparse
import os
import re
import unicodedata

import numpy as np
import pandas as pd

SHUTTLE_NETWORK = "in_OneFluxShuttle"
# Releases of the FLUXCOM v0.07 site cube, each with a site list.
V07_NETWORKS = ["in_FLUXNET2015", "in_AmeriFlux", "in_ICOS", "in_OzFLUX",
                "in_ICOSDrought2018", "in_WarmWinter2020", "in_CCBY"]
# Networks in the order their coordinates and UTC are preferred: the shuttle, the
# frozen fluxcom releases and LaThuile2007, then the networks of the input csv.
SOURCE_PRIORITY = ([SHUTTLE_NETWORK] + V07_NETWORKS + ["in_LaThuile2007"]
                   + ["in_FLUXCOM", "in_FLUXCOM_Z", "in_SINDBAD_BRK5",
                      "in_SINDBAD_LTL07", "in_RU_FLUX", "in_SAPFLUXNET", "in_GeoTrees",
                      "in_GeoTrees_candidates", "in_SOTR"])
# Networks in the order their PFT is preferred.
PFT_PRIORITY = [SHUTTLE_NETWORK] + [n for n in SOURCE_PRIORITY if n != SHUTTLE_NETWORK]
# Input networks that the source lists replace, with the networks that replace them.
REPLACED = {"in_FLUXCOM": V07_NETWORKS, "in_SINDBAD_BRK5": ["in_FLUXNET2015"],
            "in_SINDBAD_LTL07": ["in_LaThuile2007"], SHUTTLE_NETWORK: [SHUTTLE_NETWORK]}
COUNTRIES = os.path.expanduser("~/research/RnD/data/sindbad_ancillary_data/reference/"
                               "ne_10m_admin_0_map_units.zip")
SHUTTLE_DIR = os.path.expanduser("~/research/RnD/data/sindbad_ancillary_data/prefetch/"
                                 "fluxnet_shuttle")
FLUXNET_ID = re.compile(r"^[A-Z]{2}-[A-Za-z0-9]{3}$")
SAPFLUXNET_ID = re.compile(r"^[A-Z]{3}(_[A-Za-z0-9]+)+$")
# Codes of FLUXNET identifiers that are not the ISO code of their country.
ID_COUNTRY_ALIAS = {"UK": "GB"}
# Distances above which a site is reported.
SPREAD_WARN_KM = 1.0
COUNTRY_WARN_KM = 10.0
MOVE_WARN_KM = 0.1


def id_style(site_id):
    if FLUXNET_ID.match(site_id):
        return "fluxnet"
    if SAPFLUXNET_ID.match(site_id):
        return "sapfluxnet"
    return "name"


def clean_name(name):
    """Name in ASCII with every run of other characters replaced by an underscore."""
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^A-Za-z0-9]+", "_", ascii_name).strip("_")


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def country_codes(lats, lons, countries_path=COUNTRIES):
    """
    Two letter country code of each point and its distance in km to that country,
    0 inside it. Points outside every polygon take the nearest one.
    """
    import geopandas as gpd
    path = countries_path if countries_path.startswith("zip://") or \
        not countries_path.endswith(".zip") else "zip://" + countries_path
    units = gpd.read_file(path)[["ISO_A2_EH", "ISO_A3_EH", "geometry"]]
    pts = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lons, lats), crs="EPSG:4326")
    hit = gpd.sjoin(pts, units, how="left", predicate="within")
    hit = hit[~hit.index.duplicated()]
    codes = hit["ISO_A2_EH"].reindex(pts.index)
    dist = pd.Series(0.0, index=pts.index)
    out = codes.isna()
    if out.any():
        # Nearest polygon in an equal area projection, with its distance in km.
        eq = "EPSG:6933"
        near = gpd.sjoin_nearest(pts[out].to_crs(eq), units.to_crs(eq), how="left",
                                 distance_col="dist")
        near = near[~near.index.duplicated()]
        codes[out] = near["ISO_A2_EH"]
        dist[out] = near["dist"] / 1000
    codes = codes.where(codes != "-99")
    return codes.tolist(), dist.round(1).tolist(), units


def rank_rows(df, priority):
    """
    Rank of each row by the first network of priority that holds it. Networks of df
    missing from priority come after it, and rows without any network come last.
    """
    networks = [c for c in priority if c in df] + \
        [c for c in df if c.startswith("in_") and c not in priority]
    flags = df[networks].fillna(0).astype(int).values
    return np.where(flags.any(axis=1), flags.argmax(axis=1), len(networks)), networks


def standard_utc_offsets(lats, lons):
    """
    IANA time zone of each point and its standard offset from UTC in hours, the
    smaller of its offsets in January and July of the current year.
    """
    import datetime
    import zoneinfo
    from timezonefinder import TimezoneFinder
    finder = TimezoneFinder()
    year = datetime.date.today().year
    zones, offsets = [], []
    for lat, lon in zip(lats, lons):
        name = finder.timezone_at(lat=float(lat), lng=float(lon))
        zones.append(name)
        if name is None:
            offsets.append(np.nan)
            continue
        zone = zoneinfo.ZoneInfo(name)
        offsets.append(min(datetime.datetime(year, m, 15, 12, tzinfo=zone).utcoffset()
                           .total_seconds() / 3600 for m in (1, 7)))
    return zones, offsets


def shuttle_listing(shuttle_dir=SHUTTLE_DIR, snapshot_path=""):
    """
    Path of the FLUXNET Shuttle listing: snapshot_path when given, else the snapshot
    of shuttle_dir, made again when it is older than a day.
    """
    if snapshot_path:
        return snapshot_path
    from fluxnet_shuttle.shuttle import listall
    from extractors.fluxnet_shuttle import snapshot
    os.makedirs(shuttle_dir, exist_ok=True)
    return snapshot({"cache_path": os.path.join(shuttle_dir, "index.csv")}, listall)


def shuttle_rows(listing_path, columns):
    """
    Rows in the layout of the input for every site of a shuttle listing, flagged in
    SHUTTLE_NETWORK, with the columns of the input and the newest product of a site
    when it is listed more than once.
    """
    snap = pd.read_csv(listing_path)
    snap = snap.sort_values("last_year").drop_duplicates("site_id", keep="last")
    rows = pd.DataFrame({"Site_ID": snap["site_id"].astype(str).str.strip(),
                         "Lat": snap["location_lat"], "Lon": snap["location_long"],
                         "UTC": np.nan, "PFT": snap["igbp"]})
    for c in columns:
        if c.startswith("in_"):
            rows[c] = int(c == SHUTTLE_NETWORK)
    if SHUTTLE_NETWORK not in rows:
        rows[SHUTTLE_NETWORK] = 1
    hubs = snap["data_hub"].value_counts().to_dict()
    print(f"{len(rows)} shuttle sites from {listing_path}: {hubs}")
    return rows


def add_rows(df, rows, label):
    """The input with the rows of a source added, reporting the new sites."""
    known = set(df["Site_ID"].astype(str).str.strip().str.upper())
    new = rows[~rows["Site_ID"].str.upper().isin(known)]
    print(f"{len(new)} {label} sites are not in the input: " + " ".join(new["Site_ID"]))
    # Identifiers of a source keep the case of the input, as US-xYE.
    case = {s.upper(): s for s in df["Site_ID"].astype(str).str.strip()}
    rows = rows.copy()
    rows["Site_ID"] = rows["Site_ID"].map(lambda s: case.get(s.upper(), s))
    return pd.concat([df, rows], ignore_index=True)


def add_shuttle(df, listing_path):
    """The input with the rows of the shuttle listing, reporting the new sites."""
    return add_rows(df, shuttle_rows(listing_path, df.columns), "shuttle")


def parse_site_list(arg):
    """Network and path of a --site-list NAME=csv argument."""
    name, sep, path = arg.partition("=")
    if not sep or not name or not path:
        raise ValueError(f"--site-list takes NAME=csv, got {arg}")
    return (name if name.startswith("in_") else f"in_{name}"), path


def list_rows(path, network, columns):
    """
    Rows in the layout of the input for the sites of a source list, flagged in
    network. The list has the columns Site_ID, Lat, Lon and UTC, and may have PFT.
    """
    lst = pd.read_csv(path)
    rows = pd.DataFrame({"Site_ID": lst["Site_ID"].astype(str).str.strip(),
                         "Lat": lst["Lat"].astype(float),
                         "Lon": lst["Lon"].astype(float),
                         "UTC": lst["UTC"].astype(float) if "UTC" in lst else np.nan,
                         "PFT": lst["PFT"] if "PFT" in lst else np.nan})
    for c in columns:
        if c.startswith("in_"):
            rows[c] = 0
    rows[network] = 1
    print(f"{len(rows)} sites of {network} from {path}")
    return rows


def drop_replaced(df, given):
    """The input without the networks that the given sources replace."""
    drop = [c for c, by in REPLACED.items() if c in df and set(by) <= set(given)]
    kept = [c for c in REPLACED if c in df and c not in drop]
    if drop:
        print("input networks replaced by the source lists: " + ", ".join(drop))
    if kept:
        print("input networks kept, as no source list replaces them: " + ", ".join(kept))
    return df.drop(columns=drop)


def consolidate_sites(df, countries_path=COUNTRIES, priority=SOURCE_PRIORITY,
                      pft_priority=PFT_PRIORITY, overrides=None):
    """
    One row per site from a table with one row per site and network, as described in
    the module docstring. overrides holds coordinates checked by hand, indexed by
    site, with the columns latitude and longitude. Returns the consolidated table
    sorted by site.
    """
    df = df.copy()
    df["Site_ID"] = df["Site_ID"].astype(str).str.strip()
    df["rank"], networks = rank_rows(df, priority)
    df["pft_rank"], _ = rank_rows(df, pft_priority)
    df["id_style"] = df["Site_ID"].map(id_style)

    codes, dist, units = country_codes(df["Lat"].values, df["Lon"].values,
                                       countries_path)
    df["country_xy"], df["country_km"] = codes, dist
    a3_to_a2 = dict(zip(units["ISO_A3_EH"], units["ISO_A2_EH"]))

    # Names get their identifier from the country of the most preferred row, so
    # every row of the same name ends up with the same identifier.
    best = df.sort_values("rank", kind="stable").drop_duplicates("Site_ID")
    best_country = dict(zip(best["Site_ID"], best["country_xy"]))
    site = []
    for raw, style in zip(df["Site_ID"], df["id_style"]):
        if style == "name":
            cc = best_country.get(raw)
            site.append(f"{cc if isinstance(cc, str) else 'XX'}-{clean_name(raw)}")
        else:
            site.append(raw)
    df["site"] = site

    rows = []
    for sid, g in df.sort_values("rank", kind="stable").groupby("site", sort=True):
        style = g["id_style"].iloc[0]
        if style == "fluxnet":
            id_country = ID_COUNTRY_ALIAS.get(sid[:2], sid[:2])
        elif style == "sapfluxnet":
            id_country = a3_to_a2.get(sid[:3], np.nan)
        else:
            id_country = g["country_xy"].iloc[0]
        inside = g[g["country_xy"] == id_country]
        top = inside.iloc[0] if len(inside) else g.iloc[0]
        lat, lon = top["Lat"], top["Lon"]
        coord_source = networks[top["rank"]] if top["rank"] < len(networks) else ""
        if overrides is not None and sid in overrides.index:
            lat = float(overrides.loc[sid, "latitude"])
            lon = float(overrides.loc[sid, "longitude"])
            coord_source = "override"
        by_pft = g.sort_values("pft_rank", kind="stable")["PFT"]
        pfts = list(dict.fromkeys(by_pft.dropna().astype(str)))
        dist_km = haversine_km(lat, lon, g["Lat"].values, g["Lon"].values)
        spread = dist_km.max()
        # Only rows at the chosen coordinates give the offset, as a row of another
        # place may lie in another time zone.
        utc = g["UTC"][dist_km <= SPREAD_WARN_KM].dropna()
        row = {"site": sid, "latitude": lat, "longitude": lon,
               "UTC": np.nan, "UTC_input": utc.iloc[0] if len(utc) else np.nan,
               "UTC_source": "", "UTC_zone": np.nan, "time_zone": "",
               "PFT": pfts[0] if pfts else np.nan,
               "PFT_all": ";".join(pfts) if len(pfts) > 1 else "",
               "country": top["country_xy"], "id_country": id_country,
               "country_km": top["country_km"], "id_style": style,
               "source_ID": ";".join(dict.fromkeys(g["Site_ID"])),
               "n_rows": len(g), "coord_spread_km": round(float(spread), 3),
               "coord_source": coord_source}
        row.update({c: int(g[c].fillna(0).max()) for c in networks})
        rows.append(row)
    sites = pd.DataFrame(rows)
    fixed = sites["coord_source"] == "override"
    if fixed.any():
        codes, dist, _ = country_codes(sites.loc[fixed, "latitude"].values,
                                       sites.loc[fixed, "longitude"].values,
                                       countries_path)
        sites.loc[fixed, "country"] = codes
        sites.loc[fixed, "country_km"] = dist
    sites["time_zone"], std = standard_utc_offsets(sites["latitude"],
                                                   sites["longitude"])
    std = pd.Series(std, index=sites.index)
    summer = sites["UTC_input"] - std == 1
    use_zone = sites["UTC_input"].isna() | summer
    sites["UTC"] = sites["UTC_input"].where(~use_zone, std)
    sites["UTC_source"] = np.where(sites["UTC_input"].isna(), "time zone",
                                   np.where(summer, "time zone, input was summer time",
                                            "input"))
    sites["UTC_zone"] = std
    return sites


def report(sites):
    """Print the sites whose rows or coordinates need a look."""
    print(f"{len(sites)} sites: "
          + ", ".join(f"{n} {k}" for k, n in sites["id_style"].value_counts().items()))
    print("coordinates from: " + str(sites["coord_source"].value_counts().to_dict()))
    print(f"{sites['PFT'].isna().sum()} sites without PFT")
    print(f"{sites['UTC'].isna().sum()} sites without UTC, UTC_source: "
          + str(sites["UTC_source"].value_counts().to_dict()))
    odd = sites[(sites["UTC_source"] == "input") & (sites["UTC"] != sites["UTC_zone"])]
    print(f"{len(odd)} sites keep an input offset that differs from the standard offset "
          "of their time zone today:")
    for r in odd.itertuples():
        print(f"    {r.site}: input {r.UTC_input}, {r.time_zone} {r.UTC_zone}")
    spread = sites[sites["coord_spread_km"] > SPREAD_WARN_KM]
    print(f"{len(spread)} sites whose rows are more than {SPREAD_WARN_KM} km apart:")
    for r in spread.itertuples():
        print(f"    {r.site}: {r.coord_spread_km:.1f} km, kept {r.latitude}, "
              f"{r.longitude}")
    pft = sites[sites["PFT_all"] != ""]
    print(f"{len(pft)} sites with different PFTs, the first is kept:")
    print("    " + ", ".join(f"{r.site} {r.PFT_all}" for r in pft.itertuples()))
    other = sites[sites["country"] != sites["id_country"]]
    print(f"{len(other)} sites whose coordinates lie in another country than their "
          "identifier names:")
    for r in other.itertuples():
        print(f"    {r.site}: {r.country} at {r.latitude}, {r.longitude}")
    far = sites[sites["country_km"] > COUNTRY_WARN_KM]
    print(f"{len(far)} sites more than {COUNTRY_WARN_KM} km outside the country "
          "polygons:")
    for r in far.itertuples():
        print(f"    {r.site}: {r.country_km} km from {r.country}")


def drop_orphans(sites):
    """The sites that keep a network, printing those left without any."""
    flags = [c for c in sites if c.startswith("in_")]
    orphan = sites[flags].sum(axis=1) == 0
    if orphan.any():
        print(f"{orphan.sum()} sites have no network left after the replaced networks "
              "were dropped, and are removed: " + " ".join(sites.loc[orphan, "site"]))
    return sites[~orphan].reset_index(drop=True)


def moved_sites(sites, previous_path):
    """
    Print the sites whose coordinates moved more than MOVE_WARN_KM from those of an
    earlier output or whose UTC changed, and the sites that are new or gone.
    """
    old = pd.read_csv(previous_path)[["site", "latitude", "longitude", "UTC"]]
    both = sites.merge(old, on="site", suffixes=("", "_old"))
    utc = both[both["UTC"] != both["UTC_old"]]
    print(f"{len(utc)} sites changed their UTC from {previous_path}:")
    for r in utc.itertuples():
        print(f"    {r.site}: {r.UTC_old} to {r.UTC} ({r.UTC_source}, {r.time_zone})")
    both["moved_km"] = np.round(haversine_km(both["latitude_old"], both["longitude_old"],
                                             both["latitude"], both["longitude"]), 3)
    moved = both[both["moved_km"] > MOVE_WARN_KM].sort_values("moved_km",
                                                               ascending=False)
    new = sorted(set(sites["site"]) - set(old["site"]))
    gone = sorted(set(old["site"]) - set(sites["site"]))
    print(f"compared with {previous_path}: {len(new)} new sites, {len(gone)} sites "
          f"gone, {len(moved)} sites moved more than {MOVE_WARN_KM} km:")
    for r in moved.itertuples():
        print(f"    {r.site}: {r.moved_km} km, from {r.latitude_old}, "
              f"{r.longitude_old} to {r.latitude}, {r.longitude} ({r.coord_source})")
    if new:
        print("    new: " + " ".join(new))
    if gone:
        print("    gone: " + " ".join(gone))


def pft_table(sites, fallback_path="", overrides_path=""):
    """
    Site_ID, PFT and PFT_source of every site with a PFT, numbered from 1 like the
    index column of PFT.csv, with the PFT of fallback_path where the sites have none
    and the PFT of overrides_path replacing any other.
    """
    pft = sites[["site", "PFT"]].dropna().rename(columns={"site": "Site_ID"})
    pft["PFT_source"] = "networks"
    if fallback_path:
        old = pd.read_csv(fallback_path)[["Site_ID", "PFT"]].dropna()
        old["PFT"] = old["PFT"].astype(str).str.strip()
        old = old[~old["PFT"].isin(["", "TBD"]) & ~old["Site_ID"].isin(pft["Site_ID"])]
        old["PFT_source"] = os.path.basename(fallback_path)
        print(f"{len(old)} sites take their PFT from {fallback_path}: "
              + " ".join(old["Site_ID"]))
        pft = pd.concat([pft, old])
    if overrides_path:
        fix = pd.read_csv(overrides_path)[["Site_ID", "PFT"]].dropna()
        fix["PFT_source"] = os.path.basename(overrides_path)
        print(f"{len(fix)} sites take their PFT from {overrides_path}: "
              + " ".join(fix["Site_ID"]))
        pft = pd.concat([pft[~pft["Site_ID"].isin(fix["Site_ID"])], fix])
    pft = pft.sort_values("Site_ID")
    pft.index = range(1, len(pft) + 1)
    return pft


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--pft-csv", default="")
    p.add_argument("--fallback-pft", default="")
    p.add_argument("--pft-overrides", default="")
    p.add_argument("--countries", default=COUNTRIES)
    p.add_argument("--shuttle", action="store_true")
    p.add_argument("--shuttle-dir", default=SHUTTLE_DIR)
    p.add_argument("--shuttle-snapshot", default="")
    p.add_argument("--site-list", action="append", default=[],
                   help="NAME=csv, the sites of a source as network in_NAME")
    p.add_argument("--coord-overrides", default="")
    p.add_argument("--previous", default="",
                   help="earlier output to compare the coordinates with")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    lists = dict(parse_site_list(s) for s in args.site_list)
    shuttle = args.shuttle or args.shuttle_snapshot
    df = drop_replaced(df, list(lists) + ([SHUTTLE_NETWORK] if shuttle else []))
    if shuttle:
        df = add_shuttle(df, shuttle_listing(args.shuttle_dir, args.shuttle_snapshot))
    for network, path in lists.items():
        df = add_rows(df, list_rows(path, network, df.columns), network)
    overrides = None
    if args.coord_overrides:
        overrides = pd.read_csv(args.coord_overrides).dropna(
            subset=["latitude", "longitude"])
        overrides["Site_ID"] = overrides["Site_ID"].astype(str).str.strip()
        overrides = overrides.set_index("Site_ID")
        print(f"{len(overrides)} sites take their coordinates from "
              f"{args.coord_overrides}: " + " ".join(overrides.index))
    sites = drop_orphans(consolidate_sites(df, args.countries, overrides=overrides))
    report(sites)
    previous = args.previous or (args.output if os.path.exists(args.output) else "")
    if previous:
        moved_sites(sites, previous)
    sites.to_csv(args.output, index=False)
    print(f"written {args.output}")
    if args.pft_csv:
        pft_table(sites, args.fallback_pft, args.pft_overrides).to_csv(args.pft_csv)
        print(f"written {args.pft_csv}")


if __name__ == '__main__':
    main()
