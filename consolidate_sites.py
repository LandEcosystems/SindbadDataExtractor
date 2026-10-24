#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidate a site list collected from several networks into one row per site, with
FLUXNET style identifiers, coordinates, PFT and the networks that hold the site.

usage:
    python consolidate_sites.py <input csv> <output csv> [--pft-csv csv]
                                [--countries zip]

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

Rows of the same site are merged in the order of SOURCE_PRIORITY. Coordinates and
UTC_input come from the first network in that order that has the site. When the coordinates of
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
exception, as in Ireland. UTC_input keeps the offset of the networks. UTC takes the
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
"""
import argparse
import os
import re
import unicodedata

import numpy as np
import pandas as pd

# Networks in the order their coordinates and UTC are preferred.
SOURCE_PRIORITY = ["in_FLUXCOM", "in_FLUXCOM_Z", "in_SINDBAD_BRK5", "in_SINDBAD_LTL07",
                   "in_OneFluxShuttle", "in_RU_FLUX", "in_SAPFLUXNET", "in_GeoTrees",
                   "in_GeoTrees_candidates", "in_SOTR"]
# Networks in the order their PFT is preferred.
PFT_PRIORITY = ["in_OneFluxShuttle"] + [n for n in SOURCE_PRIORITY
                                        if n != "in_OneFluxShuttle"]
COUNTRIES = os.path.expanduser("~/research/RnD/data/sindbad_ancillary_data/reference/"
                               "ne_10m_admin_0_map_units.zip")
FLUXNET_ID = re.compile(r"^[A-Z]{2}-[A-Za-z0-9]{3}$")
SAPFLUXNET_ID = re.compile(r"^[A-Z]{3}(_[A-Za-z0-9]+)+$")
# Codes of FLUXNET identifiers that are not the ISO code of their country.
ID_COUNTRY_ALIAS = {"UK": "GB"}
# Distances above which a site is reported.
SPREAD_WARN_KM = 1.0
COUNTRY_WARN_KM = 10.0


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


def consolidate_sites(df, countries_path=COUNTRIES, priority=SOURCE_PRIORITY,
                      pft_priority=PFT_PRIORITY):
    """
    One row per site from a table with one row per site and network, as described in
    the module docstring. Returns the consolidated table sorted by site.
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
        by_pft = g.sort_values("pft_rank", kind="stable")["PFT"]
        pfts = list(dict.fromkeys(by_pft.dropna().astype(str)))
        utc = g["UTC"].dropna()
        spread = haversine_km(top["Lat"], top["Lon"], g["Lat"].values,
                              g["Lon"].values).max()
        row = {"site": sid, "latitude": top["Lat"], "longitude": top["Lon"],
               "UTC": np.nan, "UTC_input": utc.iloc[0] if len(utc) else np.nan,
               "UTC_source": "", "UTC_zone": np.nan, "time_zone": "",
               "PFT": pfts[0] if pfts else np.nan,
               "PFT_all": ";".join(pfts) if len(pfts) > 1 else "",
               "country": top["country_xy"], "id_country": id_country,
               "country_km": top["country_km"], "id_style": style,
               "source_ID": ";".join(dict.fromkeys(g["Site_ID"])),
               "n_rows": len(g), "coord_spread_km": round(float(spread), 3)}
        row.update({c: int(g[c].fillna(0).max()) for c in networks})
        rows.append(row)
    sites = pd.DataFrame(rows)
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
    args = p.parse_args()

    sites = consolidate_sites(pd.read_csv(args.input), args.countries)
    report(sites)
    sites.to_csv(args.output, index=False)
    print(f"written {args.output}")
    if args.pft_csv:
        pft_table(sites, args.fallback_pft, args.pft_overrides).to_csv(args.pft_csv)
        print(f"written {args.pft_csv}")


if __name__ == '__main__':
    main()
