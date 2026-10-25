"""
Paths and setup shared by the scripts that compare the FLUXNET Shuttle extraction
with the FLUXNET2015 SINDBAD site files.

DATA is the folder holding fn/ with the daily FLUXNET2015 files, hourly/ with the
hourly ones and sindbad_ancillary_data/ with the shuttle cache. Set the environment
variable SINDBAD_VALIDATION_DATA to use another folder. When the fluxcom package is
not installed, the import stub in fluxcom_stub is put on the path, which is enough
for the shuttle extractor since it never reads the FLUXCOM cube.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DATA = os.environ.get("SINDBAD_VALIDATION_DATA",
                      os.path.expanduser("~/research/RnD/data"))
CACHE = os.path.join(DATA, "sindbad_ancillary_data", "prefetch", "fluxnet_shuttle",
                     "index_fluxnet_shuttle.csv")
REF = {"daily": os.path.join(DATA, "fn", "{site}.1979.2017.daily.nc"),
       "hourly": os.path.join(DATA, "hourly", "{site}.1979.2025.hourly.nc")}

sys.path.insert(0, REPO)
try:
    import fluxcom  # noqa: F401
except ImportError:
    sys.path.insert(0, os.path.join(HERE, "fluxcom_stub"))
