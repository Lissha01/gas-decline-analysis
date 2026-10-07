"""
Re-download the monthly gas production in data/kgs_gas_monthly.csv from the
Kansas Geological Survey (KGS) oil and gas production database.

    python scripts/fetch_kgs.py              # refresh all wells in data/wells.csv
    python scripts/fetch_kgs.py 1042372714   # one lease, printed to screen

Each KGS lease page has a "monthly production" link that serves a comma-separated
text file. Columns used: MONTH-YEAR (e.g. "3-2014"), PRODUCT ("G" = gas, in Mcf) and
PRODUCTION. KGS rows with month 0 (yearly totals) or -1 (starting cumulative) are skipped.

Data terms: KGS states that production data before 1987 is licensed from IHS and may
not be redistributed. All wells in this project started producing after 2000.
"""
import csv
import io
import re
import sys
import urllib.request
from pathlib import Path

BASE = "https://chasm.kgs.ku.edu/ords"
DATA = Path(__file__).resolve().parent.parent / "data"


def _get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def file_name_for_lease(lease_kid):
    """The name of the monthly-data text file that KGS generates for a lease."""
    page = _get(f"{BASE}/oil.ogl5.MonthSave?f_lc={lease_kid}")
    m = re.search(r"p_file_name=(lp\d+\.txt)", page)
    if not m:
        raise RuntimeError(f"No monthly data file found for lease {lease_kid}")
    return m.group(1)


def parse_monthly_gas(text):
    """Parse a KGS monthly text file into [(YYYY-MM, gas_mcf), ...] for product G."""
    rows = list(csv.DictReader(io.StringIO(text)))
    out = []
    for r in rows:
        if r.get("PRODUCT") != "G":
            continue
        month, year = r["MONTH-YEAR"].split("-")
        if int(month) < 1:
            continue
        out.append((f"{int(year)}-{int(month):02d}", int(float(r["PRODUCTION"]))))
    out.sort()
    return out


def fetch_lease(lease_kid):
    text = _get(f"{BASE}/qualified.anon_blobber.download?p_file_name={file_name_for_lease(lease_kid)}")
    return parse_monthly_gas(text)


def main(argv):
    if len(argv) > 1:
        for month, gas in fetch_lease(argv[1]):
            print(month, gas)
        return
    wells = list(csv.DictReader(open(DATA / "wells.csv")))
    with open(DATA / "kgs_gas_monthly.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lease_kid", "well", "month", "gas_mcf"])
        for well in wells:
            rows = fetch_lease(well["lease_kid"])
            print(f"{well['well']}: {len(rows)} months, {sum(g for _, g in rows):,} Mcf")
            for month, gas in rows:
                w.writerow([well["lease_kid"], well["well"], month, gas])


if __name__ == "__main__":
    main(sys.argv)
