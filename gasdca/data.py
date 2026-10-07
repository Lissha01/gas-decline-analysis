"""Loading and cleaning monthly gas production for decline analysis."""
import calendar
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_monthly(path=None):
    """Monthly gas volumes, one row per well-month: lease_kid, well, month, gas_mcf."""
    df = pd.read_csv(path or DATA_DIR / "kgs_gas_monthly.csv", dtype={"lease_kid": str})
    df["month"] = pd.to_datetime(df["month"] + "-01")
    return df.sort_values(["lease_kid", "month"]).reset_index(drop=True)


def load_wells(path=None):
    return pd.read_csv(path or DATA_DIR / "wells.csv", dtype={"lease_kid": str})


def well_series(df, lease_kid):
    """
    One well as a complete monthly series (missing months filled with 0), with:
        days      calendar days in the month
        rate      calendar-day average rate, Mcf/d = monthly volume / days in month
        t_days    days from the first month's midpoint
    KGS reports volumes only, not producing days, so 'rate' is a calendar-day rate:
    a month with a week of downtime shows up as a low rate.
    """
    w = df[df["lease_kid"] == lease_kid].set_index("month")["gas_mcf"]
    idx = pd.date_range(w.index.min(), w.index.max(), freq="MS")
    w = w.reindex(idx, fill_value=0)
    out = pd.DataFrame({"gas_mcf": w.values}, index=idx)
    out["days"] = [calendar.monthrange(d.year, d.month)[1] for d in idx]
    out["rate"] = out["gas_mcf"] / out["days"]
    mid = idx + pd.to_timedelta(out["days"].values / 2.0, unit="D")
    out["t_days"] = (mid - mid[0]).total_seconds().values / 86400.0
    out["cum_mcf"] = out["gas_mcf"].cumsum()
    return out


def fit_window(series, max_gap_months=3, peak_search_months=6, max_months=None):
    """
    Choose the months a decline curve is fitted to.

    1. Start at the peak rate within the first `peak_search_months` producing months
       (the clean-up and ramp-up before the peak is not decline).
    2. Stop before the first shut-in of `max_gap_months` or more consecutive zero months.
       After a long shut-in or a workover the well is on a different curve.
    3. Optionally keep only the first `max_months` (used for the hindcast test).

    Returns a boolean mask over the series' rows.
    """
    q = series["gas_mcf"].values
    n = len(q)
    producing = np.flatnonzero(q > 0)
    first = producing[0]
    search = producing[producing < first + peak_search_months]
    start = search[np.argmax(q[search])]
    end = n
    zeros = 0
    for i in range(start, n):
        zeros = zeros + 1 if q[i] == 0 else 0
        if zeros >= max_gap_months:
            end = i - zeros + 1
            break
    if max_months is not None:
        end = min(end, start + max_months)
    mask = np.zeros(n, bool)
    mask[start:end] = True
    mask &= q > 0
    return mask


def downtime_months(series, window=5, low=0.5, high=2.0, prior_months=12, prior_low=0.3):
    """
    Flag months whose rate does not reflect the reservoir's decline:

    - partial downtime: rate below `low` x the centred rolling median of the
      surrounding `window` months (a well shut in for 10 days still reports a
      full month, so its calendar-day rate dips);
    - reporting catch-up: rate above `high` x that median (volume from earlier
      months reported late);
    - end-of-life collapse: rate below `prior_low` x the median of the previous
      `prior_months` months (the centred median misses this at the very end).

    These months are left out of the fits; their gas still counts in cumulative production.
    """
    r = series["rate"].where(series["gas_mcf"] > 0)
    med = r.rolling(window, center=True, min_periods=3).median()
    prior = r.shift(1).rolling(prior_months, min_periods=6).median()
    flag = (r < low * med) | (r > high * med) | (r < prior_low * prior)
    return flag.fillna(False).values


def segment_window(series, start_date, end_date):
    """Mask for months between two dates (inclusive), producing months only."""
    m = (series.index >= pd.Timestamp(start_date)) & (series.index <= pd.Timestamp(end_date))
    return m & (series["gas_mcf"].values > 0)
