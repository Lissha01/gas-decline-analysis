"""
The analysis steps applied to each well:

1. full_history_fits   fit exponential, harmonic and hyperbolic to the whole first decline
2. hindcast            fit only the first 3 years, forecast the next 5, compare with what happened
3. forward_forecast    fit the most recent 5 years, forecast to the economic limit (P90/P50/P10)
"""
import numpy as np
import pandas as pd

from . import arps
from .data import anomalous_months, fit_window
from .fit import MODELS, bootstrap_fits, fit_all, fit_arps, p90_p50_p10

DMIN_ANNUAL = 0.06     # terminal decline for the modified hyperbolic, 6 %/year
Q_LIMIT = 1.0          # economic limit, Mcf/d (an assumption; see README)
END_YEAR = 2050        # forecasts stop at the end of this year at the latest


def _prepare(series, filter_anomalies=True):
    full = fit_window(series)
    flag = anomalous_months(series) if filter_anomalies else np.zeros(len(series), bool)
    start = int(np.flatnonzero(full)[0])
    return full, flag, start


def full_history_fits(series, filter_anomalies=True):
    """Fit all three models to the first decline segment (flagged anomalous months left out)."""
    full, down, start = _prepare(series, filter_anomalies)
    use = full & ~down
    t0 = series["t_days"].values[start]
    fits = fit_all(series["t_days"].values[use] - t0, series["rate"].values[use])
    for f in fits.values():
        f.t0_days = t0
    return fits, {"window": full, "anomalous": down, "start": start}


def predicted_volumes(fit, series, mask, Dmin=None):
    """Monthly volumes (Mcf) a fit predicts for the months in `mask`: rate at mid-month x days."""
    return fit.rate(series["t_days"].values[mask], Dmin=Dmin) * series["days"].values[mask]


def hindcast(series, train_months=36, test_months=60, n_boot=200, seed=0, filter_anomalies=True):
    """
    Blind test of forecasting skill.

    Cutoff = `train_months` months after the peak. Everything used to build the
    forecast (fit window, anomaly flags, fits, bootstrap) is computed from the
    history *before* the cutoff only. The test period is the fixed calendar interval
    of the next `test_months` months, every month included: months with zero or no
    reported production count as zero actual gas, and the curves are still asked to
    predict them. The number of such months is reported, because KGS lists no zero
    months, so "produced nothing" and "record missing" cannot be told apart.
    """
    first = int(np.flatnonzero(series["gas_mcf"].values > 0)[0])
    start = int(np.flatnonzero(fit_window(series.iloc[: first + 6]))[0])   # peak in first 6 months
    cut = start + train_months
    if cut + test_months > len(series):
        raise ValueError(f"series too short for a {train_months}+{test_months}-month hindcast")
    hist = series.iloc[:cut]                       # only what was known at the cutoff
    window = fit_window(hist)
    flag = anomalous_months(hist) if filter_anomalies else np.zeros(cut, bool)
    train = np.zeros(len(series), bool)
    train[:cut] = window & ~flag
    test = np.zeros(len(series), bool)
    test[cut: cut + test_months] = True
    t0 = series["t_days"].values[start]
    t_tr = series["t_days"].values[train] - t0
    q_tr = series["rate"].values[train]
    fits = fit_all(t_tr, q_tr)
    for f in fits.values():
        f.t0_days = t0
    gas = series["gas_mcf"].values[test]
    actual = float(gas.sum())
    out = {"train": train, "test": test, "fits": fits, "actual_mcf": actual, "test_months": int(test.sum()),
           "test_zero_months": int((gas == 0).sum()),
           "test_missing_records": int((~series["reported"].values[test]).sum()),
           "train_months_fitted": int(train.sum()), "train_flagged": int((window & flag).sum())}
    for m, f in fits.items():
        pred = float(predicted_volumes(f, series, test).sum())
        out[f"pred_{m}"] = pred
        out[f"err_{m}"] = (pred - actual) / actual
    boots = bootstrap_fits(t_tr, q_tr, "hyperbolic", n_boot=n_boot, seed=seed)
    preds = []
    for f in boots:
        f.t0_days = t0
        preds.append(predicted_volumes(f, series, test).sum())
    p90, p50, p10 = p90_p50_p10(preds)
    out.update({"boot_fits": boots, "hyp_p90": p90, "hyp_p50": p50, "hyp_p10": p10,
                "covered": bool(p90 <= actual <= p10)})
    return out


def _future_months(series, end_year):
    last = series.index[-1]
    idx = pd.date_range(last + pd.offsets.MonthBegin(1), f"{end_year}-12-01", freq="MS")
    days = idx.days_in_month.values.astype(float)
    mid0 = series.index[0] + pd.Timedelta(days=series["days"].values[0] / 2.0)
    mids = idx + pd.to_timedelta(days / 2.0, unit="D")
    t = (mids - mid0).total_seconds().values / 86400.0
    return idx, t, days


def forward_forecast(series, as_of, recent_months=60, q_limit=Q_LIMIT, end_year=END_YEAR,
                     Dmin_annual=DMIN_ANNUAL, n_boot=200, seed=0):
    """
    Forecast from today for a well that is still producing.

    Fits a hyperbolic (b limited to 0-1, the usual range for a well in late life) to the
    last `recent_months` months, extends it month by month as a modified hyperbolic
    (terminal decline Dmin_annual) until the rate falls below q_limit or the end of
    end_year. Uncertainty comes from a bootstrap ensemble of n_boot refits; the P50
    case is the median of that ensemble.

    The ensemble percentiles (P90 = 10th percentile of remaining gas, P10 = 90th) only
    describe spread under this model and this bootstrap; they are not real-world
    probabilities. If the P50 case is still above the limit at the end of end_year,
    `p50_reaches_limit` is False and produced + remaining is recovery to that date,
    not ultimate recovery (EUR).

    Returns None for a well with no production in the 3 months up to `as_of`.
    """
    as_of = pd.Timestamp(as_of)
    q = series["gas_mcf"].values
    if series.index[-1] < as_of - pd.DateOffset(months=2) or (q[-3:] == 0).all():
        return None
    down = anomalous_months(series)        # fine here: the forecast starts after the last month
    idx = np.arange(len(series))
    win = (idx >= len(series) - recent_months) & (q > 0) & ~down
    t_all = series["t_days"].values
    t0 = t_all[win][0]
    t_w, q_w = t_all[win] - t0, series["rate"].values[win]
    Dmin = arps.annual_effective_to_nominal(Dmin_annual)
    months, t_f, days = _future_months(series, end_year)

    boots = bootstrap_fits(t_w, q_w, "hyperbolic", n_boot=n_boot, seed=seed, b_bounds=(0.0, 1.0))
    rates = np.array([arps.rate_modified(t_f - t0, f.qi, f.Di, f.b, Dmin) for f in boots])
    alive = np.cumprod(rates >= q_limit, axis=1).astype(bool)
    vols = rates * days * alive
    remaining = vols.sum(axis=1)
    p90, p50, p10 = p90_p50_p10(remaining)
    i50 = int(np.argmin(np.abs(remaining - p50)))          # the ensemble member closest to P50
    base = fit_arps(t_w, q_w, "hyperbolic", b_bounds=(0.0, 1.0))
    base.t0_days = t0
    life = alive[i50]
    masked = np.where(alive, rates, np.nan)
    some = alive.any(axis=0)
    lo = np.full(len(months), np.nan)
    hi = np.full(len(months), np.nan)
    lo[some] = np.nanpercentile(masked[:, some], 10, axis=0)
    hi[some] = np.nanpercentile(masked[:, some], 90, axis=0)
    return {"fit": base, "window": win, "months": months,
            "rate_p50": np.where(life, rates[i50], np.nan),
            "rate_band": (lo, hi),
            "rem_p90": p90, "rem_p50": p50, "rem_p10": p10,
            "cum_to_date": float(q.sum()),
            "end_date_p50": months[life][-1] if life.any() else series.index[-1],
            "p50_reaches_limit": bool(not life[-1]),
            "b_values": np.array([f.b for f in boots])}
