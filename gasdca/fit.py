"""
Fitting Arps curves to production data, forecasting, EUR and uncertainty.

Why fit on log(rate):
    A decline spans one or two orders of magnitude (1,000 down to 50 Mcf/d). In log
    space a 10% miss counts the same early and late in life, so the early high
    rates do not dominate the fit.

Why a robust loss:
    A month with a week of downtime shows as a low calendar-day rate. With a plain
    least-squares fit those months drag the curve down. 'soft_l1' treats small
    misses like least squares but grows only linearly for big misses, so
    one-off bad months have limited pull.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from . import arps

MODELS = ("exponential", "harmonic", "hyperbolic")


@dataclass
class ArpsFit:
    model: str
    qi: float          # Mcf/d at the start of the fit window
    Di: float          # 1/day, nominal
    b: float
    rmse_log: float    # RMS misfit of ln(q) over the fit window (robust-weighted not applied)
    n: int             # number of months fitted
    aic: float
    t0_days: float = 0.0   # time origin of the fit window, in the well's own t_days

    @property
    def first_year_decline(self):
        return arps.nominal_to_annual_effective(self.Di, self.b)

    def rate(self, t_days, Dmin=None):
        t = np.asarray(t_days, float) - self.t0_days
        if Dmin is None:
            return arps.rate(t, self.qi, self.Di, self.b)
        return arps.rate_modified(t, self.qi, self.Di, self.b, Dmin)

    def cumulative(self, t_days, Dmin=None):
        t = np.asarray(t_days, float) - self.t0_days
        if Dmin is None:
            return arps.cumulative(t, self.qi, self.Di, self.b)
        return arps.cumulative_modified(t, self.qi, self.Di, self.b, Dmin)


def _initial_guess(t, q):
    slope = np.polyfit(t, np.log(q), 1)[0]
    Di0 = float(np.clip(-slope, 1e-5, 0.05))
    qi0 = float(np.exp(np.median(np.log(q[:3]))))
    return qi0, Di0


def fit_arps(t, q, model="hyperbolic", b_bounds=(0.0, 2.0), f_scale=0.15):
    """
    Fit one Arps model to rates q (Mcf/d) at times t (days, starting at 0).
    Returns an ArpsFit.
    """
    t = np.asarray(t, float)
    q = np.asarray(q, float)
    t = t - t[0]
    qi0, Di0 = _initial_guess(t, q)
    lq = np.log(q)

    if model == "exponential":
        b_fixed = 0.0
    elif model == "harmonic":
        b_fixed = 1.0
    else:
        b_fixed = None

    def unpack(p):
        qi, Di = np.exp(p[0]), np.exp(p[1])
        b = b_fixed if b_fixed is not None else p[2]
        return qi, Di, b

    def resid(p):
        qi, Di, b = unpack(p)
        return arps.log_rate(t, qi, Di, b) - lq

    p0 = [np.log(qi0), np.log(Di0)]
    lo = [np.log(qi0) - 3, np.log(1e-6)]
    hi = [np.log(qi0) + 3, np.log(0.5)]
    if b_fixed is None:
        p0.append(0.5)
        lo.append(b_bounds[0])
        hi.append(b_bounds[1])
    sol = least_squares(resid, p0, bounds=(lo, hi), loss="soft_l1", f_scale=f_scale)
    qi, Di, b = unpack(sol.x)
    r = resid(sol.x)
    n, k = len(q), len(p0)
    rmse = float(np.sqrt(np.mean(r ** 2)))
    aic = n * np.log(np.mean(r ** 2)) + 2 * k
    return ArpsFit(model, float(qi), float(Di), float(b), rmse, n, float(aic))


def fit_all(t, q, **kw):
    """Fit all three models; return {model: ArpsFit}."""
    return {m: fit_arps(t, q, m, **kw) for m in MODELS}


def bootstrap_fits(t, q, model="hyperbolic", n_boot=200, block=6, seed=0, **kw):
    """
    An ensemble of plausible fits by residual block bootstrap.

    1. Fit the model once and keep its log-residuals (the month-to-month noise).
    2. Build a new synthetic history = fitted curve x exp(residuals re-sampled in
       blocks of `block` months, so runs of high or low months stay together).
    3. Refit. Repeat n_boot times.
    Each refit is one equally plausible reading of the same noisy history; the spread
    of their forecasts is the forecast uncertainty that comes from that noise.
    Returns a list of ArpsFit (time measured from t[0]).
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(t, float)
    q = np.asarray(q, float)
    base = fit_arps(t, q, model, **kw)
    tt = t - t[0]
    fitted = arps.rate(tt, base.qi, base.Di, base.b)
    res = np.log(q) - np.log(fitted)
    n = len(q)
    starts = np.arange(max(n - block + 1, 1))
    fits = []
    for _ in range(n_boot):
        pieces, have = [], 0
        while have < n:
            s = rng.choice(starts)
            pieces.append(res[s:s + block])
            have += len(pieces[-1])
        r = np.concatenate(pieces)[:n]
        fits.append(fit_arps(t, fitted * np.exp(r), model, **kw))
    return fits


def p90_p50_p10(values):
    """Industry labels: P90 = low case (10th percentile), P50 = median, P10 = high case (90th percentile).
    Here they are percentiles of a model-and-bootstrap ensemble, not real-world probabilities."""
    v = np.asarray(values, float)
    return float(np.percentile(v, 10)), float(np.percentile(v, 50)), float(np.percentile(v, 90))
