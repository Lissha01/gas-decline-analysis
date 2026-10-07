"""
Arps decline-curve equations (Arps, 1945) for a gas well.

Units used everywhere in this package:
    t   time since the start of the decline, in DAYS
    q   gas rate, in Mcf/d (thousand standard cubic feet per day)
    Di  initial nominal decline rate, in 1/day
    b   Arps hyperbolic exponent (dimensionless)
    Q   cumulative production, in Mcf

The three classic cases:
    b = 0      exponential   q = qi * exp(-Di t)
    0 < b < 1  hyperbolic    q = qi / (1 + b Di t)^(1/b)
    b = 1      harmonic      q = qi / (1 + Di t)

The "modified hyperbolic" model runs the hyperbolic curve until its instantaneous
decline rate D(t) = Di / (1 + b Di t) falls to a floor Dmin, then switches to an
exponential decline at Dmin. That stops a hyperbolic curve with b close to 1 from
producing an unrealistically long tail and an inflated EUR.
"""
import numpy as np

DAYS_PER_YEAR = 365.25
B_EXP = 1e-3   # below this b, use the exponential formula (the hyperbolic one overflows as b -> 0)


def rate(t, qi, Di, b):
    """Arps rate q(t) for any b >= 0."""
    t = np.asarray(t, float)
    if b < B_EXP:
        return qi * np.exp(-Di * t)
    return qi / (1.0 + b * Di * t) ** (1.0 / b)


def log_rate(t, qi, Di, b):
    """ln q(t), computed directly so very small rates do not underflow to log(0)."""
    t = np.asarray(t, float)
    if b < B_EXP:
        return np.log(qi) - Di * t
    return np.log(qi) - np.log1p(b * Di * t) / b


def cumulative(t, qi, Di, b):
    """Arps cumulative Q(t) = integral of q from 0 to t (Mcf)."""
    t = np.asarray(t, float)
    if b < B_EXP:
        return qi / Di * (1.0 - np.exp(-Di * t))
    if abs(b - 1.0) < 1e-6:
        return qi / Di * np.log1p(Di * t)
    return qi / ((1.0 - b) * Di) * (1.0 - (1.0 + b * Di * t) ** (1.0 - 1.0 / b))


def switch_time(Di, b, Dmin):
    """
    Time (days) when the decline rate D(t) = Di / (1 + b Di t) falls to Dmin.
    0 if the decline already starts at or below Dmin; infinity for an exponential
    (b = 0) whose constant decline stays above Dmin.
    """
    if Di <= Dmin:
        return 0.0
    if b < B_EXP:
        return np.inf
    return (Di / Dmin - 1.0) / (b * Di)


def rate_modified(t, qi, Di, b, Dmin):
    """Modified hyperbolic rate: the decline rate never falls below Dmin."""
    t = np.asarray(t, float)
    ts = switch_time(Di, b, Dmin)
    if np.isinf(ts):
        return rate(t, qi, Di, b)
    qs = rate(ts, qi, Di, b)
    return np.where(t <= ts, rate(t, qi, Di, b), qs * np.exp(-Dmin * np.maximum(t - ts, 0.0)))


def cumulative_modified(t, qi, Di, b, Dmin):
    """Modified hyperbolic cumulative (Mcf)."""
    t = np.asarray(t, float)
    ts = switch_time(Di, b, Dmin)
    if np.isinf(ts):
        return cumulative(t, qi, Di, b)
    qs = rate(ts, qi, Di, b)
    Qs = cumulative(ts, qi, Di, b)
    tail = qs / Dmin * (1.0 - np.exp(-Dmin * np.maximum(t - ts, 0.0)))
    return np.where(t <= ts, cumulative(t, qi, Di, b), Qs + tail)


def nominal_to_annual_effective(Di, b):
    """First-year effective decline: the fraction of rate lost in the first year, 1 - q(1 yr)/qi."""
    return 1.0 - rate(DAYS_PER_YEAR, 1.0, Di, b)


def annual_effective_to_nominal(De):
    """Exponential nominal decline (1/day) equal to an annual effective decline De."""
    return -np.log(1.0 - De) / DAYS_PER_YEAR
