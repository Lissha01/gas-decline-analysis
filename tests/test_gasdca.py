import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from gasdca import arps
from gasdca.data import downtime_months, load_monthly, well_series
from gasdca.fit import fit_arps
from fetch_kgs import parse_monthly_gas


@pytest.mark.parametrize("b", [0.0, 0.5, 1.0, 1.6])
def test_cumulative_is_integral_of_rate(b):
    """The analytic Arps cumulative must equal the numerical integral of the rate."""
    qi, Di = 500.0, 0.004
    t = np.linspace(0, 3650, 200001)
    q = arps.rate(t, qi, Di, b)
    numeric = np.sum((q[1:] + q[:-1]) / 2 * np.diff(t))
    assert arps.cumulative(3650, qi, Di, b) == pytest.approx(numeric, rel=1e-6)


def test_modified_hyperbolic_switches_at_dmin():
    qi, Di, b = 500.0, 0.004, 1.2
    Dmin = arps.annual_effective_to_nominal(0.06)
    ts = arps.switch_time(Di, b, Dmin)
    assert Di / (1 + b * Di * ts) == pytest.approx(Dmin)
    # rate and cumulative are continuous at the switch
    eps = 1e-6
    assert arps.rate_modified(ts - eps, qi, Di, b, Dmin) == pytest.approx(arps.rate_modified(ts + eps, qi, Di, b, Dmin))
    assert arps.cumulative_modified(ts - eps, qi, Di, b, Dmin) == pytest.approx(
        arps.cumulative_modified(ts + eps, qi, Di, b, Dmin))
    # an exponential faster than Dmin never switches
    assert np.isinf(arps.switch_time(0.001, 0.0, Dmin))


def test_fit_recovers_known_curve():
    """Fit noisy synthetic monthly data generated from a known hyperbolic."""
    rng = np.random.default_rng(1)
    t = np.arange(120) * 30.4375
    qi, Di, b = 800.0, 0.003, 0.7
    q = arps.rate(t, qi, Di, b) * np.exp(rng.normal(0, 0.05, t.size))
    f = fit_arps(t, q, "hyperbolic")
    assert f.qi == pytest.approx(qi, rel=0.08)
    assert f.Di == pytest.approx(Di, rel=0.25)
    assert f.b == pytest.approx(b, abs=0.15)


def test_downtime_month_is_flagged():
    df = load_monthly()
    s = well_series(df, "1042372714").copy()
    s.iloc[30, s.columns.get_loc("rate")] *= 0.2
    assert downtime_months(s)[30]


def test_data_matches_kgs_totals():
    """Totals checked against the KGS files when the data were downloaded (Oct 2026)."""
    totals = load_monthly().groupby("lease_kid")["gas_mcf"].sum().to_dict()
    assert totals == {"1028540731": 220333, "1030420805": 423506, "1031510232": 284522,
                      "1031936848": 170601, "1040927054": 55633, "1042372714": 300309}


def test_parse_kgs_file():
    text = ('"LEASE_KID","LEASE","MONTH-YEAR","PRODUCT","WELLS","PRODUCTION"\n'
            '"1","X","0-2013","G","1","73000"\n'
            '"1","X","2-2013","G","1","5247"\n'
            '"1","X","1-2013","G","1","6318"\n'
            '"1","X","1-2013","O","1","12"\n')
    assert parse_monthly_gas(text) == [("2013-01", 6318), ("2013-02", 5247)]
