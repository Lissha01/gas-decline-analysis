"""
Decline-curve analysis of six Kansas gas wells (public KGS data).

    python run_analysis.py

Writes results/*.csv and figures/*.png.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gasdca import arps
from gasdca.data import load_monthly, load_wells, well_series
from gasdca.fit import MODELS
from gasdca.workflow import (DMIN_ANNUAL, END_YEAR, Q_LIMIT, forward_forecast, full_history_fits, hindcast,
                             predicted_volumes)

# ---- colours (one per model, fixed order) and chart style --------------------------
COL = {"exponential": "#2a78d6", "harmonic": "#eb6834", "hyperbolic": "#1baf7a"}
INK, INK2, GRID, MUTED = "#0b0b0b", "#52514e", "#ecebe7", "#b9b8b2"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": INK2, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "figure.dpi": 120, "savefig.bbox": "tight",
    "axes.titlesize": 10, "axes.titlelocation": "left", "axes.titlecolor": INK,
})


def mmcf(x):
    return x / 1000.0


def main():
    os.makedirs("results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)
    df = load_monthly()
    wells = load_wells()
    as_of = df["month"].max()

    series, fits_full, hind, fore = {}, {}, {}, {}
    rows_fit, rows_h, rows_f = [], [], []
    for _, w in wells.iterrows():
        k, name = w["lease_kid"], w["well"]
        s = well_series(df, k)
        series[k] = s
        fits, info = full_history_fits(s)
        fits_full[k] = (fits, info)
        for m, f in fits.items():
            rows_fit.append({"well": name, "field": w["field"], "model": m, "qi_mcfd": round(f.qi, 1),
                             "Di_nominal_per_year": round(f.Di * arps.DAYS_PER_YEAR, 4), "b": round(f.b, 3),
                             "first_year_decline_pct": round(100 * f.first_year_decline, 1),
                             "rmse_ln_rate": round(f.rmse_log, 3), "aic": round(f.aic, 1), "months_fitted": f.n})
        h = hindcast(s)
        hind[k] = h
        rows_h.append({"well": name, "test_months": h["test_months"], "actual_mmcf": round(mmcf(h["actual_mcf"]), 2),
                       **{f"{m}_error_pct": round(100 * h[f"err_{m}"], 1) for m in MODELS},
                       "hyperbolic_p90_mmcf": round(mmcf(h["hyp_p90"]), 2),
                       "hyperbolic_p10_mmcf": round(mmcf(h["hyp_p10"]), 2),
                       "actual_inside_p90_p10": h["covered"]})
        f = forward_forecast(s, as_of)
        if f is not None:
            fore[k] = f
            rows_f.append({"well": name, "rate_now_mcfd": round(float(s["rate"].values[-6:].mean()), 1),
                           "produced_to_date_mmcf": round(mmcf(f["cum_to_date"]), 1),
                           "remaining_p90_mmcf": round(mmcf(f["rem_p90"]), 1),
                           "remaining_p50_mmcf": round(mmcf(f["rem_p50"]), 1),
                           "remaining_p10_mmcf": round(mmcf(f["rem_p10"]), 1),
                           "eur_p50_mmcf": round(mmcf(f["cum_to_date"] + f["rem_p50"]), 1),
                           "p50_end_date": f["end_date_p50"].strftime("%Y-%m")})

    pd.DataFrame(rows_fit).to_csv("results/fit_parameters.csv", index=False)
    hdf = pd.DataFrame(rows_h)
    hdf.to_csv("results/hindcast.csv", index=False)
    pd.DataFrame(rows_f).to_csv("results/forecast.csv", index=False)

    names = dict(zip(wells["lease_kid"], wells["well"]))
    fields = dict(zip(wells["lease_kid"], wells["field"]))
    plot_rate_fits(series, fits_full, names, fields)
    plot_hindcast(series, hind, names)
    plot_hindcast_errors(hdf)
    plot_forecast(series, fore, names)

    print(pd.DataFrame(rows_fit).query("model == 'hyperbolic'").to_string(index=False))
    print(hdf.to_string(index=False))
    print(pd.DataFrame(rows_f).to_string(index=False))
    mae = {m: float(np.mean(np.abs(hdf[f"{m}_error_pct"]))) for m in MODELS}
    print("Mean absolute 5-year hindcast error (%):", {m: round(v, 1) for m, v in mae.items()})


def _model_legend(fig, extra=()):
    handles = [plt.Line2D([], [], color=COL[m], lw=2, label=m.capitalize()) for m in MODELS]
    fig.legend(handles=list(extra) + handles, loc="upper left", ncol=len(handles) + len(extra),
               frameon=False, bbox_to_anchor=(0.0, 1.02))


def plot_rate_fits(series, fits_full, names, fields):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6), sharey=False)
    for ax, (k, s) in zip(axes.ravel(), series.items()):
        fits, info = fits_full[k]
        win, down = info["window"], info["downtime"]
        prod = s["gas_mcf"].values > 0
        ax.scatter(s.index[prod & ~win], s["rate"][prod & ~win], s=7, color=MUTED, lw=0)
        ax.scatter(s.index[win & ~down], s["rate"][win & ~down], s=8, color=INK2, lw=0)
        ax.scatter(s.index[win & down], s["rate"][win & down], s=14, facecolors="none", edgecolors=INK2, lw=0.8)
        tw = s["t_days"].values[win]
        tt = np.linspace(tw.min(), tw.max(), 300)
        dates = s.index[0] + pd.to_timedelta(tt, unit="D")
        for m in MODELS:
            ax.plot(dates, fits[m].rate(tt), color=COL[m], lw=1.8)
        f = fits["hyperbolic"]
        ax.set_yscale("log")
        ax.set_title(f"{names[k]}  ({fields[k]})")
        ax.text(0.98, 0.95, f"hyperbolic b = {f.b:.2f}\nfirst-year decline {100 * f.first_year_decline:.0f}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
        ax.set_ylabel("Gas rate (Mcf/d)")
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    extra = [plt.Line2D([], [], color=INK2, marker="o", ls="", ms=4, label="Month used in fit"),
             plt.Line2D([], [], color=INK2, marker="o", ls="", ms=5, mfc="none", label="Anomalous month (left out)"),
             plt.Line2D([], [], color=MUTED, marker="o", ls="", ms=4, label="After first long shut-in")]
    _model_legend(fig, extra)
    fig.suptitle("Arps fits to the first decline of each well (log scale)", x=0.0, ha="left", y=1.07, fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/rate_fits.png")
    plt.close(fig)


def plot_hindcast(series, hind, names):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.6))
    for ax, (k, s) in zip(axes.ravel(), series.items()):
        h = hind[k]
        tr, te = h["train"], h["test"]
        i_tr = np.flatnonzero(tr)
        i_te = np.flatnonzero(te)
        span = slice(i_tr[0], i_te[-1] + 1)
        cum = s["gas_mcf"].cumsum()
        ax.plot(s.index[span], mmcf(cum.values[span]), color=INK, lw=2, label="Actual")
        ax.axvspan(s.index[i_tr[0]], s.index[i_tr[-1]], color="#ecebe7", zorder=0)
        start_cum = cum.values[i_te[0] - 1]
        months = s.index[i_te]
        for m in MODELS:
            v = predicted_volumes(h["fits"][m], s, te)
            ax.plot(months, mmcf(start_cum + np.cumsum(v)), color=COL[m], lw=1.8)
        band = np.array([np.cumsum(predicted_volumes(f, s, te)) for f in h["boot_fits"]])
        lo, hi = np.percentile(band, 10, axis=0), np.percentile(band, 90, axis=0)
        ax.fill_between(months, mmcf(start_cum + lo), mmcf(start_cum + hi), color=COL["hyperbolic"], alpha=0.18, lw=0)
        ax.set_title(names[k])
        ax.text(0.02, 0.95, "fit on these\n3 years", transform=ax.transAxes, va="top", fontsize=8, color=INK2)
        ax.set_ylabel("Cumulative gas (MMcf)")
        ax.xaxis.set_major_locator(mdates.YearLocator(2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    extra = [plt.Line2D([], [], color=INK, lw=2, label="Actual"),
             plt.Rectangle((0, 0), 1, 1, color=COL["hyperbolic"], alpha=0.25, label="Hyperbolic P90 to P10")]
    _model_legend(fig, extra)
    fig.suptitle("Blind test: fit 3 years, forecast the next 5, compare with what happened",
                 x=0.0, ha="left", y=1.07, fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/hindcast.png")
    plt.close(fig)


def plot_hindcast_errors(hdf):
    fig, ax = plt.subplots(figsize=(8, 3.8))
    x = np.arange(len(hdf))
    wdt = 0.26
    for i, m in enumerate(MODELS):
        vals = hdf[f"{m}_error_pct"].values
        mae = np.mean(np.abs(vals))
        ax.bar(x + (i - 1) * wdt, vals, width=wdt - 0.03, color=COL[m],
               label=f"{m.capitalize()}  (mean absolute error {mae:.0f}%)")
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xticks(x, hdf["well"], fontsize=8.5)
    ax.set_ylabel("Forecast error, % of actual\n(+ = over-predicted)")
    ax.set_title("5-year forecast error after fitting 3 years: exponential under-predicts every well")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=3)
    ax.grid(axis="x", visible=False)
    fig.savefig("figures/hindcast_errors.png")
    plt.close(fig)


def plot_forecast(series, fore, names):
    ks = list(fore)
    fig, axes = plt.subplots(1, len(ks), figsize=(4.2 * len(ks), 3.8))
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, ks):
        s, f = series[k], fore[k]
        hist = s.index >= s.index[-1] - pd.DateOffset(years=10)
        prod = (s["gas_mcf"].values > 0) & hist
        ax.scatter(s.index[prod], s["rate"][prod], s=7, color=INK2, lw=0)
        ax.scatter(s.index[f["window"]], s["rate"][f["window"]], s=7, color=INK, lw=0)
        lo, hi = f["rate_band"]
        ax.fill_between(f["months"], lo, hi, color=COL["hyperbolic"], alpha=0.18, lw=0)
        ax.plot(f["months"], f["rate_p50"], color=COL["hyperbolic"], lw=2)
        ax.axhline(Q_LIMIT, color=INK2, lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set_title(names[k])
        ax.text(0.98, 0.95, f"Remaining to {END_YEAR} (MMcf)\nP90 {mmcf(f['rem_p90']):.1f} / P50 {mmcf(f['rem_p50']):.1f} / "
                f"P10 {mmcf(f['rem_p10']):.1f}", transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
        ax.set_ylabel("Gas rate (Mcf/d)")
        ax.xaxis.set_major_locator(mdates.YearLocator(10))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    handles = [plt.Line2D([], [], color=INK, marker="o", ls="", ms=4, label="Last 5 years (fitted)"),
               plt.Line2D([], [], color=COL["hyperbolic"], lw=2, label="P50 forecast (modified hyperbolic)"),
               plt.Rectangle((0, 0), 1, 1, color=COL["hyperbolic"], alpha=0.25, label="P90 to P10"),
               plt.Line2D([], [], color=INK2, lw=0.8, ls="--", label=f"Economic limit {Q_LIMIT:g} Mcf/d")]
    fig.legend(handles=handles, loc="upper left", ncol=4, frameon=False, bbox_to_anchor=(0.0, 1.06))
    fig.suptitle(f"Forecast for the wells still producing (terminal decline {100 * DMIN_ANNUAL:.0f}%/yr)",
                 x=0.0, ha="left", y=1.12, fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/forecast.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
