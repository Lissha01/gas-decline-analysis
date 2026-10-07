# Gas-well decline analysis and production forecasting

Decline-curve analysis (DCA) of six real Kansas gas wells using public monthly production data from the Kansas Geological Survey (KGS). The project:

1. Fits three Arps decline models to each well's history.
2. Tests how well each model forecasts in a blind 5-year hindcast.
3. Forecasts the wells that are still producing, with a P90 / P50 / P10 range.

New to this? Open [`notebooks/gas_dca_walkthrough.ipynb`](notebooks/gas_dca_walkthrough.ipynb) in Google Colab. It walks through every step in plain language with worked numbers.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Lissha01/gas-decline-analysis/blob/main/notebooks/gas_dca_walkthrough.ipynb)

## Units and terms

| Term | Meaning |
|---|---|
| Mcf / MMcf | thousand / million standard cubic feet of gas |
| Mcf/d | gas rate, Mcf per day |
| qi | initial rate (Mcf/d) |
| Di | initial nominal decline rate (per year in tables, per day in code) |
| b | Arps exponent: 0 = exponential, 1 = harmonic, in between = hyperbolic |
| EUR | estimated ultimate recovery = produced to date + remaining |
| P90 / P50 / P10 | low / middle / high case: 90% / 50% / 10% chance of doing at least this well |

## Data

Six single-well leases from five Kansas gas fields, all first producing after 2000. Monthly gas volumes run to June 2026.

| Well | Field | County | First month | Months | Total produced (MMcf) |
|---|---|---|---|---|---|
| KUTTLER 'H' 1 | Bradshaw Gas Area | Greeley | 2002-08 | 173 | 220.3 |
| ROBBINS 2-27 | Glick | Kiowa | 2003-07 | 275 | 423.5 |
| BAIER #3 | Aetna Gas Area | Barber | 2004-05 | 264 | 284.5 |
| CLIFT 'A' 3-33 | Bradshaw Gas Area | Greeley | 2004-09 | 163 | 170.6 |
| Jones 24-13 | Cherry Creek Niobrara Gas Area | Cheyenne | 2009-04 | 200 | 55.6 |
| Logan 'X' 9 | Stranathan | Barber | 2010-01 | 189 | 300.3 |

- `data/wells.csv` lists each well and links to its KGS lease page.
- `data/kgs_gas_monthly.csv` holds the monthly volumes.
- `python scripts/fetch_kgs.py` downloads them again from KGS.

Source: [Kansas Geological Survey oil and gas production database](https://www.kgs.ku.edu/PRS/petroDB.html). KGS licenses production data before 1987 from IHS and does not allow it to be redistributed. All the data here is from after 2000.

## Method

**1. Rate from volume.** KGS reports gas sold per month, not producing days. So the rate is a calendar-day rate:

rate = monthly volume ÷ days in the month

For example, Logan 'X' 9 produced 3,423 Mcf in January 2010, so its rate was 3,423 ÷ 31 = 110.4 Mcf/d.

**2. Fit window.**
- The window starts at the peak month within the first 6 producing months.
- It ends at the first shut-in of 3 or more months.
- Months distorted by downtime or late reporting are left out of the fit. A month is flagged if its rate is below half, or above twice, the median of the 5 months around it, or below 30% of the median of the previous year. Their gas still counts in cumulative production.

**3. Arps models.**
- Exponential: q = qi·e^(−Di·t)
- Hyperbolic: q = qi / (1 + b·Di·t)^(1/b)
- Harmonic: the hyperbolic with b = 1

Each model is fitted by robust least squares on ln(rate). Working in logs means a 10% miss counts the same at any rate.

**4. Blind test (hindcast).** For each model:
- Fit the first 36 months after the peak.
- Predict the next 60 months.
- Compare the prediction with the gas actually produced.

**5. Forecast.**
- Applies to the three wells still producing in 2026.
- Fits a hyperbolic (b between 0 and 1) to the last 60 months.
- Switches to an exponential when the annual decline reaches 6% (the *modified hyperbolic*). This stops b from forecasting gas that never runs out.
- Each forecast ends at 1 Mcf/d or at the end of 2050, whichever comes first.

**6. Uncertainty.** A residual block bootstrap builds 200 alternative histories from the fitted trend. Each one is a different reshuffle of 6-month blocks of the fit's misfit. Every history is refitted and forecast. The 10th, 50th and 90th percentiles of remaining gas give P90, P50 and P10.

## Results

### Full-history fits

![Rate history and fits](figures/rate_fits.png)

Fitted hyperbolic b ranges from 0.52 to 1.84. Five of the six wells have b above 0.7, so their decline slows down a lot over time. That is common in tight, low-permeability gas reservoirs. All parameters are in `results/fit_parameters.csv`.

### Blind 5-year test

![Hindcast](figures/hindcast.png)
![Hindcast errors](figures/hindcast_errors.png)

| Well | Actual (MMcf) | Exponential | Harmonic | Hyperbolic |
|---|---|---|---|---|
| KUTTLER 'H' 1 | 71.0 | −7% | +17% | −7% |
| ROBBINS 2-27 | 137.2 | −13% | +6% | −13% |
| BAIER #3 | 67.0 | −50% | +10% | −26% |
| CLIFT 'A' 3-33 | 55.4 | −28% | −3% | +16% |
| Jones 24-13 | 14.5 | −46% | −2% | +23% |
| Logan 'X' 9 | 73.6 | −24% | +35% | −24% |
| **Mean absolute error** | | **28%** | **12%** | **18%** |

What the test shows:
- **The exponential model under-predicted every well.** Assuming a constant decline rate is too pessimistic for these wells.
- **Harmonic did best on average,** but over-predicted Logan by 35%. No single model wins on every well.
- **The best fit to history is not the best forecast.** Hyperbolic always fits the past best because it has an extra parameter. It still came second here.
- **The bootstrap range is too narrow.** The hyperbolic P90 to P10 range contained the actual volume for only 3 of the 6 wells. Noise-based ranges miss the bigger uncertainty, which is picking the wrong model or a change in how the well is operated.

### Forecast for producing wells

![Forecast](figures/forecast.png)

| Well | Rate now (Mcf/d) | Produced (MMcf) | Remaining P90 / P50 / P10 (MMcf) | EUR P50 (MMcf) | P50 end |
|---|---|---|---|---|---|
| ROBBINS 2-27 | 12.9 | 423.5 | 44.2 / 55.6 / 59.0 | 479.1 | still producing end-2050 |
| BAIER #3 | 6.7 | 284.5 | 26.2 / 30.5 / 31.7 | 315.0 | still producing end-2050 |
| Jones 24-13 | 1.5 | 55.6 | 0.3 / 0.5 / 1.0 | 56.2 | Oct 2027 |

Rate now is the average of the last 6 months. From the hindcast, the true range is likely wider than P90 to P10 shows.

## Limitations

- **Calendar-day rates.** Without producing days, downtime can only be filtered out approximately.
- **Sales volumes, not well tests.** Line pressure, compression and curtailment all show up in the data as if they were reservoir decline.
- **Arps is empirical.** It assumes the same operating conditions continue in the future. It does not model the reservoir, and it does not account for refracs, workovers or new wells nearby.
- **Fixed economic inputs.** The 1 Mcf/d limit and 6% terminal decline are fixed assumptions, not based on prices or costs.
- **Small sample.** Six wells is a demonstration, not a statistical study.

## How to run

```bash
pip install -r requirements.txt
python -m pytest -q        # 9 tests
python run_analysis.py     # writes results/*.csv and figures/*.png
```

GitHub Actions runs the tests and refreshes the figures and results on every push to the code or data.

## Layout

```
gasdca/arps.py        Arps rate and cumulative equations, modified hyperbolic
gasdca/data.py        loading, rate conversion, fit window, downtime filter
gasdca/fit.py         robust fitting, block bootstrap, P90/P50/P10
gasdca/workflow.py    full-history fits, hindcast, forward forecast
run_analysis.py       runs everything and writes figures and tables
scripts/fetch_kgs.py  re-downloads the data from KGS
tests/                checks on the equations, fitting, filters and data
notebooks/            plain-language walkthrough for Google Colab
```

## Reference

Arps, J.J. (1945). Analysis of decline curves. *Transactions of the AIME*, 160(1), 228–247.

## License

MIT
