# Data (not included: to be filled)

This repository does not ship any data. Place the files described below in this folder, using exactly these names. The scripts read them from these locations.

```
data/
  panels/
    ringpair_yellow.parquet
    ringpair_uber.parquet
    ringpair_lyft.parquet
    ringpair_citibike.parquet
    ringpair_subway.parquet
  covariates/
    rainsnow_cov.csv
    gusts_cov.csv
```

The panels were derived from public records:

- NYC TLC trip records (Yellow Taxi, Uber, Lyft)
- MTA Subway origin-destination ridership estimates
- Citi Bike trip data
- Open-Meteo historical hourly weather
- NYC DCP 2020 Neighborhood Tabulation Areas (NTAs)

The ring-construction and panel-building code is not part of this repository. The definitions below are the ones used in the paper.

## Rings and keys

- **Rings, per mode.**
  - Ring 0 is the Congestion Relief Zone (CRZ): the 17 Manhattan NTAs at or below 60th Street.
  - All non-CRZ NTAs served by the mode are clustered (Ward, two clusters) on z-scored travel-linkage, cost, accessibility and land-use features (paper, Section 3.2).
  - The cluster closer to the CRZ is Ring 1 (near). The other cluster is Ring 2 (far).
- **Workdays.** Monday-Friday, excluding U.S. federal and NY State holidays (see the Subway note under `avg_daily_rides`).
- **Seasons.** These are calendar-year meteorological seasons:
  - winter = Jan, Feb, Dec of the same calendar year
  - spring = Mar-May
  - summer = Jun-Aug
  - fall = Sep-Nov
- **Peak windows.**
  - `am_peak` = hourly buckets 7, 8, 9.
  - `pm_peak` = hourly buckets 17, 18, 19.
  - `off_peak` = the other 18 hours.

## Panels: `panels/ringpair_<mode>.parquet`

One file per mode: `yellow`, `uber`, `lyft`, `citibike`, `subway`.

- Each file has **1,080 rows** = 9 ring pairs x 2 years x 4 seasons x 5 weekdays x 3 peak windows.
- Rows are unique on `(ring_flow, year, season, dow, peak)`.
- Row order does not matter.
- Category labels must match exactly, because the scripts match on these strings.

"Req." marks columns that at least one script reads.

| column | type | req. | definition |
|---|---|---|---|
| `mode` | string | | Mode name, constant within a file. |
| `year` | int64 | yes | `2024` or `2025`. |
| `season` | string | yes | `winter`, `spring`, `summer`, `fall`. |
| `dow` | string | yes | `Mon`, `Tue`, `Wed`, `Thu`, `Fri`. |
| `peak` | string | yes | `am_peak`, `off_peak`, `pm_peak`. |
| `origin_ring` | int32 | | Origin ring: 0 = CRZ, 1 = near, 2 = far. |
| `dest_ring` | int32 | | Destination ring, coded the same way. |
| `ring_flow` | string | yes | `"<origin_ring>-><dest_ring>"`. 9 levels: `0->0`, `0->1`, `0->2`, `1->0`, `1->1`, `1->2`, `2->0`, `2->1`, `2->2`. |
| `avg_daily_rides` | float64 | yes | Average daily rides on the ring pair in the peak window, over the workdays of `(year, season, dow)`. Yellow Taxi, Uber, Lyft and Citi Bike: total trips divided by `n_days`. Subway: the source is published as monthly averages by weekday and hour, so the monthly per-day totals are averaged over the months of the season, weighted by each month's count of that weekday. Holidays cannot be separated in that source. |
| `log1p_avg_daily_rides` | float64 | yes | Response: `log1p(avg_daily_rides)`. |
| `n_days` | int64 (float64 in the subway file) | | Number of workdays behind the average. |
| `n_od_cells` | int64 (float64 in the subway file) | | Number of NTA-level OD x hour source records aggregated into the row. |
| `post` | bool | yes | `year >= 2025`. Must be boolean. |
| `treated` | bool | yes | `origin_ring == 0 or dest_ring == 0`, i.e. 5 treated and 4 control ring pairs. Must be boolean. |
| `temperature_2m` | float64 | yes | Air temperature at 2 m, °C. See *Weather*. |
| `precipitation` | float64 | | Precipitation, mm. Not used; replaced by `rain` and `snowfall` in the covariate files. |
| `wind_speed_10m` | float64 | | Wind speed at 10 m, km/h. Not used. |
| `relative_humidity_2m` | float64 | yes | Relative humidity at 2 m, %. See *Weather*. |
| `population` | float64 | | Mean of the origin-ring and destination-ring total population, from NTA population. Not used. |
| `log_population` | float64 | | `log(population)`. Not used. |
| `income` | float64 | | Mean of the ride-weighted median worker earnings of the origin and destination NTAs. Not used. |
| `log_income` | float64 | | `log(income)`. Not used. |
| `ammi` | float64 | yes | AMMI covariate, constant within `(ring_flow, peak)`. See *AMMI*. |

## Covariates: `covariates/*.csv`

Each file covers all five modes in one CSV:

- 1,080 rows = 5 modes x 9 ring pairs x 2 years x 4 seasons x 3 peak windows.
- Rows are unique on `(mode, ring_flow, year, season, peak)`.
- There is no `dow` column, because these covariates are constant across weekdays.
- The scripts merge them onto each panel by `(ring_flow, year, season, peak)` after filtering on `mode`.
- Every panel row must find a match. The scripts stop on missing values.

`rainsnow_cov.csv` (all scripts merge it; only the weather-selection sweep uses its values):

| column | type | definition |
|---|---|---|
| `mode` | string | `yellow`, `uber`, `lyft`, `citibike`, `subway`. |
| `ring_flow` | string | As in the panels. |
| `year` | int64 | `2024`, `2025`. |
| `season` | string | As in the panels. |
| `peak` | string | As in the panels. |
| `rain` | float64 | Rain (past hour), mm. |
| `snowfall` | float64 | Snowfall (past hour), cm. |

`gusts_cov.csv`:

| column | type | definition |
|---|---|---|
| `mode`, `ring_flow`, `year`, `season`, `peak` | | Same keys as above. |
| `wind_gusts_10m` | float64 | Wind gusts at 10 m, km/h. |

## Weather

All weather variables, in the panels and in the covariate files, are built the same way:

1. Start from Open-Meteo historical hourly data on a grid of 17 cells covering NYC.
2. For each ring of a mode, average over a regular grid of points inside the ring's NTA polygons. Each point takes its nearest weather cell.
3. Keep Monday-Friday hours only.
4. Average over the hours of the peak window within `(year, season)`.
5. A ring pair's value is the mean of its origin-ring and destination-ring values.

As a result, weather varies by `(mode, ring_flow, year, season, peak)` and is constant across `dow`. The scripts z-score each variable within mode (population SD) before fitting.

## AMMI

The AMMI covariate is computed per mode from the pre-period (2024) panel rows:

1. Form the 9 x 3 matrix of means of `log1p_avg_daily_rides` over `ring_flow` x `peak`, averaging over seasons and weekdays. Use rows in `ring_flow` order and columns `am_peak`, `off_peak`, `pm_peak`.
2. Double-centre the matrix.
3. Take the leading singular vectors `g`, `h`, with unit norm and no singular-value scaling.
4. Flip both signs if needed, so that the largest-magnitude entry of `h` is positive.
5. Set `ammi = g[ring_flow] * h[peak]`.

The following reproduces the stored column to about 1e-15:

```python
import numpy as np
pre = df[df.year == 2024]
M  = pre.groupby(["ring_flow", "peak"])["log1p_avg_daily_rides"].mean().unstack("peak")
Mc = M.values - M.values.mean(1, keepdims=True) - M.values.mean(0, keepdims=True) + M.values.mean()
U, d, Vt = np.linalg.svd(Mc, full_matrices=False)
g, h = U[:, 0], Vt[0]
if h[np.argmax(np.abs(h))] < 0:
    g, h = -g, -h
df["ammi"] = df.ring_flow.map(dict(zip(M.index, g))) * df.peak.map(dict(zip(M.columns, h)))
```
