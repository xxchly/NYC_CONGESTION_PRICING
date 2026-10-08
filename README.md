# Cross-Classified Random-Effects Difference-in-Differences: NYC Congestion Pricing

Model-fitting code for

> Han, L. & Di, X. **Cross-Classified Random-Effects Difference-in-Differences Model: A Case Study for New York City's Congestion Pricing.** Transportation Research Board (TRB) Annual Meeting 2027 (accepted for presentation); under review at *Transportation Research Record*.

```bibtex
@misc{han_di_ccre_did_2027,
  author = {Han, Luyang and Di, Xuan},
  title  = {Cross-Classified Random-Effects Difference-in-Differences Model:
            A Case Study for New York City's Congestion Pricing},
  year   = {2027},
  note   = {Accepted for presentation at the TRB Annual Meeting 2027;
            under review at Transportation Research Record}
}
```

Authors: Luyang Han (Department of Statistics, Columbia University) and Xuan (Sharon) Di (Department of Civil Engineering and Engineering Mechanics and Data Science Institute, Columbia University).

## Summary

On January 5, 2025 New York City activated the Central Business District Tolling Program, a cordon charge over the Congestion Relief Zone (CRZ) in Manhattan at or below 60th Street. Yellow Taxi, Uber and Lyft trips are tolled per trip; Subway and Citi Bike are not. The paper estimates the first-year effect on all five modes. Each mode's neighborhoods outside the CRZ are clustered into a near ring (Ring 1) and a far ring (Ring 2), with Ring 0 being the CRZ. Trips are aggregated to directed origin-destination (OD) ring pairs crossed with morning-peak, evening-peak and off-peak windows. A **cross-classified random-effects DiD (M2)** is fit once per mode. One fit returns the pooled treatment effect and 15 flow-level effects (5 treated ring pairs x 3 peak windows). Parallel trends is tested at the pooled level and flow by flow with a residual-difference slope test (RDST).

Headline results: pooled ridership fell significantly for Yellow Taxi (-20.3%) and Uber (-6.2%) but not for Lyft. Subway and Citi Bike show no pooled effect. At the flow level, substitution toward the untolled modes is confined to the commute direction within Ring 1: Subway carries the morning inbound flow and Citi Bike the evening outbound flow. Suppression of the tolled modes appears across all rings and windows.

## Model

Notation (per mode $m$, superscript dropped):

- $r$: directed OD ring pair $i\to j$, $i,j\in\{0,1,2\}$ (9 levels). The pair is treated iff it touches Ring 0, so 5 pairs are treated and 4 are controls.
- $t=(y,s,d,p)$: year x season x weekday x peak window (2 x 4 x 5 x 3 = 120 levels).
- $c=(r,p)$: a *flow*, i.e. one ring pair in one peak window. There are 15 treated flows.
- Response: $Y_{r,t}=\log\bigl(1+\overline{\text{rides}}_{r,t}\bigr)$, where $\overline{\text{rides}}$ is the average daily ridership over workdays.
- $\mathrm{Post}_y=\mathbf 1(y\ge 2025)$. The DiD regressor is $\mathrm{Post}_y\,\mathrm{Treated}_r$.
- Covariates: z-scored temperature, wind gusts and relative humidity (the adopted weather set), plus the AMMI covariate. AMMI is the rank-1 interaction score $g_r h_p$ from an SVD of the double-centred 2024 (pre-period) mean-response matrix over ring pair x peak window (Gauch 1988).

**M0, two-way fixed-effects DiD** (`lm`):
$Y_{r,t}=\lambda_r+\lambda_t+\beta_1\,\mathrm{Post}_y\mathrm{Treated}_r+\textstyle\sum_k\beta_{k+1}W^{(k)}_{r,t}+\varepsilon_{r,t}$

**M1, interacted TWFE DiD** (`lm`): the same model, but with one fixed DiD coefficient $\beta_{1,c}$ per treated flow.

**M2, cross-classified random-effects DiD** (adopted; `lme4::lmer`, REML):

$$
Y_{r,t}\sim\mathcal N(\mu_{r,t},\sigma^2),\qquad
\mu_{r,t}=\beta_0+u_r+u_t+\bigl(\beta_1+u_{r,p}\bigr)\mathrm{Post}_y\mathrm{Treated}_r
+\beta_2 z(\text{temp})+\beta_3 z(\text{gusts})+\beta_4 z(\text{RH})+\beta_5\,\mathrm{AMMI}_{r,p}
$$

$$
u_r\sim\mathcal N(0,\sigma_r^2),\quad u_t\sim\mathcal N(0,\sigma_t^2),\quad u_{r,p}\sim\mathcal N(0,\sigma_{rp}^2)\ \text{(independent)}.
$$

The ring-pair and time intercepts are crossed, not nested. The DiD slope has a random effect on the flow. The pooled effect is $\beta_1$. The flow-level effects are the BLUPs $\beta_1+u_{r,p}$. Effects are reported as $\exp(\beta)-1$. In lme4 syntax (`code/fit_m2_full.R`):

```r
Y ~ PostTreated + (0 + PostTreated | ring_flow:peak)
  + z_temperature_2m + z_wind_gusts_10m + z_relative_humidity_2m
  + (1 | ring_flow) + (1 | year:season:dow:peak) + ammi
# lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)); Satterthwaite df via lmerTest
```

**RDST, the parallel-trends test.** For each treated flow, take the M2 residuals. Subtract the mean residual of the 4 control ring pairs in the same $(y,s,d,p)$ to get $dR$. Then fit `dR ~ tau + Post + (1 | year:season)`, where `tau` is the calendar-season index 1..8. The flow passes if the Satterthwaite p-value of the `tau` slope is at least 0.05. The pooled test applies the same regression to the treated flows' $dR$ averaged with weights proportional to 2024 rides per window hour.

## Repository layout

```
code/                         analysis scripts (run from the repo root, see below)
  fit_m2_full.R               M2 per mode: fixed effects, flow BLUPs, variance components, flow RDST
  gate_table3__tw79.R         parallel-trends table: RDST + pre-period mean test, flow and pooled
  integration_pvals_tw79.R    model selection: M2 vs M0 (pooled) and M1 (15 flows), Welch tests
  pt_by_model__tw79.R         RDST run on M0 / M1 / M2 residuals
  weather_sweep__tw79.R       weather selection, all 32 subsets on M2 (REML est., ML AIC)
  weather_sweep_m0m1__tw79.R  the same sweep on M0 and M1
  ammi_dow_m2__tw79.R         AMMI selection on M2, plus 225-cell DiD grid on M2
  ammi_necessity_tw79.py      AMMI selection on M1 (OLS), plus pooled M0
  decomp_dow_m1__tw79.py      225-cell DiD grid on M1 (OLS)
  fig_style.py                shared figure style and paths
  fig_cordon.py / fig_ringtriple.py / fig_covariates.py / fig_ammi_necessity.py   figures
  gen_tables.py / gen_selection_tables.py                                         LaTeX tables
data/                         NOT INCLUDED, see data/README.md
  panels/                     ringpair_<mode>.parquet  (5 files)
  covariates/                 rainsnow_cov.csv, gusts_cov.csv
results/  figures/  tables/   created and filled by the scripts
```

The suffix `tw79` in file names refers to the adopted peak windows: morning 07-09 and evening 17-19.

## Data

**The data are not included in this repository.** The panels were built from public records:

- NYC Taxi & Limousine Commission trip records (Yellow Taxi; high-volume for-hire vehicles: Uber, Lyft)
- MTA Subway origin-destination ridership estimates
- Citi Bike trip data
- Open-Meteo historical hourly weather
- NYC Department of City Planning 2020 NTA boundaries

The ring construction and panel-building code is also not included. `data/README.md` documents the exact file names and the full column schema the scripts expect.

## Requirements

- **R** with `arrow`, `data.table`, `lme4`, `lmerTest`.
- **Python 3** with `numpy`, `pandas`, `pyarrow`, `scipy`, `matplotlib`.

Tested environment: R 4.5.1 (arrow 23.0.1.2, data.table 1.17.6, lme4 2.0.1, lmerTest 3.2.1) and Python 3.11 (numpy 2.4.4, pandas 3.0.2, pyarrow 23.0.1, scipy 1.17.1, matplotlib 3.10.9). With the original panels, this environment reproduces the manuscript's result tables exactly.

```r
install.packages(c("arrow", "data.table", "lme4", "lmerTest"))
```
```bash
pip install numpy pandas pyarrow scipy matplotlib
```

## How to run

Put the data files in `data/` (see `data/README.md`), then run from the repository root. R scripts take the repo root as an optional first argument or from the environment variable `CCRE_DID_ROOT`; the default is the working directory. Python scripts locate the repo from their own path.

```bash
# 1. Main model, parallel trends, model selection
Rscript code/fit_m2_full.R               # must run before gate_table3__tw79.R
Rscript code/gate_table3__tw79.R
Rscript code/integration_pvals_tw79.R
Rscript code/pt_by_model__tw79.R

# 2. Variable-selection experiments
Rscript code/weather_sweep__tw79.R
Rscript code/weather_sweep_m0m1__tw79.R
Rscript code/ammi_dow_m2__tw79.R
python  code/ammi_necessity_tw79.py
python  code/decomp_dow_m1__tw79.py

# 3. Figures and LaTeX tables
python code/fig_cordon.py
python code/fig_ringtriple.py
python code/fig_covariates.py
python code/fig_ammi_necessity.py
python code/gen_tables.py
python code/gen_selection_tables.py
```

The whole pipeline takes a few minutes on a laptop; the 160-fit weather sweep is the slowest step. lme4 prints `boundary (singular) fit` messages for some of the small RDST regressions, where the year-season random-intercept variance is estimated at zero. These messages are expected.

## Outputs

| Paper exhibit (TRB 2027 manuscript) | Script(s) | Output |
|---|---|---|
| Table 5, weather covariate selection | `weather_sweep__tw79.R`, `weather_sweep_m0m1__tw79.R`, then `gen_selection_tables.py` | `tables/tab_weather_sel.tex` |
| Table 6, DiD interaction-dimension selection | `decomp_dow_m1__tw79.py`, `ammi_dow_m2__tw79.R`, then `gen_selection_tables.py` | `tables/tab_grid_sel.tex` |
| Figure 3, AMMI covariate selection | `ammi_necessity_tw79.py`, `ammi_dow_m2__tw79.R`, then `fig_ammi_necessity.py` | `figures/fig_ammi_necessity.{pdf,png}` |
| Table 7, M2 reproduces M0 and M1 | `integration_pvals_tw79.R`, then `gen_tables.py` | `tables/tab_m2_reproduces.tex` |
| Table 8, full M2 results | `fit_m2_full.R`, then `gen_tables.py` | `tables/tab_m2_full.tex` |
| Table 9, parallel-trends tests | `fit_m2_full.R`, `gate_table3__tw79.R`, then `gen_tables.py` | `tables/tab_parallel_trends.tex` |
| Figure 4, flow-level effects on the ring cordon | `fit_m2_full.R`, then `fig_cordon.py` | `figures/fig_cordon_peakdirection.{pdf,png}` |
| Figure 5, effects by ring distance | `fit_m2_full.R`, then `fig_ringtriple.py` | `figures/fig_ringtriple.{pdf,png}` |
| Figure 6, weather and AMMI coefficients | `fit_m2_full.R`, then `fig_covariates.py` | `figures/fig_covariates.{pdf,png}` |
| Text: RDST under M0 / M1 / M2 (63 / 65 / 65 of 75) | `pt_by_model__tw79.R` | `results/pt_by_model__tw79.csv` |

Intermediate CSVs in `results/`:

- `fit_m2_full.R` writes the following, all keyed by `mode`:
  - `fixed__tw79.csv`: fixed effects with Satterthwaite df and p.
  - `did_blup__tw79.csv`: the 15 flow effects $\beta_1+u_{r,p}$ with conditional SD and 95% interval.
  - `fitstats__tw79.csv`: n, logLik, AIC, BIC, sigma, pooled $\beta_1$ with CI, and convergence/singularity flags.
  - `varcomp__tw79.csv`: $\sigma^2_{rp},\sigma^2_r,\sigma^2_t,\sigma^2$.
  - `gate_rst_ptr__tw79.csv`: the RDST for each flow.
  - `gate_rst_pooled__tw79.csv`: an equal-weight pooled RDST. The paper reports the intensity-weighted pooled test in `table3__tw79.csv`.
- The remaining scripts write `table3`, `integration_pvals`, `pt_by_model`, `weather_sweep`, `weather_sweep_m0m1`, `ammi_necessity`, `ammi_necessity_strip`, `ammi_necessity_strip_m2`, `m0_pooled`, `decomp_m1_cells_dow` and `decomp_m2_blups_dow` (all `__tw79.csv`).

`gen_tables.py` writes `tab_m2_full.tex` as a single table. The manuscript shows its fixed-effect and variance-component rows as Table 8 and the flow-level rows in Figure 4. Table and figure numbers follow the TRB 2027 manuscript and may change in the journal version.

## License

MIT, see `LICENSE`.
