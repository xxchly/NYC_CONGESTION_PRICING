# -*- coding: utf-8 -*-
"""E2 grid decomposition with day-of-week, M1 side.
Fit M1 = 225 FIXED ring_flow x season3 x peak x dow DiD cells (+ ring_flow FE
+ ysdp FE + 3 z-weather + composite rank-1 AMMI), save the 225 cell effects
per mode. Pure OLS. Output: results/decomp_m1_cells_dow__tw79.csv
"""
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
PAN = HERE.parent / "data" / "panels"
COV = HERE.parent / "data" / "covariates"
OUT = HERE.parent / "results"; OUT.mkdir(parents=True, exist_ok=True)
MODES = ["yellow", "uber", "lyft", "citibike", "subway"]
WVARS = ["temperature_2m", "wind_gusts_10m", "relative_humidity_2m"]
FLOWS = ["0->0", "0->1", "0->2", "1->0", "2->0"]; PEAKS = ["am_peak", "off_peak", "pm_peak"]
SEA3 = ["rest", "summer", "winter"]
SEASON3 = {"winter": "winter", "summer": "summer", "spring": "rest", "fall": "rest"}


def zscore(x):
    m = x.mean(); s = np.sqrt(np.mean((x - m) ** 2)); return np.zeros_like(x) if s == 0 else (x - m) / s


def rank1_ammi(pre, dt, coldims):
    kp = pre[list(coldims)].astype(str).agg("|".join, axis=1); kd = dt[list(coldims)].astype(str).agg("|".join, axis=1)
    cm = pre.assign(_k=kp).groupby(["ring_flow", "_k"])["Y"].mean().unstack("_k")
    M = cm.to_numpy(); M = np.where(np.isnan(M), np.nanmean(M), M)
    Mc = M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + M.mean()
    U, d, Vt = np.linalg.svd(Mc, full_matrices=False); sr, sp = U[:, 0] * np.sqrt(d[0]), Vt[0, :] * np.sqrt(d[0])
    if sp[np.argmax(np.abs(sp))] < 0: sr, sp = -sr, -sp
    return dt["ring_flow"].map(dict(zip(cm.index, sr))).to_numpy() * kd.map(dict(zip(cm.columns, sp))).to_numpy()


gu_cov = pd.read_csv(COV / "gusts_cov.csv")
rows = []
for mode in MODES:
    dt = pd.read_parquet(PAN / f"ringpair_{mode}.parquet")
    dt = dt.merge(gu_cov[gu_cov["mode"] == mode].drop(columns="mode"),
                  on=["ring_flow", "year", "season", "peak"], how="left")
    assert dt["wind_gusts_10m"].notna().all()
    dt["season3"] = dt["season"].map(SEASON3); dt["Y"] = dt["log1p_avg_daily_rides"].astype(float)
    dt["PostTreated"] = ((dt["post"]) & (dt["treated"])).astype(int)
    dt["ysdp"] = dt["year"].astype(str) + "_" + dt["season"].astype(str) + "_" + dt["dow"].astype(str) + "_" + dt["peak"].astype(str)
    pre = dt[dt["year"] == 2024]; n = len(dt)
    DOWS = sorted(dt["dow"].unique())
    cells = [(f, s, p, d) for f in FLOWS for s in SEA3 for p in PEAKS for d in DOWS]
    base = [np.ones(n)] + [zscore(dt[w].to_numpy()) for w in WVARS]
    base.append(pd.get_dummies(dt["ring_flow"], drop_first=True).to_numpy(float))
    base.append(pd.get_dummies(dt["ysdp"], drop_first=True).to_numpy(float))
    ammi = rank1_ammi(pre, dt, ("season3", "peak", "dow"))                 # composite AMMI nuisance
    didX = np.column_stack([((dt["PostTreated"] == 1) & (dt["ring_flow"] == f) & (dt["season3"] == s)
                             & (dt["peak"] == p) & (dt["dow"] == d)).to_numpy(float)
                            for (f, s, p, d) in cells])
    assert (didX.sum(0) > 0).all(), "empty DiD cell"
    X = np.column_stack(base + [ammi, didX]); d0 = X.shape[1] - didX.shape[1]
    b, _, rank, _ = np.linalg.lstsq(X, dt["Y"].to_numpy(), rcond=None)
    assert rank == X.shape[1], f"rank deficient: {rank} < {X.shape[1]}"
    for (f, s, p, d), v in zip(cells, b[d0:]):
        rows.append(dict(mode=mode, ring_flow=f, season3=s, peak=p, dow=d, coef=v))
    print(mode, "done")

pd.DataFrame(rows).to_csv(OUT / "decomp_m1_cells_dow__tw79.csv", index=False)
print("DONE")
