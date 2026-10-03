# -*- coding: utf-8 -*-
"""AMMI-necessity grid (M1.1.0-M1.1.4) + M0 pooled DiD, on the am7-9/pm17-19 window panel.
For each mode fit M1 = ring_flow FE + ysdp FE + 4 z-weather + 15 ring x peak DiD cells + AMMI_variant,
for each AMMI variant; take the 15 treated-cell DiD coefficients and report their spread in DiD%
(sd, min, max) and the count beyond +/-60%. Only the matched ring x peak AMMI should collapse them.
Also fit M0 = pooled TWFE (single PostTreated, no cells, no AMMI) per mode. Pure OLS.
"""
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
PAN = HERE.parent / "data" / "panels"
COV = HERE.parent / "data" / "covariates"
OUT = HERE.parent / "results"; OUT.mkdir(parents=True, exist_ok=True)
MODES = ["yellow", "uber", "lyft", "citibike", "subway"]
# adopted 3-variable weather set: panel-level temp/humidity + cell-level gusts
WVARS = ["temperature_2m", "wind_gusts_10m", "relative_humidity_2m"]
TR_FLOWS = ["0->0", "0->1", "0->2", "1->0", "2->0"]; PEAKS = ["am_peak", "off_peak", "pm_peak"]
SEASON3 = {"winter": "winter", "summer": "summer", "spring": "rest", "fall": "rest"}
RPCELLS = [(f, p) for f in TR_FLOWS for p in PEAKS]                    # 15 ring x peak treated cells
# AMMI variants: (label, rowdim, coldims) ; None = no AMMI
VARIANTS = {"M1.1.0_none": None,
            "M1.1.1_ring_season": ("ring_flow", ("season3",)),
            "M1.1.2_ring_dow": ("ring_flow", ("dow",)),
            "M1.1.3_peak_season": ("peak", ("season3",)),
            "M1.1.4_ring_peak": ("ring_flow", ("peak",))}


def zscore(x):
    m = x.mean(); s = np.sqrt(np.mean((x - m) ** 2)); return np.zeros_like(x) if s == 0 else (x - m) / s


def rank1_ammi(pre, dt, rowdim, coldims):
    kr_p = pre[rowdim].astype(str); kc_p = pre[list(coldims)].astype(str).agg("|".join, axis=1)
    kr_d = dt[rowdim].astype(str); kc_d = dt[list(coldims)].astype(str).agg("|".join, axis=1)
    cm = pre.assign(_r=kr_p, _c=kc_p).groupby(["_r", "_c"])["Y"].mean().unstack("_c")
    M = cm.to_numpy(); M = np.where(np.isnan(M), np.nanmean(M), M)
    Mc = M - M.mean(1, keepdims=True) - M.mean(0, keepdims=True) + M.mean()
    U, d, Vt = np.linalg.svd(Mc, full_matrices=False); sr, sp = U[:, 0] * np.sqrt(d[0]), Vt[0, :] * np.sqrt(d[0])
    if sp[np.argmax(np.abs(sp))] < 0: sr, sp = -sr, -sp
    return kr_d.map(dict(zip(cm.index, sr))).to_numpy() * kc_d.map(dict(zip(cm.columns, sp))).to_numpy()


def fit_cells(dt, ammi_col):
    n = len(dt)
    base = [np.ones(n)] + [zscore(dt[w].to_numpy()) for w in WVARS]
    base.append(pd.get_dummies(dt["ring_flow"], drop_first=True).to_numpy(float))
    base.append(pd.get_dummies(dt["ysdp"], drop_first=True).to_numpy(float))
    if ammi_col is not None: base.append(ammi_col.reshape(-1, 1))
    didX = np.column_stack([((dt["PostTreated"] == 1) & (dt["ring_flow"] == f) & (dt["peak"] == p)).to_numpy(float)
                            for (f, p) in RPCELLS])
    X = np.column_stack(base + [didX]); d0 = X.shape[1] - didX.shape[1]
    b, _, _, _ = np.linalg.lstsq(X, dt["Y"].to_numpy(), rcond=None)
    return b[d0:]                                                       # 15 ring x peak DiD coefs


rs_cov = pd.read_csv(COV / "rainsnow_cov.csv")
gu_cov = pd.read_csv(COV / "gusts_cov.csv")
nec_rows = []; m0_rows = []; strip_rows = []
for mode in MODES:
    dt = pd.read_parquet(PAN / f"ringpair_{mode}.parquet")
    dt = dt.merge(rs_cov[rs_cov["mode"] == mode].drop(columns="mode"),
                  on=["ring_flow", "year", "season", "peak"], how="left")
    dt = dt.merge(gu_cov[gu_cov["mode"] == mode].drop(columns="mode"),
                  on=["ring_flow", "year", "season", "peak"], how="left")
    assert dt[["rain", "snowfall", "wind_gusts_10m"]].notna().all().all()
    dt["season3"] = dt["season"].map(SEASON3); dt["Y"] = dt["log1p_avg_daily_rides"].astype(float)
    dt["PostTreated"] = ((dt["post"]) & (dt["treated"])).astype(int)
    dt["ysdp"] = dt["year"].astype(str) + "_" + dt["season"].astype(str) + "_" + dt["dow"].astype(str) + "_" + dt["peak"].astype(str)
    pre = dt[dt["year"] == 2024]
    # M0 pooled TWFE
    n = len(dt); Xb = [np.ones(n)] + [zscore(dt[w].to_numpy()) for w in WVARS]
    Xb.append(pd.get_dummies(dt["ring_flow"], drop_first=True).to_numpy(float))
    Xb.append(pd.get_dummies(dt["ysdp"], drop_first=True).to_numpy(float))
    Xb.append(dt["PostTreated"].to_numpy(float).reshape(-1, 1))
    Xm0 = np.column_stack(Xb); bm0, _, _, _ = np.linalg.lstsq(Xm0, dt["Y"].to_numpy(), rcond=None)
    m0_rows.append(dict(mode=mode, m0_beta=bm0[-1], m0_did_pct=100 * (np.exp(bm0[-1]) - 1)))
    # AMMI necessity
    for lab, spec in VARIANTS.items():
        ac = None if spec is None else rank1_ammi(pre, dt, spec[0], spec[1])
        cells = fit_cells(dt, ac); pct = 100 * (np.exp(cells) - 1)
        nec_rows.append(dict(mode=mode, variant=lab, cell_sd=np.std(pct), cell_min=pct.min(),
                             cell_max=pct.max(), n_beyond60=int((np.abs(pct) > 60).sum())))
        for (f, p), v in zip(RPCELLS, pct):
            strip_rows.append(dict(mode=mode, variant=lab, ring_flow=f, peak=p, pct=v))

nec = pd.DataFrame(nec_rows)
m0 = pd.DataFrame(m0_rows)
pd.DataFrame(strip_rows).to_csv(OUT / "ammi_necessity_strip__tw79.csv", index=False)
# summary: mean cell_sd over modes per variant, and total n_beyond60
summ = nec.groupby("variant").agg(mean_cell_sd=("cell_sd", "mean"), tot_beyond60=("n_beyond60", "sum")).round(2)
nec.round(2).to_csv(OUT / "ammi_necessity__tw79.csv", index=False)
m0.round(4).to_csv(OUT / "m0_pooled__tw79.csv", index=False)
pd.set_option("display.width", 200)
print("=== AMMI necessity: mean cell-effect sd (DiD%) over 5 modes, and #cells beyond +/-60% (of 75) ===")
print(summ.to_string())
print("\n=== M0 pooled TWFE DiD by mode ===")
print(m0.round(2).to_string(index=False))
