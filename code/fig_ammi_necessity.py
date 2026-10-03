"""AMMI covariate selection: the 75 flow-level DiD estimates (15 treated flows
x 5 modes) under the five candidate AMMI covariates, one horizontal strip each,
on the interacted fixed-effect model (left) and the random-slope model (right);
per-row sd annotated. The model numbering is not used: the rows are named by
the covariate itself. Inputs: results/ammi_necessity_strip__tw79.csv,
results/ammi_necessity_strip_m2__tw79.csv."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fig_style as st

st.setup()
m1 = pd.read_csv(os.path.join(st.RESULTS, "ammi_necessity_strip__tw79.csv"))
m2 = pd.read_csv(os.path.join(st.RESULTS, "ammi_necessity_strip_m2__tw79.csv"))
# bottom -> top; the matched AMMI (OD ring pair x peak window) is the adopted one
VAR = [("M1.1.3_peak_season", "peak_season", "AMMI $p{\\times}s$"),
       ("M1.1.2_ring_dow",    "ring_dow",    "AMMI $r{\\times}d$"),
       ("M1.1.1_ring_season", "ring_season", "AMMI $r{\\times}s$"),
       ("M1.1.0_none",        "none",        "no AMMI"),
       ("M1.1.4_ring_peak",   "ring_peak",   "AMMI $r{\\times}p$ (matched)")]
ADOPT_ROW = 4
rng = np.random.default_rng(7)

fig, axes = plt.subplots(1, 2, figsize=(st.FULLW, 2.3), sharey=True)
for ax, d, key, title in [(axes[0], m1, 0, "interacted fixed-effect model"),
                          (axes[1], m2, 1, "random-slope model")]:
    ax.axvline(0, color=st.GRID, lw=0.8, zorder=1)
    for y, spec in enumerate(VAR):
        vals = d[d["variant"] == spec[key]]["pct"].to_numpy()
        adopted = (y == ADOPT_ROW)
        col = st.CAT3[0] if adopted else st.MUTED
        jy = y + rng.uniform(-0.14, 0.14, len(vals))
        ax.plot(vals, jy, "o", ms=2.4, mfc=col, mec="none",
                alpha=0.85 if adopted else 0.5, zorder=3)
        ax.text(0.99, y + 0.30, f"sd = {vals.std(ddof=1):.1f}%",
                transform=ax.get_yaxis_transform(), fontsize=6,
                color=st.INK if adopted else st.MUTED, ha="right", va="center")
    ax.set_title(title, fontsize=7)
    ax.set_ylim(-0.5, len(VAR) - 0.5)
    ax.tick_params(axis="y", length=0)
axes[0].set_yticks(range(len(VAR)))
axes[0].set_yticklabels([spec[2] for spec in VAR], fontsize=6.5)
fig.supxlabel("treatment effect per flow (% change), 15 treated flows $\\times$ 5 modes",
              fontsize=7, y=0.02)
fig.subplots_adjust(left=0.185, right=0.995, top=0.90, bottom=0.22, wspace=0.06)
st.save(fig, "fig_ammi_necessity")
