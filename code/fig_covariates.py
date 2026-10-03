"""F7 -- weather and AMMI coefficients of M2. Weather (top row) = % change
exp(beta)-1 per +1 SD, shared axis across modes; AMMI (bottom row) = the raw
coefficient on the log scale, per +1 unit of the peak-direction score, on its
own axis. The AMMI score carries an arbitrary unit, so exp(beta)-1 would run
past +800% and could not be read against the weather terms (the two
scales differ by an order of magnitude). Solid = 95% CI excludes zero; hollow =
includes zero. Input: results/fixed__tw79.csv"""
import os
import pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_style as st

st.setup()
fx = pd.read_csv(os.path.join(st.RESULTS, "fixed__tw79.csv"))
WTERMS = [("z_temperature_2m",       "temperature"),
          ("z_wind_gusts_10m",       "wind gusts"),
          ("z_relative_humidity_2m", "relative humidity")]
fxi = fx.set_index(["mode", "term"])

def ci_pct(mode, term, pct=True):
    """percentage change for the weather terms, raw log coefficient for AMMI"""
    r = fxi.loc[(mode, term)]
    half = stats.t.ppf(0.975, r.df) * r.se
    vals = [r.est - half, r.est, r.est + half]
    lo, mid, hi = st.pct(vals) if pct else vals
    return lo, mid, hi, (lo > 0 or hi < 0)

def draw(ax, y, lo, mid, hi, sig):
    col = st.NEG if mid < 0 else st.POS
    ax.plot([lo, hi], [y, y], color=col, lw=1.1,
            alpha=1.0 if sig else st.FADE + 0.25, zorder=2,
            solid_capstyle="butt")
    ax.plot(mid, y, "o", ms=4.2, mec=col, mew=1.1,
            mfc=col if sig else "white", alpha=1.0 if sig else 0.85, zorder=3)

fig, axes = plt.subplots(2, 5, figsize=(st.FULLW, 2.4), sharex="row",
                         gridspec_kw={"height_ratios": [2.8, 1.0]})
ys = range(len(WTERMS))[::-1]
for j, mode in enumerate(st.MODES):
    axw, axa = axes[0, j], axes[1, j]
    for ax in (axw, axa):
        ax.axvline(0, color=st.GRID, lw=0.8, zorder=1)
    for y, (term, _) in zip(ys, WTERMS):
        draw(axw, y, *ci_pct(mode, term))
    draw(axa, 0, *ci_pct(mode, "ammi", pct=False))
    axw.set_title(st.MODE_LABEL[mode], fontsize=8, pad=2)
    axw.set_ylim(-0.6, len(WTERMS) - 0.4); axw.tick_params(labelbottom=True)
    axa.set_ylim(-0.6, 0.6); axa.set_yticks([])
    axw.tick_params(axis="x", labelsize=6.5)
    axa.tick_params(axis="x", labelsize=6.5)
    if j > 0:
        axw.set_yticks([])
axes[0, 0].set_yticks(list(ys))
axes[0, 0].set_yticklabels([lab for _, lab in WTERMS])
axes[0, 0].tick_params(axis="y", length=0)
axes[1, 0].set_yticks([0]); axes[1, 0].set_yticklabels(["AMMI"])
axes[1, 0].tick_params(axis="y", length=0)
axes[0, 2].set_xlabel("% change", fontsize=7)
axes[1, 2].set_xlabel("coefficient (log scale)", fontsize=7)
leg = [Line2D([0], [0], marker="o", ls="none", mfc=st.NEG, mec=st.NEG, ms=5,
              label="negative"),
       Line2D([0], [0], marker="o", ls="none", mfc=st.POS, mec=st.POS, ms=5,
              label="positive"),
       Line2D([0], [0], marker="o", ls="none", mfc="white", mec=st.MUTED, ms=5,
              label="CI includes 0")]
fig.legend(handles=leg, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.05), handlelength=1.0, columnspacing=1.4)
fig.subplots_adjust(left=0.13, right=0.995, top=0.90, bottom=0.24,
                    wspace=0.10, hspace=0.75)
st.save(fig, "fig_covariates")
