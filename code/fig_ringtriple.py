"""F9 -- flow-level DiD (%) with 95% CI along Ring 0 -> Ring 1 -> Ring 2,
arranged as three rows of commute-direction groups (peak-direction / off-peak /
reverse peak-direction) x five modes. Colour and marker encode the direction of
the flow (inbound / outbound / within-CRZ); solid = 95% CI excludes zero, hollow
= includes zero. RDST-failing flows are omitted (as in the cordon figure).
Inputs: results/{did_blup,gate_rst_ptr}__tw79.csv"""
import os, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import fig_style as st

st.setup()
did = pd.read_csv(os.path.join(st.RESULTS, "did_blup__tw79.csv"))
gate = pd.read_csv(os.path.join(st.RESULTS, "gate_rst_ptr__tw79.csv"))
gate["rst_pass"] = gate["rst_pass"].astype(str).str.upper().isin(["TRUE", "1", "T"])
d = did.merge(gate[["mode", "ring_flow", "peak", "rst_pass"]],
              on=["mode", "ring_flow", "peak"], how="left")
for c in ["blup", "ci_lo", "ci_hi"]:
    d[c + "_p"] = st.pct(d[c])
D = d.set_index(["mode", "ring_flow", "peak"])

# rows: commute-direction groups (flows given as (ring_flow, peak))
ROWS = [("peak-direction",
         [("0->0", "am_peak"), ("0->0", "pm_peak"),      # within-CRZ
          ("1->0", "am_peak"), ("0->1", "pm_peak"),      # Ring 1
          ("2->0", "am_peak"), ("0->2", "pm_peak")]),    # Ring 2
        ("off-peak",
         [("0->0", "off_peak"),
          ("1->0", "off_peak"), ("0->1", "off_peak"),
          ("2->0", "off_peak"), ("0->2", "off_peak")]),
        ("reverse peak-direction",
         [("1->0", "pm_peak"), ("0->1", "am_peak"),
          ("2->0", "pm_peak"), ("0->2", "am_peak")])]

POS_X = {"0->0": 0, "0->1": 1, "1->0": 1, "0->2": 2, "2->0": 2}
DIR_OF = {"0->0": "within", "0->1": "out", "0->2": "out",
          "1->0": "in", "2->0": "in"}
# green is reserved in this figure for the "fades outward" highlight, so the
# direction trio uses sky / orange / purple instead of the default CAT3 green
DIR_COL = {"in": "#56B4E9", "out": "#E69F00", "within": "#CC79A7"}
MARK = {"in": "^", "out": "v", "within": "o"}
DODGE = {"in": -0.17, "out": 0.17, "within": 0.0}

TOLLED = {"yellow", "uber", "lyft"}
FIT = "#444444"                          # non-significant ring slope
HL, HL_BG = "#1A9850", "#EAF7EC"         # effect fades outward
OPP, OPP_BG = "#C0392B", "#FDECEA"       # effect intensifies outward

def wls_slope(xs, ys, sds):
    """weighted least squares of y on x (w = 1/sd^2) -> (intercept, slope, p)"""
    x = np.asarray(xs, float); y = np.asarray(ys, float)
    w = 1.0 / np.asarray(sds, float) ** 2
    if len(x) < 2 or len(np.unique(x)) < 2:
        return None
    xb = np.sum(w * x) / np.sum(w); yb = np.sum(w * y) / np.sum(w)
    sxx = np.sum(w * (x - xb) ** 2)
    b = np.sum(w * (x - xb) * (y - yb)) / sxx
    se = np.sqrt(1.0 / sxx)
    z = b / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / np.sqrt(2))))
    return yb - b * xb, b, p

# height trimmed 4.6 -> 3.4 in to fit the TRB page budget; the width stays at
# FULLW so the figure is still embedded 1:1 and every label keeps its point size
fig, axes = plt.subplots(3, 5, figsize=(st.FULLW, 2.95), sharex=True)
for i, (gname, flows) in enumerate(ROWS):
    # rings the group is defined on (reverse peak-direction has no within-CRZ flow)
    rings_expected = {POS_X[f] for f, _ in flows}
    for j, mode in enumerate(st.MODES):
        ax = axes[i, j]
        ax.axhline(0, color=st.GRID, lw=0.8, zorder=1)
        fx, fy, fsd = [], [], []
        for f, pk in flows:
            r = D.loc[(mode, f, pk)]
            if not bool(r["rst_pass"]):        # RDST-failing flows omitted
                continue
            dirn = DIR_OF[f]
            x = POS_X[f] + DODGE[dirn]
            col = DIR_COL[dirn]
            sig = r["ci_lo"] > 0 or r["ci_hi"] < 0
            ax.plot([x, x], [r["ci_lo_p"], r["ci_hi_p"]], color=col, lw=1.0,
                    alpha=1.0 if sig else st.FADE + 0.25, zorder=2)
            ax.plot(x, r["blup_p"], MARK[dirn], ms=4.0, mec=col, mew=1.0,
                    mfc=col if sig else "white", zorder=3)
            fx.append(POS_X[f]); fy.append(r["blup_p"])
            fsd.append(100 * np.exp(r["blup"]) * r["sd"])   # delta-method SD in %
        res = wls_slope(fx, fy, fsd)
        if res is not None:
            a, b, p = res
            partial = set(fx) != rings_expected   # a ring lost to the RDST
            sig = p < 0.05 and not partial    # highlight only if no ring is missing
            # green when the effect fades outward: the tolled estimates are
            # negative, so fading is b > 0; the untolled are positive near the
            # zone, so fading is b < 0. Red marks the reverse.
            fades = (b > 0) if mode in TOLLED else (b < 0)
            if not sig:
                lc, bg, ls = FIT, None, (0, (4, 2))
            elif fades:
                lc, bg, ls = HL, HL_BG, "-"
            else:
                lc, bg, ls = OPP, OPP_BG, (0, (5, 1.2, 1, 1.2))
            xx = np.array([-0.35, 2.35]) if not partial else \
                np.array([min(fx) - 0.35, max(fx) + 0.35])
            ax.plot(xx, a + b * xx, color=lc, lw=1.6, alpha=0.35,
                    ls=ls, zorder=0)                            # behind the marks
            if bg is not None:
                ax.set_facecolor(bg)
                for sp in ax.spines.values():
                    sp.set_edgecolor(lc); sp.set_linewidth(1.4)
            if partial:
                ax.text(0.5, 0.06, f"{len(fx)} flows only",
                        transform=ax.transAxes, fontsize=6, color=st.MUTED,
                        ha="center", va="bottom")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["Ring 0", "Ring 1", "Ring 2"], fontsize=6.5)
        ax.set_xlim(-0.5, 2.5)
        ax.tick_params(axis="x", length=0)
        if i == 0:
            ax.set_title(st.MODE_LABEL[mode], fontsize=8, pad=3)
        if j == 0:
            ax.set_ylabel(gname, fontsize=7)
# plain-language label, matching the cordon figure: no formula in a figure
fig.supylabel("treatment effect per flow (% change)", fontsize=8, x=0.005)

leg = [Line2D([0], [0], marker=MARK[k], ls="none", mfc=DIR_COL[k],
              mec=DIR_COL[k], ms=4.5, label=lab)
       for k, lab in [("in", "inbound"), ("out", "outbound"),
                      ("within", "within-CRZ")]]
leg += [Line2D([0], [0], marker="o", ls="none", mfc="white", mec=st.MUTED,
               ms=4.5, label="CI includes 0")]
fig.legend(handles=leg, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, 0.028), handlelength=1.0, columnspacing=1.6)

leg2 = [Line2D([0], [0], color=HL, lw=1.6, alpha=0.35,
               label="suppression (tolled) / substitution (untolled) fades outward"),
        Line2D([0], [0], color=OPP, lw=1.6, alpha=0.35,
               ls=(0, (5, 1.2, 1, 1.2)), label="the reverse: intensifies outward"),
        Line2D([0], [0], color=FIT, lw=1.6, alpha=0.35, ls=(0, (4, 2)),
               label="not significant, or fewer than three flows")]
fig.legend(handles=leg2, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.032), handlelength=1.4, columnspacing=1.4,
           fontsize=6.2)
fig.subplots_adjust(left=0.115, right=0.995, top=0.94, bottom=0.145,
                    wspace=0.34, hspace=0.22)
st.save(fig, "fig_ringtriple")
