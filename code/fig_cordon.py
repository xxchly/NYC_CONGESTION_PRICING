"""F6 -- flow-level DiD (%) on the ring cordon, grouped by commute direction.
Faithful to the original paper_tw79 cordon figure (plot_cordon_tw.py): discrete
diverging colour bins, white-haloed arrows, per-mode counts line, horizontal
colourbar. Adapted: no map panel (moved to fig_ringmaps), three glyphs scaled
up proportionally, fonts raised to print-readable size, S->Sub / CB->Citi,
retired terms replaced by canonical ones (peak-direction / flow / RDST), and
the within-CRZ value lists sit directly below the core circle (they no longer
fit inside it at readable font size).
Added: pooled DiD markers exp(beta1)-1 per mode on the colourbar axis (solid =
95% CI excludes 0), so flow-level arrows can be read against the mode average.
Inputs: results/{did_blup,gate_rst_ptr,fitstats}__tw79.csv"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import TwoSlopeNorm, to_rgb, to_hex, ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch, Circle
import fig_style as st

st.setup()
ABBR = st.MODE_ABBR
W3 = np.ones(3)
INSIG = "#E9E9E6"

def tint(hexc, f):
    return to_hex(W3 - f * (W3 - np.array(to_rgb(hexc))))

# ---------- data ----------
did = pd.read_csv(os.path.join(st.RESULTS, "did_blup__tw79.csv"))
gate = pd.read_csv(os.path.join(st.RESULTS, "gate_rst_ptr__tw79.csv"))
gate["rst_pass"] = gate["rst_pass"].astype(str).str.upper().isin(["TRUE", "1", "T"])
did = did.merge(gate[["mode", "ring_flow", "peak", "rst_pass"]],
                on=["mode", "ring_flow", "peak"], how="left")
did["pct"] = st.pct(did["blup"])
did["sig"] = (did["ci_lo"] > 0) | (did["ci_hi"] < 0)
D = did.set_index(["mode", "ring_flow", "peak"])


# ---------- discrete bins (original) ----------
old = TwoSlopeNorm(vcenter=0.0, vmin=-40, vmax=40); rdbu = plt.get_cmap("RdBu")
def anchor(v): return to_hex(rdbu(old(v)))
NEG_BINS = [(-np.inf, -30, anchor(-40)), (-30, -25, anchor(-33)), (-25, -20, anchor(-28)),
            (-20, -15, anchor(-24)), (-15, -10, anchor(-20)), (-10, -5, anchor(-15)),
            (-5, 0, anchor(-10))]
POS_BINS = [(0, 5, anchor(15)), (5, 10, anchor(22)), (10, 15, anchor(30)), (15, np.inf, anchor(38))]
def cell_color(row):
    if not row["sig"]: return INSIG
    for lo, hi, c in (NEG_BINS + POS_BINS):
        if lo < row["pct"] <= hi: return c
    return INSIG
def fmt(row):
    v = f"{row['pct']:+.1f}"
    return f"{v}*" if not row["sig"] else v

# ---------- glyph geometry (original layout; Ring 0 enlarged so the within-CRZ
# value lists fit beside the loops inside the core at readable font size) ----------
R0, R1, R2 = 4.00, 5.15, 6.30
A1 = A2 = 0.63; SP = 0.74
mid1, mid2 = (R0 + R1) / 2, (R1 + R2) / 2
def pol(r, a): return np.array([r * np.cos(a), r * np.sin(a)])

def arrow(ax, row, p0, p1, rad, loop=False):
    kw = dict(connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
              mutation_scale=7 if loop else 7.5, lw=1.1, color=cell_color(row),
              zorder=6, shrinkA=0, shrinkB=0,
              path_effects=[pe.Stroke(linewidth=2.5, foreground="white"), pe.Normal()])
    ax.add_patch(FancyArrowPatch(tuple(p0), tuple(p1), **kw))

def fan(ax, flow, peak, ang, r_out, kind):
    u = np.array([np.cos(ang), np.sin(ang)]); p = np.array([-u[1], u[0]])
    if p[1] < 0: p = -p
    ends, lines = [], []
    for k, mode in enumerate(st.MODES):
        off_v = (2 - k) * SP * p
        row = D.loc[(mode, flow, peak)]
        if not bool(row["rst_pass"]):
            continue
        t0, t1 = pol(R0 * 0.92, ang) + off_v, pol(r_out, ang) + off_v
        if kind == "out": arrow(ax, row, t0, t1, 0.10)
        else: arrow(ax, row, t1, t0, 0.10)
        ends.append(t1); lines.append(f"{ABBR[mode]} {fmt(row)}")
    if u[1] > 0:      # upper fans: one label per arrow, at its outer end
        for (ex, ey), txt in zip(ends, lines):
            lx = ex + 0.25 if u[0] > 0 else ex - 0.25
            ax.text(lx, ey, txt, fontsize=6.0, color=st.INK,
                    ha="left" if u[0] > 0 else "right", va="center", zorder=7)
    else:             # lower fans: block hangs off the outermost arrowhead
        ex, ey = ends[-1] if u[0] > 0 else ends[0]
        ax.text(ex + (0.20 if u[0] > 0 else -0.20), ey - 0.22, "\n".join(lines),
                fontsize=6.0, color=st.INK, linespacing=1.2,
                ha="left" if u[0] > 0 else "right", va="top", zorder=7)

def loops(ax, peak, x0, header=None):
    """within-CRZ loop arrows with value labels beside them (original layout)"""
    if header:
        # two lines: the spelled-out window names are too wide for one line and
        # the two core columns would collide
        ax.text(x0 + 1.10, 2.62, header, fontsize=6.0, fontweight="bold",
                color=st.MUTED, ha="center", va="center", zorder=7,
                linespacing=1.05)
    for k, mode in enumerate(st.MODES):
        row = D.loc[(mode, "0->0", peak)]
        if not bool(row["rst_pass"]):
            continue
        ly = 1.60 - k * 0.80
        arrow(ax, row, (x0, ly + 0.22), (x0, ly - 0.22), 2.1, loop=True)
        ax.text(x0 + 0.32, ly, f"{ABBR[mode]} {fmt(row)}", fontsize=6.0,
                color=st.INK, ha="left", va="center", zorder=7)

def glyph(ax, spec):
    for r, hue, f in [(R2, st.RING_COLOR[2], 0.14), (R1, st.RING_COLOR[1], 0.13),
                      (R0, st.RING_COLOR[0], 0.15)]:
        ax.add_patch(Circle((0, 0), r, fc=tint(hue, f), ec="#B9B8B2", lw=0.7, zorder=1))
    ax.add_patch(Circle((0, 0), R0, fc="none", ec=st.INK, lw=1.0, zorder=3))
    (fL_up, fL_lo), (fR_up, fR_lo) = spec["left"], spec["right"]
    fan(ax, *fL_up, np.pi - A1, R1 - 0.10, "in")      # stop just inside Ring 1
    fan(ax, *fL_lo, np.pi + A2, mid2 + 0.13, "in")
    fan(ax, *fR_up, A1, R1 - 0.10, "out")
    fan(ax, *fR_lo, -A2, mid2 + 0.13, "out")
    if spec["core"] == "both":
        loops(ax, "am_peak", -2.75, header="morning\npeak")
        loops(ax, "pm_peak", 0.45, header="evening\npeak")
    elif spec["core"] == "off":
        loops(ax, "off_peak", -1.55)
    ax.text(-4.3, 5.75, spec["labL"], fontsize=6.0, fontweight="bold",
            color=st.MUTED, ha="center", va="center", zorder=7)
    ax.text(4.3, 5.75, spec["labR"], fontsize=6.0, fontweight="bold",
            color=st.MUTED, ha="center", va="center", zorder=7)
    ax.set_xlim(-7.6, 7.6); ax.set_ylim(-8.9, 6.8)
    ax.set_aspect("equal"); ax.set_axis_off()

GLYPHS = [
    dict(title="peak-direction group",
         left=(("1->0", "am_peak"), ("2->0", "am_peak")),
         right=(("0->1", "pm_peak"), ("0->2", "pm_peak")), core="both",
         labL="morning inbound", labR="evening outbound",
         within=[("morning peak", "am_peak"), ("evening peak", "pm_peak")]),
    dict(title="reverse peak-direction group",
         left=(("1->0", "pm_peak"), ("2->0", "pm_peak")),
         right=(("0->1", "am_peak"), ("0->2", "am_peak")), core="none",
         labL="evening inbound", labR="morning outbound", within=[]),
    dict(title="off-peak group",
         left=(("1->0", "off_peak"), ("2->0", "off_peak")),
         right=(("0->1", "off_peak"), ("0->2", "off_peak")), core="off",
         labL="off-peak inbound", labR="off-peak outbound",
         within=[("off-peak", "off_peak")]),
]

# ---------- assemble ----------
# height trimmed for the TRB page budget; width stays at FULLW so the figure
# is embedded 1:1 and every label keeps its point size
fig = plt.figure(figsize=(st.FULLW, 2.85))
gs = fig.add_gridspec(1, 3, wspace=0.02, left=0.005, right=0.995,
                      top=0.955, bottom=0.20)
for j, spec in enumerate(GLYPHS):
    ax = fig.add_subplot(gs[0, j])
    glyph(ax, spec)
    ax.set_title(spec["title"], fontsize=8, pad=2)

bin_colors = [c for _, _, c in NEG_BINS] + [c for _, _, c in POS_BINS]
bounds = [-35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20]
cmapD = ListedColormap(bin_colors); normD = BoundaryNorm(bounds, cmapD.N)
smD = ScalarMappable(norm=normD, cmap=cmapD)
cax = fig.add_axes([0.28, 0.105, 0.36, 0.034])
cb = fig.colorbar(smD, cax=cax, orientation="horizontal", spacing="uniform",
                  ticks=bounds)
cb.ax.set_xticklabels(["<−30", "−30", "−25", "−20", "−15",
                       "−10", "−5", "0", "+5", "+10", "+15", ">15"],
                      fontsize=5.5)
cb.ax.tick_params(length=2)
cb.outline.set_edgecolor("#B9B8B2")
cb.set_label("treatment effect per flow (% change)", fontsize=6)
axp = fig.add_axes([0.665, 0.104, 0.010, 0.034])
axp.add_patch(plt.Rectangle((0, 0), 1, 1, fc=INSIG, ec="#B9B8B2"))
axp.set_xlim(0, 1); axp.set_ylim(0, 1); axp.set_axis_off()
fig.text(0.680, 0.121, "*  CI includes 0", fontsize=5.5,
         color=st.MUTED, va="center", ha="left")

# ---------- pooled DiD markers on the colourbar axis ----------
# The colourbar has uniform bin spacing, so map % -> bin-fraction position.
fs = pd.read_csv(os.path.join(st.RESULTS, "fitstats__tw79.csv")).set_index("mode")
def cpos(v):
    v = max(min(v, bounds[-1] - 1e-9), bounds[0] + 1e-9)
    for i in range(len(bounds) - 1):
        if bounds[i] <= v <= bounds[i + 1]:
            return (i + (v - bounds[i]) / (bounds[i + 1] - bounds[i])) / (len(bounds) - 1)
    return 0.0
pool = []
for m in st.MODES:
    v = 100 * (np.exp(fs.loc[m, "beta2_pooled"]) - 1)
    sig = fs.loc[m, "b2_hi"] < 0 or fs.loc[m, "b2_lo"] > 0
    pool.append((m, v, cpos(v), sig))
pool.sort(key=lambda r: r[2])
# dodge labels: sweep left->right enforcing a minimum horizontal gap
MS = 0.16
lx = [p[2] for p in pool]
for i in range(1, len(lx)):
    lx[i] = max(lx[i], lx[i - 1] + MS)
lx = [x - (np.mean(lx) - np.mean([p[2] for p in pool])) for x in lx]  # recentre
axo = fig.add_axes([0.28, 0.140, 0.36, 0.082]); axo.set_xlim(0, 1)
axo.set_ylim(0, 1); axo.set_axis_off()
for (m, v, x, sig), xl in zip(pool, lx):
    axo.plot(x, 0.16, marker="v", ms=3.4, mec=st.INK, mew=0.8,
             mfc=st.INK if sig else "white", zorder=3, clip_on=False)
    if abs(xl - x) > 0.01:   # leader from dodged label to its marker
        axo.plot([xl, x], [0.52, 0.30], color=st.MUTED, lw=0.5, alpha=0.7,
                 zorder=2, clip_on=False)
    txt = f"{ABBR[m]} {v:+.1f}" if sig else f"{ABBR[m]} {v:+.1f}*"
    axo.text(xl, 0.58, txt, fontsize=5.2, color=st.INK, ha="center",
             va="bottom", zorder=4)
fig.text(0.272, 0.156, "pooled treatment effect", fontsize=5.2,
         color=st.MUTED, ha="right", va="center")

st.save(fig, "fig_cordon_peakdirection")   # band / arrow key lives in the caption
