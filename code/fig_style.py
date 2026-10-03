"""
Shared figure style for paper_tw79 deliverables.
Role-based palette (fixed for the whole paper):
  NEG  blue      = negative DiD (suppression)        } the only diverging pair
  POS  vermilion = positive DiD (substitution)       }
  CAT3           = per-figure categorical trio (legend defines meaning)
  ring grays     = Ring 0/1/2 (ordered -> lightness, not hue)
Significance is encoded by fill (solid = 95% CI excludes zero, hollow/faded =
includes zero), never by colour alone. All text uses ink tokens, sans-serif.
Terminology follows the paper (canonical mode/peak/flow names).
"""
import os
import numpy as np
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
RESULTS = os.path.join(PAPER, "results")
FIGDIR = os.path.join(PAPER, "figures")
TABDIR = os.path.join(PAPER, "tables")
os.makedirs(TABDIR, exist_ok=True)

# ---------------- canonical names & orders (as in the paper) ----------------
MODES = ["yellow", "uber", "lyft", "subway", "citibike"]          # display order
MODE_LABEL = {"yellow": "Yellow Taxi", "uber": "Uber", "lyft": "Lyft",
              "subway": "Subway", "citibike": "Citi Bike"}
# Mode-name rule: within one figure use EITHER full names (MODE_LABEL) OR the
# abbreviation set Y / U / L / Sub / Citi -- never a mix.
MODE_ABBR = {"yellow": "Y", "uber": "U", "lyft": "L",
             "subway": "Sub", "citibike": "Citi"}
PEAKS = ["am_peak", "off_peak", "pm_peak"]
PEAK_LABEL = {"am_peak": "morning peak", "off_peak": "off-peak",
              "pm_peak": "evening peak"}
# spelled out everywhere: the manuscript does not use AM / PM / off anywhere,
# in text, tables or figures
PEAK_SHORT = {"am_peak": "morning peak", "off_peak": "off-peak",
              "pm_peak": "evening peak"}
TR_FLOWS = ["0->0", "0->1", "0->2", "1->0", "2->0"]

def rp(s):
    """ring-pair display: '0->1' -> '0→1' (arrows denote OD ring pairs only)"""
    return s.replace("->", "→")

# commute-direction grouping of the 15 treated flows (paper, Figure 4)
GROUPS = ["peak-direction", "reverse peak-direction", "off-peak"]
GROUP_OF = {}
for _f in TR_FLOWS:
    GROUP_OF[(_f, "off_peak")] = "off-peak"
for _f, _p in [("1->0", "am_peak"), ("2->0", "am_peak"),          # morning-inbound
               ("0->1", "pm_peak"), ("0->2", "pm_peak"),          # evening-outbound
               ("0->0", "am_peak"), ("0->0", "pm_peak")]:         # within-CRZ
    GROUP_OF[(_f, _p)] = "peak-direction"
for _f, _p in [("0->1", "am_peak"), ("0->2", "am_peak"),          # morning-outbound
               ("1->0", "pm_peak"), ("2->0", "pm_peak")]:         # evening-inbound
    GROUP_OF[(_f, _p)] = "reverse peak-direction"

# ---------------- palette (Okabe-Ito, colour-blind safe) ----------------
# sign convention follows the cordon figure's diverging scale: red-family =
# negative (suppression), blue-family = positive (substitution)
NEG = "#D55E00"      # suppression / negative DiD
POS = "#0072B2"      # substitution / positive DiD
CAT3 = ["#009E73", "#E69F00", "#CC79A7"]   # categorical trio (green/orange/purple)
CAT4 = CAT3 + ["#999999"]                  # trio + neutral (e.g. VPC 4th slot)
RING_GRAY = {0: "#555555", 1: "#9C9C9C", 2: "#D8D8D8"}
RING_COLOR = {0: "#5E3C99", 1: "#2C7FB8", 2: "#E08214"}  # ring identity (maps + cordon bands)
INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
FADE = 0.30          # alpha for CI-includes-zero marks

FULLW = 6.5          # in — full text width (8.5x11, 1in margins)

def pct(beta):
    """percentage-change transform, exp(beta)-1, in %"""
    return (np.exp(np.asarray(beta, dtype=float)) - 1.0) * 100.0

def setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
        "lines.linewidth": 1.2, "patch.linewidth": 0.5,
        "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
        "figure.dpi": 120, "savefig.dpi": 300,
    })

def save(fig, name):
    """write figures/<name>.pdf (vector, for LaTeX) + .png (300 dpi proof)"""
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"),
                    bbox_inches="tight", pad_inches=0.02)
    print(f"saved figures/{name}.pdf/.png")

# quick grayscale-separation self check (run: python fig_style.py)
if __name__ == "__main__":
    def lum(hx):
        r, g, b = (int(hx[i:i+2], 16)/255 for i in (1, 3, 5))
        return 0.2126*r + 0.7152*g + 0.0722*b
    for name, cols in [("diverging", [NEG, POS]), ("cat3", CAT3),
                       ("rings", list(RING_GRAY.values()))]:
        L = [round(lum(c), 3) for c in cols]
        print(f"{name}: {cols} luminance {L}")
