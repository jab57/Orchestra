"""
Generate the two figures for Orchestra_research_paper.tex:
  Figure 1 (figure_pipeline_schematic.pdf) -- the validation pipeline architecture.
  Figure 2 (figure_corroboration_summary.pdf) -- corroborated vs uncorroborated
    OncoKB rate by threshold, focal genes vs negative controls (Table 1 data).

Run once; regenerate by re-running this script if the underlying numbers change.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9.5,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})

# ---------------------------------------------------------------------------
# Figure 1: pipeline schematic
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.6, 4.5))
ax.set_xlim(-1.95, 10.2)
ax.set_ylim(-0.15, 10)
ax.axis("off")

BOX_STYLE = dict(boxstyle="round,pad=0.35", linewidth=1.0)

def box(x, y, w, h, text, fc, ec="black", fontsize=8.3, weight="normal", tc="black"):
    b = FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                        boxstyle="round,pad=0.12", linewidth=1.0,
                        facecolor=fc, edgecolor=ec, zorder=2)
    ax.add_patch(b)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
             weight=weight, zorder=3, linespacing=1.3, color=tc)

def arrow(x1, y1, x2, y2, label=None, label_dx=0.0, label_dy=0.15, fontsize=7):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                         linewidth=1.1, color="#333333", zorder=1)
    ax.add_patch(a)
    if label:
        ax.text((x1 + x2) / 2 + label_dx, (y1 + y2) / 2 + label_dy, label,
                 ha="center", va="center", fontsize=fontsize, color="#333333",
                 style="italic")

# Row 1: input
box(1.4, 9.15, 2.4, 0.75, "Focal gene +\nTCGA cancer type", "#eef2f6")

# Row 2: RegNetAgents
box(1.4, 8.0, 2.6, 0.9, "RegNetAgents\ncompare_network_contexts", "#dbe7f2")
arrow(1.4, 8.775, 1.4, 8.45)

# Split into the three tiers RegNetAgents returns (disjoint sets).
# This paper analyzes only tumor-acquired; the other two are dead-ends here.
box(-0.35, 6.5, 1.95, 0.9, "Population-averaged-\nonly", "#f0f0f0",
    ec="#999999", tc="#777777", fontsize=7.7)
box(1.8, 6.5, 1.95, 0.9, "Conserved\nregulators", "#f0f0f0",
    ec="#999999", tc="#777777", fontsize=7.9)
box(4.0, 6.5, 2.3, 0.9, "Tumor-acquired\nregulators (this paper)", "#ffe9cc", fontsize=7.9)
ax.text(0.7, 5.92, "not analyzed in this paper", ha="center", va="top", fontsize=6.6,
        style="italic", color="#8a8a8a")
arrow(0.7, 7.55, -0.35, 6.95)
arrow(1.4, 7.55, 1.8, 6.95)
arrow(2.1, 7.55, 4.0, 6.95)

# Rank by MI edge weight, cap at 10 (this paper's selection mechanism)
box(4.0, 5.0, 4.0, 0.9, "Rank by ARACNe MI edge weight\n(query_network), cap at 10", "#ffe9cc")
arrow(4.0, 6.05, 4.0, 5.45)

# CASCADE 4 lightweight sources
box(4.0, 3.4, 4.6, 1.4,
    "CASCADE (4 independent sources)\nLINCS  •  super-enhancer\nDoRothEA TF  •  DepMap",
    "#e3f0e6")
arrow(4.0, 4.55, 4.0, 4.1)

# Corroboration count
box(4.0, 2.0, 2.8, 0.85, "corroboration_count\n(0-4 per candidate)", "#fff3cd")
arrow(4.0, 2.7, 4.0, 2.42)

# Right column: negative controls, run through the same pipeline (including ranking)
box(8.1, 6.5, 2.6, 1.3, "Negative controls\n(housekeeping,\nnon-driver genes)", "#f6e3e6")
ax.text(8.1, 5.65, "(same pipeline,\nincl. MI ranking)", ha="center", va="top", fontsize=6.8,
         style="italic", color="#555555")

# OncoKB ground truth + stats
box(6.8, 0.6, 6.2, 1.0,
    "OncoKB ground truth  $\\rightarrow$  Fisher's exact + permutation test\n"
    "(deduplicated, focal vs. negative controls, thresholds 1-3)",
    "#eef2f6", fontsize=8)
arrow(4.0, 1.575, 5.2, 0.95)
arrow(8.1, 5.85, 7.3, 1.1)

ax.text(5.0, 9.85, "Orchestra corroboration-validation pipeline",
        ha="center", fontsize=11, weight="bold")

fig.tight_layout()
fig.savefig("figure_pipeline_schematic.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 2: corroboration summary bar chart (Table 1 data, MI-weight-ranked
# candidate selection)
# ---------------------------------------------------------------------------
thresholds = ["$\\geq$1 / 4", "$\\geq$2 / 4", "$\\geq$3 / 4"]
focal_corrob = [27.8, 44.2, 0.0]
focal_uncorrob = [18.9, 21.5, 26.1]
focal_p = ["p=0.1281", "p=0.0028", "p=1.0000"]
focal_sig = [False, True, False]
control_corrob = [18.7, 24.4, 50.0]
control_uncorrob = [8.3, 13.7, 15.8]
control_p = ["p=0.0658", "p=0.0721", "p=0.2978"]
control_sig = [False, False, False]

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)

x = np.arange(3)
width = 0.32

for ax, corrob, uncorrob, pvals, sigs, title in (
    (axes[0], focal_corrob, focal_uncorrob, focal_p, focal_sig, "Focal genes (BRCA/COAD)"),
    (axes[1], control_corrob, control_uncorrob, control_p, control_sig, "Negative controls"),
):
    corrob_plot = [v if v is not None else 0 for v in corrob]
    uncorrob_plot = [v if v is not None else 0 for v in uncorrob]
    b1 = ax.bar(x - width / 2, corrob_plot, width, label="Corroborated",
                color="#3a7d5c", edgecolor="black", linewidth=0.6, zorder=3)
    b2 = ax.bar(x + width / 2, uncorrob_plot, width, label="Uncorroborated",
                color="#b8b3a3", edgecolor="black", linewidth=0.6, zorder=3)

    for i, (c, u, p, sig) in enumerate(zip(corrob, uncorrob, pvals, sigs)):
        if c is None:
            ax.text(i, 3, "insufficient\ndata", ha="center", va="bottom",
                     fontsize=7, style="italic", color="#666666")
            continue
        top = max(c, u) + 2.5
        color = "#1f5c3d" if sig else "#666666"
        weight = "bold" if sig else "normal"
        ax.text(i, top, p, ha="center", va="bottom", fontsize=7.3,
                 color=color, weight=weight)

    ax.set_xticks(x)
    ax.set_xticklabels(thresholds, fontsize=9)
    ax.set_title(title, fontsize=9.5, weight="bold")
    ax.set_ylim(0, 58)
    ax.grid(axis="y", linewidth=0.4, alpha=0.4, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].set_ylabel("OncoKB rate (\\%)", fontsize=9.5)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 1.04), fontsize=9)
fig.suptitle("")
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig("figure_corroboration_summary.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 3: BRCA/COAD vs STAD rate-gap comparison at the >=2 threshold, across
# the two panels (BRCA/COAD from RegNetAgents' published panel; STAD newly
# constructed on a separate patient cohort), focal genes and negative controls,
# MI-weight-ranked candidate selection throughout.
# ---------------------------------------------------------------------------
groups = ["BRCA/COAD\nfocal", "BRCA/COAD\ncontrols", "STAD\nfocal", "STAD\ncontrols"]
gaps = [20.3, 10.7, 38.7, 15.2]
sig = [True, False, True, False]  # focal genes significant in both panels; controls in neither
colors = ["#3a7d5c" if s else "#b8b3a3" for s in sig]
annotations = ["p=0.0091\n(sig.)", "p=0.0721\n(n.s.)", "p=0.0021\n(sig.)", "p=0.0930\n(n.s.)"]

fig, ax = plt.subplots(figsize=(6.5, 3.6))
x = np.arange(4)
bars = ax.bar(x, gaps, width=0.55, color=colors, edgecolor="black", linewidth=0.7, zorder=3)
ax.axhline(0, color="black", linewidth=0.9, zorder=2)

for i, (g, a) in enumerate(zip(gaps, annotations)):
    va = "bottom" if g >= 0 else "top"
    offset = 1.3 if g >= 0 else -1.3
    ax.text(i, g + offset, a, ha="center", va=va, fontsize=7.6,
             color="#1f5c3d" if sig[i] else "#555555",
             weight="bold" if sig[i] else "normal", linespacing=1.3)

ax.set_xticks(x)
ax.set_xticklabels(groups, fontsize=9)
ax.set_ylabel("Corroborated $-$ uncorroborated\nOncoKB rate (percentage points)", fontsize=9)
ax.set_ylim(0, 44)
ax.grid(axis="y", linewidth=0.4, alpha=0.4, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_title("Corroboration signal replicates across two panels (BRCA/COAD and STAD)",
             fontsize=9.8, weight="bold")

fig.tight_layout()
fig.savefig("figure_stad_comparison.pdf", bbox_inches="tight")
plt.close(fig)

print("All three figures written.")
