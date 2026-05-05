#!/usr/bin/env python3
"""Generate the Orchestra architecture figure for the JOSS paper."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 8)
ax.axis('off')

# Color palette — mirrors CASCADE palette where roles overlap
C_CLIENT   = '#4A90D9'   # blue  — Claude Desktop
C_SERVER   = '#2C3E50'   # dark navy — server border/title
C_ROUTE    = '#F39C12'   # orange — decision layer
C_RNA      = '#8E44AD'   # purple — RegNetAgents
C_CASCADE  = '#E74C3C'   # red    — CASCADE child server
C_PARALLEL = '#EAFAF1'   # light green bg — evidence layer
C_SYNTH    = '#1ABC9C'   # teal   — synthesis
C_REPORT   = '#27AE60'   # green  — report
C_DARK     = '#2C3E50'


def box(x, y, w, h, color, label, fontsize: float = 8, textcolor='white',
        alpha=1.0, style='round,pad=0.1', edgecolor='#34495E',
        linestyle='solid', lw=1.2):
    fancy = FancyBboxPatch((x, y), w, h, boxstyle=style,
                           facecolor=color, edgecolor=edgecolor,
                           linewidth=lw, alpha=alpha, zorder=2,
                           linestyle=linestyle)
    ax.add_patch(fancy)
    ax.text(x + w / 2, y + h / 2, label, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', color=textcolor, zorder=3)


def arrow(x1, y1, x2, y2, color='#7F8C8D', style='-|>', lw=1.5,
          linestyle='solid'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw,
                                shrinkA=2, shrinkB=2, linestyle=linestyle),
                zorder=4)


# ── Row 0: Claude Desktop (outside Orchestra box) ──────────────────────────
box(3.0, 7.35, 4.0, 0.5, C_CLIENT, 'Claude Desktop  (MCP Client)', fontsize=9)
arrow(5.0, 7.35, 5.0, 7.02, color=C_DARK, lw=2.0)

# ── Orchestra outer box ─────────────────────────────────────────────────────
outer = FancyBboxPatch((0.2, 0.4), 9.6, 6.62, boxstyle='round,pad=0.15',
                       facecolor='#F8F9FA', edgecolor=C_SERVER,
                       linewidth=2.0, alpha=0.9, zorder=0)
ax.add_patch(outer)
ax.text(5.0, 6.88, 'Orchestra  —  LangGraph MCP Server', ha='center',
        va='center', fontsize=11, fontweight='bold', color=C_DARK, zorder=3)
ax.text(5.0, 6.63, 'MCP Server to Claude Desktop  ·  MCP Client to '
        'RegNetAgents & CASCADE  ·  3 Composite Tools',
        ha='center', va='center', fontsize=7, color='#5D6D7E',
        zorder=3, style='italic')

# ── Decision Layer ───────────────────────────────────────────────────────────
dec_bg = FancyBboxPatch((0.4, 5.45), 9.2, 1.0, boxstyle='round,pad=0.08',
                         facecolor='#FEF9E7', edgecolor='#F39C12',
                         linewidth=1.2, alpha=0.7, zorder=1)
ax.add_patch(dec_bg)
ax.text(0.65, 6.32, 'Decision Layer', ha='left', va='center',
        fontsize=8, fontweight='bold', color='#D35400', zorder=3)

box(0.6, 5.55, 2.5, 0.55, C_ROUTE, 'Classify Gene\n(CASCADE metadata)', fontsize=7)
arrow(3.1, 5.825, 3.45, 5.825, color=C_ROUTE, lw=1.8)
box(3.45, 5.55, 2.1, 0.55, C_ROUTE, 'Route Analysis\n(gene role)', fontsize=7)

# Three routing outputs
for yi, label in zip([6.17, 5.88, 5.60],
                     ['TF path', 'effector path', 'validation path']):
    ax.text(5.85, yi, label, fontsize=6.5, color='#D35400',
            fontweight='bold', zorder=3)
arrow(5.55, 5.825, 5.83, 6.12, color=C_ROUTE, lw=1.2)
arrow(5.55, 5.825, 5.83, 5.85, color=C_ROUTE, lw=1.2)
arrow(5.55, 5.825, 5.83, 5.58, color=C_ROUTE, lw=1.2)

# ── Evidence Layer ───────────────────────────────────────────────────────────
ev_bg = FancyBboxPatch((0.4, 3.75), 9.2, 1.55, boxstyle='round,pad=0.08',
                        facecolor=C_PARALLEL, edgecolor='#27AE60',
                        linewidth=1.2, alpha=0.6, zorder=1, linestyle='dashed')
ax.add_patch(ev_bg)
ax.text(0.65, 5.22, 'Evidence Layer  —  Parallel MCP Calls', ha='left',
        va='center', fontsize=8, fontweight='bold', color='#1E8449', zorder=3)

# Column headers
ax.text(2.35, 5.1, 'RegNetAgents', ha='center', va='center',
        fontsize=8, fontweight='bold', color=C_RNA, zorder=3)
ax.text(7.65, 5.1, 'CASCADE', ha='center', va='center',
        fontsize=8, fontweight='bold', color=C_CASCADE, zorder=3)

# RegNetAgents call boxes (left column)
box(0.5, 4.45, 3.7, 0.45, C_RNA, 'comprehensive_gene_analysis(gene, cell_type)', fontsize=6.5)
box(0.5, 3.88, 3.7, 0.45, C_RNA,
    'pathway_focused_analysis  ·  query_network(gene, cell_type)', fontsize=6.2)

# CASCADE call boxes (right column)
box(5.8, 4.45, 3.7, 0.45, C_CASCADE,
    'comprehensive_perturbation_analysis(gene, cell_type)', fontsize=6.2)
box(5.8, 3.88, 3.7, 0.45, C_CASCADE,
    'get_gene_metadata  ·  get_protein_interactions  ·  therapeutic_target_discovery',
    fontsize=5.7)

# Fork arrows: Route → each evidence column
arrow(5.0, 5.45, 2.35, 4.9, color=C_ROUTE, lw=1.5)
arrow(5.0, 5.45, 7.65, 4.9, color=C_ROUTE, lw=1.5)

# Evidence → Synthesize: both columns feed into Synthesize (teal), which then
# passes to Generate Report via the horizontal arrow already drawn.
arrow(2.35, 3.88, 2.45, 2.92, color='#27AE60', lw=1.5)
arrow(7.65, 3.88, 3.5,  2.92, color='#27AE60', lw=1.5)

# ── Explanation Layer ────────────────────────────────────────────────────────
ex_bg = FancyBboxPatch((0.4, 2.12), 9.2, 1.45, boxstyle='round,pad=0.08',
                        facecolor='#E8F8F5', edgecolor='#1ABC9C',
                        linewidth=1.2, alpha=0.7, zorder=1)
ax.add_patch(ex_bg)
ax.text(0.65, 3.45, 'Explanation Layer', ha='left', va='center',
        fontsize=8, fontweight='bold', color='#0E6655', zorder=3)

box(0.7, 2.22, 3.5, 0.7, C_SYNTH,
    'Synthesize\ncross-system corroboration count', fontsize=7.5)
arrow(4.2, 2.57, 4.6, 2.57, color=C_SYNTH, lw=1.8)
box(4.6, 2.22, 4.55, 0.7, C_REPORT,
    'Generate Report\nstructured evidence table  +  optional LLM narrative',
    fontsize=7)

# ── Child MCP Servers (bottom, inside outer box, dashed) ─────────────────────
# RegNetAgents server — purple dashed
box(0.35, 0.5, 3.8, 1.45, C_RNA,
    'RegNetAgents MCP Server\n\nnetwork analysis  ·  pathway enrichment\n'
    'PageRank centrality  ·  domain agents\n(Bird et al. 2026a)',
    fontsize=6.5, linestyle='dashed', edgecolor=C_RNA, lw=1.5)

# CASCADE server — red dashed
box(5.85, 0.5, 3.8, 1.45, C_CASCADE,
    'CASCADE MCP Server\n\nperturbation simulation  ·  STRING PPI\n'
    'LINCS L1000  ·  DepMap CRISPR  ·  super-enhancers\n'
    'DoRothEA  ·  cBioPortal  (Bird et al. 2026b)',
    fontsize=6.0, linestyle='dashed', edgecolor=C_CASCADE, lw=1.5)

# Arrows from child servers UP to evidence layer.
# x=0.4 is left of Synthesize (x=0.7); x=9.25 is right of Report (ends at 9.15).
# Both paths clear the Explanation layer boxes.
arrow(0.42, 1.95, 0.42, 3.75, color=C_RNA, lw=1.5, linestyle='dashed')
arrow(9.25, 1.95, 9.25, 3.75, color=C_CASCADE, lw=1.5)

# ── Legend ───────────────────────────────────────────────────────────────────
legend_items = [
    mpatches.Patch(color=C_ROUTE,   label='Decision'),
    mpatches.Patch(color=C_RNA,     label='RegNetAgents'),
    mpatches.Patch(color=C_CASCADE, label='CASCADE'),
    mpatches.Patch(color=C_SYNTH,   label='Synthesis'),
    mpatches.Patch(color=C_REPORT,  label='Report'),
    Line2D([0], [0], color=C_RNA,     lw=1.5, linestyle='dashed',
           label='MCP call (RegNetAgents)'),
    Line2D([0], [0], color=C_CASCADE, lw=1.5, linestyle='solid',
           label='MCP call (CASCADE)'),
]
ax.legend(handles=legend_items, loc='lower center', ncol=7, fontsize=6.8,
          frameon=True, fancybox=True, framealpha=0.9,
          bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.savefig('figure_architecture.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('figure_architecture.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("Saved figure_architecture.png and figure_architecture.pdf")
