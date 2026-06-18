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
C_PUBMED   = '#117A65'   # dark teal — PubMed/NCBI
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


# ── Row 0: Claude Desktop (outside Orchestra box) ───────────────────────────
box(3.0, 7.35, 4.0, 0.5, C_CLIENT, 'Claude Desktop  (MCP Client)', fontsize=9)
arrow(5.0, 7.35, 5.0, 7.02, color=C_DARK, lw=2.0)

# ── Orchestra outer box ──────────────────────────────────────────────────────
outer = FancyBboxPatch((0.2, 0.4), 9.6, 6.62, boxstyle='round,pad=0.15',
                       facecolor='#F8F9FA', edgecolor=C_SERVER,
                       linewidth=2.0, alpha=0.9, zorder=0)
ax.add_patch(outer)
ax.text(5.0, 6.88, 'Orchestra  —  LangGraph MCP Server', ha='center',
        va='center', fontsize=11, fontweight='bold', color=C_DARK, zorder=3)
ax.text(5.0, 6.63,
        'MCP Server to Claude Desktop  ·  MCP Client to RegNetAgents & CASCADE'
        '  ·  NCBI E-utilities  ·  7 Composite Tools',
        ha='center', va='center', fontsize=7, color='#5D6D7E',
        zorder=3, style='italic')

# ── Decision Layer ───────────────────────────────────────────────────────────
dec_bg = FancyBboxPatch((0.4, 5.35), 9.2, 1.10, boxstyle='round,pad=0.08',
                         facecolor='#FEF9E7', edgecolor='#F39C12',
                         linewidth=1.2, alpha=0.7, zorder=1)
ax.add_patch(dec_bg)
ax.text(0.65, 6.32, 'Decision Layer', ha='left', va='center',
        fontsize=8, fontweight='bold', color='#D35400', zorder=3)

box(0.6, 5.48, 2.5, 0.55, C_ROUTE, 'Classify Gene\n(CASCADE metadata)', fontsize=7)
arrow(3.1, 5.755, 3.45, 5.755, color=C_ROUTE, lw=1.8)
box(3.45, 5.48, 2.1, 0.55, C_ROUTE, 'Route Analysis\n(gene role)', fontsize=7)

# Four routing outputs
for yi, label in zip([6.22, 5.97, 5.72, 5.47],
                     ['TF path', 'effector path', 'validation path', 'novelty path']):
    fc = C_PUBMED if label == 'novelty path' else '#D35400'
    ax.text(5.85, yi, label, fontsize=6.5, color=fc, fontweight='bold', zorder=3)
arrow(5.55, 5.755, 5.83, 6.22, color=C_ROUTE,  lw=1.2)
arrow(5.55, 5.755, 5.83, 5.97, color=C_ROUTE,  lw=1.2)
arrow(5.55, 5.755, 5.83, 5.72, color=C_ROUTE,  lw=1.2)
arrow(5.55, 5.755, 5.83, 5.47, color=C_PUBMED, lw=1.2, linestyle='dashed')

# ── Evidence Layer ───────────────────────────────────────────────────────────
ev_bg = FancyBboxPatch((0.4, 3.70), 9.2, 1.50, boxstyle='round,pad=0.08',
                        facecolor=C_PARALLEL, edgecolor='#27AE60',
                        linewidth=1.2, alpha=0.6, zorder=1, linestyle='dashed')
ax.add_patch(ev_bg)
ax.text(0.65, 5.12, 'Evidence Layer  —  Parallel MCP & API Calls', ha='left',
        va='center', fontsize=8, fontweight='bold', color='#1E8449', zorder=3)

# Column headers
ax.text(1.9, 5.00, 'RegNetAgents', ha='center', va='center',
        fontsize=8, fontweight='bold', color=C_RNA, zorder=3)
ax.text(5.0, 5.00, 'NCBI E-utilities', ha='center', va='center',
        fontsize=7.5, fontweight='bold', color=C_PUBMED, zorder=3)
ax.text(8.1, 5.00, 'CASCADE', ha='center', va='center',
        fontsize=8, fontweight='bold', color=C_CASCADE, zorder=3)

# RegNetAgents call boxes (left column)
box(0.45, 4.38, 2.9, 0.42, C_RNA,
    'comprehensive_gene_analysis(gene, cell_type)', fontsize=6.2)
box(0.45, 3.84, 2.9, 0.42, C_RNA,
    'pathway_focused_analysis  ·  query_network', fontsize=6.0)

# PubMed call box (centre column — novelty path only)
box(3.65, 4.10, 2.7, 0.68, C_PUBMED,
    'esearch.fcgi  (hit count + recency)\nefetch.fcgi  (publication year)',
    fontsize=6.2)

# CASCADE call boxes (right column)
box(6.65, 4.38, 2.9, 0.42, C_CASCADE,
    'comprehensive_perturbation_analysis(gene, cell_type)', fontsize=6.0)
box(6.65, 3.84, 2.9, 0.42, C_CASCADE,
    'get_gene_metadata  ·  get_protein_interactions  ·  therapeutic_target_discovery',
    fontsize=5.5)

# Routing → evidence layer
arrow(5.0, 5.35, 1.9,  4.80, color=C_ROUTE,  lw=1.5)
arrow(5.0, 5.35, 5.0,  4.78, color=C_PUBMED, lw=1.5, linestyle='dashed')
arrow(5.0, 5.35, 8.1,  4.80, color=C_ROUTE,  lw=1.5)

# Evidence → Synthesis
arrow(1.9, 3.84, 1.9,  2.92, color='#27AE60', lw=1.5)
arrow(5.0, 3.84, 3.4,  2.92, color=C_PUBMED,  lw=1.5, linestyle='dashed')
arrow(8.1, 3.84, 3.55, 2.92, color='#27AE60', lw=1.5)

# ── Explanation Layer ────────────────────────────────────────────────────────
ex_bg = FancyBboxPatch((0.4, 2.07), 9.2, 1.50, boxstyle='round,pad=0.08',
                        facecolor='#E8F8F5', edgecolor='#1ABC9C',
                        linewidth=1.2, alpha=0.7, zorder=1)
ax.add_patch(ex_bg)
ax.text(0.65, 3.45, 'Explanation Layer', ha='left', va='center',
        fontsize=8, fontweight='bold', color='#0E6655', zorder=3)

box(0.7, 2.18, 3.5, 0.70, C_SYNTH,
    'Synthesize\ncross-system corroboration count', fontsize=7.5)
arrow(4.2, 2.53, 4.6, 2.53, color=C_SYNTH, lw=1.8)
box(4.6, 2.18, 4.55, 0.70, C_REPORT,
    'Generate Report\nstructured evidence table  +  optional LLM narrative',
    fontsize=7)

# ── Child servers / external APIs (bottom, dashed) ───────────────────────────
box(0.35, 0.5, 2.85, 1.45, C_RNA,
    'RegNetAgents MCP Server\n\nnetwork analysis  ·  pathway enrichment\n'
    'PageRank centrality  ·  domain agents\n(Bird et al. 2026a)',
    fontsize=6.0, linestyle='dashed', edgecolor=C_RNA, lw=1.5)

box(3.4, 0.5, 3.0, 1.45, C_PUBMED,
    'NCBI E-utilities  (PubMed)\n\nesearch.fcgi  ·  efetch.fcgi\n'
    'NCBI_API_KEY  ·  10 req/s\n(novelty_assessment path)',
    fontsize=6.0, linestyle='dashed', edgecolor=C_PUBMED, lw=1.5)

box(6.6, 0.5, 3.0, 1.45, C_CASCADE,
    'CASCADE MCP Server\n\nperturbation sim  ·  STRING PPI\n'
    'LINCS L1000  ·  DepMap CRISPR  ·  super-enhancers\n'
    'DoRothEA  ·  cBioPortal  (Bird et al. 2026b)',
    fontsize=5.8, linestyle='dashed', edgecolor=C_CASCADE, lw=1.5)

# Arrows from child services UP to evidence layer
arrow(0.5,  1.95, 0.5,  3.70, color=C_RNA,    lw=1.5, linestyle='dashed')
arrow(4.9,  1.95, 4.9,  3.70, color=C_PUBMED, lw=1.5, linestyle='dashed')
arrow(9.1,  1.95, 9.1,  3.70, color=C_CASCADE, lw=1.5)

# ── Legend ───────────────────────────────────────────────────────────────────
legend_items = [
    mpatches.Patch(color=C_ROUTE,   label='Decision'),
    mpatches.Patch(color=C_RNA,     label='RegNetAgents'),
    mpatches.Patch(color=C_PUBMED,  label='NCBI / PubMed'),
    mpatches.Patch(color=C_CASCADE, label='CASCADE'),
    mpatches.Patch(color=C_SYNTH,   label='Synthesis'),
    mpatches.Patch(color=C_REPORT,  label='Report'),
    Line2D([0], [0], color=C_RNA,     lw=1.5, linestyle='dashed',
           label='MCP call (RegNetAgents)'),
    Line2D([0], [0], color=C_PUBMED,  lw=1.5, linestyle='dashed',
           label='API call (NCBI)'),
    Line2D([0], [0], color=C_CASCADE, lw=1.5, linestyle='solid',
           label='MCP call (CASCADE)'),
]
ax.legend(handles=legend_items, loc='lower center', ncol=5, fontsize=6.2,
          frameon=True, fancybox=True, framealpha=0.9,
          bbox_to_anchor=(0.5, -0.02))

plt.tight_layout()
plt.savefig('figure_architecture.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig('figure_architecture.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("Saved figure_architecture.png and figure_architecture.pdf")
