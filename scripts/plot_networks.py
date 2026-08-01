import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

# ─────────────────────────────────────────────────────────────────
# 1. ALL edges extracted from the alignment output
#    Every [A1::A2] and [B1::B2] pair is an edge
# ─────────────────────────────────────────────────────────────────

edges_A = [
    # conserved (cost 100) edges
    ("ATM",   "NF1"),
    ("EGFR",  "ERBIN"),
    ("ERBIN", "APC"),
    ("APC",   "CDC42"),
    ("APC",   "GSK3B"),
    ("APC",   "SMAD4"),
    ("APC",   "TCF7L2"),
    ("APC",   "TNKS"),
    ("CDC42", "PTEN"),
    ("GSK3B", "MDM2"),
    ("GSK3B", "PTEN"),
    ("PTEN",  "SMAD4"),
    ("PTEN",  "TCF7L2"),
    ("PTEN",  "TNKS"),
    # substituted (cost 217) edges
    ("BRCA1", "GSK3B"),
    ("MDM2",  "NF1"),
    ("NF1",   "EGFR"),
    ("NF1",   "PTEN"),
    ("SMAD4", "ATM"),
]

edges_B = [
    # conserved (cost 100) edges
    ("FBXW7", "KRAS"),
    ("EGFR",  "ERBIN"),
    ("ERBIN", "APC"),
    ("APC",   "CDC42"),
    ("APC",   "GSK3B"),
    ("APC",   "SMAD4"),
    ("APC",   "TCF7L2"),
    ("APC",   "TNKS"),
    ("CDC42", "PTEN"),
    ("GSK3B", "MDM2"),
    ("GSK3B", "PTEN"),
    ("PTEN",  "SMAD4"),
    ("PTEN",  "TCF7L2"),
    ("PTEN",  "TNKS"),
    # substituted (cost 217) edges
    ("FBXW7", "GSK3B"),
    ("MDM2",  "KRAS"),
    ("KRAS",  "EGFR"),
    ("KRAS",  "PTEN"),
    ("SMAD4", "BCL2"),
]

# ─────────────────────────────────────────────────────────────────
# 2. Node alignment map:  A_node → B_node
#    Derived by reading A1↔B1 and A2↔B2 from every aligned pair
# ─────────────────────────────────────────────────────────────────

alignment_map = {
    "ATM":    "FBXW7",   # ATM::NF1  ↔  FBXW7::KRAS
    "NF1":    "KRAS",
    "EGFR":   "EGFR",
    "ERBIN":  "ERBIN",
    "APC":    "APC",
    "CDC42":  "CDC42",
    "GSK3B":  "GSK3B",
    "MDM2":   "MDM2",
    "PTEN":   "PTEN",
    "SMAD4":  "SMAD4",
    "TCF7L2": "TCF7L2",
    "TNKS":   "TNKS",
    "BRCA1":  "FBXW7",   # BRCA1::GSK3B  ↔  FBXW7::GSK3B
    # SMAD4::ATM ↔ SMAD4::BCL2  → ATM maps to BCL2
    # (SMAD4 already mapped to SMAD4 above)
}
# ATM appears as A2 in SMAD4::ATM ↔ SMAD4::BCL2, so ATM→BCL2 as well
# ATM is already mapped to FBXW7 as A1; BCL2 is a separate B2 target
# We keep both as separate alignment edges (one node can align to multiple)
extra_alignment = [("ATM", "BCL2")]   # from SMAD4::ATM ↔ SMAD4::BCL2

# Build full alignment edge list
alignment_edges = list(alignment_map.items()) + extra_alignment

# ─────────────────────────────────────────────────────────────────
# 3. Classify edges: conserved vs substituted
# ─────────────────────────────────────────────────────────────────

# An aligned edge pair is CONSERVED (cost=100) when A-nodes == B-nodes
# It is SUBSTITUTED (cost=217) when at least one node differs

conserved_edges_A = []
substituted_edges_A = []

conserved_edges_B = []
substituted_edges_B = []

aligned_edge_pairs = [
    (("ATM",   "NF1"),    ("FBXW7", "KRAS")),
    (("EGFR",  "ERBIN"),  ("EGFR",  "ERBIN")),
    (("ERBIN", "APC"),    ("ERBIN", "APC")),
    (("APC",   "CDC42"),  ("APC",   "CDC42")),
    (("APC",   "GSK3B"),  ("APC",   "GSK3B")),
    (("APC",   "SMAD4"),  ("APC",   "SMAD4")),
    (("APC",   "TCF7L2"), ("APC",   "TCF7L2")),
    (("APC",   "TNKS"),   ("APC",   "TNKS")),
    (("CDC42", "PTEN"),   ("CDC42", "PTEN")),
    (("GSK3B", "MDM2"),   ("GSK3B", "MDM2")),
    (("GSK3B", "PTEN"),   ("GSK3B", "PTEN")),
    (("PTEN",  "SMAD4"),  ("PTEN",  "SMAD4")),
    (("PTEN",  "TCF7L2"), ("PTEN",  "TCF7L2")),
    (("PTEN",  "TNKS"),   ("PTEN",  "TNKS")),
    (("BRCA1", "GSK3B"),  ("FBXW7", "GSK3B")),
    (("MDM2",  "NF1"),    ("MDM2",  "KRAS")),
    (("NF1",   "EGFR"),   ("KRAS",  "EGFR")),
    (("NF1",   "PTEN"),   ("KRAS",  "PTEN")),
    (("SMAD4", "ATM"),    ("SMAD4", "BCL2")),
]

for (eA, eB) in aligned_edge_pairs:
    if set(eA) == set(eB):
        conserved_edges_A.append(eA)
        conserved_edges_B.append(eB)
    else:
        substituted_edges_A.append(eA)
        substituted_edges_B.append(eB)

# ─────────────────────────────────────────────────────────────────
# 4. Build graphs
# ─────────────────────────────────────────────────────────────────

GA = nx.Graph()
GA.add_edges_from(edges_A)

GB = nx.Graph()
GB.add_edges_from(edges_B)

# ─────────────────────────────────────────────────────────────────
# 5. Shared layout: nodes with same name get mirrored Y positions
#    Network A on the left, Network B on the right
# ─────────────────────────────────────────────────────────────────

OFFSET = 7.0

# Use a combined supergraph to get a coherent spring layout
# then split positions
shared_nodes = set(GA.nodes()) & set(GB.nodes())

# Compute layout on GA as the "reference"
pos_ref = nx.spring_layout(GA, seed=7, k=2.2)

# For nodes only in GB, compute separately
only_B = set(GB.nodes()) - set(GA.nodes())
GB_sub = GB.subgraph(only_B)
pos_only_B = nx.spring_layout(GB_sub, seed=7, k=2.2) if only_B else {}

# Build final positions
pos_A = {n: (xy[0] - OFFSET / 2, xy[1]) for n, xy in pos_ref.items()}

pos_B = {}
for n in GB.nodes():
    if n in pos_ref:                    # shared name → mirror position
        pos_B[n] = (pos_ref[n][0] + OFFSET / 2, pos_ref[n][1])
    else:                               # B-only node
        raw = pos_only_B.get(n, (0, 0))
        pos_B[n] = (raw[0] + OFFSET / 2, raw[1])

# ─────────────────────────────────────────────────────────────────
# 6. Node colours
# ─────────────────────────────────────────────────────────────────

substituted_A_nodes = {"ATM", "NF1", "BRCA1"}          # differ from B
substituted_B_nodes = {"FBXW7", "KRAS", "BCL2"}

def node_colors(graph, sub_nodes, base_color, sub_color, conserved_color):
    colors = []
    for n in graph.nodes():
        if n in sub_nodes:
            colors.append(sub_color)
        else:
            colors.append(conserved_color)
    return colors

nc_A = node_colors(GA, substituted_A_nodes, "#4C9BE8", "#FF8C42", "#4C9BE8")
nc_B = node_colors(GB, substituted_B_nodes, "#E8694C", "#FF8C42", "#E8694C")

# ─────────────────────────────────────────────────────────────────
# 7. Draw
# ─────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(24, 12))
ax.set_facecolor("#1C1C2E")
fig.patch.set_facecolor("#1C1C2E")

# ── Network A edges ──────────────────────────────────────────────
nx.draw_networkx_edges(GA, pos_A, ax=ax,
                       edgelist=conserved_edges_A,
                       edge_color="#5BA4F5", width=2.5, alpha=0.9)
nx.draw_networkx_edges(GA, pos_A, ax=ax,
                       edgelist=substituted_edges_A,
                       edge_color="#FF8C42", width=2.5,
                       alpha=0.9, style="solid")

# ── Network B edges ──────────────────────────────────────────────
nx.draw_networkx_edges(GB, pos_B, ax=ax,
                       edgelist=conserved_edges_B,
                       edge_color="#F5615A", width=2.5, alpha=0.9)
nx.draw_networkx_edges(GB, pos_B, ax=ax,
                       edgelist=substituted_edges_B,
                       edge_color="#FF8C42", width=2.5,
                       alpha=0.9, style="solid")

# ── Nodes ────────────────────────────────────────────────────────
nx.draw_networkx_nodes(GA, pos_A, ax=ax,
                       node_color=nc_A, node_size=1100,
                       edgecolors="white", linewidths=1.5)
nx.draw_networkx_nodes(GB, pos_B, ax=ax,
                       node_color=nc_B, node_size=1100,
                       edgecolors="white", linewidths=1.5)

# ── Labels ───────────────────────────────────────────────────────
nx.draw_networkx_labels(GA, pos_A, ax=ax,
                        font_size=7, font_weight="bold", font_color="white")
nx.draw_networkx_labels(GB, pos_B, ax=ax,
                        font_size=7, font_weight="bold", font_color="white")

# ── Alignment dashed edges ───────────────────────────────────────
for (a_node, b_node) in alignment_edges:
    if a_node in pos_A and b_node in pos_B:
        xA, yA = pos_A[a_node]
        xB, yB = pos_B[b_node]
        same = (a_node == b_node)
        ax.annotate("",
                    xy=(xB, yB), xytext=(xA, yA),
                    arrowprops=dict(
                        arrowstyle="-",
                        color="white" if same else "#FFD700",
                        lw=1.5 if same else 2.0,
                        linestyle="dashed",
                        alpha=0.55 if same else 0.85,
                        connectionstyle="arc3,rad=0.08"
                    ))

# ─────────────────────────────────────────────────────────────────
# 8. Network labels, legend, title
# ─────────────────────────────────────────────────────────────────

ax.text(-OFFSET / 2, 1.75, "Network A",
        ha="center", fontsize=16, fontweight="bold", color="#5BA4F5")
ax.text(+OFFSET / 2, 1.75, "Network B",
        ha="center", fontsize=16, fontweight="bold", color="#F5615A")

legend_elements = [
    Line2D([0], [0], color="#5BA4F5", lw=2.5,
           label="Network A — conserved edge (cost 100)"),
    Line2D([0], [0], color="#F5615A", lw=2.5,
           label="Network B — conserved edge (cost 100)"),
    Line2D([0], [0], color="#FF8C42", lw=2.5,
           label="Substituted edge (cost 217)"),
    Line2D([0], [0], color="white",  lw=1.5, linestyle="--",
           label="Alignment — conserved node (same name)"),
    Line2D([0], [0], color="#FFD700", lw=2.0, linestyle="--",
           label="Alignment — substituted node (different name)"),
    mpatches.Patch(facecolor="#FF8C42", edgecolor="white",
                   label="Substituted node"),
]
ax.legend(handles=legend_elements, loc="lower center",
          ncol=3, fontsize=9, framealpha=0.25,
          facecolor="#2C2C3E", labelcolor="white",
          bbox_to_anchor=(0.5, -0.06))

ax.set_title(
    "Network Alignment  ·  Run 554  ·  Score 0.620096",
    fontsize=15, fontweight="bold", color="white", pad=16
)
ax.axis("off")
plt.tight_layout()
plt.savefig("network_alignment_v2.png", dpi=180, bbox_inches="tight",
            facecolor=fig.get_facecolor())
plt.show()
