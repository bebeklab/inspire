#!/usr/bin/env python3
"""Create an INSPIRE two-network alignment plot from a structured TSV log.

Expected records are emitted by writeFinalAlignmentBlock() in main.cpp:
    $ALIGNMENT_BEGIN$
    $ALIGNMENT_META$\tpetal_a=...\tpetal_b=...\t...
    $NETWORK_A_COLUMNS$\tedge_id\tnode1\tnode2
    $NETWORK_A_EDGE$\t...
    $NETWORK_B_COLUMNS$\tedge_id\tnode1\tnode2
    $NETWORK_B_EDGE$\t...
    $ALIGNMENT_COLUMNS$\trank\tage_index\tparent_age_index\ttransition_cost\ta1\ta2\tb1\tb2
    $ALIGNMENT_EDGE$\t...
    $ALIGNMENT_END$

The parser independently verifies:
  * one-to-one node mapping in both directions;
  * no repeated undirected edge in either network;
  * number of aligned edge-pair records equals tree_size + 1.

Names are used only to label and annotate the completed alignment plot. They do
not participate in INSPIRE alignment selection or scoring.
"""

import argparse
from collections import Counter
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import networkx as nx


def canonical_edge(first, second):
    """Return a stable representation of an undirected edge."""
    return tuple(sorted((first, second)))


def parse_alignment_log(log_path):
    """Parse and validate one structured INSPIRE alignment log."""
    log_path = Path(log_path)
    metadata = {}
    records = []
    network_a_edges = []
    network_b_edges = []
    inside = False
    found_begin = False
    found_end = False

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")

            if line == "$ALIGNMENT_BEGIN$":
                if found_begin:
                    raise ValueError(
                        f"{log_path}: multiple $ALIGNMENT_BEGIN$ blocks found"
                    )
                found_begin = True
                inside = True
                continue

            if line == "$ALIGNMENT_END$":
                if not inside:
                    raise ValueError(
                        f"{log_path}:{line_number}: unexpected $ALIGNMENT_END$"
                    )
                found_end = True
                inside = False
                break

            if not inside:
                continue

            fields = line.split("\t")
            record_type = fields[0]

            if record_type == "$ALIGNMENT_META$":
                for item in fields[1:]:
                    if "=" not in item:
                        raise ValueError(
                            f"{log_path}:{line_number}: malformed metadata item {item!r}"
                        )
                    key, value = item.split("=", 1)
                    metadata[key] = value

            elif record_type in ("$NETWORK_A_COLUMNS$", "$NETWORK_B_COLUMNS$"):
                expected = [record_type, "edge_id", "node1", "node2"]
                if fields != expected:
                    raise ValueError(
                        f"{log_path}:{line_number}: unexpected network columns: {fields}"
                    )

            elif record_type in ("$NETWORK_A_EDGE$", "$NETWORK_B_EDGE$"):
                if len(fields) != 4:
                    raise ValueError(
                        f"{log_path}:{line_number}: expected 4 fields for full-network edge"
                    )
                try:
                    edge_id = int(fields[1])
                except ValueError as exc:
                    raise ValueError(
                        f"{log_path}:{line_number}: invalid network edge_id"
                    ) from exc
                edge = (fields[2], fields[3])
                target = network_a_edges if record_type == "$NETWORK_A_EDGE$" else network_b_edges
                if edge_id != len(target):
                    raise ValueError(
                        f"{log_path}:{line_number}: non-consecutive network edge_id {edge_id}"
                    )
                target.append(edge)

            elif record_type == "$ALIGNMENT_COLUMNS$":
                expected = [
                    "$ALIGNMENT_COLUMNS$", "rank", "age_index",
                    "parent_age_index", "transition_cost",
                    "a1", "a2", "b1", "b2",
                ]
                if fields != expected:
                    raise ValueError(
                        f"{log_path}:{line_number}: unexpected columns: {fields}"
                    )

            elif record_type == "$ALIGNMENT_EDGE$":
                if len(fields) != 9:
                    raise ValueError(
                        f"{log_path}:{line_number}: expected 9 tab-separated "
                        f"fields, found {len(fields)}"
                    )
                try:
                    record = {
                        "rank": int(fields[1]),
                        "age_index": int(fields[2]),
                        "parent_age_index": int(fields[3]),
                        "transition_cost": int(fields[4]),
                        "a1": fields[5],
                        "a2": fields[6],
                        "b1": fields[7],
                        "b2": fields[8],
                    }
                except ValueError as exc:
                    raise ValueError(
                        f"{log_path}:{line_number}: invalid numeric field"
                    ) from exc
                records.append(record)

    if not found_begin:
        raise ValueError(f"{log_path}: no $ALIGNMENT_BEGIN$ block found")
    if not found_end:
        raise ValueError(f"{log_path}: no $ALIGNMENT_END$ marker found")
    if not metadata:
        raise ValueError(f"{log_path}: no $ALIGNMENT_META$ record found")
    if not network_a_edges:
        raise ValueError(f"{log_path}: no $NETWORK_A_EDGE$ records found")
    if not network_b_edges:
        raise ValueError(f"{log_path}: no $NETWORK_B_EDGE$ records found")
    if not records:
        raise ValueError(f"{log_path}: no $ALIGNMENT_EDGE$ records found")

    required_metadata = {
        "petal_a", "petal_b", "algorithm", "source_age_index",
        "normalized_score", "tree_cost", "tree_size",
    }
    missing = required_metadata - metadata.keys()
    if missing:
        raise ValueError(
            f"{log_path}: missing metadata keys: {', '.join(sorted(missing))}"
        )

    for key in ("source_age_index", "tree_cost", "tree_size"):
        metadata[key] = int(metadata[key])
    metadata["normalized_score"] = float(metadata["normalized_score"])

    validate_alignment(metadata, records, log_path)
    validate_full_networks(metadata, records, network_a_edges, network_b_edges, log_path)
    return metadata, records, network_a_edges, network_b_edges


def validate_alignment(metadata, records, source="alignment"):
    """Raise ValueError if the structured alignment violates core rules."""
    expected_records = metadata["tree_size"] + 1
    if len(records) != expected_records:
        raise ValueError(
            f"{source}: found {len(records)} aligned edge-pair records, "
            f"expected tree_size + 1 = {expected_records}"
        )

    ranks = [record["rank"] for record in records]
    if ranks != list(range(len(records))):
        raise ValueError(
            f"{source}: alignment ranks must be consecutive from 0; found {ranks}"
        )

    age_indices = [record["age_index"] for record in records]
    duplicates = [index for index, count in Counter(age_indices).items() if count > 1]
    if duplicates:
        raise ValueError(f"{source}: repeated AGE indices: {duplicates}")

    if records[0]["age_index"] != metadata["source_age_index"]:
        raise ValueError(
            f"{source}: first AGE index does not match source_age_index"
        )
    if records[0]["parent_age_index"] != -1:
        raise ValueError(f"{source}: source parent_age_index must be -1")

    a_to_b = {}
    b_to_a = {}
    seen_a_edges = set()
    seen_b_edges = set()

    def add_mapping(a_node, b_node, rank):
        previous_b = a_to_b.get(a_node)
        if previous_b is not None and previous_b != b_node:
            raise ValueError(
                f"{source}: one-to-many mapping at rank {rank}: "
                f"{a_node} -> {previous_b} and {b_node}"
            )
        previous_a = b_to_a.get(b_node)
        if previous_a is not None and previous_a != a_node:
            raise ValueError(
                f"{source}: many-to-one mapping at rank {rank}: "
                f"{previous_a} and {a_node} -> {b_node}"
            )
        a_to_b[a_node] = b_node
        b_to_a[b_node] = a_node

    for record in records:
        edge_a = canonical_edge(record["a1"], record["a2"])
        edge_b = canonical_edge(record["b1"], record["b2"])
        if edge_a in seen_a_edges:
            raise ValueError(
                f"{source}: repeated undirected Network A edge at rank "
                f"{record['rank']}: {edge_a}"
            )
        if edge_b in seen_b_edges:
            raise ValueError(
                f"{source}: repeated undirected Network B edge at rank "
                f"{record['rank']}: {edge_b}"
            )
        seen_a_edges.add(edge_a)
        seen_b_edges.add(edge_b)
        add_mapping(record["a1"], record["b1"], record["rank"])
        add_mapping(record["a2"], record["b2"], record["rank"])

    return a_to_b, b_to_a


def validate_full_networks(metadata, records, network_a_edges, network_b_edges, source):
    """Confirm full networks are unique and contain every aligned edge."""
    full_a = [canonical_edge(*edge) for edge in network_a_edges]
    full_b = [canonical_edge(*edge) for edge in network_b_edges]
    if len(full_a) != len(set(full_a)):
        raise ValueError(f"{source}: duplicate undirected edge in full Network A")
    if len(full_b) != len(set(full_b)):
        raise ValueError(f"{source}: duplicate undirected edge in full Network B")

    aligned_a = {canonical_edge(r["a1"], r["a2"]) for r in records}
    aligned_b = {canonical_edge(r["b1"], r["b2"]) for r in records}
    missing_a = aligned_a - set(full_a)
    missing_b = aligned_b - set(full_b)
    if missing_a:
        raise ValueError(f"{source}: aligned A edges absent from full network: {sorted(missing_a)}")
    if missing_b:
        raise ValueError(f"{source}: aligned B edges absent from full network: {sorted(missing_b)}")


def build_alignment_data(metadata, records, network_a_edges, network_b_edges):
    """Convert records into graph, edge-class, and mapping structures."""
    graph_a = nx.Graph()
    graph_b = nx.Graph()
    graph_a.add_edges_from(network_a_edges)
    graph_b.add_edges_from(network_b_edges)
    exact_edges_a = []
    exact_edges_b = []
    substituted_edges_a = []
    substituted_edges_b = []
    a_to_b = {}

    for record in records:
        edge_a = (record["a1"], record["a2"])
        edge_b = (record["b1"], record["b2"])
        # Post hoc display classification only. Names never affect alignment.
        exact_pair = (
            record["a1"] == record["b1"]
            and record["a2"] == record["b2"]
        )
        if exact_pair:
            exact_edges_a.append(edge_a)
            exact_edges_b.append(edge_b)
        else:
            substituted_edges_a.append(edge_a)
            substituted_edges_b.append(edge_b)

        a_to_b[record["a1"]] = record["b1"]
        a_to_b[record["a2"]] = record["b2"]

    aligned_a_set = {canonical_edge(r["a1"], r["a2"]) for r in records}
    aligned_b_set = {canonical_edge(r["b1"], r["b2"]) for r in records}
    unaligned_edges_a = [e for e in network_a_edges if canonical_edge(*e) not in aligned_a_set]
    unaligned_edges_b = [e for e in network_b_edges if canonical_edge(*e) not in aligned_b_set]

    return {
        "graph_a": graph_a,
        "graph_b": graph_b,
        "exact_edges_a": exact_edges_a,
        "exact_edges_b": exact_edges_b,
        "substituted_edges_a": substituted_edges_a,
        "substituted_edges_b": substituted_edges_b,
        "unaligned_edges_a": unaligned_edges_a,
        "unaligned_edges_b": unaligned_edges_b,
        "a_to_b": a_to_b,
    }


def alignment_aware_layout(graph_a, graph_b, a_to_b, seed=7, offset=7.5):
    """Lay out A, then mirror mapped B nodes at matching vertical positions."""
    pos_a_raw = nx.spring_layout(graph_a, seed=seed, k=1.5)
    pos_b_raw = nx.spring_layout(graph_b, seed=seed, k=1.5)

    pos_a = {
        node: (float(x) - offset / 2.0, float(y))
        for node, (x, y) in pos_a_raw.items()
    }

    reverse_mapping = {b: a for a, b in a_to_b.items()}
    pos_b = {}
    for node, (x, y) in pos_b_raw.items():
        if node in reverse_mapping and reverse_mapping[node] in pos_a_raw:
            a_node = reverse_mapping[node]
            a_x, a_y = pos_a_raw[a_node]
            pos_b[node] = (float(a_x) + offset / 2.0, float(a_y))
        else:
            pos_b[node] = (float(x) + offset / 2.0, float(y))

    return pos_a, pos_b


def draw_alignment(metadata, records, network_a_edges, network_b_edges,
                   output_path, dpi=180, seed=7, show=False):
    NETWORK_GAP = 3.6
    data = build_alignment_data(metadata, records, network_a_edges, network_b_edges)
    graph_a = data["graph_a"]
    graph_b = data["graph_b"]
    a_to_b = data["a_to_b"]
    pos_a, pos_b = alignment_aware_layout(
        graph_a, graph_b, a_to_b, seed=seed, offset= NETWORK_GAP
    )

    substituted_a_nodes = {a for a, b in a_to_b.items() if a != b}
    substituted_b_nodes = {b for a, b in a_to_b.items() if a != b}

    color_a = [
        "#FF8C42" if node in substituted_a_nodes else "#4C9BE8"
        for node in graph_a.nodes()
    ]
    color_b = [
        "#FF8C42" if node in substituted_b_nodes else "#E8694C"
        for node in graph_b.nodes()
    ]

    # fig, ax = plt.subplots(figsize=(24, 12))
    fig, ax = plt.subplots(figsize=(16, 10))
    background = "#FFFFFF"
    ax.set_facecolor(background)
    fig.patch.set_facecolor(background)

    # Full-network context first: unaligned interactions are thin and subdued.
    nx.draw_networkx_edges(
        graph_a,
        pos_a,
        ax=ax,
        edgelist=data["unaligned_edges_a"],
        edge_color="#93A4B8",
        width=1.0,
        alpha=0.65
    )
    
    nx.draw_networkx_edges(
        graph_b,
        pos_b,
        ax=ax,
        edgelist=data["unaligned_edges_b"],
        edge_color="#C49A9A",
        width=1.0,
        alpha=0.65
    )
    # Aligned interactions are emphasized.
    nx.draw_networkx_edges(
        graph_a, pos_a, ax=ax, edgelist=data["exact_edges_a"],
        edge_color="#5BA4F5", width=2.5, alpha=0.90
    )
    nx.draw_networkx_edges(
        graph_a, pos_a, ax=ax, edgelist=data["substituted_edges_a"],
        edge_color="#FF8C42", width=2.5, alpha=0.90
    )
    nx.draw_networkx_edges(
        graph_b, pos_b, ax=ax, edgelist=data["exact_edges_b"],
        edge_color="#F5615A", width=2.5, alpha=0.90
    )
    nx.draw_networkx_edges(
        graph_b, pos_b, ax=ax, edgelist=data["substituted_edges_b"],
        edge_color="#FF8C42", width=2.5, alpha=0.90
    )

    nx.draw_networkx_nodes(
        graph_a, pos_a, ax=ax, node_color=color_a, node_size=1100,
        edgecolors="#374151", linewidths=1.5
    )
    nx.draw_networkx_nodes(
        graph_b, pos_b, ax=ax, node_color=color_b, node_size=1100,
        edgecolors="#374151", linewidths=1.5
    )
    nx.draw_networkx_labels(
        graph_a, pos_a, ax=ax, font_size=7,
        font_weight="bold", font_color="#111827"
    )
    nx.draw_networkx_labels(
        graph_b, pos_b, ax=ax, font_size=7,
        font_weight="bold", font_color="#111827"
    )

    for a_node, b_node in sorted(a_to_b.items()):
        if a_node not in pos_a or b_node not in pos_b:
            continue
        same_label = a_node == b_node
        ax.annotate(
            "",
            xy=pos_b[b_node],
            xytext=pos_a[a_node],
            arrowprops={
                "arrowstyle": "-",
                "color": "#6B7280" if same_label else "#C2410C",
                "lw": 1.5 if same_label else 2.0,
                "linestyle": "dashed",
                "alpha": 0.75 if same_label else 0.90,
                "connectionstyle": "arc3,rad=0.08",
            },
        )

    offset = NETWORK_GAP
    ax.text(
        -offset / 2.0, 1.35, f"Network A: {metadata['petal_a']}",
        ha="center", fontsize=16, fontweight="bold", color="#5BA4F5"
    )
    ax.text(
        offset / 2.0, 1.35, f"Network B: {metadata['petal_b']}",
        ha="center", fontsize=16, fontweight="bold", color="#F5615A"
    )

    legend_elements = [
        Line2D([0], [0], color="#355A7A", lw=1.0, alpha=0.7,
               label="Network A: unaligned interaction"),
        Line2D([0], [0], color="#7A3F45", lw=1.0, alpha=0.7,
               label="Network B: unaligned interaction"),
        Line2D([0], [0], color="#5BA4F5", lw=2.5,
               label="Network A: exact displayed edge pair"),
        Line2D([0], [0], color="#F5615A", lw=2.5,
               label="Network B: exact displayed edge pair"),
        Line2D([0], [0], color="#FF8C42", lw=2.5,
               label="Non-identical aligned edge pair"),
        Line2D([0], [0], color="#6B7280", lw=1.5, linestyle="--",
               label="Mapping with same displayed label"),
        Line2D([0], [0], color="#C2410C", lw=2.0, linestyle="--",
               label="Mapping with different displayed labels"),
        mpatches.Patch(facecolor="#FF8C42", edgecolor="#111827",
                       label="Node mapped to a different displayed label"),
    ]
    legend = ax.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=4,
        fontsize=9,
        framealpha=1.0,
        facecolor="#FFFFFF",
        edgecolor="#9CA3AF",
        labelcolor="#111827",
        bbox_to_anchor=(0.5, -0.07)
    )
    for text in legend.get_texts():
        text.set_color("#111827")
    ax.set_title(
        f"Network Alignment · {metadata['petal_a']} vs {metadata['petal_b']} · "
        f"Score {metadata['normalized_score']:.6f} · "
        f"Aligned edge pairs {len(records)} · "
        f"Full edges {len(network_a_edges)}/{len(network_b_edges)}",
        fontsize=15, fontweight="bold", color="#111827", pad=16
    )
    ax.axis("off")
    plt.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path, dpi=dpi, bbox_inches="tight",
        facecolor=fig.get_facecolor()
    )
    if show:
        plt.show()
    plt.close(fig)


def default_output_path(log_path, output_dir=None):
    log_path = Path(log_path)
    filename = f"{log_path.stem}__alignment.png"
    if output_dir is None:
        return log_path.with_name(filename)
    return Path(output_dir) / filename


def main():
    parser = argparse.ArgumentParser(
        description="Plot a validated INSPIRE alignment from its structured TSV log."
    )
    parser.add_argument("log", nargs="+", help="Structured alignment TSV file(s)")
    parser.add_argument(
        "-o", "--output",
        help="Output PNG path. Allowed only when one log is supplied."
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for generated PNG files."
    )
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--show", action="store_true")
    parser.add_argument(
        "--validate-only", action="store_true",
        help="Validate logs without generating plots."
    )
    args = parser.parse_args()

    if args.output and len(args.log) != 1:
        parser.error("--output can be used only with one input log")
    if args.output and args.output_dir:
        parser.error("use either --output or --output-dir, not both")

    failed = False
    for log_name in args.log:
        try:
            metadata, records, network_a_edges, network_b_edges = parse_alignment_log(log_name)
            print(
                f"VALID {log_name}: {metadata['petal_a']} vs "
                f"{metadata['petal_b']}, score={metadata['normalized_score']:.6f}, "
                f"aligned_edge_pairs={len(records)}, "
                f"full_edges={len(network_a_edges)}/{len(network_b_edges)}"
            )
            if args.validate_only:
                continue

            output_path = (
                Path(args.output)
                if args.output
                else default_output_path(log_name, args.output_dir)
            )
            draw_alignment(
                metadata, records, network_a_edges, network_b_edges, output_path,
                dpi=args.dpi, seed=args.seed, show=args.show
            )
            print(f"WROTE {output_path}")
        except Exception as exc:
            failed = True
            print(f"ERROR {log_name}: {exc}", file=sys.stderr)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
