"""Phase 87 — NetworkX Bipartite Graph: Compounds ↔ KRAS target."""
import sys
import os

if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite


FAMILY_COLORS = {
    "benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
    "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860",
    "other": "#808080",
}


def get_family(compound_name: str) -> str:
    """Extract scaffold family from compound name."""
    prefix = compound_name.split("_")[0]
    return prefix if prefix in FAMILY_COLORS else "other"


def build_bipartite_graph(df: pd.DataFrame) -> nx.Graph:
    """Build bipartite graph: compounds (set 0) ↔ KRAS (set 1)."""
    G = nx.Graph()
    # Add target node
    G.add_node("KRAS", bipartite=1)
    # Add compound nodes and edges
    for _, row in df.iterrows():
        name = row["compound_name"]
        G.add_node(name, bipartite=0, family=get_family(name), pic50=row["pic50"])
        G.add_edge(name, "KRAS", weight=row["pic50"])
    return G


def visualize(G: nx.Graph, df: pd.DataFrame, out_path: str) -> None:
    """Visualize bipartite graph with scaffold family coloring."""
    compounds = [n for n, d in G.nodes(data=True) if d.get("bipartite") == 0]
    targets = [n for n, d in G.nodes(data=True) if d.get("bipartite") == 1]

    # Bipartite layout
    pos = {}
    for i, c in enumerate(sorted(compounds)):
        pos[c] = (0, -i)
    pos["KRAS"] = (2, -len(compounds) / 2)

    fig, ax = plt.subplots(figsize=(14, max(10, len(compounds) * 0.3)))

    # Draw edges with alpha proportional to pIC50
    pic50_vals = [G.edges[e]["weight"] for e in G.edges()]
    min_p, max_p = min(pic50_vals), max(pic50_vals)
    for u, v, d in G.edges(data=True):
        alpha = 0.3 + 0.7 * (d["weight"] - min_p) / (max_p - min_p) if max_p > min_p else 0.5
        ax.plot(
            [pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
            color="#999999", alpha=alpha, linewidth=0.8,
        )

    # Draw compound nodes colored by family
    for c in compounds:
        family = G.nodes[c].get("family", "other")
        color = FAMILY_COLORS.get(family, "#808080")
        ax.scatter(pos[c][0], pos[c][1], c=color, s=60, zorder=3)
        ax.text(pos[c][0] - 0.05, pos[c][1], c, fontsize=6, ha="right", va="center")

    # Draw KRAS node
    ax.scatter(pos["KRAS"][0], pos["KRAS"][1], c="red", s=200, zorder=3, marker="s")
    ax.text(pos["KRAS"][0] + 0.1, pos["KRAS"][1], "KRAS", fontsize=12, fontweight="bold", va="center")

    # Legend
    for family, color in FAMILY_COLORS.items():
        ax.scatter([], [], c=color, label=family, s=60)
    ax.legend(title="Scaffold Family", loc="upper right", fontsize=8)

    ax.set_title(f"Bipartite Graph: {len(compounds)} Compounds ↔ KRAS", fontsize=14)
    ax.set_xlim(-1.5, 3)
    ax.axis("off")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Graph saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build bipartite graph of compounds ↔ KRAS target.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--compounds", default="data/compounds.csv", help="Path to compounds CSV")
    args = parser.parse_args()

    os.makedirs("output", exist_ok=True)

    df = pd.read_csv(args.compounds)
    print(f"Loaded {len(df)} compounds from {args.compounds}")

    G = build_bipartite_graph(df)
    compounds = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
    targets = {n for n, d in G.nodes(data=True) if d["bipartite"] == 1}

    print(f"Graph: {len(compounds)} compounds, {len(targets)} targets, {G.number_of_edges()} edges")
    print(f"Is bipartite: {bipartite.is_bipartite(G)}")

    # Family breakdown
    families = df["compound_name"].apply(get_family).value_counts()
    print(f"\nScaffold families:")
    for fam, count in families.items():
        print(f"  {fam}: {count}")

    visualize(G, df, "output/bipartite_graph.png")


if __name__ == "__main__":
    main()
