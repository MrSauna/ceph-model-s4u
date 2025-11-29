import json
import argparse
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt

#!/usr/bin/env python3
"""
crushmap2dot.py

Read a Ceph CRUSH map JSON (e.g. `ceph osd crush dump -f json`) and build a NetworkX DiGraph
representing the bucket/device hierarchy.

Usage:
    python crushmap2dot.py -i crushmap.json -o crushmap.graphml
"""



def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def index_items(crush):
    """
    Build an index mapping CRUSH item id -> item dict with at least 'id', 'name', 'type' (if present).
    Supports common keys: 'buckets', 'devices', 'nodes', 'osds'.
    """
    idx = {}

    # buckets
    for key in ("buckets",):
        for b in crush.get(key, []):
            idx[b["id"]] = dict(b)  # copy so we can annotate
            idx[b["id"]].setdefault("type", b.get("type", "bucket"))
            idx[b["id"]].setdefault("name", b.get("name", f"bucket-{b['id']}"))

    # devices / osds / nodes
    for key in ("devices", "osds", "nodes"):
        for d in crush.get(key, []):
            idx[d["id"]] = dict(d)
            # normalize name/type fields
            idx[d["id"]].setdefault("name", d.get("name", d.get("device", f"dev-{d['id']}")))
            if "type" not in idx[d["id"]]:
                # typical osd entries might have 'type' or not; set sensible default
                idx[d["id"]].setdefault("type", d.get("type", "device"))

    return idx


def build_graph(crush):
    """
    Return a NetworkX DiGraph where edges point from bucket -> child (bucket or device).
    Node attributes include: id, name, type, weight (if available), and raw (original dict).
    """
    G = nx.DiGraph()
    idx = index_items(crush)

    # connect bucket items to their children
    for bucket in crush.get("buckets", []):
        if "~" in bucket["name"]:
            continue
        if bucket["name"] == "default":
            continue

        parent_id = bucket["id"]
        G.add_node(parent_id, weight=bucket.get("weight"), type=bucket.get("type_name"))

        for item in bucket.get("items", []):
            child_id = item.get("id")

            # osds
            u, v = parent_id, child_id
            if child_id >= 0: #devices
                G.add_node(child_id, weight=item.get("weight"), type="device")
                G.add_edge(u, v, width=0.1)

            else:
                G.add_edge(u, v, width=1)


    return G


def main():
    p = argparse.ArgumentParser(description="Build NetworkX graph from CRUSH map JSON.")
    p.add_argument("-i", "--input", default="crushmap.json", help="CRUSH JSON file (default: crushmap.json)")
    p.add_argument("-o", "--output", help="Optional output graphml filename to write the graph")
    p.add_argument("--show-summary", action="store_true", help="Print a short summary of the graph")
    args = p.parse_args()

    path = Path(args.input)
    if not path.exists():
        raise SystemExit(f"Input file not found: {path}")

    crush = load_json(path)
    G = build_graph(crush)
    # print the graph to console in plain text (for debugging)
    nx.write_network_text(G)

    # figure out hierarchy labels
    roots = [n for n, d in G.in_degree() if d == 0]
    levels = []
    for level_nodes in nx.bfs_layers(G, roots):
        level = G.nodes[level_nodes[0]].get("type")
        levels.append(level)

    pos = nx.nx_agraph.graphviz_layout(G,
                                       prog="dot",
                                       args="-Gnodesep=0.5 -Granksep=1.0 -Gdpi=300")
    # pos = nx.nx_agraph.graphviz_layout(G, prog="twopi")

    # edge widths
    widths = [G[u][v]['width'] for u,v in G.edges()]

    # node sizes
    node_sizes = [node.get("weight", 1) / 1e6 *2 for _, node in G.nodes(data=True)]

    fig, ax = plt.subplots(figsize=(16, 8))
    nx.draw(G, pos,
            ax=ax,
            node_size=node_sizes,
            width=widths,
            arrows=False,)
    ax.set_axis_on()

    uniq_y = sorted(list(set(y for x, y in pos.values())), reverse=True)
    ax.set_yticks(uniq_y, labels=levels)

    # try to hack the yticks to show
    ax.tick_params(axis='y', left=True, labelleft=True, direction='in')
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # nx.draw_networkx_labels(G, pos)
    # plt.show()    
    plt.savefig("crushmap-graph.svg", format="svg", bbox_inches="tight")

    if args.show_summary:
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        types = {}
        for _, d in G.nodes(data=True):
            t = d.get("type", "unknown")
            types[t] = types.get(t, 0) + 1
        print(f"Graph: {num_nodes} nodes, {num_edges} edges")
        print("Node types:", types)

    if args.output:
        out = Path(args.output)
        # choose GraphML for richer attribute support; NetworkX writes .graphml
        nx.write_graphml(G, out)
        print(f"Wrote graph to {out}")

    # default behaviour: if neither output nor summary requested, print summary
    if not args.output and not args.show_summary:
        print(f"Loaded CRUSH map with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

if __name__ == "__main__":
    main()