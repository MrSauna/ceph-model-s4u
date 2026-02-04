import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import json

try:
    from tools.sim_analysis import SimulationRun
except ImportError:
    # Handle case where tools is not in pythonpath (e.g. running from root)
    import sys
    sys.path.append(os.getcwd())
    from tools.sim_analysis import SimulationRun


if "snakemake" not in globals():
    mon_metrics_file = "results/vscode/mon_metrics.csv"
    net_metrics_file = "results/vscode/net_metrics.csv"
    client_metrics_file = "results/vscode/client_metrics.csv"
    output_dir = "results/vscode/figures"
    git_hash = "unknown"
else:
    mon_metrics_file = snakemake.input.mon_metrics
    net_metrics_file = snakemake.input.net_metrics
    client_metrics_file = snakemake.input.client_metrics
    output_dir = snakemake.params.output_dir
    git_hash = snakemake.params.get("git_hash", "unknown")


def add_git_hash(fig):
    fig.text(0.99, 0.01, f'Commit: {git_hash}', 
             ha='right', va='bottom', 
             fontsize=8, color='gray', alpha=0.5)


def viz_mon_metrics(sim_run):
    mon_metrics_df = sim_run.mon_df
    if mon_metrics_df is None:
        return

    fig = plt.figure(figsize=(10, 6))
    
    # Plot all columns except 'time'
    for column in mon_metrics_df.columns:
        if column == "time":
            continue
        plt.plot(mon_metrics_df["time"], mon_metrics_df[column], label=column, marker='o', linestyle='-', markersize=2)

    plt.xlabel("Time (s)")
    plt.ylabel("Number of PGs")
    plt.title("PG States Over Time")
    plt.legend()
    plt.grid(True)
    add_git_hash(fig)
    
    output_path = os.path.join(output_dir, "pgmap_plot.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()


def viz_client_metrics(sim_run):
    throughput = sim_run.get_client_throughput()
    
    if throughput.empty:
        print("Client metrics empty, creating empty plot")
        fig = plt.figure(figsize=(10, 6))
        plt.xlabel("Time (s)")
        plt.ylabel("Throughput (MiB/s)")
        plt.title("Client Throughput (No Data)")
        plt.grid(True)
        add_git_hash(fig)
        output_path = os.path.join(output_dir, "client_throughput.svg")
        plt.savefig(output_path)
        plt.close()
        return

    fig = plt.figure(figsize=(10, 6))
    plot_index = throughput.index
    
    window_size = 10 # This is hardcoded in sim_analysis too, maybe should expose it?

    if "read" in throughput.columns:
        plt.plot(plot_index, throughput["read"], label=f"Read ({window_size}s avg)", linestyle='-', linewidth=1.5)
    
    if "write" in throughput.columns:
        plt.plot(plot_index, throughput["write"], label=f"Write ({window_size}s avg)", linestyle='--', linewidth=1.5)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (MiB/s)")
    plt.title(f"Client Throughput Over Time (Smoothed {window_size}s)")
    plt.legend()
    plt.grid(True)
    add_git_hash(fig)
    
    output_path = os.path.join(output_dir, "client_throughput.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()

def viz_net_metrics(sim_run):
    net_df = sim_run.net_df
    if net_df is None or net_df.empty:
        print("Net metrics empty, skipping viz_net_metrics")
        return

    topology = sim_run.topology
    if not topology:
        print(f"Topology not loaded")
        return

    # 1. Parse Topology and calculate Tree Layout
    # Simple recursive layout:
    # x: centered above children
    # y: based on depth
    
    node_pos = {}  # id -> (x, y)
    node_names = {} # id -> name
    leaf_nodes = set() # ids of leaf nodes (hosts)
    edges = []     # (parent_id, child_id, uplink_name, bandwidth)

    def format_bw(bps):
        if bps >= 1e9:
            return f"{bps/1e9:.0f}G"
        elif bps >= 1e6:
            return f"{bps/1e6:.0f}M"
        elif bps >= 1e3:
            return f"{bps/1e3:.0f}K"
        return f"{bps:.0f}"

    def parse_and_layout(node, depth, x_start):
        """
        Returns (width, x_center) of the subtree
        """
        raw_children = node.get("children", [])
        node_id = node["id"]
        
        # Separate children: actors vs structural
        structural_children = []
        actor_children = []
        for child in raw_children:
            if child.get("type") == "actor":
                actor_children.append(child)
            else:
                structural_children.append(child)
        
        # Construct Label
        name_str = node.get("name", node_id)
        if actor_children:
            # Sort for consistency
            actor_names = []
            for c in actor_children:
                raw_name = c.get("name", c["id"])
                # Remove ", disk..." if present
                if "," in raw_name:
                    clean_name = raw_name.split(",")[0].strip()
                else:
                    clean_name = raw_name
                
                if clean_name.startswith("osd."):
                    clean_name = clean_name[4:]
                    
                actor_names.append(clean_name)

            formatted_actors = [f"({n})" for n in actor_names]
            name_str += "\\n" + "\\n".join(formatted_actors)
            
        if node_id != "_world_":
            node_names[node_id] = name_str
        
        if not structural_children:
            # Leaf (topology-wise)
            leaf_nodes.add(node_id)
            width = 1.0
            x_center = x_start + 0.5
            if node_id != "_world_":
                node_pos[node_id] = (x_center, -float(depth))
            return width, x_center
        
        # Internal node
        current_x = x_start
        child_xs = []
        
        for child in structural_children:
            w, cx = parse_and_layout(child, depth + 1, current_x)
            current_x += w
            child_xs.append(cx)
            
            # Store edge info
            uplink = child.get("uplink")
            bw = child.get("bandwidth", 0)
            if uplink:
                edges.append((node_id, child["id"], uplink, bw))
                
        width = current_x - x_start
        x_center = sum(child_xs) / len(child_xs)
        if node_id != "_world_":
            node_pos[node_id] = (x_center, -float(depth))
        
        return width, x_center

    parse_and_layout(topology, 0, 0)

    # 2. Process Metrics
    # Calculate average utilization for each link
    avg_usage = sim_run.get_average_link_utilization()

    # 3. Plot
    fig = plt.figure(figsize=(14, 10))
    ax = plt.gca()

    # Draw Edges
    for u, v, name, bw in edges:
        if u not in node_pos or v not in node_pos:
            continue
            
        ux, uy = node_pos[u]
        vx, vy = node_pos[v]
        
        # Determine Usage
        up_bps = avg_usage.get(f"{name}_UP", 0.0)
        down_bps = avg_usage.get(f"{name}_DOWN", 0.0)
        
        # Fallback if specific UP/DOWN not found but base name exists
        if up_bps == 0 and down_bps == 0 and name in avg_usage:
            up_bps = avg_usage[name]
        
        # metric is in bytes/s, bandwidth is in bits/s
        # Convert average usage to bits/s
        up_bps *= 8.0
        down_bps *= 8.0

        up_pct = (up_bps / bw * 100.0) if bw > 0 else 0.0
        down_pct = (down_bps / bw * 100.0) if bw > 0 else 0.0
        util_pct = max(up_pct, down_pct)
        
        bw_str = format_bw(bw)
        label_text = f"{bw_str}\\nU: {up_pct:.0f}% / D: {down_pct:.0f}%"
        
        # Style
        color = "#aaaaaa" # default gray
        width = 1.0
        font_weight = 'normal'
        
        if util_pct > 80:
            color = "#d62728" # red
            width = 2.5
            font_weight = 'bold'
        elif util_pct > 50:
            color = "#ff7f0e" # orange
            width = 2.0
            font_weight = 'bold'
        elif util_pct > 10:
            color = "#1f77b4" # blue
            width = 1.5

        plt.plot([ux, vx], [uy, vy], color=color, linewidth=width, zorder=1)
        
        # Label with utilization
        mid_x = (ux + vx) / 2
        mid_y = (uy + vy) / 2
        
        plt.text(mid_x, mid_y, label_text, 
                 color='black', fontsize=7, ha='center', va='center',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))

    # Draw Nodes
    for nid, (x, y) in node_pos.items():
        plt.scatter(x, y, s=150, c='lightgray', zorder=2, edgecolors='gray')
        # Offset label to avoid overlap with dot
        label = node_names.get(nid, nid)
        
        if nid in leaf_nodes:
            # Place label BELOW the node
            plt.text(x, y - 0.15, label, fontsize=9, ha='center', va='top', fontweight='bold')
        else:
            # Place label ABOVE the node (default)
            plt.text(x, y + 0.15, label, fontsize=9, ha='center', va='bottom', fontweight='bold')

    plt.title("Network Topology & Link Utilization")
    plt.axis('off')
    add_git_hash(fig)
    
    # Save
    output_path = os.path.join(output_dir, "network_topology.svg")
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Saved plot to {output_path}")
    plt.close()

def viz_star_link_utilization(sim_run):
    links = sim_run.get_star_links()
    if not links:
        print("No star links found or topology missing")
        return

    # If exactly two data centers (links) are connected, only draw the first one.
    if len(links) == 2:
        print("Exactly 2 star links found, visualizing only the first one.")
        links = links[:1]
        
    # Get all resource names we care about
    resource_names = []
    for link_name, _, _ in links:
        resource_names.append(f"{link_name}_UP")
        resource_names.append(f"{link_name}_DOWN")

    data = sim_run.get_link_utilization_over_time(resource_names)
    
    if data.empty:
        print("No metrics found for star links")
        return

    fig = plt.figure(figsize=(12, 8))
    
    time_vals = data.index
    
    # Needs reindexing explicitly if we want to ensure continuity, but pandas plot usually handles line gaps okay-ish
    # or we can reindex like before if strictly needed.
    
    for link_name, bw_bps, display_name in links:
        for direction in ["UP", "DOWN"]:
            resource_name = f"{link_name}_{direction}"
            if resource_name not in data.columns:
                continue
                
            s = data[resource_name]
            # Convert values (bytes/s) to utilization % of bandwidth (bits/s)
            utilization = (s * 8.0 / bw_bps) * 100.0
            
            label = f"{display_name} {direction}"
            linestyle = '-' if direction == 'UP' else '--'
            
            plt.plot(time_vals, utilization, label=label, linestyle=linestyle)

    plt.xlabel("Time (s)")
    plt.ylabel("Bandwidth Utilization (%)")
    plt.title("Star Link Utilization Over Time")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    add_git_hash(fig)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, "star_link_utilization.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()

def main():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # infer topology path relative to net_metrics_file if possible
    # logic copied from old viz_net_metrics but applied centrally
    metrics_dir = os.path.dirname(net_metrics_file)
    topo_path = os.path.join(metrics_dir, "topology.json")

    sim_run = SimulationRun(
        mon_path=mon_metrics_file,
        net_path=net_metrics_file,
        client_path=client_metrics_file,
        topology_path=topo_path
    )

    viz_mon_metrics(sim_run)
    viz_client_metrics(sim_run)
    viz_net_metrics(sim_run)
    viz_star_link_utilization(sim_run)

if __name__ == "__main__":
    main()
