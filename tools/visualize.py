import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import json


if "snakemake" not in globals():
    mon_metrics_file = "results/vscode/mon_metrics.csv"
    net_metrics_file = "results/vscode/net_metrics.csv"
    client_metrics_file = "results/vscode/client_metrics.csv"
    output_dir = "results/vscode/figures"
else:
    mon_metrics_file = snakemake.input.mon_metrics
    net_metrics_file = snakemake.input.net_metrics
    client_metrics_file = snakemake.input.client_metrics
    output_dir = snakemake.params.output_dir


def viz_mon_metrics(mon_metrics_df):
    # try/except block removed as reading happens in main
    if mon_metrics_df is None:
        return

    plt.figure(figsize=(10, 6))
    
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
    
    output_path = os.path.join(output_dir, "pgmap_plot.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()


def viz_client_metrics(client_metrics_df):
    if client_metrics_df.empty:
        print("Client metrics empty, creating empty plot")
        plt.figure(figsize=(10, 6))
        plt.xlabel("Time (s)")
        plt.ylabel("Throughput (MiB/s)")
        plt.title("Client Throughput (No Data)")
        plt.grid(True)
        output_path = os.path.join(output_dir, "client_throughput.svg")
        plt.savefig(output_path)
        plt.close()
        return

    # Sort by time
    client_metrics_df = client_metrics_df.sort_values(by="time")
    
    # Convert bytes to MiB
    client_metrics_df["size_mib"] = client_metrics_df["op_size"] / (1024 * 1024)
    
    # Bin by second (floor) and calculate rolling average
    # We use a 1s resample periodicity, but smooth over a window (e.g., 5s)
    
    # Create a datetime-like index for resampling
    client_metrics_df["time_pd"] = pd.to_timedelta(client_metrics_df["time"], unit="s")
    client_metrics_df = client_metrics_df.set_index("time_pd")
    
    # Group by op_type and resample to 1s, summing size_mib
    # We unstack 'op_type' to get columns for 'read', 'write'
    grouped = client_metrics_df.groupby("op_type")["size_mib"].resample("1s").sum().unstack(level=0, fill_value=0)
    
    # Apply rolling mean to smooth out the "quantization" effect of large ops
    # Window size of 5s - 10s is usually good for readability
    window_size = 5
    smoothed = grouped.rolling(window=window_size, min_periods=1).mean()
    
    # Convert index back to seconds (float) for plotting
    plot_index = smoothed.index.total_seconds()
    
    plt.figure(figsize=(10, 6))
    
    if "read" in smoothed.columns:
        plt.plot(plot_index, smoothed["read"], label=f"Read ({window_size}s avg)", linestyle='-', linewidth=1.5)
    
    if "write" in smoothed.columns:
        plt.plot(plot_index, smoothed["write"], label=f"Write ({window_size}s avg)", linestyle='--', linewidth=1.5)

    plt.xlabel("Time (s)")
    plt.ylabel("Throughput (MiB/s)")
    plt.title(f"Client Throughput Over Time (Smoothed {window_size}s)")
    plt.legend()
    plt.grid(True)
    
    output_path = os.path.join(output_dir, "client_throughput.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()

def viz_net_metrics(net_metrics_df):
    if net_metrics_df.empty:
        print("Net metrics empty, skipping viz_net_metrics")
        return

    # Infer topology path relative to net_metrics_file or default
    metrics_dir = os.path.dirname(net_metrics_file)
    topo_path = os.path.join(metrics_dir, "topology.json")

    if not os.path.exists(topo_path):
        print(f"Topology file not found: {topo_path}")
        return

    try:
        with open(topo_path, "r") as f:
            topology = json.load(f)
    except Exception as e:
        print(f"Error loading topology: {e}")
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
            name_str += "\n" + "\n".join(formatted_actors)
            
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
    # avg_val = mean(value)
    
    avg_usage = net_metrics_df.groupby("resource_name")["value"].mean()

    # 3. Plot
    plt.figure(figsize=(14, 10))
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
        label_text = f"{bw_str}\nU: {up_pct:.0f}% / D: {down_pct:.0f}%"
        
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
    
    # Save
    output_path = os.path.join(output_dir, "network_topology.svg")
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Saved plot to {output_path}")
    plt.close()

def viz_star_link_utilization(net_metrics_df):
    if net_metrics_df.empty:
        print("Net metrics empty, skipping viz_star_link_utilization")
        return

    # Infer topology path relative to net_metrics_file or default
    metrics_dir = os.path.dirname(net_metrics_file)
    topo_path = os.path.join(metrics_dir, "topology.json")

    if not os.path.exists(topo_path):
        print(f"Topology file not found: {topo_path}")
        return

    try:
        with open(topo_path, "r") as f:
            topology = json.load(f)
    except Exception as e:
        print(f"Error loading topology: {e}")
        return

    # Find star node and its children (data centers)
    # The topology structure is usually _world_ -> star -> dc-X
    
    star_node = None
    
    # Helper to find node by id recursively
    def find_node(root, target_id):
        if root.get("id") == target_id:
            return root
        for child in root.get("children", []):
            res = find_node(child, target_id)
            if res:
                return res
        return None

    star_node = find_node(topology, "star")
    
    if not star_node:
        print("Star node not found in topology, skipping star link viz")
        return
        
    # Identify links to monitor
    # child["uplink"] is the link name. child["bandwidth"] is the link capacity (bits/s)
    links = [] # (link_name, bandwidth, display_name)
    
    for child in star_node.get("children", []):
        uplink = child.get("uplink")
        bw = child.get("bandwidth", 0)
        name = child.get("name", child.get("id"))
        
        if uplink and bw > 0:
            links.append((uplink, bw, name))
            
    if not links:
        print("No star links found")
        return

    # If exactly two data centers (links) are connected, only draw the first one.
    if len(links) == 2:
        print("Exactly 2 star links found, visualizing only the first one.")
        links = links[:1]

    plt.figure(figsize=(12, 8))
    
    # Prepare data for plotting
    # net_metrics_df has columns: timestamp, resource_type, resource_name, value
    # value is in bytes/s (utilization)
    
    # Sort by time
    df = net_metrics_df.sort_values(by="timestamp")
    time_vals = df["timestamp"].unique()
    time_vals.sort()
    
    has_data = False
    
    # Define colors/linestyles for clarity
    # We might have UP and DOWN for each link
    
    for link_name, bw_bps, display_name in links:
        # Filter for this link
        # Look for _UP and _DOWN suffixes
        
        for direction in ["UP", "DOWN"]:
            resource_name = f"{link_name}_{direction}"
            link_data = df[df["resource_name"] == resource_name]
            
            if link_data.empty:
                continue
                
            has_data = True
            
            # Reindex to ensure we have values for all timestamps (fill missing with 0)
            # Use pivot or set_index to align with time_vals
            s = link_data.set_index("timestamp")["value"]
            s = s.reindex(time_vals, fill_value=0)
            
            # Convert values (bytes/s) to utilization % of bandwidth (bits/s)
            # value * 8 / bw_bps * 100
            utilization = (s * 8.0 / bw_bps) * 100.0
            
            label = f"{display_name} {direction}"
            linestyle = '-' if direction == 'UP' else '--'
            
            plt.plot(time_vals, utilization, label=label, linestyle=linestyle)

    if not has_data:
        print("No metrics found for star links")
        plt.close()
        return

    plt.xlabel("Time (s)")
    plt.ylabel("Bandwidth Utilization (%)")
    plt.title("Star Link Utilization Over Time")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, "star_link_utilization.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()

def main():
    mon_metrics_df = pd.read_csv(mon_metrics_file)
    client_metrics_df = pd.read_csv(client_metrics_file)
    net_metrics_df = pd.read_csv(net_metrics_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    viz_mon_metrics(mon_metrics_df)
    viz_client_metrics(client_metrics_df)
    viz_net_metrics(net_metrics_df)
    viz_star_link_utilization(net_metrics_df)

if __name__ == "__main__":
    main()