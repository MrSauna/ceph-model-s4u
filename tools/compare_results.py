import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys

# Ensure tools is in path
try:
    from tools.sim_analysis import SimulationRun
except ImportError:
    sys.path.append(os.getcwd())
    from tools.sim_analysis import SimulationRun

if "snakemake" not in globals():
    # Standalone mode - mostly for testing or if manually invoked
    # This might need argument parsing if used manually, but for now we assume snakemake
    print("This script is intended to be run via Snakemake.")
    sys.exit(1)
    
client_metrics_files = snakemake.input.client_metrics
mon_metrics_files = snakemake.input.mon_metrics
net_metrics_files = snakemake.input.net_metrics
output_dir = snakemake.params.output_dir

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

def main():
    runs = []
    
    # Zip together the file lists
    # We assume they are strictly ordered and correspond to each other
    for c_file, m_file, n_file in zip(client_metrics_files, mon_metrics_files, net_metrics_files):
        # Use directory name as label, or part of it
        dirname = os.path.dirname(c_file)
        # simplistic label: last component of directory
        label = os.path.basename(dirname) 
        if label == "vscode": # fallback if run in vscode dev setup
             label = "current_run"
             
        # Try to find topology
        topo_file = os.path.join(dirname, "topology.json")
        
        run = SimulationRun(
            mon_path=m_file,
            net_path=n_file,
            client_path=c_file,
            topology_path=topo_file
        )
        runs.append((label, run))

    # Plot 1: Compare Total Throughput
    fig = plt.figure(figsize=(12, 6))
    
    for label, run in runs:
        tp = run.get_client_throughput(window_size=10)
        if tp.empty:
            continue
            
        # Calculate total throughput (read + write)
        total_tp = tp.sum(axis=1)
        
        line, = plt.plot(total_tp.index, total_tp.values, label=label, linewidth=1.5)
        color = line.get_color()
        
        # Add recovery markers
        start_time, end_time = run.get_recovery_times()
        if start_time is not None:
            plt.axvline(x=start_time, color=color, linestyle=':', alpha=0.8)
        if end_time is not None:
            plt.axvline(x=end_time, color=color, linestyle='--', alpha=0.8)
        
    plt.xlabel("Time (s)")
    plt.ylabel("Total Throughput (MiB/s)")
    plt.title("Client Throughput Comparison")
    plt.legend()
    plt.grid(True)
    
    output_path = os.path.join(output_dir, "client_throughput.svg")
    plt.savefig(output_path)
    print(f"Saved comparison to {output_path}")
    plt.close()

    # Plot 2: Recovery Detailed Progress
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax2 = ax1.twinx()
    
    # Calculate global recovery window
    global_start = float('inf')
    global_end = float('-inf')
    has_recovery_data = False
    
    for label, run in runs:
        start_time, end_time = run.get_recovery_times()
        if start_time is None:
            continue
            
        has_recovery_data = True
        global_start = min(global_start, start_time)
        # If end_time is None (ongoing), distinct logic might be needed, 
        # but get_recovery_times returns max time if not finished.
        if end_time is not None:
             global_end = max(global_end, end_time)

        # Prepare Data
        mon_df = run.mon_df
        # Filter for relevant range (optional, but good for plotting)
        # We assume monotonic time
        df_rec = mon_df[mon_df["time"] >= start_time].copy()
        
        if df_rec.empty:
            continue

        # Calculate Backlog
        df_rec["backlog"] = df_rec["backfill"] + df_rec["backfill_wait"]
        
        # Plot Primary (Backlog) on Left Axis
        # We need to maintain color cycle manually or retrieve it
        # Since we are iterating, we can trust matplotlib's cycle if we plot on same ax, 
        # but here we are plotting on two axes.
        # Let's verify throughput plot order or just rely on default cycle order matching loop order
        
        # Plot Backlog
        l1, = ax1.plot(df_rec["time"], df_rec["backlog"], label=f"{label} (Backlog)", linewidth=1.5)
        color = l1.get_color()
        
        # Plot Active Backfill on Right Axis (Stairs)
        # Match color to the backlog line
        ax2.plot(df_rec["time"], df_rec["backfill"], label=f"{label} (Active)", 
                 color=color, linestyle=':', drawstyle='steps-post', linewidth=1.5, alpha=0.8)

    if has_recovery_data:
        # Set X-Axis Limits
        # Add a small buffer?
        margin = (global_end - global_start) * 0.05
        ax1.set_xlim(global_start - margin, global_end + margin)
        
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Total Backlog (PGs)")
        ax2.set_ylabel("Active Backfill (PGs)")
        
        # Scale ax2 so that the max value is at 10% of height
        # Get max data value plotted on ax2
        # We didn't track it, so let's get it from the lines if possible, or recalculate
        max_backfill = 0
        for line in ax2.get_lines():
             y_data = line.get_ydata()
             if len(y_data) > 0:
                 max_backfill = max(max_backfill, max(y_data))
        
        if max_backfill > 0:
            ax2.set_ylim(0, max_backfill * 10)
        
        # Scale ax1 so that 0 is at 10% of height (to avoid overlap with ax2)
        # 0 = bottom + 0.1 * (top - bottom)
        # => -bottom = 0.1 * top - 0.1 * bottom
        # => -0.9 * bottom = 0.1 * top
        # => bottom = -top / 9
        
        max_backlog = 0
        for line in ax1.get_lines():
             y_data = line.get_ydata()
             if len(y_data) > 0:
                 max_backlog = max(max_backlog, max(y_data))
        
        if max_backlog > 0:
            # Add a bit of top margin (e.g. 5%)
            top_limit = max_backlog * 1.05
            bottom_limit = -top_limit / 9.0
            ax1.set_ylim(bottom_limit, top_limit)
        
        ax1.set_title("Recovery Progress: Backlog vs Active Backfill")
        ax1.grid(True)
        
        # Combine legends?
        # It's tricky with dual axes. Let's just put them separately or let user figure it out by color.
        # Or distinct locations.
        ax1.legend(loc='upper left')
        ax2.legend(loc='upper right')
        
        output_path_rec = os.path.join(output_dir, "recovery_progress.svg")
        plt.savefig(output_path_rec)
        print(f"Saved recovery visualization to {output_path_rec}")
    else:
        print("No recovery data found, skipping recovery plot.")
    
    plt.close()

    # Plot 3: Latency Comparison
    latency_data = []
    labels = []
    
    for label, run in runs:
        stats = run.get_latency_stats()
        if not stats:
            continue
            
        stats["label"] = label
        latency_data.append(stats)
        labels.append(label)
        
    if latency_data:
        df_lat = pd.DataFrame(latency_data)
        df_lat = df_lat.set_index("label")
        
        # Reorder columns for logical progression
        cols = ["avg", "p95", "p99", "p99.5"]
        # Filter strictly those we put in
        cols = [c for c in cols if c in df_lat.columns]
        df_lat = df_lat[cols]
        
        fig = plt.figure(figsize=(12, 6))
        ax = plt.gca()
        
        df_lat.plot(kind='bar', ax=ax, rot=0)
        
        plt.xlabel("Scenario")
        plt.ylabel("Latency (s)")
        plt.title("Client Latency Statistics")
        plt.grid(True, axis='y')
        plt.legend(title="Metric")
        
        # Rotate x labels if too many or long
        plt.xticks(rotation=45, ha='right')
        
        plt.tight_layout()
        
        output_path_lat = os.path.join(output_dir, "latency_comparison.svg")
        plt.savefig(output_path_lat)
        print(f"Saved latency visualization to {output_path_lat}")
        plt.close()
    else:
        print("No latency data found.")

if __name__ == "__main__":
    main()