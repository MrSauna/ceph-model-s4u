import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
import textwrap

# Ensure tools is in path
try:
    from tools.sim_analysis import SimulationRun
except ImportError:
    sys.path.append(os.getcwd())
    from tools.sim_analysis import SimulationRun


def main():
    if "snakemake" in globals():
        client_metrics_files = snakemake.input.client_metrics
        mon_metrics_files = snakemake.input.mon_metrics
        net_metrics_files = snakemake.input.net_metrics
        output_dir = snakemake.params.output_dir
        git_hash = snakemake.params.get("git_hash", "unknown")
    else:
        # Fallback for manual testing (if set by wrapper) or sys.argv could be added here
        try:
             # Check if we have manually injected variables in the module scope (from wrapper)
             # This is a bit hacky but keeps the wrapper working if we update it
             client_metrics_files = globals().get("client_metrics_files")
             mon_metrics_files = globals().get("mon_metrics_files")
             net_metrics_files = globals().get("net_metrics_files")
             output_dir = globals().get("output_dir")
             git_hash = globals().get("git_hash", "unknown")
             
             if not (client_metrics_files and mon_metrics_files and net_metrics_files and output_dir):
                 raise ValueError("Missing inputs")
                 
        except (NameError, ValueError):
            print("This script is intended to be run via Snakemake or with proper globals injected.")
            sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # Helper for git hash
    def add_git_hash(fig):
        fig.text(0.99, 0.01, f'commit: {git_hash}', 
                 ha='right', va='bottom', 
                 fontsize=8, color='gray', alpha=0.5)

    original_labels = []
    raw_labels = []
    run_objects = []

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
        original_labels.append(label)
        raw_labels.append(label)
        run_objects.append(run)

    import re
    parsed_runs = []
    profiles = set()
    
    for orig_label, label, run in zip(original_labels, raw_labels, run_objects):
        sat_time = run.get_saturation_time(["dc-0_uplink_DOWN", "dc-0_uplink_UP"], threshold=0.8)
        
        m = re.match(r'big_cluster_add_(\d+)_(.*)', orig_label)
        if m:
            n_val = int(m.group(1))
            prof = m.group(2)
            parsed_runs.append({
                "type": "grid",
                "n": n_val,
                "profile": prof,
                "orig_label": orig_label,
                "label": label,
                "run": run,
                "sat_time": sat_time
            })
            profiles.add(prof)
        else:
            parsed_runs.append({
                "type": "other",
                "orig_label": orig_label,
                "label": label,
                "run": run,
                "sat_time": sat_time
            })
            
    # New Grid Visualization
    grid_runs = [r for r in parsed_runs if r["type"] == "grid"]
    other_runs = [r for r in parsed_runs if r["type"] == "other"]
    
    if grid_runs or other_runs:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [3, 1]})
        ax_grid, ax_others = axes
        
        if grid_runs:
            df_grid = pd.DataFrame(grid_runs)
            pivot = df_grid.pivot(index="n", columns="profile", values="sat_time")
            pivot.plot(kind="line", marker="o", ax=ax_grid)
            ax_grid.set_xlabel("N (Added Scale)")
            ax_grid.set_ylabel("Saturation Time (s)")
            ax_grid.set_title("Cross-DC Link Saturation (>80%)")
            ax_grid.grid(True)
            
        if other_runs:
            df_others = pd.DataFrame(other_runs)
            # Remove prefix for display
            df_others["display_label"] = df_others["orig_label"].str.replace("big_cluster_", "")
            df_others.plot(kind="bar", x="display_label", y="sat_time", ax=ax_others, color="gray", legend=False)
            ax_others.set_xlabel("Other Scenarios")
            ax_others.set_ylabel("Saturation Time (s)")
            ax_others.set_title("Other Scenarios Saturation")
            for tick in ax_others.get_xticklabels():
                tick.set_rotation(45)
                tick.set_ha('right')
            ax_others.grid(True, axis='y')
            
        plt.tight_layout()
        add_git_hash(fig)
        out_path = os.path.join(output_dir, "link_saturation.svg")
        plt.savefig(out_path)
        print(f"Saved link saturation visualization to {out_path}")
        plt.close()

    # Filter for standard plots
    selected_runs = []
    if profiles:
        for prof in profiles:
            prof_runs = [r for r in parsed_runs if r["type"] == "grid" and r["profile"] == prof]
            prof_runs.sort(key=lambda x: x["n"])
            valid_runs = [r for r in prof_runs if r["sat_time"] <= 0]
            if valid_runs:
                selected_runs.append(valid_runs[-1])
            elif prof_runs:
                selected_runs.append(prof_runs[0])
                
        for r in parsed_runs:
            if r["type"] == "other":
                selected_runs.append(r)
                
        selected_runs.sort(key=lambda x: x["orig_label"])
        raw_labels = [r["label"] for r in selected_runs]
        run_objects = [r["run"] for r in selected_runs]

    # Shorten labels by removing common prefix
    if raw_labels:
        # Pass 1: Global common prefix
        prefix = os.path.commonprefix(raw_labels)
        if len(prefix) > 4:
             raw_labels = [l[len(prefix):].lstrip(" _-") for l in raw_labels]
        
        # Pass 2: Common prefix of non-empty labels (e.g. 'clients_')
        non_empty = [l for l in raw_labels if l]
        if non_empty:
             prefix2 = os.path.commonprefix(non_empty)
             # Only strip if it covers a significant portion or is a specific marker
             if len(prefix2) > 2:
                 raw_labels = [l[len(prefix2):].lstrip(" _-") if l else l for l in raw_labels]
                 
        # Finalize
        final_labels = []
        for l in raw_labels:
            if not l:
                l = "baseline"
            final_labels.append(l)
        raw_labels = final_labels

    # Combine with runs and wrap text
    runs = []
    for label, run in zip(raw_labels, run_objects):
        # Wrap label
        wrapped_label = "\n".join(textwrap.wrap(label, width=12, break_long_words=False))
        runs.append((wrapped_label, run))

    # Plot 1: Compare Client Throughput (Read and Write)
    fig, (ax_read, ax_write) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    
    for label, run in runs:
        tp = run.get_client_throughput(window_size=10)
        if tp.empty:
            continue
            
        start_time, end_time = run.get_recovery_times()
        
        # Plot Read Throughput
        if 'read' in tp.columns:
            line_r, = ax_read.plot(tp.index, tp['read'].values, label=label, linewidth=1.5)
            color_r = line_r.get_color()
            
            if start_time is not None:
                ax_read.axvline(x=start_time, color=color_r, linestyle=':', alpha=0.8)
            if end_time is not None:
                ax_read.axvline(x=end_time, color=color_r, linestyle='--', alpha=0.8)
                
        # Plot Write Throughput
        if 'write' in tp.columns:
            line_w, = ax_write.plot(tp.index, tp['write'].values, label=label, linewidth=1.5)
            color_w = line_w.get_color()
            
            if start_time is not None:
                ax_write.axvline(x=start_time, color=color_w, linestyle=':', alpha=0.8)
            if end_time is not None:
                ax_write.axvline(x=end_time, color=color_w, linestyle='--', alpha=0.8)
        
    ax_read.set_xlabel("Time (s)")
    ax_read.set_ylabel("Read Throughput (MiB/s)")
    ax_read.set_title("Client Read Throughput")
    ax_read.legend()
    ax_read.grid(True)
    
    ax_write.set_xlabel("Time (s)")
    ax_write.set_ylabel("Write Throughput (MiB/s)")
    ax_write.set_title("Client Write Throughput")
    ax_write.legend()
    ax_write.grid(True)

    plt.tight_layout()
    add_git_hash(fig)
    
    output_path = os.path.join(output_dir, "client_throughput.svg")
    plt.savefig(output_path)
    print(f"Saved comparison to {output_path}")
    plt.close()

    # Plot 2: Recovery Detailed Progress
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 1]})
    
    # Calculate global recovery window
    global_start = float('inf')
    global_end = float('-inf')
    has_recovery_data = False
    
    # Helper to collect data first so we can calculate offsets
    recovery_series = []
    max_concurrency = 0

    for label, run in runs:
        start_time, end_time = run.get_recovery_times()
        if start_time is None:
            continue
            
        has_recovery_data = True
        global_start = min(global_start, start_time)
        if end_time is not None:
             global_end = max(global_end, end_time)

        # Prepare Data
        mon_df = run.mon_df
        # Filter for relevant range 
        df_rec = mon_df[mon_df["time"] >= start_time].copy()
        
        if df_rec.empty:
            continue

        # Calculate Backlog
        df_rec["backlog"] = df_rec["backfill"] + df_rec["backfill_wait"]
        
        # Track max concurrency for offset calculation
        if not df_rec["backfill"].empty:
             max_concurrency = max(max_concurrency, df_rec["backfill"].max())
        
        recovery_series.append({
            "label": label,
            "df": df_rec
        })

    if has_recovery_data:
        # Determine offset for seismic plot
        # Add a bit of padding, e.g., 20% or at least 1
        offset_step = max(1, max_concurrency * 1.2)
        
        for i, series in enumerate(recovery_series):
            label = series["label"]
            df_rec = series["df"]
            
            # Plot Backlog on Top Axis (ax1)
            l1, = ax1.plot(df_rec["time"], df_rec["backlog"], label=label, linewidth=1.5)
            color = l1.get_color()
            
            # Plot Active Backfill on Bottom Axis (ax2)
            # Offset = i * offset_step
            base_y = i * offset_step
            # We want the baseline to be at base_y, and values adding to it
            # But "seismic" usually means centered or just displaced. 
            # The user asked for "starting from different height".
            # Let's plot (y + base_y)
            
            y_values = df_rec["backfill"] + base_y
            
            # Add a baseline for reference? Or just the trace?
            # User said "discrete behavior... quantized to integers". 
            # Step plot is good.
            ax2.plot(df_rec["time"], y_values, label=label, 
                     color=color, linestyle='-', drawstyle='steps-post', linewidth=1.5)
            
            # Maybe add a faint zero line for this series?
            ax2.axhline(y=base_y, color=color, linestyle=':', alpha=0.3, linewidth=0.5)

        # Set X-Axis Limits
        margin = (global_end - global_start) * 0.05
        ax1.set_xlim(global_start - margin, global_end + margin)
        # ax2 shares x, so it will update automatically
        
        ax1.set_ylabel("Total Backlog (PGs)")
        ax1.set_title("Recovery Progress: Total Backlog")
        ax1.grid(True)
        ax1.legend(loc='upper right')
        
        ax2.set_xlabel("Time (s)")
        ax2.set_ylabel("Active Backfill (Stacked)")
        ax2.set_title("Recovery Progress: Active Backfill (Offset per Scenario)")
        # Remove Y ticks on the second plot as they are arbitrary due to offset
        ax2.set_yticks([])
        
        # Maybe add text labels for each series on the left or right of the plot?
        for i, series in enumerate(recovery_series):
             label = series["label"]
             y_pos = i * offset_step + (offset_step / 2) # approximate center of the band
             # Or just at the base?
             y_pos = i * offset_step
             # Place label at the start time?
             # Let's just rely on the legend for color matching, 
             # OR put text at the far right.
             ax2.text(global_end + margin, y_pos, label, va='bottom', ha='right', 
                      fontsize=8, color='black', alpha=0.7)

        add_git_hash(fig)
        
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
        
        # Flatten nested stats
        flat_stats = {}
        for op_type, type_stats in stats.items():
            if isinstance(type_stats, dict):
                for metric, val in type_stats.items():
                    key = f"{metric} ({op_type})" if op_type != 'all' else metric
                    flat_stats[key] = val
            else:
                flat_stats[op_type] = type_stats

        # Skip cases with no actual client traffic (all values zero)
        values = [v for k, v in flat_stats.items() if k != "label"]
        if not values or all(v == 0 for v in values):
            continue
            
        flat_stats["label"] = label
        latency_data.append(flat_stats)
        labels.append(label)
        
    if latency_data:
        df_lat = pd.DataFrame(latency_data)
        df_lat = df_lat.set_index("label")
        
        # Reorder columns for logical progression
        cols = [
            "avg (read)", "avg (write)", 
            "p95 (read)", "p95 (write)",
            "p99 (read)", "p99 (write)",
            "p99.5 (read)", "p99.5 (write)",
            "avg", "p95", "p99", "p99.5"
        ]
        # Filter strictly those we put in
        cols = [c for c in cols if c in df_lat.columns]
        df_lat = df_lat[cols]
        
        # Determine if we should display in milliseconds based on max value
        max_latency = df_lat.max().max()
        use_ms = max_latency < 1.0  # Use ms if all values are under 1 second
        
        if use_ms:
            # Convert to milliseconds for display
            df_lat = df_lat * 1000
            latency_unit = "ms"
        else:
            latency_unit = "s"
        
        def format_latency(val):
            """Format latency value with appropriate precision."""
            if val >= 100:
                return f'{val:.0f}'
            elif val >= 10:
                return f'{val:.1f}'
            elif val >= 1:
                return f'{val:.2f}'
            else:
                return f'{val:.3f}'
        
        # Split into Read and Write dataframes if columns exist
        read_cols = [c for c in cols if "(read)" in c]
        write_cols = [c for c in cols if "(write)" in c]
        
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6), sharey=True)
        
        # Plot Reads
        if read_cols:
            df_read = df_lat[read_cols]
            # Rename columns to remove " (read)" for cleaner legend
            df_read.columns = [c.replace(" (read)", "") for c in df_read.columns]
            df_read.plot(kind='bar', ax=axes[0], rot=0)
            
            for p in axes[0].patches:
                height = p.get_height()
                if height > 0:
                    axes[0].annotate(format_latency(height), 
                                (p.get_x() + p.get_width() / 2., height), 
                                ha='center', va='bottom', 
                                fontsize=9, xytext=(0, 2), 
                                textcoords='offset points',
                                rotation=90)
            
            axes[0].set_xlabel("Scenarios")
            axes[0].set_ylabel(f"Latency ({latency_unit})")
            axes[0].set_title("Read Latency")
            axes[0].grid(True, axis='y')
            axes[0].legend(title="Metric")
        
        # Plot Writes
        if write_cols:
            df_write = df_lat[write_cols]
            df_write.columns = [c.replace(" (write)", "") for c in df_write.columns]
            df_write.plot(kind='bar', ax=axes[1], rot=0)
            
            for p in axes[1].patches:
                height = p.get_height()
                if height > 0:
                    axes[1].annotate(format_latency(height), 
                                (p.get_x() + p.get_width() / 2., height), 
                                ha='center', va='bottom', 
                                fontsize=9, xytext=(0, 2), 
                                textcoords='offset points',
                                rotation=90)
            
            axes[1].set_xlabel("Scenarios")
            axes[1].set_title("Write Latency")
            axes[1].grid(True, axis='y')
            axes[1].legend(title="Metric")

        # Rotate x labels if too many or long
        for ax in axes:
            ax.tick_params(axis='x', rotation=0)

        plt.suptitle("Client Latency Statistics")
        plt.tight_layout()
        add_git_hash(fig)
        
        output_path_lat = os.path.join(output_dir, "latency_comparison.svg")
        plt.savefig(output_path_lat)
        print(f"Saved latency visualization to {output_path_lat}")
        plt.close()
    else:
        print("No latency data found.")

    # Plot 4: IOPS Comparison
    iops_data = []
    iops_labels = []
    
    for label, run in runs:
        stats = run.get_iops_stats()
        if not stats:
            continue
            
        stats["label"] = label
        iops_data.append(stats)
        iops_labels.append(label)
        
    if iops_data:
        df_iops = pd.DataFrame(iops_data)
        df_iops = df_iops.set_index("label")
        
        # Only show read and write IOPS
        cols = ["read", "write"]
        cols = [c for c in cols if c in df_iops.columns]
        df_iops = df_iops[cols]
        
        # Drop columns where all values are zero (no traffic of that type)
        df_iops = df_iops.loc[:, (df_iops != 0).any()]
        
        if df_iops.empty:
            print("No non-zero IOPS data found, skipping IOPS plot.")
            return
        
        fig = plt.figure(figsize=(12, 6))
        ax = plt.gca()
        
        df_iops.plot(kind='bar', ax=ax, rot=0)
        
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.annotate(f'{height:.0f}', 
                            (p.get_x() + p.get_width() / 2., height), 
                            ha='center', va='bottom', 
                            fontsize=9, xytext=(0, 2), 
                            textcoords='offset points')

        plt.xlabel("Scenarios")
        plt.ylabel("IOPS")
        plt.title("Client IOPS Comparison")
        plt.grid(True, axis='y')
        plt.legend(title="Metric")
        
        # Horizontal x labels
        plt.xticks(rotation=0, ha='center')
        
        plt.tight_layout()
        add_git_hash(fig)
        
        output_path_iops = os.path.join(output_dir, "iops_comparison.svg")
        plt.savefig(output_path_iops)
        print(f"Saved IOPS visualization to {output_path_iops}")
        plt.close()
    else:
        print("No IOPS data found.")

if __name__ == "__main__":
    main()