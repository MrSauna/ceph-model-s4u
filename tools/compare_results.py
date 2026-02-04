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
        
        plt.plot(total_tp.index, total_tp.values, label=label, linewidth=1.5)
        
    plt.xlabel("Time (s)")
    plt.ylabel("Total Throughput (MiB/s)")
    plt.title("Throughput Comparison")
    plt.legend()
    plt.grid(True)
    
    output_path = os.path.join(output_dir, "comparison_throughput.svg")
    plt.savefig(output_path)
    print(f"Saved comparison to {output_path}")
    plt.close()

if __name__ == "__main__":
    main()