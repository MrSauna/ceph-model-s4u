import os
import matplotlib.pyplot as plt
import pandas as pd


if "snakemake" not in globals():
    mon_metrics_file = "results/vscode/mon_metrics.csv"
    net_metrics_file = "results/vscode/net_metrics.csv"
    client_metrics_file = "results/vscode/client_metrics.csv"
    output_dir = "results/vscode/figures"
else:
    mon_metrics_file = snakemake.input.mon_metrics
    net_metrics_file = snakemake.input.net_metrics
    client_metrics_file = snakemake.input.client_metrics
    output_dir = snakemake.output


def viz_mon_metrics(mon_metrics_df):
    # try/except block removed as reading happens in main
    if mon_metrics_df is None:
        return

    total_pgs = len(mon_metrics_df)
    
    times = [0] + mon_metrics_df["time"].tolist()
    pgs_left = [total_pgs]
    
    current_pgs = total_pgs
    for _ in range(len(mon_metrics_df)):
        current_pgs -= 1
        pgs_left.append(current_pgs)
        
    plt.figure(figsize=(10, 6))
    plt.plot(times, pgs_left, marker='o', linestyle='-', markersize=2)
    plt.xlabel("Time (s)")
    plt.ylabel("PGs Left")
    plt.title("PG Completion Over Time")
    plt.grid(True)
    
    output_path = os.path.join(output_dir, "mon_plot.svg")
    plt.savefig(output_path)
    print(f"Saved plot to {output_path}")
    plt.close()


def viz_client_metrics(client_metrics_df):
    if client_metrics_df.empty:
        print("Client metrics empty, skipping viz_client_metrics")
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

def main():
    mon_metrics_df = pd.read_csv(mon_metrics_file)
    client_metrics_df = pd.read_csv(client_metrics_file)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    viz_mon_metrics(mon_metrics_df)
    viz_client_metrics(client_metrics_df)


if __name__ == "__main__":
    main()