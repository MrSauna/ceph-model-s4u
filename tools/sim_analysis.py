import os
import pandas as pd
import json

class SimulationRun:
    """
    Represents a single simulation run with lazy-loaded metrics.
    
    Uses pre-aggregated files from C++ simulation for memory efficiency:
    - client_metrics_agg.csv: per-second throughput buckets
    - client_latency_summary.csv: pre-computed percentiles and stats
    - client_latency_digest.csv: T-Digest centroids for custom percentile queries
    """
    
    def __init__(self, mon_path, net_path, client_path, topology_path=None):
        self.mon_path = mon_path
        self.net_path = net_path
        self.client_path = client_path
        self.topology_path = topology_path
        
        # Derive aggregated file paths from client_path directory
        self.dir_path = os.path.dirname(client_path)
        self.agg_path = os.path.join(self.dir_path, "client_metrics_agg.csv")
        self.summary_path = os.path.join(self.dir_path, "client_latency_summary.csv")
        self.digest_path = os.path.join(self.dir_path, "client_latency_digest.csv")
        
        self._mon_df = None
        self._net_df = None
        self._agg_df = None
        self._topology = None
        self._latency_stats_cache = None
        self._iops_stats_cache = None

    @classmethod
    def from_directory(cls, dir_path):
        """
        Initialize from a standard results directory structure.
        """
        return cls(
            mon_path=os.path.join(dir_path, "mon_metrics.csv"),
            net_path=os.path.join(dir_path, "net_metrics.csv"),
            client_path=os.path.join(dir_path, "client_metrics.csv"),
            topology_path=os.path.join(dir_path, "topology.json")
        )

    @property
    def mon_df(self):
        if self._mon_df is None and os.path.exists(self.mon_path):
            self._mon_df = pd.read_csv(self.mon_path)
        return self._mon_df

    @property
    def net_df(self):
        if self._net_df is None and os.path.exists(self.net_path):
            self._net_df = pd.read_csv(self.net_path)
        return self._net_df
    
    @property
    def agg_df(self):
        """Load aggregated client metrics (tiny file, ~KB)."""
        if self._agg_df is None and os.path.exists(self.agg_path):
            self._agg_df = pd.read_csv(self.agg_path)
        return self._agg_df

    @property
    def topology(self):
        if self._topology is None and self.topology_path and os.path.exists(self.topology_path):
            with open(self.topology_path, "r") as f:
                self._topology = json.load(f)
        return self._topology

    def get_client_throughput(self, resample_interval="1s", window_size=10):
        """
        Returns a smoothed DataFrame of client read/write throughput in MiB/s.
        Uses pre-aggregated data from C++ simulation.
        """
        df = self.agg_df
        if df is None or df.empty:
            return pd.DataFrame()
        
        # Convert bytes to MiB/s (already per-second buckets)
        result = pd.DataFrame(index=df['time'])
        
        if 'read_bytes' in df.columns:
            result['read'] = df['read_bytes'].values / (1024 * 1024)
        if 'write_bytes' in df.columns:
            result['write'] = df['write_bytes'].values / (1024 * 1024)
        
        result.index = result.index.astype(float)
        result = result.sort_index()
        
        # Apply rolling mean for smoothing
        smoothed = result.rolling(window=window_size, min_periods=1).mean()
        
        return smoothed

    def get_average_link_utilization(self):
        """Returns a Series mapping resource_name -> average bytes/s."""
        df = self.net_df
        if df is None or df.empty:
            return pd.Series(dtype=float)
        return df.groupby("resource_name")["value"].mean()

    def get_link_utilization_over_time(self, link_names):
        """Returns a DataFrame of utilization over time for specified links."""
        df = self.net_df
        if df is None or df.empty:
            return pd.DataFrame()
        
        mask = df["resource_name"].isin(link_names)
        filtered = df[mask]
        
        if filtered.empty:
            return pd.DataFrame()

        pivoted = filtered.pivot(index="timestamp", columns="resource_name", values="value")
        return pivoted

    def get_star_links(self):
        """Identify links connected to the 'star' node."""
        topology = self.topology
        if not topology:
            return []

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
            return []
            
        links = []
        for child in star_node.get("children", []):
            uplink = child.get("uplink")
            bw = child.get("bandwidth", 0)
            name = child.get("name", child.get("id"))
            
            if uplink and bw > 0:
                links.append((uplink, bw, name))
        return links

    def get_recovery_times(self):
        """
        Returns a tuple (start_time, end_time) of the recovery process.
        """
        df = self.mon_df
        if df is None or df.empty:
            return None, None
            
        df["total_pgs"] = df["active_clean"] + df["backfill"] + df["backfill_wait"]
        max_pgs = df["total_pgs"].max()
        
        not_clean = df[df["active_clean"] < max_pgs]
        
        if not_clean.empty:
            return None, None
            
        start_time = not_clean["time"].min()
        
        clean_after_start = df[(df["time"] > start_time) & (df["active_clean"] == max_pgs)]
        
        if clean_after_start.empty:
            end_time = df["time"].max()
        else:
            end_time = clean_after_start["time"].min()
            
        return start_time, end_time

    def get_latency_stats(self):
        """
        Returns a dictionary with latency statistics from pre-computed summary.
        """
        if self._latency_stats_cache is not None:
            return self._latency_stats_cache
            
        if not os.path.exists(self.summary_path):
            return {}
            
        try:
            df = pd.read_csv(self.summary_path)
            stats = {}
            for _, row in df.iterrows():
                metric = row['metric']
                value = row['value']
                if metric == 'avg':
                    stats['avg'] = value
                elif metric == 'p95':
                    stats['p95'] = value
                elif metric == 'p99':
                    stats['p99'] = value
                elif metric == 'p99.5':
                    stats['p99.5'] = value
            self._latency_stats_cache = stats
            return stats
        except Exception:
            return {}

    def get_iops_stats(self):
        """
        Returns IOPS statistics computed from aggregated data.
        """
        if self._iops_stats_cache is not None:
            return self._iops_stats_cache
            
        df = self.agg_df
        if df is None or df.empty:
            return {}
        
        duration = df['time'].max() - df['time'].min()
        if duration <= 0:
            duration = 1
        
        total_read_ops = df['read_ops'].sum()
        total_write_ops = df['write_ops'].sum()
        total_ops = total_read_ops + total_write_ops
        
        stats = {
            'avg': total_ops / duration,
            'read': total_read_ops / duration,
            'write': total_write_ops / duration
        }
        self._iops_stats_cache = stats
        return stats
