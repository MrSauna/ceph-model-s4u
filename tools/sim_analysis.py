import os
import pandas as pd
import json

class SimulationRun:
    def __init__(self, mon_path, net_path, client_path, topology_path=None):
        self.mon_path = mon_path
        self.net_path = net_path
        self.client_path = client_path
        self.topology_path = topology_path
        
        self._mon_df = None
        self._net_df = None
        self._client_df = None
        self._topology = None

    @classmethod
    def from_directory(cls, dir_path):
        """
        Initialize from a standard results directory structure.
        Assumes standard filenames: mon_metrics.csv, net_metrics.csv, client_metrics.csv, topology.json
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
    def client_df(self):
        if self._client_df is None and os.path.exists(self.client_path):
            self._client_df = pd.read_csv(self.client_path)
        return self._client_df

    @property
    def topology(self):
        if self._topology is None and self.topology_path and os.path.exists(self.topology_path):
            with open(self.topology_path, "r") as f:
                self._topology = json.load(f)
        return self._topology

    def get_client_throughput(self, resample_interval="1s", window_size=10):
        """
        Returns a smoothed DataFrame of client read/write throughput in MiB/s.
        Columns: 'read', 'write' (if present)
        Index: Time (float seconds)
        """
        df = self.client_df
        if df is None or df.empty:
            return pd.DataFrame()

        # Work on a copy to avoid mutating cache
        df = df.copy()
        
        # Convert bytes to MiB
        df["size_mib"] = df["op_size"] / (1024 * 1024)
        
        # Create a datetime-like index for resampling
        df["time_pd"] = pd.to_timedelta(df["time"], unit="s")
        df = df.set_index("time_pd")
        
        # Group by op_type and resample
        grouped = df.groupby("op_type")["size_mib"].resample(resample_interval).sum().unstack(level=0, fill_value=0)
        
        # Apply rolling mean
        smoothed = grouped.rolling(window=window_size, min_periods=1).mean()
        
        # Convert index back to seconds (float)
        smoothed.index = smoothed.index.total_seconds()
        
        return smoothed

    def get_average_link_utilization(self):
        """
        Returns a Series mapping resource_name -> average bytes/s
        """
        df = self.net_df
        if df is None or df.empty:
            return pd.Series(dtype=float)
        
        return df.groupby("resource_name")["value"].mean()

    def get_link_utilization_over_time(self, link_names):
        """
        Returns a DataFrame of utilization (bytes/s) over time for specified links.
        Index: timestamp
        Columns: resource_name
        """
        df = self.net_df
        if df is None or df.empty:
            return pd.DataFrame()
        
        # Filter for relevant links
        mask = df["resource_name"].isin(link_names)
        filtered = df[mask]
        
        if filtered.empty:
            return pd.DataFrame()

        # Pivot
        pivoted = filtered.pivot(index="timestamp", columns="resource_name", values="value")
        # Fill missing intermediate timestamps with 0? 
        # For now, let's just return what we have, resampling might be better but let's stick to what visualize.py did (it mostly just plotted what was there, but visualization 4 did reindexing)
        return pivoted

    def get_star_links(self):
        """
        Identify links connected to the 'star' node.
        Returns a list of tuples: (link_name, bandwidth, display_name)
        """
        topology = self.topology
        if not topology:
            return []

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
            return []
            
        links = []
        for child in star_node.get("children", []):
            uplink = child.get("uplink")
            bw = child.get("bandwidth", 0)
            name = child.get("name", child.get("id"))
            
            if uplink and bw > 0:
                links.append((uplink, bw, name))
        return links

