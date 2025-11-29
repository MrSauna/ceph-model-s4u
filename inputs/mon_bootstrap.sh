#!/bin/bash

# export local IP (or localhost)
export MON_IP=127.0.0.1

# 1. Generate a minimal config
cat <<EOF > /etc/ceph/ceph.conf
[global]
fsid = $(uuidgen)
mon_host = $MON_IP
auth_cluster_required = none
auth_service_required = none
auth_client_required = none
osd_pool_default_size = 1
EOF

# 2. Create a monmap for this single temporary monitor
monmaptool --create --add a $MON_IP:6789 --fsid $(grep fsid /etc/ceph/ceph.conf | awk '{print $3}') /tmp/monmap

# 3. Initialize the monitor data directory
mkdir -p /var/lib/ceph/mon/ceph-a
ceph-mon --mkfs -i a --monmap /tmp/monmap

# 4. Run the monitor in the background
ceph-mon -i a --public-addr $MON_IP:6789
