#!/bin/bash
set -e # Fail if any command fails

# 1. Create a basic hierarchy
ceph osd crush add-bucket default root
ceph osd crush add-bucket datacenter1 datacenter
ceph osd crush move datacenter1 root=default

# 2. Add some racks and hosts
for rack in {1..2}; do
    ceph osd crush add-bucket rack$rack rack
    ceph osd crush move rack$rack datacenter=datacenter1
    
    for host in {1..4}; do
        host_name="r${rack}h${host}"
        ceph osd crush add-bucket $host_name host
        ceph osd crush move $host_name rack=rack$rack
        
        # 3. Add OSDs (Mocking them up)
        # We generate a UUID and an ID for them
        ceph osd create $uuid $id
        # Add them to the crush map (weight 1.0)
        # We need the last ID created, let's assume sequential for this mock
        osd_id=$(( ($rack - 1) * 4 + ($host - 1) ))
        ceph osd crush add osd.$osd_id 1.0 host=$host_name
    done
done

# 4. Create a pool so the map isn't empty
ceph osd pool create mypool 32