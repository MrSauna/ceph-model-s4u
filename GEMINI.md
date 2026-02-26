# Instructions

- This is a simulation of a Ceph cluster. The simulator focuses on modeling Ceph's backfill and client traffic. A lot of details of real Ceph are omitted.

- We are using simgrid as the simulation framework. The model code is written in c++ using the s4u api.

- The first significant simplification to real ceph, is that pgdumps from Ceph are used directly to get the data placement (pg placement). Objects are backfilled from the primary OSDs to up osds that previously didn't have the object. Objects are tracked so that the amount of disk read/write is simulated, network send/receive and per operation type possibly receiving end's disk write is simulated. Also when the 3-replicated writes are done each osd needs to ack the primary osd before ack is sent to a client. 

- The Osds are "Actors" in the sense that they share the same code but have different data. Each of them act individually and autonomously. 

- There is a monitor actor which has the pgdump shared to everybody via memory. Acceptable due to simulation framework technicalities.

- There is a Client actor which sends client requests to the cluster. 


## Building the c++ code

- I'm using cmake and directory `build` is pre-configured. It can be built using `cmake --build build`.


## Running the simulator

Here is a quick guide on how to run the simulator. Remember to wait at least 10 seconds for the simulation to finish.

```
./build/ceph-sim --cfg=network/model:CM02 --dc-shape=2:10:24 --dc-speed=100:100:25 --pgdump=results/big_cluster_no_backfill2/0/pgdump.txt --pgdump=results/big_cluster_no_backfill2/1/pgdump.txt --start-up-delay=0 --shut-down-delay=300 --profile=balanced --pg-objects=250 --object-size=4194304 --pool-id=1 --disk-write-bandwidth=268435456 --iops=80 --dc-clients=1 --client-op-size=8388608 --client-read-queue-depth=128 --client-write-queue-depth=32 --client-read-bandwidth=229780184 --client-write-bandwidth=70988595 --output-dir=results/big_cluster_no_backfill2
```

- Use --output-dir vscode to output the results to the vscode directory.

## Reading outputs (if necessary)

Most of the artifacts produced by the simulation are files. The files are in the `results/vscode` directory which is gitignored.

- `topology.json` describes the topology of the cluster. including network links and their properties.
- `mon_metrics.csv` describes the mon metrics.
- `net_metrics.csv` describes the network metrics.
- `client_metrics.csv` describes the client metrics.


## Running python scripts (mandatory prerequisites)

- You must run python scripts with `.venv/bin/python` otherwise the command will fail.
- The python scripts are in the `tools` directory. The `visualize.py` script can be run using `.venv/bin/python tools/visualize.py`. The output will go in `results/vscode/figures/`.
- If running other python commands do it in the `.venv` virtual environment.
