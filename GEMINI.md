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
./build/ceph-sim \
    --dc-shape=2:2:2 \
    --dc-speed=10:10:10 \
    --dc-weight=10:: \
    --pgdump=results/small_cluster/0/pgdump.txt \
    --pgdump=results/small_cluster/1/pgdump.txt \
    --output-dir=results/vscode \
    --pg-objects=1000 \
    --object-size=4194304 \
    --pool-id=1 \
    --disk-read-bandwidth=125829120 \
    --disk-write-bandwidth=83886080 \
    --clients=1 \
    --client-read-queue-depth=4 \
    --client-write-queue-depth=4
```

- The --dc-shape is a toplogoy descriptor, you don't need to modify unless specifically asked to. The --dc-speed is the speed of the links, you don't need to modify unless specifically asked to. The --dc-weight is the weight of the links, you don't need to modify unless specifically asked to.

- The --pgdump defines the data that needs to be moved in the simulation but you can not create a pgdump file yourself. 

- The --pg-objects=1000 says that 1000 objects exist per pg. The --object-size=4194304 says that each object is 4MB. This practically means the whole simulation that should run in about 5 seconds wall time.

- Use --output-dir vscode to output the results to the vscode directory.

## Reading outputs (if necessary)

Most of the artifacts produced by the simulation are files. The files are in the `results/vscode` directory which is gitignored.

- `topology.json` describes the topology of the cluster. including network links and their properties.
- `mon_metrics.csv` describes the mon metrics.
- `net_metrics.csv` describes the network metrics.
- `client_metrics.csv` describes the client metrics.


## Running python scripts

- The python scripts are in the `tools` directory. The `visualize.py` script can be run using `.venv/bin/python tools/visualize.py`. The output will go in `results/vscode/figures/`.
- If running other python commands do it in the `.venv` virtual environment.

