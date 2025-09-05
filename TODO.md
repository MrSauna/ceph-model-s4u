# TODO

todo stuff 


## pg map approach

I need to decide wether I'll implement crush algorithm myself or do I just import pg map. 

Personally, I would like to implement crush algorithm just of curiosity. However, that is probably a big task. 
=> Really look into easier alternaties.
Is there a offline tool in ceph which yields a pg map based on a crush map.

What is osdmaptool capable of?
- --test-map-pgs-dump should be interesting
    - doesn't take into account reweights eg. osd out
    - thus crushdiff is incapable of taking that into account as well
- modify stuff
- can extract and embed the embedded crush map into the osdmap
    - a finding in itself - osdmap contains crush map

crushtool:
- For manipulating the crush map itsef. 
- Is capable of simulating mapping of input values to PGs

crushdiff:
- wrapper around osdmaptool's --test-map-pgs-dump
- needs an osdmap (includes crushmap) and a pg dump. By default takes them by connecting to cluster.
- explicitly tells me for each pg.id the {from set} -> {to set}

With these lessons, I can input only the PG map to the model. It includes the acting set and up set; or the current set vs the target set. This should be all the information, I need to accurately model client traffic to the acting set and backfill to the up/target set. 

Quick look at ceph's pgmap and pg class code tells me including the code would import half of ceph codebase. Seems excessive and quaranteed to be difficult. It seems that it should be simpler to implement my own parser and data structures for handling and updating the PGMap. 

=> In practice: modify osdmap so that 1. it has osdmap representing the initial state 2. modify its crush map so that it reflects the target state. 
=> crushdiff vs osdmaptool


## simulate rw ops

At first skip the scheduler completely. 


## simulate backfill

- Osd sends stuff to new location
- Receiving end is capable to receive it
- Keeps track of progress
    - How much of pg is sent
    - which pg is backfilling (limit parallel backfills)
- When ready, notifies mon
