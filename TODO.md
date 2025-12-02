# TODO

todo stuff


## Import the pgdump.txt from osdmaptool --test-map-pgs-dump-all to the simulation

Do I parse the pgdump.txt fully in c++?

intermediate formats would be kinda dumb, parse directly in c++.


## How do I generate the platform?

Preferably it should be generated dynamically from the same parameters as as the crushmap; or from the crushmap.


## Create deployment

Defining it in c++ code would remove the requirement for deployment.xml. I can always use snakemake for defining different experiments using different cli arguments.


## add visualization of crushmaps to Snakemake

use the map-visualizer as base. Implement colors for unique DCs and make it part of the pipeline so it gets generated automatically