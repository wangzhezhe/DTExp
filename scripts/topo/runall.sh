#!/bin/bash

rm -rf topo_*

while [ ! -d ./topo_pr ]
do
    mkdir topo_pr
done
cd topo_pr
sbatch ~/cworkspace/src/DTExp/scripts/topo/pr.scripts
cd ..



while [ ! -d ./topo_cr ]
do
    mkdir topo_cr
done
cd topo_cr
sbatch ~/cworkspace/src/DTExp/scripts/topo/cr.scripts
cd ..


while [ ! -d ./topo_mr ]
do
    mkdir topo_mr
done
cd topo_mr
sbatch ~/cworkspace/src/DTExp/scripts/topo/mr.scripts
cd ..


while [ ! -d ./topo_nr ]
do
    mkdir topo_nr
done
cd topo_nr
sbatch ~/cworkspace/src/DTExp/scripts/topo/nr.scripts
cd ..