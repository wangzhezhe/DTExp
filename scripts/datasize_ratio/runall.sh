#!/bin/bash

rm -rf datasize_ratio_*

while [ ! -d ./datasize_ratio_pr ]
do
    mkdir datasize_ratio_pr
done
cd datasize_ratio_pr
sbatch ~/cworkspace/src/DTExp/scripts/datasize_ratio/pr.scripts
cd ..



while [ ! -d ./datasize_ratio_cr ]
do
    mkdir datasize_ratio_cr
done
cd datasize_ratio_cr
sbatch ~/cworkspace/src/DTExp/scripts/datasize_ratio/cr.scripts
cd ..


while [ ! -d ./datasize_ratio_mr ]
do
    mkdir datasize_ratio_mr
done
cd datasize_ratio_mr
sbatch ~/cworkspace/src/DTExp/scripts/datasize_ratio/mr.scripts
cd ..


while [ ! -d ./datasize_ratio_nr ]
do
    mkdir datasize_ratio_nr
done
cd datasize_ratio_nr
sbatch ~/cworkspace/src/DTExp/scripts/datasize_ratio/nr.scripts
cd ..