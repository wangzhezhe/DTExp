#!/bin/bash

# mpi is loaded
# dir is the build dir
# workflow server exist at current dir
# metaserver exist at current dir

rm *.log
rm -rf gs.bp*
rm -rf pdf.bp*
rm -rf iso.bp*
rm -rf Metaserver multinodeip

# start meta server for recording workfow running time
./metaserver & 

# start workflow server for event matching
rm -rf multinodeip          
export NUM_SERVER=2
export NUM_GROUPSIZE=2
export NUM_GROUP=1
export THREAD_POOL=4
export NOTIFY_SLEEP=1000
export COORDINATOR_THRESHOLD=1024
export DIR=.

./workflowserver 1500 eth0 $NUM_GROUPSIZE $NUM_GROUP $THREAD_POOL $NOTIFY_SLEEP $COORDINATOR_THRESHOLD 1 >& workflowserver.log &

sleep 0.5

# start sim
mpirun -n 1 ./gray-scott ./settings.json &> sim.log &

while [ ! -d ./gs.bp ]
do
    sleep 0.01
    echo "dir not exist"
done

# start pdfpushevent 
./pdfpushevent gs.bp &>pdfpushevent.log &