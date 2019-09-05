

#!/bin/bash

# mpi is loaded
# dir is the build dir

rm *.log
rm -rf gs.bp*
rm -rf pdf.bp*
rm -rf iso.bp*
rm -rf Metaserver

# run metaserer

./metaserver &> metaserver.log &


while [ ! -d ./Metaserver ]
do
    sleep 0.01
    echo "Metaserver dir not exist"
done

echo "---ok to start metaserver"


# start original gray-scott
mpirun -n 1 ./gray-scott settings.json &>./gray-scott.log &

while [ ! -d ./gs.bp ]
do
    sleep 0.01
    echo "gs.bp dir not exist"
done

# start checking data
TASKNUM=32
#x=0
#while [ $x -lt $TASKNUM ]
#do
#  echo $x
  ./anapullmeta 1 $TASKNUM &>anapullmeta.log &
#  x=$(( $x + 1 ))
#done

# 2 4 8 16 32 64
# run pdfpushmeta
./pdfpushmeta gs.bp $TASKNUM &>pdfpushmeta.log &







