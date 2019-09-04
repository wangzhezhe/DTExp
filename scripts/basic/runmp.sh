

#!/bin/bash

# mpi is loaded
# dir is the build dir

rm *.log
rm -rf gs.bp*
rm -rf pdf.bp*
rm -rf iso.bp*
rm -rf Metaserver

# run metaserer

./metaserver & 


while [ ! -d ./Metaserver ]
do
    sleep 0.01
    echo "Metaserver dir not exist"
done

echo "---ok to start metaserver"


# run sim
mpirun -n 1 ./gray-scott ./settings.json &> sim.log &

echo "ok to start sim"

# check bp file
# use BP4, there is only dir
while [ ! -d ./gs.bp ]
do
    sleep 0.01
    echo "dir not exist"
done


if [ ! -d ./vtkdata ]; then
  mkdir ./vtkdata;
fi

./anapullmeta 4 &> anapullmeta.log & 
#echo "ok to start anapullmeta"

# run pdfpushmeta
./pdfpushmeta gs.bp &>pdfpushmeta.log &

echo "ok to start pdfpushmeta"







