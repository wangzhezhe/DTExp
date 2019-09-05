
#!/bin/sh

module load openmpi/2.1.1
PJPATH_SCRATCH=/scratch/zw241/cworkspace/build/build_DTExp
PJPATH_HOME=/home/zw241/cworkspace/src/DTExp


cd $PJPATH_SCRATCH
cp $PJPATH_HOME/adios2.xml .
cp $PJPATH_HOME/simulation/settings.json .

# clean the previous adios file
rm *.log
rm -rf gs.bp*
rm -rf pdf.bp*
rm -rf iso.bp*

# check vtkdata
if [ ! -d ./vtkdata ]; then
  mkdir ./vtkdata;
fi


# start metaserver for testing
rm -rf Metaserver
cp /home/zw241/cworkspace/build/build_MMServer/metaserver .
sleep 1
./metaserver &> metaserver.log &

# check dir
while [ ! -d ./Metaserver ]
do
    sleep 0.01
    echo "Metaserver dir not exist"
done

cp /home/zw241/cworkspace/build/build_ewfs/workflowserver .
rm -rf multinodeip