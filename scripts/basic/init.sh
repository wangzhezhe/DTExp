
#!/bin/bash


#module load openmpi/2.1.1

module use /projects/community/modulefiles/
module load gcc/5.4/openmpi/3.1.2-kholodvl

#set by env
#PJPATH_SCRATCH=/scratch/zw241/cworkspace/build/build_DTExp
#PJPATH_HOME=/home/zw241/cworkspace/src/DTExp

PJPATH_HOME=/home/zw241/cworkspace/src/DTExp

#cd $PJPATH_SCRATCH
cp $PJPATH_HOME/adios2.xml .
cp $PJPATH_HOME/simulation/settings.json .
cp ../anastartbyevent .

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

srun --mem-per-cpu=2000 --ntasks=1 --cpus-per-task=16 ./metaserver eno1 &> metaserver.log &

# check dir
while [ ! -d ./Metaserver ]
do
    sleep 0.01
done

echo "start meta server"


rm -rf multinodeip
cp /home/zw241/cworkspace/build/build_ewfs/workflowserver .

PARTITION=8