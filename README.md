# tasks

The simulation code comes from the https://github.com/pnorbert/adiosvm

# scripts for testing

Explain four patterns in detail

# other middleware used in this experiments

**The metadata service**

**The event matching service**


**temp**

module use /projects/community/modulefiles/
module load gcc/5.4/openmpi/3.1.2-kholodvl


on amarel
cmake ~/cworkspace/src/DTExp/ -DADIOS2_DIR=~/cworkspace/build/build_ADIOS2/ -DVTK_DIR=~/cworkspace/build/build_vtk -DUSE_TIMERS=ON -DVTK=ON