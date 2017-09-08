#!/bin/csh

cd Intel/AMPS

module purge
module load comp-intel/2016.2.181 
module load mpi-sgi
module load tecplot/2017r2
module load boost/1.62
make $1 MPIRUN="mpiexec -n 8 dplace -s1 -c 0-7 " TESTMPIRUN1="mpiexec -n 1 omplace " >>& test_amps.log

cd ../..
echo Done > AmpsTestIntelComplete 
