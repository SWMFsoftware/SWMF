#!/bin/csh

cd Intel/AMPS

module purge
module load comp-intel/2016.2.181 
module load mpi-sgi
module load tecplot/2017r2
make $1 MPIRUN=mpiexec >>& test_amps.log

cd ../..
echo Done > AmpsTestIntelComplete 
