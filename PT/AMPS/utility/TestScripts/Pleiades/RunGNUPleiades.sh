#!/bin/csh

cd GNU/AMPS
module purge
module load gcc/6.2 
module load mpi-sgi
module load tecplot/2017r2
make $1 "MPIRUN=mpiexec -n 8" >>& test_amps.log

cd ../..
echo Done > AmpsTestGNUComplete 

