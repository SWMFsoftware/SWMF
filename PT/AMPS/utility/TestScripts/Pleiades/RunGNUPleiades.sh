#!/bin/csh

cd GNU/AMPS
module purge
module load gcc/4.9.4
module load mpi-sgi
module load tecplot/2017r2
make $1 MPIRUN=mpiexec >>& test_amps.log

cd ../..
echo Done > AmpsTestGNUComplete 

