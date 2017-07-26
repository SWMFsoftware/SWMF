#!/bin/csh

cd PGI/AMPS
module purge
module load comp-pgi/17.1
module load mpi-sgi/mpt
module load tecplot/2017r2
make $1 MPIRUN="mpiexec -n 8" >>& test_amps.log

cd ../..
echo Done > AmpsTestPGIComplete 
