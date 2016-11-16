#!/bin/bash

# Environment configuration-----------------------
# The following commands (without '#') can be added to ~/.modules:
# module swap PrgEnv-cray PrgEnv-intel
# module load cray-hdf5-parallel
# module unload darshan
# module load idl

# Environment configuration-----------------------

# To run on XE node
#PBS -l nodes=1:ppn=32:xe

# To run on XK node
### PBS -l nodes=1:ppn=16:xk

#PBS -l walltime=48:00:00

# Specify queue priority: normal, high, low, debug
# Low priority job can potentially be killed after 30 minutes by higher
# priority jobs and restarted.
#PBS -q high

#PBS -N PostProc

# Send email if something happens
#PBS -m abe

cd $PBS_O_WORKDIR

echo "PBS_NUM_NODES = "$PBS_NUM_NODES
echo "PBS_NUM_PPN = "$PBS_NUM_PPN
echo "PBS_NP = "$PBS_NP

# This is needed to make the stacksize large
# This can/should also be put into the ~/.bashrc
ulimit -s unlimited

# Suppress some error messages of Blue Waters
export PMI_NO_FORK=1
export PMI_NO_PREINITIALIZE=1

aprun ./PostProc.pl -n=$PBS_NP -r=60 -s=7  >& PostProc.log

