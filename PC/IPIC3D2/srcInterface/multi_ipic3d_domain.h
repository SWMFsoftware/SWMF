#ifndef MULTI_IPIC3D_DOMAIN_H
#define MULTI_IPIC3D_DOMAIN_H

#include "iPic3D.h"
#include <mpi.h>

using namespace iPic3D;

static int nDim = 0;
static int myrank = 0;
static int numProc = 0;
static MPI_Comm ipic3dComm = MPI_COMM_NULL;
 
const int nmaxIPIC = 100;
int nIPIC;
int iIPIC;
int *iSimCycle;

// first time for region
bool firstcall[nmaxIPIC];

// pointer to all IPIC3D 
iPic3D::c_Solver **SimRun;

// pointer to the stringstream object that contains the param.in
std::stringstream *param;

//store start time
double starttime[nmaxIPIC];
static double timenow;



// Store variables asosiated with comunicating the grid to and from GM
// When we recive data form BATSRUS we need to split it up and give it to
// the respective SimRun object. This to arrays will be used for this purpus.
int nGridPntSim[nmaxIPIC];  //number of grid pints *ndim used for each sim box
int nShiftGridPntSim[nmaxIPIC]; //fist index for each sims grid ponts

#endif
