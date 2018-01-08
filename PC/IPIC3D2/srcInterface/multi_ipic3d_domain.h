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
extern int nIPIC;
extern int iIPIC; 
extern int *iSimCycle;

// pointer to all IPIC3D 
extern iPic3D::c_Solver **SimRun;

// pointer to the stringstream object that contains the param.in
extern std::stringstream *param;

#endif
