//=================================================================
//$Id$
//=================================================================
//define the random generator

#include "mpi.h"

#include "rnd.h"


unsigned long int RandomNumberGenerator::rndLastSeed=0;

void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&thread);

  if (seed==-1) seed=thread;

  RandomNumberGenerator::rndLastSeed=seed;
}
