//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=================================================================
//$Id$
//=================================================================
//define the random generator

#include "mpi.h"

#include "rnd.h"

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_

unsigned long int RandomNumberGenerator::rndLastSeed=0;
unsigned long int *RandomNumberGenerator::rndLastSeedArray=NULL;

void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&thread);

  if (seed==-1) seed=thread;

  RandomNumberGenerator::rndLastSeed=seed;

  //init the seed array in case OpenMP is used
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    int nThreadsOpenMP,i;

    #pragma omp parallel private(i,nThreadsOpenMP)
    {
      #pragma omp single
      {
        nThreadsOpenMP=omp_get_num_threads();

        if (RandomNumberGenerator::rndLastSeedArray==NULL) RandomNumberGenerator::rndLastSeedArray=new unsigned long int[nThreadsOpenMP];
        for (i=0;i<nThreadsOpenMP;i++) RandomNumberGenerator::rndLastSeedArray[i]=abs(seed)+i+thread*nThreadsOpenMP;
      }
    }
  #endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

}
