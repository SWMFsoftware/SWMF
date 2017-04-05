//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=================================================================
//$Id$
//=================================================================
//define the random generator

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef _RND_
#define _RND_

#ifndef _DO_NOT_LOAD_GLOBAL_H_
#include "global.h"
#endif 

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_

namespace RandomNumberGenerator {
  extern unsigned long int rndLastSeed;
  extern unsigned long int *rndLastSeedArray;

}

void rnd_seed(int seed=-1);

inline double rnd() {
  double res;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  RandomNumberGenerator::rndLastSeed*=48828125;
  RandomNumberGenerator::rndLastSeed&=2147483647; // pow(2,31) - 1
  if (RandomNumberGenerator::rndLastSeed==0) RandomNumberGenerator::rndLastSeed=1;
  res=double(RandomNumberGenerator::rndLastSeed/2147483648.0); //(pow(2,31) - 1) + 1

  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread=omp_get_thread_num();
  RandomNumberGenerator::rndLastSeedArray[thread]*=48828125;
  RandomNumberGenerator::rndLastSeedArray[thread]&=2147483647; // pow(2,31) - 1
  if (RandomNumberGenerator::rndLastSeedArray[thread]==0) RandomNumberGenerator::rndLastSeedArray[thread]=1;
  res=double(RandomNumberGenerator::rndLastSeedArray[thread]/2147483648.0); //(pow(2,31) - 1) + 1

  #else
  #error Unknown option
  #endif //_COMPILATION_MODE_

  return res;
}


#endif
