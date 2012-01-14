//=================================================================
//$Id$
//=================================================================
//define the random generator

#ifndef _RND_
#define _RND_

namespace RandomNumberGenerator {
  extern int rndLastSeed;
}

void rnd_seed(int seed=-1);

inline double rnd() {
  RandomNumberGenerator::rndLastSeed*=48828125;
  if (RandomNumberGenerator::rndLastSeed<0) RandomNumberGenerator::rndLastSeed=(RandomNumberGenerator::rndLastSeed+2147483647)+1;
  if (RandomNumberGenerator::rndLastSeed==0) RandomNumberGenerator::rndLastSeed=1;

  return RandomNumberGenerator::rndLastSeed/2147483647.0;
}


#endif
