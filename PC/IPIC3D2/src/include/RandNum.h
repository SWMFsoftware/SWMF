#ifndef RANDNUM_H
#define RANDNUM_H

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double PseudoRand(long *idum){

  // "minimal" random number generator of Park and Miller. Taken from
  // "numerical recipes in C"
  // Return a uniform random deviate between 0.0 and 1.0. Set or reset idum 
  // to any integer value (except the unlikely value MASK) it initialize 
  // the sequence; idum must not be alterd between calls or successive 
  // deviates in a sequence.
  
  long k;
  double ans;
  
  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum <0) *idum +=IM;
  ans=AM*(*idum);
  *idum ^=MASK;
  return ans;
  
};

// cleaning
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

#endif
