#include <time.h>
#include <math.h>
#include <stdlib.h>

#ifdef MPI_ON
  #include "mpi.h"
#endif

long int nint(double a)
{
   long int n;
   n=(long int) a;
   if (a-n>0.5) n++;
   return n;
}

void rnd_seed() {
  static char random_state_buffer[256];

#ifdef MPI_ON
  int thread;
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);
  initstate((unsigned) (thread+time(NULL)),random_state_buffer,256);
#else  
  initstate((unsigned) time(NULL),random_state_buffer,256);
#endif
}

//===================================================
double rnd() {
  return ((double)(random())+1.0)/2147483649.0;
}

//===================================================
double erf(double s) { 
  double b,d,c,t;

  b=fabs(s);
  if (b>4.0)
    d=1.0;
  else { 
    c=exp(-b*b);
    t=1.0/(1.0+0.3275911*b);
    d=1.0-(0.254829592*t-0.284496736*t*t+1.421413741*t*t*t-
      1.453152027*t*t*t*t+1.061405429*t*t*t*t*t)*c;
  }

  if (s<0.0) d=-d;
  return d;
}  

//===================================================
double gam(double x) {
  double a,y;

  a=1.0;
  y=x;

  if (y<1.0) 
    a/=y;
  else {
    y--;
    while (y>1.0) {
      a*=y;
      y--;
    }
  } 

  return a*(1.0-0.5748646*y+0.9512363*y*y-0.6998588*y*y*y+
    0.4245549*y*y*y*y-0.1010678*y*y*y*y*y); 
}
