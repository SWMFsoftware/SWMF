

#ifdef MPI_ON
#include "mpi.h"
#endif

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern int ThisThread;
extern int TotalThreadsNumber;

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
  initstate((unsigned) (thread /*+time(NULL)*/ ),random_state_buffer,256);
#else  
  initstate((unsigned) 0 /*time(NULL)*/ ,random_state_buffer,256);
#endif
}

//===================================================
double rnd() {
  return ((double)(random())+1.0)/2147483649.0;
}

//===================================================
/*
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
*/

//===================================================
/*
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
*/
//===================================================
void PrintErrorLog(const char* message) {
  FILE* errorlog=fopen("error.log","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"%s\n\n",message);

  fclose(errorlog);
}

//use: PrintErrorLog(__LINE__,__FILE__, "mesage")
void PrintErrorLog(long int nline, const char* fname, const char* message) {

return;

  FILE* errorlog=fopen("error.log","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"file=%s, line=%ld\n",fname,nline);
  fprintf(errorlog,"%s\n\n",message);

  fclose(errorlog);
}


//===================================================
void StampSignature(char* message) {
  double *buffer=new double[TotalThreadsNumber];
  double sign=0.0;
  int thread;

  buffer[0]=rnd();

#ifdef MPI_ON
  MPI_Gather(buffer,1,MPI_DOUBLE,buffer,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

  if (ThisThread==0) {
    for (thread=0;thread<TotalThreadsNumber;thread++) sign+=buffer[thread];
     
    printf("Signature=%e (msg: %s)\n",sign,message); 
  }

  delete [] buffer;
}

//===================================================
//use: exit(__LINE__,__FILE__, "mesage")
void exit(long int nline, const char* fname, const char* msg) { 
  char str[1000];

  if (msg==NULL) sprintf(str," exit: line=%ld, file=%s\n",nline,fname);
  else sprintf(str," exit: line=%ld, file=%s, message=%s\n",nline,fname,msg); 

  PrintErrorLog(str);
  exit(0);
}

void PrintLineMark(long int nline ,char* fname ,char* msg) {
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (ThisThread==0) {
    if (msg==NULL) printf("linemark: line=%ld, file=%s\n",nline,fname);
    else printf("linemark: line=%ld, file=%s, message=%s\n",nline,fname,msg);
  }
}






