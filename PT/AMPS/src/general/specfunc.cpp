//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf


#ifdef MPI_ON
#include "mpi.h"
#endif

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "rnd.h"
#include "specfunc.h"


int ThisThread;
int TotalThreadsNumber;
int ExitErrorCode=0;

long int nint(double a)
{
   long int n;
   n=(long int) a;
   if (a-n>0.5) n++;
   return n;
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
*/
//===================================================
void PrintErrorLog(const char* message) {
  FILE* errorlog=fopen("$ERRORLOG","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"%s\n\n",message);

  fclose(errorlog);
}

//use: PrintErrorLog(__LINE__,__FILE__, "mesage")
void PrintErrorLog(long int nline, const char* fname, const char* message) {
  FILE* errorlog=fopen("$ERRORLOG","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"file=%s, line=%ld\n",fname,nline);
  fprintf(errorlog,"%s\n\n",message);

#if _STDOUT_ERRORLOG_MODE_ == _STDOUT_ERRORLOG_MODE__ON_
  printf("$PREFIX:Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  printf("$PREFIX:file=%s, line=%ld\n",fname,nline);
  printf("$PREFIX:%s\n\n",message);
#endif

  fclose(errorlog);
}


//===================================================
void StampSignature(char* message) {
  double *buffer=new double[TotalThreadsNumber];
  double sign=0.0;
  int thread;

  buffer[0]=rnd();

#ifdef MPI_ON
  double bufferRecv[TotalThreadsNumber];

  MPI_Gather(buffer,1,MPI_DOUBLE,bufferRecv,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(buffer,bufferRecv,TotalThreadsNumber*sizeof(double));
#endif

  if (ThisThread==0) {
    for (thread=0;thread<TotalThreadsNumber;thread++) sign+=buffer[thread];
     
    printf("$PREFIX:Signature=%e (msg: %s)\n",sign,message);
  }

  delete [] buffer;
}

//===================================================
//use: exit(__LINE__,__FILE__, "mesage")
void exit(long int nline, const char* fname, const char* msg) {
  char str[1000];

  if (msg==NULL) sprintf(str," exit: line=%ld, file=%s (error code=%i)\n",nline,fname,ExitErrorCode);
  else sprintf(str," exit: line=%ld, file=%s, message=%s (error code=%i)\n",nline,fname,msg,ExitErrorCode);

  printf("$PREFIX:%s",str);

  PrintErrorLog(str);
  exit(0);
}

void PrintLineMark(long int nline ,char* fname ,char* msg) {
#ifdef MPI_ON
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
#endif

  if (ThisThread==0) {
    if (msg==NULL) printf("$PREFIX:linemark: line=%ld, file=%s\n",nline,fname);
    else printf("$PREFIX:linemark: line=%ld, file=%s, message=%s\n",nline,fname,msg);
  }
}

//=============================================================
//print a message into the debugger stream
void Debugger::SaveDataIntoStream(void* data,int length,const char* msg) {

  struct cStreamBuffer {
    int CollCounter;
    unsigned long CheckSum;
    char CallPoint[200];
  };

  const int CheckSumBufferLength=500;
  static cStreamBuffer StreamBuffer[CheckSumBufferLength];

  static int BufferPointer=0;
  static int CallCounter=0;
  static CRC32 CheckSum;

  //create a new copy of the dbuffer stream file at the first call of this function
  if (CallCounter==0) {
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",ThisThread);
    fout=fopen(fn,"w");
    fclose(fout);
  }

  //increment the call counter
  CallCounter++;

  //update the check sum
  CheckSum.add((char*)data,length);

  //save the checksum in the buffer
  StreamBuffer[BufferPointer].CollCounter=CallCounter;
  StreamBuffer[BufferPointer].CheckSum=CheckSum.checksum();
  sprintf(StreamBuffer[BufferPointer].CallPoint,"%s",msg);
  BufferPointer++;

  if (BufferPointer>=CheckSumBufferLength) {
    //save the accumulated checksums into a file
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",ThisThread);
    fout=fopen(fn,"a");

    for (int i=0;i<BufferPointer;i++) fprintf(fout,"%i: 0x%lx\t%s\n",StreamBuffer[i].CollCounter,StreamBuffer[i].CheckSum,StreamBuffer[i].CallPoint);

    BufferPointer=0;
    fclose(fout);
  }
}

template <class T>
void Debugger::SaveDataIntoStream(T data,const char* msg) {
  Debugger::SaveDataIntoStream(&data,sizeof(T),msg);
}



