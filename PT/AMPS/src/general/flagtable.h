//$Id$
//the class contaning the bitwise flag table

//returned values: FlagTable==0 -> false; FlagTable!=0 -> true
//the default value of the table is 'false'


/*
 * flagtable.h
 *
 *  Created on: Apr 12, 2017
 *      Author: vtenishe
 */

#ifndef _FLAGTABLE_H_
#define _FLAGTABLE_H_

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "specfunc.h"


class cBitwiseFlagTable {
public:
  unsigned char **FlagTable;
  int *FlagTableLength;
  int nThreadsOpenMP;


  //the default value of the FlagTable when initialized
  int DefaultFlagTableLength;
  double FlagTableLengthIncrement; //increment of the table length (persents) when re-allocated

  void AllocateTable(int size,int iThreadOpenMP) {
//    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
//    #pragma omp critical
//    {
//    #endif

      if (size>10000000) printf("asdasdlkahjds\n");
      if (iThreadOpenMP>=16) printf("zXCMfgaksjhdfb\n");

      if (FlagTableLength==NULL) printf("xkljghfdslkghjs\n");
      if (FlagTable==NULL)  printf("xzxvdhbkjhv\n");

    if (size>FlagTableLength[iThreadOpenMP]) {
      //allocate the new table
      unsigned char *tmpTable=new unsigned char [size];
      int i;

      for (i=0;i<size;i++) tmpTable[i]=0;

      //copy the old table
      if (FlagTable[iThreadOpenMP]!=NULL) {
        memcpy(tmpTable,FlagTable[iThreadOpenMP],FlagTableLength[iThreadOpenMP]*sizeof(unsigned char));
        delete [] FlagTable[iThreadOpenMP];
      }

      //set up the new table
      FlagTable[iThreadOpenMP]=tmpTable;
      FlagTableLength[iThreadOpenMP]=size;
    }

//    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
//    }
//    #endif
  }

  void SetDefaultParameters() {
    //allocate the tables
    nThreadsOpenMP=1;

   #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
   #pragma omp parallel default(none)
   {
     #pragma omp single
     {
     nThreadsOpenMP=omp_get_num_threads();
     }
   }
   #endif

    FlagTable=new unsigned char *[nThreadsOpenMP];
    FlagTableLength=new int [nThreadsOpenMP];

    for (int thread=0;thread<nThreadsOpenMP;thread++) FlagTable[thread]=NULL,FlagTableLength[thread]=0;

    //default values of the memory allocation parameters
    DefaultFlagTableLength=1000;
    FlagTableLengthIncrement=1.20; //increment of the table length (persents) when re-allocated
  }

  cBitwiseFlagTable() {
    SetDefaultParameters();
  }

  cBitwiseFlagTable(int size) {
    SetDefaultParameters();
    for (int thread=0;thread<nThreadsOpenMP;thread++) AllocateTable(1+(int)(size/8),thread);
  }

  ~cBitwiseFlagTable() {
    if (FlagTable!=NULL) {
      for (int thread=0;thread<nThreadsOpenMP;thread++) if (FlagTable[thread]!=NULL) delete [] FlagTable[thread];

      delete [] FlagTable;
      delete [] FlagTableLength;

      FlagTable=NULL;
    }
  }

  bool Test(int GlobalBitePosition) {
    int iBit,iByte;
    unsigned char mask;

    iByte=GlobalBitePosition/8;
    iBit=GlobalBitePosition%8;

    if (iByte>=FlagTableLength[0]) exit(__LINE__,__FILE__,"Error: iByte exeeds the size of the data bufer");

    mask=(unsigned char)(1<<iBit);
    return ((FlagTable[0][iByte]&mask)!=0) ? true : false;
  }

  void SetFlag(bool flag,int GlobalBitePosition) {
    int iBit,iByte,iThreadOpenMP=0;
    unsigned char mask;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    iThreadOpenMP=omp_get_thread_num();
    #endif

    iByte=GlobalBitePosition/8;
    iBit=GlobalBitePosition%8;

//    printf("%i %i %i %i %i \n",iThreadOpenMP,GlobalBitePosition,iByte,iBit,FlagTableLength[iThreadOpenMP]);

    if (iByte>=FlagTableLength[iThreadOpenMP]) {
      //re-allocate the table
      int NewFlagTableSize;

      NewFlagTableSize=(FlagTableLength[iThreadOpenMP]==0) ? DefaultFlagTableLength : (int)(FlagTableLengthIncrement*FlagTableLength[iThreadOpenMP]);
      if (NewFlagTableSize<=iByte) NewFlagTableSize=(int)(FlagTableLengthIncrement*iByte);


//      printf("!!!!!!  1: %i %i %i %i\n",iThreadOpenMP,iByte,iBit,FlagTableLength[iThreadOpenMP]);
      AllocateTable(NewFlagTableSize,iThreadOpenMP);
//      printf("!!!!!!  2: %i %i %i %i\n",iThreadOpenMP,iByte,iBit,FlagTableLength[iThreadOpenMP]);
    }



    if (flag==true) {
      mask=(unsigned char)(1<<iBit);
      FlagTable[iThreadOpenMP][iByte]|=mask;
    }
    else {
      mask=(unsigned char)(0xff^(1<<iBit));
      FlagTable[iThreadOpenMP][iByte]&=mask;
    }
  }

  //collect flags from all processors (strided)
  void Gather() {
    int rank,size,i,thread;
    int *LengthTable,MaxTableLength=0;

    //get the maximum size of the flag buffer
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&size);

    LengthTable=new int[size];

    MaxTableLength=FlagTableLength[0];
    for (thread=0;thread<nThreadsOpenMP;thread++) if (MaxTableLength<FlagTableLength[thread]) MaxTableLength=FlagTableLength[thread];

    MPI_Allgather(&MaxTableLength,1,MPI_INT,LengthTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

    for (i=0;i<size;i++) if (MaxTableLength<LengthTable[i]) MaxTableLength=LengthTable[i];

    delete [] LengthTable;

    //increase the length of the table if needed
    for (thread=0;thread<nThreadsOpenMP;thread++) if (MaxTableLength>FlagTableLength[thread]) AllocateTable(MaxTableLength,thread);

    //combine flags from all OpenMP threads
    for (thread=1;thread<nThreadsOpenMP;thread++) for (i=0;i<FlagTableLength[0];i++) FlagTable[0][i]|=FlagTable[thread][i];

    //gather the flag table on the root processor
    if (rank==0) {
      MPI_Status status;
      unsigned char ExchangeBuffer[MaxTableLength];
      bool flag;

      for (int thread=1;thread<size; thread++) {
        MPI_Recv(ExchangeBuffer,MaxTableLength,MPI_UNSIGNED_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (i=0;i<MaxTableLength;i++) for (i=0;i<FlagTableLength[0];i++) FlagTable[0][i]|=ExchangeBuffer[i];
      }
    }
    else {
      MPI_Send(FlagTable[0],MaxTableLength,MPI_UNSIGNED_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
    }

    //broadcast the table
    MPI_Bcast(FlagTable[0],MaxTableLength,MPI_UNSIGNED_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  }

};



#endif /* _FLAGTABLE_H_ */
