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


class cBitwiseFlagTable {
public:
  unsigned char *FlagTable;
  int FlagTableLength;

  //the default value of the FlagTable when initialized
  int DefaultFlagTableLength;
  double FlagTableLengthIncrement; //increment of the table length (persents) when re-allocated

  void AllocateTable(int size) {
    if (size<=FlagTableLength) return; //there is nothing to do

    //allocate the new table
    unsigned char *tmpTable=new unsigned char [size];
    int i;

    for (i=0;i<size;i++) tmpTable[i]=0;

    //copy the old table
    if (FlagTable!=NULL) {
      memcpy(tmpTable,FlagTable,FlagTableLength*sizeof(unsigned char));
      delete [] FlagTable;
    }

    //set up the new table
    FlagTable=tmpTable;
    FlagTableLength=size;
  }

  void SetDefaultParameters() {
    FlagTable=NULL, FlagTableLength=0;
    DefaultFlagTableLength=1000;
    FlagTableLengthIncrement=1.20; //increment of the table length (persents) when re-allocated
  }

  cBitwiseFlagTable() {
    SetDefaultParameters();
  }

  cBitwiseFlagTable(int size) {
    SetDefaultParameters();
    AllocateTable(1+(int)(size/8));
  }

  ~cBitwiseFlagTable() {
    if (FlagTable!=NULL) delete [] FlagTable;
  }

  bool Test(int GlobalBitePosition) {
    int iBit,iByte;
    unsigned char mask;

    iByte=GlobalBitePosition/8;
    iBit=GlobalBitePosition%8;

    mask=(unsigned char)(1<<iBit);
    return ((FlagTable[iByte]&mask)!=0) ? true : false;
  }

  void SetFlag(bool flag,int GlobalBitePosition) {
    int iBit,iByte;
    unsigned char mask;

    iByte=GlobalBitePosition/8;
    iBit=GlobalBitePosition%8;

    if (iByte>=FlagTableLength) {
      //re-allocate the table
      int NewFlagTableSize;

      NewFlagTableSize=(FlagTableLength==0) ? DefaultFlagTableLength : (int)(FlagTableLengthIncrement*FlagTableLength);
      AllocateTable(NewFlagTableSize);
    }

    mask=(unsigned char)(1<<iBit);

    if (flag==true) {
      FlagTable[iByte]|=mask;
    }
    else {
      FlagTable[iByte]^=mask;
    }
  }

  //collect flags from all processors (strided)
  void GatherStrided() {
    int rank,size,i,thread;
    int *LengthTable,MaxTableLength;

    //get the maximum size of the flag buffer
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&size);

    LengthTable=new int[size];
    MPI_Allgather(&FlagTableLength,1,MPI_INT,LengthTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

    for (MaxTableLength=FlagTableLength,i=0;i<size;i++) if (MaxTableLength<LengthTable[i]) MaxTableLength=LengthTable[i];

    delete [] LengthTable;

    //increase the length of the table if needed
    if (MaxTableLength>FlagTableLength) AllocateTable(MaxTableLength);

    //gather the flag table on the root processor
    if (rank==0) {
      MPI_Status status;
      unsigned char ExchangeBuffer[MaxTableLength];
      bool flag;

      for (int thread=1;thread<size; thread++) {
        MPI_Recv(ExchangeBuffer,MaxTableLength,MPI_UNSIGNED_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (i=0;i<MaxTableLength;i++) if (i%size==thread) {
          int iBit,iByte;
          unsigned char mask;

          iByte=i/8;
          iBit=i%8;

          mask=(unsigned char)(1<<iBit);
          flag=((ExchangeBuffer[iByte]&mask)!=0) ? true : false;

          SetFlag(flag,i);
        }
      }
    }
    else {
      MPI_Send(FlagTable,MaxTableLength,MPI_UNSIGNED_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
    }

    //broadcast the table
    MPI_Bcast(FlagTable,MaxTableLength,MPI_UNSIGNED_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  }
};



#endif /* _FLAGTABLE_H_ */
