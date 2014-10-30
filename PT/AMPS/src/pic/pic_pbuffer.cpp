//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//==========================================================
//$Id$
//==========================================================
//the particle data buffer


#include "pic.h"

long int PIC::ParticleBuffer::ParticleDataLength=_PIC_PARTICLE_DATA__DEFAULT_DATA_LENGTH_;
PIC::ParticleBuffer::byte *PIC::ParticleBuffer::ParticleDataBuffer=NULL;
long int PIC::ParticleBuffer::MaxNPart=0;
long int PIC::ParticleBuffer::NAllPart=0;
long int PIC::ParticleBuffer::FirstPBufferParticle=-1;

//==========================================================
//init the buffer
void PIC::ParticleBuffer::Init(long int BufrerLength) {

  if ((ParticleDataBuffer!=NULL)||(MaxNPart!=0)) exit(__LINE__,__FILE__,"Reallocation of the particle data buffer");
  if (sizeof(byte)!=1) exit(__LINE__,__FILE__,"The size of 'byte' is diferent from 1");
  if (BufrerLength<=0) exit(__LINE__,__FILE__,"BufrerLength is less that zero");

  //reserve the space for additional 'particle's variables'

  //allocate the memory for the buffer
  MaxNPart=BufrerLength;
  ParticleDataBuffer=(PIC::ParticleBuffer::byte*) malloc(ParticleDataLength*MaxNPart);

  //init the list of particles in the buffer
  for (long int ptr=0;ptr<MaxNPart-1;ptr++) SetNext(ptr+1,ptr);
  SetNext(-1,MaxNPart-1);
  FirstPBufferParticle=0;

}

//==========================================================
//Request additional data for a particle
void PIC::ParticleBuffer::RequestDataStorage(long int &offset,int TotalDataLength) {
  if (ParticleDataBuffer!=NULL) exit(__LINE__,__FILE__,"Error: the particle data buffer is already initialized. Request the particle data storage before the initialization of the particle data buffer");

  offset=ParticleDataLength;
  ParticleDataLength+=TotalDataLength;
}

//==========================================================
//the basic data access functions for a particle
PIC::ParticleBuffer::byte *PIC::ParticleBuffer::GetParticleDataPointer(long int ptr) {
  return ParticleDataBuffer+ptr*ParticleDataLength;
}


//==========================================================
//the functions that controls the particle buffer
long int PIC::ParticleBuffer::GetMaxNPart() {return MaxNPart;}
long int PIC::ParticleBuffer::GetAllPartNum() {return NAllPart;}
long int PIC::ParticleBuffer::GetParticleDataLength() {return ParticleDataLength;}

long int PIC::ParticleBuffer::GetNewParticle() {
  long int newptr;
  byte *pdataptr;

  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(pdataptr);
  #endif

  FirstPBufferParticle=GetNext(pdataptr);
  SetPrev(-1,pdataptr);
  SetNext(-1,pdataptr);

//#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (IsParticleAllocated(pdataptr)==true) exit(__LINE__,__FILE__,"Error: the particle is re-allocated");
  SetParticleAllocated(pdataptr);
//#endif

  return newptr;
}

long int PIC::ParticleBuffer::GetNewParticle(long int &ListFirstParticle) {
  long int newptr;
  byte *pdataptr;

  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(pdataptr);
  #endif

  FirstPBufferParticle=GetNext(pdataptr);

  if (ListFirstParticle>=0) {
    byte *listFirstPData=GetParticleDataPointer(ListFirstParticle);

    SetPrev(GetPrev(listFirstPData),pdataptr);
    SetNext(ListFirstParticle,pdataptr);
    SetPrev(newptr,listFirstPData);
  }
  else {
    SetPrev(ListFirstParticle,pdataptr);
    SetNext(ListFirstParticle,pdataptr);
  }

//#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (IsParticleAllocated(pdataptr)==true) exit(__LINE__,__FILE__,"Error: the particle is re-allocated");
  SetParticleAllocated(pdataptr);
//#endif

  ListFirstParticle=newptr;
  return newptr;
}

void PIC::ParticleBuffer::ExcludeParticleFromList(long int ptr,long int& ListFirstParticle) {
  byte *pdataptr=GetParticleDataPointer(ptr);
  long int prev,next;

  //exclude the particle from the list
  prev=GetPrev(pdataptr);
  next=GetNext(pdataptr);

  if (ptr==ListFirstParticle) {
    SetPrev(prev,next);
    ListFirstParticle=next;
  }
  else {
    SetNext(next,prev);
    SetPrev(prev,next);
  }
}


void PIC::ParticleBuffer::DeleteParticle(long int ptr) {
  //terminate the particle trajectory sampling
  #if _PIC_PARTICLE_TRACKER_MODE_  == _PIC_MODE_ON_
  byte *ParticleData=GetParticleDataPointer(ptr);
  PIC::ParticleTracker::FinilazeParticleRecord(ParticleData);
  #endif

  DeleteParticle_withoutTrajectoryTermination(ptr);
}

void PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(long int ptr) {

//#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (IsParticleAllocated(ptr)==false) exit(__LINE__,__FILE__,"Error: the particle is re-deleted");
  SetParticleDeleted(ptr);
//#endif

  NAllPart--;
  SetNext(FirstPBufferParticle,ptr);
  FirstPBufferParticle=ptr;
}


void PIC::ParticleBuffer::DeleteParticle(long int ptr,long int& ListFirstParticle) {
  ExcludeParticleFromList(ptr,ListFirstParticle);
  DeleteParticle(ptr);
}



void PIC::ParticleBuffer::CloneParticle(long int copy,long int source) {
  byte *SourceData,*CopyData;
  long int next,prev;

  SourceData=GetParticleDataPointer(source);
  CopyData=GetParticleDataPointer(copy);

  prev=GetPrev(CopyData);
  next=GetNext(CopyData);

//  for (i=0;i<ParticleDataLength;i++) CopyData[i]=SourceData[i];
  memcpy(CopyData,SourceData,ParticleDataLength*sizeof(byte));

  SetPrev(prev,CopyData);
  SetNext(next,CopyData);
}

//==========================================================
//save the particle buffer in a restart file
void PIC::ParticleBuffer::SaveImageFile(int fd) {
  exit(__LINE__,__FILE__,"Not implemented");
}

void PIC::ParticleBuffer::LoadImageFile(int fd) {
  exit(__LINE__,__FILE__,"not implemented");
}


//==========================================================
//pack the particle data

void PIC::ParticleBuffer::PackParticleData(char* buffer,long int ptr) {
  byte *SourceData=GetParticleDataPointer(ptr);
//  long int i;

//  for (int i=0;i<ParticleDataLength;i++) buffer[i]=SourceData[i];

  memcpy(buffer,SourceData,ParticleDataLength*sizeof(byte));
}


void PIC::ParticleBuffer::UnPackParticleData(char* buffer,long int ptr) {
  byte *pdata;
  long int next,prev;

  pdata=GetParticleDataPointer(ptr);
  prev=GetPrev(pdata);
  next=GetNext(pdata);

//  for (int i=0;i<ParticleDataLength;i++) pdata[i]=buffer[i];

  memcpy(pdata,buffer,ParticleDataLength*sizeof(byte));

  SetPrev(prev,pdata);
  SetNext(next,pdata);
}

//==========================================================
//get the checksum of the particle buffer
unsigned long PIC::ParticleBuffer::GetChecksum() {
  CRC32 sum;

  //save the particle's buffer internal data
  sum.add(&ParticleDataLength,1);
  sum.add(&MaxNPart,1);
  sum.add(&NAllPart,1);
  sum.add(&FirstPBufferParticle,1);

  //save the particle's data
  sum.add(ParticleDataBuffer,MaxNPart*ParticleDataLength);

  unsigned long int *buffer=new unsigned long int[TotalThreadsNumber];
  char str[10*_MAX_STRING_LENGTH_PIC_];

  buffer[0]=sum.checksum();

  unsigned long int bufferRecv[TotalThreadsNumber];
  MPI_Gather(buffer,1,MPI_UNSIGNED_LONG,bufferRecv,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(buffer,bufferRecv,TotalThreadsNumber*sizeof(unsigned long int));

  if (ThisThread==0) {
    sprintf(str,"Cdsmc::pbuffer CRC32 checksum: ");

    for (long int thread=0;thread<TotalThreadsNumber;thread++) sprintf(str,"%s 0x%lx ",str,buffer[thread]);

    printf("$PREFIX:%s\n",str);
    PrintErrorLog(str);
  }

  delete [] buffer;
  return sum.checksum();
}

//==========================================================
//check particle list -> calculate the number of particles stored in the lists and compare with the total number of particles stored in the particle buffer
void PIC::ParticleBuffer::CheckParticleList() {
  long int nTotalListParticles=0;
  int i,j,k; //,LocalCellNumber;
  long int ParticleList;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
//  PIC::Mesh::cDataCenterNode *cell;


  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

//  int centerNodeIndexCounter;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  PIC::Mesh::cDataBlockAMR block;

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));


  //the tables of the first particles in the cells
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],tempParticleMovingListTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];


  for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
    node=(thread==PIC::Mesh::mesh.ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    if (node==NULL) continue;

    //sample the processor load
  #if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    double EndTime,StartTime=MPI_Wtime();
  #endif


    while (node!=NULL) {

      memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
      memcpy(tempParticleMovingListTable,node->block->tempParticleMovingListTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));


      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              /*LocalCellNumber=*/   PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

/*
      memcpy(&block,node->block,sizeof(PIC::Mesh::cDataBlockAMR));

      {
        {
          for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
              LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
              cell=block.GetCenterNode(LocalCellNumber);

//              cell=node->block->GetCenterNode(LocalCellNumber);
 */

//              if (cell!=NULL) {
                if (tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]!=-1) exit(__LINE__,__FILE__,"Error: the temp list is not empty");

                ParticleList=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

                while (ParticleList!=-1) {
                  ++nTotalListParticles;
                  ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                  if (nTotalListParticles>NAllPart) exit(__LINE__,__FILE__,"Error: the list particles' number exeeds the total number of particles stored in the buffer");
                }
//              }
            }

            if (DIM==1) break;
         }

         if ((DIM==1)||(DIM==2)) break;
      }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

      node=node->nextNodeThisThread;
    }
  }

  if (nTotalListParticles!=NAllPart) exit(__LINE__,__FILE__,"Error: the total number of particles stored in the lists is different from that stored in the particle buffer");
}
