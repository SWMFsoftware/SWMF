//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//==========================================================
//$Id$
//==========================================================
//the particle data buffer


#include "pic.h"

long int PIC::ParticleBuffer::ParticleDataLength=_PIC_PARTICLE_DATA__FULL_DATA_LENGTH_;
PIC::ParticleBuffer::byte *PIC::ParticleBuffer::ParticleDataBuffer=NULL;
long int PIC::ParticleBuffer::MaxNPart=0;
long int PIC::ParticleBuffer::NAllPart=0;
long int PIC::ParticleBuffer::FirstPBufferParticle=-1;

int PIC::ParticleBuffer::Thread::NTotalThreads=0;
long int *PIC::ParticleBuffer::Thread::AvailableParticleListLength=NULL,*PIC::ParticleBuffer::Thread::FirstPBufferParticle=NULL;


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

  if (ParticleDataBuffer==NULL) {
    char msg[500];

    sprintf(msg,"Error: cannot allocate the particle data buffer (%ld byte). Decrease the total number of the reserved particles.",ParticleDataLength*MaxNPart);
    exit(__LINE__,__FILE__,msg);
  }  

  if (PIC::ThisThread==0) printf("$PREFIX: The total particle buffer length=%i\n",BufrerLength);

  //init the list of particles in the buffer
  for (long int ptr=0;ptr<MaxNPart-1;ptr++) {
    SetNext(ptr+1,ptr);
    SetParticleDeleted(ptr);
  }

  SetNext(-1,MaxNPart-1);
  SetParticleDeleted(MaxNPart-1);

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  FirstPBufferParticle=0;

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  //allocate the Thread::FirstPBufferParticle array and distribute particles between processes
  #pragma omp parallel
  {
    #pragma omp single
    {
      long int thread,nParticlePerThread;

      Thread::NTotalThreads=omp_get_num_threads();

      Thread::AvailableParticleListLength=new long int [Thread::NTotalThreads];
      Thread::FirstPBufferParticle=new long int [Thread::NTotalThreads];

      nParticlePerThread=MaxNPart/Thread::NTotalThreads;

      for (thread=0;thread<Thread::NTotalThreads;thread++) {
        long int nStartPart,ListLength;

        nStartPart=nParticlePerThread*thread;
        ListLength=(thread!=Thread::NTotalThreads-1) ? nParticlePerThread : MaxNPart-nParticlePerThread*(Thread::NTotalThreads-1);

        SetPrev(-1,nStartPart);
        if (nStartPart!=0) SetNext(-1,nStartPart-1);

        Thread::AvailableParticleListLength[thread]=ListLength;
        Thread::FirstPBufferParticle[thread]=nStartPart;
      }
    }
  }
#else
  #error The mode is unknown
#endif //_COMPILATION_MODE_

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

long int PIC::ParticleBuffer::GetAllPartNum() {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  return NAllPart;
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int nTotalParticles=0;

  #pragma omp parallel default(none) shared(Thread::AvailableParticleListLength,MaxNPart,nTotalParticles,PIC::nTotalThreadsOpenMP)
  {
    #pragma omp single
    {
      int threadOpenMP,nTotalAvailableParticles=0;

      for (threadOpenMP=0;threadOpenMP<PIC::nTotalThreadsOpenMP;threadOpenMP++) nTotalAvailableParticles+=Thread::AvailableParticleListLength[threadOpenMP];

      nTotalParticles=MaxNPart-nTotalAvailableParticles;
    }
  }

  return nTotalParticles;
#else
  #error The mode is unknown
#endif //_COMPILATION_MODE_
}

long int PIC::ParticleBuffer::GetParticleDataLength() {return ParticleDataLength;}

//option RandomThreadOpenMP==true can be used ONLY when the code is outside of any OpenMP sections
long int PIC::ParticleBuffer::GetNewParticle(bool RandomThreadOpenMP) {
  long int newptr;
  byte *pdataptr;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);
  FirstPBufferParticle=GetNext(pdataptr);

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread; //=omp_get_thread_num();

  thread=(RandomThreadOpenMP==false) ? omp_get_thread_num() : (int)(rnd()*Thread::NTotalThreads);

  if (Thread::AvailableParticleListLength[thread]==0) {
    printf("$PREFIX: The particle buffer is full (thread=%i)\nParticle allocation report:\nOpenMP thread\tAvailable Particles\n",PIC::ThisThread);

    for (int tt=0;tt<Thread::NTotalThreads;tt++) {
      printf("$PREFIX: %i\t%i\n",tt,Thread::AvailableParticleListLength[tt]); 
    }

    exit(__LINE__,__FILE__,"The particle buffer is full");
  }

  Thread::AvailableParticleListLength[thread]--;
  newptr=Thread::FirstPBufferParticle[thread];
  pdataptr=GetParticleDataPointer(newptr);
  Thread::FirstPBufferParticle[thread]=GetNext(pdataptr);

#else
  #error The mode is unknown
#endif //_COMPILATION_MODE_

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(pdataptr);
  #endif


  SetPrev(-1,pdataptr);
  SetNext(-1,pdataptr);

//#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (IsParticleAllocated(pdataptr)==true) exit(__LINE__,__FILE__,"Error: the particle is re-allocated");
  SetParticleAllocated(pdataptr);
//#endif

  return newptr;
}

//option RandomThreadOpenMP==true can be used ONLY when the code is outside of any OpenMP sections
long int PIC::ParticleBuffer::GetNewParticle(long int &ListFirstParticle,bool RandomThreadOpenMP) {
  long int newptr;
  byte *pdataptr;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);
  FirstPBufferParticle=GetNext(pdataptr);

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread; //=omp_get_thread_num();

  if (RandomThreadOpenMP==false) {
    thread=omp_get_thread_num();
  }
  else {
    bool found=false;
    
    //random search for a thread
    for (int i=0;i<Thread::NTotalThreads;i++) {  
      thread=(RandomThreadOpenMP==false) ? omp_get_thread_num() : (int)(rnd()*Thread::NTotalThreads);

      if (Thread::AvailableParticleListLength[thread]!=0) {
        found=true;
        break;
      }
    }

    //go through each thread
    if (found==false) for (thread=0;thread<Thread::NTotalThreads;thread++) if (Thread::AvailableParticleListLength[thread]!=0) {
      found=true;
      break;  
    }
 
    if (found==false) {
      thread=omp_get_thread_num();
      printf("$PREFIX: The particle buffer is full [MPI process=%i,OpenMP thread=%i]\n",PIC::ThisThread,thread);

      if (RandomThreadOpenMP==true) {
        printf("$PREFIX:RandomThreadOpenMP==true\nThread\tThe number of the available particles\n");
      }
      else {
        printf("$PREFIX:RandomThreadOpenMP==false\nThread\tThe number of the available particles\n");
      }

      for (thread=0;thread<Thread::NTotalThreads;thread++) printf("$PREFIX: %i\t%i\n",thread,Thread::AvailableParticleListLength[thread]);

      exit(__LINE__,__FILE__,"The particle buffer is full");
    }
  }

  if (Thread::AvailableParticleListLength[thread]==0) {
    thread=omp_get_thread_num();

    if (RandomThreadOpenMP==true) {
      printf("$PREFIX:RandomThreadOpenMP==true\nThread\tThe number of the available particles\n");
    }
    else {
      printf("$PREFIX:RandomThreadOpenMP==false\nThread\tThe number of the available particles\n");
    }

    printf("$PREFIX: The particle buffer is full [MPI process=%i,OpenMP thread=%i]\nThread\tThe number of the available particles\n",PIC::ThisThread,thread);
    for (thread=0;thread<Thread::NTotalThreads;thread++) printf("$PREFIX: %i\t%i\n",thread,Thread::AvailableParticleListLength[thread]);

    exit(__LINE__,__FILE__,"The particle buffer is full");
  }

  Thread::AvailableParticleListLength[thread]--;
  newptr=Thread::FirstPBufferParticle[thread];
  pdataptr=GetParticleDataPointer(newptr);
  Thread::FirstPBufferParticle[thread]=GetNext(pdataptr);

#else
  #error The mode is unknown
#endif //_COMPILATION_MODE_


  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(pdataptr);
  #endif

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
    if (next>=0) SetPrev(prev,next);
    ListFirstParticle=next;
  }
  else {
    if (prev>=0) SetNext(next,prev);
    if (next>=0) SetPrev(prev,next);
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


//option RandomThreadOpenMP==true can be used ONLY when the code is outside of any OpenMP sections
void PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(long int ptr,bool RandomThreadOpenMP) {

//#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (IsParticleAllocated(ptr)==false) exit(__LINE__,__FILE__,"Error: the particle is re-deleted");
  SetParticleDeleted(ptr);
//#endif

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  NAllPart--;
  SetNext(FirstPBufferParticle,ptr);
  FirstPBufferParticle=ptr;
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread; //=omp_get_thread_num();

  thread=(RandomThreadOpenMP==false) ? omp_get_thread_num() : (int)(rnd()*Thread::NTotalThreads);
  Thread::AvailableParticleListLength[thread]++;
  SetNext(Thread::FirstPBufferParticle[thread],ptr);
  Thread::FirstPBufferParticle[thread]=ptr;
#else
  #error The mode is unknown
#endif //_COMPILATION_MODE_
}


void PIC::ParticleBuffer::DeleteParticle(long int ptr,long int& ListFirstParticle) {
  ExcludeParticleFromList(ptr,ListFirstParticle);
  DeleteParticle(ptr);
}



void PIC::ParticleBuffer::CloneParticle(byte* CopyData,byte* SourceData) {
  long int next,prev;

  prev=GetPrev(CopyData);
  next=GetNext(CopyData);

  memcpy(CopyData,SourceData,ParticleDataLength*sizeof(byte));

  SetPrev(prev,CopyData);
  SetNext(next,CopyData);
}

void PIC::ParticleBuffer::CloneParticle(long int copy,long int source) {
  byte *SourceData,*CopyData;

  SourceData=GetParticleDataPointer(source);
  CopyData=GetParticleDataPointer(copy);

  CloneParticle(CopyData,SourceData);
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
unsigned long PIC::ParticleBuffer::GetChecksum(const char *msg) {
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
    if (msg==NULL) {
      sprintf(str,"Cdsmc::pbuffer CRC32 checksum: ");
    }
    else {
      sprintf(str,"Cdsmc::pbuffer CRC32 checksum (msg=%s): ",msg);
    }

    for (long int thread=0;thread<TotalThreadsNumber;thread++) sprintf(str,"%s 0x%lx ",str,buffer[thread]);

    printf("$PREFIX:%s\n",str);
    PrintErrorLog(str);
  }

  delete [] buffer;
  return sum.checksum();
}

unsigned long PIC::ParticleBuffer::GetChecksum() {
  return GetChecksum(NULL);
}

unsigned long PIC::ParticleBuffer::GetChecksum(int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"[line=%i,file=%s]",nline,fname);

  return GetChecksum(msg);
}

unsigned long PIC::ParticleBuffer::GetChecksum(int code,int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"[code=%i, line=%i,file=%s]",code,nline,fname);

  return GetChecksum(msg);
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
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_]; //,tempParticleMovingListTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];


  int nTotalThreads_OpenMP=1,thread_OpenMP;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel default(none) shared(nTotalThreads_OpenMP)
  {
    #pragma omp single
    {
      nTotalThreads_OpenMP=omp_get_num_threads();
    }
  }
#endif


  for (int thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
    node=(thread==PIC::Mesh::mesh.ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];

    if (node==NULL) continue;

    //sample the processor load
  #if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    double EndTime,StartTime=MPI_Wtime();
  #endif

    while (node!=NULL) {

      memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
//      memcpy(tempParticleMovingListTable,node->block->tempParticleMovingListTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));


      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

                //check the tempoparyly particle lists
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
                if (node->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]!=-1) exit(__LINE__,__FILE__,"Error: the temp list is not empty");
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
                for (thread_OpenMP=0;thread_OpenMP<nTotalThreads_OpenMP;thread_OpenMP++) {
                  if (*(node->block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k))!=-1) exit(__LINE__,__FILE__,"Error: the temp list is not empty");
                }
#else
#error The option is unknown
#endif

                //count the number and allocation flag of the active particles
                ParticleList=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

                while (ParticleList!=-1) {
                  double *x=PIC::ParticleBuffer::GetX(ParticleList);
                  int idim;

                  for (idim=0;idim<3;idim++) if ((x[idim]<node->xmin[idim])||(node->xmax[idim]<x[idim])) {
                    exit(__LINE__,__FILE__,"Error: a particle is outside of the block limit");
                  }

                  if (PIC::ParticleBuffer::IsParticleAllocated(ParticleList)==false) {
                    exit(__LINE__,__FILE__,"Error: a particle in the list is not allocated");
                  }

                  ++nTotalListParticles;
                  ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                  if (nTotalListParticles>MaxNPart) exit(__LINE__,__FILE__,"Error: the list particles' number exeeds the maximum number of particles that can be stored in the buffer");
                }

/*
                //check the allocation flags of the not-active particles
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
                ParticleList=FirstPBufferParticle;
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel default(none) private(ParticleList) shared(Thread::FirstPBufferParticle)
                {
                int thread;

                thread=omp_get_thread_num();
                ParticleList=Thread::FirstPBufferParticle[thread];

#endif //_COMPILATION_MODE_

                while (ParticleList!=-1) {
                  if (PIC::ParticleBuffer::IsParticleAllocated(ParticleList)==true) {
                     exit(__LINE__,__FILE__,"Error: a NON-active particle is marked as allocated");
                  }

                  ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);
                }

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
                }
#endif //_COMPILATION_MODE_
*/
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

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  if (nTotalListParticles!=NAllPart) exit(__LINE__,__FILE__,"Error: the total number of particles stored in the lists is different from that stored in the particle buffer");
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  //calculate the total number of the particles available in the particle buffer
  long int nTotalAvialableParticles=0;

  for (thread_OpenMP=0;thread_OpenMP<nTotalThreads_OpenMP;thread_OpenMP++) nTotalAvialableParticles+=Thread::AvailableParticleListLength[thread_OpenMP];

  if (nTotalAvialableParticles+nTotalListParticles!=MaxNPart) exit(__LINE__,__FILE__,"Error: the total number of particles stored in the lists is different from that stored in the particle buffer");
#else
#error The option is unknown
#endif
}

//==========================================================
//initiate the new particle
int PIC::ParticleBuffer::InitiateParticle(double *x,double *v,double *WeightCorrectionFactor,int *spec,PIC::ParticleBuffer::byte* ParticleData,int InitMode,void *node) {
  int ptr,ptrSpec;
  byte* ptrData;

  ptr=PIC::ParticleBuffer::GetNewParticle();
  ptrData=GetParticleDataPointer(ptr);

  //default settings
  SetIndividualStatWeightCorrection(1.0,ptrData);

  //set up the fields of the new particle with the user-defined data
  if (ParticleData!=NULL) {
    memcpy((void*)ptrData,(void*)ParticleData,ParticleDataLength);
    SetParticleAllocated(ptrData);
  }

  if (x!=NULL) SetX(x,ptrData);
  if (v!=NULL) SetV(v,ptrData);
  if (spec!=NULL) SetI(*spec,ptrData);
  if (WeightCorrectionFactor!=NULL) SetIndividualStatWeightCorrection(*WeightCorrectionFactor,ptrData);

  //determine the species number
  ptrSpec=(spec!=NULL) ? *spec : GetI(ptrData);

  //apply the particle tracking condition
  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(ptrData);
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,ptrSpec,ptrData,node);
  #endif

  //check if mover guiding center motion integration is used
#if _PIC_MOVER_INTEGRATOR_MODE_ ==_PIC_MOVER_INTEGRATOR_MODE__GUIDING_CENTER_
  PIC::Mover::GuidingCenter::InitiateMagneticMoment(GetI(ptr),
						    GetX(ptr),GetV(ptr),
						    ptr, node);
#endif //_PIC_MOVER_INTEGRATOR_MODE_
  


  //add the paticle to the cell's particle list
  long int FirstCellParticle;
  int iCell,jCell,kCell;

  switch (InitMode) {
  case _PIC_INIT_PARTICLE_MODE__ADD2LIST_:
    PIC::Mesh::mesh.fingCellIndex(x,iCell,jCell,kCell,(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)node);
    FirstCellParticle=((cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)node)->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

    SetNext(FirstCellParticle,ptr);
    SetPrev(-1,ptr);

    if (FirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstCellParticle);
    ((cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)node)->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=ptr;

    break;
  case _PIC_INIT_PARTICLE_MODE__MOVE_:
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,((cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)node)->block->GetLocalTimeStep(ptrSpec)*rnd(),(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)node);

    break;
  default:
    exit(__LINE__,__FILE__,"Error: the opiton is not found");
  }

  return ptr;
}


//==========================================================
//rebalance the particle lists
void PIC::ParticleBuffer::Thread::RebalanceParticleList() {
  int SendThread=0,RecvThread=0,ExchangeListLength=0;
  long int ExchangeListBegin=-1,ExchangeListEnd=-1;
  int thread,nParticlePerThread,i;

  for (thread=0,nParticlePerThread=0;thread<NTotalThreads;thread++) nParticlePerThread+=AvailableParticleListLength[thread];
  nParticlePerThread/=NTotalThreads;

  //redistribute available particle slots netween threads on the same node
  do {
    if ((AvailableParticleListLength[SendThread]>nParticlePerThread+NTotalThreads)&&(nParticlePerThread-NTotalThreads>AvailableParticleListLength[RecvThread])) {
      //exchange particles between threads SendThread and RecvThread
      ExchangeListLength=min(AvailableParticleListLength[SendThread]-(nParticlePerThread+NTotalThreads),
          nParticlePerThread-NTotalThreads-AvailableParticleListLength[RecvThread]);

      ExchangeListBegin=FirstPBufferParticle[SendThread];
      ExchangeListEnd=FirstPBufferParticle[SendThread];

      for (i=0;i<ExchangeListLength-1;i++) ExchangeListEnd=PIC::ParticleBuffer::GetNext(ExchangeListEnd);

      //remove the particle list from SendThread
      FirstPBufferParticle[SendThread]=PIC::ParticleBuffer::GetNext(ExchangeListEnd);
      AvailableParticleListLength[SendThread]-=ExchangeListLength;

      //add the particle list to RecvThread
      PIC::ParticleBuffer::SetNext(FirstPBufferParticle[RecvThread],ExchangeListEnd);
      FirstPBufferParticle[RecvThread]=ExchangeListBegin;
      AvailableParticleListLength[RecvThread]+=ExchangeListLength;
    }

    //advance SendThread and RecvThread numbers if needed
    if (AvailableParticleListLength[SendThread]<=nParticlePerThread+NTotalThreads) SendThread++;
    if (AvailableParticleListLength[RecvThread]>=nParticlePerThread-NTotalThreads) RecvThread++;
  }
  while ((SendThread<NTotalThreads)&&(RecvThread<NTotalThreads));

  //check that all threads has particles
  for (thread=0;thread<NTotalThreads;thread++) if (AvailableParticleListLength[thread]<=0) {
    char msg[100];

    sprintf(msg,"Error: AvailableParticleListLength[thread]==0 for thread=%i",thread);
    exit(__LINE__,__FILE__,msg);
  }

}
