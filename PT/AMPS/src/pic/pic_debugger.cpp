/*
 * pic_debugger.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: vtenishe
 */
//$Id$
//contains functions used for debugging AMPS

#include "pic.h"


//Save particle data into a debugger data stream
void PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(void* data,int length,int nline,const char* fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"%s, line %i",fname,nline);
  PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(data,length,msg);
}


void PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(void* data,int length,const char* msg) {

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

    sprintf(fn,"DebuggerStream.thread=%i.dbg",PIC::ThisThread);
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

    sprintf(fn,"DebuggerStream.thread=%i.dbg",PIC::ThisThread);
    fout=fopen(fn,"a");

    for (int i=0;i<BufferPointer;i++) fprintf(fout,"%i: 0x%lx\t%s\n",StreamBuffer[i].CollCounter,StreamBuffer[i].CheckSum,StreamBuffer[i].CallPoint);

    BufferPointer=0;
    fclose(fout);
  }
}


//InfiniteLoop==false ==> no problem found; InfiniteLoop==true => the actual number of particles does not consider the that in teh particle buffer
bool PIC::Debugger::InfiniteLoop(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  bool res=false;
  long int ptr;
  static long int nAllCountedParticles=0;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;
  if (startNode==PIC::Mesh::mesh.rootTree) nAllCountedParticles=0;


  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InfiniteLoop(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block
    int nTotalThreads_OpenMP=1,thread_OpenMP;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    nTotalThreads_OpenMP=omp_get_thread_num();
#endif


    for (thread_OpenMP=0;thread_OpenMP<nTotalThreads_OpenMP;thread_OpenMP++) for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
      ptr=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      ptr=*(startNode->block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k));
#else
#error The option is unknown
#endif

      while (ptr!=-1) {
        ptr=PIC::ParticleBuffer::GetNext(ptr);

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          FindDoubleReferencedParticle();

          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //collect the information from all processors;
    int Flag,FlagTable[PIC::nTotalThreads];

    Flag=(nAllCountedParticles==PIC::ParticleBuffer::NAllPart) ? true : false;
    MPI_Gather(FlagTable,1,MPI_INT,&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      for (int thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==false) {
        printf("$PREFIX: the actual number of particles does not coniside with that in the particle buffer (Thread=%i,line=%i,file=%s)\n",thread,__LINE__,__FILE__);
        res=true;
      }
    }

    MPI_Bcast(&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    nAllCountedParticles=0;
  }

  return res;
}



void PIC::Debugger::FindDoubleReferencedParticle(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  long int ptr,ptrNext;
  static long int nAllCountedParticles=0;

  static char* ParticleAllocationTable=NULL;
  static long int* ptrPrevTable=NULL;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;

  //create table of allocated/not-allocated flag table
  if (startNode==PIC::Mesh::mesh.rootTree) {
    long int cnt=0;

    nAllCountedParticles=0;
    ParticleAllocationTable=new char [PIC::ParticleBuffer::MaxNPart];
    ptrPrevTable=new long int [PIC::ParticleBuffer::MaxNPart];

    //the definition of the bits of ParticleAllocationTable[]
    //ParticleAllocationTable[] & 1: == 0 --> particle is not allocated; == 1 --> the particle is allocated
    //ParticleAllocationTable[] & 2: == 0 --> the particle is not lonked yet; == 2 --> the particle is already linked

    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) ptrPrevTable[ptr]=-1,ParticleAllocationTable[ptr]=1; //particle is allocated

    if (PIC::ParticleBuffer::FirstPBufferParticle!=-1) ParticleAllocationTable[PIC::ParticleBuffer::FirstPBufferParticle]=2; //particle is links and not allocated

    for (cnt=0,ptr=PIC::ParticleBuffer::FirstPBufferParticle;ptr!=-1;ptr=ptrNext) {
      ptrNext=PIC::ParticleBuffer::GetNext(ptr);

      if (ptrNext!=-1) {
        if ((ParticleAllocationTable[ptrNext]&2)==2) {
          printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);

          exit(__LINE__,__FILE__,"Error: an un-allocated particle is double-referenced");
        }
        else ParticleAllocationTable[ptrNext]=2,ptrPrevTable[ptrNext]=ptr;
      }

      if (++cnt>PIC::ParticleBuffer::MaxNPart) {
        exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
      }

      ptr=ptrNext;
    }
  }

  //check the particle lists
  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) FindDoubleReferencedParticle(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block

    int nTotalThreads_OpenMP=1,thread_OpenMP;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    nTotalThreads_OpenMP=omp_get_thread_num();
  #endif

    for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (int npass=0;npass<2;npass++)  for (thread_OpenMP=0;thread_OpenMP<((npass==0) ? 1 : nTotalThreads_OpenMP);thread_OpenMP++) {
      if (npass==0) ptr=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
      else {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
        ptr=startNode->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        ptr=(*startNode->block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k));
#endif
      }

//      ptr=(npass==0) ? startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)] : startNode->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      if (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        if ((ParticleAllocationTable[ptr]&2)==2) {
          printf("Error: the first particle in the list is referenced: %ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptr],ptr);
          exit(__LINE__,__FILE__,"Error: the first particle in the list is referenced");
        }
        else ParticleAllocationTable[ptr]|=2;

        if (PIC::ParticleBuffer::GetPrev(ptr)>=0) {
          exit(__LINE__,__FILE__,"Error: the first particle in the list has non-negarive value of GetPrev");
        }
      }


      while (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        ptrNext=PIC::ParticleBuffer::GetNext(ptr);

        if (ptrNext!=-1) {
          if ((ParticleAllocationTable[ptrNext]&1)==0) {
            exit(__LINE__,__FILE__,"Error: reference to an un-allocated particles is found in the particle list");
          }

          if ((ParticleAllocationTable[ptrNext]&2)==2) {
            printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);
            exit(__LINE__,__FILE__,"Error: have found double-referenced particle in the list");
          }

          if (PIC::ParticleBuffer::GetPrev(ptrNext)!=ptr) {
            exit(__LINE__,__FILE__,"Error: PIC::ParticleBuffer::GetPrev(ptrNext) != ptr");
          }

          ParticleAllocationTable[ptrNext]|=2;
          ptrPrevTable[ptrNext]=ptr;
        }

        ptr=ptrNext;

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) {
      if ((ParticleAllocationTable[ptr]&2)==0) {
        exit(__LINE__,__FILE__,"Error: found particles that is not referenced at all");
      }
    }

    delete [] ParticleAllocationTable;
    delete [] ptrPrevTable;
  }
}


//catch the out of limit value in the sample buffer (check only the base quantaty)
void PIC::Sampling::CatchOutLimitSampledValue() {
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int s,i,j,k;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber;
  char *SamplingBuffer;



  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
      LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
      cell=block->GetCenterNode(LocalCellNumber);

      SamplingBuffer=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;


      for (s=0;s<PIC::nTotalSpecies;s++) {
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleWeghtRelativeOffset)),1,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberRelativeOffset)),1,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberDensityRelativeOffset)),1,__LINE__,__FILE__);

        PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocityRelativeOffset)),DIM,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocity2RelativeOffset)),DIM,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleSpeedRelativeOffset)),1,__LINE__,__FILE__);
      }

    }
  }

#endif
#endif
}


//==========================================================================================
//get checksum of the corner and center node associated data
unsigned long int PIC::Debugger::SaveCornerNodeAssociatedDataSignature(long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int i,j,k;
  PIC::Mesh::cDataCornerNode *CornerNode;
  PIC::Mesh::cDataBlockAMR *block;

  static int EntryCounter;  //used to ensure the same allocation of the blocks and nodes
  static CRC32 CheckSum,SingleVectorCheckSum;
  static CMPI_channel pipe(1000000);
  static int initflag=false;
  static FILE *fout=NULL;

  //coundate of the fucntion calls
  static int nCallCounter=0;

  if (startNode==NULL) {
    startNode=PIC::Mesh::mesh.rootTree;
    EntryCounter=0;
    CheckSum.clear();

    nCallCounter++;

    if ((fnameOutput!=NULL)&&(PIC::ThisThread==0)) fout=fopen(fnameOutput,"w");

    if (initflag==false) {
      initflag=true;

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else {
        pipe.openSend(0);
      }
    }
  }

  //add the associated node data
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((startNode->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      EntryCounter++;
      CheckSum.add(EntryCounter);

      block=startNode->block;

      for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_+1;k++) {
        for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_+1;j++) {
          for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_+1;i++) {
            EntryCounter++;
            CheckSum.add(EntryCounter);
            SingleVectorCheckSum.clear();

            if (startNode->Thread==0) {
              //the associated data is located the the root
              if (block!=NULL) if ((CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                CheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

                if (fnameOutput!=NULL) {
                  SingleVectorCheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
                }
              }

              if (fnameOutput!=NULL) {
                fprintf(fout,"node: id=%i, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
              }
            }
            else {
              if (PIC::ThisThread==0) {
                //recieve the data vector
                bool DataSendMode;

                pipe.recv(DataSendMode,startNode->Thread);

                if (DataSendMode==true) {
                  char *ptr;

                  ptr=pipe.recvPointer<char>(PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,startNode->Thread);
                  CheckSum.add(ptr,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

                  if (fnameOutput!=NULL) {
                    SingleVectorCheckSum.add(ptr,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
                  }
                }

                if (fnameOutput!=NULL) {
                  fprintf(fout,"node: id=%i, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
                }
              }
              else {
                //send the data vector
                bool DataSendMode;

                if (block!=NULL) {
                  if ((CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                    DataSendMode=true;
                    pipe.send(DataSendMode);

                    pipe.send(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
                  }
                  else {
                    DataSendMode=false;
                    pipe.send(DataSendMode);
                  }
                }
                else {
                  DataSendMode=false;
                  pipe.send(DataSendMode);
                }

              }
            }

          }

        }
      }
    }
  }
  else {
    //add daugher blocks
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) SaveCornerNodeAssociatedDataSignature(nline,fnameSource,fnameOutput,downNode);
  }

  //output the checksum
  if (startNode==PIC::Mesh::mesh.rootTree) {
    pipe.flush();

    if (PIC::ThisThread==0) {
      char msg[500];

      sprintf(msg," line=%ld, file=%s (Call Counter=%i)",nline,fnameSource,nCallCounter);
      CheckSum.PrintChecksumSingleThread(msg);

      if (fnameOutput!=NULL) {
        fclose(fout);
        fout=NULL;
      }
    }
  }

  return CheckSum.checksum();
}

unsigned long int PIC::Debugger::GetCornerNodeAssociatedDataSignature(long int nline,const char* fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return SaveCornerNodeAssociatedDataSignature(nline,fname,NULL,startNode);
}

unsigned long int PIC::Debugger::GetCenterNodeAssociatedDataSignature(long int nline,const char* fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return SaveCenterNodeAssociatedDataSignature(nline,fname,NULL,startNode);
}

unsigned long int PIC::Debugger::SaveCenterNodeAssociatedDataSignature(long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int i,j,k;
  PIC::Mesh::cDataCenterNode *CenterNode;
  PIC::Mesh::cDataBlockAMR *block;

  static int EntryCounter;  //used to ensure the same allocation of the blocks and nodes
  static CRC32 CheckSum,SingleVectorCheckSum;
  static CMPI_channel pipe(1000000);
  static int initflag=false;
  static FILE *fout=NULL;

  //coundate of the fucntion calls
  static int nCallCounter=0;

  if (startNode==NULL) {
    startNode=PIC::Mesh::mesh.rootTree;
    EntryCounter=0;
    CheckSum.clear();

    if ((fnameOutput!=NULL)&&(PIC::ThisThread==0)) fout=fopen(fnameOutput,"w");

    nCallCounter++;

    if (initflag==false) {
      initflag=true;

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else {
        pipe.openSend(0);
      }
    }
  }

  //add the associated node data
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((startNode->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      EntryCounter++;
      CheckSum.add(EntryCounter);

      block=startNode->block;

      for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
        for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) {
          for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
            EntryCounter++;
            CheckSum.add(EntryCounter);
            SingleVectorCheckSum.clear();

            if (startNode->Thread==0) {
              //the associated data is located the the root
              if (block!=NULL) if ((CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k)))!=NULL) {
                CheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);

                if (fnameOutput!=NULL) {
                  SingleVectorCheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
                }
              }

              if (fnameOutput!=NULL) {
                fprintf(fout,"node: id=%i, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
              }

            }
            else {
              if (PIC::ThisThread==0) {
                //recieve the data vector
                bool DataSendMode;

                pipe.recv(DataSendMode,startNode->Thread);

                if (DataSendMode==true) {
                  char *ptr;

                  ptr=pipe.recvPointer<char>(PIC::Mesh::cDataCenterNode::totalAssociatedDataLength,startNode->Thread);
                  CheckSum.add(ptr,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);

                  if (fnameOutput!=NULL) {
                    SingleVectorCheckSum.add(ptr,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
                  }
                }

                if (fnameOutput!=NULL) {
                  fprintf(fout,"node: id=%i, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
                }
              }
              else {
                //send the data vector
                bool DataSendMode;

                if (block!=NULL) {
                  if ((CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k)))!=NULL) {
                    DataSendMode=true;
                    pipe.send(DataSendMode);

                    pipe.send(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
                  }
                  else {
                    DataSendMode=false;
                    pipe.send(DataSendMode);
                  }
                }
                else {
                  DataSendMode=false;
                  pipe.send(DataSendMode);
                }

              }
            }

          }

        }
      }
    }
  }
  else {
    //add daugher blocks
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) SaveCenterNodeAssociatedDataSignature(nline,fnameSource,fnameOutput,downNode);
  }

  //output the checksum
  if (startNode==PIC::Mesh::mesh.rootTree) {
    pipe.flush();

    if (PIC::ThisThread==0) {
      char msg[500];

      sprintf(msg," line=%ld, file=%s (Call Counter=%i)",nline,fnameSource,nCallCounter);
      CheckSum.PrintChecksumSingleThread(msg);

      if (fnameOutput!=NULL) {
        fclose(fout);
        fout=NULL;
      }
    }
  }

  return CheckSum.checksum();
}


//=====================================================================================
//get signature describe the particle population
unsigned long int PIC::Debugger::GetParticlePopulationSignature(long int nline,const char* fname) {
  CRC32 Checksum;
  PIC::ParticleBuffer::byte *ParticleDataPtr,ParticleBuffer[PIC::ParticleBuffer::ParticleDataLength];
  int i,j,k,ptr;

  CMPI_channel pipe;
  pipe.init(1000000);

  //init the particle buffer
  for (i=0;i<PIC::ParticleBuffer::ParticleDataLength;i++) ParticleBuffer[i]=0;

  const int CommunicationCompleted_SIGNAL=0;
  const int ParticleDataSend_SIGNAL=1;
  const int BockDataStarted_SIGNAL=2;

  if (PIC::ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);

  //loop through all blocks
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    if ((node->Thread==0)&&(PIC::ThisThread==0)) {
      //the block belongs to the root

      if (node->block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
        ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        //collect signature
        while (ptr!=-1) {
          //copy the state vector of the particle without 'next' and 'prev'
          ParticleDataPtr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          PIC::ParticleBuffer::CloneParticle(ParticleBuffer,ParticleDataPtr);

          //add signature of the particle
          Checksum.add(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength);
          ptr=PIC::ParticleBuffer::GetNext(ptr);
        }
      }
    }
    else if (PIC::ThisThread==0) {
      //this is the root BUT the block belongs to another MPI process
      unsigned long t;
      int Signal;

      pipe.recv(Signal,node->Thread);

      switch (Signal) {
      case BockDataStarted_SIGNAL:
        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          pipe.recv(Signal,node->Thread);

          //collect signatures
          while (Signal!=CommunicationCompleted_SIGNAL) {
            pipe.recv(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength,node->Thread);
            Checksum.add(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength);

            pipe.recv(Signal,node->Thread);
          }
        }

        break;
      case CommunicationCompleted_SIGNAL:
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the sigmal is not recognized");
      }


    }
    else if (node->Thread==PIC::ThisThread) {
      //this is NOT the root BUT the block belongs to the current MPI process
      //loop through all cells and particles

      if (node->block!=NULL) {
        pipe.send(BockDataStarted_SIGNAL);

        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            pipe.send(ParticleDataSend_SIGNAL);

            //copy the state vector of the particle without 'next' and 'prev'
            ParticleDataPtr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            PIC::ParticleBuffer::CloneParticle(ParticleBuffer,ParticleDataPtr);
            pipe.send(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength);

            ptr=PIC::ParticleBuffer::GetNext(ptr);
          }

          pipe.send(CommunicationCompleted_SIGNAL);
        }

      }
      else {
        pipe.send(CommunicationCompleted_SIGNAL);
      }

    }
  }

  //output the checksum
  if (PIC::ThisThread==0) {
    char msg[500];

    pipe.closeRecvAll();

    sprintf(msg," line=%ld, file=%s",nline,fname);
    Checksum.PrintChecksumSingleThread(msg);
  }
  else pipe.closeSend();

  return Checksum.checksum();
}


//=========================================================================================================
//save the map of the domain decomposition
void PIC::Debugger::SaveDomainDecompositionMap(long int nline,const char* fname,int Index) {
  FILE *fout;
  char FullFileName[100];
  int id,i,j,iface,iedge,icorner;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,*neibNode;

  sprintf(FullFileName,"DomainDecompositionMap.thread=%i(line=%i,file=%s,Index=%i).dat",PIC::ThisThread,nline,fname,Index);
  fout=fopen(FullFileName,"w");

  fprintf(fout,"VARIABLES=\"NodeTempID\", \"Thread\", \"Face Neib\", \"Edge Neib\", \"Corner Neib\"\n");

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    fprintf(fout,"node->Temp_ID=%i, thread=%i\n",node->Temp_ID,node->Thread);

    //face neib
    for (iface=0;iface<6;iface++) for (i=0;i<2;i++) for (j=0;j<2;j++) {
      if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"iface=%i,i=%i,j=%i,neib=%i\n",iface,i,j,id);
    }

    //edge neib
    for (iedge=0;iedge<12;iedge++) for (i=0;i<2;i++) {
      if ((neibNode=node->GetNeibEdge(iedge,i))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"iedge=%i,i=%i,neib=%i\n",iface,i,id);
    }

    //corner
    for (icorner=0;icorner<8;icorner++) {
      if ((neibNode=node->GetNeibCorner(icorner))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"icorner=%i,neib=%i\n",icorner,id);
    }

    //end the line
    fprintf(fout,"\n");
  }

  //close the file
  fclose(fout);
}


//====================================================================================================
//output the debug debug particle data
list<PIC::Debugger::ParticleDebugData::cDebugData> PIC::Debugger::ParticleDebugData::DebugParticleData;

//method for sorting the debug particle list
bool PIC::Debugger::ParticleDebugData::CompareParticleDebugData(const PIC::Debugger::ParticleDebugData::cDebugData& first, const PIC::Debugger::ParticleDebugData::cDebugData& second) {
  return (first.initCheckSum < second.initCheckSum);
}

//accumulate the partilce debug data
void PIC::Debugger::ParticleDebugData::AddParticleDebugData(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,bool InitChckSumMode) {
  int i,j,k,di,dj,dk;
  double *x;
  cDebugData p;

  //add the particle data
  CRC32 CheckSum;

  if (InitChckSumMode==true) {
    p.initCheckSum=PIC::ParticleBuffer::GetParticleSignature(ptr);

    x=PIC::ParticleBuffer::GetX(ptr);
    PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false);

    p.i=i;
    p.j=j;
    p.k=k;
    p.nodeid=node->Temp_ID;
    p.ptr=ptr;
  }
  else {
    //look for the particle
    list<cDebugData>::iterator pp;

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      if (pp->ptr==ptr) {
        //the particle has been found
        pp->finalCheckSum=PIC::ParticleBuffer::GetParticleSignature(ptr);
        return;
      }
    }

    //when come to this point --> the particle has not beed found
    exit(__LINE__,__FILE__,"the particle has not been found");
  }

  //add checksum of the corner and center nodes
  for (dk=0;dk<2;dk++) {
    for (dj=0;dj<2;dj++)  {
      for (di=0;di<2;di++) {
        PIC::Mesh::cDataCornerNode *CornerNode;

        if (node->Thread==PIC::ThisThread) {
          if (node->block!=NULL) {
            CheckSum.clear();

            CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(p.i+di,p.j+dj,p.k+dk));
            CheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

            p.CornerNodeChecksum[di][dj][dk]=CheckSum.checksum();
          }
        }
      }
    }
  }


  for (dk=-1;dk<1;dk++) {
    for (dj=-1;dj<1;dj++)  {
      for (di=-1;di<1;di++) {
        PIC::Mesh::cDataCenterNode *CenterNode;

        if (node->Thread==PIC::ThisThread) {
          if (node->block!=NULL) {
            CheckSum.clear();

            CenterNode=node->block->GetCenterNode(_getCenterNodeLocalNumber(p.i+di,p.j+dj,p.k+dk));
            CheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);

            p.CenterNodeChecksum[1+di][1+dj][1+dk]=CheckSum.checksum();
          }
        }
      }
    }
  }


  //add the particle data to the list
  DebugParticleData.push_front(p);
}

//============================================================================================
//output particle debug data
void PIC::Debugger::ParticleDebugData::OutputParticleDebugData(int nLineSource,const char *fNameSource,int Index) {
  int size;
  cDebugData p;
  list<cDebugData>::iterator pp;

  CMPI_channel pipe;
  pipe.init(1000000);

  if (PIC::ThisThread==0) {
      pipe.openRecvAll();
  }
  else pipe.openSend(0);

  //collect the debug data from all MPI processes
  if (PIC::ThisThread==0) {
    for (int thread=1;thread<PIC::nTotalThreads;thread++) {
      pipe.recv(size,thread);

      for (int i=0;i<size;i++) {
        pipe.recv(p,thread);
        DebugParticleData.push_front(p);
      }
    }


    //sort and output the particle data
    DebugParticleData.sort(CompareParticleDebugData);

    //output the debug data
    FILE *fout;
    char fname[100];
    int di,dj,dk;

    sprintf(fname,"ParticleDebugData(nline=%i,file=%s).Index=%i.dat",nLineSource,fNameSource,Index);
    fout=fopen(fname,"w");

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      fprintf(fout,"particle(init)=0x%lx, particle(final)=0x%lx, nodeid=%i, i=%i, j=%i, k=%i\n",
        pp->initCheckSum,pp->finalCheckSum,pp->nodeid,pp->i,pp->j,pp->k);

      fprintf(fout,"CenterNodeChecksum:\n");

      for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
        fprintf(fout,"di=%i, dj=%i, dk=%i, CenterNodeChecksum=0x%lx\n",di-1,dj-1,dk-1,pp->CenterNodeChecksum[di][dj][dk]);
      }


      fprintf(fout,"CornerNodeChecksum:\n");

      for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
        fprintf(fout,"di=%i, dj=%i, dk=%i, CornerNodeChecksum=0x%lx\n",di,dj,dk,pp->CornerNodeChecksum[di][dj][dk]);
      }

      fprintf(fout,"\n");
    }

    fclose(fout);
  }
  else {
    size=DebugParticleData.size();
    pipe.send(size);

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      p=*pp;
      pipe.send(p);
    }
  }


   //cloe the pipe
   if (PIC::ThisThread==0) {
     pipe.closeRecvAll();
   }
   else pipe.closeSend();

   //remove the particle debug data list
   DebugParticleData.clear();
}

//========================================================================================================================
//save signatures of the nodes
void PIC::Debugger::SaveNodeSignature(int nline,const char *fname) {
  int i,j,k;
  static FILE *fout=NULL;
  static int cnt;
  static int ncall=0;

  ncall++;

  if ((fout==NULL)&&(PIC::ThisThread==0))  {
    fout=fopen("SaveNodeSignature.dat","w");
  }

  unsigned long s=PIC::Debugger::GetCornerNodeAssociatedDataSignature(nline,fname);
  unsigned long int p=PIC::Debugger::GetParticlePopulationSignature(nline,fname);

  //'corner' data
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_+1;k++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_+1;j++)  {
        for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_+1;i++) {
          unsigned long int cs;
          CRC32 CheckSum;
          PIC::Mesh::cDataCornerNode *CornerNode;

          if (node->Thread==PIC::ThisThread) {
            if (node->block!=NULL) {
              CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
              CheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
            }

            cs=CheckSum.checksum();

            if (PIC::ThisThread!=0) {
              //send the checksum to the root
              MPI_Send(&cs,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
            }
          }
          else if (PIC::ThisThread==0) {
            //recieve the checksum
            MPI_Status status;

            MPI_Recv(&cs,1,MPI_UNSIGNED_LONG,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          }

          //print the checksum to a file
          if (PIC::ThisThread==0) {
            fprintf(fout,"Corner CheckSum=0x%lx, i=%i, j=%i, k=%i,id=%i, ncall=%i, line=%i,file=%s \n",cs,i,j,k,node->Temp_ID,ncall,nline,fname);
          }
        }
      }
    }
  }

  //'center' data
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)  {
         for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
           unsigned long int cs;
           CRC32 CheckSum;
           PIC::Mesh::cDataCenterNode *CenterNode;

           if (node->Thread==PIC::ThisThread) {
             if (node->block!=NULL) {
               CenterNode=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
               CheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
             }

             cs=CheckSum.checksum();

             if (PIC::ThisThread!=0) {
               //send the checksum to the root
               MPI_Send(&cs,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
             }
           }
           else if (PIC::ThisThread==0) {
             //recieve the checksum
             MPI_Status status;

             MPI_Recv(&cs,1,MPI_UNSIGNED_LONG,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
           }

           //print the checksum to a file
           if (PIC::ThisThread==0) {
             fprintf(fout,"Center CheckSum=0x%lx, i=%i, j=%i, k=%i,id=%i, ncall=%i, line=%i,file=%s \n",cs,i,j,k,node->Temp_ID,ncall,nline,fname);
           }
        }
      }
    }
  }

  if (fout!=NULL) {
    fflush(fout);
  }
}





















