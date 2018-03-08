//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that control the interprocessor communication of the code

#include "pic.h"

long int PIC::Parallel::sendParticleCounter=0,PIC::Parallel::recvParticleCounter=0,PIC::Parallel::IterationNumberAfterRebalancing=0;
double PIC::Parallel::RebalancingTime=0.0,PIC::Parallel::CumulativeLatency=0.0;
double PIC::Parallel::EmergencyLoadRebalancingFactor=3.0;
double PIC::Parallel::Latency=0.0;

//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::fUserDefiendProcessNodeAssociatedData PIC::Parallel::ProcessCenterNodeAssociatedData=NULL,PIC::Parallel::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::fUserDefiendProcessNodeAssociatedData PIC::Parallel::CopyCenterNodeAssociatedData=NULL,PIC::Parallel::CopyCornerNodeAssociatedData=NULL;

//====================================================
//Exchange particles between Processors
void PIC::Parallel::ExchangeParticleData() {
  int From,To;
  long int Particle,NextParticle,newParticle,LocalCellNumber=-1;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *sendNode=NULL,*recvNode=NULL;


#if DIM == 3
//  cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
  cAMRnodeID nodeid;
#elif DIM == 2
  cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
#else
  cMeshAMR1d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
#endif


  //The signals
  int Signal;
  const int _NEW_BLOCK_ID_SIGNAL_=       0;
  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;

#if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
#endif




  int sendProcList[PIC::Mesh::mesh.nTotalThreads],nSendProc=0;
  int recvProcList[PIC::Mesh::mesh.nTotalThreads],nRecvProc=0;
  int sendProcVector[PIC::Mesh::mesh.nTotalThreads],exchangeProcMatrix[PIC::Mesh::mesh.nTotalThreads*PIC::Mesh::mesh.nTotalThreads];
  int thread;

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) sendProcList[thread]=-1,recvProcList[thread]=-1,sendProcVector[thread]=0;


  //local copy of the block's cells
//  int cellListLength=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();

  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  //PIC::Mesh::cDataCenterNode *cellList[cellListLength],*cell;

  //calculate the number of bytes that will be send
  for (To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
      bool CommunicationInitialed_BLOCK_;
      int iCell,jCell,kCell;

      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;
        memcpy(FirstCellParticleTable,sendNode->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));


        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

          if  (Particle!=-1) {
            if (CommunicationInitialed_BLOCK_==false) {
              sendProcVector[To]+=sizeof(nodeid)+sizeof(int);

              CommunicationInitialed_BLOCK_=true;
            }

            sendProcVector[To]+=sizeof(int)+sizeof(LocalCellNumber);

            while (Particle!=-1) {
              sendProcVector[To]+=sizeof(int)+PIC::ParticleBuffer::ParticleDataLength;

              Particle=PIC::ParticleBuffer::GetNext(Particle);
            }


          }
        }
      }


      //end the part of the sender
      if (sendProcVector[To]!=0) {
        sendProcVector[To]+=sizeof(int);
        sendProcList[nSendProc++]=To;
      }
   }


  //collect the data exchenge matrix
  Latency=MPI_Wtime();
  MPI_Allgather(sendProcVector,PIC::Mesh::mesh.nTotalThreads,MPI_INT,exchangeProcMatrix,PIC::Mesh::mesh.nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
  Latency=MPI_Wtime()-Latency;

  //extract the list of the processors that will send information to 'ThisThread'
  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) if (exchangeProcMatrix[thread*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread]!=0) recvProcList[nRecvProc++]=thread;

  //Allocate all nessesary data exchange buffers
  int nproc;
  char *recvDataBuffer[PIC::Mesh::mesh.nTotalThreads],*sendDataBuffer[PIC::Mesh::mesh.nTotalThreads];

  for (nproc=0;nproc<nRecvProc;nproc++) {
    thread=recvProcList[nproc];
    recvDataBuffer[nproc]=new char[exchangeProcMatrix[thread*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread]];
  }

  for (nproc=0;nproc<nSendProc;nproc++) {
    thread=sendProcList[nproc];
    sendDataBuffer[nproc]=new char[sendProcVector[thread]];
  }

  //collect the send data
  int offset;
  char *buffer;
  MPI_Request SendRequest[PIC::Mesh::mesh.nTotalThreads],RecvRequest[PIC::Mesh::mesh.nTotalThreads];

  sendParticleCounter=0;
  recvParticleCounter=0;

  //the part of the sender
  for (nproc=0;nproc<nSendProc;nproc++) {
      bool CommunicationInitialed_BLOCK_;
      int iCell,jCell,kCell;
      bool CellParticleTableModified;

      To=sendProcList[nproc];
      offset=0;
      buffer=sendDataBuffer[nproc];

      //reset the proceesed flaf for the blocks to be send
      //send the nodes' data
      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;


        //memcpy(cellList,sendNode->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));
        memcpy(FirstCellParticleTable,sendNode->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
        CellParticleTableModified=false;

        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          Particle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

          if  (Particle!=-1) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);

            if (CommunicationInitialed_BLOCK_==false) {
              nodeid=sendNode->AMRnodeID;

              //pipe.send(_NEW_BLOCK_ID_SIGNAL_);
              *((int*)(buffer+offset))=_NEW_BLOCK_ID_SIGNAL_;
              offset+=sizeof(int);

              #if DIM == 3
              *((cAMRnodeID*)(buffer+offset))=nodeid;
              #else
              exit(__LINE__,__FILE__,"Error: not implemetned");
              #endif

              offset+=sizeof(nodeid);

              CommunicationInitialed_BLOCK_=true;
            }

            *((int*)(buffer+offset))=_CENTRAL_NODE_NUMBER_SIGNAL_;
            offset+=sizeof(int);

            *((long int*)(buffer+offset))=LocalCellNumber;
            offset+=sizeof(long int);

            while (Particle!=-1) {
              *((int*)(buffer+offset))=_NEW_PARTICLE_SIGNAL_;
              offset+=sizeof(int);

              PIC::ParticleBuffer::PackParticleData(buffer+offset,Particle);
              offset+=PIC::ParticleBuffer::ParticleDataLength;
              sendParticleCounter++;

              NextParticle=PIC::ParticleBuffer::GetNext(Particle);
              PIC::ParticleBuffer::DeleteParticle_withoutTrajectoryTermination(Particle,true);
              Particle=NextParticle;
            }

            FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=-1;
            CellParticleTableModified=true;
          }
        }

        if (CellParticleTableModified==true) memcpy(sendNode->block->FirstCellParticleTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
      }

      *((int*)(buffer+offset))=_END_COMMUNICATION_SIGNAL_;
      offset+=sizeof(int);

      if ((offset!=sendProcVector[To])||(offset!=exchangeProcMatrix[PIC::Mesh::mesh.ThisThread*PIC::Mesh::mesh.nTotalThreads+To])) exit(__LINE__,__FILE__,"Error: the data anount to be send is not consistent");

      //end the part of the sender - initiate non-blocks send
      MPI_Isend(buffer,offset,MPI_CHAR,To,0,MPI_GLOBAL_COMMUNICATOR,SendRequest+nproc);
   }


  //initiate reciving of the data
  for (nproc=0;nproc<nRecvProc;nproc++) {
    thread=recvProcList[nproc];
    MPI_Irecv(recvDataBuffer[nproc],exchangeProcMatrix[thread*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread],MPI_CHAR,thread,0,MPI_GLOBAL_COMMUNICATOR,RecvRequest+nproc);
  }


  //recieve particles
  int DataStillToRecieve=nRecvProc;
  MPI_Status status;
  int flag;
  int iCell=-10,jCell=-10,kCell=-10;

  while (DataStillToRecieve!=0) {

    //determine the processor that has finished sending the data
    do {
      for (nproc=0;nproc<nRecvProc;nproc++) if (recvDataBuffer[nproc]!=NULL) {
        MPI_Test(RecvRequest+nproc,&flag,&status);
        if (flag==true) break;

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
        break;
#endif
      }
    }
    while (flag==false);

    buffer=recvDataBuffer[nproc];
    offset=0;
    From=recvProcList[nproc];

    //recieve the data
    Signal=*((int*)(buffer+offset));
    offset+=sizeof(int);

     while (Signal!=_END_COMMUNICATION_SIGNAL_) {

       switch (Signal) {
       case _NEW_BLOCK_ID_SIGNAL_ :
         #if DIM == 3
         nodeid=*((cAMRnodeID*)(buffer+offset));
         #else
         exit(__LINE__,__FILE__,"Error: not implemetned");
         #endif

         offset+=sizeof(nodeid);
         recvNode=PIC::Mesh::mesh.findAMRnodeWithID(nodeid);

         memcpy(FirstCellParticleTable,recvNode->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

         if (recvNode->block==NULL) exit(__LINE__,__FILE__,"Error: the node is not allocated");
         break;
       case _CENTRAL_NODE_NUMBER_SIGNAL_ :
         //pipe.recv(LocalCellNumber,From);
         LocalCellNumber=*((long int*)(buffer+offset));

         PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(LocalCellNumber,iCell,jCell,kCell);
         offset+=sizeof(long int);

         break;
       case _NEW_PARTICLE_SIGNAL_ :
         newParticle=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)],true);

         PIC::ParticleBuffer::UnPackParticleData(buffer+offset,newParticle);
         recvParticleCounter++;

         offset+=PIC::ParticleBuffer::ParticleDataLength;
         break;
       default:
         exit(__LINE__,__FILE__,"Error: the option is not recognized");
       }

       Signal=*((int*)(buffer+offset));
       offset+=sizeof(int);

       if ((Signal==_NEW_BLOCK_ID_SIGNAL_)||(Signal==_END_COMMUNICATION_SIGNAL_)) {
         memcpy(recvNode->block->FirstCellParticleTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
       }

     }

      //end the part of the receiver
     if (offset!=exchangeProcMatrix[From*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread]) exit(__LINE__,__FILE__,"Error: the amount of recieved data is not consistent");

     DataStillToRecieve--;
     delete [] recvDataBuffer[nproc];
     recvDataBuffer[nproc]=NULL;
  }

  //finish the send operations
  for (nproc=0;nproc<nSendProc;nproc++) {
    MPI_Wait(SendRequest+nproc,&status);
    delete [] sendDataBuffer[nproc];
  }
}

void PIC::Parallel::ProcessCornerBlockBoundaryNodes() {
  int thread,iThread,i,j,k,iface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  PIC::Mesh::cDataCornerNode *CornerNode;
  char *CornerNodeAssociatedData;
  PIC::Mesh::cDataBlockAMR *block;
  MPI_Status status;

  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};

  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};

  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};

  struct cStencilElement {
    int StencilLength;
    int StencilThreadTable[8];
    int iCornerNode,jCornerNode,kCornerNode;
    cAMRnodeID nodeid;
  };

  int iStencil;
  static int StencilTableLength=0;
  static cStencilElement *StencilTable=NULL;


  //generate a new stencil table
  static int nMeshModificationCounter=-1;

  if (nMeshModificationCounter!=PIC::Mesh::mesh.nMeshModificationCounter) {
    int NewTableLength=0;

    //update the coundater
    nMeshModificationCounter=PIC::Mesh::mesh.nMeshModificationCounter;

    //reset the 'processed' flag
    for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      block=node->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) {
        for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++)  for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)  {
          CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

          if (CornerNode!=NULL) CornerNode->SetProcessedFlag(false);
        }
      }
    }

    //determine the new length of the table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      int flag;

      block=node->block;

      for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++)  for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (block!=NULL) {
          flag=1;

          CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

          if (CornerNode!=NULL) {
            CornerNodeAssociatedData=CornerNode->GetAssociatedDataBufferPointer();
            if (CornerNode->TestProcessedFlag()==true) flag=0;
            CornerNode->SetProcessedFlag(true);
          }
          else flag=0;
        }
        else flag=0;

        //combine the array of flags
        int FlagTable[PIC::nTotalThreads],FlagSum;

        MPI_Allgather(&flag,1,MPI_INT,FlagTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
        for (thread=0,FlagSum=0;thread<PIC::nTotalThreads;thread++) FlagSum+=FlagTable[thread];

        if ((flag==1)&&(FlagSum!=1)) {
          NewTableLength++;
        }
      }
    }

    //allocate the new Stencile Table
    if (StencilTableLength!=0) delete [] StencilTable;

    StencilTableLength=NewTableLength;
    StencilTable=new cStencilElement[NewTableLength];

    //populate the Stencil Table
    //reset the 'processed' flag
    for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      block=node->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) {
        for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++)  for (k=kFaceMin[iface];k<=kFaceMax[iface];k++)  {
          CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

          if (CornerNode!=NULL) CornerNode->SetProcessedFlag(false);
        }
      }
    }

    //populate the table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      int flag;

      block=node->block;

      for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++)  for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (block!=NULL) {
          flag=1;

          CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

          if (CornerNode!=NULL) {
            CornerNodeAssociatedData=CornerNode->GetAssociatedDataBufferPointer();
            if (CornerNode->TestProcessedFlag()==true) flag=0;
            CornerNode->SetProcessedFlag(true);
          }
          else flag=0;
        }
        else flag=0;

        //combine the array of flags
        int FlagTable[PIC::nTotalThreads],FlagSum;

        MPI_Allgather(&flag,1,MPI_INT,FlagTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
        for (thread=0,FlagSum=0;thread<PIC::nTotalThreads;thread++) FlagSum+=FlagTable[thread];

        if ((flag==1)&&(FlagSum!=1)) {
          //the thread will participate in the data exchage
          StencilTable[iStencil].StencilLength=0;
          StencilTable[iStencil].iCornerNode=i,StencilTable[iStencil].jCornerNode=j,StencilTable[iStencil].kCornerNode=k;
          StencilTable[iStencil].nodeid=node->AMRnodeID;

          for (thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==1) StencilTable[iStencil].StencilThreadTable[StencilTable[iStencil].StencilLength++]=thread;

          iStencil++;
        }
      }
    }

  }


  //combine the 'corner' associated data vectors from the 'corner' nodes at the boundary of the blocks
  if (ProcessCornerNodeAssociatedData!=NULL) for (iStencil=0;iStencil<StencilTableLength;iStencil++) {
    node=PIC::Mesh::mesh.findAMRnodeWithID(StencilTable[iStencil].nodeid);
    CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(StencilTable[iStencil].iCornerNode,StencilTable[iStencil].jCornerNode,StencilTable[iStencil].kCornerNode));
    CornerNodeAssociatedData=CornerNode->GetAssociatedDataBufferPointer();

    char tempCornerNodeAssociatedData[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];
    char RecvBuffer[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

    if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[0]) {
      memcpy(tempCornerNodeAssociatedData,CornerNodeAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

      //process the associated data
      for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
        MPI_Recv(RecvBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_CHAR,StencilTable[iStencil].StencilThreadTable[iThread],0,MPI_GLOBAL_COMMUNICATOR,&status);
        ProcessCornerNodeAssociatedData(tempCornerNodeAssociatedData,RecvBuffer);
      }

      //save the processes data vector
      if (CopyCornerNodeAssociatedData!=NULL) CopyCornerNodeAssociatedData(CornerNodeAssociatedData,tempCornerNodeAssociatedData);
      else memcpy(CornerNodeAssociatedData,tempCornerNodeAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

      //send out the associated data
      for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
        MPI_Send(tempCornerNodeAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_CHAR,StencilTable[iStencil].StencilThreadTable[iThread],0,MPI_GLOBAL_COMMUNICATOR);
      }
    }
    else {
      //send the original associated data vector
      MPI_Send(CornerNodeAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_CHAR,StencilTable[iStencil].StencilThreadTable[0],0,MPI_GLOBAL_COMMUNICATOR);

      //recieve the associated data vector combined with that from other subdomain
      MPI_Recv(RecvBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_CHAR,StencilTable[iStencil].StencilThreadTable[0],0,MPI_GLOBAL_COMMUNICATOR,&status);


      if (CopyCornerNodeAssociatedData!=NULL) CopyCornerNodeAssociatedData(CornerNodeAssociatedData,RecvBuffer);
      else memcpy(CornerNodeAssociatedData,RecvBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
    }
  }
}

