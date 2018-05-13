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
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData=NULL;
PIC::Parallel::CornerBlockBoundaryNodes::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData=NULL; 
bool PIC::Parallel::CornerBlockBoundaryNodes::ActiveFlag=false;
void PIC::Parallel::CornerBlockBoundaryNodes::SetActiveFlag(bool flag) {ActiveFlag=flag;}

//default function forprocessing of the corner node associated data
void PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
}

void PIC::Parallel::CornerBlockBoundaryNodes::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
}

//-----------------------------------------------------------------


/*
//processing 'corner' and 'center' node associated data vectors when perform syncronization
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::ProcessCenterNodeAssociatedData=NULL,PIC::Parallel::ProcessCornerNodeAssociatedData=NULL;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCenterNodeAssociatedData=PIC::Parallel::CopyCenterNodeAssociatedData_default;
PIC::Parallel::fUserDefinedProcessNodeAssociatedData PIC::Parallel::CopyCornerNodeAssociatedData=PIC::Parallel::CopyCornerNodeAssociatedData_default;

//default function forprocessing of the corner node associated data
void PIC::Parallel::CopyCornerNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
}

void PIC::Parallel::CopyCenterNodeAssociatedData_default(char *TargetBlockAssociatedData,char *SourceBlockAssociatedData) {
  memcpy(TargetBlockAssociatedData,SourceBlockAssociatedData,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
}
*/ 

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
      MPI_Isend(buffer,offset,MPI_BYTE,To,0,MPI_GLOBAL_COMMUNICATOR,SendRequest+nproc);
   }


  //initiate reciving of the data
  for (nproc=0;nproc<nRecvProc;nproc++) {
    thread=recvProcList[nproc];
    MPI_Irecv(recvDataBuffer[nproc],exchangeProcMatrix[thread*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread],MPI_BYTE,thread,0,MPI_GLOBAL_COMMUNICATOR,RecvRequest+nproc);
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

  if (CornerBlockBoundaryNodes::ActiveFlag==false) return;

  const int iFaceMin[6]={0,_BLOCK_CELLS_X_,0,0,0,0};
  const int iFaceMax[6]={0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_};

  const int jFaceMin[6]={0,0,0,_BLOCK_CELLS_Y_,0,0};
  const int jFaceMax[6]={_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_};

  const int kFaceMin[6]={0,0,0,0,0,_BLOCK_CELLS_Z_};
  const int kFaceMax[6]={_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,0,_BLOCK_CELLS_Z_};

  static const int BlockCornerTable[8][3]={
      {0,0,0}, {_BLOCK_CELLS_X_,0,0}, {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0}, {0,_BLOCK_CELLS_Y_,0},
      {0,0,_BLOCK_CELLS_Z_}, {_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_}, {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}, {0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}
  };

  static const int BlockCornerOffsetTable[8][3]={
      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_}
  };

  static const int BlockEdgeOffsetTable[12][3]={
      {0,-_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {0,+_BLOCK_CELLS_Y_,-_BLOCK_CELLS_Z_},
      {0,+_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},
      {0,-_BLOCK_CELLS_Y_,+_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,0,-_BLOCK_CELLS_Z_},
      {+_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},
      {-_BLOCK_CELLS_X_,0,+_BLOCK_CELLS_Z_},

      {-_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},
      {+_BLOCK_CELLS_X_,-_BLOCK_CELLS_Y_,0},
      {+_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0},
      {-_BLOCK_CELLS_X_,+_BLOCK_CELLS_Y_,0}
  };

  static const int BlockFaceOffsetTable[6][3]={
      {-_BLOCK_CELLS_X_,0,0},{+_BLOCK_CELLS_X_,0,0},
      {0,-_BLOCK_CELLS_Y_,0},{0,+_BLOCK_CELLS_Y_,0},
      {0,0,-_BLOCK_CELLS_Z_},{0,0,+_BLOCK_CELLS_Z_}
  };

  static const int iMinEdgeTable[12]={
      0,0,0,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0
  };

  static const int iMaxEdgeTable[12]={
      _BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0,
      0,_BLOCK_CELLS_X_,_BLOCK_CELLS_X_,0
  };

  static const int jMinEdgeTable[12]={
      0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,
      0,0,0,0,
      0,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_
  };

  static const int jMaxEdgeTable[12]={
      0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,0,
      _BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_,
      0,0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Y_
  };

  static const int kMinEdgeTable[12]={
    0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
    0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
    0,0,0,0
  };

  static const int kMaxEdgeTable[12]={
      0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
      0,0,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,
      _BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_,_BLOCK_CELLS_Z_
  };

  struct cStencilElement {
    int StencilLength;
    int StencilThreadTable[80];
    int iCornerNode,jCornerNode,kCornerNode;
    cAMRnodeID nodeid;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
    PIC::Mesh::cDataCornerNode *CornerNode;
    char *AssociatedDataPointer;
  };

  int iStencil;
  static int StencilTableLength=0;
  static cStencilElement *StencilTable=NULL;

  //SendTablePBC,RecvTablePBC,SendTableLengthPBC,RecvTableLengthPBC used for data exchange in case the periodic boundary conditions are enforsed
  struct cSendRecvCornerNodeDataElement {
    int i,j,k;
    cAMRnodeID SendBlockID,RecvBlockID;
    PIC::Mesh::cDataCornerNode *CornerNode;
    char *AssociatedDataPointer;
    int OperationID;
  };

  //pocesseing of the associated data of the nodes located at the boundary of the "real" computational domain in case the periodic boundary conditions are in use
  struct cNodeSetElement {
    cAMRnodeID nodeid;
    int i,j,k;
  };

  struct cNodeSet {
    cNodeSetElement NodeTable[80];
    int NodeTableLength;
  };

  list<cNodeSet> NodeSetList;


  struct cStencilPoint {
    int ThreadTable[80];
    cNodeSetElement NodeSetTable[80];
    int ThreadTableLength;

    char *AssociatedDataPointer;
  };

  struct cStencilPBC {
    cStencilPoint PointTable[80];
    int nStencilPoints;
    int ProcessingThread;
    int SourceThreadTable[80];

    //threads to which the processed associated vector will be send out back
    int InvolvedThreadTable[80],InvolvedThreadTableLength;
    bool InvolvedFlag;
  };

  struct cBlockTable {
    cAMRnodeID GhostBloks[80],RealBlockPair[80];
    int iCorner[80],jCorner[80],kCorner[80];
    int BlockTableLength;
  };

  static int StencilTablePBCLength=0;
  static cStencilPBC *StencilTablePBC=NULL;




  //generate a new stencil table
  static int nMeshModificationCounter=-1;


  //function to reset the corner node 'processed' flag
  auto ResetProcessedFlag = [&] () {
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      PIC::Mesh::cDataBlockAMR *block=node->block;

      if (block!=NULL) for (i=0;i<_BLOCK_CELLS_X_+1;i++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        PIC::Mesh::cDataCornerNode *CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

        if (CornerNode!=NULL) CornerNode->SetProcessedFlag(false);
      }
    }
  };



  auto VerifyAllocatedNeighbourBlock = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool FoundFlag=false;
    int iNeighbour;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if ((neibNode=node->neibNodeFace[iNeighbour])!=NULL) if (neibNode->block!=NULL)  {
      FoundFlag=true;
      break;
    }

    //check corners
    if (FoundFlag==false) for (iNeighbour=0;iNeighbour<8;iNeighbour++) if ((neibNode=node->neibNodeCorner[iNeighbour])!=NULL) if (neibNode->block!=NULL) {
      FoundFlag=true;
      break;
    }

    //check edges
    if (FoundFlag==false) for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if ((neibNode=node->neibNodeEdge[iNeighbour])!=NULL) if (neibNode->block!=NULL) {
      FoundFlag=true;
      break;
    }

    return FoundFlag;
  };


  //find a corner node in neighbour blocks
  auto FindCornerNbN = [&] (int &iNode,int &jNode,int &kNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &NeibNode,PIC::Mesh::cDataCornerNode *CornerNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool res=false;
    int i,j,k,iface,iNeighbour;
    PIC::Mesh::cDataBlockAMR *block;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if (node->neibNodeFace[iNeighbour]!=NULL) {
      block=node->neibNodeFace[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeFace[iNeighbour];
          return res;
        }
      }
    }

    //check corners
    for (iNeighbour=0;iNeighbour<8;iNeighbour++) if (node->neibNodeCorner[iNeighbour]!=NULL) {
      block=node->neibNodeCorner[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeCorner[iNeighbour];
          return res;
        }
      }
    }

    //check edges
    for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if (node->neibNodeEdge[iNeighbour]!=NULL) {
      block=node->neibNodeEdge[iNeighbour]->block;

      if (block!=NULL) for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
        if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
          res=true;
          iNode=i,jNode=j,kNode=k;
          NeibNode=node->neibNodeEdge[iNeighbour];
          return res;
        }
      }
    }

    NeibNode=NULL;
    return res;
  };

  auto GetCornerNode = [&] (int &iNode,int &jNode,int &kNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &NodeOut,cAMRnodeID NodeidIn) {
    PIC::Mesh::cDataCornerNode *res=NULL;
    int iface,iedge,icorner,i,j,k,iNeighbour;
    PIC::Mesh::cDataBlockAMR *block;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode,*node=PIC::Mesh::mesh.findAMRnodeWithID(NodeidIn);


    //check self
    if ((block=node->block)!=NULL) {
      res=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(iNode,jNode,kNode));
      NodeOut=node;
      return res;
    }
    else {
      //search connections through faces
      for (iface=0;iface<6;iface++) if ((neibNode=node->GetNeibFace(iface,0,0))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockFaceOffsetTable[iface][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockFaceOffsetTable[iface][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockFaceOffsetTable[iface][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res==neibNode->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

      //search connections through edges
      for (iedge=0;iedge<12;iedge++) if ((neibNode=node->GetNeibEdge(iedge,0))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockEdgeOffsetTable[iedge][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockEdgeOffsetTable[iedge][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockEdgeOffsetTable[iedge][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res==neibNode->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

      //search connections through corners
      for (icorner=0;icorner<8;icorner++) if ((neibNode=node->GetNeibCorner(icorner))!=NULL) if (neibNode->block!=NULL) {
        if (node->RefinmentLevel!=neibNode->RefinmentLevel) exit(__LINE__,__FILE__,"Error: not implemented for the case when neibours has different resolution levels");

        i=iNode-BlockCornerOffsetTable[icorner][0];

        if ((0<=i)&&(i<=_BLOCK_CELLS_X_)) {
          j=jNode-BlockCornerOffsetTable[icorner][1];

          if ((0<=j)&&(j<=_BLOCK_CELLS_Y_)) {
            k=kNode-BlockCornerOffsetTable[icorner][2];

            if ((0<=k)&&(k<=_BLOCK_CELLS_Z_)) {
              if ((res==neibNode->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                NodeOut=neibNode;
                iNode=i,jNode=j,kNode=k;
                return res;
              }
            }
          }
        }
      }

    }

    res=NULL;
    return res;
  };


  auto VerifyBoundaryCornerNode = [&] (PIC::Mesh::cDataCornerNode *CornerNode,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    bool res=false;
    int iNeighbour,i,j,k,iface;
    bool EnternalBlockFlag,StartBlockExternalFlag;
    PIC::Mesh::cDataBlockAMR *block;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode,*neibNode_LastProcessed=NULL;

    StartBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

    //check faces
    for (iNeighbour=0;iNeighbour<6*4;iNeighbour++) if ((neibNode=node->neibNodeFace[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    //check corners
    for (iNeighbour=0;iNeighbour<8;iNeighbour++) if ((neibNode=node->neibNodeCorner[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    //check edges
    for (iNeighbour=0;iNeighbour<12*2;iNeighbour++) if ((neibNode=node->neibNodeEdge[iNeighbour])!=NULL) if (neibNode!=neibNode_LastProcessed) if ((block=neibNode->block)!=NULL) {
      neibNode_LastProcessed=neibNode;
      EnternalBlockFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if ( ((EnternalBlockFlag==true)&&(StartBlockExternalFlag==false)) || ((EnternalBlockFlag==false)&&(StartBlockExternalFlag==true)) ) {
        for (iface=0;iface<6;iface++) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          if (CornerNode==block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))) {
            //the 'corner' node is at the boundary of the 'real' domain
            res=true;
            return res;
          }
        }
      }
    }

    return res;
  };

  auto TestCheckFace = [&] (int iface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;
    int i,j;

    //check face neighbors
    for (i=0;i<2;i++) for (j=0;j<2;j++) {
      if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) {
        if (node->Thread!=neibNode->Thread) {
         //the check need to be performed
         return true;
        }

        if (node->RefinmentLevel==neibNode->RefinmentLevel) break;
      }
    }

    //check edge neighbors
    static const int EdgeFaceTable[6][4]={{4,11,7,8},{5,10,6,9}, {0,9,3,8},{1,10,2,11}, {0,5,1,4},{3,6,2,7}};

    for (int iedge=0;iedge<4;iedge++) {
      for (i=0;i<2;i++) {
        if ((neibNode=node->GetNeibEdge(EdgeFaceTable[iface][iedge],i))!=NULL) {
          if (node->Thread!=neibNode->Thread) {
            //the check need to be performed
            return true;
          }

          if (node->RefinmentLevel==neibNode->RefinmentLevel) break;
        }
      }
    }

    //check corner neighbors
    static const int CornerFaceTable[6][4]={{0,3,7,4},{1,2,6,5}, {0,1,5,4},{3,2,6,7}, {0,1,2,3},{4,5,6,7}};

    for (i=0;i<4;i++) if ((neibNode=node->GetNeibCorner(CornerFaceTable[iface][i]))!=NULL) if (node->Thread!=neibNode->Thread) {
      //the check need to be performed
      return true;
    }

    return false;
  };

  auto TestRealBoundaryFace = [&] (int iface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    int i,j;
    PIC::Mesh::cDataCornerNode *CornerNode;
    bool StartBlockExternalFlag,NeibBlockExternalFlag;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;

    StartBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

    for (i=0;i<2;i++) for (j=0;j<2;j++) if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) {
      NeibBlockExternalFlag=(PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) ? true : false;

      if (((StartBlockExternalFlag==true)&&(NeibBlockExternalFlag==false)) || ((StartBlockExternalFlag==false)&&(NeibBlockExternalFlag==true))) {
        return true;
      }
      else {
        return false;
      }
    }

    return false;
  };


  //determine whether the mesh/domain decomposition have been changed
  int localMeshChangeFlag,globalMeshChangeFlag;

  localMeshChangeFlag=(nMeshModificationCounter==PIC::Mesh::mesh.nMeshModificationCounter) ? 0 : 1;
  MPI_Allreduce(&localMeshChangeFlag,&globalMeshChangeFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);


  if (globalMeshChangeFlag!=0) {
    //the mesh or the domain decomposition has been modified. Need to create a new communucation table
    int NewTableLength=0;
    int iNode,jNode,kNode,FlagSum;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ActualNode;
    bool meshModifiedFlag_CountMeshElements=PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements;
    double StartTime=MPI_Wtime();

    //reset the flags
    ResetProcessedFlag();

    //determine the new length of the table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      if ((_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_)||(PIC::Mesh::mesh.ExternalBoundaryBlock(node)!=_EXTERNAL_BOUNDARY_BLOCK_)) {
        int flag;
        bool TemporaralyAllocatedBlock=false;
        cAMRnodeID nodeid=node->AMRnodeID;

        //for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_+1;j++)  for (int k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        for (iface=0;iface<6;iface++) if (TestCheckFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          //determine whether the corner node has been already processes
          if (node->Thread==PIC::ThisThread) {
            CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
            flag=(CornerNode->TestProcessedFlag()==true) ? 0 : 1;
          }

          MPI_Bcast(&flag,1,MPI_INT,node->Thread,MPI_GLOBAL_COMMUNICATOR);

          if (flag==0) continue;

          //the corner node has not been processed yet
          iNode=i,jNode=j,kNode=k;

          if ((CornerNode=GetCornerNode(iNode,jNode,kNode,ActualNode,nodeid))!=NULL) {
            //the corner node exists at the current MPI process
            CornerNode->SetProcessedFlag(true);
            flag=1;
          }
          else flag=0;

          MPI_Allreduce(&flag,&FlagSum,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

          if (FlagSum>1) {
            NewTableLength++;
          }
        }
      }
    }

    //allocate the new Stencile Table
    if (StencilTableLength!=0) delete [] StencilTable;

    StencilTableLength=NewTableLength;
    StencilTable=new cStencilElement[NewTableLength];

    //reset the 'processed' flag
    ResetProcessedFlag();

    //populate the Stencil Table
    for (iStencil=0,node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
      if ((_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_)||(PIC::Mesh::mesh.ExternalBoundaryBlock(node)!=_EXTERNAL_BOUNDARY_BLOCK_)) {
        int flag;
        bool TemporaralyAllocatedBlock=false;
        cAMRnodeID nodeid=node->AMRnodeID;

        //for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_+1;j++)  for (int k=0;k<_BLOCK_CELLS_Z_+1;k++) {
        for (iface=0;iface<6;iface++) if (TestCheckFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
          //determine whether the corner node has been already processes
          if (node->Thread==PIC::ThisThread) {
            CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
            flag=(CornerNode->TestProcessedFlag()==true) ? 0 : 1;
          }

          MPI_Bcast(&flag,1,MPI_INT,node->Thread,MPI_GLOBAL_COMMUNICATOR);

          if (flag==0) continue;

          //the corner node has not been processed yet
          iNode=i,jNode=j,kNode=k;

          if ((CornerNode=GetCornerNode(iNode,jNode,kNode,ActualNode,nodeid))!=NULL) {
            //the corner node exists at the current MPI process
            CornerNode->SetProcessedFlag(true);
            flag=1;
          }
          else flag=0;

          MPI_Allreduce(&flag,&FlagSum,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

          if (FlagSum>1) {
            StencilTable[iStencil].StencilLength=0;

            StencilTable[iStencil].iCornerNode=iNode,StencilTable[iStencil].jCornerNode=jNode,StencilTable[iStencil].kCornerNode=kNode;
            StencilTable[iStencil].nodeid=ActualNode->AMRnodeID;
            StencilTable[iStencil].node=ActualNode;
            StencilTable[iStencil].CornerNode=CornerNode;
            StencilTable[iStencil].AssociatedDataPointer=(CornerNode!=NULL) ? CornerNode->GetAssociatedDataBufferPointer() : NULL;

            //combine information of which MPI processes has a copy of the current corner node
            int FlagTable[PIC::nTotalThreads];

            MPI_Allgather(&flag,1,MPI_INT,FlagTable,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

            for (thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==1) {
              StencilTable[iStencil].StencilThreadTable[StencilTable[iStencil].StencilLength++]=thread;
            }

            iStencil++;
          }
        }
      }
    }

    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      //In case when the periodic boundary conditions are inforced
      //additional information exchange table need to be generated to link points located at the boundary of the "real" domain
      list<cStencilPBC> StencilListPBC;
      cBlockTable BlockTable; //contains information of all blocks that containes a given corner node

      ResetProcessedFlag();

      for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode)  {
        if (PIC::Mesh::mesh.ExternalBoundaryBlock(node)==_EXTERNAL_BOUNDARY_BLOCK_)  {
          for (iface=0;iface<6;iface++) if (TestRealBoundaryFace(iface,node)==true) for (i=iFaceMin[iface];i<=iFaceMax[iface];i++) for (j=jFaceMin[iface];j<=jFaceMax[iface];j++) for (k=kFaceMin[iface];k<=kFaceMax[iface];k++) {
            //the analysis will be performed by the MPI process that the block belongs to
            int ii,jj,kk,ff,iNeighbour,BoundaryNodeFlag=false;
            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
            PIC::Mesh::cDataBlockAMR *neibBlock;
            bool found;

            if (node->Thread==PIC::ThisThread) {
              BlockTable.BlockTableLength=0;
              CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));

              //determine whether the point is at the boundary
              BoundaryNodeFlag=VerifyBoundaryCornerNode(CornerNode,node);

              if ((BoundaryNodeFlag==true)&&(CornerNode->TestProcessedFlag()==false)) {
                //the corner node has not been used yet. The set of the corner nodes that are connected due to enforcing the periodic boundary conditions is not determined yet

                BlockTable.iCorner[0]=i;
                BlockTable.jCorner[0]=j;
                BlockTable.kCorner[0]=k;
                BlockTable.GhostBloks[0]=node->AMRnodeID;
                BlockTable.BlockTableLength=1;


                //find all 'ghost' nodes that contains that corner node
                //check faces
                for (int iiface=0;iiface<6;iiface++) if ((neibNode=node->GetNeibFace(iiface,0,0))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                  neibBlock=neibNode->block;

                  if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                    exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                  }

                  for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                    found=true;
                    break;
                  }

                  if (found==true) continue;

                  if (neibBlock!=NULL) {
                    switch (iiface) {
                    case 0:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k;
                      break;
                    case 1:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k;
                      break;
                    case 2:
                      ii=i;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    case 3:
                      ii=i;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    case 4:
                      ii=i;
                      jj=j;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 5:
                      ii=i;
                      jj=j;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    }

                    if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)&&(0<=jj)&&(jj<=_BLOCK_CELLS_Y_)&&(0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                      if (CornerNode==neibBlock->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(ii,jj,kk))) {
                        //a new block has been found that contains tested 'corner' node
                        BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                        BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                        BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                        BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                        BlockTable.BlockTableLength++;
                      }
                    }
                  }
                }

                //check corners
                for (int icorner=0;icorner<8;icorner++) if ((neibNode=node->GetNeibCorner(icorner))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                  neibBlock=neibNode->block;
                    
                  if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                    exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                  }
                    
                  for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                    found=true;
                    break;
                  }
                      
                  if (found==true) continue;

                  if (neibBlock!=NULL) {
                    switch (icorner) {
                    case 0:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 1:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 2:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 3:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;

                    case 4:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    case 5:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    case 6:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    case 7:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    }
                          
                    if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)&&(0<=jj)&&(jj<=_BLOCK_CELLS_Y_)&&(0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                      if (CornerNode==neibBlock->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(ii,jj,kk))) {
                        //a new block has been found that contains tested 'corner' node
                        BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                        BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                        BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                        BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                        BlockTable.BlockTableLength++;
                      }
                    }
                  }
                }
                  

                //check edges
                for (int iedge=0;iedge<12;iedge++) if ((neibNode=node->GetNeibEdge(iedge,0))!=NULL) if (PIC::Mesh::mesh.ExternalBoundaryBlock(neibNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
                  neibBlock=neibNode->block;
                     
                  if (node->RefinmentLevel!=neibNode->RefinmentLevel) {
                    exit(__LINE__,__FILE__,"Error: blocks has different refinement levels");
                  }
                     
                  for (found=false,ii=0;ii<BlockTable.BlockTableLength;ii++) if (BlockTable.GhostBloks[ii]==neibNode->AMRnodeID) {
                    found=true;
                    break;
                  }
                      
                  if (found==true) continue;
                      
                  if (neibBlock!=NULL) {
                    switch (iedge) {
                    case 0:
                      ii=i;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 1:
                      ii=i;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 2:
                      ii=i;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    case 3:
                      ii=i;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;

                    case 4:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 5:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k+_BLOCK_CELLS_Z_;
                      break;
                    case 6:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;
                    case 7:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j;
                      kk=k-_BLOCK_CELLS_Z_;
                      break;

                    case 8:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    case 9:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j+_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    case 10:
                      ii=i-_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    case 11:
                      ii=i+_BLOCK_CELLS_X_;
                      jj=j-_BLOCK_CELLS_Y_;
                      kk=k;
                      break;
                    }


                    if ((0<=ii)&&(ii<=_BLOCK_CELLS_X_)&&(0<=jj)&&(jj<=_BLOCK_CELLS_Y_)&&(0<=kk)&&(kk<=_BLOCK_CELLS_Z_)) {
                      if (CornerNode==neibBlock->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(ii,jj,kk))) {
                        //a new block has been found that contains tested 'corner' node
                        BlockTable.iCorner[BlockTable.BlockTableLength]=ii;
                        BlockTable.jCorner[BlockTable.BlockTableLength]=jj;
                        BlockTable.kCorner[BlockTable.BlockTableLength]=kk;
                        BlockTable.GhostBloks[BlockTable.BlockTableLength]=neibNode->AMRnodeID;
                        BlockTable.BlockTableLength++;
                          
                        found=true;
                      }
                    }
                  }
                }
                  
                //determine 'real' blocks tha correspond to those in 'BlockTable'
                for (ii=0;ii<BlockTable.BlockTableLength;ii++) {
                  BlockTable.RealBlockPair[ii]=PIC::BC::ExternalBoundary::Periodic::findCorrespondingRealBlock(PIC::Mesh::mesh.findAMRnodeWithID(BlockTable.GhostBloks[ii]))->AMRnodeID;
                }
              }
            }

            //Broadcast the 'BlockTable' to other MPI processes
            int iBlock;
            int DataRequestFlagTable[PIC::nTotalThreads];
            cNodeSetElement Set;
            cStencilPBC NewStencilElementPBC;

            MPI_Bcast(&BlockTable,sizeof(cBlockTable),MPI_BYTE,node->Thread,MPI_GLOBAL_COMMUNICATOR);

            NewStencilElementPBC.nStencilPoints=0;

            //loop through the 'BlockTable'
            for (int ipass=0;ipass<2;ipass++) for (iBlock=0;iBlock<((ipass==0) ? 1: BlockTable.BlockTableLength);iBlock++) {
              switch (ipass) {
              case 0:
                Set.i=BlockTable.iCorner[iBlock];
                Set.j=BlockTable.jCorner[iBlock];
                Set.k=BlockTable.kCorner[iBlock];
                Set.nodeid=BlockTable.GhostBloks[iBlock];
                break;
              default:
                Set.i=BlockTable.iCorner[iBlock];
                Set.j=BlockTable.jCorner[iBlock];
                Set.k=BlockTable.kCorner[iBlock];
                Set.nodeid=BlockTable.RealBlockPair[iBlock];
              }

              //determine whether a corner node descrived by 'Set' exist in the current MPI process
              //if the corner node exists -> set the processed flag state "true"
              int CornerNodeExistFlag;
              int CornerNodeRecvFlagTable[PIC::nTotalThreads],CornerNodeSendFlagTable[PIC::nTotalThreads];
              PIC::Mesh::cDataCornerNode *CornerNode;
              cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* NodeOut;

              CornerNode=GetCornerNode(Set.i,Set.j,Set.k,NodeOut,Set.nodeid);

              if (CornerNode!=NULL) {
                //the corner exists in the currect MPI process
                if (CornerNode->TestProcessedFlag()==false) {
                  CornerNode->SetProcessedFlag(true);
                  CornerNodeExistFlag=true;

                  Set.nodeid=NodeOut->AMRnodeID;
                }
                else CornerNodeExistFlag=false;
              }
              else {
                CornerNodeExistFlag=false;
              }

              //gather information from all MPI processes at the root MPI process
              MPI_Gather(&CornerNodeExistFlag,1,MPI_INT,CornerNodeRecvFlagTable,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

              if (PIC::ThisThread==0) {
                cStencilPoint StencilPoint;

                StencilPoint.ThreadTableLength=0;

                if (CornerNodeRecvFlagTable[0]==true) {
                  //root thread has the corner node
                  StencilPoint.NodeSetTable[0]=Set;
                  StencilPoint.ThreadTable[0]=PIC::ThisThread;
                  StencilPoint.ThreadTableLength=1;
                }

                for (thread=1;thread<PIC::nTotalThreads;thread++) if (CornerNodeRecvFlagTable[thread]==true) {
                  MPI_Status status;

                  MPI_Recv(&Set,sizeof(cNodeSetElement),MPI_BYTE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
                  StencilPoint.NodeSetTable[StencilPoint.ThreadTableLength]=Set;
                  StencilPoint.ThreadTable[StencilPoint.ThreadTableLength]=thread;
                  StencilPoint.ThreadTableLength++;
                }

                if (StencilPoint.ThreadTableLength!=0) NewStencilElementPBC.PointTable[NewStencilElementPBC.nStencilPoints++]=StencilPoint;
              }
              else if (CornerNodeExistFlag==true) {
                MPI_Send(&Set,sizeof(cNodeSetElement),MPI_BYTE,0,0,MPI_GLOBAL_COMMUNICATOR);
              }
            }

            //create a new entry to the stencil table
            if ((PIC::ThisThread==0)&&(NewStencilElementPBC.nStencilPoints!=0)) {
              //determine thread the will process the data
              NewStencilElementPBC.ProcessingThread=NewStencilElementPBC.PointTable[0].ThreadTable[0];

              //determine thread that is the source of the data for each point of the stencil
              for (ii=0;ii<NewStencilElementPBC.nStencilPoints;ii++) NewStencilElementPBC.SourceThreadTable[ii]=NewStencilElementPBC.PointTable[ii].ThreadTable[0];

              //determine all threads that are involved into the communications
              NewStencilElementPBC.InvolvedThreadTableLength=0;

              for (ii=0;ii<NewStencilElementPBC.nStencilPoints;ii++) {
                for (jj=0;jj<NewStencilElementPBC.PointTable[ii].ThreadTableLength;jj++) {
                  bool found=false;

                  for (kk=0;kk<NewStencilElementPBC.InvolvedThreadTableLength;kk++) if (NewStencilElementPBC.PointTable[ii].ThreadTable[jj]==NewStencilElementPBC.InvolvedThreadTable[kk]) {
                    found=true;
                    break;
                  }

                  if (found==false) NewStencilElementPBC.InvolvedThreadTable[NewStencilElementPBC.InvolvedThreadTableLength++]=NewStencilElementPBC.PointTable[ii].ThreadTable[jj];
                }
              }

              //add the new stencil to the stencil list
              StencilListPBC.push_front(NewStencilElementPBC);
            }
          }
        }
      }


    //now StencilPBC is complete. Convert the list into an array with the communication rules and broadcasr it to other MPI processes
    if (StencilTablePBCLength!=0) {
      delete [] StencilTablePBC;
      StencilTablePBCLength=0;
    }

    if (PIC::ThisThread==0) {
      //create the stencil table
      int i;
      list<cStencilPBC>::iterator p;

      StencilTablePBCLength=StencilListPBC.size();
      StencilTablePBC=new cStencilPBC[StencilTablePBCLength];

      for (i=0,p=StencilListPBC.begin();i<StencilTablePBCLength;i++,p++) StencilTablePBC[i]=*p;

      MPI_Bcast(&StencilTablePBCLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);



      CMPI_channel pipe;

      pipe.init(5000000);
      pipe.openBcast(0);

      for (i=0;i<StencilTablePBCLength;i++) pipe.send(StencilTablePBC[i]);

      pipe.closeBcast();



 //     MPI_Bcast(StencilTablePBC,StencilTablePBCLength*sizeof(cStencilPBC),MPI_BYTE,0,MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Bcast(&StencilTablePBCLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

      StencilTablePBC=new cStencilPBC[StencilTablePBCLength];


      CMPI_channel pipe;

      pipe.init(5000000);
      pipe.openBcast(0);

      for (i=0;i<StencilTablePBCLength;i++) pipe.recv(StencilTablePBC[i],0);

      pipe.closeBcast();


//      MPI_Bcast(StencilTablePBC,StencilTablePBCLength*sizeof(cStencilPBC),MPI_BYTE,0,MPI_GLOBAL_COMMUNICATOR);
    }

    //prepare data that is needed for a 'fast' performing of the data exchange operation
    for (iStencil=0;iStencil<StencilTablePBCLength;iStencil++) {
      //determine whether the current MPI processes in involved into the communication
      StencilTablePBC[iStencil].InvolvedFlag=false;

      for (i=0;i<StencilTablePBC[iStencil].InvolvedThreadTableLength;i++) if (StencilTablePBC[iStencil].InvolvedThreadTable[i]==PIC::ThisThread) {
        StencilTablePBC[iStencil].InvolvedFlag=true;
        break;
      }

      //if the current thread is involved into the communication -> initialize pointer to the associated data vector
      if (StencilTablePBC[iStencil].InvolvedFlag==true) {
        //loop through all points in the stencil
        for (i=0;i<StencilTablePBC[iStencil].nStencilPoints;i++) {
          //loop through all threads that have this point
          for (j=0;j<StencilTablePBC[iStencil].PointTable[i].ThreadTableLength;j++) {
            if (StencilTablePBC[iStencil].PointTable[i].ThreadTable[j]==PIC::ThisThread) {
              //the current MPI process has point 'j' of the stencil
              cNodeSetElement Set=StencilTablePBC[iStencil].PointTable[i].NodeSetTable[j];
              cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.findAMRnodeWithID(Set.nodeid);


              if (node->block!=NULL) {
                StencilTablePBC[iStencil].PointTable[i].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(Set.i,Set.j,Set.k))->GetAssociatedDataBufferPointer();
              }
              else {
                //the block is not allocated => the corner associated data vectors are not definied
                exit(__LINE__,__FILE__,"Error: the block is not allocated");
              }
            }
          }
        }
      }

    }
  }

    //update the coundater
    nMeshModificationCounter=PIC::Mesh::mesh.nMeshModificationCounter;
    PIC::Mesh::mesh.meshModifiedFlag_CountMeshElements=meshModifiedFlag_CountMeshElements;
    PIC::Parallel::RebalancingTime+=MPI_Wtime()-StartTime;
  }


  //performe the data exchange session
  char recvAssociatedDataBuffer[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

  //1. combine 'corner' node data
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) for (iStencil=0;iStencil<StencilTableLength;iStencil++) {
    int iThread;

    if (StencilTable[iStencil].StencilLength>1) {
      //there are more that one MPI processes that contributed to the state vector of the corner node
      if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[0]) {
        //this thread will combine data from all other MPI processes
        //recieve the data
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
          MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

          if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
            PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer);
          }
          else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
        }

        //send out the state vector to MPI processes that have contributed to it
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) {
          MPI_Send(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[iThread],iStencil,MPI_GLOBAL_COMMUNICATOR);
        }
      }

      else {
        for (iThread=1;iThread<StencilTable[iStencil].StencilLength;iThread++) if (PIC::ThisThread==StencilTable[iStencil].StencilThreadTable[iThread]) {
          //this thread will contribute to the colelcted corner node associated data
          MPI_Send(StencilTable[iStencil].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR);

          //recieve the associated data
          MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTable[iStencil].StencilThreadTable[0],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

          if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
            PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer);
          }
          else{
            memcpy(StencilTable[iStencil].AssociatedDataPointer,recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
          }

          break;
        }
      }
    }
  }

  //2. processes the 'real' domain boundary in case periodic boundary conditions are in use
  if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) for (iStencil=0;iStencil<StencilTablePBCLength;iStencil++) {
     int iThread,thread,ipoint;

     if (StencilTablePBC[iStencil].InvolvedFlag==true) {
       char *AssociateDataVector=StencilTablePBC[iStencil].PointTable[0].AssociatedDataPointer;

       for (ipoint=1;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) { //ipoint=1 is correct because the processing node is selected such that point=0 is already accounted for
         if (StencilTablePBC[iStencil].ProcessingThread==PIC::ThisThread) {
           //the current thread is the processing manager of the stencil
           if (StencilTablePBC[iStencil].SourceThreadTable[ipoint]==PIC::ThisThread) {
             //the data are located with the same MPI process
             if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
               PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(AssociateDataVector,StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer);
             }
             else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
           }
           else {
             //the data need to be recieved before processing
             MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].SourceThreadTable[ipoint],iStencil,MPI_GLOBAL_COMMUNICATOR,&status);

             if (PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData!=NULL) {
               PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData(AssociateDataVector,recvAssociatedDataBuffer);
             }
             else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::ProcessCornerNodeAssociatedData is not defined");
           }

         }
         else if (StencilTablePBC[iStencil].SourceThreadTable[ipoint]==PIC::ThisThread) {
           //the MPI process is not "Processing" but serves as a source of the associated data vector for the 'ipoint' of the stencil
           MPI_Send(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].ProcessingThread,iStencil,MPI_GLOBAL_COMMUNICATOR);
         }
       }

       //processing of the data is complete -> send it out to the involved MPI processes
       if (StencilTablePBC[iStencil].ProcessingThread==PIC::ThisThread) {
         int cnt=0;

         //copy processes associated data vector into a temporary buffer
         if (PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData!=NULL) {
           PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(recvAssociatedDataBuffer,AssociateDataVector);
         }
         else exit(__LINE__,__FILE__,"Error: PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData is not defined");


         //copy/send thecontent of the processed associated data buffer to all MPI processes that are involved in the stencil
         for (ipoint=0;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) for (iThread=0;iThread<StencilTablePBC[iStencil].PointTable[ipoint].ThreadTableLength;iThread++) {
           thread=StencilTablePBC[iStencil].PointTable[ipoint].ThreadTable[iThread];

           if (thread==PIC::ThisThread) {
             //copy the processes associated data
             PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,recvAssociatedDataBuffer);
           }
           else {
             //send processes associated data vector
             MPI_Send(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,thread,0,MPI_GLOBAL_COMMUNICATOR);
           }
         }
       }
       else {
         //recieve the processes assoviated vactor
         for (ipoint=0;ipoint<StencilTablePBC[iStencil].nStencilPoints;ipoint++) for (iThread=0;iThread<StencilTablePBC[iStencil].PointTable[ipoint].ThreadTableLength;iThread++) {
           thread=StencilTablePBC[iStencil].PointTable[ipoint].ThreadTable[iThread];

           if (thread==PIC::ThisThread) {
             //recieve the processed associated data
             MPI_Status status;

             MPI_Recv(recvAssociatedDataBuffer,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,MPI_BYTE,StencilTablePBC[iStencil].ProcessingThread,0,MPI_GLOBAL_COMMUNICATOR,&status);
             PIC::Parallel::CornerBlockBoundaryNodes::CopyCornerNodeAssociatedData(StencilTablePBC[iStencil].PointTable[ipoint].AssociatedDataPointer,recvAssociatedDataBuffer);
           }
         }

       }
     }
  }
}

