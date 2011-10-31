//====================================================
//$Id$
//====================================================
//the functions that control the interprocessor communication of the code

#include "pic.h"

long int PIC::Parallel::sendParticleCounter=0,PIC::Parallel::recvParticleCounter=0,PIC::Parallel::IterationNumberAfterRebalancing=0;
double PIC::Parallel::RebalancingTime=0.0,PIC::Parallel::CumulativeLatency=0.0;
double PIC::Parallel::EmergencyLoadRebalancingFactor=3.0;
double PIC::Parallel::Latency=0.0;

//====================================================
//Exchange particles between Processors
/*
void PIC::Parallel::ExchangeParticleData() {
  int From,To,pipeLastRecvThread=-1;
  long int Particle,NextParticle,newParticle,LocalCellNumber=-1;
  CMPI_channel pipe(100000);

  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *sendNode=NULL,*recvNode=NULL;


#if DIM == 3
  cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
#elif DIM == 2
  cMeshAMR2d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
#else
  cMeshAMR1d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
#endif


  pipe.openSend(0);
  pipe.openRecv(0);
  pipeLastRecvThread=0;

  sendParticleCounter=0,recvParticleCounter=0;

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


  //the data exchange loop
  for (From=0;From<PIC::Mesh::mesh.nTotalThreads;From++) for (To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((From!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[From][To]==true)) {

    //the part of the sender
    if (PIC::ThisThread==From) {
      bool CommunicationInitialed_BLOCK_;
      int iCell,jCell,kCell;

      //redirect the send pipe buffers
      pipe.RedirectSendBuffer(To);

      //reset the proceesed flaf for the blocks to be send
      //send the nodes' data
      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;

        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
          Particle=sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;


          if  (Particle!=-1) {
            if (CommunicationInitialed_BLOCK_==false) {
              PIC::Mesh::mesh.GetAMRnodeID(nodeid,sendNode);
              pipe.send(_NEW_BLOCK_ID_SIGNAL_);
              pipe.send((char*)(&nodeid),sizeof(nodeid));
              CommunicationInitialed_BLOCK_=true;
            }


            pipe.send(_CENTRAL_NODE_NUMBER_SIGNAL_);
            pipe.send(LocalCellNumber);

            while (Particle!=-1) {
              PIC::ParticleBuffer::PackParticleData(buffer,Particle);
              pipe.send(_NEW_PARTICLE_SIGNAL_);
              pipe.send(buffer,PIC::ParticleBuffer::ParticleDataLength);
              sendParticleCounter++;

              NextParticle=PIC::ParticleBuffer::GetNext(Particle);
              PIC::ParticleBuffer::DeleteParticle(Particle);
              Particle=NextParticle;
            }

            sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle=-1;
          }
        }
      }

      pipe.send(_END_COMMUNICATION_SIGNAL_);
      pipe.flush();
     //end the part of the sender
   }
   else if (PIC::ThisThread==To) {
     //the part of the receiver

     //redirect the recv's pipe buffers
     pipe.RedirectRecvBuffer(From);
     pipeLastRecvThread=From;
     pipe.recv(Signal,From);

     while (Signal!=_END_COMMUNICATION_SIGNAL_) {

       switch (Signal) {
       case _NEW_BLOCK_ID_SIGNAL_ :
         pipe.recv((char*)(&nodeid),sizeof(nodeid),From);
         recvNode=PIC::Mesh::mesh.findAMRnodeWithID(nodeid);

         if (recvNode->block==NULL) exit(__LINE__,__FILE__,"Error: the node is not allocated");
         break;
       case _CENTRAL_NODE_NUMBER_SIGNAL_ :
         pipe.recv(LocalCellNumber,From);
         break;
       case _NEW_PARTICLE_SIGNAL_ :
         pipe.recv(buffer,PIC::ParticleBuffer::ParticleDataLength,From);
         recvParticleCounter++;

         newParticle=PIC::ParticleBuffer::GetNewParticle(recvNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle);
         PIC::ParticleBuffer::UnPackParticleData(buffer,newParticle);
         break;
       default:
         exit(__LINE__,__FILE__,"Error: the option is not recognized");
       }

       pipe.recv(Signal,From);
     }

      //end the part of the receiver
    }

  }

  delete [] buffer;

  pipe.closeSend();
  pipe.closeRecv(pipeLastRecvThread);
}
*/

void PIC::Parallel::ExchangeParticleData() {
  int From,To;
  long int Particle,NextParticle,newParticle,LocalCellNumber=-1;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *sendNode=NULL,*recvNode=NULL;


#if DIM == 3
  cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID nodeid;
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

  //calculate the number of bytes that will be send
  for (To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
      bool CommunicationInitialed_BLOCK_;
      int iCell,jCell,kCell;

      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;

        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
          Particle=sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;


          if  (Particle!=-1) {
            if (CommunicationInitialed_BLOCK_==false) {
              //pipe.send(_NEW_BLOCK_ID_SIGNAL_);
              //pipe.send((char*)(&nodeid),sizeof(nodeid));
              sendProcVector[To]+=sizeof(nodeid)+sizeof(int);

              CommunicationInitialed_BLOCK_=true;
            }


            //pipe.send(_CENTRAL_NODE_NUMBER_SIGNAL_);
            //pipe.send(LocalCellNumber);
            sendProcVector[To]+=sizeof(int)+sizeof(LocalCellNumber);

            while (Particle!=-1) {
              //pipe.send(_NEW_PARTICLE_SIGNAL_);
              //pipe.send(buffer,PIC::ParticleBuffer::ParticleDataLength);
              sendProcVector[To]+=sizeof(int)+PIC::ParticleBuffer::ParticleDataLength;

              Particle=PIC::ParticleBuffer::GetNext(Particle);
            }


          }
        }
      }


      //end the part of the sender
      //pipe.send(_END_COMMUNICATION_SIGNAL_);
      if (sendProcVector[To]!=0) {
        sendProcVector[To]+=sizeof(int);
        sendProcList[nSendProc++]=To;
      }
   }


  //collect the data exchenge matrix
  Latency=MPI_Wtime();
  MPI_Allgather(sendProcVector,PIC::Mesh::mesh.nTotalThreads,MPI_INT,exchangeProcMatrix,PIC::Mesh::mesh.nTotalThreads,MPI_INT,MPI_COMM_WORLD);
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

      To=sendProcList[nproc];
      offset=0;
      buffer=sendDataBuffer[nproc];

      //reset the proceesed flaf for the blocks to be send
      //send the nodes' data
      for (sendNode=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];sendNode!=NULL;sendNode=sendNode->nextNodeThisThread) {
        CommunicationInitialed_BLOCK_=false;

        for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
          Particle=sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;


          if  (Particle!=-1) {
            if (CommunicationInitialed_BLOCK_==false) {
              PIC::Mesh::mesh.GetAMRnodeID(nodeid,sendNode);

              //pipe.send(_NEW_BLOCK_ID_SIGNAL_);
              *((int*)(buffer+offset))=_NEW_BLOCK_ID_SIGNAL_;
              offset+=sizeof(int);

              //pipe.send((char*)(&nodeid),sizeof(nodeid));
              #if DIM == 3
              *((cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID*)(buffer+offset))=nodeid;
              #else
              exit(__LINE__,__FILE__,"Error: not implemetned");
              #endif

              offset+=sizeof(nodeid);

              CommunicationInitialed_BLOCK_=true;
            }


            //pipe.send(_CENTRAL_NODE_NUMBER_SIGNAL_);
            *((int*)(buffer+offset))=_CENTRAL_NODE_NUMBER_SIGNAL_;
            offset+=sizeof(int);

            //pipe.send(LocalCellNumber);
            *((long int*)(buffer+offset))=LocalCellNumber;
            offset+=sizeof(long int);

            while (Particle!=-1) {
              //pipe.send(_NEW_PARTICLE_SIGNAL_);
              *((int*)(buffer+offset))=_NEW_PARTICLE_SIGNAL_;
              offset+=sizeof(int);

              //pipe.send(buffer,PIC::ParticleBuffer::ParticleDataLength);
              PIC::ParticleBuffer::PackParticleData(buffer+offset,Particle);
              offset+=PIC::ParticleBuffer::ParticleDataLength;
              sendParticleCounter++;

              NextParticle=PIC::ParticleBuffer::GetNext(Particle);
              PIC::ParticleBuffer::DeleteParticle(Particle);
              Particle=NextParticle;
            }

            sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle=-1;
          }
        }
      }

      //pipe.send(_END_COMMUNICATION_SIGNAL_);
      *((int*)(buffer+offset))=_END_COMMUNICATION_SIGNAL_;
      offset+=sizeof(int);

      if ((offset!=sendProcVector[To])||(offset!=exchangeProcMatrix[PIC::Mesh::mesh.ThisThread*PIC::Mesh::mesh.nTotalThreads+To])) exit(__LINE__,__FILE__,"Error: the data anount to be send is not consistent");

      //end the part of the sender - initiate non-blocks send
      MPI_Isend(buffer,offset,MPI_CHAR,To,0,MPI_COMM_WORLD,SendRequest+nproc);
   }


  //initiate reciving of the data
  for (nproc=0;nproc<nRecvProc;nproc++) {
    thread=recvProcList[nproc];
    MPI_Irecv(recvDataBuffer[nproc],exchangeProcMatrix[thread*PIC::Mesh::mesh.nTotalThreads+PIC::Mesh::mesh.ThisThread],MPI_CHAR,thread,0,MPI_COMM_WORLD,RecvRequest+nproc);
  }


  //recieve particles
  int DataStillToRecieve=nRecvProc;
  MPI_Status status;
  int flag;

  while (DataStillToRecieve!=0) {

    //determine the processor that has finished sending the data
    do {
      for (nproc=0;nproc<nRecvProc;nproc++) if (recvDataBuffer[nproc]!=NULL) {
        MPI_Test(RecvRequest+nproc,&flag,&status);
        if (flag==true) break;
      }
    }
    while (flag==false);

    buffer=recvDataBuffer[nproc];
    offset=0;
    From=recvProcList[nproc];

    //recieve the data
    //pipe.recv(Signal,From);
    Signal=*((int*)(buffer+offset));
    offset+=sizeof(int);

     while (Signal!=_END_COMMUNICATION_SIGNAL_) {

       switch (Signal) {
       case _NEW_BLOCK_ID_SIGNAL_ :
         //pipe.recv((char*)(&nodeid),sizeof(nodeid),From);
         #if DIM == 3
         nodeid=*((cMeshAMR3d<PIC::Mesh::cDataCornerNode,PIC::Mesh::cDataCenterNode,PIC::Mesh::cDataBlockAMR > :: cAMRnodeID*)(buffer+offset));
         #else
         exit(__LINE__,__FILE__,"Error: not implemetned");
         #endif

         offset+=sizeof(nodeid);
         recvNode=PIC::Mesh::mesh.findAMRnodeWithID(nodeid);

         if (recvNode->block==NULL) exit(__LINE__,__FILE__,"Error: the node is not allocated");
         break;
       case _CENTRAL_NODE_NUMBER_SIGNAL_ :
         //pipe.recv(LocalCellNumber,From);
         LocalCellNumber=*((long int*)(buffer+offset));
         offset+=sizeof(long int);

         break;
       case _NEW_PARTICLE_SIGNAL_ :
         //pipe.recv(buffer,PIC::ParticleBuffer::ParticleDataLength,From);

         newParticle=PIC::ParticleBuffer::GetNewParticle(recvNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle);
         PIC::ParticleBuffer::UnPackParticleData(buffer+offset,newParticle);
         recvParticleCounter++;

         offset+=PIC::ParticleBuffer::ParticleDataLength;
         break;
       default:
         exit(__LINE__,__FILE__,"Error: the option is not recognized");
       }

//       pipe.recv(Signal,From);
       Signal=*((int*)(buffer+offset));
       offset+=sizeof(int);
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

