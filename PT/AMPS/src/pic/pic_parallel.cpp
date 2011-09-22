//====================================================
//$Id$
//====================================================
//the functions that control the interprocessor communication of the code

#include "pic.h"

long int PIC::Parallel::sendParticleCounter=0,PIC::Parallel::recvParticleCounter=0,PIC::Parallel::IterationNumberAfterRebalancing=0;
double PIC::Parallel::RebalancingTime=0.0,PIC::Parallel::CumulativeLatency=0.0;
double PIC::Parallel::EmergencyLoadRebalancingFactor=3.0;

//====================================================
//Exchange particles between Processors
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

