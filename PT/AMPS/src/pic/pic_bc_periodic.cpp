//$Id$
//functionality of the periodic boundary condition manager


/*
 * pic_bc_periodic.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: vtenishe
 */


#include "pic.h"

//originaly requsted limits of the computational domain
double PIC::BC::ExternalBoundary::Periodic::xminOriginal[3]={0.0,0.0,0.0};
double PIC::BC::ExternalBoundary::Periodic::xmaxOriginal[3]={0.0,0.0,0.0};

//the actual limits of the computational part (without "ghost" block) of the dimain
double PIC::BC::ExternalBoundary::Periodic::xminDomain[3]={0.0,0.0,0.0};
double PIC::BC::ExternalBoundary::Periodic::xmaxDomain[3]={0.0,0.0,0.0};

//the period length in each dimention
double PIC::BC::ExternalBoundary::Periodic::L[3]={0.0,0.0,0.0};

//the table of the correspondence of the 'ghost' and 'real' blocks
PIC::BC::ExternalBoundary::Periodic::cBlockPairTable *PIC::BC::ExternalBoundary::Periodic::BlockPairTable=NULL;
int PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength=0;

//thecommunication channel for message exchange between MPI processes
CMPI_channel PIC::BC::ExternalBoundary::Periodic::pipe;

//update data between the 'real' and 'ghost' blocks
void PIC::BC::ExternalBoundary::Periodic::UpdateData() {
  int i,j,k,iBlockPair,RealBlockThread,GhostBlockThread;

  //loop through all block pairs
  for (iBlockPair=0;iBlockPair<BlockPairTableLength;iBlockPair++) {
    if ( ((GhostBlockThread=BlockPairTable[iBlockPair].GhostBlock->Thread)==PIC::ThisThread) || ((RealBlockThread=BlockPairTable[iBlockPair].RealBlock->Thread)==PIC::ThisThread) ) {
      if (GhostBlockThread==RealBlockThread) {
        ExchangeBlockDataLocal(BlockPairTable[iBlockPair]);
      }
      else {
        ExchangeBlockDataMPI(BlockPairTable[iBlockPair]);
      }
    }
  }
}

//process the data update for the 'ghost' block
void PIC::BC::ExternalBoundary::Periodic::ExchangeBlockDataLocal(cBlockPairTable& BlockPair) {
  int i,j,k,idim;
  long int ptr,NextPtr;
  double *x;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *GhostBlock=BlockPair.GhostBlock;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *RealBlock=BlockPair.RealBlock;

  //attach particle list from the 'ghost' block to the 'real block'
  for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) if ((ptr=GhostBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)])!=-1) {
    //find the last particle in the block
    NextPtr=PIC::ParticleBuffer::GetNext(ptr);

    if (NextPtr!=-1) {
      do {
        ptr=NextPtr;

        //shift location of the particle
        x=PIC::ParticleBuffer::GetX(ptr);
        for (idim=0;idim<3;idim++) x[idim]=(x[idim]+L[idim]<xmaxOriginal[idim]) ? x[idim]+L[idim] : x[idim]-L[idim];

        NextPtr=PIC::ParticleBuffer::GetNext(ptr);
      }
      while (NextPtr!=-1);
    }
    else {
      //shift location of the particle
      x=PIC::ParticleBuffer::GetX(ptr);
      for (idim=0;idim<3;idim++) x[idim]=(x[idim]+L[idim]<xmaxOriginal[idim]) ? x[idim]+L[idim] : x[idim]-L[idim];
    }

    //reconnect the lists
    PIC::ParticleBuffer::SetNext(RealBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)],ptr);

    if (RealBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]!=-1) {
      PIC::ParticleBuffer::SetPrev(ptr,RealBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
    }

    RealBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=GhostBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    GhostBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
  }
}


void PIC::BC::ExternalBoundary::Periodic::ExchangeBlockDataMPI(cBlockPairTable& BlockPair) {
  int i,j,k;
  long int ptr,NewParticle,NextPtr;

  int GhostBlockThread=BlockPair.GhostBlock->Thread;
  int RealBlockThread=BlockPair.RealBlock->Thread;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *GhostBlock=BlockPair.GhostBlock;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *RealBlock=BlockPair.RealBlock;

  pipe.openSend(0);
  pipe.openRecv(0);
  int pipeLastRecvThread=0;

  //the list of signals
  const int NewCellCoordinatesSignal=0;
  const int NewParticelDataSignal=1;
  const int ParticleSendCompleteSignal=2;
  int Signal;

  //send the particles associated with the block
  if (GhostBlockThread==PIC::ThisThread) {
    //send out all particles
    //redirect the send pipe buffers
    pipe.RedirectSendBuffer(RealBlockThread);

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) if ((ptr=GhostBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)])!=-1) {
      pipe.send(NewCellCoordinatesSignal);
      pipe.send(i);
      pipe.send(j);
      pipe.send(k);

      while (ptr!=-1) {
        pipe.send(NewParticelDataSignal);
        pipe.send((char*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);


        NextPtr=PIC::ParticleBuffer::GetNext(ptr);
        PIC::ParticleBuffer::DeleteParticle(ptr);
        ptr=NextPtr;
      }
   
      GhostBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
    }

    pipe.send(ParticleSendCompleteSignal);
  }
  else {
    //recieve particles
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
    double *x;

    pipe.RedirectRecvBuffer(GhostBlockThread);
    Signal=pipe.recv<int>(GhostBlockThread);

    if (Signal!=ParticleSendCompleteSignal) {
      do {
        switch (Signal) {
        case NewCellCoordinatesSignal:
          i=pipe.recv<int>(GhostBlockThread);
          j=pipe.recv<int>(GhostBlockThread);
          k=pipe.recv<int>(GhostBlockThread);
          break;

        case NewParticelDataSignal:
          pipe.recv(tempParticleData,PIC::ParticleBuffer::ParticleDataLength,GhostBlockThread);

          //shift location of the particle
          x=PIC::ParticleBuffer::GetX((PIC::ParticleBuffer::byte *)tempParticleData);
          for (int idim=0;idim<3;idim++) x[idim]=(x[idim]+L[idim]<xmaxOriginal[idim]) ? x[idim]+L[idim] : x[idim]-L[idim];

          //generate a new particle
          NewParticle=PIC::ParticleBuffer::GetNewParticle(RealBlock->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
          PIC::ParticleBuffer::CloneParticle(PIC::ParticleBuffer::GetParticleDataPointer(NewParticle),(PIC::ParticleBuffer::byte *)tempParticleData);
          break;

        default:
          exit(__LINE__,__FILE__,"Error: unknown signal");
        }

        Signal=pipe.recv<int>(GhostBlockThread);
      }
      while (Signal!=ParticleSendCompleteSignal);
    }
  }

  //flush the pipe
  pipe.flush();
}

//Initialize the periodic boundary manager
void PIC::BC::ExternalBoundary::Periodic::Init(double* xmin,double* xmax,double (*localResuestedResolutionFunction)(double*)) {
  int idim;

  //init the pipe for message exchange between MPI processes
  pipe.init(1000000);

  //save the initial and modified domain limits
  for (idim=0;idim<3;idim++) {
    xminOriginal[idim]=xmin[idim],xmaxOriginal[idim]=xmax[idim];
    L[idim]=xmax[idim]-xmin[idim];
  }

  //initiate the mesh

}










