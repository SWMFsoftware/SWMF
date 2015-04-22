/*
 * pic_initial_conditions.cpp
 *
 *  Created on: Apr 17, 2015
 *      Author: vtenishe
 */

//$Id$
//the functionality of setting the initial distribution of the particles within the computational domain

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "pic.h"

long int PIC::InitialCondition::PrepopulateDomain(int spec,double NumberDensity,double *Velocity,double Temperature) {
  int iCell,jCell,kCell,nInjectionProcess;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;

  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];

  //particle ejection parameters
  double beta=PIC::MolecularData::GetMass(spec)/(2*Kbol*Temperature);

  //particle stat weight
#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  double ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

  //the boundaries of the block and middle point of the cell
  double *xmin,*xmax,*xMiddle;
  double x[3],v[3],anpart;
  int npart,idim;
  long int ptr,FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];


  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    xmin=node->xmin,xmax=node->xmax;
    memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      cell=cellList[nd];
      xMiddle=cell->GetX();

      //inject particles into the cell
      anpart=NumberDensity*cell->Measure/ParticleWeight;
      npart=(int)(anpart);
      if (rnd()<anpart-npart) npart++;
      nLocalInjectedParticles+=npart;

      while (npart-->0) {
        x[0]=xMiddle[0]+(xmax[0]-xmin[0])/_BLOCK_CELLS_X_*(rnd()-0.5);
        x[1]=xMiddle[1]+(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_*(rnd()-0.5);
        x[2]=xMiddle[2]+(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_*(rnd()-0.5);

        for (idim=0;idim<3;idim++) v[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())/beta)+Velocity[idim];

        //init the new particle
        ptr=PIC::ParticleBuffer::GetNewParticle();
        PIC::ParticleBuffer::SetX(x,ptr);
        PIC::ParticleBuffer::SetV(v,ptr);
        PIC::ParticleBuffer::SetI(spec,ptr);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,ptr);

        //apply the particle tracking condition
        #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
        PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

        PIC::ParticleTracker::InitParticleID(ParticleData);
        PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,ParticleData);
        #endif

        //add the paticle to the cell's particle list
        long int FirstCellParticle=FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

        PIC::ParticleBuffer::SetNext(FirstCellParticle,ptr);
        PIC::ParticleBuffer::SetPrev(-1,ptr);

        if (FirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstCellParticle);
        FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=ptr;
      }
      //end of the particle injection block

    }

    //update the particle list table
    memcpy(node->block->FirstCellParticleTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
  }

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return nGlobalInjectedParticles;
}


