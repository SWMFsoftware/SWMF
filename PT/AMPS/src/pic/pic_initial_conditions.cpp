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

long int PIC::InitialCondition::PutParticle(int spec, double *x, double *v){
  // the function places a single (for each processor) particle of spec species
  // to a given location with given velocity vector
  // this function is meant for DEBUG purposes
  int iCell,jCell,kCell,nInjectionProcess;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;

  node=PIC::Mesh::mesh.findTreeNode(x);

  //initiate the new particle
  PIC::ParticleBuffer::InitiateParticle(x,v,NULL,&spec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
  
  return 1;  
}

long int PIC::InitialCondition::PrepopulateDomain(int spec,double NumberDensity,double *Velocity,double Temperature) {
  int iCell,jCell,kCell,nInjectionProcess;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;

  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];

  //particle ejection parameters
  double ParticleWeight,beta=PIC::MolecularData::GetMass(spec)/(2*Kbol*Temperature);

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

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    xmin=node->xmin,xmax=node->xmax;

    //particle stat weight
    #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    #else
    ParticleWeight=node->block->GetLocalParticleWeight(spec);
    #endif

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

        //initiate the new particle
        PIC::ParticleBuffer::InitiateParticle(x,v,NULL,&spec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
      }
      //end of the particle injection block
    }
  }

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return nGlobalInjectedParticles;
}


