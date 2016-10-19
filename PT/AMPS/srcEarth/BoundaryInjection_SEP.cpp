//Physical model of SEP at the boundary of the computational domain
//$Id$

/*
 * BoundaryInjection_SEP.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include "Earth.h"

double Earth::BoundingBoxInjection::SEP::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  exit(__LINE__,__FILE__,"Error: not implemented");

  return 1.0;
}

void Earth::BoundingBoxInjection::SEP::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {

  exit(__LINE__,__FILE__,"Error: not implemented");
}


