//$Id$

/*
 * Surface_particle_interaction.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 */

#include "Surface.h"

//the accommodation coefficient table
double Surface::MaxwellReflection::AccommodationCoefficient[PIC::nTotalSpecies];


//particle reflection in the model of the specular reflection
int Surface::SpecularReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double c;

  c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];
  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//diffuse reflection
int Surface::DiffuseReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double BulkFlowVelocity[3]={0.0,0.0,0.0},temp,*ExternalNormal;
  int spec;

  spec=PIC::ParticleBuffer::GetI(ptr);
  temp=Surface::GetSurfaceTemeprature(TriangleCutFace,startNode);
  ExternalNormal=TriangleCutFace->ExternalNormal;

  PIC::Distribution::InjectMaxwellianDistribution(vInit,BulkFlowVelocity,temp,ExternalNormal,spec,_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_);

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//Maxwell model
int Surface::MaxwellReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int res,spec;

  spec=PIC::ParticleBuffer::GetI(ptr);

  if (rnd()<AccommodationCoefficient[spec]) res=Surface::DiffuseReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
  else res=Surface::SpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);

  return res;
}
