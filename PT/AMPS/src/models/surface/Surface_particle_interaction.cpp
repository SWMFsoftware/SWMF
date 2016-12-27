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
  double BulkFlowVelocity[3]={0.0,0.0,0.0},temp,*e0,*e1,*e2,ksi_n_r,ksi_t1_r,ksi_t2_r,beta,phi,c;
  int spec;

  spec=PIC::ParticleBuffer::GetI(ptr);
  temp=Surface::GetSurfaceTemeprature(TriangleCutFace,startNode);
  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*temp));

  //get the internal frame of reference related to the surface element
  e0=TriangleCutFace->e0Orthogonal;
  e1=TriangleCutFace->e1Orthogonal;
  e2=TriangleCutFace->ExternalNormal;

  //the normal component of the reflected velocity
  ksi_n_r=sqrt(-log(rnd()))/beta;

  //tangenetial component of thereflected velocity
  phi=PiTimes2*rnd();
  c=sqrt(-log(rnd()))/beta;

  ksi_t1_r=c*cos(phi);
  ksi_t2_r=c*sin(phi);

  //combine the reflection velocity vector
  for (int idim=0;idim<3;idim++) vInit[idim]=ksi_t1_r*e0[idim]+ksi_t2_r*e1[idim]+ksi_n_r*e2[idim];

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

//CLL model. Implements from Padilla-2009-JTHT, Padilla-2008-PHD (as in MONACO)
int Surface::CLL::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double *e0,*e1,*e2;
  int idim,spec;
  double ksi_n_i=0.0,ksi_t1_i=0.0;
  double r1,phi2,ksi_n_m,ksi_mp_w,ksi_n_r;
  double r3,phi4,ksi_t_m,ksi_t1_r;
  double r5,phi6,ksi_t2_r;
  double alpha_n,alpha_t;

  //get the internal frame of reference related to the surface element
  e0=TriangleCutFace->e0Orthogonal;
  e1=TriangleCutFace->e1Orthogonal;
  e2=TriangleCutFace->ExternalNormal;

  //get thenorml and tangential components of the incident velocities
  for (idim=0;idim<3;idim++) ksi_n_i+=vInit[idim]*e2[idim],ksi_t1_i+=vInit[idim]*e1[idim];

  spec=PIC::ParticleBuffer::GetI(ptr);
  ksi_mp_w=sqrt(2.0*Kbol/PIC::MolecularData::GetMass(spec)*Surface::GetSurfaceTemeprature(TriangleCutFace,startNode));

  //get the normal component of the refrected velocity
  r1=sqrt(-alpha_n*log(rnd()));
  phi2=PiTimes2*rnd();
  ksi_n_m=fabs(ksi_n_i/ksi_mp_w)*sqrt(1.0-alpha_n);
  ksi_n_r=ksi_mp_w*sqrt(r1*r1+ksi_n_m*ksi_n_m+2.0*r1*ksi_n_m*cos(phi2));

  //get the first tangential component of the refrected velocity
  r3=sqrt(-alpha_t*log(rnd()));
  phi4=PiTimes2*rnd();
  ksi_t_m=fabs(ksi_t1_i/ksi_mp_w)*sqrt(1-alpha_t);
  ksi_t1_r=ksi_mp_w*(ksi_t_m+r3*cos(phi4));

  //get the second tangential component of the refrected velocity
  r5=sqrt(-alpha_t*log(rnd()));
  phi6=PiTimes2*rnd();
  ksi_t2_r=ksi_mp_w*r5*cos(phi6);

  //combine the reflection velocity vector
  for (idim=0;idim<3;idim++) vInit[idim]=ksi_t1_r*e0[idim]+ksi_t2_r*e1[idim]+ksi_n_r*e2[idim];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


