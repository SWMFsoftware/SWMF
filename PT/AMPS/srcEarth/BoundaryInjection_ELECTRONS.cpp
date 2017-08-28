//Physical model of electrons at the boundary of the computational domain
//$Id$

#include "pic.h"
#include "Earth.h"

double Earth::BoundingBoxInjection::Electrons::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  return 1.0;
}

void Earth::BoundingBoxInjection::Electrons::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {

  double Energy=100.0*KeV2J;
  double v[3];

  //generate velocity of the new particle
  PIC::Distribution::InjectRingDistribution(v,Energy,ExternalNormal,spec);


  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::Mover::GuidingCenter::InitiateMagneticMoment(spec,x,v,ParticleData,startNode);
}

//init the Electrons injection model
void Earth::BoundingBoxInjection::Electrons::Init() {
  //do nothing for now
}
