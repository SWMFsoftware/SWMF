//$Id$
//The physical model of GCR at the boundary of the computational domain

/*
 * BoundaryInjection_GCR.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include "Earth.h"

//double **Earth::BoundingBoxInjection::GCR::BoundaryFaceLocalInjectionRate=NULL;
//double *Earth::BoundingBoxInjection::GCR::maxBoundaryFaceLocalInjectionRate=NULL;
//double *Earth::BoundingBoxInjection::GCR::BoundaryFaceTotalInjectionRate=NULL;
//cBoundaryFaceDescriptor *Earth::BoundingBoxInjection::BoundaryFaceDescriptor=NULL;
//int Earth::BoundingBoxInjection::nTotalBoundaryInjectionFaces=0;

//get the total injection rate
double Earth::BoundingBoxInjection::GCR::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

//  exit(__LINE__,__FILE__,"Error: not implemented");

  return 1.0;
}



//generate properties of a new particle
void Earth::BoundingBoxInjection::GCR::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {
  double v[3],Energy,mass,WeightCorrectionFactor;

  //power law energy distribution will be used for testing
  //later more realisting models of GCR (Nymmic, etc) will be added to the library
  //of GCR boundary box injection models

  const double EnergyDistributionPowerIndex=3;

  //Distribute energy
  Energy=minEnergy*(1.0+rnd()*(maxEnergy-minEnergy));
  WeightCorrectionFactor=pow(Energy,-EnergyDistributionPowerIndex)/pow(minEnergy,-EnergyDistributionPowerIndex);


  //generate velocity of the new particle
  PIC::Distribution::InjectRingDistribution(v,Energy,ExternalNormal,spec);

  //save parameters of the new particle
  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,ParticleData);

}

//init the SEP injection model
void Earth::BoundingBoxInjection::GCR::Init() {
/*

  //Init the particle injection tables
  Earth::BoundingBoxInjection::InitBoundingBoxInjectionTable(BoundaryFaceLocalInjectionRate,maxBoundaryFaceLocalInjectionRate,
      BoundaryFaceTotalInjectionRate,InjectionRate);
*/
}
