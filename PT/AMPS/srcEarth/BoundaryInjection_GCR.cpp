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


double *Earth::BoundingBoxInjection::GCR::InjectionRateTable=NULL;
double *Earth::BoundingBoxInjection::GCR::maxEnergySpectrumValue=NULL;

//get the total injection rate
double Earth::BoundingBoxInjection::GCR::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

//  exit(__LINE__,__FILE__,"Error: not implemented");

  //init the injection rate tables if needed
  if (InjectionRateTable==NULL) {
    int s,iCompositionGroup,iInterval,i;
    double minEnergy,maxEnergy,minVelocity,maxVelocity,mass,de,f0,f1;

    const int nTotalIntegrationIntervals=1000;

    InjectionRateTable=new double [PIC::nTotalSpecies];
    maxEnergySpectrumValue=new double [PIC::nTotalSpecies];

    for (s=0;s<PIC::nTotalSpecies;s++) InjectionRateTable[s]=0.0,maxEnergySpectrumValue[s]=0.0;

    //integrate the energy spectrum
    for (s=0;s<PIC::nTotalSpecies;s++) {
      iCompositionGroup=Earth::CompositionGroupTableIndex[s];
      mass=PIC::MolecularData::GetMass(s);

      minVelocity=Earth::CompositionGroupTable[iCompositionGroup].GetMinVelocity(s);
      maxVelocity=Earth::CompositionGroupTable[iCompositionGroup].GetMaxVelocity(s);

      minEnergy=Relativistic::Speed2E(minVelocity,mass);
      maxEnergy=Relativistic::Speed2E(maxVelocity,mass);

      //integrate the energy spectrum, and determine its maximum value in this energy interval
      de=(maxEnergy-minEnergy)/nTotalIntegrationIntervals;
      f0=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(minEnergy);
      maxEnergySpectrumValue[s]=f0;

      for (i=0;i<nTotalIntegrationIntervals;i++) {
        f1=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(minEnergy+de*(i+1));
        InjectionRateTable[s]+=(f0+f1)*de/2.0;

        if (f1>maxEnergySpectrumValue[s]) maxEnergySpectrumValue[s]=f1;

        f0=f1;
      }
    }
  }

  return InjectionRateTable[spec];
}



//generate properties of a new particle
void Earth::BoundingBoxInjection::GCR::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {
  double v[3],mass,WeightCorrectionFactor;
  double minV,maxV,minE,maxE,E;
  int CompositionGraoupIndex;

  mass=PIC::MolecularData::GetMass(spec);

  //get velocity limits
  CompositionGraoupIndex=Earth::CompositionGroupTableIndex[spec];
  minV=Earth::CompositionGroupTable[CompositionGraoupIndex].GetMinVelocity(spec);
  maxV=Earth::CompositionGroupTable[CompositionGraoupIndex].GetMaxVelocity(spec);

  //convert velocity into energy and distribute energy of a new particles
  minE=Relativistic::Speed2E(minV,mass);
  maxE=Relativistic::Speed2E(maxV,mass);

  do {
    E=minE+rnd()*(maxE-minE);
    WeightCorrectionFactor=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(E)/Earth::BoundingBoxInjection::GCR::maxEnergySpectrumValue[spec];
  }
  while (rnd()>WeightCorrectionFactor);

  WeightCorrectionFactor=1.0;


  //generate velocity of the new particle
  PIC::Distribution::InjectRingDistribution(v,E,ExternalNormal,spec);

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
