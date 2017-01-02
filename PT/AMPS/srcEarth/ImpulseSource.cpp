//$Id$
//functionality of the impulse source of the energetic particles in the magnetosphere

/*
 * ImpulseSource.cpp
 *
 *  Created on: Dec 26, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "Earth.h"

Earth::ImpulseSource::cImpulseSourceData Earth::ImpulseSource::ImpulseSourceData[]={0,0.0,false,{0.0,0.0,0.0},0.0,0.0};
int Earth::ImpulseSource::nTotalSourceLocations=0;
double Earth::ImpulseSource::TimeCounter=0.0;
double Earth::ImpulseSource::EnergySpectrum::Constant::e=1.0*MeV2J;
int Earth::ImpulseSource::EnergySpectrum::Mode=Earth::ImpulseSource::EnergySpectrum::Mode_Constatant;
bool Earth::ImpulseSource::Mode=false;

//inject the energetic particles
long int Earth::ImpulseSource::InjectParticles() {
  int nTotalInjectedParticles=0;
  int idim,spec,iTotalSourceLocations;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  double mass,a,v[3];
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  for (iTotalSourceLocations=0;iTotalSourceLocations<nTotalSourceLocations;iTotalSourceLocations++) {
    if ((TimeCounter>=ImpulseSourceData[iTotalSourceLocations].time)&&(ImpulseSourceData[iTotalSourceLocations].ProcessedFlag==false)) {
      ImpulseSourceData[iTotalSourceLocations].ProcessedFlag=true;

      //inject energetic particles
      startNode=PIC::Mesh::mesh.findTreeNode(ImpulseSourceData[iTotalSourceLocations].x);
      if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

      //determine the number of the injected particles
      spec=ImpulseSourceData[iTotalSourceLocations].spec;
      a=ImpulseSourceData[iTotalSourceLocations].Source/startNode->block->GetLocalParticleWeight(spec);
      nTotalInjectedParticles=(int)a;
      a-=nTotalInjectedParticles;
      if (rnd()<a) nTotalInjectedParticles++;

      mass=PIC::MolecularData::GetMass(spec);

      for (int iPart=0;iPart<nTotalInjectedParticles;iPart++) {
        //generate new particle velocity
        double v[3],speed;

        switch (EnergySpectrum::Mode) {
        case EnergySpectrum::Mode_Constatant:
          speed=Relativistic::E2Speed(EnergySpectrum::Constant::e,mass);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is not found");
        }

        Vector3D::Distribution::Uniform(v);
        for (idim=0;idim<3;idim++) v[idim]*=speed;

        //generate particles' velocity
        newParticle=PIC::ParticleBuffer::GetNewParticle();
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

        PIC::ParticleBuffer::SetX(ImpulseSourceData[iTotalSourceLocations].x,newParticleData);
        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetI(spec,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

        //inject the particle into the system
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,TimeCounter-ImpulseSourceData[iTotalSourceLocations].time,startNode);
      }

    }

    //increment the species dependent time counter

  }

}

//set weights of the model particles
void Earth::ImpulseSource::InitParticleWeight() {
  int spec,iTotalSourceLocations;
  double WeightTable[PIC::nTotalSpecies],TotalSourceTable[PIC::nTotalSpecies],TotalRequestedModelParticlesTable[PIC::nTotalSpecies];
  double MinParticleWeight=-1.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) WeightTable[spec]=-1.0,TotalSourceTable[spec]=0.0,TotalRequestedModelParticlesTable[spec]=0.0;

  //loop through all injection events and collect the injection data
  for (iTotalSourceLocations=0;iTotalSourceLocations<nTotalSourceLocations;iTotalSourceLocations++) {
    TotalSourceTable[ImpulseSourceData[iTotalSourceLocations].spec]+=ImpulseSourceData[iTotalSourceLocations].Source;
    TotalRequestedModelParticlesTable[ImpulseSourceData[iTotalSourceLocations].spec]+=ImpulseSourceData[iTotalSourceLocations].NumberInjectedParticles;
  }

  //calculate the weight
  for (spec=0;spec<PIC::nTotalSpecies;spec++) if (TotalRequestedModelParticlesTable[spec]>0.0) {
    WeightTable[spec]=TotalSourceTable[spec]/TotalRequestedModelParticlesTable[spec];
    if ((MinParticleWeight<0.0)||(WeightTable[spec]<MinParticleWeight)) MinParticleWeight=WeightTable[spec];
  }

  //1) check is weights os all model spaces are defined. if not -> used the min weight
  //2) set the particle weight in the system
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (WeightTable[spec]<0.0) WeightTable[spec]=MinParticleWeight;
    PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(spec,WeightTable[spec]);
  }
}

