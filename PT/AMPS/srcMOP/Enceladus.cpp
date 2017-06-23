//$Id$
//source of volatiles in the Jovian sytem


/*
 * Enceladus.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "MOP.h"

double MOP::SaturninanSystem::Enceladus::TotalSourceRate=1.0E23;
double MOP::SaturninanSystem::Enceladus::SourceTemperature=300.0;

double MOP::SaturninanSystem::Enceladus::vEnceladus[3]={0.0,0.0,0.0};
double MOP::SaturninanSystem::Enceladus::xEnceladus[3]={0.0,0.0,0.0};

//the volatile source
double MOP::SaturninanSystem::Enceladus::SourceRate(int spec) {
  double res=0.0;

  switch (spec) {
  case _H2O_SPEC_:
    res=TotalSourceRate;
    break;
  }

  return res;
}


long int MOP::SaturninanSystem::Enceladus::InjectParticles() {
  int nTotalInjectedParticles=0;
  int idim,spec;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  double mass,a,v[3],TimeCounter,TimeIncrement,ModelParticleInjectionRate,LocalTimeStep,ExternalNormal[3],InjectedFlowBulkVelocity[3]={0.0,0.0,0.0};
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  //determine location of Enceladus
  SpiceDouble lt,et,StateEnceladus[6];
  utc2et_c(Exosphere::SimulationStartTimeString,&et);

  et+=PIC::SimulationTime::Get();
  spkezr_c("Enceladus",et,"SSO","none","Saturn",StateEnceladus,&lt);

  for (idim=0;idim<3;idim++) {
    xEnceladus[idim]=1.0E3*StateEnceladus[idim];
    vEnceladus[idim]=1.0E3*StateEnceladus[3+idim];
  }

  //only water will be injected at this stage
  spec=_H2O_SPEC_;

  //etermine the current location of the source
  //inject energetic particles
  startNode=PIC::Mesh::mesh.findTreeNode(xEnceladus);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return nTotalInjectedParticles;


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel
  {
  #pragma omp single
  {
#endif

  //the number of the OpenMP threads
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int nThreadsOpenMP=omp_get_num_threads();
  #else
  int nThreadsOpenMP=1;
  #endif

  ModelParticleInjectionRate=SourceRate(_H2O_SPEC_)/startNode->block->GetLocalParticleWeight(spec);
  LocalTimeStep=startNode->block->GetLocalTimeStep(spec);

  if (ModelParticleInjectionRate>0.0) {
    TimeCounter=0.0;
    TimeIncrement=-log(rnd())/ModelParticleInjectionRate *rnd(); //<- *rnd() is to account for the injection of the first particle in the curent interaction

    while (TimeCounter+TimeIncrement<LocalTimeStep) {
      TimeCounter+=TimeIncrement;
      TimeIncrement=-log(rnd())/ModelParticleInjectionRate;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      #pragma omp task default (none) firstprivate (newParticle) private (idim,newParticleData)  \
        shared (SourceLocationB,GyroFrequencySample,GyroRadiiSample,SampleCounter,TimeCounter,iSource,nTotalInjectedParticles,startNode,spec,mass,ElectricCharge,EnergySpectrum::Mode,EnergySpectrum::Mode_Constatant,EnergySpectrum::Constant::e,ImpulseSourceData)
        {
#endif

        //generate the initial particle velocity
        Vector3D::Distribution::Uniform(ExternalNormal);
        PIC::Distribution::InjectMaxwellianDistribution(v,InjectedFlowBulkVelocity,SourceTemperature,ExternalNormal,spec,-1);
        for (idim=0;idim<3;idim++) v[idim]+=vEnceladus[idim];

        newParticle=PIC::ParticleBuffer::GetNewParticle();
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

        PIC::ParticleBuffer::SetX(xEnceladus,newParticleData);
        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetI(spec,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

        //apply the partile tracking condition if needed
        //apply condition of tracking the particle
        #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
        PIC::ParticleTracker::InitParticleID(newParticleData);
        PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xEnceladus,v,spec,newParticleData,(void*)startNode);
        #endif

        //inject the particle into the system
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        }
#endif

      nTotalInjectedParticles++;
    }

    //account for the particle productions rate
    PIC::BC::nTotalInjectedParticles+=nTotalInjectedParticles;
    PIC::BC::nInjectedParticles[spec]+=nTotalInjectedParticles;

    PIC::BC::ParticleProductionRate[spec]+=nTotalInjectedParticles/LocalTimeStep;
    PIC::BC::ParticleMassProductionRate[spec]+=nTotalInjectedParticles/LocalTimeStep*PIC::MolecularData::GetMass(spec);
  }
}
