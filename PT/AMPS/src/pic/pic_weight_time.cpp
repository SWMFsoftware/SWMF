//====================================================
//$Id$
//====================================================
//the time step and particle weight adaptation procedures



#include "pic.h"

//set particle wight and local time step
PIC::ParticleWeightTimeStep::fSetFunction PIC::ParticleWeightTimeStep::LocalParticleWeight=NULL,PIC::ParticleWeightTimeStep::LocalTimeStep=NULL;
PIC::ParticleWeightTimeStep::fSetFunction PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=NULL;
double PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=400;

//the global particle weight/time step
double *PIC::ParticleWeightTimeStep::GlobalParticleWeight=NULL,*PIC::ParticleWeightTimeStep::GlobalTimeStep=NULL;


//====================================================
//get the maximum Block injection rate across the computational domain
double PIC::ParticleWeightTimeStep::GetMaximumBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double res=0.0;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    res=LocalBlockInjectionRate(spec,startNode)*startNode->block->GetLocalTimeStep(spec);
  }
  else {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;
    double c;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
       c=GetMaximumBlockInjectionRate(spec,downNode);
       if (res<c) res=c;
    }
  }

  return res;
}


//====================================================
//set particle's local weight (the weight is constant  across the domain)
void PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static double GlobalParticleWeight=0.0;
  double GlobalTimeStep,ParticleInjectionRate;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //injection rate from the boundariees of the box
    ParticleInjectionRate=GetMaximumBlockInjectionRate(spec);


    //injection rate from the internal surfaces
    list<cInternalBoundaryConditionsDescriptor>::iterator descriptor;
    cInternalSphericalData *Sphere;

    for (descriptor=PIC::Mesh::mesh.InternalBoundaryList.begin();descriptor!=PIC::Mesh::mesh.InternalBoundaryList.end();descriptor++) {
      switch (descriptor->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
        Sphere=(cInternalSphericalData*)descriptor->BoundaryElement;
        ParticleInjectionRate+=Sphere->InjectionRate(spec,(void*)Sphere);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
      }
    }

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    GlobalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
exit(__LINE__,__FILE__,"The mode is not implemented");
#endif

    GlobalParticleWeight=ParticleInjectionRate*GlobalTimeStep/maxReferenceInjectedParticleNumber;
  }



  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    startNode->block->SetLocalParticleWeight(GlobalParticleWeight,spec);
  }
  else {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) initParticleWeight_ConstantWeight(spec,downNode);
  }
}

void PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight() {
  int s;

  for (s=0;s<PIC::nTotalSpecies;s++) initParticleWeight_ConstantWeight(s);
}


//====================================================
//set particle's local time step

void PIC::ParticleWeightTimeStep::initTimeStep(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double blockTimeStep=0.0;
    int s;

    for (s=0;s<PIC::nTotalSpecies;s++) {
      blockTimeStep=PIC::ParticleWeightTimeStep::LocalTimeStep(s,startNode);
      startNode->block->SetLocalTimeStep(blockTimeStep,s);
    }

  }
  else {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) initTimeStep(downNode);
  }
}

//====================================================
//copy local time step and weight distribution from one spece to another
void PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(int specTarget,int specSource,double ProportionaltyCoefficient) {
#if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  GlobalParticleWeight[specTarget]=ProportionaltyCoefficient*GlobalParticleWeight[specSource];
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif
}

void PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(int specTarget,int specSource,double ProportionaltyCoefficient) {
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  GlobalTimeStep[specTarget]=ProportionaltyCoefficient*GlobalTimeStep[specSource];
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

}
