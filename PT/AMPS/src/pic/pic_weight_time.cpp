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
  double ParticleInjection=0.0;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //injection rate from the boundariees of the box
//    ParticleInjection=GetMaximumBlockInjectionRate(spec);


    //injection rate from the internal surfaces
    list<cInternalBoundaryConditionsDescriptor>::iterator descriptor;
    cInternalSphericalData *Sphere;
    cInternalCircleData *Circle;
    cInternalSphere1DData *Sphere1D;

    for (descriptor=PIC::Mesh::mesh.InternalBoundaryList.begin();descriptor!=PIC::Mesh::mesh.InternalBoundaryList.end();descriptor++) {
      switch (descriptor->BondaryType) {
      case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
        if (DIM!=3) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 3D simulations");

        Sphere=(cInternalSphericalData*)descriptor->BoundaryElement;
        ParticleInjection+=Sphere->InjectionRate(spec,(void*)Sphere)*Sphere->maxIntersectedNodeTimeStep[spec];
        break;
      case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
        if (DIM!=2) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 2D simulations");

        Circle=(cInternalCircleData*)descriptor->BoundaryElement;
        ParticleInjection+=Circle->InjectionRate(spec,(void*)Circle)*Circle->maxIntersectedNodeTimeStep[spec];
        break;
      case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
        if (DIM!=1) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 1D simulations");

        Sphere1D=(cInternalSphere1DData*)descriptor->BoundaryElement;
        ParticleInjection+=Sphere1D->InjectionRate(spec,(void*)Sphere1D)*Sphere1D->maxIntersectedNodeTimeStep[spec];
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
      }
    }


    GlobalParticleWeight=ParticleInjection/maxReferenceInjectedParticleNumber;

    //exchange the particle weight
    double WeightArray[PIC::nTotalThreads];
    MPI_Gather(&GlobalParticleWeight,1,MPI_DOUBLE,WeightArray,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (PIC::ThisThread==0) for (int thread=0;thread<PIC::nTotalThreads;thread++) if (GlobalParticleWeight<WeightArray[thread]) GlobalParticleWeight=WeightArray[thread];
    MPI_Bcast(&GlobalParticleWeight,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if (GlobalParticleWeight<=0.0) exit(__LINE__,__FILE__,"Error: ParticleInjection has zero value");
  }



  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block!=NULL) startNode->block->SetLocalParticleWeight(GlobalParticleWeight,spec);
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

    if (startNode->block!=NULL) for (s=0;s<PIC::nTotalSpecies;s++) {
      blockTimeStep=PIC::ParticleWeightTimeStep::LocalTimeStep(s,startNode);
      startNode->block->SetLocalTimeStep(blockTimeStep,s);

      //injection rate from the internal surfaces
 //     list<cInternalBoundaryConditionsDescriptor>::iterator descriptor;

      cInternalBoundaryConditionsDescriptor *descriptor;
      cInternalSphericalData *Sphere;
      cInternalCircleData *Circle;
      cInternalSphere1DData *Sphere1D;

//      for (descriptor=PIC::Mesh::mesh.InternalBoundaryList.begin();descriptor!=PIC::Mesh::mesh.InternalBoundaryList.end();descriptor++) {
      for (descriptor=startNode->InternalBoundaryDescriptorList;descriptor!=NULL;descriptor=descriptor->nextInternalBCelement) {

        switch (descriptor->BondaryType) {
        case _INTERNAL_BOUNDARY_TYPE_SPHERE_:
          if (DIM!=3) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 3D simulations");

          Sphere=(cInternalSphericalData*)descriptor->BoundaryElement;
          if (Sphere->maxIntersectedNodeTimeStep[s]<blockTimeStep) Sphere->maxIntersectedNodeTimeStep[s]=blockTimeStep;
          break;
        case _INTERNAL_BOUNDARY_TYPE_CIRCLE_:
          if (DIM!=2) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 2D simulations");

          Circle=(cInternalCircleData*)descriptor->BoundaryElement;
          if (Circle->maxIntersectedNodeTimeStep[s]<blockTimeStep) Circle->maxIntersectedNodeTimeStep[s]=blockTimeStep;
          break;
        case _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_:
          if (DIM!=1) exit(__LINE__,__FILE__,"Error: cInternalSphericalData can be used ONLY for 1D simulations");

          Sphere1D=(cInternalSphere1DData*)descriptor->BoundaryElement;
          if (Sphere1D->maxIntersectedNodeTimeStep[s]<blockTimeStep) Sphere1D->maxIntersectedNodeTimeStep[s]=blockTimeStep;
          break;
        default:
          exit(__LINE__,__FILE__,"Error: the boundary type is not recognized");
        }
      }
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
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_

  class cCopyTimeStepDistribution {
  public:
    double dt;

    void CopyTimeStep(int s1,int s0,double t,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
      if (startNode->block!=NULL) {
        dt=startNode->block->GetLocalTimeStep(s0)*t;
        startNode->block->SetLocalTimeStep(dt,s1);
      }

      for (int nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) CopyTimeStep(s1,s0,t,startNode->downNode[nDownNode]);
    }
  } CopyTimeStepDistribution;

  CopyTimeStepDistribution.CopyTimeStep(specTarget,specSource,ProportionaltyCoefficient,PIC::Mesh::mesh.rootTree);
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

}
