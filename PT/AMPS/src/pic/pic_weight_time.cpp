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
  static double ConstantWeightValue;

  if (startNode==PIC::Mesh::mesh.rootTree) ConstantWeightValue=GetMaximumBlockInjectionRate(spec)/maxReferenceInjectedParticleNumber;



  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    startNode->block->SetLocalParticleWeight(ConstantWeightValue,spec);
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
