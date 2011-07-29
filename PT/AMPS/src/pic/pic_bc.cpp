//====================================================
//$Id$
//====================================================
//the functions that controls execution of the boundary conditions

#include "pic.h"


//the list of blocks where the injection BCs are applied

PIC::BC::fBlockInjectionIndicator PIC::BC::BlockInjectionBCindicatior=NULL;
PIC::BC::fBlockInjectionBC PIC::BC::userDefinedBoundingBlockInjectionFunction=NULL;
list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > PIC::BC::boundingBoxInjectionBlocksList;
long int PIC::BC::nInjectedParticles=0;

//====================================================
//create the list of blocks where the injection BCs are applied
void PIC::BC::InitBoundingBoxInjectionBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  if (startNode==PIC::Mesh::mesh.rootTree) {
    if (boundingBoxInjectionBlocksList.size()!=0) exit(__LINE__,__FILE__,"Error: reinitialization of the 'boundingBoxInjectionBlocksList' list");
    if (BlockInjectionBCindicatior==NULL) exit(__LINE__,__FILE__,"Error: 'BlockInjectionBCindicatior' is not defined");
  }


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (BlockInjectionBCindicatior(startNode)==true) boundingBoxInjectionBlocksList.push_back(startNode);
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (int i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) InitBoundingBoxInjectionBlockList(downNode);
  }
}
//====================================================
//the function controls the overall execution of the injection boundary conditions
void PIC::BC::InjectionBoundaryConditions() {

  nInjectedParticles=0;

  //model the particle injection through the face of the bounding box
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  for (nodeptr=boundingBoxInjectionBlocksList.begin(),end=boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (node->Thread==PIC::Mesh::mesh.ThisThread) {
      nInjectedParticles+=userDefinedBoundingBlockInjectionFunction(node);

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
      EndTime=MPI_Wtime();
      node->ParallelLoadMeasure+=EndTime-StartTime;
      StartTime=EndTime;
#endif
    }
  }

}



