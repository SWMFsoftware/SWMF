//====================================================
//$Id$
//====================================================
//the functions that controls execution of the boundary conditions

#include "pic.h"


//the list of blocks where the injection BCs are applied

PIC::BC::fBlockInjectionBC PIC::BC::BlockInjectionBCindicatior=NULL;
PIC::BC::fBlockInjectionBC PIC::BC::userDefinedBoundingBlockInjectionFunction=NULL;
list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > PIC::BC::boundingBoxInjectionBlocksList;

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

  //model the particle injection through the face of the bounding box
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;

  for (nodeptr=boundingBoxInjectionBlocksList.begin(),end=boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    if ((*nodeptr)->Thread==PIC::Mesh::mesh.ThisThread) userDefinedBoundingBlockInjectionFunction(*nodeptr);
  }

}



