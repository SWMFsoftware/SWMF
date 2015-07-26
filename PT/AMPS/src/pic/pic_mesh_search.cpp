
//$Id$
//funstion that accelerates block/cell search on the AMPS's mesh


/*
 * pic_mesh_search.cpp
 *
 *  Created on: Jul 25, 2015
 *      Author: vtenishe
 */


#include "pic.h"

cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ***PIC::Mesh::Search::HashTable=NULL;
double PIC::Mesh::Search::dx[3],PIC::Mesh::Search::xMinDomain[3];

//======================================================================================================
//init the hash table
void PIC::Mesh::Search::ScanTree(int imin,int jmin,int kmin,int IndexInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k;


  if (IndexInterval==1) HashTable[imin][jmin][kmin]=node;
  else {
    if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      //the node is the last in the tree structure -> save it in the table
      for (i=0;i<IndexInterval;i++) for (j=0;j<IndexInterval;j++) for (k=0;k<IndexInterval;k++) HashTable[i+imin][j+jmin][k+kmin]=node;
    }
    else {
      //move down along the tree
      for (i=0;i<2;i++) for (j=0;j<2;j++) for (k=0;k<2;k++) {
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* DownNode=node->GetDownNode(i,j,k);

        if (DownNode!=NULL) ScanTree(imin+i*IndexInterval/2,jmin+j*IndexInterval/2,kmin+k*IndexInterval/2,IndexInterval/2,DownNode);
      }
    }
  }

}


void PIC::Mesh::Search::Init() {
  int i,j,k,cnt;
  int nElements=1<<HashTableLevel;

  //allocate the table
  if (HashTable==NULL) {
    HashTable=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* **[nElements];

    HashTable[0]=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* *[nElements*nElements];
    for (i=1;i<nElements;i++) HashTable[i]=HashTable[i-1]+nElements;

    HashTable[0][0]=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nElements*nElements*nElements];

    for (cnt=0,i=0;i<nElements;i++) for (j=0;j<nElements;j++) {
      HashTable[i][j]=HashTable[0][0]+cnt;
      cnt+=nElements;
    }
  }

  //reset the table
  for (i=0;i<nElements;i++) for (j=0;j<nElements;j++) for (k=0;k<nElements;k++) HashTable[i][j][k]=NULL;

  //reset the geometry information
  memcpy(xMinDomain,PIC::Mesh::mesh.rootTree->xmin,3*sizeof(double));
  for (i=0;i<3;i++) dx[i]=(PIC::Mesh::mesh.rootTree->xmax[i]-PIC::Mesh::mesh.rootTree->xmin[i])/nElements;


  //init the table with the pointers to the blocks
  ScanTree(0,0,0,nElements,PIC::Mesh::mesh.rootTree);
}

//======================================================================================================
//Find block
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::Mesh::Search::FindBlock(double *x) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* StartBlock;

  if (HashTable==NULL) Init();

  //determine the start block from the hash table
  int i,j,k;
  int nElements=1<<HashTableLevel;

  i=(x[0]-xMinDomain[0])/dx[0];

  if ((i<0)||(i>nElements-1)) StartBlock=NULL;
  else {
    j=(x[1]-xMinDomain[1])/dx[1];

    if ((j<0)||(j>nElements-1)) StartBlock=NULL;
    else {
      k=(x[2]-xMinDomain[2])/dx[2];

      if ((k<0)||(k>nElements-1)) StartBlock=NULL;
      else
        StartBlock=HashTable[i][j][k];
    }
  }

  //determine the actual block using procedures from the mesh obsect
  return PIC::Mesh::mesh.findTreeNode(x,StartBlock);
}

//======================================================================================================
//find cell
PIC::Mesh::cDataCenterNode *PIC::Mesh::Search::FindCell(double *x) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  int i,j,k,nd;

  node=FindBlock(x);
  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);

  if (node->block==NULL) return NULL;
  return node->block->GetCenterNode(nd);
}





