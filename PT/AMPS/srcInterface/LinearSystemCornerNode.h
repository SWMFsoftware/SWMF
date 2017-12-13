//$Id$
//interface to SWMF's GMRES solver

/*
 * LinearSystemCornerNode.h
 *
 *  Created on: Nov 20, 2017
 *      Author: vtenishe
 */

#ifndef _LINEARSYSTEMCORNERNODE_H_
#define _LINEARSYSTEMCORNERNODE_H_

#include "pic.h"
#include "linear_solver_wrapper_c.h"

class cLinearSystemCornerNodeDataRequestListElement {
public:
  int CornerNodeID;
  cAMRnodeID NodeID;
};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
class cLinearSystemCornerNode {
public:
  double *SubdomainPartialRHS,*SubdomainPartialUnknownsVector;
  int SubdomainPartialUnknownsVectorLength;

  class cStencilElement {
  public:
    cCornerNode* CornerNode;
    int CornerNodeID;
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    int UnknownVectorIndex,iVar,Thread;

    cStencilElement() {
      node=NULL,CornerNode=NULL;
      MatrixElementValue=0.0;
      UnknownVectorIndex=-1,iVar=-1,Thread=-1,CornerNodeID=-1;
    }
  };

  class cRhsSupportTable {
  public:
    double Coefficient;
    char *CornerNodeAssociatedDataPointer;
  };

  class cMatrixRow {
  public:
    cMatrixRow* next;
    double Rhs;

    cStencilElement Elements[MaxStencilLength];
    int nNonZeroElements;

    //pointed to the beginitg of the associated data vectros used in calculating of the right-hand side vectors
    cRhsSupportTable RhsSupportTable[MaxRhsSupportLength];
    int RhsSupportLength;

    int i,j,k;
    int iVar;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    cCornerNode* CornerNode;

    cMatrixRow() {
      next=NULL,Rhs=0.0,node=NULL,CornerNode=NULL;
      i=0,j=0,k=0,iVar=0,nNonZeroElements=0,RhsSupportLength=0;
    }
  };

  cStack<cMatrixRow> MatrixRowStack;

  //Table of the rows local to the current MPI process
  cMatrixRow *MatrixRowTable,*MatrixRowLast;

  //exchange buffers
  double** RecvExchangeBuffer;  //the actual recv data buffer
  int* RecvExchangeBufferLength;  //the number of elementf in the recv list

  double** SendExchangeBuffer;
  int* SendExchangeBufferLength;
  int **SendExchangeBufferElementIndex;

  //add new row to the matrix
  //for that a user function is called that returns stenal (elements of the matrix), corresponding rhs, and the number of the elements in the amtrix' line
  class cMatrixRowNonZeroElementTable {
  public:
    int i,j,k; //coordintes of the corner block used in the equation
    int iVar; //the index of the variable used in the stencil
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node;

    cMatrixRowNonZeroElementTable() {
      i=0,j=0,k=0,iVar=0,MatrixElementValue=0.0,Node=NULL;
    }
  };

  cMatrixRowNonZeroElementTable MatrixRowNonZeroElementTable[MaxStencilLength];

  //provcess the right boundary
  int nVirtualBlocks,nRealBlocks;
  void ProcessRightDomainBoundary(int *RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));


  //constructor
  cLinearSystemCornerNode() {
    MatrixRowTable=NULL,MatrixRowLast=NULL;

    //exchange buffers
    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL;
    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    //partial data of the linear system to be solved
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
    SubdomainPartialUnknownsVectorLength=0;

    nVirtualBlocks=0,nRealBlocks=0;
  }

  //void reset the data of the obsect to the default state
  void Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
  void Reset();


  //reset indexing of the nodes
  void ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //build the matrix
  void BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));

  //exchange the data
  void ExchangeIntermediateUnknownsData(double *x);

  //matrix/vector multiplication
  void MultiplyVector(double *p,double *x,int length);

  //call the linear system solver, and unpack the solution afterward
  void Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode));

  //update the RHS vector
  void UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int));

  //destructor
  ~cLinearSystemCornerNode() {
    Reset();

    MatrixRowStack.~cStack();
  }


};

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=node->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (node->downNode[i]!=NULL) ResetUnknownVectorIndex(node->downNode[i]);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::ProcessRightDomainBoundary(int* RecvDataPointCounter,void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  int i,j,k,iface,iblock;
  cMatrixRow* NewRow;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  cStencilElement* el;
//  int Index=PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;
  int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

//  int Debug=0;

  struct cOffset {
    int di,dj,dk;
  };


  nVirtualBlocks=0;

//  return ;

  cOffset xOffset[1]={{_BLOCK_CELLS_X_,0,0}};
  cOffset yOffset[1]={{0,_BLOCK_CELLS_Y_,0}};
  cOffset zOffset[1]={{0,0,_BLOCK_CELLS_Z_}};

  cOffset xyOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0}};
  cOffset xzOffset[3]={{_BLOCK_CELLS_X_,0,0},{0,0,_BLOCK_CELLS_Z_},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_}};
  cOffset yzOffset[3]={{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  cOffset xyzOffset[7]={{_BLOCK_CELLS_X_,0,0},{0,_BLOCK_CELLS_Y_,0},{0,0,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,0},{_BLOCK_CELLS_X_,0,_BLOCK_CELLS_Z_},{0,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_},
      {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_}};

  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    bool flag=false;

    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    if (node->Thread==PIC::ThisThread) for (iface=1;iface<6;iface+=2) if (node->GetNeibFace(iface,0,0)==NULL) {
      flag=true;
      break;
    }

    if (flag==false) continue;

/*Debug++;

//debug -> count the number of rows and blocks
int ntotbl=0;
for ( cMatrixRow* Row=MatrixRowTable;Row!=NULL;Row=Row->next) ntotbl++;*/

    //the block has at least one face at the 'right' boundary of the domain
    //determine the combination of the faces that are at the boundary
    int xBoundaryFace=0,yBoundaryFace=0,zBoundaryFace=0;

    if (node->GetNeibFace(1,0,0)==NULL) xBoundaryFace=10;
    if (node->GetNeibFace(3,0,0)==NULL) yBoundaryFace=100;
    if (node->GetNeibFace(5,0,0)==NULL) zBoundaryFace=1000;

    //build "virtual" blocks for the combination of the boundaries
    cOffset *OffsetTable;
    int OffsetTableLength;

    switch (xBoundaryFace+yBoundaryFace+zBoundaryFace) {
     case 10+0+0:
      OffsetTable=xOffset;
      OffsetTableLength=1;
      break;

     case 0+100+0:
      OffsetTable=yOffset;
      OffsetTableLength=1;
      break;

     case 0+0+1000:
      OffsetTable=zOffset;
      OffsetTableLength=1;
      break;

     case 10+100+0:
      OffsetTable=xyOffset;
      OffsetTableLength=3;
      break;

     case 10+0+1000:
      OffsetTable=xzOffset;
      OffsetTableLength=3;
      break;

     case 0+100+1000:
      OffsetTable=yzOffset;
      OffsetTableLength=3;
      break;

     case 10+100+1000:
      OffsetTable=xyzOffset;
      OffsetTableLength=7;
      break;
    }

    nVirtualBlocks+=OffsetTableLength;

    //loop through all 'virtual' blocks
    for (iblock=0;iblock<OffsetTableLength;iblock++) {
      //loop through all 'internal' corner nodes of the block
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        if ( (i+OffsetTable[iblock].di<=_BLOCK_CELLS_X_) && (j+OffsetTable[iblock].dj<=_BLOCK_CELLS_Y_) && (k+OffsetTable[iblock].dk<=_BLOCK_CELLS_Z_) ) {
          //the point is located at the boundary of the domain and the user-defined funtion is needed to be called to get the interpolation stencil
          int NonZeroElementsFound;
          double rhs;

          NewRow=MatrixRowStack.newElement();

          //call the user defined function to determine the non-zero elements of the matrix
          NewRow->RhsSupportLength=0;
          rhs=0.0;

          f(i+OffsetTable[iblock].di,j+OffsetTable[iblock].dj,k+OffsetTable[iblock].dk,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable,NewRow->RhsSupportLength,node);
          if (NonZeroElementsFound>MaxStencilLength) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

          //populate the new row
          NewRow->i=i+OffsetTable[iblock].di,NewRow->j=j+OffsetTable[iblock].dj,NewRow->k=k+OffsetTable[iblock].dk;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=rhs;
          NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(NewRow->i,NewRow->j,NewRow->k));
          NewRow->nNonZeroElements=NonZeroElementsFound;

          //add the new row to the matrix
          if (MatrixRowTable==NULL) {
            MatrixRowLast=NewRow;
            MatrixRowTable=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowLast->next=NewRow;
            MatrixRowLast=NewRow;
            NewRow->next=NULL;
          }

          //add to the row non-zero elements
          for (int ii=0;ii<NonZeroElementsFound;ii++) {
            MatrixRowNonZeroElementTable[ii].Node=node;

            if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0);
            }
            else if (MatrixRowNonZeroElementTable[ii].i<0) {
              MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0);
            }

            if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0);
            }
            else if (MatrixRowNonZeroElementTable[ii].j<0) {
              MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0);
            }

            if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0)!=NULL)) {
              MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0);
            }
            else if (MatrixRowNonZeroElementTable[ii].k<0) {
              MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
              MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0);
            }

            if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
              exit(__LINE__,__FILE__,"Error: the block is not found");
            }

            cCornerNode* CornerNode=MatrixRowNonZeroElementTable[ii].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k));

            if (CornerNode==NULL) exit(__LINE__,__FILE__,"Error: something is wrong");

            el=NewRow->Elements+ii;

            el->CornerNode=CornerNode;
            el->CornerNodeID=PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[ii].i,MatrixRowNonZeroElementTable[ii].j,MatrixRowNonZeroElementTable[ii].k);
            el->UnknownVectorIndex=-1;
            el->MatrixElementValue=MatrixRowNonZeroElementTable[ii].MatrixElementValue;
            el->node=MatrixRowNonZeroElementTable[ii].Node;
            el->iVar=iVar;

            //count the number of the element that are needed
            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[ii].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }
        else {
          //the point is outside of the domain
          NewRow=MatrixRowStack.newElement();

          NewRow->i=-1000000,NewRow->j=-1000000,NewRow->k=-1000000;
          NewRow->iVar=iVar;
          NewRow->node=node;
          NewRow->Rhs=1.0;
          NewRow->CornerNode=NULL;
          NewRow->nNonZeroElements=1;
          NewRow->RhsSupportLength=0;

          NewRow->Elements[0].CornerNode=NULL;
          NewRow->Elements[0].CornerNodeID=-1;
          NewRow->Elements[0].UnknownVectorIndex=-1; //Index;
          NewRow->Elements[0].MatrixElementValue=1.0;
          NewRow->Elements[0].node=node;
          NewRow->Elements[0].iVar=iVar;

          //add the new row to the matrix
          if (MatrixRowTable==NULL) {
            MatrixRowLast=NewRow;
            MatrixRowTable=NewRow;
            NewRow->next=NULL;
          }
          else {
            MatrixRowLast->next=NewRow;
            MatrixRowLast=NewRow;
            NewRow->next=NULL;
          }

//          if (iVar==0) Index++;
        }
      }
    }
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& Rhs,cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k,thread,nLocalNode;
  int iRow=0;
  cMatrixRow* Row;
  cStencilElement* el;
  cCornerNode* CornerNode;

  int Debug=0;

  nRealBlocks=0;

  //reset indexing of the nodes
  Reset();

  //allocate the counter of the data points to be recieve
  int RecvDataPointCounter[PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) RecvDataPointCounter[thread]=0;

  list<cLinearSystemCornerNodeDataRequestListElement> DataRequestList[PIC::nTotalThreads];

  //build the matrix
  for (nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    //in case of the periodic boundary condition it is only the points that are inside the "real" computational domain that are considered
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;

      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
        //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
        BoundaryBlock=true;
        break;
      }

      if (BoundaryBlock==true) continue;
    }

    nRealBlocks++;

    //the limits are correct: the point i==_BLOCK_CELLS_X_ belongs to the onother block
    int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        int NonZeroElementsFound;
        double rhs;

        //create the new Row entry
        cMatrixRow* NewRow=MatrixRowStack.newElement();

        //call the user defined function to determine the non-zero elements of the matrix
        NewRow->RhsSupportLength=0;
        rhs=0.0;

        f(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,NewRow->RhsSupportTable,NewRow->RhsSupportLength,node);
        if (NonZeroElementsFound>MaxStencilLength) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

        //scan through the found stencil and correct blocks and indexing is needed
        for (int ii=0;ii<NonZeroElementsFound;ii++) {
          MatrixRowNonZeroElementTable[ii].Node=node;

          if ((MatrixRowNonZeroElementTable[ii].i>=iMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].i<0) {
            MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0);
          }

          if ((MatrixRowNonZeroElementTable[ii].j>=jMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].j<0) {
            MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0);
          }

          if ((MatrixRowNonZeroElementTable[ii].k>=kMax) && (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0)!=NULL)) {
            MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].k<0) {
            MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0);
          }

          if (MatrixRowNonZeroElementTable[ii].Node==NULL) {
            exit(__LINE__,__FILE__,"Error: the block is not found");
          }

          //check whether the periodic boundary conditions are in use, and the new block is on the boundary of the domain
          if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
            bool BoundaryBlock=false;

            for (int iface=0;iface<6;iface++) if (MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(iface,0,0)==NULL) {
              //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
              BoundaryBlock=true;
              break;
            }

            if (BoundaryBlock==true) {
              //the block is at the domain boundary -> find the point located in the "real" part of the computational domain
              MatrixRowNonZeroElementTable[ii].Node=PIC::BC::ExternalBoundary::Periodic::findCorrespondingRealBlock(MatrixRowNonZeroElementTable[ii].Node);
            }
          }

        }

        //populate the new row
        NewRow->i=i,NewRow->j=j,NewRow->k=k;
        NewRow->iVar=iVar;
        NewRow->node=node;
        NewRow->Rhs=rhs;
        NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
        NewRow->nNonZeroElements=NonZeroElementsFound;


        if (MatrixRowTable==NULL) {
          MatrixRowLast=NewRow;
          MatrixRowTable=NewRow;
          NewRow->next=NULL;
        }
        else {
          MatrixRowLast->next=NewRow;
          MatrixRowLast=NewRow;
          NewRow->next=NULL;
        }

        //add to the row non-zero elements
        for (int iElement=0;iElement<NonZeroElementsFound;iElement++) {
          el=NewRow->Elements+iElement;
          CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));

          el->CornerNode=CornerNode;
          el->CornerNodeID=PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k);
          el->UnknownVectorIndex=-1;
          el->MatrixElementValue=MatrixRowNonZeroElementTable[iElement].MatrixElementValue;
          el->node=MatrixRowNonZeroElementTable[iElement].Node;
          el->iVar=iVar;

          //count the number of the element that are needed
          if (MatrixRowNonZeroElementTable[iElement].Node->Thread==PIC::ThisThread) {
            //the point is located at the current MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[PIC::ThisThread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
          else {
            //the data point is on another MPI process

            if (CornerNode->LinearSolverUnknownVectorIndex==-1) {
              RecvDataPointCounter[MatrixRowNonZeroElementTable[iElement].Node->Thread]++;
              CornerNode->LinearSolverUnknownVectorIndex=-2;
            }
          }
        }


      }
    }
  }


  //add 'virtual' blocks at the right boundary of the domain
  if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_OFF_) {
    ProcessRightDomainBoundary(RecvDataPointCounter,f);
  }


  //allocate the data buffers for the partial vectors and link the pointers points that are inside the subdomian
  //allocate the exchange buffer for the data the needs to be recieved from other MPI processess
  int cnt=0,iElement;

  int *DataExchangeTableCounter=new int [PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) DataExchangeTableCounter[thread]=0;

  for (iRow=0,Row=MatrixRowTable;Row!=NULL;iRow++,Row=Row->next) {

    for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
      el=Row->Elements+iElement;

      if (el->CornerNode!=NULL) {
        if (el->CornerNode->LinearSolverUnknownVectorIndex==-2) {
          //the node is still not linked to the 'partial data vector'
          el->CornerNode->LinearSolverUnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;

          if (el->node->Thread!=PIC::ThisThread) {
            //add the point to the data request list
            cLinearSystemCornerNodeDataRequestListElement DataRequestListElement;
            cAMRnodeID NodeID;

            PIC::Mesh::mesh.GetAMRnodeID(NodeID,el->node);
            DataRequestListElement.CornerNodeID=el->CornerNodeID;
            DataRequestListElement.NodeID=NodeID;

            DataRequestList[el->node->Thread].push_back(DataRequestListElement);
          }
        }

        el->UnknownVectorIndex=el->CornerNode->LinearSolverUnknownVectorIndex;
        el->Thread=el->node->Thread;
      }
      else {
        //the point is within the 'virtual' block
        el->UnknownVectorIndex=DataExchangeTableCounter[el->node->Thread]++;
        el->Thread=el->node->Thread;
      }
    }
  }

  //re-count the 'internal' corner nodes so they follow the i,j,k order as well as rows
  for (Row=MatrixRowTable;Row!=NULL;Row=Row->next) if (Row->CornerNode!=NULL) Row->CornerNode->LinearSolverUnknownVectorIndex=-1;

  for (iRow=0,Row=MatrixRowTable;Row!=NULL;Row=Row->next) {
    if (Row->CornerNode!=NULL) {
      if (Row->iVar==0) {
        if (Row->CornerNode->LinearSolverUnknownVectorIndex>=0) exit(__LINE__,__FILE__,"Error: a courner node is double counted");
        Row->CornerNode->LinearSolverUnknownVectorIndex=iRow;
      }
    }
    else Row->Elements[0].UnknownVectorIndex=iRow; //rows corresponding to the virtual blocks has only one non-zero element

    if (Row->iVar==NodeUnknownVariableVectorLength-1) iRow++; //all rows pointing to the same 'CornerNode' have the same 'iRow'
  }

  //verify that all all corner nodes have been counted
  for (Row=MatrixRowTable;Row!=NULL;iRow++,Row=Row->next) for (iElement=0;iElement<Row->nNonZeroElements;iElement++) {
    int n;

    if (Row->Elements[iElement].CornerNode!=NULL) {
      if ((n=Row->Elements[iElement].CornerNode->LinearSolverUnknownVectorIndex)<0) {
        exit(__LINE__,__FILE__,"Error: an uncounted corner node has been found");
      }

      Row->Elements[iElement].UnknownVectorIndex=n;
    }
  }


  //exchange the request lists
  int From,To;
  int *nGlobalDataPointTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
  cLinearSystemCornerNodeDataRequestListElement* ExchangeList;

  MPI_Allgather(DataExchangeTableCounter,PIC::nTotalThreads,MPI_INT,nGlobalDataPointTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

  SendExchangeBuffer=new double* [PIC::nTotalThreads];
  SendExchangeBufferLength=new int [PIC::nTotalThreads];
  SendExchangeBufferElementIndex=new int* [PIC::nTotalThreads];

  RecvExchangeBuffer=new double* [PIC::nTotalThreads];
  RecvExchangeBufferLength=new int [PIC::nTotalThreads];

  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    SendExchangeBuffer[thread]=NULL,SendExchangeBufferLength[thread]=0,SendExchangeBufferElementIndex[thread]=0;
    RecvExchangeBuffer[thread]=NULL,RecvExchangeBufferLength[thread]=0;
  }


  //create the lists
  for (To=0;To<PIC::nTotalThreads;To++) for (From=0;From<PIC::nTotalThreads;From++) if ( (To!=From) && (nGlobalDataPointTable[From+To*PIC::nTotalThreads]!=0) && ((To==PIC::ThisThread)||(From==PIC::ThisThread)) ) {
    ExchangeList=new cLinearSystemCornerNodeDataRequestListElement [nGlobalDataPointTable[From+To*PIC::nTotalThreads]];

    if (PIC::ThisThread==To) {
      list<cLinearSystemCornerNodeDataRequestListElement>::iterator ptr;

      for (i=0,ptr=DataRequestList[From].begin();ptr!=DataRequestList[From].end();ptr++,i++) ExchangeList[i]=*ptr;
      MPI_Send(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,From,0,MPI_GLOBAL_COMMUNICATOR);

      //create the recv list
      RecvExchangeBufferLength[From]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
      RecvExchangeBuffer[From]=new double[NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From]];
    }
    else {
      MPI_Status status;

      MPI_Recv(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,To,0,MPI_GLOBAL_COMMUNICATOR,&status);

      //unpack the SendDataList
      SendExchangeBufferLength[To]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
      SendExchangeBuffer[To]=new double[NodeUnknownVariableVectorLength*SendExchangeBufferLength[To]];
      SendExchangeBufferElementIndex[To]=new int[SendExchangeBufferLength[To]];

      for (int ii=0;ii<SendExchangeBufferLength[To];ii++) {
        node=PIC::Mesh::mesh.findAMRnodeWithID(ExchangeList[ii].NodeID);
        CornerNode=node->block->GetCornerNode(ExchangeList[ii].CornerNodeID);

        SendExchangeBufferElementIndex[To][ii]=CornerNode->LinearSolverUnknownVectorIndex;
      }
    }

    delete [] ExchangeList;
  }

  //deallocate 'nGlobalDataPointTable'
  delete [] nGlobalDataPointTable;
  delete [] DataExchangeTableCounter;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::ExchangeIntermediateUnknownsData(double *x) {
  int To,From;
  MPI_Request SendRequest[PIC::nTotalThreads],RecvRequest[PIC::nTotalThreads];
  MPI_Status SendStatus[PIC::nTotalThreads],RecvStatus[PIC::nTotalThreads];
  int RecvThreadCounter=0,SendThreadCounter=0;


  //initiate the non-blocked recieve
  for (From=0;From<PIC::nTotalThreads;From++) if (RecvExchangeBufferLength[From]!=0) {
    MPI_Irecv(RecvExchangeBuffer[From],NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From],MPI_DOUBLE,From,0,MPI_GLOBAL_COMMUNICATOR,RecvRequest+RecvThreadCounter);
    RecvThreadCounter++;
  }

  //prepare data to send and initiate the non-blocked send
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=PIC::ThisThread)&&(SendExchangeBufferLength[To]!=0)) {
    int i,offset=0;
    double *Buffer=SendExchangeBuffer[To];

    for (i=0;i<SendExchangeBufferLength[To];i++) {
      memcpy(Buffer+offset,x+NodeUnknownVariableVectorLength*SendExchangeBufferElementIndex[To][i],NodeUnknownVariableVectorLength*sizeof(double));
      offset+=NodeUnknownVariableVectorLength;
    }

    MPI_Isend(Buffer,NodeUnknownVariableVectorLength*SendExchangeBufferLength[To],MPI_DOUBLE,To,0,MPI_GLOBAL_COMMUNICATOR,SendRequest+SendThreadCounter);
    SendThreadCounter++;
  }

  //finalize send and recieve
  MPI_Waitall(RecvThreadCounter,RecvRequest,RecvStatus);
  MPI_Waitall(SendThreadCounter,SendRequest,SendStatus);
}


template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::Reset() {
  Reset(PIC::Mesh::mesh.rootTree);
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  if ((startNode==PIC::Mesh::mesh.rootTree)&&(RecvExchangeBuffer!=NULL)) {
    //clear the row stack
    MatrixRowStack.resetStack();

    //deallocate the allocated data buffers
    for (int thread=0;thread<PIC::nTotalThreads;thread++) {
      if (RecvExchangeBuffer!=NULL) if (RecvExchangeBuffer[thread]!=NULL) {
        delete [] RecvExchangeBuffer[thread];
      }

      if (SendExchangeBufferElementIndex!=NULL) if (SendExchangeBufferElementIndex[thread]!=NULL) {
        delete [] SendExchangeBufferElementIndex[thread];
        delete [] SendExchangeBuffer[thread];
      }
    }

    if (RecvExchangeBuffer!=NULL) {
      delete [] RecvExchangeBuffer;
      delete [] RecvExchangeBufferLength;
    }

    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL;

    if (SendExchangeBufferElementIndex!=NULL) {
      delete [] SendExchangeBuffer;
      delete [] SendExchangeBufferLength;
      delete [] SendExchangeBufferElementIndex;
    }

    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    if (SubdomainPartialRHS!=NULL) {
      delete [] SubdomainPartialRHS;
      delete [] SubdomainPartialUnknownsVector;
    }

    MatrixRowTable=NULL,MatrixRowLast=NULL;
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
  }


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=NULL;
    PIC::Mesh::cDataCornerNode *CornerNode=NULL;

    if ((block=startNode->block)!=NULL) {
      for (int i=0;i<_BLOCK_CELLS_X_+1;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
        if ((CornerNode=block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k)))!=NULL) {
          CornerNode->LinearSolverUnknownVectorIndex=-1;
        }
      }
    }

  }
  else for (int i=0;i<(1<<DIM);i++) if (startNode->downNode[i]!=NULL) Reset(startNode->downNode[i]);
}


template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::UpdateRhs(double (*fSetRhs)(int,cRhsSupportTable*,int)) {
  cMatrixRow* row;

  for (row=MatrixRowTable;row!=NULL;row=row->next) if (row->RhsSupportLength!=0) {
    row->Rhs=fSetRhs(row->iVar,row->RhsSupportTable,row->RhsSupportLength);
  }
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::MultiplyVector(double *p,double *x,int length) {
  cMatrixRow* row;
  cStencilElement StencilElement,*Elements;
  int cnt,iElement,iElementMax;
  double res;

  ExchangeIntermediateUnknownsData(x);

  RecvExchangeBuffer[PIC::ThisThread]=x;

  for (row=MatrixRowTable,cnt=0;row!=NULL;row=row->next,cnt++) {
    iElementMax=row->nNonZeroElements;
    Elements=row->Elements;

    for (res=0.0,iElement=0;iElement<iElementMax;iElement++) {
      memcpy(&StencilElement,Elements+iElement,sizeof(cStencilElement));

      res+=StencilElement.MatrixElementValue*RecvExchangeBuffer[StencilElement.Thread][StencilElement.iVar+NodeUnknownVariableVectorLength*StencilElement.UnknownVectorIndex];
    }

     p[cnt]=res;
  }

  RecvExchangeBuffer[PIC::ThisThread]=NULL;
}

template <class cCornerNode, int NodeUnknownVariableVectorLength,int MaxStencilLength, int MaxRhsSupportLength>
void cLinearSystemCornerNode<cCornerNode, NodeUnknownVariableVectorLength,MaxStencilLength, MaxRhsSupportLength>::Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode)) {
  cMatrixRow* row;
  int GlobalVariableIndex,cnt;

  //count the total number of the variables, and create the data buffers
  if (SubdomainPartialRHS==NULL) {
    for (SubdomainPartialUnknownsVectorLength=0,row=MatrixRowTable;row!=NULL;row=row->next) SubdomainPartialUnknownsVectorLength+=NodeUnknownVariableVectorLength;

    SubdomainPartialRHS=new double [SubdomainPartialUnknownsVectorLength];
    SubdomainPartialUnknownsVector=new double [SubdomainPartialUnknownsVectorLength];
  }

  //calculate the mean value of the Rhs that could be used later for the point located in the 'virtual' blocks
  int i;
  double MeanRhs[NodeUnknownVariableVectorLength];
  int MeanRhsCounter[NodeUnknownVariableVectorLength];

  for (i=0;i<NodeUnknownVariableVectorLength;i++) MeanRhs[i]=0.0,MeanRhsCounter[i]=0;

  for (row=MatrixRowTable;row!=NULL;cnt++,row=row->next) if (row->CornerNode!=NULL) {
    MeanRhs[row->iVar]+=row->Rhs;
    MeanRhsCounter[row->iVar]++;
  }

  for (i=0;i<NodeUnknownVariableVectorLength;i++) {
    if (MeanRhsCounter[i]==0) MeanRhs[i]=0.0;
    else MeanRhs[i]/=MeanRhsCounter[i];
  }


  //populate the buffers
  for (cnt=0,row=MatrixRowTable;row!=NULL;cnt++,row=row->next) {
    if (row->CornerNode!=NULL) {
      if (cnt%NodeUnknownVariableVectorLength==0) {
        fInitialUnknownValues(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
      }

      SubdomainPartialRHS[cnt]=row->Rhs;
    }
    else {
      //the point is located in the 'virtual' block
      SubdomainPartialUnknownsVector[row->iVar+NodeUnknownVariableVectorLength*row->Elements->UnknownVectorIndex]=MeanRhs[row->iVar];
      SubdomainPartialRHS[cnt]=MeanRhs[row->iVar];
    }


  }


  //call the iterative solver
  double Tol=1e-5;// the max iteration error allowed
  int nIter=100; //iter number
  int nVar=NodeUnknownVariableVectorLength; //variable number
  int nDim = 3; //dimension
  int nI=_BLOCK_CELLS_X_;
  int nJ=_BLOCK_CELLS_Y_;
  int nK=_BLOCK_CELLS_Z_;
  int nBlock = nRealBlocks+nVirtualBlocks;

  MPI_Fint iComm = MPI_Comm_c2f(MPI_GLOBAL_COMMUNICATOR);

  //double Rhs_I(nVar*nVar*nI*nJ*nK*nBlock)// RHS of the equation
  //double Sol_I(nVar*nVar*nI*nJ*nK*nBlock)// vector for solution
  double * Rhs_I=SubdomainPartialRHS;
  double * Sol_I=SubdomainPartialUnknownsVector;
  double PrecondParam=0; // not use preconditioner
  double ** precond_matrix_II;// pointer to precondition matrix; us  null if no preconditioner
  int lTest=1;//1: need test output; 0: no test statement

/*  printf("nBlock:%d\n",nBlock);
  for (int ii=0; ii<cnt;ii++){
    printf("No.%d, RHS_I:%f\n", ii, Rhs_I[ii]);
  }*/

  linear_solver_wrapper("GMRES", &Tol,&nIter, &nVar, &nDim,&nI, &nJ, &nK, &nBlock, &iComm, Rhs_I,Sol_I, &PrecondParam, NULL, &lTest);

/*  for (int ii=0; ii<cnt;ii++){
    printf("No.%d, Sol_I:%f\n", ii, Sol_I[ii]);
  }*/


  //unpack the solution
  for (row=MatrixRowTable;row!=NULL;row=row->next) if (row->CornerNode!=NULL)  {
    fUnpackSolution(SubdomainPartialUnknownsVector+NodeUnknownVariableVectorLength*row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
  }

  //execute the data exchange between 'ghost' blocks
  PIC::Mesh::mesh.ParallelBlockDataExchange();
}

#endif /* _LINEARSYSTEMCORNERNODE_H_ */
