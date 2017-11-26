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

class cLinearSystemCornerNodeDataRequestListElement {
public:
  int CornerNodeID;
  cAMRnodeID NodeID;
};

template <class cCornerNode>
class cLinearSystemCornerNode {
public:
  int NodeUnknownVariableVectorLength;
  int NodeUnknownVariableVectorOffset;

  double *SubdomainPartialRHS,*SubdomainPartialUnknownsVector;
  int SubdomainPartialUnknownsVectorLength;

  class cStencilElement {
  public:
    cStencilElement* next;
    cCornerNode* CornerNode;
    int CornerNodeID;
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    int UnknownVectorIndex,iVar,Thread;

    cStencilElement() {
      next=NULL,node=NULL,CornerNode=NULL;
      MatrixElementValue=0.0;
      UnknownVectorIndex=-1,iVar=-1,Thread=-1,CornerNodeID=-1;
    }
  };

  class cMatrixRow {
  public:
    cMatrixRow* next;
    cStencilElement* FirstElement;
    double rhs;

    int i,j,k;
    int iVar;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    cCornerNode* CornerNode;

    cMatrixRow() {
      next=NULL,FirstElement=NULL,rhs=0.0,node=NULL,CornerNode=NULL;
      i=0,j=0,k=0,iVar=0;
    }
  };

  cStack<cStencilElement> StencilElementStack;
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
  static const int nMaxMatrixNonzeroElementsDefault=1000;
  int nMaxMatrixNonzeroElements;

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

  cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable;

  //constructor
  void ExplicitConstructor(int nUnknownsPerNode,int UnknownsVectorOffset,int MaxNonZeroElementNumber) {
    nMaxMatrixNonzeroElements=MaxNonZeroElementNumber;
    NodeUnknownVariableVectorLength=nUnknownsPerNode;
    NodeUnknownVariableVectorOffset=UnknownsVectorOffset;

    MatrixRowTable=NULL,MatrixRowLast=NULL;

    //exchange buffers
    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL;
    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;

    MatrixRowNonZeroElementTable=new cMatrixRowNonZeroElementTable[nMaxMatrixNonzeroElements];

    //partial data of the linear system to be solved
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
    SubdomainPartialUnknownsVectorLength=0;
  }

  cLinearSystemCornerNode(int nUnknownsPerNode,int UnknownsVectorOffset, int MaxNonZeroElementNumber) {ExplicitConstructor(nUnknownsPerNode,UnknownsVectorOffset,MaxNonZeroElementNumber);}
  cLinearSystemCornerNode(int nUnknownsPerNode,int UnknownsVectorOffset) {ExplicitConstructor(nUnknownsPerNode,UnknownsVectorOffset,nMaxMatrixNonzeroElementsDefault);}

  //void reset the data of the obsect to the default state
  void Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
  void Reset();

  //reset indexing of the nodes
  void ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

  //build the matrix
  void BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& rhs,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node));

  //exchange the data
  void ExchageIntermediateUnknownsData(double *x);

  //matrix/vector multiplication
  void MultiplyVector(double *p,double *x,int length);

  //call the linear system solver, and unpack the solution afterward
  void Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode));

  //destructor
  ~cLinearSystemCornerNode() {
    Reset();

    StencilElementStack.~cStack();
    MatrixRowStack.~cStack();
  }
};

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::ResetUnknownVectorIndex(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
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

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& rhs,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k,thread;
  int iRow=0;
  cMatrixRow* Row;
  cStencilElement* el;
  cCornerNode* CornerNode;

  //reset indexing of the nodes
  Reset();

  //allocate the counter of the data points to be recieve
  int RecvDataPointCounter[PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) RecvDataPointCounter[thread]=0;

  list<cLinearSystemCornerNodeDataRequestListElement> DataRequestList[PIC::nTotalThreads];

  //build the matrix
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    //in case of the periodic boundary condition it is only the points that are inside the "real" computational domain that are considered

    //the limits are correct: the point i==_BLOCK_CELLS_X_ belongs to the onother block
    int kMax=_BLOCK_CELLS_Z_,jMax=_BLOCK_CELLS_Y_,iMax=_BLOCK_CELLS_X_;

    if (node->GetNeibFace(1,0,0)==NULL) iMax=_BLOCK_CELLS_X_+1;
    if (node->GetNeibFace(3,0,0)==NULL) jMax=_BLOCK_CELLS_Y_+1;
    if (node->GetNeibFace(5,0,0)==NULL) kMax=_BLOCK_CELLS_Z_+1;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        int NonZeroElementsFound;
        double rhs;

        //call the user defined function to determine the non-zero elements of the matrix
        f(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs,node);
        if (NonZeroElementsFound>=nMaxMatrixNonzeroElements) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

        //scan through the found stencil and correct blocks and indexing is needed
        for (int ii=0;ii<NonZeroElementsFound;ii++) {
          MatrixRowNonZeroElementTable[ii].Node=node;

          if (MatrixRowNonZeroElementTable[ii].i>=iMax) {
            MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].i<0) {
            MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0);
          }

          if (MatrixRowNonZeroElementTable[ii].j>=jMax) {
            MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].j<0) {
            MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0);
          }

          if (MatrixRowNonZeroElementTable[ii].k>=kMax) {
            MatrixRowNonZeroElementTable[ii].k-=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(5,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].k<0) {
            MatrixRowNonZeroElementTable[ii].k+=_BLOCK_CELLS_Z_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(4,0,0);
          }

          //check whether the periodic boundary conditions are in use, and the new block is on the boundary of the domain
        }

        //create the new Row entry
        cMatrixRow* NewRow=MatrixRowStack.newElement();

        NewRow->i=i,NewRow->j=j,NewRow->k=k;
        NewRow->iVar=iVar;
        NewRow->node=node;
        NewRow->rhs=rhs;
        NewRow->FirstElement=NULL;
        NewRow->CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));


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
          el=StencilElementStack.newElement();
          CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));

          el->CornerNode=CornerNode;
          el->CornerNodeID=PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k);
          el->UnknownVectorIndex=-1;
          el->MatrixElementValue=MatrixRowNonZeroElementTable[iElement].MatrixElementValue;
          el->node=MatrixRowNonZeroElementTable[iElement].Node;
          el->iVar=iVar;

          el->next=NewRow->FirstElement;
          NewRow->FirstElement=el;

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

  //allocate the data buffers for the partial vectors and link the pointers points that are inside the subdomian
  //allocate the exchange buffer for the data the needs to be recieved from other MPI processess
  int cnt=0;

  int DataExchangeTableCounter[PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) DataExchangeTableCounter[thread]=0;

  for (iRow=0,Row=MatrixRowTable;Row!=NULL;iRow++,Row=Row->next) {
    for (el=Row->FirstElement;el!=NULL;el=el->next) {
      if (el->CornerNode->LinearSolverUnknownVectorIndex==-2) {
        //the node is still not linked to the 'partial data vector'
        el->CornerNode->LinearSolverUnknownVectorIndex=DataExchangeTableCounter[PIC::ThisThread]++;

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
  }

  //exchange the request lists
  int From,To;
  int *nGlobalDataPointTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
  cLinearSystemCornerNodeDataRequestListElement* ExchangeList;

  MPI_Allgather(DataExchangeTableCounter,PIC::nTotalThreads,MPI_INT,nGlobalDataPointTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

  SendExchangeBuffer=new double* [PIC::nTotalThreads];
  SendExchangeBufferLength=new int [PIC::nTotalThreads];

  RecvExchangeBuffer=new double* [PIC::nTotalThreads];
  RecvExchangeBufferLength=new int [PIC::nTotalThreads];

  for (thread=0;thread<PIC::nTotalThreads;thread++) {
    SendExchangeBuffer[thread]=NULL,SendExchangeBufferLength[thread]=0;
    RecvExchangeBuffer[thread]=NULL,RecvExchangeBufferLength[thread]=0;
  }


  //create the lists
  for (To=0;To<PIC::nTotalThreads;To++) for (From=0;From<PIC::nTotalThreads;From++) if ( (To!=From) && (nGlobalDataPointTable[From+To*PIC::nTotalThreads]!=0) && ((To==PIC::ThisThread)||(From==PIC::ThisThread)) ) {
    ExchangeList=new cLinearSystemCornerNodeDataRequestListElement [nGlobalDataPointTable[From+To*PIC::nTotalThreads]];

    if (PIC::ThisThread==To) {
      list<cLinearSystemCornerNodeDataRequestListElement>::iterator ptr;

      for (i=0,ptr=DataRequestList[From].begin();ptr!=DataRequestList[From].end();ptr++,i++) ExchangeList[i]=*ptr;
      MPI_Send(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cLinearSystemCornerNodeDataRequestListElement),MPI_CHAR,0,From,MPI_GLOBAL_COMMUNICATOR);

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
}

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::ExchageIntermediateUnknownsData(double *x) {
  int To,From;
  MPI_Request SendRequest[PIC::nTotalThreads],RecvRequest[PIC::nTotalThreads];
  MPI_Status SendStatus[PIC::nTotalThreads],RecvStatus[PIC::nTotalThreads];
  int RecvThreadCounter=0,SendThreadCounter=0;


  //initiate the non-blocked recieve
  for (From=0;From<PIC::nTotalThreads;From++) if (RecvExchangeBufferLength[From]!=0) {
    MPI_Irecv(RecvExchangeBuffer[From],NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From],MPI_DOUBLE,From,0,MPI_GLOBAL_COMMUNICATOR,RecvRequest+RecvThreadCounter);
    RecvThreadCounter++;
  }

  //prepare data to send and initial send {
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=PIC::ThisThread)&&(SendExchangeBufferLength[To]!=0)) {
    int i,j,cnt=0;

    for (i=0;i<SendExchangeBufferLength[To];i++) for (j=0;j<NodeUnknownVariableVectorLength;j++) {
      SendExchangeBuffer[To][cnt++]=x[j+NodeUnknownVariableVectorLength*SendExchangeBufferElementIndex[To][i]];
    }

    MPI_Isend(SendExchangeBuffer[To],NodeUnknownVariableVectorLength*SendExchangeBufferLength[To],MPI_DOUBLE,To,0,MPI_GLOBAL_COMMUNICATOR,SendRequest+SendThreadCounter);
    SendThreadCounter++;
  }

  //finalize send and recieve
  MPI_Waitall(RecvThreadCounter,RecvRequest,RecvStatus);
  MPI_Waitall(SendThreadCounter,SendRequest,SendStatus);

}


template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::Reset() {
  Reset(PIC::Mesh::mesh.rootTree);
}

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  if ((startNode==PIC::Mesh::mesh.rootTree)&&(RecvExchangeBuffer!=NULL)) {
    //clear the stacks
    StencilElementStack.resetStack();
    MatrixRowStack.resetStack();

    //deallocate the allocated data buffers
    for (int thread=0;thread<PIC::nTotalThreads;thread++) {
      if (RecvExchangeBuffer[thread]!=NULL) delete [] RecvExchangeBuffer[thread];

      if (SendExchangeBufferElementIndex[thread]!=NULL) {
        delete [] SendExchangeBufferElementIndex;
        delete [] SendExchangeBuffer;
      }
    }

    delete [] RecvExchangeBuffer;
    delete [] RecvExchangeBufferLength;

    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL;

    delete [] SendExchangeBuffer;
    delete [] SendExchangeBufferLength;
    delete [] SendExchangeBufferElementIndex;

    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL,SendExchangeBufferElementIndex=NULL;


    delete [] SubdomainPartialRHS;
    delete [] SubdomainPartialUnknownsVector;

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

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::MultiplyVector(double *p,double *x,int length) {
  cMatrixRow* row;
  cStencilElement *el,StencilElement;
  int cnt;
  double res;

  ExchageIntermediateUnknownsData(x);

  for (row=MatrixRowTable,cnt=0;row!=NULL;row=row->next,cnt++) {
    for (res=0.0,el=row->FirstElement;el!=NULL;el++) {
      memcpy(&StencilElement,el,sizeof(cStencilElement));

      res+=StencilElement.MatrixElementValue*((StencilElement.Thread==PIC::ThisThread) ? x[StencilElement.iVar+NodeUnknownVariableVectorLength*StencilElement.UnknownVectorIndex] : RecvExchangeBuffer[StencilElement.Thread][StencilElement.iVar+NodeUnknownVariableVectorLength*StencilElement.UnknownVectorIndex]);
    }

     p[cnt]=res;
  }
}

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::Solve(void (*fInitialUnknownValues)(double* x,cCornerNode* CornerNode),void (*fUnpackSolution)(double* x,cCornerNode* CornerNode)) {
  cMatrixRow* row;
  int GlobalVariableIndex;

  //count the total number of the variables, and create the data buffers
  if (SubdomainPartialRHS==NULL) {
    for (row=MatrixRowTable;row!=NULL;row=row->next) SubdomainPartialUnknownsVectorLength+=NodeUnknownVariableVectorLength;

    SubdomainPartialRHS=new double [SubdomainPartialUnknownsVectorLength];
    SubdomainPartialUnknownsVector=new double [SubdomainPartialUnknownsVectorLength];
  }

  //populate the buffers
  for (row=MatrixRowTable;row!=NULL;row=row->next) fInitialUnknownValues(SubdomainPartialUnknownsVector+row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);
/*
  //call the iterative solver
  double Tol=1e-5;// the max iteration error allowed
  int nIter=200; //iter number
  int nVar=1; //variable number
  int nDim =1; //dimension
  int nI, nJ,nK,nBlock;

  MPI_Fint iComm = MPI_Comm_c2f(MPI_COMM_WORLD);

  double Rhs_I(nVar*nVar*nI*nJ*nK*nBlock)// RHS of the equation
  double Sol_I(nVar*nVar*nI*nJ*nK*nBlock)// vector for solution
  double PrecondParam=0; // not use preconditioner
  double ** precond_matrix_II;// pointer to precondition matrix; us  null if no preconditioner
  int lTest=1;//1: need test output; 0: no test statement

  linear_solver_wrapper("GMRES", &Tol,&nIter, &nVar, &nDim,
                         &nI, &nJ, &nK, &nBlock, &iComm, Rhs_I,
  Sol_I, &PrecondParam, precond_matrix_II[0],
                         &lTest);*/


  //unpack the solution
  for (row=MatrixRowTable;row!=NULL;row=row->next) fUnpackSolution(SubdomainPartialUnknownsVector+row->CornerNode->LinearSolverUnknownVectorIndex,row->CornerNode);

  //execute the data exchange between 'ghost' blocks
  PIC::Mesh::mesh.ParallelBlockDataExchange();
}

#endif /* _LINEARSYSTEMCORNERNODE_H_ */
