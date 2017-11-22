
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

template <class cCornerNode>
class cLinearSystemCornerNode {
public:
  int NodeUnknownVariableVectorLength;
  int NodeUnknownVariableVectorOffset;

  double *SubdomainPartialRHS,*SubdomainPartialUnknownsVector;

  class cStencilElement {
  public:
    cStencilElement* next;
    cCornerNode* CornerNode;
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    int UnknownVectorIndex,Thread;

    cStencilElement() {
      next=NULL,node=NULL,CornerNode=NULL;
      MatrixElementValue=0.0;
      UnknownVectorIndex=-1,Thread=-1;
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

    cMatrixRow() {
      next=NULL,FirstElement=NULL,rhs=0.0,node=NULL;
      i=0,j=0,k=0,iVar=0;
    }
  };

  cStack<cStencilElement> StencilElementStack;
  cStack<cMatrixRow> MatrixRowStack;

  //Table of the rows local to the current MPI process
  cMatrixRow* MatrixRowTable;

  //exchange buffers
  int** RecvExchangeBuffer;
  int* RecvExchangeBufferLength;

  int** SendExchangeBuffer;
  int* SendExchangeBufferLength;

  int *RecvDataPointCounter;

  //add new row to the matrix
  //for that a user function is called that returns stenal (elements of the matrix), corresponding rhs, and the number of the elements in the amtrix' line
  static const int nMaxMatrixNonzeroElementsDefault=1000;
  int nMaxMatrixNonzeroElements;

  class cMatrixRowNonZeroElementTable {
  public:
    int i,j,k; //coordintes of the corner block used in the equation
    double MatrixElementValue;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node;

    cMatrixRowNonZeroElementTable() {
      i=0,j=0,k=0,MatrixElementValue=0.0,Node=NULL;
    }
  };

  cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable;

  //constructor
  void ExplicitConstructor(int nUnknownsPerNode,int UnknownsVectorOffset,int MaxNonZeroElementNumber) {
    nMaxMatrixNonzeroElements=MaxNonZeroElementNumber;
    NodeUnknownVariableVectorLength=nUnknownsPerNode;
    NodeUnknownVariableVectorOffset=UnknownsVectorOffset;

    MatrixRowTable=NULL;

    //exchange buffers
    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL,RecvDataPointCounter=NULL;
    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL;

    MatrixRowNonZeroElementTable=new cMatrixRowNonZeroElementTable[nMaxMatrixNonzeroElements];

    //partial data of the linear system to be solved
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;
  }

  cLinearSystemCornerNode(int nUnknownsPerNode,int UnknownsVectorOffset, int MaxNonZeroElementNumber) {ExplicitConstructor(nUnknownsPerNode,UnknownsVectorOffset,MaxNonZeroElementNumber);}
  cLinearSystemCornerNode(int nUnknownsPerNode,int UnknownsVectorOffset) {ExplicitConstructor(nUnknownsPerNode,UnknownsVectorOffset,nMaxMatrixNonzeroElementsDefault);}

  //void reset the data of the obsect to the default state
  void Reset(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
  void Reset();

  //build the matrix
  struct cDataRequestListElement {
    int i,j,k;
    cAMRnodeID NodeID;
  };

  void BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* Set,int& NonZeroElementsFound,double& rhs));

  //exchange the data
  void ExchageIntermediateUnknownsData();

  //matrix/vector multiplication
  void MultiplyVector(double *p,double *x,int length);

  //destructor
  ~cLinearSystemCornerNode() {
    Reset();

    StencilElementStack.~cStack();
    MatrixRowStack.~cStack();
  }
};

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::BuildMatrix(void(*f)(int i,int j,int k,int iVar,cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs)) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k,thread;
  int iRow=0;
  cMatrixRow* Row;
  cStencilElement* el;
  cCornerNode* CornerNode;

  //allocate the counter of the data points to be recieve
  RecvDataPointCounter=new int [PIC::nTotalThreads];
  for (thread=0;thread<PIC::nTotalThreads;thread++) RecvDataPointCounter[thread]=0;

  list<cDataRequestListElement> DataRequestList[PIC::ThisThread];

  //build the matrix
  for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

    //in case of the periodic boundary condition it is only the points that are inside the "real" computational domain that are considered

    //the limits are correct: the point i==_BLOCK_CELLS_X_ belongs to the onother block
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      for (int iVar=0;iVar<NodeUnknownVariableVectorLength;iVar++) {
        int NonZeroElementsFound;
        double rhs;

        //call the user defined function to determine the non-zero elements of the matrix
        f(i,j,k,iVar,MatrixRowNonZeroElementTable,NonZeroElementsFound,rhs);
        if (NonZeroElementsFound>=nMaxMatrixNonzeroElements) exit(__LINE__,__FILE__,"Error: NonZeroElementsFound>=nMaxMatrixNonzeroElement; Need to increase the value of nMaxMatrixNonzeroElement");

        //scan through the found stencil and correct blocks and indexing is needed
        for (int ii=0;ii<NonZeroElementsFound;ii++) {
          MatrixRowNonZeroElementTable[ii].Node=node;

          if (MatrixRowNonZeroElementTable[ii].i>=_BLOCK_CELLS_X_) {
            MatrixRowNonZeroElementTable[ii].i-=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(1,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].i<0) {
            MatrixRowNonZeroElementTable[ii].i+=_BLOCK_CELLS_X_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(0,0,0);
          }

          if (MatrixRowNonZeroElementTable[ii].j>=_BLOCK_CELLS_Y_) {
            MatrixRowNonZeroElementTable[ii].j-=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(3,0,0);
          }
          else if (MatrixRowNonZeroElementTable[ii].j<0) {
            MatrixRowNonZeroElementTable[ii].j+=_BLOCK_CELLS_Y_;
            MatrixRowNonZeroElementTable[ii].Node=MatrixRowNonZeroElementTable[ii].Node->GetNeibFace(2,0,0);
          }

          if (MatrixRowNonZeroElementTable[ii].k>=_BLOCK_CELLS_Z_) {
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

        NewRow->next=MatrixRowTable;
        MatrixRowTable=NewRow;

        //add to the row non-zero elements
        for (int iElement=0;iElement<NonZeroElementsFound;iElement++) {
          el=StencilElementStack.newElement();
          CornerNode=MatrixRowNonZeroElementTable[iElement].Node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(MatrixRowNonZeroElementTable[iElement].i,MatrixRowNonZeroElementTable[iElement].j,MatrixRowNonZeroElementTable[iElement].k));

          el->CornerNode=CornerNode;
          el->UnknownVectorIndex=-1;
          el->MatrixElementValue=MatrixRowNonZeroElementTable[iElement].MatrixElementValue;
          el->node=MatrixRowNonZeroElementTable[iElement].Node;

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

              //add into the data requst list
              cDataRequestListElement DataRequestListElement;
              cAMRnodeID NodeID;

              PIC::Mesh::mesh.GetAMRnodeID(NodeID,MatrixRowNonZeroElementTable[iElement].Node);

              DataRequestListElement.i=MatrixRowNonZeroElementTable[iElement].i;
              DataRequestListElement.j=MatrixRowNonZeroElementTable[iElement].j;
              DataRequestListElement.k=MatrixRowNonZeroElementTable[iElement].k;
              DataRequestListElement.NodeID=NodeID;

              DataRequestList[MatrixRowNonZeroElementTable[iElement].Node->Thread].push_front(DataRequestListElement);
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

  for (iRow=0,Row=MatrixRowTable;Row!=NULL;iRow,Row=Row->next) {
    for (el=Row->FirstElement;el!=NULL;el=el->next) {
      if (el->CornerNode->LinearSolverUnknownVectorIndex==-2) {
        //the node is still not linked to the 'partial data vector'
        el->CornerNode->LinearSolverUnknownVectorIndex=DataExchangeTableCounter[PIC::ThisThread]++;
      }

      el->UnknownVectorIndex=el->CornerNode->LinearSolverUnknownVectorIndex;
      el->Thread=el->node->Thread;
    }
  }

  //exchange the request lists
  int From,To;
  int *nGlobalDataPointTable=new int [PIC::nTotalThreads*PIC::nTotalThreads];
  cDataRequestListElement* ExchangeList;

  MPI_Allgather(DataExchangeTableCounter,PIC::nTotalThreads,MPI_INT,nGlobalDataPointTable,PIC::nTotalThreads,MPI_INT,MPI_GLOBAL_COMMUNICATOR);

  SendExchangeBuffer=new int* [PIC::nTotalThreads];
  SendExchangeBufferLength=new int [PIC::nTotalThreads];

  for (thread=0;thread<PIC::nTotalThreads;thread++) SendExchangeBuffer[thread]=NULL,SendExchangeBufferLength[thread]=0;


  //create the lists
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=From)&&(nGlobalDataPointTable[From+To*PIC::nTotalThreads]!=0)) {
    ExchangeList=new cDataRequestListElement [nGlobalDataPointTable[From+To*PIC::nTotalThreads]];

    if (PIC::ThisThread==To) {
//      list<cDataRequestListElement>::iterator ptr;

//      for (i=0,ptr=DataRequestList[From].begin();ptr!=DataRequestList[From].end();ptr++,i++) ExchangeList[i]=*ptr;
//      MPI_Send(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cDataRequestListElement),MPI_CHAR,0,From,MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Status status;

      MPI_Recv(ExchangeList,nGlobalDataPointTable[From+To*PIC::nTotalThreads]*sizeof(cDataRequestListElement),MPI_CHAR,To,0,MPI_GLOBAL_COMMUNICATOR,&status);

      //unpack the SendDataList
      SendExchangeBufferLength[To]=nGlobalDataPointTable[From+To*PIC::nTotalThreads];
      SendExchangeBuffer[To]=new int[SendExchangeBufferLength[To]];


      for (int ii=0;ii<SendExchangeBufferLength[To];ii++) {
        node=PIC::Mesh::mesh.findAMRnodeWithID(ExchangeList[ii].NodeID);
        CornerNode=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(ExchangeList[ii].i,ExchangeList[ii].j,ExchangeList[ii].k));

        SendExchangeBuffer[To][ii]=CornerNode->LinearSolverUnknownVectorIndex;
      }

      delete [] ExchangeList;
    }
  }

  //deallocate 'nGlobalDataPointTable'
  delete [] nGlobalDataPointTable;
}

template <class cCornerNode>
void cLinearSystemCornerNode<cCornerNode>::ExchageIntermediateUnknownsData() {
  int To,From;
  MPI_Request SendRequest[PIC::nTotalThreads],RecvRequest[PIC::nTotalThreads];
  MPI_Status SendStatus[PIC::nTotalThreads],RecvStatus[PIC::nTotalThreads];
  int RecvThreadCounter=0,SendThreadCounter=0;


  //initiate the non-blocked recieve
  for (From=0;From<PIC::nTotalThreads;From++) if (RecvExchangeBufferLength[From]!=0) {
    MPI_Irecv(RecvExchangeBuffer[From],NodeUnknownVariableVectorLength*RecvExchangeBufferLength[From],MPI_DOUBLE,From,0,MPI_GLOBAL_COMMUNICATOR,RecvRequest+nRecvThreads);
    RecvThreadCounter++;
  }

  //prepare data to send and initial send {
  for (To=0;To<PIC::nTotalThreads;To++) if ((To!=PIC::ThisThread)&&(SendExchangeBufferLength[To]!=0)) {
    int i,j,cnt=0;
    double *p;

    for (i=0;i<SendExchangeBufferLength[To];i++) for (p=SendExchangeDataPointerList[To],j=0;j<NodeUnknownVariableVectorLength;j++) {
      SendExchangeBuffer[To][cnt]=p[j];
    }

    MPI_Isend(SendExchangeBufferLength[To],NodeUnknownVariableVectorLength*SendExchangeBufferLength[To],MPI_DOUBLE,To,0,MPI_GLOBAL_COMMUNICATOR,SendRequest+To);
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

  if ((startNode==PIC::Mesh::mesh.rootTree)&&(RecvDataPointCounter!=NULL)) {
    //clear the stacks
    StencilElementStack.resetStack();
    MatrixRowStack.resetStack();

    //deallocate the allocated data buffers
    for (int thread=0;thread<PIC::nTotalThreads;thread++) if (RecvExchangeBuffer[thread]!=NULL) delete [] RecvExchangeBuffer[thread];

    delete [] RecvExchangeBuffer;
    delete [] RecvExchangeBufferLength;

    RecvExchangeBuffer=NULL,RecvExchangeBufferLength=NULL;

    for (int thread=0;thread<PIC::nTotalThreads;thread++) if (SendExchangeBuffer[thread]!=NULL) {
      delete [] SendExchangeBuffer[thread];
    }

    delete [] SendExchangeBuffer;
    delete [] SendExchangeBufferLength;

    delete [] SubdomainPartialRHS;
    delete [] SubdomainPartialUnknownsVector;

    SendExchangeBuffer=NULL,SendExchangeBufferLength=NULL;
    MatrixRowTable=NULL;
    SubdomainPartialRHS=NULL,SubdomainPartialUnknownsVector=NULL;

    delete [] RecvDataPointCounter;
    RecvDataPointCounter=NULL;
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

  ExchageIntermediateUnknownsData();

  for (row=MatrixRowTable,cnt=0,row!=NULL;row=row->next,cnt++) {
    for (res=0.0,el=row->FirstElement;el!=NULL;el++) {
      memcpy(&StencilElement,el,sizeof(cStencilElement));

      res+=StencilElement.MatrixElementValue*((StencilElement.Thread==PIC::ThisThread) ? x[StencilElement.UnknownVectorIndex] : RecvExchangeBuffer[StencilElement.Thread][StencilElement.UnknownVectorIndex];
    }

     p[cnt]=res;
  }
}



#endif /* _LINEARSYSTEMCORNERNODE_H_ */
