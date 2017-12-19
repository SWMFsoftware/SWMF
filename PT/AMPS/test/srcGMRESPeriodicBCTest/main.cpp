/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and vtenishe
 */

//$Id$


#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "meshAMRgeneric.h"

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"

int nVars=1; //number of variables in center associated data
double Background[3]={100.0,-20.0,10.0};


//#define _UNIFORM_MESH_ 1
//#define _NONUNIFORM_MESH_ 2

#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif


double xmin[3]={-1.0,-3.0,-4.0};
double xmax[3]={3.0,2.0,3.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;


//the solver used for solving the implicit wave transport equation
namespace TransportEquation {
  cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1> Solver;

  double c[3]={-0.4,0.0,0.0};

  void matvec(double* VecIn, double * VecOut, int n){
    Solver.MultiplyVector(VecOut,VecIn,n);
  }

  //create the stencil for the linear equaton solver test
  void GetTestStencilTimeDependent(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int& RhsSupportLength_CornerNodes,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int& RhsSupportLength_CenterNodes, 
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {


    MatrixRowNonZeroElementTable[0].i=i-1;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=-0.001;
    MatrixRowNonZeroElementTable[0].iVar=0;

    MatrixRowNonZeroElementTable[1].i=i+1;
    MatrixRowNonZeroElementTable[1].j=j;
    MatrixRowNonZeroElementTable[1].k=k;
    MatrixRowNonZeroElementTable[1].MatrixElementValue=0.001;
    MatrixRowNonZeroElementTable[1].iVar=0;


    MatrixRowNonZeroElementTable[2].i=i;
    MatrixRowNonZeroElementTable[2].j=j;
    MatrixRowNonZeroElementTable[2].k=k;
    MatrixRowNonZeroElementTable[2].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[2].iVar=0;


    RhsSupportTable_CornerNodes[0].Coefficient=1.0;
    RhsSupportTable_CornerNodes[0].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();
    RhsSupportLength_CornerNodes=1;

    RhsSupportLength_CenterNodes=0;
    NonZeroElementsFound=3;

   if (RhsSupportLength_CornerNodes!=0) for (i=0,rhs=0.0;i<RhsSupportLength_CornerNodes;i++) rhs+=*((double*)(RhsSupportTable_CornerNodes[i].AssociatedDataPointer+CurrentCornerNodeOffset))*RhsSupportTable_CornerNodes[i].Coefficient;

return ;




////////////////////////////////////

    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=0;

    RhsSupportTable_CornerNodes[0].Coefficient=1.0;
    RhsSupportTable_CornerNodes[0].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength_CornerNodes=1;
    RhsSupportLength_CenterNodes=0; 
    NonZeroElementsFound=1;

   if (RhsSupportLength_CornerNodes!=0) for (i=0,rhs=0.0;i<RhsSupportLength_CornerNodes;i++) rhs+=*((double*)(RhsSupportTable_CornerNodes[i].AssociatedDataPointer+CurrentCornerNodeOffset))*RhsSupportTable_CornerNodes[i].Coefficient;

return ;



///////////////////

     //'x'-direction
    MatrixRowNonZeroElementTable[0].i=i+1;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=c[0];
    MatrixRowNonZeroElementTable[0].iVar=0;

    MatrixRowNonZeroElementTable[1].i=i-1;
    MatrixRowNonZeroElementTable[1].j=j;
    MatrixRowNonZeroElementTable[1].k=k;
    MatrixRowNonZeroElementTable[1].MatrixElementValue=-c[0];
    MatrixRowNonZeroElementTable[1].iVar=0;

    //'y'-direction
    MatrixRowNonZeroElementTable[2].i=i;
    MatrixRowNonZeroElementTable[2].j=j+1;
    MatrixRowNonZeroElementTable[2].k=k;
    MatrixRowNonZeroElementTable[2].MatrixElementValue=c[1];
    MatrixRowNonZeroElementTable[2].iVar=0;

    MatrixRowNonZeroElementTable[3].i=i;
    MatrixRowNonZeroElementTable[3].j=j-1;
    MatrixRowNonZeroElementTable[3].k=k;
    MatrixRowNonZeroElementTable[3].MatrixElementValue=-c[1];
    MatrixRowNonZeroElementTable[3].iVar=0;

    //'z'-direction
    MatrixRowNonZeroElementTable[4].i=i;
    MatrixRowNonZeroElementTable[4].j=j;
    MatrixRowNonZeroElementTable[4].k=k+1;
    MatrixRowNonZeroElementTable[4].MatrixElementValue=c[2];
    MatrixRowNonZeroElementTable[4].iVar=0;

    MatrixRowNonZeroElementTable[5].i=i;
    MatrixRowNonZeroElementTable[5].j=j;
    MatrixRowNonZeroElementTable[5].k=k-1;
    MatrixRowNonZeroElementTable[5].MatrixElementValue=-c[2];
    MatrixRowNonZeroElementTable[5].iVar=0;

    //self
    MatrixRowNonZeroElementTable[6].i=i;
    MatrixRowNonZeroElementTable[6].j=j;
    MatrixRowNonZeroElementTable[6].k=k;
    MatrixRowNonZeroElementTable[6].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[6].iVar=0;

    //right-hand side
    RhsSupportTable_CornerNodes[0].Coefficient=c[0];
    RhsSupportTable_CornerNodes[0].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportTable_CornerNodes[1].Coefficient=-c[0];
    RhsSupportTable_CornerNodes[1].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();


    RhsSupportTable_CornerNodes[2].Coefficient=c[1];
    RhsSupportTable_CornerNodes[2].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();

    RhsSupportTable_CornerNodes[3].Coefficient=-c[1];
    RhsSupportTable_CornerNodes[3].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();


    RhsSupportTable_CornerNodes[4].Coefficient=c[2];
    RhsSupportTable_CornerNodes[4].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();

    RhsSupportTable_CornerNodes[5].Coefficient=-c[2];
    RhsSupportTable_CornerNodes[5].AssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();


    RhsSupportLength_CornerNodes=6;
    NonZeroElementsFound=7;

  if (RhsSupportLength_CornerNodes!=0) for (i=0,rhs=0.0;i<RhsSupportLength_CornerNodes;i++) rhs+=*((double*)(RhsSupportTable_CornerNodes[i].AssociatedDataPointer+CurrentCenterNodeOffset))*RhsSupportTable_CornerNodes[i].Coefficient;

  }




  //update the RHS vector
  double UpdateRhs(int iVar,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1>::cRhsSupportTable* RhsSupportTable_CornerNodes,int RhsSupportLength_CornerNodes,
    cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,6,1>::cRhsSupportTable* RhsSupportTable_CenterNodes,int RhsSupportLength_CenterNodes) {
    int i;
    double res=0.0;

    for (i=0;i<RhsSupportLength_CornerNodes;i++) res+=((double*)(RhsSupportTable_CornerNodes[i].AssociatedDataPointer+CurrentCornerNodeOffset))[iVar]*RhsSupportTable_CornerNodes[i].Coefficient;

    return res;
  }

  //set the initial guess
  void SetInitialGuess(double* x,PIC::Mesh::cDataCornerNode* CornerNode) {
   x[0]=*((double*)(CornerNode->GetAssociatedDataBufferPointer()+CurrentCornerNodeOffset));
 //   x[0]=0.0;
  }

  //process the solution vector
  void ProcessFinalSolution(double* x,PIC::Mesh::cDataCornerNode* CornerNode) {
    char *offset=CornerNode->GetAssociatedDataBufferPointer();

    ((double*)(offset+NextCornerNodeOffset))[0]=x[0]; ///+((double*)(offset+CurrentCornerNodeOffset))[0];
  }


  void BuildMatrix() {
    Solver.Reset();
    Solver.BuildMatrix(GetTestStencilTimeDependent);
  }

  void TimeStep() {

    Solver.UpdateRhs(UpdateRhs);

    linear_solver_matvec_c = matvec;
    Solver.Solve(SetInitialGuess,ProcessFinalSolution);
  }
}




void PrintCenterNodeVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"test wave (center node)\"");
}

void InterpolateCenterNode(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {

  double Wave=0.0;
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+CurrentCenterNodeOffset;
    Wave += (*((double*)SamplingBuffer))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+CurrentCenterNodeOffset,&Wave,sizeof(double));

}

void PrintCenterNodeData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+CurrentCenterNodeOffset));
  }
  
  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);
    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

}

void PrintCornerNodeVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"test wave (corner node)\"");
}

void PrintCornerNodeData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CornerNodeThread,PIC::Mesh::cDataCornerNode *CornerNode) {
  int idim;
  double t;

  if (pipe->ThisThread==CornerNodeThread) {
    t= *((double*)(CornerNode->GetAssociatedDataBufferPointer()+CurrentCornerNodeOffset));
  }

  if (pipe->ThisThread==0) {
    if (CornerNodeThread!=0) pipe->recv(t,CornerNodeThread);
    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
}


double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      static const double CharacteristicSpeed=1.0E-2;

      ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode); //CharacteristicSpeed=3.0;
      //return 0.3*startNode->GetCharacteristicCellSize()/CharacteristicSpeed;
      return dt;
    } else {
      CharacteristicSpeed=5.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
    }
#else
    //CharacteristicSpeed=5.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
    CharacteristicSpeed=100.0;
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    if (spec==_H_SPEC_) CharacteristicSpeed*=30.0;
    if (spec==_O_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_H2_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_OH_SPEC_) CharacteristicSpeed*=5.0;
#endif
#endif

    CellSize=startNode->GetCharacteristicCellSize();
    //return 0.3*CellSize/CharacteristicSpeed;

    return 1.0E-2;
}


double BulletLocalResolution(double *x) {                                                                                           
  double dist = xmax[0]-xmin[0];

#if _TEST_MESH_MODE_==_UNIFORM_MESH_  
  double res =4;
#endif

#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  double highRes = dist/32.0, lowRes= dist/2.0;     
  double res =(5-1)/dist*(x[0]-xmin[0])+1;  
#endif

  res=dist/pow(2,res);
  
  return res;
}
                       


void SetIC(int nVars) {
  int i,j,k;
  char *offset;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  double x[3];
  double x0All[3][3]={{0.0,0.0,0.0},{0.0,-1.5,-1.5},{0.0,1.5,1.5}};
  double x0l[3];

  for (int s=0;s<nVars;s++) {
    for (i=0;i<3;i++) x0l[i]=x0All[s][i];

    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      for (k=0;k<_BLOCK_CELLS_Z_+1;k++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (i=0;i<_BLOCK_CELLS_X_+1;i++) {
        node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

        offset=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

        x[0]=node->xmin[0]+(i*(node->xmax[0]-node->xmin[0]))/_BLOCK_CELLS_X_;
        x[1]=node->xmin[1]+(j*(node->xmax[1]-node->xmin[1]))/_BLOCK_CELLS_Y_;
        x[2]=node->xmin[2]+(k*(node->xmax[2]-node->xmin[2]))/_BLOCK_CELLS_Z_;

/*
        double r=Background[s],x0=sqrt(pow(x[0]-x0l[0],2)+pow(x[1]-x0l[1],2)+pow(x[2]-x0l[2],2));

        if (x0<1.0) r=Background[s]*(1.0-pow(x0-1.0,2));
*/

double r=Background[s];
if (x[1]*x[1]+x[2]*x[2]<1.0) if (fabs(x[0]-1.0)<0.5) {
r-=(0.5-fabs(x[0]-1.0))/1.0*Background[s];
}

r=0.0;
if ((0.0<=x[0])&&(x[0]<=2.0)) r=sin(x[0]*Pi/2.0);

        ((double*)(offset+CurrentCornerNodeOffset))[s]=r;
        ((double*)(offset+NextCornerNodeOffset))[s]=r;

      }
    }
  }

  //exchange data in the 'ghost' blocks
//  PIC::BC::ExternalBoundary::Periodic::UpdateData();
}

int main(int argc,char **argv) {
  PIC::InitMPI();
  PIC::Init_BeforeParser();


  int RelativeOffset=0;
  
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  printf("non-uniform mesh!\n");
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_
  printf("uniform mesh!\n");
#endif


  CurrentCenterNodeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=nVars*sizeof(double);
  NextCenterNodeOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=nVars*sizeof(double);

  CurrentCornerNodeOffset=PIC::Mesh::cDataCornerNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=nVars*sizeof(double);
  NextCornerNodeOffset=PIC::Mesh::cDataCornerNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCornerNode::totalAssociatedDataLength+=nVars*sizeof(double);
  
  PIC::Mesh::mesh.GetCenterNodesInterpolationCoefficients=PIC::Mesh::GetCenterNodesInterpolationCoefficients;
  //seed the random number generator
  rnd_seed(100);

  //generate mesh or read from file
  char mesh[_MAX_STRING_LENGTH_PIC_]="none";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  sprintf(mesh,"amr.sig=%s.mesh.bin","test_mesh");

  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::BC::ExternalBoundary::Periodic::Init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.memoryAllocationReport();

  //generate mesh or read from file
  bool NewMeshGeneratedFlag=false;

  char fullname[STRING_LENGTH];
  sprintf(fullname,"%s/%s",PIC::UserModelInputDataPath,mesh);

  FILE *fmesh=NULL;

  fmesh=fopen(fullname,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(fullname);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh.ThisThread==0) {
       PIC::Mesh::mesh.buildMesh();
       PIC::Mesh::mesh.saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }
  }


  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  

  //PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();
  PIC::Mesh::mesh.InitCellMeasure();

  PIC::Init_AfterParser();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  if (PIC::ThisThread==0) printf("test1\n");
  PIC::Mesh::mesh.outputMeshTECPLOT("mesh_test.dat");

  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable();
  

  double v[10][3]={{-1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,-1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,1.0},{0.0,0.0,-1.0},{-1,-1,0},{-1,-1,-1},{1,1,1},{0.8,0.8,-1}};
  double xparticle[10][3]={{-0.8,0,0},{0.0,1.9,0.0},{0.0,-2.5,0.0},{2.9,0.0,0.0},{0.0,0.0,2.9},{0.0,0.0,-3.5},{-0.9,-2.9,-3},{-0.9,-2.9,-3.9},{2.9,1.9,2.9},{2.9,1.9,2.9}};
  int s,i,j,k;


  int parSize=10;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  long int newParticle;

  if (PIC::ThisThread==0) printf("test2\n");
 
  // PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(0);
  //PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(1);
  
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1.0);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(1,1.0);
  
  
 
  
  PIC::Mesh::AddVaraibleListFunction(PrintCenterNodeVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintCenterNodeData);
  PIC::Mesh::InterpolateCenterNode.push_back(InterpolateCenterNode);

  PIC::Mesh::PrintVariableListCornerNode.push_back(PrintCornerNodeVariableList);
  PIC::Mesh::PrintDataCornerNode.push_back(PrintCornerNodeData);

  PIC::DomainBlockDecomposition::UpdateBlockTable();

  //solve the transport equation
  //set the initial conditions for the transport equation
  SetIC(1);
  PIC::Mesh::mesh.outputMeshDataTECPLOT("ic.dat",0);



  //transport test
  if (false) for (int niter=0;niter<100;niter++) {
      for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
        bool flag=false;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
        int i,j,k;

        for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
         flag=true;
         break;
        }

        if (flag==true) continue;

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          char *source,*target;

           source=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();  
           target=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

           *(double*)(target+NextCornerNodeOffset)=*(double*)(source+CurrentCornerNodeOffset);
        } 
      }

PIC::BC::ExternalBoundary::Periodic::UpdateData();

      int t=CurrentCornerNodeOffset;
      CurrentCornerNodeOffset=NextCornerNodeOffset;
      NextCornerNodeOffset=t;

      char fname[100];
      sprintf(fname,"test.out=%i.dat",niter);
      if (niter%5==0) PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,0);
    }



 PIC::BC::ExternalBoundary::Periodic::UpdateData();

  TransportEquation::BuildMatrix();

  for (int niter=0;niter<200;niter++) {

//PIC::Mesh::mesh.outputMeshDataTECPLOT("1.dat",0);

    TransportEquation::TimeStep();

//PIC::Mesh::mesh.outputMeshDataTECPLOT("2.dat",0);

    int t=CurrentCornerNodeOffset;
    CurrentCornerNodeOffset=NextCornerNodeOffset;
    NextCornerNodeOffset=t;

///PIC::Mesh::mesh.outputMeshDataTECPLOT("3.dat",0);

 //   PIC::BC::ExternalBoundary::Periodic::UpdateData();
//    SolverTransportEquation.UpdateRhs(UpdateRhsTransportEquation);

    char fname[100];
    sprintf(fname,"Transport.out=%i.dat",niter);
    if (niter%10==0) PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,0);
  }
 

  // countNumbers();

  
  for (int iPar=0;iPar<parSize; iPar++ ){
    newNode=PIC::Mesh::mesh.findTreeNode(xparticle[iPar]);
    
    if (newNode->Thread==PIC::ThisThread) {
      PIC::Mesh::mesh.fingCellIndex(xparticle[iPar],i,j,k,newNode);
      
      newParticle=PIC::ParticleBuffer::GetNewParticle(newNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      
      PIC::ParticleBuffer::SetV(v[iPar],newParticle);
      PIC::ParticleBuffer::SetX(xparticle[iPar],newParticle);
      PIC::ParticleBuffer::SetI(0,newParticle);
    }
  }
  

  MPI_Finalize();
  cout << "End of the run" << endl;

  return EXIT_SUCCESS;
}
