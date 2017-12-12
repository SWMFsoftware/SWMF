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

#include "stencil.h"

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"

#define _UNIFORM_MESH_ 1
#define _NONUNIFORM_MESH_ 2

#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif

double xmin[3]={-1.0,-3.0,-4.0};
double xmax[3]={3.0,2.0,3.0};

int nVars=3; //number of variables in center associated data
double Background[3]={100.0,-20.0,10.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,1> Solver1;
cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,7,7> SolverTimeDependent;

void solver1_matvec(double* VecIn, double * VecOut, int n){  
  Solver1.MultiplyVector(VecOut,VecIn,n);
}

void SolverTimeDependent_matvec(double* VecIn, double * VecOut, int n){
  SolverTimeDependent.MultiplyVector(VecOut,VecIn,n);
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
  for (int i=0;i<nVars; i++) fprintf(fout,", \"test wave (corner node)[i=%i]\"",i);
}

void PrintCornerNodeData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CornerNodeThread,PIC::Mesh::cDataCornerNode *CornerNode) {
  int idim;
  double t;


  for (int i=0;i<nVars;i++) {
    if (pipe->ThisThread==CornerNodeThread) {
      t= ((double*)(CornerNode->GetAssociatedDataBufferPointer()+CurrentCornerNodeOffset))[i];
    }

    if (pipe->ThisThread==0) {
      if (CornerNodeThread!=0) pipe->recv(t,CornerNodeThread);
      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
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
  double res =2;
#endif

#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  double highRes = dist/32.0, lowRes= dist/2.0;     
  double res =(5-1)/dist*(x[0]-xmin[0])+1;  
#endif

  res=dist/pow(2,res);
  
  return res/2.0;
}
                       
double waveFunction(double * xLocation, double * waveCenter, double * k, double waveWidth){
  // k is normalized wave number
  double sum=0;
  for (int iDim=0; iDim<3; iDim++){
    sum+=k[iDim]*(xLocation[iDim]-waveCenter[iDim]);
  }
  if (abs(sum)>waveWidth/2.) return 0;

  return 10*cos(sum*Pi/waveWidth);
}

int nExternalBlock=0;


// draw line and count numbers of ghost cells and true cells along the line 
void countNumbers(){
 
  int dataPoints=1024;
  double startPoint[3]={PIC::BC::ExternalBoundary::Periodic::xminDomain[0],0.0,0.0};
  double endPoint[3]={PIC::BC::ExternalBoundary::Periodic::xmaxDomain[0],0.0,0.0};
  double dl[3];
  dl[0]=(PIC::BC::ExternalBoundary::Periodic::xmaxDomain[0]-PIC::BC::ExternalBoundary::Periodic::xminDomain[0])/((double)dataPoints-1)/2.;
  dl[1]=0.0;
  dl[2]=0.0;
  double x[3];
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;  
  int positiveTen[]={0,0}, negativeTen[]={0,0};
  
  for (int i=0;i<3;i++) x[i]=startPoint[i]+1e-3;

  for (int iPoint=0;iPoint<dataPoints;iPoint++){
    printf("x:%f,%f,%f\n",x[0],x[1],x[2]);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
    newNode=PIC::Mesh::mesh.findTreeNode(x);
    
    int i,j,k;
    PIC::Mesh::mesh.fingCellIndex(x,i,j,k,newNode);
    printf("i,j,k:%d,%d,%d\n",i,j,k);
    int nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    if ((CenterNode=newNode->block->GetCenterNode(nd))==NULL) exit(__LINE__,__FILE__,"Error: not in the domain");

    offset=CenterNode->GetAssociatedDataBufferPointer();
    
    if (abs(((double*)(offset+CurrentCenterNodeOffset))[0]-10)<1e-3) positiveTen[0]++;
    if (abs(((double*)(offset+CurrentCenterNodeOffset))[0]+10)<1e-3) negativeTen[0]++;

    for (int iDim=0;iDim<3;iDim++) x[iDim] += dl[iDim];
  }

  printf("finished half\n");
   
  for (int i=0;i<3;i++) x[i]=endPoint[i]-1e-3;
  
  for (int iPoint=0;iPoint<dataPoints;iPoint++){
    printf("x:%f,%f,%f\n",x[0],x[1],x[2]);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
    newNode=PIC::Mesh::mesh.findTreeNode(x);
    
    int i,j,k;
    PIC::Mesh::mesh.fingCellIndex(x,i,j,k,newNode);

    int nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    if ((CenterNode=newNode->block->GetCenterNode(nd))==NULL) exit(__LINE__,__FILE__,"Error: not in the domain");

    offset=CenterNode->GetAssociatedDataBufferPointer();
    
    if (abs(((double*)(offset+CurrentCenterNodeOffset))[0]-10)<1e-3) positiveTen[1]++;
    if (abs(((double*)(offset+CurrentCenterNodeOffset))[0]+10)<1e-3) negativeTen[1]++;

    for (int iDim=0;iDim<3;iDim++) x[iDim] -= dl[iDim];
  }

  printf("positiveTen:%d,%d\n",positiveTen[0],positiveTen[1]);
  printf("negativeTen:%d,%d\n",negativeTen[0],negativeTen[1]);

}

void InitCenterData(int ipass,int nVars, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode, double * waveCenter, double * waveNumber, double waveWidth) {
 
  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd,idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  PIC::Mesh::cDataCornerNode *CornerNode;

  char *offset;
  /*
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
  */
  const int iMin=0,iMax=_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=_BLOCK_CELLS_Y_-1;                                                 
  const int kMin=0,kMax=_BLOCK_CELLS_Z_-1; 

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double xLOCAL[3];
    int ii,jj;
    /*
    double waveCenter[3]={1,0,0};
    double waveNumber[3]={1,0,0};
    double waveWidth=2.0;
    */
    bool externalBlockFlag=false;
    if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
      externalBlockFlag=true; 
      nExternalBlock++;
    }
    
    if ((startNode->Thread==PIC::ThisThread)&&(startNode->block!=NULL)) {
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        //the interpolation location
        xLOCAL[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
        xLOCAL[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
        xLOCAL[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);
	    
        //locate the cell
        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
        //  offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;
        offset=CenterNode->GetAssociatedDataBufferPointer();
	    
        if (externalBlockFlag)  for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCenterNodeOffset))[iVar] = -10;
	    
        if (!externalBlockFlag) for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCenterNodeOffset))[iVar] = waveFunction(xLOCAL, waveCenter, waveNumber, waveWidth); 
        // if (!externalBlockFlag) for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCenterNodeOffset))[iVar] = 10;
															   
        // *((double*)(offset+dataOffset+iVar*sizeof(double))) = waveFunction(xLOCAL, waveCenter, waveNumber, waveWidth);
      }//if (startNode->block!=NULL) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) 

      for (k=kMin;k<=kMax+1;k++) for (j=jMin;j<=jMax+1;j++) for (i=iMin;i<=iMax+1;i++) {
        //the interpolation location
        xLOCAL[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*i;
        xLOCAL[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*j;
        xLOCAL[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*k;

        //locate the cell
        nd=PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k);
        if ((CornerNode=startNode->block->GetCornerNode(nd))==NULL) continue;
        //  offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;
        offset=CornerNode->GetAssociatedDataBufferPointer();

        if (ipass=0) if (externalBlockFlag)  for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCornerNodeOffset))[iVar] = -10;

        if (!externalBlockFlag) for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCornerNodeOffset))[iVar] = waveFunction(xLOCAL, waveCenter, waveNumber, waveWidth);
        // if (!externalBlockFlag) for (int iVar=0;iVar<nVars;iVar++) ((double*)(offset+CurrentCornerNodeOffset))[iVar] = 10;

        // *((double*)(offset+dataOffset+iVar*sizeof(double))) = waveFunction(xLOCAL, waveCenter, waveNumber, waveWidth);
      }//if (startNode->block!=NULL) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)
    }
  }//if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) 
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitCenterData(ipass,nVars,startNode->downNode[nDownNode],waveCenter, waveNumber, waveWidth);
  }
  
}

void PropagateCenterData(double * v, int nVars, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
 
  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd,idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  PIC::Mesh::cDataCornerNode *CornerNode;
  char *offset;
  
  /*
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
  */
  const int iMin=0,iMax=_BLOCK_CELLS_X_-1;                                                 
  const int jMin=0,jMax=_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=_BLOCK_CELLS_Z_-1; 

  
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    // only propagate waves in the user defined region
    if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode)==_EXTERNAL_BOUNDARY_BLOCK_) return;
    
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double xLOCAL[3];
    int ii,jj;
   

    if ((startNode->Thread==PIC::ThisThread)&&(startNode->block!=NULL)) {
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        //the interpolation location
        xLOCAL[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
        xLOCAL[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
        xLOCAL[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);
	    
        //locate the cell
        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;

        //  offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;
        double xInterpolate[3];
        //find the location of the corresponding point at previous time step
        for (int iDim=0; iDim<3; iDim++) xInterpolate[iDim]=xLOCAL[iDim]-v[iDim]*PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
	
	      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *tempNode=PIC::Mesh::mesh.findTreeNode(xInterpolate,startNode);

	   
        //PIC::CPLR::InitInterpolationStencil(xInterpolate,startNode);
        //PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(xInterpolate,startNode);
        PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInterpolate,tempNode); 
        int iStencil;
        double currentInterpolation[nVars];
        for (int iVar=0; iVar<nVars; iVar++) currentInterpolation[iVar]=0.0;
        PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

        memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
	    
        for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
          for (int iVar=0; iVar<nVars; iVar++){
            currentInterpolation[iVar]+=((double*)(CurrentCenterNodeOffset+Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()))[iVar]*Stencil.Weight[iStencil];
          }
        }
	    
        offset=CenterNode->GetAssociatedDataBufferPointer();
	    
        for (int iVar=0;iVar<nVars;iVar++)
          ((double*)(offset+NextCenterNodeOffset))[iVar] = currentInterpolation[iVar];

       }//if (startNode->block!=NULL) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) 


      for (k=kMin;k<=kMax+1;k++) for (j=jMin;j<=jMax+1;j++) for (i=iMin;i<=iMax+1;i++) {
        //the interpolation location
        xLOCAL[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*i;
        xLOCAL[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*j;
        xLOCAL[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*k;

        //locate the cell
        nd=PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k);
        if ((CornerNode=startNode->block->GetCornerNode(nd))==NULL) continue;

        //  offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;
        double xInterpolate[3];
        //find the location of the corresponding point at previous time step
        for (int iDim=0; iDim<3; iDim++) xInterpolate[iDim]=xLOCAL[iDim]-v[iDim]*PIC::ParticleWeightTimeStep::GlobalTimeStep[0];


        //PIC::CPLR::InitInterpolationStencil(xInterpolate,startNode);
        //PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(xInterpolate,startNode);
        PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInterpolate,startNode);
        int iStencil;
        double currentInterpolation[nVars];
        for (int iVar=0; iVar<nVars; iVar++) currentInterpolation[iVar]=0.0;
        PIC::InterpolationRoutines::CellCentered::cStencil Stencil;

        memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));

        for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
          for (int iVar=0; iVar<nVars; iVar++){
            currentInterpolation[iVar]+=((double*)(CurrentCornerNodeOffset+Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()))[iVar]*Stencil.Weight[iStencil];
          }
        }

        offset=CornerNode->GetAssociatedDataBufferPointer();

        for (int iVar=0;iVar<nVars;iVar++)
          ((double*)(offset+NextCornerNodeOffset))[iVar] = currentInterpolation[iVar];

       }//if (startNode->block!=NULL) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)

    }
  }//if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) 
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) PropagateCenterData(v,nVars,startNode->downNode[nDownNode]);
  }

 
}


void test_wave(int iTest, int nVars, double * waveCenter, double * waveNumber, double waveWidth, double *vProp){
    //init center wave data
  InitCenterData(0,nVars,PIC::Mesh::mesh.rootTree, waveCenter, waveNumber, waveWidth);
  InitCenterData(1,nVars,PIC::Mesh::mesh.rootTree, waveCenter, waveNumber, waveWidth);
 
  

  PIC::Mesh::mesh.ParallelBlockDataExchange();
   
  char wave_init_fname[STRING_LENGTH];
  sprintf(wave_init_fname,"test%d-wave-init.dat",iTest);   

  PIC::Mesh::mesh.outputMeshDataTECPLOT(wave_init_fname,0);

  PIC::BC::ExternalBoundary::Periodic::UpdateData();
 
  sprintf(wave_init_fname,"test%d-wave-init-updated.dat",iTest);
  PIC::Mesh::mesh.outputMeshDataTECPLOT(wave_init_fname,0);
  
  
  for (int iter=0; iter<1; iter++) {
    PropagateCenterData(vProp,nVars,PIC::Mesh::mesh.rootTree);

    int temp=CurrentCenterNodeOffset;
    CurrentCenterNodeOffset=NextCenterNodeOffset;
    NextCenterNodeOffset=temp;

    temp=CurrentCornerNodeOffset;
    CurrentCornerNodeOffset=NextCornerNodeOffset;
    NextCornerNodeOffset=temp;

    PIC::BC::ExternalBoundary::Periodic::UpdateData();

    if (iter%5==0){
      char wave_fname[STRING_LENGTH];
      sprintf(wave_fname,"test%d-wave%d.dat",iTest,iter);
      PIC::Mesh::mesh.outputMeshDataTECPLOT(wave_fname,0);
    }
    
    PIC::TimeStep();
  }

  


}


//create the stencil for the linear equaton solver test
void GetTestStencil(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,1>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,1>::cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  /*
 MatrixRowNonZeroElementTable[0].i=i,MatrixRowNonZeroElementTable[0].j=j,MatrixRowNonZeroElementTable[0].k=k,MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
NonZeroElementsFound=1;
rhs=100.0;
MatrixRowNonZeroElementTable[0].iVar=0;
return;
  */
  
  //check whether the point is located at the boundary
  //boundary cell at the larger end is absorbed
  int index[3] = {i,j,k};
  int nCell[3] = {_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  // set  bc for the lower boundary
  /*
  for (int iDim=0;iDim<3;iDim++) if ((node->GetNeibFace(iDim*2,0,0)==NULL) && index[iDim]==0 ) {
    //the point is located at the boundary of the domain
    MatrixRowNonZeroElementTable[0].i=i,MatrixRowNonZeroElementTable[0].j=j,MatrixRowNonZeroElementTable[0].k=k,MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=0;

    rhs=100.0;   // boundary value
    NonZeroElementsFound=1;

    return;
  }
  */
 
  int iElement = 0;
  int nDirBoundary =0;
  //the point is within the domain, fill matrix for neighbor cells
  rhs =0.0;
  for (int ii=0;ii<3; ii++){
    int jMax =2;
    int sign=1;
    // set bc for the neihbor cell of higher boundary
    if (index[ii]== nCell[ii]-1 && node->GetNeibFace(2*ii+1,0,0)==NULL) {
      jMax = 1;
      sign = 1; 
      nDirBoundary++;
      rhs -= 100; //boundary condition
    }

    if (index[ii]== 0 && node->GetNeibFace(2*ii,0,0)==NULL) {
      jMax = 1;
      sign = -1; 
      nDirBoundary++;
      rhs -= 100; //boundary condition
    }
   
    int addition[2]={-1,1};
    for (int jj=0;jj<jMax;jj++){             
      MatrixRowNonZeroElementTable[iElement].i=(ii!=0)?i:i+sign*addition[jj];
      MatrixRowNonZeroElementTable[iElement].j=(ii!=1)?j:j+sign*addition[jj];
      MatrixRowNonZeroElementTable[iElement].k=(ii!=2)?k:k+sign*addition[jj];
      MatrixRowNonZeroElementTable[iElement].MatrixElementValue=1.0;
      MatrixRowNonZeroElementTable[iElement].iVar=0;
      iElement++;
    }
  }
  // for the cell itself
  MatrixRowNonZeroElementTable[iElement].i=i;
  MatrixRowNonZeroElementTable[iElement].j=j;
  MatrixRowNonZeroElementTable[iElement].k=k;
  MatrixRowNonZeroElementTable[iElement].MatrixElementValue=-6.0;
  MatrixRowNonZeroElementTable[iElement].iVar=0;
  iElement++;
  
  NonZeroElementsFound=iElement;
}


//create the stencil for the linear equaton solver test
void GetTestStencilTimeDependent(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,7>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,1,7,7>::cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  /*
 MatrixRowNonZeroElementTable[0].i=i,MatrixRowNonZeroElementTable[0].j=j,MatrixRowNonZeroElementTable[0].k=k,MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
NonZeroElementsFound=1;
rhs=100.0;
MatrixRowNonZeroElementTable[0].iVar=0;
return;
  */

  static char *StateVectorBS=NULL;
  static double ValueBC=100;

  double c=0.05;

  if (StateVectorBS==NULL) {
    StateVectorBS=new char[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

    *((double*)(StateVectorBS+CurrentCornerNodeOffset))=ValueBC;
    *((double*)(StateVectorBS+NextCornerNodeOffset))=ValueBC;
  }

  //populate the spencil
  if ( ((i==0)&&(node->GetNeibFace(0,0,0)==NULL)) || ((j==0)&&(node->GetNeibFace(2,0,0)==NULL)) || ((k==0)&&(node->GetNeibFace(4,0,0)==NULL)) ) {
    //the point is at the left boundary
    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=0;

    rhs=0.0;
    NonZeroElementsFound=1;
    RhsSupportLength=0;
  }
  else if ( ((i==_BLOCK_CELLS_X_-1)&&(node->GetNeibFace(1,0,0)==NULL)) || ((j==_BLOCK_CELLS_Y_-1)&&(node->GetNeibFace(3,0,0)==NULL)) || ((k==_BLOCK_CELLS_Z_-1)&&(node->GetNeibFace(5,0,0)==NULL)) ) {
    //the point is at the right boundary
    //x-direction

    int ElementCounter=0;

    if (i+1<_BLOCK_CELLS_X_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i+1;
      MatrixRowNonZeroElementTable[ElementCounter].j=j;
      MatrixRowNonZeroElementTable[ElementCounter].k=k;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[0].Coefficient=-c;
      RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[0].Coefficient=-c;
      RhsSupportTable[0].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i-1;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[1].Coefficient=-c;
    RhsSupportTable[1].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    //y-direction
    if (j+1<_BLOCK_CELLS_Y_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i;
      MatrixRowNonZeroElementTable[ElementCounter].j=j+1;
      MatrixRowNonZeroElementTable[ElementCounter].k=k;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[2].Coefficient=-c;
      RhsSupportTable[2].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[2].Coefficient=-c;
      RhsSupportTable[2].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j-1;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[3].Coefficient=-c;
    RhsSupportTable[3].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();

    //z-direction
    if (k+1<_BLOCK_CELLS_Z_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i;
      MatrixRowNonZeroElementTable[ElementCounter].j=j;
      MatrixRowNonZeroElementTable[ElementCounter].k=k+1;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[4].Coefficient=-c;
      RhsSupportTable[4].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[4].Coefficient=-c;
      RhsSupportTable[4].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k-1;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[5].Coefficient=-c;
    RhsSupportTable[5].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();

    //self
    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=(1.0+6.0*c);
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[6].Coefficient=6.0*c;
    RhsSupportTable[6].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=7;
    NonZeroElementsFound=ElementCounter;
  }
  else {
    //the point is inside the domain
    //x-direction
    MatrixRowNonZeroElementTable[0].i=i+1;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[0].iVar=0;

    RhsSupportTable[0].Coefficient=-c;
    RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[1].i=i-1;
    MatrixRowNonZeroElementTable[1].j=j;
    MatrixRowNonZeroElementTable[1].k=k;
    MatrixRowNonZeroElementTable[1].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[1].iVar=0;

    RhsSupportTable[1].Coefficient=-c;
    RhsSupportTable[1].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    //y-direction
    MatrixRowNonZeroElementTable[2].i=i;
    MatrixRowNonZeroElementTable[2].j=j+1;
    MatrixRowNonZeroElementTable[2].k=k;
    MatrixRowNonZeroElementTable[2].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[2].iVar=0;

    RhsSupportTable[2].Coefficient=-c;
    RhsSupportTable[2].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[3].i=i;
    MatrixRowNonZeroElementTable[3].j=j-1;
    MatrixRowNonZeroElementTable[3].k=k;
    MatrixRowNonZeroElementTable[3].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[3].iVar=0;

    RhsSupportTable[3].Coefficient=-c;
    RhsSupportTable[3].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();

    //z-direction
    MatrixRowNonZeroElementTable[4].i=i;
    MatrixRowNonZeroElementTable[4].j=j;
    MatrixRowNonZeroElementTable[4].k=k+1;
    MatrixRowNonZeroElementTable[4].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[4].iVar=0;

    RhsSupportTable[4].Coefficient=-c;
    RhsSupportTable[4].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[5].i=i;
    MatrixRowNonZeroElementTable[5].j=j;
    MatrixRowNonZeroElementTable[5].k=k-1;
    MatrixRowNonZeroElementTable[5].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[5].iVar=0;

    RhsSupportTable[5].Coefficient=-c;
    RhsSupportTable[5].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();

    //self
    MatrixRowNonZeroElementTable[6].i=i;
    MatrixRowNonZeroElementTable[6].j=j;
    MatrixRowNonZeroElementTable[6].k=k;
    MatrixRowNonZeroElementTable[6].MatrixElementValue=1.0+6.0*c;
    MatrixRowNonZeroElementTable[6].iVar=0;

    RhsSupportTable[6].Coefficient=6.0*c;
    RhsSupportTable[6].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=7;
    NonZeroElementsFound=7;
  }

  //set the RHS value
  if (RhsSupportLength!=0) for (i=0,rhs=0.0;i<RhsSupportLength;i++) rhs+=*((double*)(RhsSupportTable[i].CornerNodeAssociatedDataPointer+CurrentCenterNodeOffset))*RhsSupportTable[i].Coefficient;
}


void GetTestStencilTimeDependent1(int i,int j,int k,int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,7,7>::cMatrixRowNonZeroElementTable* MatrixRowNonZeroElementTable,int& NonZeroElementsFound,double& rhs,
cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,7,7>::cRhsSupportTable* RhsSupportTable,int &RhsSupportLength,
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  /*
 MatrixRowNonZeroElementTable[0].i=i,MatrixRowNonZeroElementTable[0].j=j,MatrixRowNonZeroElementTable[0].k=k,MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
NonZeroElementsFound=1;
rhs=100.0;
MatrixRowNonZeroElementTable[0].iVar=0;
return;
  */

  static char *StateVectorBS=NULL;
  static double ValueBC=Background[iVar];

  double c=0.05;

  if (StateVectorBS==NULL) {
    StateVectorBS=new char[PIC::Mesh::cDataCornerNode::totalAssociatedDataLength];

    *((double*)(StateVectorBS+CurrentCornerNodeOffset))=ValueBC;
    *((double*)(StateVectorBS+NextCornerNodeOffset))=ValueBC;
  }

  //build the stencil
  static bool InitFlag=false;
  static cStencil Stencil;
  static cStencil::cStencilElement* StencilTable;
  static int StencilTableLength=0;

  if (InitFlag==false) {
    InitFlag=true;

    cStencil dx,dy,dz;

    dx.add(1.0,1,0,0);
    dx.add(-1.0,0,0,0);

    dy.add(1.0,0,1,0);
    dy.add(-1.0,0,0,0);

    dz.add(1.0,0,0,1);
    dz.add(-1.0,0,0,0);

    Stencil=dx;
    Stencil.SubstractShifted(dx,-1,0,0);

    Stencil+=dy;
    Stencil.SubstractShifted(dy,0,-1,0);

    Stencil+=dz;
    Stencil.SubstractShifted(dz,0,0,-1);

    Stencil*=-c;
    Stencil.add(1.0,0,0,0);

    StencilTable=Stencil.GetStencil();
    StencilTableLength=Stencil.StencilLength;
  }


  //boundary
  if ( ((i==0)&&(node->GetNeibFace(0,0,0)==NULL)) ||
      ((j==0)&&(node->GetNeibFace(2,0,0)==NULL)) ||
      ((k==0)&&(node->GetNeibFace(4,0,0)==NULL)) ||
      ((i==_BLOCK_CELLS_X_)&&(node->GetNeibFace(1,0,0)==NULL)) ||
      ((j==_BLOCK_CELLS_Y_)&&(node->GetNeibFace(3,0,0)==NULL)) ||
      ((k==_BLOCK_CELLS_Z_)&&(node->GetNeibFace(5,0,0)==NULL)) ) {
    //the point is at the boundary
    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=iVar;

    rhs=Background[iVar];
    NonZeroElementsFound=1;
    RhsSupportLength=0;
  }
  else {
    for (int ii=0;ii<StencilTableLength;ii++) {
      MatrixRowNonZeroElementTable[ii].i=i+StencilTable[ii].i;
      MatrixRowNonZeroElementTable[ii].j=j+StencilTable[ii].j;
      MatrixRowNonZeroElementTable[ii].k=k+StencilTable[ii].k;
      MatrixRowNonZeroElementTable[ii].MatrixElementValue=StencilTable[ii].a;
      MatrixRowNonZeroElementTable[ii].iVar=iVar;
    }

    RhsSupportTable[0].Coefficient=1.0;
    RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=1;
    NonZeroElementsFound=StencilTableLength;

/*
    //the point is inside the domain
    //x-direction
    MatrixRowNonZeroElementTable[0].i=i+1;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[0].iVar=0;

    RhsSupportTable[0].Coefficient=0.0;
    RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[1].i=i-1;
    MatrixRowNonZeroElementTable[1].j=j;
    MatrixRowNonZeroElementTable[1].k=k;
    MatrixRowNonZeroElementTable[1].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[1].iVar=0;

    RhsSupportTable[1].Coefficient=0.0;
    RhsSupportTable[1].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    //y-direction
    MatrixRowNonZeroElementTable[2].i=i;
    MatrixRowNonZeroElementTable[2].j=j+1;
    MatrixRowNonZeroElementTable[2].k=k;
    MatrixRowNonZeroElementTable[2].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[2].iVar=0;

    RhsSupportTable[2].Coefficient=0.0;
    RhsSupportTable[2].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[3].i=i;
    MatrixRowNonZeroElementTable[3].j=j-1;
    MatrixRowNonZeroElementTable[3].k=k;
    MatrixRowNonZeroElementTable[3].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[3].iVar=0;

    RhsSupportTable[3].Coefficient=0.0;
    RhsSupportTable[3].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();

    //z-direction
    MatrixRowNonZeroElementTable[4].i=i;
    MatrixRowNonZeroElementTable[4].j=j;
    MatrixRowNonZeroElementTable[4].k=k+1;
    MatrixRowNonZeroElementTable[4].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[4].iVar=0;

    RhsSupportTable[4].Coefficient=0.0;
    RhsSupportTable[4].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[5].i=i;
    MatrixRowNonZeroElementTable[5].j=j;
    MatrixRowNonZeroElementTable[5].k=k-1;
    MatrixRowNonZeroElementTable[5].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[5].iVar=0;

    RhsSupportTable[5].Coefficient=0.0;
    RhsSupportTable[5].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();

    //self
    MatrixRowNonZeroElementTable[6].i=i;
    MatrixRowNonZeroElementTable[6].j=j;
    MatrixRowNonZeroElementTable[6].k=k;
    MatrixRowNonZeroElementTable[6].MatrixElementValue=1.0+6.0*c;
    MatrixRowNonZeroElementTable[6].iVar=0;

    RhsSupportTable[6].Coefficient=1.0;
    RhsSupportTable[6].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=7;
    NonZeroElementsFound=7;*/
  }

  //set the RHS value
  if (RhsSupportLength!=0) for (i=0,rhs=0.0;i<RhsSupportLength;i++) rhs+=((double*)(RhsSupportTable[i].CornerNodeAssociatedDataPointer+CurrentCenterNodeOffset))[iVar]*RhsSupportTable[i].Coefficient;


  return;

  /*
  //populate the spencil
  if ( ((i==0)&&(node->GetNeibFace(0,0,0)==NULL)) || ((j==0)&&(node->GetNeibFace(2,0,0)==NULL)) || ((k==0)&&(node->GetNeibFace(4,0,0)==NULL)) ) {
    //the point is at the left boundary
    MatrixRowNonZeroElementTable[0].i=i;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=1.0;
    MatrixRowNonZeroElementTable[0].iVar=0;

    rhs=ValueBC;
    NonZeroElementsFound=1;
    RhsSupportLength=0;
  }
  else if ( ((i==_BLOCK_CELLS_X_-1)&&(node->GetNeibFace(1,0,0)==NULL)) || ((j==_BLOCK_CELLS_Y_-1)&&(node->GetNeibFace(3,0,0)==NULL)) || ((k==_BLOCK_CELLS_Z_-1)&&(node->GetNeibFace(5,0,0)==NULL)) ) {
    //the point is at the right boundary
    //x-direction

    int ElementCounter=0;

    if (i+1<_BLOCK_CELLS_X_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i+1;
      MatrixRowNonZeroElementTable[ElementCounter].j=j;
      MatrixRowNonZeroElementTable[ElementCounter].k=k;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[0].Coefficient=0.0;
      RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[0].Coefficient=c;
      RhsSupportTable[0].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i-1;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[1].Coefficient=0;
    RhsSupportTable[1].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    //y-direction
    if (j+1<_BLOCK_CELLS_Y_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i;
      MatrixRowNonZeroElementTable[ElementCounter].j=j+1;
      MatrixRowNonZeroElementTable[ElementCounter].k=k;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[2].Coefficient=0.0;
      RhsSupportTable[2].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[2].Coefficient=c;
      RhsSupportTable[2].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j-1;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[3].Coefficient=0.0;
    RhsSupportTable[3].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();

    //z-direction
    if (k+1<_BLOCK_CELLS_Z_) {
      MatrixRowNonZeroElementTable[ElementCounter].i=i;
      MatrixRowNonZeroElementTable[ElementCounter].j=j;
      MatrixRowNonZeroElementTable[ElementCounter].k=k+1;
      MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
      MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

      ElementCounter++;

      RhsSupportTable[4].Coefficient=0.0;
      RhsSupportTable[4].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();
    }
    else {
      RhsSupportTable[4].Coefficient=c;
      RhsSupportTable[4].CornerNodeAssociatedDataPointer=StateVectorBS;
    }

    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k-1;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[5].Coefficient=0.0;
    RhsSupportTable[5].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();

    //self
    MatrixRowNonZeroElementTable[ElementCounter].i=i;
    MatrixRowNonZeroElementTable[ElementCounter].j=j;
    MatrixRowNonZeroElementTable[ElementCounter].k=k;
    MatrixRowNonZeroElementTable[ElementCounter].MatrixElementValue=(1.0+6.0*c);
    MatrixRowNonZeroElementTable[ElementCounter].iVar=0;

    ElementCounter++;

    RhsSupportTable[6].Coefficient=1.0;
    RhsSupportTable[6].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=7;
    NonZeroElementsFound=ElementCounter;
  }
  else {
    //the point is inside the domain
    //x-direction
    MatrixRowNonZeroElementTable[0].i=i+1;
    MatrixRowNonZeroElementTable[0].j=j;
    MatrixRowNonZeroElementTable[0].k=k;
    MatrixRowNonZeroElementTable[0].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[0].iVar=0;

    RhsSupportTable[0].Coefficient=0.0;
    RhsSupportTable[0].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i+1,j,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[1].i=i-1;
    MatrixRowNonZeroElementTable[1].j=j;
    MatrixRowNonZeroElementTable[1].k=k;
    MatrixRowNonZeroElementTable[1].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[1].iVar=0;

    RhsSupportTable[1].Coefficient=0.0;
    RhsSupportTable[1].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i-1,j,k))->GetAssociatedDataBufferPointer();

    //y-direction
    MatrixRowNonZeroElementTable[2].i=i;
    MatrixRowNonZeroElementTable[2].j=j+1;
    MatrixRowNonZeroElementTable[2].k=k;
    MatrixRowNonZeroElementTable[2].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[2].iVar=0;

    RhsSupportTable[2].Coefficient=0.0;
    RhsSupportTable[2].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j+1,k))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[3].i=i;
    MatrixRowNonZeroElementTable[3].j=j-1;
    MatrixRowNonZeroElementTable[3].k=k;
    MatrixRowNonZeroElementTable[3].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[3].iVar=0;

    RhsSupportTable[3].Coefficient=0.0;
    RhsSupportTable[3].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j-1,k))->GetAssociatedDataBufferPointer();

    //z-direction
    MatrixRowNonZeroElementTable[4].i=i;
    MatrixRowNonZeroElementTable[4].j=j;
    MatrixRowNonZeroElementTable[4].k=k+1;
    MatrixRowNonZeroElementTable[4].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[4].iVar=0;

    RhsSupportTable[4].Coefficient=0.0;
    RhsSupportTable[4].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k+1))->GetAssociatedDataBufferPointer();

    MatrixRowNonZeroElementTable[5].i=i;
    MatrixRowNonZeroElementTable[5].j=j;
    MatrixRowNonZeroElementTable[5].k=k-1;
    MatrixRowNonZeroElementTable[5].MatrixElementValue=-c;
    MatrixRowNonZeroElementTable[5].iVar=0;

    RhsSupportTable[5].Coefficient=0.0;
    RhsSupportTable[5].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k-1))->GetAssociatedDataBufferPointer();

    //self
    MatrixRowNonZeroElementTable[6].i=i;
    MatrixRowNonZeroElementTable[6].j=j;
    MatrixRowNonZeroElementTable[6].k=k;
    MatrixRowNonZeroElementTable[6].MatrixElementValue=1.0+6.0*c;
    MatrixRowNonZeroElementTable[6].iVar=0;

    RhsSupportTable[6].Coefficient=1.0;
    RhsSupportTable[6].CornerNodeAssociatedDataPointer=node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer();

    RhsSupportLength=7;
    NonZeroElementsFound=7;
  }

  //set the RHS value
  if (RhsSupportLength!=0) for (i=0,rhs=0.0;i<RhsSupportLength;i++) rhs+=*((double*)(RhsSupportTable[i].CornerNodeAssociatedDataPointer+CurrentCenterNodeOffset))*RhsSupportTable[i].Coefficient;*/
}

//update RHS vector
double UpdateRhs(int iVar,cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode,3,7,7>::cRhsSupportTable* RhsSupportTable,int RhsSupportLength) {
  int i;
  double res=0.0;

  for (i=0;i<RhsSupportLength;i++) res+=((double*)(RhsSupportTable[i].CornerNodeAssociatedDataPointer+CurrentCenterNodeOffset))[iVar]*RhsSupportTable[i].Coefficient;

  return res;
}


//set the initial guess for the unknowns' values
void SetInitialGuess(double *x,PIC::Mesh::cDataCornerNode* CornerNode) {
  for (int i=0;i<nVars;i++) x[i]=Background[i];
}

void ProcessFinalSolution(double *x,PIC::Mesh::cDataCornerNode* CornerNode) {
  char* offset=CornerNode->GetAssociatedDataBufferPointer();

  for (int i=0;i<nVars;i++) ((double*)(offset+NextCornerNodeOffset))[i]=x[i];
}

void ProcessFinalSolutionTimeDependent(double *x,PIC::Mesh::cDataCornerNode* CornerNode) {
  char *offset=CornerNode->GetAssociatedDataBufferPointer();

  ((double*)(offset+NextCornerNodeOffset))[0]+=x[0];
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

        double r=Background[s],x0=sqrt(pow(x[0]-x0l[0],2)+pow(x[1]-x0l[1],2)+pow(x[2]-x0l[2],2));

        if (x0<1.0) r=Background[s]*(1.0-pow(x0-1.0,2));

        ((double*)(offset+CurrentCornerNodeOffset))[s]=r;
        ((double*)(offset+NextCornerNodeOffset))[s]=r;

      }
    }
  }

  //exchange data in the 'ghost' blocks
  PIC::Mesh::mesh.ParallelBlockDataExchange();
}

int main(int argc,char **argv) {
  PIC::InitMPI();
  PIC::Init_BeforeParser();

//  int nVars=3; //number of variables in center associated data
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

  //set the file output
  PIC::Mesh::AddVaraibleListFunction(PrintCenterNodeVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintCenterNodeData);
  PIC::Mesh::InterpolateCenterNode.push_back(InterpolateCenterNode);

  PIC::Mesh::PrintVariableListCornerNode.push_back(PrintCornerNodeVariableList);
  PIC::Mesh::PrintDataCornerNode.push_back(PrintCornerNodeData);

  //generate mesh or read from file
  char mesh[_MAX_STRING_LENGTH_PIC_]="none";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  sprintf(mesh,"amr.sig=%s.mesh.bin","test_mesh");

  PIC::Mesh::mesh.AllowBlockAllocation=false;
  //PIC::BC::ExternalBoundary::Periodic::Init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
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

 
//  cLinearSystemCornerNode<PIC::Mesh::cDataCornerNode> Solver2(2,CurrentCornerNodeOffset+sizeof(double));

  PIC::DomainBlockDecomposition::UpdateBlockTable();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,1.0);



  SetIC(3);
  SolverTimeDependent.Reset();
//  SolverTimeDependent.BuildMatrix(GetTestStencilTimeDependent);

  SolverTimeDependent.BuildMatrix(GetTestStencilTimeDependent1);

  linear_solver_matvec_c = SolverTimeDependent_matvec;

  PIC::Mesh::mesh.outputMeshDataTECPLOT("ic.dat",0);

  for (int niter=0;niter<20;niter++) {
//    SolverTimeDependent.Solve(SetInitialGuess,ProcessFinalSolutionTimeDependent);

    SolverTimeDependent.Solve(SetInitialGuess,ProcessFinalSolution);

    int t=CurrentCornerNodeOffset;
    CurrentCornerNodeOffset=NextCornerNodeOffset;
    NextCornerNodeOffset=t;

    PIC::Mesh::mesh.ParallelBlockDataExchange();
    SolverTimeDependent.UpdateRhs(UpdateRhs);

    char fname[100];
    sprintf(fname,"out=%i.dat",niter);
    PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,0);
  }


/*  Solver1.Reset();
  Solver1.BuildMatrix(GetTestStencil);
  linear_solver_matvec_c = solver1_matvec;
  Solver1.Solve(SetInitialGuess,ProcessFinalSolution);*/

  //double *xx;
  //Solver1.ExchageIntermediateUnknownsData(xx);

  PIC::Init_AfterParser();
  PIC::Mover::Init();



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
  
//  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1.0);

  
  
 
  

 

  // countNumbers();
  PIC::Mesh::mesh.ParallelBlockDataExchange();
  
  char wave_init_fname[STRING_LENGTH];
  sprintf(wave_init_fname,"test-solver-init.dat");   

  PIC::Mesh::mesh.outputMeshDataTECPLOT(wave_init_fname,0);

 
  MPI_Finalize();
  cout << "End of the run" << endl;

  return EXIT_SUCCESS;
}
