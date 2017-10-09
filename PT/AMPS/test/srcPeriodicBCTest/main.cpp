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


#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif

double xmin[3]={-1.0,-3.0,-4.0};
double xmax[3]={3.0,2.0,3.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

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
    CharacteristicSpeed=1.0;
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
    if (spec==_H_SPEC_) CharacteristicSpeed*=30.0;
    if (spec==_O_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_H2_SPEC_) CharacteristicSpeed*=10.0;
    if (spec==_OH_SPEC_) CharacteristicSpeed*=5.0;
#endif
#endif

    CellSize=startNode->GetCharacteristicCellSize();
    return 0.3*CellSize/CharacteristicSpeed;
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

void InitCenterData(int ipass,int nVars, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
 
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
    double waveCenter[3]={1,0,0};
    double waveNumber[3]={1,0,0};
    double waveWidth=2.0;
    bool externalBlockFlag=false;
    if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode)==_EXTERNAL_BOUNDARY_BLOCK_) {
      externalBlockFlag=true; 
      nExternalBlock++;
    }
    
    if (startNode->block!=NULL) {
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
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitCenterData(ipass,nVars,startNode->downNode[nDownNode]);
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
    double waveCenter[3]={1,0,0};
    double waveNumber[3]={1,0,0};
    double waveWidth=2.0;

    if (startNode->block!=NULL) {
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

        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *tempNode=PIC::Mesh::mesh.findTreeNode(xInterpolate,startNode);
        PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInterpolate,tempNode);

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


int main(int argc,char **argv) {
  PIC::InitMPI();
  PIC::Init_BeforeParser();

  int nVars=1; //number of variables in center associated data
  int RelativeOffset=0;
  
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
  
  
  //init center wave data
  InitCenterData(0,nVars,PIC::Mesh::mesh.rootTree);
  InitCenterData(1,nVars,PIC::Mesh::mesh.rootTree);

  printf("nExternalBlock:%d\n",nExternalBlock);
 
  PIC::Mesh::AddVaraibleListFunction(PrintCenterNodeVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintCenterNodeData);
  PIC::Mesh::InterpolateCenterNode.push_back(InterpolateCenterNode);

  PIC::Mesh::PrintVariableListCornerNode.push_back(PrintCornerNodeVariableList);
  PIC::Mesh::PrintDataCornerNode.push_back(PrintCornerNodeData);


  // countNumbers();

  PIC::Mesh::mesh.ParallelBlockDataExchange();
  PIC::Mesh::mesh.outputMeshDataTECPLOT("wave-init.dat",0);

  PIC::BC::ExternalBoundary::Periodic::UpdateData();
 
  PIC::Mesh::mesh.outputMeshDataTECPLOT("wave-init-1.dat",0);
  
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
  
  printf("test3\n");
  
  double vPropagate[]={1.0,0.0,0.0};
  for (int iter=0; iter<51; iter++) {
    //  PIC::BC::ExternalBoundary::Periodic::UpdateData();
    
    PropagateCenterData(vPropagate,nVars,PIC::Mesh::mesh.rootTree);

    int temp=CurrentCenterNodeOffset;
    CurrentCenterNodeOffset=NextCenterNodeOffset;
    NextCenterNodeOffset=temp;

    temp=CurrentCornerNodeOffset;
    CurrentCornerNodeOffset=NextCornerNodeOffset;
    NextCornerNodeOffset=temp;

    PIC::BC::ExternalBoundary::Periodic::UpdateData();
//    PIC::Mesh::mesh.ParallelBlockDataExchange();

    if (iter%10==0){
      char wave_fname[STRING_LENGTH];
      sprintf(wave_fname,"wave%d.dat",iter);
      PIC::Mesh::mesh.outputMeshDataTECPLOT(wave_fname,0);
    }
    
    PIC::TimeStep();
  }
  
  MPI_Finalize();
  cout << "End of the run" << endl;

  return EXIT_SUCCESS;
}
