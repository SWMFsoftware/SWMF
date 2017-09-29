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



double xmin[3]={-1.0,-3.0,-4.0};
double xmax[3]={3.0,2.0,3.0};

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
  double highRes = dist/32.0, lowRes= dist/2.0;

  double res =(5-1)/dist*(x[0]-xmin[0])+1;
  
  res=dist/pow(2,res);
    
  return res;
}


int main(int argc,char **argv) {
  PIC::InitMPI();
  PIC::Init_BeforeParser();

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

  printf("test1\n");

PIC::Mesh::mesh.outputMeshTECPLOT("mesh_test.dat");


  /*
  PIC::Mesh::mesh.outputMeshTECPLOT("mesh_test.dat");
  
  int nExternalLeafBlks;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> ** GhostBlockArray=NewExternalLeafBlockArray(nExternalLeafBlks);
  
  printf("nExternalLeafBlks:%d\n",nExternalLeafBlks);

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> ** TrueBlockArray=CorrespondingTrueBlockArray(nExternalLeafBlks, xmin, xmax, extensionDx, GhostBlockArray);
  */
  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable();

  /*
  printf("GhostBlockArray[0]->xmin[0],GhostBlockArray[0]->xmin[1],GhostBlockArray[0]->xmin[2]:%f,%f,%f\n", GhostBlockArray[0]->xmin[0],GhostBlockArray[0]->xmin[1],GhostBlockArray[0]->xmin[2]);
  printf("GhostBlockArray[0]->xmax[0],GhostBlockArray[0]->xmax[1],GhostBlockArray[0]->xmax[2]:%f,%f,%f\n", GhostBlockArray[0]->xmax[0],GhostBlockArray[0]->xmax[1],GhostBlockArray[0]->xmax[2]);

  printf("TrueBlockArray[0]->xmin[0],TrueBlockArray[0]->xmin[1],TrueBlockArray[0]->xmin[2]:%f,%f,%f\n", TrueBlockArray[0]->xmin[0],TrueBlockArray[0]->xmin[1],TrueBlockArray[0]->xmin[2]);
  printf("TrueBlockArray[0]->xmax[0],TrueBlockArray[0]->xmax[1],TrueBlockArray[0]->xmax[2]:%f,%f,%f\n", TrueBlockArray[0]->xmax[0],TrueBlockArray[0]->xmax[1],TrueBlockArray[0]->xmax[2]);


  
  printf("GhostBlockArray[1001]->xmin[0],GhostBlockArray[1001]->xmin[1],GhostBlockArray[1001]->xmin[2]:%f,%f,%f\n", GhostBlockArray[1001]->xmin[0],GhostBlockArray[1001]->xmin[1],GhostBlockArray[1001]->xmin[2]);
  printf("GhostBlockArray[1001]->xmax[0],GhostBlockArray[1001]->xmax[1],GhostBlockArray[1001]->xmax[2]:%f,%f,%f\n", GhostBlockArray[1001]->xmax[0],GhostBlockArray[1001]->xmax[1],GhostBlockArray[1001]->xmax[2]);

  printf("TrueBlockArray[1001]->xmin[0],TrueBlockArray[1001]->xmin[1],TrueBlockArray[1001]->xmin[2]:%f,%f,%f\n", TrueBlockArray[1001]->xmin[0],TrueBlockArray[1001]->xmin[1],TrueBlockArray[1001]->xmin[2]);
  printf("TrueBlockArray[1001]->xmax[0],TrueBlockArray[1001]->xmax[1],TrueBlockArray[1001]->xmax[2]:%f,%f,%f\n", TrueBlockArray[1001]->xmax[0],TrueBlockArray[1001]->xmax[1],TrueBlockArray[1001]->xmax[2]);
  printf("diff:%f,%f,%f", abs(TrueBlockArray[1001]->xmax[0]-GhostBlockArray[1001]->xmax[0]), abs(TrueBlockArray[1001]->xmax[1]-GhostBlockArray[1001]->xmax[1]), abs(TrueBlockArray[1001]->xmax[2]-GhostBlockArray[1001]->xmax[2]));
  */
  //newParticle=PIC::ParticleBuffer::GetNewParticle();
  
  
  double v[10][3]={{-1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,-1.0,0.0},{1.0,0.0,0.0},{0.0,0.0,1.0},{0.0,0.0,-1.0},{-1,-1,0},{-1,-1,-1},{1,1,1},{0.8,0.8,-1}};
double xparticle[10][3]={{-0.8,0,0},{0.0,1.9,0.0},{0.0,-2.5,0.0},{2.9,0.0,0.0},{0.0,0.0,2.9},{0.0,0.0,-3.5},{-0.9,-2.9,-3},{-0.9,-2.9,-3.9},{2.9,1.9,2.9},{2.9,1.9,2.9}};
  int s,i,j,k;


  int parSize=10;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode[parSize];
  long int newParticle[parSize];

  printf("test2\n");
//  PIC::ParticleBuffer::Init(10000);
 
  for (int iPar=0;iPar<parSize; iPar++ ){
  newNode[iPar]=PIC::Mesh::mesh.findTreeNode(xparticle[iPar]);
  // if(newNode->Thread==PIC::ThisThread){    
  PIC::Mesh::mesh.fingCellIndex(xparticle[iPar],i,j,k,newNode[iPar]); 
  
  newParticle[iPar]=PIC::ParticleBuffer::GetNewParticle(newNode[iPar]->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
  
  PIC::ParticleBuffer::SetV(v[iPar],newParticle[iPar]);
  PIC::ParticleBuffer::SetX(xparticle[iPar],newParticle[iPar]);
  PIC::ParticleBuffer::SetI(0,newParticle[iPar]);
  }
    //  }


  /*  double xminUser[3], xmaxUser[3];
  for (int iDim=0; iDim<3; iDim++){
    xminUser[iDim]=xmin[iDim]+extensionDx[iDim];
    xmaxUser[iDim]=xmax[iDim]-extensionDx[iDim];
  }
  */
    int nMoved=0;
    for (int iter=0; iter<10; iter++){
      double vtemp[3],xtemp[3];

    /* if (nMoved!=0){
      findXLocation(xtemp,xminUser,xmaxUser);
      newNode=PIC::Mesh::mesh.findTreeNode(xtemp);
      PIC::Mesh::mesh.fingCellIndex(xtemp,i,j,k,newNode);
      newParticle=PIC::ParticleBuffer::GetNewParticle(newNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
    }
    */
      for (int iPar=0;iPar<parSize; iPar++ ){
	s=PIC::ParticleBuffer::GetI(newParticle[iPar]);
	PIC::ParticleBuffer::GetV(vtemp,newParticle[iPar]);
	PIC::ParticleBuffer::GetX(xtemp,newParticle[iPar]);
	
       
	printf("iPar:%d, iter:%d, v:%f,%f,%f\n",iPar,iter,vtemp[0],vtemp[1],vtemp[2]);
	printf("iPar:%d, iter:%d, x:%f,%f,%f\n",iPar,iter,xtemp[0],xtemp[1],xtemp[2]);
      }
      
  // nMoved=MoveParticleFromGhostBlock(TrueBlockArray, GhostBlockArray, nExternalLeafBlks);
  
  PIC::BC::ExternalBoundary::Periodic::UpdateData();
  PIC::TimeStep();
}
  
  MPI_Finalize();
  cout << "End of the run" << endl;

  return EXIT_SUCCESS;
}
