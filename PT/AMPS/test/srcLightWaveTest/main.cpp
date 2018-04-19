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

#include "PeriodicBCTest.dfn"


//for lapenta mover

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"


int nVars=3; //number of variables in center associated data
double Background[3]={100.0,-20.0,10.0};


//#define _UNIFORM_MESH_ 1
//#define _NONUNIFORM_MESH_ 2

#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif


double xmin[3]={-2.0,-2.0,-2.0};
double xmax[3]={2.0,2.0,2.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

int iCase;

void CleanParticles(){
  
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) if (node->block!=NULL) {
   
     long int *  FirstCellParticleTable=node->block->FirstCellParticleTable;
     if (FirstCellParticleTable==NULL) continue;
     for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
	 for (int i=0;i<_BLOCK_CELLS_X_;i++) {
	   long int * ptr=FirstCellParticleTable+(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));
	   while ((*ptr)!=-1) PIC::ParticleBuffer::DeleteParticle(*ptr,*ptr);

//////
/*
long int next,ptr=FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]; 

while (ptr!=-1) {

  next=PIC::ParticleBuffer::GetNext(ptr);
  PIC::ParticleBuffer::DeleteParticle(ptr);
  ptr=next; 
}

FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]=-1;
*/

/////	   
	 }   
       }
     }
     
  }

}

long int PrepopulateDomain(int spec,double NumberDensity,double Temperature) {
  int iCell,jCell,kCell;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;
  double Velocity[3];
  /*
  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  */
  //particle ejection parameters
  double ParticleWeight;//beta=PIC::MolecularData::GetMass(spec)/(2*Kbol*Temperature);

#ifndef DIM
#error ERROR: DIM is used but not defined
#endif
#ifndef DIM
#error ERROR: DIM is used but not defined
#endif
#ifndef DIM
#error ERROR: DIM is used but not defined
#endif
#if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

  //the boundaries of the block and middle point of the cell
  double *xmin,*xmax,*xMiddle;
  double x[3],v[3],anpart;
  int npart,idim;

  //  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
  //  {
        for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
	  //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	  BoundaryBlock=true;
	  break;
	}
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;


    // }
    //local copy of the block's cells
    int cellListLength=node->block->GetCenterNodeListLength();
    PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  
    memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    xmin=node->xmin,xmax=node->xmax;

    //particle stat weight
#ifndef _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
#error ERROR: _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_ is used but not defined
#endif
#ifndef _SIMULATION_PARTICLE_WEIGHT_MODE_
#error ERROR: _SIMULATION_PARTICLE_WEIGHT_MODE_ is used but not defined
#endif
    #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    #else
    ParticleWeight=node->block->GetLocalParticleWeight(spec);
    #endif

    for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      cell=cellList[nd];
      xMiddle=cell->GetX();

      //inject particles into the cell
      anpart=NumberDensity*cell->Measure/ParticleWeight;
      npart=(int)(anpart);
      //if (rnd()<anpart-npart) npart++;
      nLocalInjectedParticles+=npart;

      while (npart-->0) {
	/*
        x[0]=xMiddle[0]+(xmax[0]-xmin[0])/_BLOCK_CELLS_X_*(rnd()-0.5);
        x[1]=xMiddle[1]+(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_*(rnd()-0.5);
        x[2]=xMiddle[2]+(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_*(rnd()-0.5);	
	*/
	
	x[0]=xMiddle[0];
        x[1]=xMiddle[1];
        x[2]=xMiddle[2];
       

	if (spec==1){
	  Velocity[0]=0.05*sin(Pi*(x[0]+x[1])); // ocillate in x direction
	  Velocity[1]=0.05*sin(Pi*(x[0]+x[1]));
 
	}else if (spec==0){
	  Velocity[0]=0.0;
	  Velocity[1]=0.0;
	}
	
	//Velocity[1]=0.0;
	Velocity[2]=0.0;
        for (idim=0;idim<3;idim++) v[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())*(2*Kbol*Temperature)/PIC::MolecularData::GetMass(spec))+Velocity[idim];

        //initiate the new particle
        PIC::ParticleBuffer::InitiateParticle(x,v,NULL,&spec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
	
      }
      //end of the particle injection block
    }
  }

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  return nGlobalInjectedParticles;
}



void SetIC() {
  
    int i,j,k;
    char *offset;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    double cPi = 3.14159265;
    double waveNumber[3]={0.0,0.0,0.0};

    if (iCase!=1){
      waveNumber[0]=cPi/2;
    }else{
      waveNumber[0]=cPi/2;
      waveNumber[1]=cPi/2;     
    }
  
    double x[3];
   
    using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
    int nBreak=0;

    printf("User Set IC called\n");
    printf("Set IC case number:%d\n",iCase);
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      
      if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
	bool BoundaryBlock=false;
	
	for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
	    //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	    BoundaryBlock=true;
	    break;
	  }
	
	if (BoundaryBlock==true) continue;
      }
      

     
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
	    
	    PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh.getCornerNodeLocalNumber(i,j,k));
	    if (CornerNode!=NULL){
	      offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
	      
	      x[0]=node->xmin[0]+(i*(node->xmax[0]-node->xmin[0]))/_BLOCK_CELLS_X_;
	      x[1]=node->xmin[1]+(j*(node->xmax[1]-node->xmin[1]))/_BLOCK_CELLS_Y_;
	      x[2]=node->xmin[2]+(k*(node->xmax[2]-node->xmin[2]))/_BLOCK_CELLS_Z_;

	 
	      double E = 10*sin(waveNumber[0]*(x[0]-xmin[0])+waveNumber[1]*(x[1]-xmin[1])+waveNumber[2]*(x[2]-xmin[2]));
	     
	      //1e9 charge in species.input
	      if (_CURRENT_MODE_==_PIC_MODE_OFF_){ 
		if (iCase==1){
		((double*)(offset+CurrentEOffset))[ExOffsetIndex]=-E/sqrt(2);
		((double*)(offset+CurrentEOffset))[EyOffsetIndex]=E/sqrt(2);
		((double*)(offset+CurrentEOffset))[EzOffsetIndex]=0.0;
		}else{
		((double*)(offset+CurrentEOffset))[ExOffsetIndex]=0.0;                                                           
                ((double*)(offset+CurrentEOffset))[EyOffsetIndex]=E; 
		((double*)(offset+CurrentEOffset))[EzOffsetIndex]=E;      
		}
		
	      ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=0.0;
	      ((double*)(offset+OffsetE_HalfTimeStep))[EyOffsetIndex]=0.0;
	      ((double*)(offset+OffsetE_HalfTimeStep))[EzOffsetIndex]=0.0;
	      }else{
		
		((double*)(offset+CurrentEOffset))[ExOffsetIndex]=0.0;
		((double*)(offset+CurrentEOffset))[EyOffsetIndex]=0.0;
		((double*)(offset+CurrentEOffset))[EzOffsetIndex]=0.0;
		
		
		((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=0.0;
		((double*)(offset+OffsetE_HalfTimeStep))[EyOffsetIndex]=0.0;
		((double*)(offset+OffsetE_HalfTimeStep))[EzOffsetIndex]=0.0;		

	      }

	    

	    // ((double*)(offset+CurrentCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;


	    
	    //((double*)(offset+NextCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
	    }//
	  }//for (k=0;k<_BLOCK_CELLS_Z_+1;k++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (i=0;i<_BLOCK_CELLS_X_+1;i++) 
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
	     
	    PIC::Mesh::cDataCenterNode *CenterNode= node->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
	    if (CenterNode!=NULL){
	      offset=node->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

	      x[0]=node->xmin[0]+((i+0.5)*(node->xmax[0]-node->xmin[0]))/_BLOCK_CELLS_X_;
	      x[1]=node->xmin[1]+((j+0.5)*(node->xmax[1]-node->xmin[1]))/_BLOCK_CELLS_Y_;
	      x[2]=node->xmin[2]+((k+0.5)*(node->xmax[2]-node->xmin[2]))/_BLOCK_CELLS_Z_;

	      double B = 10*sin(waveNumber[0]*(x[0]-xmin[0])+waveNumber[1]*(x[1]-xmin[1])+waveNumber[2]*(x[2]-xmin[2]));
	    
	    
	      if (_CURRENT_MODE_==_PIC_MODE_OFF_){
	      
		if (iCase==1){
		  ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=0.0;
		  ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=0.0;
		  ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=B;
		}else{
		  ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=0.0;
		  ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=B;
		  ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=-B;
		}
		
		((double*)(offset+PrevBOffset))[BxOffsetIndex]=0.0;
		((double*)(offset+PrevBOffset))[ByOffsetIndex]=0.0;
		((double*)(offset+PrevBOffset))[BzOffsetIndex]=0.0;
	      
	      }else{
			      
		((double*)(offset+CurrentBOffset))[BxOffsetIndex]=0.0;
		((double*)(offset+CurrentBOffset))[ByOffsetIndex]=0.0;
		((double*)(offset+CurrentBOffset))[BzOffsetIndex]=0.0;
		
		
		((double*)(offset+PrevBOffset))[BxOffsetIndex]=0.0;
		((double*)(offset+PrevBOffset))[ByOffsetIndex]=0.0;
		((double*)(offset+PrevBOffset))[BzOffsetIndex]=0.0;
		
	      }
	      
	    }// if (CenterNode!=NULL)
	  }//for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) 
    }
   
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh.ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::Periodic::UpdateData();
      break;
    }
}


double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;


    CellSize=startNode->GetCharacteristicCellSize();
    //return 0.3*CellSize/CharacteristicSpeed;

    //return 0.05;
    return 0.2;
}


double BulletLocalResolution(double *x) {                                                                                           
  double dist = xmax[0]-xmin[0];

#ifndef _UNIFORM_MESH_
#error ERROR: _UNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_  
  double res = 3;
#endif

#ifndef _NONUNIFORM_MESH_
#error ERROR: _NONUNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  double highRes = dist/32.0, lowRes= dist/2.0;     
  double res =(5-1)/dist*(x[0]-xmin[0])+1;  
#endif

  res=dist/pow(2,res);
  
  return res;
}
                       

int main(int argc,char **argv) {
  PIC::InitMPI();
  PIC::Init_BeforeParser();


  int RelativeOffset=0;
  
#ifndef _NONUNIFORM_MESH_
#error ERROR: _NONUNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  printf("non-uniform mesh!\n");
#endif
#ifndef _UNIFORM_MESH_
#error ERROR: _UNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_
  printf("uniform mesh!\n");
#endif


#ifndef _PIC_MODE_ON_
#error ERROR: _PIC_MODE_ON_ is used but not defined
#endif
#ifndef _CURRENT_MODE_
#error ERROR: _CURRENT_MODE_ is used but not defined
#endif
#if _CURRENT_MODE_==_PIC_MODE_ON_
  printf("current on!\n");
#endif
#ifndef _PIC_MODE_OFF_
#error ERROR: _PIC_MODE_OFF_ is used but not defined
#endif
#ifndef _CURRENT_MODE_
#error ERROR: _CURRENT_MODE_ is used but not defined
#endif
#if _CURRENT_MODE_==_PIC_MODE_OFF_
  printf("current mode off!\n");
#endif



  //seed the random number generator
  rnd_seed(100);

  //generate mesh or read from file
  char mesh[_MAX_STRING_LENGTH_PIC_]="none";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  sprintf(mesh,"amr.sig=%s.mesh.bin","test_mesh");

  PIC::Mesh::mesh.AllowBlockAllocation=false;
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
  PIC::BC::ExternalBoundary::Periodic::Init(xmin,xmax,BulletLocalResolution);
  }else{
    PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
  }
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
  

  PIC::Mesh::initCellSamplingDataBuffer();

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
  
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable();
  }
  //-387.99e2
  double v[14][3]={{1.0, 1.0, 1.0},{-1.0,1.0, 1.0},{1.0,-1.0,1.0},{-1.0,-1.0,1.0},{1.0,1.0,-1.0},{-1.0,1.0,-1.0},{1.0,-1.0,-1.0},{-1.0,-1.0,-1.0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
  double xparticle[10][3]={{0.125,1.875,0.125},{0.125,1.875,0.125},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},{2.9,0.9,2.9},{2.9,0.9,2.9}};
  int s,i,j,k;
  int species[10]={0,1,0,1,0,1,0,1,0,1};
  
  int parSize;
  if (_CURRENT_MODE_==_PIC_MODE_OFF_){
    parSize=0;
  }else{
    parSize=2;
  }


  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* newNode;
  long int newParticle;

  if (PIC::ThisThread==0) printf("test2\n");
 
  // PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(0);
  //PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(1);
  
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,0.0625);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(1,0.0625);

  PIC::DomainBlockDecomposition::UpdateBlockTable();

  //solve the transport equation
  //set the initial conditions for the transport equation
  //  TransportEquation::SetIC(3);
 

  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh.ParallelBlockDataExchange();
      break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::Periodic::UpdateData();
      break;
  }
  //PIC::FieldSolver::Init(); 
  PIC::FieldSolver::Electromagnetic::ECSIM::SetIC=SetIC;
    
  int  totalIter,CaseNumber;
  if (_CURRENT_MODE_==_PIC_MODE_OFF_){
    CaseNumber=2;// one at 0 degree, one propagates at 45 degree
  }else{
    CaseNumber=2;//one particle energy test, one langmuir wave test
  }


  for (iCase=0;iCase<CaseNumber;iCase++){
    printf("iCase:%d\n",iCase);
    PIC::FieldSolver::Electromagnetic::ECSIM::Init_IC();

    if (_CURRENT_MODE_==_PIC_MODE_OFF_){
	totalIter = 4/PIC::FieldSolver::Electromagnetic::ECSIM::cDt;
    }else{
      totalIter = 7;
    }
    


     
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh.ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::Periodic::UpdateData();
      break;
    }
    PIC::Mesh::mesh.outputMeshDataTECPLOT("ic.dat",0);
  
 
    // countNumbers();

#ifndef _PIC_MODE_ON_
#error ERROR: _PIC_MODE_ON_ is used but not defined
#endif
#ifndef _CURRENT_MODE_
#error ERROR: _CURRENT_MODE_ is used but not defined
#endif
#if _CURRENT_MODE_==_PIC_MODE_ON_
    //for (int iPar=0;iPar<parSize; iPar++ ){
    if (iCase==0){
    for (int ii=0;ii<17;ii++){
      for(int jj=0;jj<17;jj++){
	for(int kk=0;kk<17;kk++){
	  double xLocation[3]={ii*0.25-2,jj*0.25-2,kk*0.25-2};
	  newNode=PIC::Mesh::mesh.findTreeNode(xLocation);
	  for (int iPar=0;iPar<14;iPar++){
	    if (newNode->Thread==PIC::ThisThread) {
	      PIC::Mesh::mesh.fingCellIndex(xLocation,i,j,k,newNode);
	    
	      newParticle=PIC::ParticleBuffer::GetNewParticle(newNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
	     
	      PIC::ParticleBuffer::SetV(v[iPar],newParticle);
	      PIC::ParticleBuffer::SetX(xLocation,newParticle);
	      PIC::ParticleBuffer::SetI(0,newParticle);
	      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticle);
	    }
	    
	  }
	  // }
	
	}
      }
    }
    }

    if (iCase==1){
      double protonNumDensity=4, antiprotonNumDensity=4;
      double Temperature=0.0;
      long int popNum1,popNum2;
      int LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      int GlobalParticleNumber;
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("Before cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
      std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;

      CleanParticles();
      LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("After cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);

      //PrepopulateDomain(int spec,double NumberDensity,double Temperature);
      popNum1=PrepopulateDomain(0,protonNumDensity,Temperature);
      popNum2=PrepopulateDomain(1,antiprotonNumDensity,Temperature);
 

      LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("After prepopulating, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
      std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;

    }
#endif

    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh.ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::Periodic::UpdateData();
      break;
    }

    for (int niter=0;niter<totalIter;niter++) {
    
      //PIC::Mesh::mesh.outputMeshDataTECPLOT("1.dat",0);
    
      //TransportEquation::TimeStep();
  
      PIC::TimeStep();
      //PIC::FieldSolver::Electromagnetic::ECSIM::TimeStep();

      //PIC::Mesh::mesh.outputMeshDataTECPLOT("2.dat",0);


      switch (_PIC_BC__PERIODIC_MODE_) {
      case _PIC_BC__PERIODIC_MODE_OFF_:
	PIC::Mesh::mesh.ParallelBlockDataExchange();
	break;

      case _PIC_BC__PERIODIC_MODE_ON_:
	PIC::BC::ExternalBoundary::Periodic::UpdateData();
	break;
      }


      char fname[100];
      if (_CURRENT_MODE_==_PIC_MODE_OFF_){
	sprintf(fname,"LightWaveCase%i.out=%i.dat",iCase,niter);
      }else{
	sprintf(fname,"PIC_particle_case%i.out=%i.dat",iCase,niter);
      }
    
      // if (niter%10==0) PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,0);
  
      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,0);
    }
  }

  MPI_Finalize();
  cout << "End of the run" << endl;

  return EXIT_SUCCESS;
}
