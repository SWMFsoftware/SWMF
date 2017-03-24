//$Id$

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
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "OH.h"




//the parameters of the domain
const double StretchCoefficient=1.0;
const double DebugRunMultiplier=4.0;

//the mesh resolution
double localResolution(double *x) {
  int idim;
  double lnR,res,r=0.0;
  
  res = OH::DomainDXMin;
  return res;
}

//set up the local time step
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed;

  if (PIC::ParticleWeightTimeStep::GlobalTimeStep) 
    if(PIC::ParticleWeightTimeStep::GlobalTimeStep[spec] > 0.0)
      return PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];

  switch (spec) {
  case _H_SPEC_: 
    if (PIC::nTotalSpecies == 1) CharacteristicSpeed=50.0E3;
    else CharacteristicSpeed=25.0E3;
    break;
  case _H_ENA_V1_SPEC_:
    CharacteristicSpeed=75.0E3;
    break;
  case _H_ENA_V2_SPEC_:
    CharacteristicSpeed=250.0E3;
    break;
  case _H_ENA_V3_SPEC_:
    CharacteristicSpeed=650.0E3;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the species is unknown");
  }

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}



double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;  
  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
   bool ExternalFaces[6];
   double ExternalNormal[3],ModelParticlesInjectionRate;
   int nface;

   if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
     for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
       startNode->GetExternalNormal(ExternalNormal,nface);
       ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(OH::InjectionNDensity,OH::InjectionTemperature,OH::InjectionVelocity,ExternalNormal,_H_SPEC_);

       if (ModelParticlesInjectionRate>0.0) return true;
     }
   }

   return false;
 }

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
   bool ExternalFaces[6];
   double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
   int nface,idim;
   long int newParticle;
   PIC::ParticleBuffer::byte *newParticleData;
   long int nInjectedParticles=0;

   if (spec!=_H_SPEC_) return 0; //inject only spec=0

   double v[3];


   double ModelParticlesInjectionRate;

   if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
     ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
     LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


     for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
       startNode->GetExternalNormal(ExternalNormal,nface);
       TimeCounter=0.0;

       //injection boundary conditions: -x -> inflow, other -> open boundary
       if (-ExternalNormal[0]<0.9) {
         nInjectedParticles+=PIC::BC::ExternalBoundary::OpenFlow::InjectBlock(spec,startNode,nface);
         continue; 
       } 


       ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(OH::InjectionNDensity,OH::InjectionTemperature,OH::InjectionVelocity,ExternalNormal,spec);


       if (ModelParticlesInjectionRate>0.0) {
	 ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

	 PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

	 while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
	   //generate the new particle position on the face
	   for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

	   //generate a particle
	   newParticle=PIC::ParticleBuffer::GetNewParticle();
	   newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
	   nInjectedParticles++;

	   //generate particles' velocity
	   PIC::Distribution::InjectMaxwellianDistribution(v,OH::InjectionVelocity,OH::InjectionTemperature,ExternalNormal,spec,-1);

	   PIC::ParticleBuffer::SetX(x,newParticleData);
	   PIC::ParticleBuffer::SetV(v,newParticleData);
	   PIC::ParticleBuffer::SetI(spec,newParticleData);
	   PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

	   // tagging the particle to the right population that it was created
	   double InjectionPressure = 0.0;
	   InjectionPressure = 2.0*OH::InjectionNDensity*Kbol*OH::InjectionTemperature;

	   OH::SetOriginTag(OH::GetEnaOrigin(OH::InjectionNDensity,InjectionPressure,OH::InjectionVelocity), newParticleData);

	   //inject the particle into the system
	   _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
	 }
       }


     }
   }

   return nInjectedParticles;
 }

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;
  
  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);
  
  return nInjectedParticles;
}

double BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  
  bool ExternalFaces[6];
  double ExternalNormal[3],BlockSurfaceArea;
  int nface;
  
  if (spec!=_H_SPEC_) return 0; //inject only spec=0
  
  double ModelParticlesInjectionRate=0.0;
  
  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
	startNode->GetExternalNormal(ExternalNormal,nface);
	BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);
	ModelParticlesInjectionRate+=BlockSurfaceArea*PIC::BC::CalculateInjectionRate_MaxwellianDistribution(OH::InjectionNDensity,OH::InjectionTemperature,OH::InjectionVelocity,ExternalNormal,spec);
      }
  }
  
  return ModelParticlesInjectionRate;
}


void amps_init_mesh(){


  PIC::InitMPI();
  
   //SetUp the alarm
  //  PIC::Alarm::SetAlarm(2000);
  
  rnd_seed();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  
  
  OH::Init_BeforeParser();
  PIC::Init_BeforeParser();
  OH::Init_AfterParser();
  
  
  
  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();
  
  //init the mesh
  if (PIC::Mesh::mesh.ThisThread==0){
    cout << "Init the mesh" << endl;}
  
  int maxBlockCellsnumber,minBlockCellsnumber,idim;
  
  maxBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);
  
  minBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);
  
  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(OH::DomainXMin,OH::DomainXMax,localResolution);
  PIC::Mesh::mesh.memoryAllocationReport();
  
  
  if (PIC::Mesh::mesh.ThisThread==0) {
    PIC::Mesh::mesh.buildMesh();
    PIC::Mesh::mesh.saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh.readMeshFile("mesh.msh");
  }
  
  // cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;
  
  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
  
  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();
  
#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif
  
  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();
  
  //initialize the blocks
  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();
  
  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();
  
#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif
  
  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();

  // allocate array with global times step and reset them
  if (! PIC::ParticleWeightTimeStep::GlobalTimeStep) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  }


}

void amps_init() {
 
  
  
  //init the PIC solver
   PIC::Init_AfterParser ();
   PIC::Mover::Init();

   
   //set up the time step
   PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
   PIC::ParticleWeightTimeStep::initTimeStep();

   //create the list of mesh nodes where the injection boundary conditinos are applied
   PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
   PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
   PIC::BC::InitBoundingBoxInjectionBlockList();
   
   //set up the particle weight
   PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=BoundingBoxInjectionRate;
   PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H_SPEC_);
   //setting the stat weight of the charge-exchange species such that wieght/dt is const.
   if (_H_ENA_V1_SPEC_>=0) {
     PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(_H_ENA_V1_SPEC_,PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H_SPEC_]*PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_ENA_V1_SPEC_]/PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_SPEC_]);
   }   
   if (_H_ENA_V2_SPEC_>=0) {
     PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(_H_ENA_V2_SPEC_,PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H_SPEC_]*PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_ENA_V2_SPEC_]/PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_SPEC_]);
   }   
   if (_H_ENA_V3_SPEC_>=0) {
     PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(_H_ENA_V3_SPEC_,PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H_SPEC_]*PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_ENA_V3_SPEC_]/PIC::ParticleWeightTimeStep::GlobalTimeStep[_H_SPEC_]);
   }   

   MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
   if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;
   
   //init the particle buffer
   PIC::ParticleBuffer::Init(10000000);

   // change global time step if it's set in the input file
   if(OH::UserGlobalTimeStep > 0.0){
     if (PIC::Mesh::mesh.ThisThread==0) 
       cout << "AMPS:: Global time steps of all species are changed to " << 
	 OH::UserGlobalTimeStep << 
	 " by user (see oh.input)\n";
     for(int spec=0; spec < PIC::nTotalSpecies; spec++)
       PIC::ParticleWeightTimeStep::GlobalTimeStep[spec] = OH::UserGlobalTimeStep;
   }
}

//time step
void amps_time_step(){

    //make the time advance
    static int LastDataOutputFileNumber=0;

    // change time after las coupling session for each species
    for(int spec=0; spec < PIC::nTotalSpecies; spec++)
      OH::Coupling::TimeAfterCoupling[spec] += 
	PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];

    //make the time advance
     PIC::TimeStep();

     // write output file
     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) 
	 cout << "AMPS: Output file " << PIC::DataOutputFileNumber<< " is done" << endl;
     }
     
}



