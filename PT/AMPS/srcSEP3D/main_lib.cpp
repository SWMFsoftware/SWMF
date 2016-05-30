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

#include "pic.h"
#include "constants.h"

#include "Exosphere.h"

#include "SEP3D.h"

//default definition of the functions for the exospehre module ==================================
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {return 0.0;}
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {return 0.0;}
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EPHIOD";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Europa";
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;}
//===============================================================================================

static double DomainDX = 2E6;
static double DomainXMin[3]={8.760E8,-0.5*DomainDX,-1.3E7};
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
static double DomainXMax[3]={9.400E8, 0.5*DomainDX, 1.3E7};
#else
static double DomainXMax[3]={9.445E8, 0.5*DomainDX, 1.3E7};
#endif

double localResolution(double *x) {

  return DomainDX;
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CharacteristicSpeed;

  switch (spec) {
  case _H_PLUS_SPEC_: case _ELECTRON_SPEC_:
    CharacteristicSpeed=1.0e7;
    break;
  default:
    exit(__LINE__,__FILE__,"unknown species");
   }

  return 0.3*startNode->GetCharacteristicCellSize()/CharacteristicSpeed;
}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
	double res=1.0;

	// for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

	return res;
}

bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  return false;
}



void amps_init_mesh() {
  PIC::InitMPI();
	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

	//init the particle solver
	PIC::Init_BeforeParser();

	PIC::Mover::Init_BeforeParser();

#if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
	// initialize field lines
	PIC::FieldLine::Init();
#endif//_PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_

	//init the solver
	PIC::Mesh::initCellSamplingDataBuffer();

	//init the mesh
	double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};
	int idim;


	//generate only the tree
	PIC::Mesh::mesh.AllowBlockAllocation=false;
	PIC::Mesh::mesh.init(DomainXMin,DomainXMax,localResolution);
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

	//read the data file
	if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_)
	  PIC::CPLR::DATAFILE::MULTIFILE::Init(true,33);
	else
	  PIC::CPLR::DATAFILE::MULTIFILE::Init(true,300);//15
}

void amps_init(){
	//init the PIC solver
	PIC::Init_AfterParser();
	PIC::Mover::Init();

  //create the list of mesh nodes where the injection boundary conditinos are applied
	/*
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=BoundingBoxInjectionRate;
	*/

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();
  
  //init particle weight
  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,3.0E+26);
  

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;
  
  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
#if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
  {  // create field lines and inject particles
    //    PIC::FieldLine::Init();
    double xStart[3] = {0.0,0.0,-1.3e+6};
    for(xStart[0] = 9.435E+8; xStart[0]> 9.40E+8; xStart[0]-=0.0025E+8)
      PIC::FieldLine::InitLoop2D(xStart,0.1,1e+5,3e+5);
    for(int i=0; i<2000000; i++)
      PIC::FieldLine::InjectParticle(0);
    SEP3D::GlobalEnergyDistribution::sample();
    SEP3D::GlobalEnergyDistribution::print("EnergyDistribution.000000",0);
    }
#else
  {  // prepopulate the domain
    double NDensity=1.0E+10, Temperature=6000, Velocity[3]={1.0E6,0.0,0.0};
    for (int s=0;s<PIC::nTotalSpecies;s++)
      PIC::InitialCondition::PrepopulateDomain(s,NDensity, Velocity, Temperature);
  }
#endif
    
  PIC::Mesh::mesh.outputMeshDataTECPLOT("plasma-data.dat",0);

}



int amps_time_step () {

#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
#else
  static int cnt=0;
  static int file=0;
  cnt++;
  if(cnt==10){
    file++;
    cnt = 0;
    char Name[100];
    sprintf(Name,"EnergyDistribution.%06d", file);
    SEP3D::GlobalEnergyDistribution::sample();
    SEP3D::GlobalEnergyDistribution::print(Name,0);
  }
#endif
  return PIC::TimeStep();
  
}
