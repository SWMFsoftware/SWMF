
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

//$Id$

#include "pic.h"
#include "constants.h"

#include "Exosphere.h"

//default definition of the functions for the exospehre module
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {return 0.0;}
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {return 0.0;}
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EPHIOD";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Europa";
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;}

const double DomainLength[3]={1.0E7,1.0E7,1.0E7},DomainCenterOffset[3]={0.0,0.0,0.0};

double localResolution(double *x) {

	return 0.025*DomainLength[0];
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	double CellSize;

//	double CharacteristicSpeed_NA=5.0E3;

	const double maxInjectionEnergy=1E5*1000*ElectronCharge;
	double CharacteristicSpeed = sqrt(2*maxInjectionEnergy/_MASS_(_O_));


	//  CharacteristicSpeed*=sqrt(PIC::MolecularData::GetMass(NA)/PIC::MolecularData::GetMass(spec));

	CellSize=startNode->GetCharacteristicCellSize();




//	CharacteristicSpeed=1.0E6;

/*if (spec==_OPLUS_HIGH_SPEC_) CharacteristicSpeed=10.0*1.6E6;
if (spec==_OPLUS_THERMAL_SPEC_) CharacteristicSpeed=10.0*9.6E4;*/

  switch (spec) {
  case _OPLUS_HIGH_SPEC_:
    CharacteristicSpeed=10.0*1.6E6;
    break;
  case _OPLUS_THERMAL_SPEC_:
    CharacteristicSpeed=10.0*9.6E4;
    break;

  case _O2_SPEC_:
    CharacteristicSpeed=1.0e4;
    break;
  case _O2PLUS_SPEC_:
    CharacteristicSpeed=10*1.0e4;
    break;
  default:
    exit(__LINE__,__FILE__,"unknown species");
   }

	return 0.3*CellSize/CharacteristicSpeed;


}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
	double res=1.0;

	// for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

	return res;
}



void amps_init() {
  PIC::InitMPI();
	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

	//init the particle solver
	PIC::Init_BeforeParser();

	//init the solver
	PIC::Mesh::initCellSamplingDataBuffer();

	//init the mesh
	double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};
	int idim;

	for (idim=0;idim<DIM;idim++) {
		xmax[idim]=DomainCenterOffset[idim]+DomainLength[idim]/2.0;
		xmin[idim]=DomainCenterOffset[idim]-DomainLength[idim]/2.0;
	}

	//generate only the tree
	PIC::Mesh::mesh.AllowBlockAllocation=false;
	PIC::Mesh::mesh.init(xmin,xmax,localResolution);
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


	//init the PIC solver
	PIC::Init_AfterParser ();
	PIC::Mover::Init();

	//set up the time step
	PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
	PIC::ParticleWeightTimeStep::initTimeStep();

      	//init particle weight
	for (int s=0;s<PIC::nTotalSpecies;s++) PIC::Mesh::mesh.rootBlock->SetLocalParticleWeight(0.0,s);

	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
	if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

	//init the particle buffer
	PIC::ParticleBuffer::Init(10000000);
}

	//time step

void amps_time_step () {
 		PIC::TimeStep();
}
