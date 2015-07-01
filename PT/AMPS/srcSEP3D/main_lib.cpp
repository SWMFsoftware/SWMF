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
static double DomainXMin[3]={8.76E8,-0.5*DomainDX,-1.3E7};
static double DomainXMax[3]={9.40E8, 0.5*DomainDX, 1.3E7};

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
/*
bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double v[3]={1.0e6,0.0,0.0}, NDensity=1.0E+13;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=-NDensity*(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);

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

  static double v[3]={1.0e6,0.0,0.0}, NDensity=1.0E+13;

  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=-NDensity*(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);


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

          //apply condition of tracking the particle
          #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
          PIC::ParticleTracker::InitParticleID(newParticleData);
          PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData);
          #endif

          //generate particles' velocity
          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

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


  double ModelParticlesInjectionRate=0.0;
  static double v[3]={1.0e6,0.0,0.0}, NDensity=1.0E+13;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]<0.0) {
        ModelParticlesInjectionRate+=-NDensity*BlockSurfaceArea*(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);
      }
    }
  }

  return ModelParticlesInjectionRate;
}
*/

bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  return false;
}

void amps_init_mesh() {
  PIC::InitMPI();
	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

	//init the particle solver
	PIC::Init_BeforeParser();

	PIC::Mover::Init_BeforeParser();

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
	PIC::CPLR::DATAFILE::MULTIFILE::Init("dataARMS",20);
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

  {  // prepopulate the domain
    double NDensity=1.0E+10, Temperature=6000, Velocity[3]={1.0E6,0.0,0.0};
    for (int s=0;s<PIC::nTotalSpecies;s++)
      PIC::InitialCondition::PrepopulateDomain(s,NDensity, Velocity, Temperature);
  }
    
  PIC::Mesh::mesh.outputMeshDataTECPLOT("plasma-data.dat",0);

}





void amps_time_step () {

  PIC::TimeStep();

}
