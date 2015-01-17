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
  double CharacteristicSpeed;

  switch (spec) {
  case _O2_SPEC_:
    CharacteristicSpeed=1.0e4;
    break;
  case _O2PLUS_SPEC_:
    CharacteristicSpeed=1.0e6;
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

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double v[3]={1.0e4,0.0,0.0};

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=-(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);

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

//  if (spec!=_O_SPEC_) return 0; //inject only spec=0

  static double v[3]={1.0e4,0.0,0.0};

  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=-(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);


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
  static double v[3]={1.0e4,0.0,0.0};

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]<0.0) {
        ModelParticlesInjectionRate+=-BlockSurfaceArea*(v[0]*ExternalNormal[0]+v[1]*ExternalNormal[1]+v[2]*ExternalNormal[2]);
      }
    }
  }

  return ModelParticlesInjectionRate;
}


bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  return false;
}

void TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

  memcpy(x_LOCAL,x,3*sizeof(double));
  memcpy(v_LOCAL,v,3*sizeof(double));

#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
    long int nd;
    char *offset;
    int i,j,k;
    PIC::Mesh::cDataCenterNode *CenterNode;
    double E[3],B[3];

    if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
      startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL,startNode);

      if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
    #endif //<-- _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_

    CenterNode=startNode->block->GetCenterNode(nd);
    offset=CenterNode->GetAssociatedDataBufferPointer();

    PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
    PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);

    double ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    double mass=PIC::MolecularData::GetMass(spec);

    accl_LOCAL[0]+=ElectricCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
    accl_LOCAL[1]+=ElectricCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
    accl_LOCAL[2]+=ElectricCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;
#endif //<-- _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

  //copy the local value of the acceleration to the global one
  memcpy(accl,accl_LOCAL,3*sizeof(double));
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

	//read the domain size from the data file
	PIC::CPLR::CCMC::LFM::GetDomainLimits(xmin,xmax,"/Users/ccmc/Example_Runs/LFM/IMFBy.LFM.Asher_Pembroke_112513_1_mhd_2000-01-01T07-07-00Z.cdf");

	for (idim=0;idim<3;idim++) {
	  double xCenter,dx;

	  xCenter=xmax[idim]+xmin[idim];
	  dx=0.8*(xmax[idim]-xmin[idim]);

	  xmin[idim]=xCenter-0.5*dx;
	  xmax[idim]=xCenter+0.5*dx;
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

	//read the data file
	PIC::CPLR::CCMC::LFM::LoadDataFile("/Users/ccmc/Example_Runs/LFM/IMFBy.LFM.Asher_Pembroke_112513_1_mhd_2000-01-01T07-07-00Z.cdf");

	//init the PIC solver
	PIC::Init_AfterParser();
	PIC::Mover::Init();

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=BoundingBoxInjectionRate;


	//set up the time step
	PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
	PIC::ParticleWeightTimeStep::initTimeStep();

  //init particle weight
	for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);

	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
	if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

	//init the particle buffer
	PIC::ParticleBuffer::Init(10000000);

  //init ICES
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
  #ifdef _ICES_CREATE_COORDINATE_LIST_
  PIC::CPLR::ICES::createCellCenterCoordinateList();
  PIC::CPLR::ICES::retriveSWMFdata(); //"Europa09"); //("RESTART_t001.52m"); //("MERCURY_RESTART_n070100");  ////("MERCURY_RESTART_n070001");
  #endif

  PIC::CPLR::ICES::readSWMFdata(1.0);
#endif

  PIC::Mesh::mesh.outputMeshDataTECPLOT("plasma-data.dat",0);
}





void amps_time_step () {
 		PIC::TimeStep();
}
