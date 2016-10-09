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
#include "Earth.h"


const double rSphere=_RADIUS_(_TARGET_);

double xMaxDomain=5; //modeling the vicinity of the planet
double yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


double dxMinSphere=0.1,dxMaxSphere=0.1;
double dxMinGlobal=1,dxMaxGlobal=1;


//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
  double res=0.0;
  return res;
}


//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double res,r,l[3] = {1.0,0.0,0.0};


  if ( (strcmp(Earth::Mesh::sign,"0x301020156361a50")!=0)) {
    //test mesh
    return 0.05*_RADIUS_(_TARGET_);
  }
  else {
    return 0.5*_RADIUS_(_TARGET_);
  }





  res=dxMinSphere;
  res/=2.1;
  return rSphere*res;
}


double localResolution(double *x) {
  double res;

  if ( (strcmp(Earth::Mesh::sign,"0x301020156361a50")!=0)) {
    //test mesh
    return 0.5*_RADIUS_(_TARGET_);
  }
  else {
    return 0.5*_RADIUS_(_TARGET_);
  }

  //  if (strcmp(Earth::Mesh::sign,"new")==0) 
{ // new mesh
    int idim;
    double r=0.0;

    for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

    r=sqrt(r);

    if (r<0.98*rSphere) return rSphere;
    else if (r<1.05*rSphere) return localSphericalSurfaceResolution(x);
    else return rSphere * dxMinGlobal;
  }
}

//set up the local time step
double localTimeStep(int spec,
		     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  CellSize = startNode->GetCharacteristicCellSize();
  
  const double maxInjectionEnergy=1E5*1000*ElectronCharge;
  double CharacteristicSpeed = sqrt(2*maxInjectionEnergy/_MASS_(_O_));
  return 0.2* 0.3*CellSize/CharacteristicSpeed;
}


double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  return 0.0;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  
  return false;
}



//injection of model particles through the faces of the bounding box
double  BoundingBoxInjection(int spec,
			       cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  return 0;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  return 0;
}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;
  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,
			      double *x,double *v, 
			      double &dtTotal, void *NodeDataPonter,
			      void *SphereDataPointer)  {
  return _PARTICLE_DELETED_ON_THE_FACE_;
}




double sphereInjectionRate(int spec,void *SphereDataPointer) {
  double res=0.0;
  return res;
}



void amps_init_mesh() {

  //if (strcmp(Earth::Mesh::sign,"new")==0) 
{ //full mesh

    dxMinGlobal=0.4/2.1,dxMaxGlobal=1;

    //modeling the vicinity of the planet
    xMaxDomain=16 *_RADIUS_(_TARGET_); 
    //the minimum size of the domain in the direction perpendicular
    //to the direction to the sun
    yMaxDomain=16 * _RADIUS_(_TARGET_); 

    dxMinSphere=20E3,dxMaxSphere=100E3;
  }
//  else exit(__LINE__,__FILE__,"Error: unknown option");


 PIC::InitMPI();
 rnd_seed();

 MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


 Earth::Init_BeforeParser();

 //init the particle solver
 PIC::Init_BeforeParser();



 const int InitialSampleLength=600;
 
 //register the sphere
 static const bool SphereInsideDomain=true;
 
 if (SphereInsideDomain) {
   double sx0[3]={0.0,0.0,0.0};
   cInternalBoundaryConditionsDescriptor SphereDescriptor;
   cInternalSphericalData *Sphere;
   
   //reserve memory for sampling of the surface balance of sticking species
   long int ReserveSamplingSpace[PIC::nTotalSpecies];
   for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=0;
   
   
   cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);
   
   
   
   PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
   SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
   Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
   Sphere->SetSphereGeometricalParameters(sx0,rSphere);
   
   
   
   
   
   
   
   
   Sphere->Radius=_RADIUS_(_TARGET_);
   Sphere->PrintSurfaceMesh("Sphere.dat");
   Sphere->PrintSurfaceData("SpheraData.dat",0);
   Sphere->localResolution=localSphericalSurfaceResolution;
   Sphere->faceat=0;
   
   Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);
 }
 
 //init the solver
 PIC::Mesh::initCellSamplingDataBuffer();
 
 //init the mesh
 if (PIC::ThisThread==0)
   cout << "Init the mesh" << endl;
 
 int maxBlockCellsnumber,minBlockCellsnumber,idim;
 
 maxBlockCellsnumber=_BLOCK_CELLS_X_;
 if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
 if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);
 
 minBlockCellsnumber=_BLOCK_CELLS_X_;
 if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
 if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);
 
 double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

 for (idim=0;idim<DIM;idim++) {
   xmax[idim]= 16 * _RADIUS_(_TARGET_);
   xmin[idim]=-16 * _RADIUS_(_TARGET_);
 }
 
 
 //generate only the tree
 PIC::Mesh::mesh.AllowBlockAllocation=false;
 PIC::Mesh::mesh.init(xmin,xmax,localResolution);
 PIC::Mesh::mesh.memoryAllocationReport();
 
 
 char mesh[200]="";
 bool NewMeshGeneratedFlag=false;
 FILE *fmesh=NULL;
 
 sprintf(mesh,"amr.sig=%s.mesh.bin",Earth::Mesh::sign);
 fmesh=fopen(mesh,"r");
 
 if (fmesh!=NULL) {
   fclose(fmesh);
   PIC::Mesh::mesh.readMeshFile(mesh);
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
 
 if (NewMeshGeneratedFlag==true) PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
 
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
 
 
 if (PIC::ThisThread==0) cout << "AMPS' Initialization is complete" << endl;
 
 MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
 
 
}
 
 void amps_init() {
   int idim;
   
   //init the PIC solver
   PIC::Init_AfterParser ();
   PIC::Mover::Init();
   
   //set up the time step
   PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
   PIC::ParticleWeightTimeStep::initTimeStep();
   
   //set up the particle weight
   PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=BoundingBoxInjection;
   PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1.0);
   
   MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
   if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;
   
   //output final data
   //create the list of mesh nodes where the injection boundary conditions are applied
   PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
   PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
   PIC::BC::InitBoundingBoxInjectionBlockList();
   
   
   //init the particle buffer
   PIC::ParticleBuffer::Init(10000000);

   int LastDataOutputFileNumber=-1;
   
   
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_

#if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_
   // BATL reader

   // initialize the reader
   PIC::CPLR::DATAFILE::BATSRUS::Init("3d__ful_2_t00000000_n00020000.idl");
   PIC::CPLR::DATAFILE::BATSRUS::LoadDataFile();


#else
   exit(__LINE__,__FILE__,"ERROR: unrecognized datafile reader mode");
   
#endif //_PIC_COUPLER_DATAFILE_READER_MODE_
   

#endif //_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
	

    //init particle weight of neutral species that primary source is sputtering



  if (_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE_==_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE__TECPLOT_ASCII_) {
    PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);
  }



}

 //time step

void amps_time_step () {
  
  PIC::TimeStep();
  
}
