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
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"

double BulletLocalResolution(double *x) {
  int idim;
  double res,r;



return  ((fabs(x[0])<100.0)||(x[1]*x[1]+x[2]*x[2]<40.0*40.0)) ? 5.0 : 100.0;

  r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

  if (r<0.5) res=10.0;
  else if (res<5.0) res=0.24;
  else res=10.0;

  return res;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double dt,c,vInitUnchanged[3],LocalParticleWeight,RateFactor,mass;
  int spec;

  memcpy(vInitUnchanged,vInit,3*sizeof(double));

  c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];
  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  //sample the energy and momentum transfer rates
  spec=PIC::ParticleBuffer::GetI(ptr);
  mass=PIC::MolecularData::GetMass(spec);

  LocalParticleWeight=startNode->block->GetLocalParticleWeight(spec);
  LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

  dt=startNode->block->GetLocalTimeStep(spec);
  RateFactor=mass*LocalParticleWeight/dt; 

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  double *xx,dp,de; 

  dp=(vInit[0]-vInitUnchanged[0])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpX[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateX;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  dp=(vInit[1]-vInitUnchanged[1])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpY[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateY;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  dp=(vInit[2]-vInitUnchanged[2])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpZ[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateZ;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  de=(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2] -
      vInitUnchanged[0]*vInitUnchanged[0]-vInitUnchanged[1]*vInitUnchanged[1]-vInitUnchanged[2]*vInitUnchanged[2])*RateFactor;
  xx=TriangleCutFace->UserData.EnergyTransferRate;
  xx[spec]+=de/TriangleCutFace->SurfaceArea;

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  double *xx,dp,de; 
  int thread=omp_get_thread_num();

  //dpX
  dp=(vInit[0]-vInitUnchanged[0])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpX[thread]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateX;
  dp/=TriangleCutFace->SurfaceArea;

  #pragma omp atomic
  xx[spec]+=dp;

  //dpY
  dp=(vInit[1]-vInitUnchanged[1])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpY[thread]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateY;
  dp/=TriangleCutFace->SurfaceArea; 

  #pragma omp atomic
  xx[spec]+=dp;

  //dpX
  dp=(vInit[2]-vInitUnchanged[2])*RateFactor;
  Orbiter::Sampling::DragCoefficient::dpZ[thread]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateZ;
  dp/=TriangleCutFace->SurfaceArea;

  #pragma omp atomic
  xx[spec]+=dp;

  //dE
  de=(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2] -
      vInitUnchanged[0]*vInitUnchanged[0]-vInitUnchanged[1]*vInitUnchanged[1]-vInitUnchanged[2]*vInitUnchanged[2])*RateFactor;
  xx=TriangleCutFace->UserData.EnergyTransferRate; 
  de/=TriangleCutFace->SurfaceArea; 

  #pragma omp atomic
  xx[spec]+=de;

#else
  exit(__LINE__,__FLIE__,"Error: the option is unknown");
#endif



  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  double res,size;

  size=t->CharacteristicSize();

  if (size<0.01) res=0.01;
  else if (size<1.0) res=0.01*pow(10.0,size);
  else res=0.25;

  //reduce the mesh resolution when run tests
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  if (res<0.20) res=0.20;
  #endif

  return res;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed=50.0E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}



double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}




/*  double Orbiter::GetTotalProduction(int spec,void *BoundaryElement) {
    return 1.0E20;
  }

  double Orbiter::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
    return GetTotalProduction(spec,BoundaryElement);
  }*/



int main(int argc,char **argv) {
  char fname[_MAX_STRING_LENGTH_PIC_];

  //init the particle solver
  PIC::InitMPI();

  //seed the random number generatore
  if (_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_) {
    rnd_seed(100);
  }


  PIC::Init_BeforeParser();
  PIC::Alarm::SetAlarm(48*3600-30*60);

  Orbiter::Init_BeforeParser();


  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  sprintf(fname,"%s/%s",PIC::UserModelInputDataPath,Orbiter::SurfaceModel::MeshFileName);

  switch (Orbiter::SurfaceModel::MeshFileFormat) {
  case Orbiter::SurfaceModel::MeshFileFormat_CEA:
    PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(fname,20.0);
    break;

  case Orbiter::SurfaceModel::MeshFileFormat_NASTRAN:
    PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat(fname);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the format is unknown");
  }


  PIC::Mesh::IrregularSurface::GetSurfaceSizeLimits(xmin,xmax);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);

  //calculate the projection are when the drag coefficient it to be calculated 
  if (Orbiter::Sampling::DragCoefficient::SamplingMode==true) {
    Orbiter::CalculateProjectionArea();
  }


//  PIC::Mesh::IrregularSurface::SmoothRefine(0.5);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }


  //set up the s/c object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor OrbiterSurfaceDescriptor;
  cInternalNastranSurfaceData *Orbiter;

//  PIC::BC::InternalBoundary::RotationBody::Init(ReserveSamplingSpace,NULL);
  OrbiterSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  Orbiter=(cInternalNastranSurfaceData*) OrbiterSurfaceDescriptor.BoundaryElement;

//  Orbiter->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
//  Orbiter->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;

//  Nucleus->localResolution=Comet::localSphericalSurfaceResolution;
//  Orbiter->InjectionRate=Orbiter::GetTotalProduction;
  Orbiter->faceat=0;
//  Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
//  Orbiter->InjectionBoundaryCondition=Orbiter::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
//  Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
  //PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;


  for (int idim=0;idim<3;idim++) {
    double OrbiterSize=xmax[idim]-xmin[idim];

    xmin[idim]-=OrbiterSize*Orbiter::DomainSize::xMinOffset[idim];
    xmax[idim]+=OrbiterSize*Orbiter::DomainSize::xMaxOffset[idim];
  }

  PIC::Mesh::mesh.CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();


  //PIC::Mesh::mesh.buildMesh();

  //generate mesh or read from file
  char mesh[400];
  bool NewMeshGeneratedFlag=false;

  FILE *fmesh=NULL;

  sprintf(mesh,"%s/amr.sig=%s.mesh.bin",PIC::UserModelInputDataPath,Orbiter::Mesh::sign);
  fmesh=fopen(mesh,"r");

  if (fmesh!=NULL) {
    if (PIC::ThisThread==0) printf("$PREFIX: mesh %s is found: loading.... ",mesh);

    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(mesh);

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) printf("done. \n");
  }
  else {
    NewMeshGeneratedFlag=true;
    if (PIC::ThisThread==0) {
      printf("$PREFIX: mesh %s is not found: generating.... ",mesh);
      fflush(stdout);
    }

    if (PIC::Mesh::mesh.ThisThread==0) {
      PIC::Mesh::mesh.buildMesh();
      PIC::Mesh::mesh.saveMeshFile("mesh.msh");
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      printf("  done.\n",mesh);
      fflush(stdout);
    }
  }


  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();



  //initialize the blocks
  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();


  //init the external normals of the cut faces
  double xLightSource[3]={6000.0e3,0,0.0};

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  PIC::Mesh::IrregularSurface::InitExternalNormalVector();

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  //test the shadow procedure
  PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,true);

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation-shadow.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh(fname);
  }

  //output the volume mesh
  sprintf(fname,"%s/VolumeMesh.dat",PIC::OutputDataFileDirectory);
  PIC::Mesh::mesh.outputMeshTECPLOT(fname);

  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      printf("$PREFIX: Renaming the mesh file to amr.sig=0x%lx.mesh.bin.....  ",MeshSignature);
      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
      printf(" done.\n");
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);



//  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=2000; //0; //00; //*10;
//  PIC::RequiredSampleLength=500; //00; //0; //0;


  PIC::Init_AfterParser ();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();



  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=Orbiter::UpstreamBC::BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=Orbiter::UpstreamBC::BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();


  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=Orbiter::UpstreamBC::BoundingBoxInjectionRate;
  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);

  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::AdjustParticleWeight_ConstantWeightOverTimeStep_KeepMinParticleWeight(s);

  //init the particle buffer
//  PIC::ParticleBuffer::Init(1000000);

  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;


  //determine the total number of the iterations to perform
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=550;

  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
    static int LastDataOutputFileNumber=-1;

    PIC::TimeStep();

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>1000) PIC::RequiredSampleLength=1000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }

    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }
  }

  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_Orbiter.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }

  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
