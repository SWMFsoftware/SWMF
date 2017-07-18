//$Id$


#include "pic.h"
#include "constants.h"
#include "Surface.h"

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



namespace ParticleSurfaceInterationModel {
  int COPS(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    double BulkFlowVelocity[3]={0.0,0.0,0.0},Twall=293.0;
    int spec;

    spec=PIC::ParticleBuffer::GetI(ptr); 

    //reflect particles from the walls with themerature of 293 K;
    do {
      PIC::Distribution::MaxwellianVelocityDistribution(vInit,BulkFlowVelocity,Twall,spec);
    }
    while (Vector3D::DotProduct(vInit,TriangleCutFace->ExternalNormal)<=0.0);

    return _PARTICLE_REJECTED_ON_THE_FACE_; 
  }

  int ParticleSurfaceInteraction_CYCNSS(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    double BulkFlowVelocity[3]={0.0,0.0,0.0},Twall=50.0+273.0;
    int spec;

    spec=PIC::ParticleBuffer::GetI(ptr);

    if (TriangleCutFace->attribute==8) Twall=100.0+273.0;

    //reflect particles from the walls with themerature of 293 K;
    do {
      PIC::Distribution::MaxwellianVelocityDistribution(vInit,BulkFlowVelocity,Twall,spec);
    }
    while (Vector3D::DotProduct(vInit,TriangleCutFace->ExternalNormal)<=0.0);

    return _PARTICLE_REJECTED_ON_THE_FACE_;
  }
}


//namespace containing function used in the nightly test 
namespace DragCoefficientTest {
  int iTestCode=0;

  const int TestCode__Adsorption=0;
  const int TestCode__DiffuseReflection=1;
  const int TestCode__MaxwellReflection=2; 
  const int TestCode__QuasiSpecularReflection=3;
  const int TestCode__SpecularReflection=4;
  const int TestCode__CLL=5; 

  const int nTotalTestCodes=6;

  int ParticleSurfaceInteractionProcessor_TEST(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
    int res=_PARTICLE_REJECTED_ON_THE_FACE_;

    switch (iTestCode) {
    case TestCode__Adsorption:
      res=_PARTICLE_DELETED_ON_THE_FACE_;
      break;

    case TestCode__DiffuseReflection:
      res=Surface::DiffuseReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode); 
      break; 

    case TestCode__MaxwellReflection:
      res=Surface::MaxwellReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
      break;

    case TestCode__QuasiSpecularReflection: 
   //   res=Surface::QuasiSpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);

      res=_PARTICLE_DELETED_ON_THE_FACE_;
      break;

    case TestCode__SpecularReflection:
      res=Surface::SpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
      break;

    case TestCode__CLL:
      //res=Surface::CLL::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);

      res=_PARTICLE_DELETED_ON_THE_FACE_;
      break;

    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown"); 
    } 

    return res;
  }

  void RemoveAllParticles() {
    int Ptr,nextPrt;

    //reset particle lists in all cells
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    PIC::Mesh::cDataBlockAMR *block;
    int i,j,k;

    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];

      if ((block=node->block)!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) if ((Ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)])!=-1) {
        block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1; 

        do {
          nextPrt=PIC::ParticleBuffer::GetNext(Ptr);
          PIC::ParticleBuffer::DeleteParticle(Ptr);
  
          Ptr=nextPrt;
        }
        while (Ptr!=-1); 

      }
    }
  }

  void ExecuteTest() {
    int LastDataOutputFileNumber=-1,nTotalIterations=600;   
    int initRequiredSampleLength=PIC::RequiredSampleLength; 

    RemoveAllParticles(); 
    PIC::DataOutputFileNumber=0;

    //time step
    for (long int niter=0;;niter++) {
      static int LastDataOutputFileNumber=-1;

      PIC::TimeStep();

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        if (niter>nTotalIterations) break;

        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>1000) PIC::RequiredSampleLength=1000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }

      if (PIC::Mesh::mesh.ThisThread==0) {
        time_t TimeValue=time(NULL);
        tm *ct=localtime(&TimeValue);

        printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current iTest=%i, current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,iTestCode,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }
  }

    //output the particle statistics for the nightly tests
    char fname[400];

    switch (iTestCode) {
    case TestCode__Adsorption:
      sprintf(fname,"%s/test_Orbiter--Adsorption.dat",PIC::OutputDataFileDirectory);
      
      if (PIC::ThisThread==0) {
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--Adsorption.dat"); 
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--Adsorption.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--Adsorption.dat"); 
      }

      break;
    case TestCode__DiffuseReflection:
      sprintf(fname,"%s/test_Orbiter--DiffuseReflection.dat",PIC::OutputDataFileDirectory);

      if (PIC::ThisThread==0) {
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--DiffuseReflection.dat");
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--DiffuseReflection.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--DiffuseReflection.dat");
      }

      break;
    case TestCode__MaxwellReflection:
      sprintf(fname,"%s/test_Orbiter--MaxwellReflection.dat",PIC::OutputDataFileDirectory);
   
      if (PIC::ThisThread==0) {
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--MaxwellReflection.dat");
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--MaxwellReflection.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--MaxwellReflection.dat");
      }

      break;
    case TestCode__QuasiSpecularReflection:
      sprintf(fname,"%s/test_Orbiter--QuasiSpecularReflection.dat",PIC::OutputDataFileDirectory);

      if (PIC::ThisThread==0) { 
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--QuasiSpecularReflection.dat");
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--QuasiSpecularReflection.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--QuasiSpecularReflection.dat");
      }

      break;
    case TestCode__SpecularReflection:
      sprintf(fname,"%s/test_Orbiter--SpecularReflection.dat",PIC::OutputDataFileDirectory);

      if (PIC::ThisThread==0) {
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--SpecularReflection.dat");
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--SpecularReflection.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--SpecularReflection.dat");
      }

      break;
    case TestCode__CLL:
      sprintf(fname,"%s/test_Orbiter--CLL.dat",PIC::OutputDataFileDirectory);
 
      if (PIC::ThisThread==0) {
        system("mv PT/plots/amps.cut-cell.surface-data.out=2.dat PT/plots/amps.cut-cell.surface-data.out=2--CLL.dat");
        system("mv PT/plots/pic.H2O.s=0.out=2.dat PT/plots/pic.H2O.s=0.out=2--CLL.dat");
        system("mv PT/plots/test_Orbiter_Drag_Coefficient.dat PT/plots/test_Orbiter_Drag_Coefficient--CLL.dat");
      }

      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }

    //output the run statistics
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

    //reset the required sample length
    PIC::RequiredSampleLength=initRequiredSampleLength;
  }


  void ExecuteAllTests() {

    //set the test conditions 
    Surface::Temeprature::Isothremal::Temp=300.0;
    Surface::MaxwellReflection::AccommodationCoefficient[0]=0.5;

    if (PIC::nTotalSpecies!=1) exit(__LINE__,__FILE__,"Error: the test is set to be exdcuted only with a single model speceis");


    //execute the tests
    for (iTestCode=0;iTestCode<nTotalTestCodes;iTestCode++) ExecuteTest();  

    MPI_Finalize();
    cout << "End of the run:" << PIC::nTotalSpecies << endl;

    exit(EXIT_SUCCESS); 
  }
}  



int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double dt,c,vInitUnchanged[3],ParticleWeight,RateFactor,mass;
  int code,spec;

  spec=PIC::ParticleBuffer::GetI(ptr);
  mass=PIC::MolecularData::GetMass(spec);
  ParticleWeight=startNode->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr); 

  memcpy(vInitUnchanged,vInit,3*sizeof(double));

  //process particle/surface interaction event
  code=_ORBITER__PARTICLE_SURFACE_INTERACTION_PROCESSOR_(ptr,xInit,vInit,TriangleCutFace,startNode);

  //check for adsorption of the particle on teh face 
  if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
    if (_ORBITER_ADSORPTION_MODEL__MODE_ == _ORBITER_ADSORPTION_MODEL__MODE__OFF_) {
      exit(__LINE__,__FILE__,"Error: with the surrent settings only refrection on the surface is permitted. Change the model settings");
    }

    //the particle has been adsorbed on the surface
    vInit[0]=0.0,vInit[1]=0.0,vInit[2]=0.0;  //this is needed to calculating of the momentum transfer to the surface

    //update the adsorption flux
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    double *p=TriangleCutFace->UserData.AdsorptionFlux;

    #pragma omp atomic
    p[spec]+=ParticleWeight;

    #else
    TriangleCutFace->UserData.AdsorptionFlux[spec]+=ParticleWeight;
    #endif
  }


  //sample the energy and momentum transfer rates
  dt=startNode->block->GetLocalTimeStep(spec);
  RateFactor=mass*ParticleWeight/dt; 

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  double *xx,dp,de; 

  Orbiter::Sampling::DragCoefficient::Flux[0]+=ParticleWeight/dt;
  Orbiter::Sampling::DragCoefficient::ModelParticleCounter[0]++;

  dp=-(vInit[0]-vInitUnchanged[0])*RateFactor;  //the rate of the momentu gain by the orbiter's surface
  Orbiter::Sampling::DragCoefficient::dpX[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateX;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  dp=-(vInit[1]-vInitUnchanged[1])*RateFactor;   //the rate of the momentu gain by the orbiter's surface
  Orbiter::Sampling::DragCoefficient::dpY[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateY;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  dp=-(vInit[2]-vInitUnchanged[2])*RateFactor;   //the rate of the momentu gain by the orbiter's surface
  Orbiter::Sampling::DragCoefficient::dpZ[0]+=dp;
  xx=TriangleCutFace->UserData.MomentumTransferRateZ;
  xx[spec]+=dp/TriangleCutFace->SurfaceArea;

  de=(vInit[0]*vInit[0]+vInit[1]*vInit[1]+vInit[2]*vInit[2] -
      vInitUnchanged[0]*vInitUnchanged[0]-vInitUnchanged[1]*vInitUnchanged[1]-vInitUnchanged[2]*vInitUnchanged[2])*RateFactor/2.0;
  xx=TriangleCutFace->UserData.EnergyTransferRate;
  xx[spec]+=de/TriangleCutFace->SurfaceArea;

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  double *xx,dp,de; 
  int thread=omp_get_thread_num();

  Orbiter::Sampling::DragCoefficient::Flux[thread]+=ParticleWeight/dt;
  Orbiter::Sampling::DragCoefficient::ModelParticleCounter[thread]++;

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
      vInitUnchanged[0]*vInitUnchanged[0]-vInitUnchanged[1]*vInitUnchanged[1]-vInitUnchanged[2]*vInitUnchanged[2])*RateFactor/1.0;
  xx=TriangleCutFace->UserData.EnergyTransferRate; 
  de/=TriangleCutFace->SurfaceArea; 

  #pragma omp atomic
  xx[spec]+=de;

#else
  exit(__LINE__,__FLIE__,"Error: the option is unknown");
#endif



  return code;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  double res,size;

  size=t->CharacteristicSize();

  if (size<0.01) res=0.01;
  else if (size<1.0) res=0.01*pow(10.0,size);
  else res=0.25;

  //reduce the mesh resolution when run tests
  if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) { 
    if (_ORBITER__NIGHTLY_TEST_REDUCE_RESOLUTION_MODE_==_PIC_MODE_ON_) {  
      res*=_PIC_NIGHTLY_TEST__GRID_RESOLUTION_MULTIPLIER_;
    }

    //limit the grid resolution with the floor value 
    if (res<_PIC_NIGHTLY_TEST__GRID_RESOLUTION_FLOOR_VALUE_) res=_PIC_NIGHTLY_TEST__GRID_RESOLUTION_FLOOR_VALUE_; 
  }

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
  list<CutCell::cSurfaceMeshFile> SurfaceMeshFileList;
  CutCell::cSurfaceMeshFile MeshFile;

  //load the surface model
  for (int iSurfaceModelPart=0;iSurfaceModelPart<Orbiter::SurfaceModel::nTotalSurfaceModelFiles;iSurfaceModelPart++) {
    MeshFile.faceat=Orbiter::SurfaceModel::SurfaceModelSet[iSurfaceModelPart].faceat;
    sprintf(MeshFile.MeshFileName,"%s/%s",PIC::UserModelInputDataPath,Orbiter::SurfaceModel::SurfaceModelSet[iSurfaceModelPart].FileName);

    SurfaceMeshFileList.push_back(MeshFile);
  }

  switch (Orbiter::SurfaceModel::MeshFileFormat) {
  case Orbiter::SurfaceModel::MeshFileFormat_CEA:
    PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(SurfaceMeshFileList,Orbiter::SurfaceModel::ScalingFactor);
    break;

  case Orbiter::SurfaceModel::MeshFileFormat_NASTRAN:
    PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat(SurfaceMeshFileList,Orbiter::SurfaceModel::ScalingFactor);
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
      printf("  done.\n");
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
  PIC::Mesh::IrregularSurface::CutFaceAccessCounter::Init();
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

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {

    if (_ORBITER_EXECUTE_NEW_TEST_==_PIC_MODE_ON_) {
      DragCoefficientTest::ExecuteAllTests();
    }

    nTotalIterations=_PIC_NIGHTLY_TEST__TOTAL_ITERATION_NUMBER_;   
  } 

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
