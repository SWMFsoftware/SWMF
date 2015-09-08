
//$Id$

/*
 * main_lib.cpp
 *
 *  Created on: May 26, 2015
 *      Author: vtenishe
 */


#include "mars-ions.h"

double localSphericalSurfaceResolution(double *x);
double localResolution(double *x);
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  return res;
}



void amps_init_mesh() {
  PIC::InitMPI();

  //SetUp the alarm
  //  PIC::Alarm::SetAlarm(2000);

  rnd_seed();

  //init the physical model
  MarsIon::Init_BeforeParser();

  //init the particle solver
  PIC::Init_BeforeParser();

  MarsIon::Init_AfterParser();



  //register the sphere
  static const bool SphereInsideDomain=true;

  if (SphereInsideDomain==true) {
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
    Sphere->SetSphereGeometricalParameters(sx0,_RADIUS_(_TARGET_));

    Sphere->Radius=_RADIUS_(_TARGET_);
    Sphere->PrintSurfaceMesh("Sphere.dat");
    Sphere->PrintSurfaceData("SpheraData.dat",0);
    Sphere->localResolution=localSphericalSurfaceResolution;
//    Sphere->InjectionRate=Europa::SourceProcesses::totalProductionRate;
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=MarsIon::ParticleSphereInteraction;

//    Sphere->PrintTitle=Europa::Sampling::OutputSurfaceDataFile::PrintTitle;
//    Sphere->PrintVariableList=Europa::Sampling::OutputSurfaceDataFile::PrintVariableList;
//    Sphere->PrintDataStateVector=Europa::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

    //set up the planet pointer in Europa model
    MarsIon::Planet=Sphere;
    Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);
  }

  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  cout << "Init the mesh" << endl;

  double xmin[3],xmax[3];

  for (int idim=0;idim<DIM;idim++) {
    xmax[idim]=4*_RADIUS_(_TARGET_);
    xmin[idim]=-4*_RADIUS_(_TARGET_);
  }

  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,localResolution);
  PIC::Mesh::mesh.memoryAllocationReport();

  //generate mesh or read from file
  char mesh[200]="amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  bool NewMeshGeneratedFlag=false;

  FILE *fmesh=NULL;

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

  //allocate the mesh data buffers
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


  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
  if (PIC::ThisThread==0) cout << "AMPS' Initialization is complete" << endl;

}

void amps_init() {
   int idim;

   //init the PIC solver
   PIC::Init_AfterParser ();
   PIC::Mover::Init();

   //create the list of mesh nodes where the injection boundary conditions are applied
//   PIC::BC::BlockInjectionBCindicatior=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxParticleInjectionIndicator;
//   PIC::BC::userDefinedBoundingBlockInjectionFunction=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjection;
//   PIC::BC::InitBoundingBoxInjectionBlockList();




   //init the particle buffer
   PIC::ParticleBuffer::Init(10000000);
   int LastDataOutputFileNumber=-1;


   //init the sampling of the particls' distribution functions
   //const int nSamplePoints=3;
   //double SampleLocations[nSamplePoints][DIM]={{2.0E6,0.0,0.0}, {0.0,2.0E6,0.0}, {-2.0E6,0.0,0.0}};

   /* THE DEFINITION OF THE SAMPLE LOCATIONS IS IN THE INPUT FILE
      PIC::DistributionFunctionSample::vMin=-40.0E3;
      PIC::DistributionFunctionSample::vMax=40.0E3;
      PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

      PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);
   */

   //also init the sampling of the particles' pitch angle distribution functions
   //PIC::PitchAngleDistributionSample::nSampledFunctionPoints=101;

   //PIC::PitchAngleDistributionSample::Init(SampleLocations,nSamplePoints);






#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
#if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_
  //TECPLOT
  //read the background data
    if (PIC::CPLR::DATAFILE::BinaryFileExists("MARS-BATSRUS")==true)  {
      PIC::CPLR::DATAFILE::LoadBinaryFile("MARS-BATSRUS");
    }
    else {
      double xminTECPLOT[3]={-5.1,-5.1,-5.1},xmaxTECPLOT[3]={5.1,5.1,5.1};

      double RotationMatrix_BATSRUS2AMPS[3][3]={ { 1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

      //  1  0  0
      //  0  1  0
      //  0  0  1

      PIC::CPLR::DATAFILE::TECPLOT::SetRotationMatrix_DATAFILE2LocalFrame(RotationMatrix_BATSRUS2AMPS);

      PIC::CPLR::DATAFILE::TECPLOT::UnitLength=_MARS__RADIUS_;
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsXYZ(xminTECPLOT,xmaxTECPLOT);
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsSPHERICAL(1.001,10.0);

      PIC::CPLR::DATAFILE::TECPLOT::DataMode=PIC::CPLR::DATAFILE::TECPLOT::DataMode_SPHERICAL;
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedVelocityVariableData(39,1.0E3);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedIonPressureVariableData(26,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedMagneticFieldVariableData(8,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedDensityVariableData(22,1.0E6/16.);
      PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=41;
      PIC::CPLR::DATAFILE::TECPLOT::ImportData("data_mhd_PERmax-SSLONG180U.plt");

      PIC::CPLR::DATAFILE::SaveBinaryFile("MARS-BATSRUS");
    }

#else
    exit(__LINE__,__FILE__,"ERROR: unrecognized datafile reader mode");
#endif //_PIC_COUPLER_DATAFILE_READER_MODE_
#endif //_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_

    //init background data based on that loaded from TECPLOT
    MarsIon::InitBackgroundData();

    //set up the time step
    PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
    PIC::ParticleWeightTimeStep::initTimeStep();

    //set up the particle weight
    PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=MarsIon::SourceProcesses::GetBlockInjectionRate;

    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
      PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(spec);
    }

  PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);
}


