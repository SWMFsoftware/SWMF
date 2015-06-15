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

//the parameters of the domain and the sphere
//setting of ICES

namespace ICES {
  const double DebugRunMultiplier=8.0;
  const double rSphere=_ENCELADUS__RADIUS_;
  const double xMaxDomain=5.0;
  const double dxMinGlobal=1.0,dxMaxGlobal=5.0;
  const double dxMinSphere=DebugRunMultiplier*2.0/100,dxMaxSphere=DebugRunMultiplier*4.0/100.0;


  //Preprocessor of the SWMF data the obtained through ICES -> rotate the vectors
  //Need to rotate: E[3],B[3],swVel[3]
  void SWMFdataPreProcessor(double *x,PIC::CPLR::DATAFILE::ICES::cDataNodeSWMF& data) {
    VectorRotation::Along_Z_direction(data.B,-270.0/180.0*Pi);
    VectorRotation::Along_Z_direction(data.E,-270.0/180.0*Pi);
    VectorRotation::Along_Z_direction(data.swVel,-270.0/180.0*Pi);
  }

  //the mesh resolution
  double localResolution(double *x) {
    int idim;
    double res=0.0,r=0.0;

    if ((x[2]<0.0)&&(asin(sqrt((x[0]*x[0]+x[1]*x[1])/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])))/Pi*180.0<45.0)) return DebugRunMultiplier*0.25*(dxMaxSphere+dxMinSphere)*rSphere;

    return DebugRunMultiplier*rSphere*dxMaxGlobal;
  }

  //get the domain limits used in the test
  void GetDomainLimit(double *xmin,double *xmax) {
    for (int idim=0;idim<DIM;idim++) {
      xmin[idim]=-xMaxDomain*rSphere;
      xmax[idim]=xMaxDomain*rSphere;
    }
  }

  void Read() {
    PIC::CPLR::DATAFILE::ICES::SWMFdataPreProcessor=SWMFdataPreProcessor;

    PIC::CPLR::DATAFILE::ICES::createCellCenterCoordinateList();
    PIC::CPLR::DATAFILE::ICES::retriveSWMFdata("Enceladus");

    PIC::CPLR::DATAFILE::ICES::readSWMFdata(1.0);
  }
}

namespace BATL {
  const double rSphere=_MERCURY__RADIUS_;
  double *xmin,*xmax;

  double localResolution(double *x) {
    double l=0.0;

    for (int idim=0;idim<3;idim++) l+=pow(xmax[idim]-xmin[idim],2);

    return sqrt(l)/20;
  }
}

namespace TECPLOT {
  const double rSphere=_ENCELADUS__RADIUS_;
  double *xmin,*xmax;

  double localResolution(double *x) {
    double l=0.0;

    for (int idim=0;idim<3;idim++) l+=pow(xmax[idim]-xmin[idim],2);

    return sqrt(l)/20;
  }
}

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {return 1.0;}

int main(int argc,char **argv) {
  PIC::InitMPI();

  rnd_seed();
  PIC::Init_BeforeParser();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //init the mesh
  char TestFileName[400];
  double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

  if (PIC::ThisThread==0) cout << "Init the mesh" << endl;

  //Init the mesh
  PIC::Mesh::mesh.AllowBlockAllocation=false;

  //init the datafile reader
  switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
  case _PIC_COUPLER_DATAFILE_READER_MODE__ICES_ :
    ICES::GetDomainLimit(xmin,xmax);
    PIC::Mesh::mesh.init(xmin,xmax,ICES::localResolution);
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
    PIC::CPLR::DATAFILE::BATSRUS::Init("3d__mhd_1_n00000001.idl");
    PIC::CPLR::DATAFILE::BATSRUS::GetDomainLimits(xmin,xmax);
    PIC::CPLR::DATAFILE::BATSRUS::UnitLength=BATL::rSphere;

    //convert the BATSRUS coordinate units into that of AMPS
    for (int idim=0;idim<DIM;idim++) xmin[idim]*=BATL::rSphere,xmax[idim]*=BATL::rSphere; 
    BATL::xmin=xmin,BATL::xmax=xmax;
 
    //init the mesh object
    PIC::Mesh::mesh.init(xmin,xmax,BATL::localResolution);
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
    {
      double xminTECPLOT[3]={-5.1,-5.1,-5.1},xmaxTECPLOT[3]={5.1,5.1,5.1};
      double RotationMatrix_BATSRUS2AMPS[3][3]={ { 1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

      //  1  0  0
      //  0  1  0
      //  0  0  1

      PIC::CPLR::DATAFILE::TECPLOT::SetRotationMatrix_DATAFILE2LocalFrame(RotationMatrix_BATSRUS2AMPS);

      PIC::CPLR::DATAFILE::TECPLOT::UnitLength=_EUROPA__RADIUS_;
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsXYZ(xminTECPLOT,xmaxTECPLOT);
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsSPHERICAL(1.001,10.0);

      PIC::CPLR::DATAFILE::TECPLOT::DataMode=PIC::CPLR::DATAFILE::TECPLOT::DataMode_SPHERICAL;
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedVelocityVariableData(5,1.0E3);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedIonPressureVariableData(11,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedMagneticFieldVariableData(8,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedDensityVariableData(4,1.0E6);
      PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=11;
    }


    //init the mesh object
    for (int idim=0;idim<3;idim++) xmin[idim]=-5.0*_EUROPA__RADIUS_,xmax[idim]=5.0*_EUROPA__RADIUS_;
    PIC::Mesh::mesh.init(xmin,xmax,TECPLOT::localResolution);

    break;

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //generate AMPS' mesh
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

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

  //initialize the blocks
  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();
  PIC::Mesh::mesh.InitCellMeasure();

  PIC::Init_AfterParser();


  switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
  case _PIC_COUPLER_DATAFILE_READER_MODE__ICES_ :
    ICES::Read();
    sprintf(TestFileName,"%s/test_ices-reader.dat",PIC::OutputDataFileDirectory);
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
    PIC::CPLR::DATAFILE::BATSRUS::LoadDataFile();
    sprintf(TestFileName,"%s/test_batl-reader.dat",PIC::OutputDataFileDirectory);
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
    PIC::CPLR::DATAFILE::TECPLOT::ImportData("3d__mhd_3_n00045039-extracted.plt");
    sprintf(TestFileName,"%s/test_tecplot-reader.dat",PIC::OutputDataFileDirectory);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //create the reference file with the extracted data
  PIC::CPLR::DATAFILE::SaveTestReferenceData(TestFileName);

  //output ascii file with the interpolated data
  PIC::Mesh::cDataBlockAMR block;

  block.SetLocalParticleWeight(0.0,0);
  block.SetLocalTimeStep(0.0,0);

  sprintf(TestFileName,"%s/loaded-data.dat",PIC::OutputDataFileDirectory);
  PIC::Mesh::mesh.outputMeshDataTECPLOT(TestFileName,0);

  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}





