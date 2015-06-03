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

  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;

  switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
  case _PIC_COUPLER_DATAFILE_READER_MODE__ICES_ :
    ICES::GetDomainLimit(xmin,xmax);
    PIC::Mesh::mesh.init(xmin,xmax,ICES::localResolution);
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
    PIC::CPLR::DATAFILE::BATSRUS::Init("3d__mhd_1_n00000001.idl");
    PIC::CPLR::DATAFILE::BATSRUS::GetDomainLimits(xmin,xmax);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

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

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //create the reference file with the extracted data
  PIC::CPLR::DATAFILE::SaveTestReferenceData(TestFileName);



  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}





