//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//the interface between AMPS and SWMF

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <sstream>

#include <sys/time.h>
#include <sys/resource.h>

#include "pic.h"
#include "amps2swmf.h"


using namespace std;


void amps_init();
void amps_init_mesh();
void amps_time_step();


extern "C" { 
  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit);
  int initamps_();
  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc);
  void amps_finalize_();

  //import magnetic field from GM onto the 'center' nodes
  void amps_get_center_point_number(int*);
  void amps_get_center_point_coordinates(double*);

  //return the number of the AMPS' mesh rebalancing operations
  void amps_mesh_id_(int* id) {
    *id=PIC::Mesh::mesh.nParallelListRedistributions;
  }

  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc) {
    PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(iComm,iProc,nProc);

    //initialize the coupler and AMPS
    PIC::CPLR::SWMF::init();
    amps_init_mesh();
  }


  void amps_get_center_point_number_(int *nCenterPoints) {
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints);
  }

  void amps_get_center_point_coordinates_(double *x) {
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x);
  }


  void amps_recieve_gm2amps_center_point_data_(char *NameVar, int *nVar, double *data,int *index) {
    PIC::CPLR::SWMF::RecieveCenterPointData(NameVar,*nVar,data,index);
  }


  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit) {
    static bool InitFlag=false;
    static double swmfTimeSimulation=-1.0;

    if (swmfTimeSimulation<0.0) swmfTimeSimulation=*TimeSimulation;

    //call AMPS only after the first coupling has occured
    if (PIC::CPLR::SWMF::FirstCouplingOccured==false) {
      *TimeSimulation=*TimeSimulationLimit;
      return;
    }

    //init AMPS
    if (InitFlag==false) {
      //initamps_();

      amps_init();
      InitFlag=true;

      //print the output file on each iteration
      //PIC::RequiredSampleLength=1;
    }

    //determine whether to proceed with the current iteraction
    if (swmfTimeSimulation+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]>*TimeSimulationLimit) {
      *TimeSimulation=*TimeSimulationLimit;
      return;
    }
    else {
      swmfTimeSimulation+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
      *TimeSimulation=swmfTimeSimulation;
    }

    //call AMPS
    static long int counter=0;
    counter++;

    amps_time_step();

    if (PIC::ModelTestRun::mode==true) if (counter==PIC::ModelTestRun::nTotalIteraction) {
      char fname[400];

      sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

      exit(0);
    }

  }

  void amps_finalize_() {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

    //save particle trajectory file
    #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    sprintf(fname,"%s/amps.TrajectoryTracking.out=Final",PIC::OutputDataFileDirectory);
    PIC::ParticleTracker::OutputTrajectory(fname);
   #endif
  }


  int amps_read_param_(char *param, int *nlines, int *ncharline, int *iProc){
    // convert character array to string stream object                                                                                                      
    std::stringstream ss;

    AMPS2SWMF::PARAMIN::char_to_stringstream(param, *nlines, *ncharline,&ss);
    AMPS2SWMF::PARAMIN::read_paramin(&ss);

    return 0;
  }


}

/*
int main () {
  return 1;
*/
