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

#include <sys/time.h>
#include <sys/resource.h>

#include "pic.h"

using namespace std;


void amps_init();
void amps_time_step();


extern "C" { 
  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit);
  int initamps_();
  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc);

  //import magnetic field from GM onto the 'center' nodes
  void amps_get_center_point_number(int*);
  void amps_get_center_point_coordinates(double*);

  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc) {
    PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(iComm,iProc,nProc);

    //initialize the coupler and AMPS
    PIC::CPLR::SWMF::init();
    amps_init();
  }


  void amps_get_center_point_number_(int *nCenterPoints) {
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints);
  }

  void amps_get_center_point_coordinates_(double *x) {
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x);
  }


  void amps_recieve_gm2amps_center_point_data_(double *data,int *index) {
    PIC::CPLR::SWMF::RecieveCenterPointData(data,index);
  }


  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit) {
    static bool InitFlag=false;

    if (InitFlag==false) {
      //initamps_();

      //amps_init();

      InitFlag=true;

      //print the output file on each iteration
      PIC::RequiredSampleLength=1;
    }


    static int counter=0;
    counter++;

    amps_time_step();

    if (counter==100) {
      char fname[400];

      sprintf(fname,"%s/amsp.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

      exit(0);
    }    

  }

}

/*
int main () {
  return 1;
*/
