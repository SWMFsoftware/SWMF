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

  //void amps_init();
  //void amps_time_step();

extern "C" { 
  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit);
  int initamps_();
  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc);

  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc) {
    MPI_GLOBAL_COMMUNICATOR=MPI_Comm_f2c(*iComm);
    PIC::InitMPI();

    if (PIC::ThisThread==0) {
      printf("AMPS: MPI Communicatior is imported from SWMF, size=%i\n",PIC::nTotalThreads);
    }
  }

  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit) {
    static bool InitFlag=false;

    if (InitFlag==false) {
      //initamps_();

      //amps_init();

      InitFlag=true;
    }


    static int counter=0;
    counter++;

    //amps_time_step();

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
