//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf



//the particle class
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

//$Id$


void amps_init();
void amps_time_step();

int main(int argc,char **argv) {
  
  clock_t runtime =-clock();
  
  amps_init();


  int niter,nTotalIterations=100000001;
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=100;

  //time step
  static int LastDataOutputFileNumber=0;

  for (niter=0;niter<nTotalIterations;niter++) {

    amps_time_step();

    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>40000) PIC::RequiredSampleLength=40000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }

  }


  char fname[400];

  sprintf(fname,"%s/test_Moon.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

//  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  MPI_Finalize();

  runtime+=clock();
  
  if(PIC::Mesh::mesh.ThisThread==0)
    cout << "Total AMPS runtime is "
	 << (double)runtime / CLOCKS_PER_SEC << endl;

  return EXIT_SUCCESS;
}
