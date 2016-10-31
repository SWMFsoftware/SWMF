
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



//the particle class
#include "pic.h"
#include "constants.h"
#include "Earth.h"


void amps_init();
void amps_init_mesh();
void amps_time_step();

int main(int argc,char **argv) {
  static int LastDataOutputFileNumber=0;

  //init Earth magnetosphere model
  Earth::Init();

  amps_init_mesh();
  amps_init();




  PIC::RequiredSampleLength=600;


#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  int nTotalIterations = 0;
#else
  int nTotalIterations = 100000001;
#endif

  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
    
    amps_time_step();
    
    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);
      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",
	     ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,
	     PIC::RequiredSampleLength,
	     PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>25000) PIC::RequiredSampleLength=25000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }

  }
  
  
  //output the particle statistics of the test run
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  char fname[300];
  sprintf(fname,"%s/test_Earth.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
#endif
  
  
  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  
  return EXIT_SUCCESS;
}
