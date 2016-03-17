

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


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "Titan.h"


void amps_init();
void amps_time_step();


int main(int argc,char **argv) {
  //      MPI_Init(&argc,&argv);



  amps_init();


  //determine the total number of the iterations to perform
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;
  static int LastDataOutputFileNumber=0;


  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=100;


  //time step
  for (long int niter=0;niter<nTotalIterations;niter++) {
    amps_time_step();

     //print the iteration number
     if (PIC::Mesh::mesh.ThisThread==0) {
       time_t TimeValue=time(NULL);
       tm *ct=localtime(&TimeValue);

       printf(": (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
     }

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>100000) PIC::RequiredSampleLength=100000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }

  }


  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_Titan.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }


  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
