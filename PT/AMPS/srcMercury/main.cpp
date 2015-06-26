

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
#include "Mercury.h"


void amps_init();
void amps_init_mesh();
void amps_time_step();


int main(int argc,char **argv) {
  amps_init_mesh();
  amps_init();


  //create the reference file with the extracted data
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  char fname[400];

  sprintf(fname,"%s/test_Mercury_background.dat",PIC::OutputDataFileDirectory);
  PIC::CPLR::DATAFILE::SaveTestReferenceData(fname);
  #endif




  //time step
  int nTotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_OFF_) ? 100000001 : 150;

  for (long int niter=0;niter<nTotalIterations;niter++) {
    amps_time_step();
  }


  //output the particle statistics of the test run
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  sprintf(fname,"%s/test_Mercury.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  #endif


  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
