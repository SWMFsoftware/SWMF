


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

  amps_init();


  int niter;
  const int niterMax=100;

  for (niter=0;niter<niterMax;niter++) {

    amps_time_step();
  }


  char fname[400];

  sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
