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


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"


//$Id$


void amps_init();
void amps_time_step();

int main(int argc,char **argv) {


  PIC::InitMPI();

  amps_init();


  for (long int niter=0;niter<100000001 /*000001*/;niter++) {

    amps_time_step();
  }


  char fname[400];

  sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
