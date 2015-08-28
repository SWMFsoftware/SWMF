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



#include "pic.h"
#include "constants.h"


void amps_init();
void amps_init_mesh();
int  amps_time_step();


int main(int argc,char **argv) {

  amps_init_mesh();
  amps_init();

  //time step
  for (long int niter=0;niter<100000001;niter++) {

    if(amps_time_step() == _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_) break;

  }
  
  //output the particle statistics for the nightly tests                        
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_SEP3D.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }

  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
