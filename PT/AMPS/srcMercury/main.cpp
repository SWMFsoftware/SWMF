

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
  //      MPI_Init(&argc,&argv);

//  PIC::InitMPI();
//  PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::Init("3d__mhd_1_n00000001.idl"); //"3d__all_3_t00000010_n0000059.idl"); //"3d__mhd_1_n00000001.idl");

cout << "Start: MERCURY" << endl;

  amps_init_mesh();
  amps_init();




  //time step
  for (long int niter=0;niter<100000001;niter++) {

    amps_time_step();




  }


  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
