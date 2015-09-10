
//$Id$

/*
 * main.cpp
 *
 *  Created on: May 26, 2015
 *      Author: vtenishe
 */



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

//the particle class
#include "pic.h"
#include "constants.h"
#include "mars-ions.h"


void amps_init();
void amps_init_mesh();
void amps_time_step();



int main(int argc,char **argv) {
//  MPI_Init(&argc,&argv);

  amps_init_mesh();
  amps_init();

  //time step
  for (long int niter=0;niter<100000001;niter++) {
    PIC::TimeStep();
  }



  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}

