
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


//#include <saturn.h>

//#define _MAIN_PROFILE_



//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"


void amps_init();
void amps_time_step();





int main(int argc,char **argv) {
//	MPI_Init(&argc,&argv);

  amps_init();




	//time step
	for (long int niter=0;niter<100000001;niter++) {

	  amps_time_step();




	}


	cout << "End of the run:" << PIC::nTotalSpecies << endl;

	return 1;
}
