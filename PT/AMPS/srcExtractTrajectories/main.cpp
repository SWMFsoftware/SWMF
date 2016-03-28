//$Id$


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


//#include "meshAMRcutcell.h"
//#include "cCutBlockSet.h"
//#include "Comet.h"
//#include "Exosphere.h"


int main(int argc,char **argv) {
  PIC::InitMPI(); 

  int TrajectoryPointBufferLength=10000000;
  char fname[500]="PT/plots/amps.TrajectoryTracking.out=2";

  printf("Extracting particle trajectories....  ");

  PIC::ParticleTracker::CreateTrajectoryOutputFiles(fname,PIC::OutputDataFileDirectory,TrajectoryPointBufferLength);

  printf("done.\n");
  MPI_Finalize();

  return EXIT_SUCCESS;
}
