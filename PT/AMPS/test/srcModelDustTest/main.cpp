//$Id$

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
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>


#include "pic.h"
#include "constants.h"
#include "Exosphere.h"
#include "Dust.h"
#include "Dust.dfn"

using namespace ElectricallyChargedDust::Charging;
using namespace std;

int main(int argc,char **argv) {
  //initialize MPI
  PIC::InitMPI();

  rnd_seed();
  PIC::Init_BeforeParser();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //set distanc to Sun
  Exosphere::xSun_SO[0]=2.0*_AU_;

  //output file
  char TestFileName[400] = "";
  sprintf(TestFileName,"%s/test_model-dust.dat",PIC::OutputDataFileDirectory);

  ofstream OutFile;
  if(PIC::ThisThread==0)
    OutFile.open(TestFileName);

  //grain parameters
  double GrainElectricChargeInit = 0.0;
  double GrainElectricCharge, EquilibriumCharge;
  const int nSize=4;
  double GrainRadii[nSize] = {1E-7, 1E-6, 1E-5, 1E-4};
  //plasma parameters
  double PlasmaVelocity[3] = {4E5,0,0};
  double PlasmaTemperature = 1E6;
  double PlasmaNumberDensity = 1E7;
  double PlasmaPressure = 1.5*Kbol*PlasmaNumberDensity*PlasmaTemperature;
  //time step
  double dt=1E+0;
  //particle velocity
  double ParticleVelocity[3] = {0};

  for(int iSize=0; iSize<nSize; iSize++){
    double GrainRadius = GrainRadii[iSize];
    
    //equlibrium charging
    EquilibriumCharge = GrainElectricChargeInit;
    UpdateGrainCharge(dt, ParticleVelocity,
		      PlasmaVelocity, PlasmaPressure,
		      PlasmaTemperature, PlasmaNumberDensity, 
		      GrainRadius, EquilibriumCharge,
		      CHARGE_INTEGRATION_MODE__EQUILIBRIUM_POTENTIAL);
    
    if(PIC::ThisThread==0){
      OutFile<<"\nCharging of a grain with radius "<<GrainRadius<<" m"<<
	"-------------------------------\n"<<
      "Equilibrium charge:   "<<EquilibriumCharge<<endl;
    }

    //time dependent charging
    GrainElectricCharge=GrainElectricChargeInit;
    int iIter=0;
    do{
      UpdateGrainCharge(dt,ParticleVelocity, 
			PlasmaVelocity, PlasmaPressure, 
			PlasmaTemperature, PlasmaNumberDensity, 
			GrainRadius, GrainElectricCharge, 
			CHARGE_INTEGRATION_MODE__TIME_DEPENDENT);
      
      if(PIC::ThisThread==0)
	OutFile<<"Charging iteration "<<++iIter<<": " <<
	  GrainElectricCharge<<endl;
    }while (fabs(GrainElectricCharge/EquilibriumCharge-1) > 1E-5);
    
  }

  if(PIC::ThisThread==0)
    OutFile.close();
  
  //finish the run
  MPI_Finalize();
  //  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}





