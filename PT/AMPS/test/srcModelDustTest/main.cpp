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

int main(int argc,char **argv) {
  //initialize MPI
  PIC::InitMPI();

  rnd_seed();
  PIC::Init_BeforeParser();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //set distanc to Sun
  Exosphere::xSun_SO[0]=1.0*_AU_;

  //grain parameters
  double GrainElectricChargeInit = 0.0;
  double GrainElectricCharge;
  double GrainRadius = 1E-6;
  //plasma parameters
  double PlasmaVelocity[3] = {4E5,0,0};
  double PlasmaTemperature = 100;
  double PlasmaNumberDensity = 1E7;
  double PlasmaPressure = 1.5*Kbol*PlasmaNumberDensity*PlasmaTemperature;
  //time step
  double dt=1E-3;
  //particle velocity
  double ParticleVelocity[3] = {0};

  //equlibrium charging
  GrainElectricCharge = GrainElectricChargeInit;
  ElectricallyChargedDust::Charging::UpdateGrainCharge(dt, ParticleVelocity,
     PlasmaVelocity, PlasmaPressure,
     PlasmaTemperature, PlasmaNumberDensity, 
     GrainRadius, GrainElectricCharge,
     ElectricallyChargedDust::Charging::CHARGE_INTEGRATION_MODE__EQUILIBRIUM_POTENTIAL);

  if(PIC::ThisThread==0)
    std::cout<<GrainElectricCharge<<std::endl;

  //time dependent charging
  GrainElectricCharge=GrainElectricChargeInit;
  for(int i=0; i<10;i++){
    ElectricallyChargedDust::Charging::UpdateGrainCharge(dt,ParticleVelocity, PlasmaVelocity, PlasmaPressure, PlasmaTemperature, PlasmaNumberDensity, GrainRadius, GrainElectricCharge, ElectricallyChargedDust::Charging::CHARGE_INTEGRATION_MODE__TIME_DEPENDENT);

  if(PIC::ThisThread==0)
    std::cout<<GrainElectricCharge<<std::endl;
  }
  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}





