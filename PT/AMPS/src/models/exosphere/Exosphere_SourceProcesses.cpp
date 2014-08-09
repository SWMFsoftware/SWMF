//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury_SourceProcesses.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"
#include "Exosphere.h"


/*----------------------------------------   PHOTON STIMULATED DESORPTION -------------------------------------------*/
cSingleVariableDistribution<int> Exosphere::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[PIC::nTotalSpecies];
cSingleVariableDiscreteDistribution<int> Exosphere::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[PIC::nTotalSpecies];

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::EnergyDistributionFunction(double e, int *spec) {
  static const double x=0.7;
  static const double U=0.052;

  e/=eV2J;

  return x*(1+x)*e*pow(U,x)/pow(e+U,2.0+x);
}

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate(int nElement,int *spec) {
  return GetSurfaceElementProductionRate(*spec,nElement,Planet);
}

/*----------------------------------------   THERMAL DESORPTION -------------------------------------------*/
cSingleVariableDiscreteDistribution<int> Exosphere::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[PIC::nTotalSpecies];

double Exosphere::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate(int nElement,int *spec) {
  return GetSurfaceElementProductionRate(*spec,nElement,Planet);
}

/*----------------------------------------   SOLAR WIND SPUTTERING -------------------------------------------*/
cSingleVariableDistribution<int> Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistribution[PIC::nTotalSpecies];
cSingleVariableDiscreteDistribution<int> Exosphere::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[PIC::nTotalSpecies];

double Exosphere::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate(int nElement,int *spec) {
  return GetSurfaceElementProductionRate(*spec,nElement,Planet);
}

double Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistributionFunction(double e,int *spec) {
  static const double Eb=2.0*eV2J;
  static const double Tm=500.0*eV2J;

  double t=e/pow(e+Eb,3)*(1.0-sqrt((e+Eb)/Tm));

  return (t>0.0) ? t : 0.0;
}

/*--------------------------------- SOURCE: Impact Vaporization -----------------------------------*/
double Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
  return ImpactVaporization_SourceRate[spec]*pow(ImpactVaporization_HeliocentricDistance/Exosphere::xObjectRadial,ImpactVaporization_SourceRatePowerIndex);
}
