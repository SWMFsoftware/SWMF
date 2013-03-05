/*
 * Mercury_SourceProcesses.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"


/*----------------------------------------   PHOTON STIMULATED DESORPTION -------------------------------------------*/
cSingleVariableDistribution Exosphere::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution;
cSingleVariableDiscreteDistribution Exosphere::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution;

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::EnergyDistributionFunction(double e) {
  const double x=0.7;
  const double U=0.052;

  e/=eV2J;

  return x*(1+x)*e*pow(U,x)/pow(e+U,2.0+x);
}

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementSodiumProductionRate(int nElement) {
  return GetSurfaceElementProductionRate(_NA_SPEC_,nElement,Planet);
}

/*----------------------------------------   THERMAL DESORPTION -------------------------------------------*/
cSingleVariableDiscreteDistribution Exosphere::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution;

double Exosphere::SourceProcesses::ThermalDesorption::GetSurfaceElementSodiumProductionRate(int nElement) {
  return GetSurfaceElementProductionRate(_NA_SPEC_,nElement,Planet);
}

/*----------------------------------------   SOLAR WIND SPUTTERING -------------------------------------------*/
cSingleVariableDistribution Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistribution;
cSingleVariableDiscreteDistribution Exosphere::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution;

double Exosphere::SourceProcesses::SolarWindSputtering::GetSurfaceElementSodiumProductionRate(int nElement) {
  return GetSurfaceElementProductionRate(_NA_SPEC_,nElement,Planet);
}

double Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistributionFunction(double e) {
  static const double Eb=2.0*eV2J;
  static const double Tm=500.0*eV2J;

  double t=e/pow(e+Eb,3)*(1.0-sqrt((e+Eb)/Tm));

  return (t>0.0) ? t : 0.0;
}

/*--------------------------------- SOURCE: Impact Vaporization -----------------------------------*/
double Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(int spec,void *SphereDataPointer) {
  return SourceRate*pow(HeliocentricDistance/Exosphere::xObjectRadial,SourceRatePowerIndex);
}
