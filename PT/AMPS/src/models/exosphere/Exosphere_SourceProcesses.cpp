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

//init the model
void Exosphere::SourceProcesses::Init() {
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  SolarWindSputtering::Init();
#endif
}

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

double Exosphere::SourceProcesses::SolarWindSputtering::TypicalIonFluxSputteringRate(int spec) {
  return SolarWindSputtering_Yield[spec]*Pi*pow(Exosphere::Planet->Radius,2)*Exosphere::swNumberDensity_Typical*
      sqrt(pow(Exosphere::swVelocity_Typical[0],2)+pow(Exosphere::swVelocity_Typical[1],2)+pow(Exosphere::swVelocity_Typical[2],2));
}

double Exosphere::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate(int nElement,int *spec) {
  return GetSurfaceElementProductionRate(*spec,nElement,Planet);
}

double Exosphere::SourceProcesses::SolarWindSputtering::DefaultEnergyDistributionFunction(double e,int *spec) {
  static const double Eb=2.0*eV2J;
  static const double Tm=500.0*eV2J;

  double t=e/pow(e+Eb,3)*(1.0-sqrt((e+Eb)/Tm));

  return (t>0.0) ? t : 0.0;
}

//calculate the total and maximum source rates
void Exosphere::SourceProcesses::SolarWindSputtering::Init() {
  double ElementSourceRate,maxRatePerM2,t;
  int el,nTotalSurfaceElements,spec;

  if (Exosphere::Planet==NULL) return;

  nTotalSurfaceElements=Planet->GetTotalSurfaceElementsNumber();

  for (spec=0;spec<PIC::nTotalSpecies;spec++) for (el=0;el<nTotalSurfaceElements;el++) {
    ElementSourceRate=GetSurfaceElementProductionRate(spec,el,Exosphere::Planet);
    t=ElementSourceRate/Planet->GetSurfaceElementArea(el);

    if (t>maxLocalSourceRate[spec]) maxLocalSourceRate[spec]=t;
    SourceRate[spec]+=ElementSourceRate;
  }
}

/*--------------------------------- SOURCE: Impact Vaporization -----------------------------------*/
double Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
  return ImpactVaporization_SourceRate[spec]*pow(ImpactVaporization_HeliocentricDistance/Exosphere::xObjectRadial,ImpactVaporization_SourceRatePowerIndex);
}
