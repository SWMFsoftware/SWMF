/*
 * Europa_Parallel.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: vtenishe
 */

//$Id$


#include "pic.h"

void Europa::ExchangeSurfaceAreaDensity() {

/*

  if (Europa::Planet==NULL) return;

#if _EUROPA_INTEGRATION_MODE_ == _EUROPA_INTEGRATION_MODE__TIME_DEPENDENT_
  int el;

  double TotalFlux_LOCAL[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
  double TotalFlux_GLOBAL[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

  cInternalSphericalData *Sphere=(cInternalSphericalData*)(PIC::Mesh::mesh.InternalBoundaryList.begin()->BoundaryElement);

  memcpy(TotalFlux_LOCAL,Sphere->O2SurfaceAreaDensity_FluxDOWN,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
  memcpy(TotalFlux_GLOBAL,Sphere->O2SurfaceAreaDensity_FluxUP,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
  for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) TotalFlux_LOCAL[el]-=TotalFlux_GLOBAL[el];

  MPI_Allreduce(TotalFlux_LOCAL,TotalFlux_GLOBAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
#endif


  //flush the counters for the UP and DOWN fluxes
    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) TotalFlux_LOCAL[el]=0.0;
    memcpy(Sphere->O2SurfaceAreaDensity_FluxUP,TotalFlux_LOCAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
    memcpy(Sphere->O2SurfaceAreaDensity_FluxDOWN,TotalFlux_LOCAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));

    //update the value of the local O2 surface density
    double TotalPhotonStimulatedDesorptionFlux=0.0,localMaxPhotonStimulatedDesorptionFlux=0.0;
    double TotalThermalDesorptionFlux=0.0,localMaxThermalDesorptionFlux=0.0;
    double TotalSolarWindSputteringFlux=0.0,localMaxSolarWindSputteringFlux=0.0;
    double SurfaceElementArea[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    memcpy(SurfaceElementArea,Sphere->SurfaceElementArea,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));

    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      register double t;

      Sphere->O2SurfaceAreaDensity[el]+=TotalFlux_GLOBAL[el];

  #if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
      t=Europa::SourceProcesses::PhotonStimulatedDesorption::GetLocalProductionRate(_NA_SPEC_,el,Sphere);

      TotalPhotonStimulatedDesorptionFlux+=t;

      t/=SurfaceElementArea[el];
      Sphere->ElementSourceRate[el].PhotoStimulatedDesorptionFlux=t;
      if (localMaxPhotonStimulatedDesorptionFlux<t) localMaxPhotonStimulatedDesorptionFlux=t;
  #endif


  #if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
      t=Europa::SourceProcesses::ThermalDesorption::GetLocalProductionRate(_NA_SPEC_,el,Sphere);

      TotalThermalDesorptionFlux+=t;

      t/=SurfaceElementArea[el];
      Sphere->ElementSourceRate[el].ThermalDesorptionFlux=t;
      if (localMaxThermalDesorptionFlux<t) localMaxThermalDesorptionFlux=t;
  #endif

  #if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
      t=Europa::SourceProcesses::SolarWindSputtering::GetLocalProductionRate(_NA_SPEC_,el,Sphere);

      TotalSolarWindSputteringFlux+=t;

      t/=SurfaceElementArea[el];
      Sphere->ElementSourceRate[el].SolarWindSputteringFlux=t;
      if (localMaxSolarWindSputteringFlux<t) localMaxSolarWindSputteringFlux=t;
  #endif


    }

    //update the total source rate values for specific source processes
    Europa::SourceProcesses::PhotonStimulatedDesorption::SourceRate=TotalPhotonStimulatedDesorptionFlux;
    Europa::SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate=localMaxPhotonStimulatedDesorptionFlux;

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
    Europa::SourceProcesses::ThermalDesorption::SourceRate=TotalThermalDesorptionFlux;
    Europa::SourceProcesses::ThermalDesorption::maxLocalSourceRate=localMaxThermalDesorptionFlux;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
    Europa::SourceProcesses::SolarWindSputtering::SourceRate=TotalSolarWindSputteringFlux;
    Europa::SourceProcesses::SolarWindSputtering::maxLocalSourceRate=localMaxSolarWindSputteringFlux;
#endif

#elif _EUROPA_INTEGRATION_MODE_ == _EUROPA_INTEGRATION_MODE__STEADY_STATE_
    //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option in not recognized");
#endif
*/

}
