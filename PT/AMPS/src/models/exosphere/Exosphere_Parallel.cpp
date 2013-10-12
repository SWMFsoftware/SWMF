//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury_Parallel.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: vtenishe
 */

//$Id$


#include "pic.h"
#include "Exosphere.h"

void Exosphere::ExchangeSurfaceAreaDensity() {

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
  #if _EXOSPHERE_INTEGRATION_MODE_ == _EXOSPHERE_INTEGRATION_MODE__TIME_DEPENDENT_
    int el;

    double TotalFlux_LOCAL[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
    double TotalFlux_GLOBAL[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

    cInternalSphericalData *Sphere=(cInternalSphericalData*)(PIC::Mesh::mesh.InternalBoundaryList.begin()->BoundaryElement);

    memcpy(TotalFlux_LOCAL,Sphere->SurfaceElementAdsorptionFluxDOWN[spec],PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
    memcpy(TotalFlux_GLOBAL,Sphere->SurfaceElementDesorptionFluxUP[spec],PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) TotalFlux_LOCAL[el]-=TotalFlux_GLOBAL[el];

    MPI_Allreduce(TotalFlux_LOCAL,TotalFlux_GLOBAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  #if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
      exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
  #endif




      //update the value of the local sodium surface density
      double TotalPhotonStimulatedDesorptionFlux=0.0,localMaxPhotonStimulatedDesorptionFlux=0.0;
      double TotalThermalDesorptionFlux=0.0,localMaxThermalDesorptionFlux=0.0;
      double TotalSolarWindSputteringFlux=0.0,localMaxSolarWindSputteringFlux=0.0;
      double SurfaceElementArea[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      memcpy(SurfaceElementArea,Sphere->SurfaceElementArea,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));

      //redistribute the surface content
  #if _EXOSPHERE__SURFACE_CONTENT_ == _EXOSPHERE__SURFACE_CONTENT__BALANCE_FLUXES_
      //do noting
  #elif _EXOSPHERE__SURFACE_CONTENT_ == _EXOSPHERE__SURFACE_CONTENT__UNIFORM_
      //uniform surface density distribution
      double TotalSurfaceContent=0.0,TotalSurfaceArea=0.0;

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        TotalSurfaceContent+=Sphere->SodiumSurfaceElementPopulation[el]+TotalFlux_GLOBAL[el];
        TotalSurfaceArea+=SurfaceElementArea[el];
      }

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        Sphere->SodiumSurfaceElementPopulation[el]=0.0;
        TotalFlux_GLOBAL[el]=TotalSurfaceContent*SurfaceElementArea[el]/TotalSurfaceArea;
      }
  #elif _EXOSPHERE__SURFACE_CONTENT_ == _EXOSPHERE__SURFACE_CONTENT__RADIAL_DISTRIBUTION_
      //uniform surface density along azimuth on the sphere
      double TotalSurfaceContent,TotalSurfaceArea;
      int iZenith,iAzimuth;

      for (iZenith=0;iZenith<Exosphere::Planet->nZenithSurfaceElements;iZenith++) {
        TotalSurfaceContent=0.0,TotalSurfaceArea=0.0;

        for (iAzimuth=0;iAzimuth<Exosphere::Planet->nAzimuthalSurfaceElements;iAzimuth++) {
          el=Exosphere::Planet->GetLocalSurfaceElementNumber(iZenith,iAzimuth);

          TotalSurfaceContent+=Sphere->SodiumSurfaceElementPopulation[el]+TotalFlux_GLOBAL[el];
          TotalSurfaceArea+=SurfaceElementArea[el];
        }

        for (iAzimuth=0;iAzimuth<Exosphere::Planet->nAzimuthalSurfaceElements;iAzimuth++) {
          el=Exosphere::Planet->GetLocalSurfaceElementNumber(iZenith,iAzimuth);

          Sphere->SodiumSurfaceElementPopulation[el]=0.0;
          TotalFlux_GLOBAL[el]=TotalSurfaceContent*SurfaceElementArea[el]/TotalSurfaceArea;
        }
      }
  #elif _EXOSPHERE__SURFACE_CONTENT_ == _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_
      //the value of the surface content will be set up in the next loop
      //do nothing here
  #else
      exit(__LINE__,__FILE__,"Error: the option is not relcognized");
  #endif

      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        register double t;


  //====================   DEBUG ========================
  /*
  if (Sphere->SodiumSurfaceElementPopulation[el]+TotalFlux_GLOBAL[el]<0.0) {
    cout << __LINE__ << "  " << __FILE__ << endl;
  }
  */
  //====================  END DEBIG =====================



        Sphere->SurfaceElementPopulation[spec][el]+=TotalFlux_GLOBAL[el];

  #if _EXOSPHERE__SURFACE_CONTENT_ == _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_
        //use the user defined function for the local surface content
        Sphere->SurfaceElementPopulation[spec][el]=SurfaceElementArea[el]*(_EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(spec,el));
  #endif




    #if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
        t=Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate(spec,el,Sphere);

        TotalPhotonStimulatedDesorptionFlux+=t;

        t/=SurfaceElementArea[el];
        Sphere->ElementSourceRate[spec][el].PhotoStimulatedDesorptionFlux=t;
        if (localMaxPhotonStimulatedDesorptionFlux<t) localMaxPhotonStimulatedDesorptionFlux=t;
    #endif


    #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
        t=Exosphere::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate(spec,el,Sphere);

        TotalThermalDesorptionFlux+=t;

        t/=SurfaceElementArea[el];
        Sphere->ElementSourceRate[spec][el].ThermalDesorptionFlux=t;
        if (localMaxThermalDesorptionFlux<t) localMaxThermalDesorptionFlux=t;
    #endif

    #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
        t=Exosphere::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate(spec,el,Sphere);

        TotalSolarWindSputteringFlux+=t;

        t/=SurfaceElementArea[el];
        Sphere->ElementSourceRate[spec][el].SolarWindSputteringFlux=t;
        if (localMaxSolarWindSputteringFlux<t) localMaxSolarWindSputteringFlux=t;
    #endif


      }

      //update the total source rate values for specific source processes
      Exosphere::SourceProcesses::PhotonStimulatedDesorption::SourceRate[spec]=TotalPhotonStimulatedDesorptionFlux;
      Exosphere::SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate[spec]=localMaxPhotonStimulatedDesorptionFlux;

  #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      Exosphere::SourceProcesses::ThermalDesorption::SourceRate[spec]=TotalThermalDesorptionFlux;
      Exosphere::SourceProcesses::ThermalDesorption::maxLocalSourceRate[spec]=localMaxThermalDesorptionFlux;
  #endif

  #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      Exosphere::SourceProcesses::SolarWindSputtering::SourceRate[spec]=TotalSolarWindSputteringFlux;
      Exosphere::SourceProcesses::SolarWindSputtering::maxLocalSourceRate[spec]=localMaxSolarWindSputteringFlux;
  #endif

  #elif _EXOSPHERE_INTEGRATION_MODE_ == _EXOSPHERE_INTEGRATION_MODE__STEADY_STATE_
      //do nothing
  #else
      exit(__LINE__,__FILE__,"Error: the option in not recognized");
  #endif


      //flush the counters for the UP and DOWN fluxes
      for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) TotalFlux_LOCAL[el]=0.0;
      memcpy(Sphere->SurfaceElementDesorptionFluxUP[spec],TotalFlux_LOCAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
      memcpy(Sphere->SurfaceElementAdsorptionFluxDOWN[spec],TotalFlux_LOCAL,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*sizeof(double));
  }
}
