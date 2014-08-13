//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//the user defined classes used by the mesher
//$Id$


#ifndef _EXOSPHERE_DEFAULT_SURFACE_SAMPLING_DATASTRUCTURE_
#define _EXOSPHERE_DEFAULT_SURFACE_SAMPLING_DATASTRUCTURE_



#include "specfunc.h"

class cInternalSphericalData_UserDefined {
public :
  int faceat;
  double *SamplingBuffer;

  //the maximum time step in the blocks that are intersected by the sphere
  double *maxIntersectedNodeTimeStep;


  //Collect surface area density on the surface
  double **SurfaceElementDesorptionFluxUP,**SurfaceElementAdsorptionFluxDOWN,**SurfaceElementPopulation;

  //surface element normal and area
  struct cSurfaceElementExternalNormal {
    double norm[3];
  };

  cSurfaceElementExternalNormal *SurfaceElementExternalNormal;
  double *SurfaceElementArea;

  //local production rate due to different soruce processes for each element
  class cElementSourceRate {
  public:
    double PhotoStimulatedDesorptionFlux,ThermalDesorptionFlux,SolarWindSputteringFlux;

    cElementSourceRate() {
      PhotoStimulatedDesorptionFlux=0.0;
      ThermalDesorptionFlux=0.0;
      SolarWindSputteringFlux=0.0;
    }
  };

  cElementSourceRate **ElementSourceRate;

  //source rates, return fluxes and surface content
  double ***SampleSpeciesSurfaceSourceRate;
  double **SampleSpeciesSurfaceAreaDensity,**SampleSpeciesSurfaceReturnFlux,**SampleSpeciesSurfaceInjectionFlux;

  //bulk velocity of returning particles
  double **SampleReturnFluxBulkSpeed;

  //the bulk velocity of injected particles
  double **SampleInjectedFluxBulkSpeed;

  //solar wind wurface flux
  double *SolarWindSurfaceFlux;

  //injection of particles from the sphere
  typedef double (*fInjectionRate)(int spec,int BoundaryElementType,void *BoundaryElement);
  fInjectionRate InjectionRate;

  //injection of particles from the sphere
  typedef long int (*fInjectionBoundaryCondition)(int BoundaryElementType,void *BoundaryElement);
  fInjectionBoundaryCondition InjectionBoundaryCondition;

  //interaction of aprticles with the sphere
  typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal,void *startNode,void *InternalSphere);
  fParticleSphereInteraction ParticleSphereInteraction;

  //the flag indicating that the boundary prosidure associated with the sphere is already executed
  bool ProcessedBCflag;

  cInternalSphericalData_UserDefined () {
    faceat=-1,SamplingBuffer=NULL,ProcessedBCflag=false;
    InjectionRate=NULL,InjectionBoundaryCondition=NULL,ParticleSphereInteraction=NULL;
    maxIntersectedNodeTimeStep=NULL;

    SurfaceElementDesorptionFluxUP=NULL,SurfaceElementAdsorptionFluxDOWN=NULL,SurfaceElementPopulation=NULL;
    SurfaceElementExternalNormal=NULL,SurfaceElementArea=NULL;

    ElementSourceRate=NULL;

    SampleSpeciesSurfaceSourceRate=NULL,SampleSpeciesSurfaceAreaDensity=NULL,SampleSpeciesSurfaceReturnFlux=NULL,SampleSpeciesSurfaceInjectionFlux=NULL,SampleReturnFluxBulkSpeed=NULL,SampleInjectedFluxBulkSpeed=NULL;
    SolarWindSurfaceFlux=NULL;
  }

  void SaveSurfaceDensity(char *fname,int nSurfaceElements,int nSpecies) {
    FILE *fout=fopen(fname,"wb");
    fwrite(SurfaceElementPopulation[0],sizeof(double),nSurfaceElements*nSpecies,fout);
    fwrite(ElementSourceRate,sizeof(cElementSourceRate),nSurfaceElements,fout);

    fclose(fout);
  }

  void LoadSurfaceDensity(char *fname,int nSurfaceElements,int nSpecies) {
    FILE *fout=fopen(fname,"r");

    if ((SurfaceElementPopulation==NULL)||(ElementSourceRate==NULL)) exit(__LINE__,__FILE__,"Error: bufefr is not initialized");

    fread(SurfaceElementPopulation[0],sizeof(double),nSurfaceElements*nSpecies,fout);
    fread(ElementSourceRate,sizeof(cElementSourceRate),nSurfaceElements,fout);

    fclose(fout);
  }
};


class cInternalRotationBodyData_UserDefined : public cInternalSphericalData_UserDefined {
public:
  cInternalRotationBodyData_UserDefined() : cInternalSphericalData_UserDefined() {

  }
};

class cInternalNastranSurfaceData_UserDefined : public cInternalSphericalData_UserDefined {
public:
  cInternalNastranSurfaceData_UserDefined() : cInternalSphericalData_UserDefined() {
  }
};
#endif
