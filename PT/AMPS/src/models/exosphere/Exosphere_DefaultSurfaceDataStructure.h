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

  template <class T>
  void Allocate(int nTotalSpecies, int TotalSurfaceElementNumber, int EXOSPHERE__SOURCE_MAX_ID_VALUE,T* Surface) {
    //allocate the buffers for collecting the sodium surface density
    SurfaceElementDesorptionFluxUP=new double*[nTotalSpecies];
    SurfaceElementAdsorptionFluxDOWN=new double*[nTotalSpecies];
    SurfaceElementPopulation=new double*[nTotalSpecies];

    SurfaceElementDesorptionFluxUP[0]=new double[nTotalSpecies*TotalSurfaceElementNumber];
    SurfaceElementAdsorptionFluxDOWN[0]=new double[nTotalSpecies*TotalSurfaceElementNumber];
    SurfaceElementPopulation[0]=new double[nTotalSpecies*TotalSurfaceElementNumber];

    for (int spec=1;spec<nTotalSpecies;spec++) {
      SurfaceElementDesorptionFluxUP[spec]=SurfaceElementDesorptionFluxUP[spec-1]+TotalSurfaceElementNumber;
      SurfaceElementAdsorptionFluxDOWN[spec]=SurfaceElementAdsorptionFluxDOWN[spec-1]+TotalSurfaceElementNumber;
      SurfaceElementPopulation[spec]=SurfaceElementPopulation[spec-1]+TotalSurfaceElementNumber;
    }

    SurfaceElementExternalNormal=new cInternalSphericalData_UserDefined::cSurfaceElementExternalNormal[TotalSurfaceElementNumber];
    SurfaceElementArea=new double[TotalSurfaceElementNumber];


    ElementSourceRate=new cInternalSphericalData_UserDefined::cElementSourceRate*[nTotalSpecies];
    ElementSourceRate[0]=new cInternalSphericalData_UserDefined::cElementSourceRate[nTotalSpecies*TotalSurfaceElementNumber];

    for (int spec=1;spec<nTotalSpecies;spec++) {
      ElementSourceRate[spec]=ElementSourceRate[spec-1]+TotalSurfaceElementNumber;
    }


    SolarWindSurfaceFlux=new double[TotalSurfaceElementNumber];

    for (int el=0;el<TotalSurfaceElementNumber;el++) {
      for (int spec=0;spec<nTotalSpecies;spec++) {
        SurfaceElementDesorptionFluxUP[spec][el]=0.0;
        SurfaceElementAdsorptionFluxDOWN[spec][el]=0.0;
        SurfaceElementPopulation[spec][el]=0.0;
      }

      SolarWindSurfaceFlux[el]=-1.0;

      SurfaceElementArea[el]=Surface->GetSurfaceElementArea(el);
      Surface->GetSurfaceElementNormal((SurfaceElementExternalNormal+el)->norm,el);
    }

    //allocate buffers for sampling surface sodium source rates and sodikum surface content
    int offsetSpecie,offsetElement,s,el,i;

    SampleSpeciesSurfaceSourceRate=new double** [nTotalSpecies];
    SampleSpeciesSurfaceSourceRate[0]=new double *[nTotalSpecies*TotalSurfaceElementNumber];
    SampleSpeciesSurfaceSourceRate[0][0]=new double [nTotalSpecies*TotalSurfaceElementNumber*(EXOSPHERE__SOURCE_MAX_ID_VALUE+1)];

    for (offsetSpecie=0,s=0,offsetElement=0;s<nTotalSpecies;s++) {
      SampleSpeciesSurfaceSourceRate[s]=SampleSpeciesSurfaceSourceRate[0]+offsetSpecie;
      offsetSpecie+=TotalSurfaceElementNumber;

      for (el=0;el<TotalSurfaceElementNumber;el++) {
        SampleSpeciesSurfaceSourceRate[s][el]=SampleSpeciesSurfaceSourceRate[0][0]+offsetElement;
        offsetElement+=EXOSPHERE__SOURCE_MAX_ID_VALUE+1;

        for (i=0;i<EXOSPHERE__SOURCE_MAX_ID_VALUE+1;i++) SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;
      }
    }

    SampleSpeciesSurfaceAreaDensity=new double* [nTotalSpecies];
    SampleSpeciesSurfaceAreaDensity[0]=new double [nTotalSpecies*TotalSurfaceElementNumber];

    for (offsetSpecie=0,s=0;s<nTotalSpecies;s++) {
      SampleSpeciesSurfaceAreaDensity[s]=SampleSpeciesSurfaceAreaDensity[0]+offsetSpecie;
      offsetSpecie+=TotalSurfaceElementNumber;

      for (el=0;el<TotalSurfaceElementNumber;el++) {
        SampleSpeciesSurfaceAreaDensity[s][el]=0.0;
      }
    }

    SampleSpeciesSurfaceReturnFlux=new double* [nTotalSpecies];
    SampleSpeciesSurfaceReturnFlux[0]=new double [nTotalSpecies*TotalSurfaceElementNumber];

    SampleSpeciesSurfaceInjectionFlux=new double* [nTotalSpecies];
    SampleSpeciesSurfaceInjectionFlux[0]=new double [nTotalSpecies*TotalSurfaceElementNumber];

    SampleReturnFluxBulkSpeed=new double* [nTotalSpecies];
    SampleReturnFluxBulkSpeed[0]=new double [nTotalSpecies*TotalSurfaceElementNumber];

    SampleInjectedFluxBulkSpeed=new double* [nTotalSpecies];
    SampleInjectedFluxBulkSpeed[0]=new double [nTotalSpecies*TotalSurfaceElementNumber];

    for (offsetSpecie=0,s=0;s<nTotalSpecies;s++) {
      SampleSpeciesSurfaceReturnFlux[s]=SampleSpeciesSurfaceReturnFlux[0]+offsetSpecie;
      SampleSpeciesSurfaceInjectionFlux[s]=SampleSpeciesSurfaceInjectionFlux[0]+offsetSpecie;

      SampleReturnFluxBulkSpeed[s]=SampleReturnFluxBulkSpeed[0]+offsetSpecie;
      SampleInjectedFluxBulkSpeed[s]=SampleInjectedFluxBulkSpeed[0]+offsetSpecie;

      offsetSpecie+=TotalSurfaceElementNumber;

      for (el=0;el<TotalSurfaceElementNumber;el++) {
        SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
        SampleSpeciesSurfaceInjectionFlux[s][el]=0.0;

        SampleReturnFluxBulkSpeed[s][el]=0.0;
        SampleInjectedFluxBulkSpeed[s][el]=0.0;
      }
    }
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
