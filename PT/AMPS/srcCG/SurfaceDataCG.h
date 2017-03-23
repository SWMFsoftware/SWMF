
#ifndef _CG_SURFACE_DATA_
#define _CG_SURFACE_DATA_

#include "mpichannel.h"

class cSurfaceDataCG {
public:
  double InjectionFlux[_TOTAL_SPECIES_NUMBER_];
  double MassInjectionFlux[_TOTAL_SPECIES_NUMBER_];

  //contribution to the Rosina nude nad ram gauges
  double NudeGaugeDensityContribution[_TOTAL_SPECIES_NUMBER_],NudeGaugeFluxContribution[_TOTAL_SPECIES_NUMBER_];
  double RamGaugeDensityContribution[_TOTAL_SPECIES_NUMBER_],RamGaugeFluxContribution[_TOTAL_SPECIES_NUMBER_];

  //save and output of the original and modified sources rates
  double OriginalSourceRate[_TOTAL_SPECIES_NUMBER_],ModifiedSourceRate[_TOTAL_SPECIES_NUMBER_];

  //save and output the scalar product of the external normal to the surface element and the vector pointing to the spacecraft location
  double CrossProduct_FaceNormal_SpacecraftLocation;

  //the marker saying that the surface element can be seen from the spacecraft location by the ram/nude gauge
  bool FieldOfView_NudeGauge,FieldOfView_RamGauge;

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int *Stencil,int StencilLength);
  void Flush();  
  void Gather(CMPI_channel* pipe);
}; 

#endif
