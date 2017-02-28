
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

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int *Stencil,int StencilLength);
  void Flush();  
  void Gather(CMPI_channel* pipe);
}; 

#endif
