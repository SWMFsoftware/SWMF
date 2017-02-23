
#ifndef _ORBITER_SURFACE_DATA_
#define _ORBITER_SURFACE_DATA_

#include "global.h"
#include "specfunc.h"
#include "mpichannel.h"

class cSurfaceDataOrbiter {
public:
  double EnergyTransferRate[_TOTAL_SPECIES_NUMBER_];
  double MomentumTransferRateX[_TOTAL_SPECIES_NUMBER_],MomentumTransferRateY[_TOTAL_SPECIES_NUMBER_],MomentumTransferRateZ[_TOTAL_SPECIES_NUMBER_];

  //the surface aboundance of the model species (if needed for modeing of adsorption/desorption)
  double SpeciesSurfaceAboundance[_TOTAL_SPECIES_NUMBER_];

  //sampling of the adsorption/desorption rates
  double AdsorptionFlux[_TOTAL_SPECIES_NUMBER_],DesorptionFlux[_TOTAL_SPECIES_NUMBER_];

  //source rate of the simulated speces ejected from the surface
  double SourceRate[_TOTAL_SPECIES_NUMBER_];

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength);
  void Flush();  
  void Gather(CMPI_channel* pipe);

  //constructor set the defauls values of the sampling buffers
  cSurfaceDataOrbiter() {
    //init the surface aboundance buffer
    for (int spec=0;spec<_TOTAL_SPECIES_NUMBER_;spec++) SpeciesSurfaceAboundance[spec]=0.0;

    //init other sampling buffers
    Flush();
  }
}; 

#endif
