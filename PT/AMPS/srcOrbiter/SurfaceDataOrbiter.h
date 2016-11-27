
#ifndef _ORBITER_SURFACE_DATA_
#define _ORBITER_SURFACE_DATA_

#include "global.h"
#include "specfunc.h"
#include "mpichannel.h"

class cSurfaceDataOrbiter {
public:
  double EnergyTransferRate[_TOTAL_SPECIES_NUMBER_];
  double MomentumTransferRateX[_TOTAL_SPECIES_NUMBER_],MomentumTransferRateY[_TOTAL_SPECIES_NUMBER_],MomentumTransferRateZ[_TOTAL_SPECIES_NUMBER_];

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength);
  void Flush();  
  void Gather(CMPI_channel* pipe);
}; 

#endif
