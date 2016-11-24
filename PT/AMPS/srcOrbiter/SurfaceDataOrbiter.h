
#ifndef _ORBITER_SURFACE_DATA_
#define _ORBITER_SURFACE_DATA_

#include "mpichannel.h"

class cSurfaceDataOrbiter {
public:
  double EnergyTransferRate[_TOTAL_SPECIES_NUMBER_];
  double MomentumTransferRate[_TOTAL_SPECIES_NUMBER_];

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength);
  void Flush();  
  void Gather(CMPI_channel* pipe);
}; 

#endif
