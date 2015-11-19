
#ifndef _CG_SURFACE_DATA_
#define _CG_SURFACE_DATA_

#include "mpichannel.h"

class cSurfaceDataCG {
public:
  double InjectionFlux[_TOTAL_SPECIES_NUMBER_];

  void PrintVarableList(FILE* fout);  
  void Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int StencilLength); 
  void Flush();  
  void Gather(CMPI_channel* pipe);
}; 

#endif
