
//$Id$
//sampling of the ejection flux at the nucleus

#include "pic.h"
#include "SurfaceDataCG.h"


void cSurfaceDataCG::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"SourceRate[s=%i]\"",spec);
}

void cSurfaceDataCG::Gather(CMPI_channel* pipe) {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (PIC::ThisThread==0) {
      for (int thread=1;thread<PIC::nTotalThreads;thread++){
        double t;

        pipe->recv(t,thread);
        InjectionFlux[spec]+=t;
      }
    }
    else pipe->send(InjectionFlux[spec]);
  }
}

void cSurfaceDataCG::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) InjectionFlux[spec]=0.0;
}

void cSurfaceDataCG::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int StencilLength) {
  int i,spec;
  double t;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->InjectionFlux[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }
}
