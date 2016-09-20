
//$Id$
//sampling of the ejection flux at the nucleus

#include "pic.h"
#include "SurfaceDataCG.h"


void cSurfaceDataCG::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"SourceRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"MassSourceRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"Probability Density[s=%i]\"",spec);
}

void cSurfaceDataCG::Gather(CMPI_channel* pipe) {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (PIC::ThisThread==0) {
      for (int thread=1;thread<PIC::nTotalThreads;thread++){
        double t;

        pipe->recv(t,thread);
        InjectionFlux[spec]+=t;

        pipe->recv(t,thread);
        MassInjectionFlux[spec]+=t;
      }
    }
    else {
      pipe->send(InjectionFlux[spec]);
      pipe->send(MassInjectionFlux[spec]);
    }
  }
}

void cSurfaceDataCG::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) InjectionFlux[spec]=0.0,MassInjectionFlux[spec]=0.0;
}

void cSurfaceDataCG::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int *Stencil,int StencilLength) {
  int i,spec,nface;
  double t;

  //output the calculated source rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->InjectionFlux[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }

  //output the calculated mass source rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->MassInjectionFlux[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }

  //output the theoretical value proportianal to the probability density of particle be generated at a particular location
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

/*
    if (Comet::BjornNASTRAN::SurfaceInjectionProbability[spec].ProbabilityDensityTable!=NULL) {
      for (i=0;i<StencilLength;i++) {
        nface=Stencil[i];
        t+=InterpolationWeightList[i]*Comet::BjornNASTRAN::SurfaceInjectionProbability[spec].ProbabilityDensityTable[nface]/CutCell::BoundaryTriangleFaces[nface].SurfaceArea;
      }
    }
*/

    fprintf(fout," %e",t);
  }
}
