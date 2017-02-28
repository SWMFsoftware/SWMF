
//$Id$
//sampling of the ejection flux at the nucleus

#include "pic.h"
#include "SurfaceDataCG.h"


void cSurfaceDataCG::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"SourceRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"MassSourceRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"Probability Density[s=%i]\"",spec);

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    fprintf(fout,", \"Nude Gauge Density Contribution[s=%i]\",  \"Nude Gauge Flux Contribution[s=%i]\", \"Ram Gauge Density Contribution[s=%i]\",  \"Ram Gauge Flux Contribution[s=%i]\"",spec,spec,spec,spec);
  }
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
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    InjectionFlux[spec]=0.0,MassInjectionFlux[spec]=0.0;

    NudeGaugeDensityContribution[spec]=0.0,NudeGaugeFluxContribution[spec]=0.0;
    RamGaugeDensityContribution[spec]=0.0,RamGaugeFluxContribution[spec]=0.0;
  }

}

void cSurfaceDataCG::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataCG** InterpolationFaceList,int *Stencil,int StencilLength) {
  int i,spec,nface;
  double t;

  //output the calculated source rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->InjectionFlux[spec];
    if (PIC::LastSampleLength!=0) t/=PIC::LastSampleLength;

    fprintf(fout," %e",t);
  }

  //output the calculated mass source rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->MassInjectionFlux[spec];
    if (PIC::LastSampleLength!=0) t/=PIC::LastSampleLength;

    fprintf(fout," %e",t);
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

  //output contribution of the oarticular surface elements to the ram nad nude gauges measuremetns
  // \"Nude Gauge Density Contribution[s=%i]\",  \"Nude Gauge Pressure Contribution[s=%i]\", \"Ram Gauge Density Contribution[s=%i]\",  \"Ram Gauge Pressure Contribution[s=%i]
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->NudeGaugeDensityContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->NudeGaugeFluxContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->RamGaugeDensityContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->RamGaugeFluxContribution[spec];
    fprintf(fout," %e",t);

  }


}
