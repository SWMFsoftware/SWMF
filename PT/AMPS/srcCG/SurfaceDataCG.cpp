
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
    fprintf(fout,", \"Nude Gauge Density Contribution[s=%i]\",  \"Nude Gauge Flux Contribution[s=%i]\", \"Ram Gauge Density Contribution[s=%i]\",  \"Ram Gauge Flux Contribution[s=%i]\", \"Original SourceRate[s=%i]\", \"Modified SourceRate[s=%i]\"",spec,spec,spec,spec,spec,spec);
  }

  fprintf(fout,", \"Scalar Product of the surface external notmal to the vector pf the spacecraft location\"");
  fprintf(fout,", \"Ram Gauge Field of View\", \"Nude Gauge Field of View\"");
  fprintf(fout,", \"ExternalNormX\", \"ExternalNormY\", \"ExternalNormZ\"");
  fprintf(fout,", \"Face ID\"");
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

        pipe->recv(t,thread);
        NudeGaugeDensityContribution[spec]+=t;

        pipe->recv(t,thread);
        NudeGaugeFluxContribution[spec]+=t;

        pipe->recv(t,thread);
        RamGaugeDensityContribution[spec]+=t;

        pipe->recv(t,thread);
        RamGaugeFluxContribution[spec]+=t;
      }
    }
    else {
      pipe->send(InjectionFlux[spec]);
      pipe->send(MassInjectionFlux[spec]);

      pipe->send(NudeGaugeDensityContribution[spec]);
      pipe->send(NudeGaugeFluxContribution[spec]);
      pipe->send(RamGaugeDensityContribution[spec]);
      pipe->send(RamGaugeFluxContribution[spec]);
    }
  }
}

void cSurfaceDataCG::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    InjectionFlux[spec]=0.0,MassInjectionFlux[spec]=0.0;

    NudeGaugeDensityContribution[spec]=0.0,NudeGaugeFluxContribution[spec]=0.0;
    RamGaugeDensityContribution[spec]=0.0,RamGaugeFluxContribution[spec]=0.0;

    OriginalSourceRate[spec]=0.0,ModifiedSourceRate[spec]=0.0;
  }

  ScalarProduct_FaceNormal_SpacecraftLocation=0.0;
  FieldOfView_NudeGauge=false,FieldOfView_RamGauge=false;
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
  // \"Nude Gauge Density Contribution[s=%i]\",  \"Nude Gauge Pressure Contribution[s=%i]\", \"Ram Gauge Density Contribution[s=%i]\",  \"Ram Gauge Pressure Contribution[s=%i], \"Original SourceRate[s=%i]\", \"Modified SourceRate[s=%i]\"
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->NudeGaugeDensityContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->NudeGaugeFluxContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->RamGaugeDensityContribution[spec];
    fprintf(fout," %e",t);

    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->RamGaugeFluxContribution[spec];
    fprintf(fout," %e",t);

    //the original source rate
    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->OriginalSourceRate[spec];
    fprintf(fout," %e",t);

    //the modified source rate
    for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->ModifiedSourceRate[spec];
    fprintf(fout," %e",t);
  }

  //the scalar product of the local external normal and the vector of the direction toward the spacecraft
  for (i=0,t=0.0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->ScalarProduct_FaceNormal_SpacecraftLocation;
  fprintf(fout," %e",t);

  //ram and hude gauges fields of view
  double RamGaugeFieldOfVewFlag=-1.0,NudeGaugeFieldOfVewFlag=-1.0;
  double ExternalNormX=0.0,ExternalNormY=0.0,ExternalNormZ=0.0,*l;

  for (i=0;i<StencilLength;i++) {
    if (InterpolationFaceList[i]->FieldOfView_NudeGauge>0.0) NudeGaugeFieldOfVewFlag=1.0;
    if (InterpolationFaceList[i]->FieldOfView_RamGauge>0.0) RamGaugeFieldOfVewFlag=1.0;

    l=CutCell::BoundaryTriangleFaces[Stencil[i]].ExternalNormal;

    ExternalNormX+=InterpolationWeightList[i]*l[0];
    ExternalNormY+=InterpolationWeightList[i]*l[1];
    ExternalNormZ+=InterpolationWeightList[i]*l[2];

  }

  fprintf(fout," %e  %e",RamGaugeFieldOfVewFlag,NudeGaugeFieldOfVewFlag);
  fprintf(fout," %e  %e  %e %i",ExternalNormX,ExternalNormY,ExternalNormZ,Stencil[0]);
}
