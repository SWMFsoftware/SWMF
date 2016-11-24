
//$Id$
//sampling of the ejection flux at the nucleus

#include "SurfaceDataOrbiter.h"

#include "pic.h"


void cSurfaceDataOrbiter::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"EnergyTransferRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"MomentumTransferRate[s=%i]\"",spec);
}

void cSurfaceDataOrbiter::Gather(CMPI_channel* pipe) {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (PIC::ThisThread==0) {
      for (int thread=1;thread<PIC::nTotalThreads;thread++){
        double t;

        pipe->recv(t,thread);
        EnergyTransferRate[spec]+=t;

        pipe->recv(t,thread);
        MomentumTransferRate[spec]+=t;
      }
    }
    else {
      pipe->send(EnergyTransferRate[spec]);
      pipe->send(MomentumTransferRate[spec]);
    }
  }
}

void cSurfaceDataOrbiter::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) EnergyTransferRate[spec]=0.0,MomentumTransferRate[spec]=0.0;
}

void cSurfaceDataOrbiter::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength) {
  int i,spec,nface;
  double t;

  //output the calculated Energy Transfer rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->EnergyTransferRate[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }

  //output the calculated Momentum Transfer rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->MomentumTransferRate[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }
}
