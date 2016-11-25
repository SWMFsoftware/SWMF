
//$Id$
//sampling of the ejection flux at the nucleus

#include "SurfaceDataOrbiter.h"

#include "pic.h"


void cSurfaceDataOrbiter::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"EnergyTransferRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"absMomentumTransferRate[s=%i]\", \"MomentumTransferRateX[s=%i]\", \"MomentumTransferRateY[s=%i]\", \"MomentumTransferRateZ[s=%i]\"",spec,spec,spec,spec);
}

void cSurfaceDataOrbiter::Gather(CMPI_channel* pipe) {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (PIC::ThisThread==0) {
      for (int thread=1;thread<PIC::nTotalThreads;thread++){
        double t;

        pipe->recv(t,thread);
        EnergyTransferRate[spec]+=t;

        pipe->recv(t,thread);
        MomentumTransferRateX[spec]+=t;

        pipe->recv(t,thread);
        MomentumTransferRateY[spec]+=t;

        pipe->recv(t,thread);
        MomentumTransferRateZ[spec]+=t;
      }
    }
    else {
      pipe->send(EnergyTransferRate[spec]);
      pipe->send(MomentumTransferRateX[spec]);
      pipe->send(MomentumTransferRateY[spec]);
      pipe->send(MomentumTransferRateZ[spec]);
    }
  }
}

void cSurfaceDataOrbiter::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    EnergyTransferRate[spec]=0.0;
    MomentumTransferRateX[spec]=0.0,MomentumTransferRateY[spec]=0.0,MomentumTransferRateZ[spec]=0.0;
  }
}

void cSurfaceDataOrbiter::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength) {
  int i,spec,nface;
  double t,tX,tY,tZ;

  //output the calculated Energy Transfer rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->EnergyTransferRate[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);
  }

  //output the calculated Momentum Transfer rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    tX=0.0,tY=0.0,tZ=0.0;

    for (i=0;i<StencilLength;i++) {
      tX+=InterpolationWeightList[i]*InterpolationFaceList[i]->MomentumTransferRateX[spec];
      tY+=InterpolationWeightList[i]*InterpolationFaceList[i]->MomentumTransferRateY[spec];
      tZ+=InterpolationWeightList[i]*InterpolationFaceList[i]->MomentumTransferRateZ[spec];
    }

    tX/=PIC::LastSampleLength;
    tY/=PIC::LastSampleLength;
    tZ/=PIC::LastSampleLength;

    fprintf(fout," %e %e %e %e",sqrt(tX*tX+tY*tY+tZ*tZ),tX,tY,tZ);
  }
}
