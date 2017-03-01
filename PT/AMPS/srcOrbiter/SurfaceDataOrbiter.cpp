
//$Id$
//sampling of the ejection flux at the nucleus

#include "SurfaceDataOrbiter.h"

#include "pic.h"


void cSurfaceDataOrbiter::PrintVarableList(FILE* fout) {
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"EnergyTransferRate[s=%i]\"",spec);
  for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout,", \"absMomentumTransferRate[s=%i]\", \"MomentumTransferRateX[s=%i]\", \"MomentumTransferRateY[s=%i]\", \"MomentumTransferRateZ[s=%i]\"",spec,spec,spec,spec);

  //surface aboundance and adsorption/desorption rates
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    fprintf(fout,", \"Surface Aboundance [m^-2][s=%i]\" ",spec);
  }

  //source rate of the injected species
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    fprintf(fout,", \"SourceRate [m^-2 s^-1][s=%i]\" ",spec);
  }
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

        pipe->recv(t,thread);
        SourceRate[spec]+=t;
      }
    }
    else {
      pipe->send(EnergyTransferRate[spec]);
      pipe->send(MomentumTransferRateX[spec]);
      pipe->send(MomentumTransferRateY[spec]);
      pipe->send(MomentumTransferRateZ[spec]);
      pipe->send(SourceRate[spec]);
    }
  }
}

void cSurfaceDataOrbiter::Flush() {
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    EnergyTransferRate[spec]=0.0;
    MomentumTransferRateX[spec]=0.0,MomentumTransferRateY[spec]=0.0,MomentumTransferRateZ[spec]=0.0;
    AdsorptionFlux[spec]=0.0,DesorptionFlux[spec]=0.0,SourceRate[spec]=0.0;
  }
}

void cSurfaceDataOrbiter::Print(FILE *fout,double* InterpolationWeightList,cSurfaceDataOrbiter** InterpolationFaceList,int *Stencil,int StencilLength) {
  int i,spec,nface;
  double t,tX,tY,tZ,tRate;

  //output the calculated Energy Transfer rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    t=0.0;

    for (i=0;i<StencilLength;i++) t+=InterpolationWeightList[i]*InterpolationFaceList[i]->EnergyTransferRate[spec];
    fprintf(fout," %e",t/PIC::LastSampleLength);

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(t,__LINE__,__FILE__);
    #endif
    #endif
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

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(tX,__LINE__,__FILE__);
    PIC::Debugger::CatchOutLimitValue(tY,__LINE__,__FILE__);
    PIC::Debugger::CatchOutLimitValue(tZ,__LINE__,__FILE__);
    #endif
    #endif
  }


  //surface aboundance and adsorption/desorption rates
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double tAboundance=0.0,TotalStencilArea=0.0;

    //average the surface aboundance weighting not with the surface area but.....
    double AboundanceInterpolationWeight[StencilLength],summAboundanceInterpolationWeight=0.0;

    for (i=0;i<StencilLength;i++) if (CutCell::BoundaryTriangleFaces[Stencil[i]].SurfaceArea>0.0) {
      double w=InterpolationFaceList[i]->SpeciesSurfaceAboundance[spec]/CutCell::BoundaryTriangleFaces[Stencil[i]].SurfaceArea;

      tAboundance+=w*InterpolationFaceList[i]->SpeciesSurfaceAboundance[spec]/CutCell::BoundaryTriangleFaces[Stencil[i]].SurfaceArea;
      summAboundanceInterpolationWeight+=w;
    }

    //if (TotalStencilArea>0.0) tAboundance/=TotalStencilArea;

    if (summAboundanceInterpolationWeight>0.0) tAboundance/=summAboundanceInterpolationWeight;


/*    for (i=0;i<StencilLength;i++) {
      tAboundance+=InterpolationFaceList[i]->SpeciesSurfaceAboundance[spec];
      TotalStencilArea+=CutCell::BoundaryTriangleFaces[Stencil[i]].SurfaceArea;
    }*/

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(tAboundance,__LINE__,__FILE__);
 //   PIC::Debugger::CatchOutLimitValue(TotalStencilArea,__LINE__,__FILE__);
    #endif
    #endif


    fprintf(fout," %e ",tAboundance);
  }

  //the source rate
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    tRate=0.0;

    for (i=0;i<StencilLength;i++) {
      tRate+=InterpolationWeightList[i]*InterpolationFaceList[i]->SourceRate[spec];
    }

    tRate/=PIC::LastSampleLength;
    fprintf(fout," %e ",tRate);

    #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    #if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    PIC::Debugger::CatchOutLimitValue(tRate,__LINE__,__FILE__);
    #endif
    #endif
  }
}
