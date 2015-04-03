//$Id$

#ifndef _OH_H_
#define _OH_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

// charge exchange physical model
#include "ChargeExchange.h"

// SPICE is not used for this application
// use empty functions to fill in
#include "SpiceEmptyDefinitions.h"

namespace OH {
  using namespace Exosphere;
  
  //  injection boundary condition
  extern double InjectionVelocity[3];
  extern double InjectionNDensity;
  extern double InjectionTemperature;
  
  // computational domain size
  extern double DomainXMin[3];
  extern double DomainXMax[3];
  extern double DomainDXMin;
  extern double DomainDXMax;
  

  void Init_BeforeParser();

  //---------------------------------------------------------------------------
  namespace Sampling{
    using namespace Exosphere::Sampling;
  }

  namespace Coupling {
    extern double TimeAfterCoupling[PIC::nTotalSpecies];
    void Send(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI);
  }

  //---------------------------------------------------------------------------
  namespace Output{

    extern int TotalDataLength;
    extern int ohSourceDensityOffset; 
    extern int ohSourceMomentumOffset;
    extern int ohSourceEnergyOffset;
    
    void Init();

    void PrintVariableList(FILE* fout,int DataSetNumber);

    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

    int RequestDataBuffer(int offset);
  }

  //---------------------------------------------------------------------------
  namespace Loss {
    double LifeTime(double *x, int spec, long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    int ReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
  }

  //---------------------------------------------------------------------------
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {

    accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

  }
  

}

#endif
