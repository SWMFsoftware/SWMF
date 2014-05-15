//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal Sphere boundary


#ifndef _AMR_INTERNAL_SURFACE_NASTRAN_SURFACE_
#define _AMR_INTERNAL_SURFACE_NASTRAN_SURFACE_

#include "math.h"


#include "meshAMRdef.h"
#include "mpichannel.h"


class cInternalNastranSurfaceData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalNastranSurfaceData_UserDefined
#endif
{
protected:



public:

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nElement,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalNastranSurfaceData *Sutface,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
  fPrintDataStateVector PrintDataStateVector;

  typedef double (*fLocalResolution)(double *);
  fLocalResolution localResolution;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  void cleanDataBuffer() {

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL,localResolution=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalNastranSurfaceData()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
  : cInternalNastranSurfaceData_UserDefined()
#endif
  {
    cleanDataBuffer();
  }




  void GetSurfaceElementIndex() {}
  long int GetTotalSurfaceElementsNumber() {return 1;}
  double GetSurfaceElementArea(int nPolarElement) { return 1.0;}

  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
  }


  void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

  //=======================================================================
  //intersection of a block with the Sphere
  int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
    return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
  }


  double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
    return 0.0;
  }


};

#endif
