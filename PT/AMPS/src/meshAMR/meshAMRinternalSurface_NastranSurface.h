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

#include "meshAMRcutcell.h"


class cInternalNastranSurfaceData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_ON_
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
#if _USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_NASTRAN_SURFACE_MODE_ON_
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
    int nt;

    //check possible intersection between the block and the surface
    for (nt=0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
      if (CutCell::BoundaryTriangleFaces[nt].BlockIntersection(xBlockMin,xBlockMax,EPS)==true) {
        return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
      }
    }

    //if no intersection is found, determin if the miggle point of the block is within the surface
    double x[3];
    int res;

    for (int i=0;i<3;i++) x[i]=0.5*(xBlockMin[i]+xBlockMax[i]);

    res=(CutCell::CheckPointInsideDomain(x,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,EPS)==true) ? _AMR_BLOCK_INSIDE_DOMAIN_ : _AMR_BLOCK_OUTSIDE_DOMAIN_;

    return res;
  }


  double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
    double res=1.0;

    for (int idim=0;idim<3;idim++) res*=fabs(xBlockMaxInit[idim]-xBlockMinInit[idim]);

    return res;
  }


};

#endif
