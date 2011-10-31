//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal Sphere boundary


#ifndef _AMR_INTERNAL_SURFACE_1D_SPHERE_
#define _AMR_INTERNAL_SURFACE_1D_SPHERE_

#include "math.h"


#include "meshAMRdef.h"
#include "mpichannel.h"


class cInternalSphere1DData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalSphericalData_UserDefined
#endif
{
protected:
  double OriginPosition[3],Radius;


public:

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
  fPrintDataStateVector PrintDataStateVector;

  typedef double (*fLocalResolution)(double *);
  fLocalResolution localResolution;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  void cleanDataBuffer() {
    for (int idim=0;idim<3;idim++) OriginPosition[idim]=0.0;
    Radius=0.0;

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL,localResolution=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalSphere1DData ()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
  : cInternalSphericalData_UserDefined()
#endif
  {
    cleanDataBuffer();
  }


  void SetSphereGeometricalParameters(double r) {
     for (int idim=0;idim<1;idim++) OriginPosition[idim]=0.0;
     Radius=r;
  }

  void GetSphereGeometricalParameters(double* &x0,double &r) {
     x0=OriginPosition;
     r=Radius;
  }

  void GetSurfaceElementIndex() {}
  long int GetTotalSurfaceElementsNumber() {return 1;}
  double GetSurfaceElementArea(int nPolarElement) { return 4.0*Pi*pow(Radius,2);}

  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
    FILE *fout=NULL;

    CMPI_channel pipe(1000000);
    int ThisThread=0,nTotalThreads=1;
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
    MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);

    if (ThisThread==0) {
      fout=fopen(fname,"w");
      pipe.openRecvAll();

      //print the output file title
      if (PrintTitle!=NULL) {
        PrintTitle(fout);
        fprintf(fout,"\n");
      }

      //print the variable list
      fprintf(fout,"VARIABLES=");
      if (PrintStateVectorFlag==true) {
        if (PrintVariableList==NULL) exit(__LINE__,__FILE__,"Error: PrintVariableList is not defined");
        PrintVariableList(fout);
      }

      fprintf(fout,"\n");
    }
    else pipe.openSend(0);

    if (PrintStateVectorFlag==true) {
      if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

      if (PrintDataStateVector!=NULL) PrintDataStateVector(fout,this,nDataSet,&pipe,ThisThread,nTotalThreads);

      if (ThisThread==0) fprintf(fout,"\n");
    }

    //close the pipe
    if (ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();

    if (ThisThread==0) fclose(fout);
  }


  void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

  //=======================================================================
  //intersection of a block with the Sphere
  int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
    if (xBlockMax[0]<Radius) return _AMR_BLOCK_OUTSIDE_DOMAIN_;
    else if (xBlockMin[0]>Radius) return _AMR_BLOCK_INSIDE_DOMAIN_;
    else return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
  }


  double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
    int BlockIntersectionCode;
    double res=0.0;

    BlockIntersectionCode=BlockIntersection(xBlockMinInit,xBlockMaxInit,EPS);
    *IntersectionStatus=BlockIntersectionCode;

    switch(BlockIntersectionCode) {
    case _AMR_BLOCK_OUTSIDE_DOMAIN_ :
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      res=1.0E-10*EPS;
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      res=pow(EPS,3);
#else
      exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

      break;
    case _AMR_BLOCK_INSIDE_DOMAIN_:
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      res=xBlockMaxInit[0]-xBlockMinInit[0];
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      res=4.0/3.0*Pi*(pow(xBlockMaxInit[0],3)-pow(xBlockMinInit[0],3));
#else
      exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

      break;
    case _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_:
#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
      res=xBlockMaxInit[0]-Radius;
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_
      res=4.0/3.0*Pi*(pow(xBlockMaxInit[0],3)-pow(Radius,3));
#else
      exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not found");
    }

    return res;
  }


};

#endif
