//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal surfaces installed into the mesh


#ifndef _AMR_INTERNAL_SURFACE_
#define _AMR_INTERNAL_SURFACE_

#include "math.h"


#include "meshAMRdef.h"

//=======================================================================
//the descriptor of the internal boundary conditions
class cInternalBoundaryConditionsDescriptor {
public:
  unsigned char BondaryType;
  void *BoundaryElementProcessor;
  cInternalBoundaryConditionsDescriptor *nextInternalBCelement;

  void cleanDataBuffer() {
    BondaryType=_INTERNAL_BOUNDARY_TYPE_UNDEFINED_;
    BoundaryElementProcessor=NULL,nextInternalBCelement=NULL;
  }

  cInternalBoundaryConditionsDescriptor() {
    cleanDataBuffer();
  }
};


class cInternalSphericalData : public cAMRexit {
protected:
  double OriginPosition[3],Radius;
  void *SurfaceData;

  static long int nZenithSurfaceElements,nAzimuthalSurfaceElements;
  static double dZenithAngle,dAzimuthalAngle;

public:

  typedef void (*fPrintVarianleList)(FILE*);
  fPrintVarianleList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength);
  fPrintDataStateVector PrintDataStateVector;


  cInternalSphericalData () {
    OriginPosition[0]=0.0,OriginPosition[1]=0.0,OriginPosition[2]=0.0;
    Radius=0.0,SurfaceData=NULL;

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL;
  }


  void SetGeneralSurfaceMeshParameters(long int nZenithElements,long int nAzimuthalElements) {
    nZenithSurfaceElements=nZenithElements,nAzimuthalSurfaceElements=nAzimuthalElements;

  }

  void SetSurfaceDataPointer(void *DataPointer) {SurfaceData=DataPointer;}

  void SetSphereGeometricalParameters(double *x0,double r) {
     for (int idim=0;idim<3;idim++) OriginPosition[idim]=x0[idim];
     Radius=r;
  }

  long int GetLocalSurfaceElementNumber(long int nZenithElement,long int nAzimuthalElement) {
    return nZenithElement+nZenithSurfaceElements*nAzimuthalElement;
  }

  long int GetTotalSurfaceElementsNumber() {return nZenithSurfaceElements*nAzimuthalSurfaceElements;}

  void GetSurfaceElementProjectionIndex(double *x,long int &nZenithElement,long int &nAzimuthalElement) {
    double r,r2,ZenithAngle,AzimuthalAngle,xNormalized[3];
    int idim;

    for (r2=0.0,idim=0;idim<3;idim++) r2+=pow(x[idim],2);
    for (r=sqrt(r2),idim=0;idim<3;idim++) xNormalized[idim]=x[idim]/r;

    ZenithAngle=acos(xNormalized[2]);

    if ((r=pow(xNormalized[0],2)+pow(xNormalized[1],2))>0.0) {
      AzimuthalAngle=acos(xNormalized[0]/sqrt(r));
      if (xNormalized[1]<0.0) AzimuthalAngle=2.0*Pi-ZenithAngle;
    }
    else AzimuthalAngle=0.0;

    nZenithElement=(long int)(ZenithAngle/dZenithAngle);
    nAzimuthalElement=(long int)(AzimuthalAngle/dAzimuthalAngle);
  }

  void GetSurfaceCoordinate(double *x,long int iZenithPoint,long int iAzimutalPoint) {
    double ZenithAngle,AzimuthalAngle;

    ZenithAngle=dZenithAngle*iZenithPoint;
    AzimuthalAngle=dAzimuthalAngle*iAzimutalPoint;

    x[0]=Radius*sin(ZenithAngle)*cos(AzimuthalAngle)+OriginPosition[0];
    x[1]=Radius*sin(ZenithAngle)*sin(AzimuthalAngle)+OriginPosition[1];
    x[2]=Radius*cos(ZenithAngle)+OriginPosition[2];
  }

  void PrintSurfaceData(const char *fname,bool PrintStateVectorFlag=false) {
    long int iZenith,iAzimuthal;
    FILE *fout;
    double x[3];

    fout=fopen(fname,"w");

    //print the output file title
    if (PrintTitle!=NULL) {
      PrintTitle(fout);
      fprintf(fout,"\n");
    }

    //print the variable list
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");
    if (PrintStateVectorFlag==true) {
      if (PrintVariableList==NULL) exit(__LINE__,__FILE__,"Error: PrintVariableList is not defined");
      PrintVariableList(fout);
    }

    //print the number of variables and blocks
    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",2+(nZenithSurfaceElements-1)*nAzimuthalSurfaceElements,nZenithSurfaceElements*nAzimuthalSurfaceElements);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      GetSurfaceCoordinate(x,iZenith,iAzimuthal);
      fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");
        PrintDataStateVector(fout,iZenith,iAzimuthal,InterpolationList,InterpolationListLength);
      }

      fprintf(fout,"\n");

      //only one point is printed for azimuthal angle of iAzimuthal==0, nAzimuthalSurfaceElements
      if ((iZenith==0)||(iZenith==nZenithSurfaceElements)) break;
    }

    //print the connectivity list
    long int iAzimuthalMax,iAzimuthalMin;
    long int nd0,nd1,nd2,nd3;

    for (iZenith=0;iZenith<nZenithSurfaceElements;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
      iAzimuthalMin=iAzimuthal;

      nd0=(iZenith!=0) ? 2+iAzimuthalMin+(iZenith-1)*nAzimuthalSurfaceElements : 1; //iZenithMin,iAzimuthalMin
      nd1=(iZenith!=nZenithSurfaceElements-1) ? 2+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements : 2+(nZenithSurfaceElements-1)*nAzimuthalSurfaceElements; //iZenithMax,iAzimuthalMin
      nd2=(iZenith!=nZenithSurfaceElements-1) ? 2+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements : 2+(nZenithSurfaceElements-1)*nAzimuthalSurfaceElements; //iZenithMax,iAzimuthalMax
      nd3=(iZenith!=0) ? 2+iAzimuthalMax+(iZenith-1)*nAzimuthalSurfaceElements : 1; //iZenithMin,iAzimuthalMax

      fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
    }

    fclose(fout);
  }


   void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,false);}

};


#endif
