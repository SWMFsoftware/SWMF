//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal surfaces installed into the mesh


#ifndef _AMR_INTERNAL_SURFACE__ROTATION_BODY_
#define _AMR_INTERNAL_SURFACE__ROTATION_BODY_

#include "math.h"

#include "meshAMRdef.h"
#include "meshAMRcutcell.h"

#include "mpichannel.h"
#include "rnd.h"


//=======================================================================
//the class describes the data that defines the spherical internal boundary; the class may contains the user defined data
class cInternalRotationBodyData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalRotationBodyData_UserDefined
#endif
{
//protected:
public:
  double OriginPosition[3],AxisOfSymmetry[3],xAxisMin,xAxisMax,rSurfaceMax;
  double e0[3],e1[3],e2[3]; //internal frame of reference
  double *SurfaceElementNormalVector;


  static int nAzimuthalSurfaceElements;
  double dAzimuthalAngle;

  static int nAxisSurfaceElements;
  double *xAxisSurfaceElement;

  typedef bool (*fSurfaceCurve)(double&,double);
  fSurfaceCurve SurfaceCurve;

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalRotationBodyData *Surface,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
  fPrintDataStateVector PrintDataStateVector;

  typedef double (*fLocalResolution)(double *);
  fLocalResolution localResolution;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif


  void cleanDataBuffer() {
    for (int i=0;i<3;i++) OriginPosition[i]=0.0,AxisOfSymmetry[i]=0.0;
    xAxisMin=0.0,xAxisMax=0.0,rSurfaceMax=0.0;

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL,localResolution=NULL,SurfaceCurve=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalRotationBodyData ()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
 : cInternalRotationBodyData_UserDefined()
#endif
  {

    SurfaceElementNormalVector=NULL;
//    nAxisSurfaceElements=0,nAzimuthalSurfaceElements=0;

    cleanDataBuffer();
//    SetGeneralSurfaceMeshParameters(nZenithSurfaceElements,nAzimuthalSurfaceElements);
  }

  static void SetGeneralSurfaceMeshParameters(long int nAxisElements,long int nAzimuthalElements) {
    nAxisSurfaceElements=nAxisElements,nAzimuthalSurfaceElements=nAzimuthalElements;
  }


  void SetGeometricalParameters(double *x0,double *l,double xmin,double xmax,fSurfaceCurve f) {
     int idim,iAxis;
     double c=0.0;

     if ((nAxisSurfaceElements==0)||(nAzimuthalSurfaceElements==0)) exit(__LINE__,__FILE__,"Error: nAxisSurfaceElements and nAzimuthalSurfaceElements has to be set up first");

     for (idim=0;idim<3;idim++) {
       OriginPosition[idim]=x0[idim],AxisOfSymmetry[idim]=l[idim],e0[idim]=l[idim];
       c+=pow(l[idim],2);
     }

     xAxisMin=xmin,xAxisMax=xmax;
     SurfaceCurve=f;

     dAzimuthalAngle=2.0*Pi/nAzimuthalSurfaceElements;
     xAxisSurfaceElement=new double [nAxisSurfaceElements+1];

     SurfaceElementNormalVector=new double [3*nAzimuthalSurfaceElements*nAxisSurfaceElements];

     //derive the local coordinate frame
     int j,imax=-1;
     double t=-1.0;

     for (idim=0,c=sqrt(c);idim<3;idim++) {
       e0[idim]/=c;
       e1[idim]=0.0;

       if (t<fabs(e0[idim])) t=fabs(e0[idim]),imax=idim;
     }

     j=imax+1;
     if (j==3) j=0;

     e1[j]=e0[imax];
     e1[imax]=-e0[j];

     for (idim=0,c=0.0;idim<3;idim++) c+=pow(e1[idim],2);
     for (idim=0,c=sqrt(c);idim<3;idim++) e1[idim]/=c;

     e2[0]=+(e0[1]*e1[2]-e0[2]*e1[1]);
     e2[1]=-(e0[0]*e1[2]-e0[2]*e1[0]);
     e2[2]=+(e0[0]*e1[1]-e0[1]*e1[0]);

     //Uniform distribution
     for (int i=0;i<nAxisSurfaceElements+1;i++) xAxisSurfaceElement[i]=xAxisMin+(xAxisMax-xAxisMin)/double(nAxisSurfaceElements)*i;


     //init the surface normals
     for (iAxis=0;iAxis<nAxisSurfaceElements;iAxis++) for (int iAzimuth=0;iAzimuth<nAzimuthalSurfaceElements;iAzimuth++) {
       int nSurfaceElementNumber=GetLocalSurfaceElementNumber(iAxis,iAzimuth);
       double x[3],norm[3];

       GetSurfaceElementMiddlePoint(x,iAxis,iAzimuth);
       GetSurfaceNormal(norm,x);

       memcpy(SurfaceElementNormalVector+3*nSurfaceElementNumber,norm,3*sizeof(double));
     }

     //calcualte the maximum value of the radius for the surface 'rSurfaceMax'
     for (iAxis=0,rSurfaceMax=-1.0;iAxis<nAxisSurfaceElements;iAxis++) {
       double r0;

       if (SurfaceCurve(r0,xAxisSurfaceElement[iAxis])==false) exit(__LINE__,__FILE__,"Error: the surface is not defined in the whole rangle (xmin,xmax)");
       if (r0>rSurfaceMax) rSurfaceMax=r0;
     }

  }

  void GetSphereGeometricalParameters(double* &x0,double* &l,double &xmin,double &xmax) {
     x0=OriginPosition,l=AxisOfSymmetry;
     xmin=xAxisMin,xmax=xAxisMax;
  }

  inline long int GetLocalSurfaceElementNumber(long int nAxisElement,long int nAzimuthalElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nAxisElement<0)||(nAxisElement>=nAxisSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    #endif

    return nAxisElement+nAxisSurfaceElements*nAzimuthalElement;
  }

  inline void GetSurfaceElementIndex(int &nAxisElement,int &nAzimuthalElement,int nSurfaceElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nSurfaceElement<0)||(nSurfaceElement>=nAxisSurfaceElements*nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nSurfaceElement' is out of range");
    #endif

    nAzimuthalElement=(int)(nSurfaceElement/nAxisSurfaceElements);
    nAxisElement=nSurfaceElement-nAxisSurfaceElements*nAzimuthalElement;
  }

  long int GetTotalSurfaceElementsNumber() {return nAxisSurfaceElements*nAzimuthalSurfaceElements;}


  double GetSurfaceElementArea(int nAxisElement,int nAzimuthalElement) {
    double theta,dx,r1,r0;

    dx=xAxisSurfaceElement[nAxisElement+1]-xAxisSurfaceElement[nAxisElement];
    if ((SurfaceCurve(r0,xAxisSurfaceElement[nAxisElement])==false)||(SurfaceCurve(r1,xAxisSurfaceElement[nAxisElement+1])==false)) exit(__LINE__,__FILE__,"Error: there is a problem somewhere :-(");
    theta=atan(fabs((r1-r0)/dx));

    return Pi/sin(theta)*fabs(r1*r1-r0*r0)/nAzimuthalSurfaceElements;
  }

  double GetSurfaceElementArea(int SurfaceElementNumber) {
    int nAxisElement,nAzimuthalElement;

    GetSurfaceElementIndex(nAxisElement,nAzimuthalElement,SurfaceElementNumber);
    return GetSurfaceElementArea(nAxisElement,nAzimuthalElement);
  } 


  inline void GetSurfaceElementNormal(double *norm,int nAxisElement,int nAzimuthalElement) {
    GetSurfaceElementNormal(norm,GetLocalSurfaceElementNumber(nAxisElement,nAzimuthalElement));
  }

  inline void GetSurfaceElementNormal(double *norm,int SurfaceElementNumber) {
    for (int i=0;i<3;i++) norm[i]=SurfaceElementNormalVector[i+3*SurfaceElementNumber];
  }



  inline void GetSurfaceElementMiddlePoint(double *x,int nAxisElement,int nAzimuthalElement) {
    double locx[3],r,AzimuthAngle;

    locx[0]=(xAxisSurfaceElement[nAxisElement]+xAxisSurfaceElement[nAxisElement+1])/2.0;
    AzimuthAngle=(nAzimuthalElement+0.5)*dAzimuthalAngle;

    if (SurfaceCurve(r,locx[0])==false) exit(__LINE__,__FILE__,"Error: there is a problem somewhere :-(");
    locx[1]=r*cos(AzimuthAngle);
    locx[2]=r*sin(AzimuthAngle);

    for (int i=0;i<3;i++) x[i]=OriginPosition[i]+locx[0]*e0[i]+locx[1]*e1[i]+locx[2]*e2[i];
  }


  inline void GetSurfaceElementMiddlePoint(double *x,int SurfaceElementNumber) {
    int nAxisElement,nAzimuthalElement;

    GetSurfaceElementIndex(nAxisElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementMiddlePoint(x,nAxisElement,nAzimuthalElement);
  }


  inline void GetSurfaceElementRandomPoint(double *x,int nAxisElement,int nAzimuthalElement) {
    double locx[3],AzimuthAngle,r;
    double alpha,r1,r0;

    SurfaceCurve(r1,xAxisSurfaceElement[nAxisElement+1]);
    SurfaceCurve(r0,xAxisSurfaceElement[nAxisElement]);

    alpha=(r1-r0)/(xAxisSurfaceElement[nAxisElement+1]-xAxisSurfaceElement[nAxisElement]);

    if (fabs(alpha)<1.0E-10) {
      //cyclinder
      locx[0]=xAxisSurfaceElement[nAxisElement]+rnd()*(xAxisSurfaceElement[nAxisElement+1]-xAxisSurfaceElement[nAxisElement]);
    }
    else {
      //cone
      double x0,x1,x0Cone;

      x0Cone= (r0>r1) ? (xAxisSurfaceElement[nAxisElement]-alpha/r0) : (xAxisSurfaceElement[nAxisElement+1]-alpha/r1);

      if (alpha>0.0) {
        x1=xAxisSurfaceElement[nAxisElement+1]-x0Cone;
        x0=xAxisSurfaceElement[nAxisElement]-x0Cone;
   
        locx[0]=sqrt(x0*x0+rnd()*(x1*x1-x0*x0))+x0Cone;
      }
      else {
        x0=x0Cone-xAxisSurfaceElement[nAxisElement+1];
        x1=x0Cone-xAxisSurfaceElement[nAxisElement];

        locx[0]=x0Cone-sqrt(x0*x0+rnd()*(x1*x1-x0*x0));
      }
    }


    SurfaceCurve(r,locx[0]);

    AzimuthAngle=2.0*Pi*rnd();
    locx[1]=r*cos(AzimuthAngle);
    locx[2]=r*sin(AzimuthAngle);

    for (int i=0;i<3;i++) x[i]=OriginPosition[i]+locx[0]*e0[i]+locx[1]*e1[i]+locx[2]*e2[i];
  }

  inline bool GetSurfaceElementProjectionIndex(double *x,long int &nAxisElement,long int &nAzimuthalElement) {
    double AzimuthAngle,locx[3]={0.0,0.0,0.0},r;
    int i;
    bool flag;

    for (i=0;i<3;i++) {
      locx[0]+=e0[i]*(x[i]-OriginPosition[i]);
      locx[1]+=e1[i]*(x[i]-OriginPosition[i]);
      locx[2]+=e2[i]*(x[i]-OriginPosition[i]);
    }

    for (i=0,flag=false;i<nAxisSurfaceElements;i++) if ((xAxisSurfaceElement[i]<=locx[0])&&(locx[0]<=xAxisSurfaceElement[i+1])) {
      nAxisElement=i;
      break;
    }

    if (flag==false) {
      nAxisElement=-1,nAzimuthalElement=-1;
      return false;
    }

    r=sqrt(locx[1]*locx[1]+locx[2]*locx[2]);
    AzimuthAngle=acos(locx[1]/r);
    if (locx[2]<0.0) AzimuthAngle=2.0*Pi-AzimuthAngle;

    nAzimuthalElement=(long int)(AzimuthAngle/dAzimuthalAngle);
    return true;
  }


  inline void GetSurfaceNormal(double *norm,double *x) {
    double AzimuthAngle,locx[3]={0.0,0.0,0.0},r;
    int i,iAxis;
    bool flag;

    //determine the local coordinates of the point
    for (i=0;i<3;i++) {
      locx[0]+=e0[i]*(x[i]-OriginPosition[i]);
      locx[1]+=e1[i]*(x[i]-OriginPosition[i]);
      locx[2]+=e2[i]*(x[i]-OriginPosition[i]);
    }

    for (i=0,flag=false;i<nAxisSurfaceElements;i++) if ((xAxisSurfaceElement[i]<=locx[0])&&(locx[0]<=xAxisSurfaceElement[i+1])) {
      iAxis=i;
      flag=true;
      break;
    }

    if (flag==false) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    r=sqrt(locx[1]*locx[1]+locx[2]*locx[2]);
    AzimuthAngle=acos(locx[1]/r);
    if (locx[2]<0.0) AzimuthAngle=2.0*Pi-AzimuthAngle;

    //calculate the normal
    double r0,r1,norm2d[2],norm3d[3];

    if ((SurfaceCurve(r0,xAxisSurfaceElement[iAxis])==false)||(SurfaceCurve(r1,xAxisSurfaceElement[iAxis+1])==false)) exit(__LINE__,__FILE__,"Error: there is a problem somewhere :-(");

    norm2d[0]=-(r1-r0),norm2d[1]=xAxisSurfaceElement[iAxis+1]-xAxisSurfaceElement[iAxis];
    r=sqrt(norm2d[0]*norm2d[0]+norm2d[1]*norm2d[1]);
    norm2d[0]/=r,norm2d[1]/=r;

    if (norm2d[1]<0.0) norm2d[0]*=-1.0,norm2d[1]*=-1.0;

    norm3d[0]=norm2d[0];
    norm3d[1]=norm2d[1]*cos(AzimuthAngle);
    norm3d[2]=norm2d[1]*sin(AzimuthAngle);

    norm[0]=norm3d[0]*e0[0]+norm3d[1]*e1[0]+norm3d[2]*e2[0];
    norm[1]=norm3d[0]*e0[1]+norm3d[1]*e1[1]+norm3d[2]*e2[1];
    norm[2]=norm3d[0]*e0[2]+norm3d[1]*e1[2]+norm3d[2]*e2[2];
  }

  inline void GetSurfaceCoordinate(double *x,int iAxis,int iAzimutalPoint) {
    double l[3],r;

    l[0]=xAxisSurfaceElement[iAxis];
    if (SurfaceCurve(r,l[0])==false) exit(__LINE__,__FILE__,"Error: there is a problem somewhere :-(");
    if (r<0.0) r=0.0;

    l[1]=r*cos(dAzimuthalAngle*iAzimutalPoint);
    l[2]=r*sin(dAzimuthalAngle*iAzimutalPoint);

    for (int i=0;i<3;i++) x[i]=OriginPosition[i]+l[0]*e0[i]+l[1]*e1[i]+l[2]*e2[i];
  }

  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
    long int iAxis,iAzimuthal;
    FILE *fout=NULL;
    double x[3];

    CMPI_channel pipe(1000000);
    int ThisThread=0,nTotalThreads=1;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

    if (ThisThread==0) {
      fout=fopen(fname,"w");
      pipe.openRecvAll();

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
      fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",(nAxisSurfaceElements+1)*nAzimuthalSurfaceElements,nAxisSurfaceElements*nAzimuthalSurfaceElements);
    }
    else pipe.openSend(0);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iAxis=0;iAxis<nAxisSurfaceElements+1;iAxis++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      GetSurfaceCoordinate(x,iAxis,iAzimuthal);
      if (ThisThread==0) fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iAxis==0) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else if (iAxis==nAxisSurfaceElements) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nAxisSurfaceElements-1,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nAxisSurfaceElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else {
          int iA,iZ,A[2],Z[2];

          Z[0]=iAxis-1,Z[1]=iAxis;

          A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1;
          A[1]=iAzimuthal;

          for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
        }


        PrintDataStateVector(fout,iAxis,iAzimuthal,InterpolationList,InterpolationListLength,this,nDataSet,&pipe,ThisThread,nTotalThreads);
      }

      if (ThisThread==0) fprintf(fout,"\n");
    }

    //close the pipe
    if (ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();

    //print the connectivity list
    long int iAzimuthalMax,iAzimuthalMin;
    long int nd0,nd1,nd2,nd3;

    if (ThisThread==0) {
      for (iAxis=0;iAxis<nAxisSurfaceElements;iAxis++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
        iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
        iAzimuthalMin=iAzimuthal;

        nd0=1+iAzimuthalMin+iAxis*nAzimuthalSurfaceElements;
        nd1=1+iAzimuthalMax+iAxis*nAzimuthalSurfaceElements;
        nd2=1+iAzimuthalMax+(iAxis+1)*nAzimuthalSurfaceElements;
        nd3=1+iAzimuthalMin+(iAxis+1)*nAzimuthalSurfaceElements;

        fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
      }

      fclose(fout);
    }

  }


   void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

   //=======================================================================
   //calculate the point of intersection between an edge and a segment of the surface
   bool LinearSegment_Surface_Intersection(double x0Surface_LocalFrame,double x1Surface_LocalFrame,double *x0LinearSegment_GlobalFrame,double *x1LinearSegment_GlobalFrame,double *IntersectionPoint_GlobalFrame) {
     double r0Surface,r1Surface;

     double x0SurfaceMin,x0SurfaceMax;


     x0SurfaceMin=min(x0Surface_LocalFrame,x1Surface_LocalFrame);
     x0SurfaceMax=max(x0Surface_LocalFrame,x1Surface_LocalFrame);

     if (x0SurfaceMin<xAxisMin) x0SurfaceMin=xAxisMin;
     if (x0SurfaceMax>xAxisMax) x0SurfaceMax=xAxisMax;

     x0Surface_LocalFrame=x0SurfaceMin;
     x1Surface_LocalFrame=x0SurfaceMax;


     if (SurfaceCurve(r0Surface,x0Surface_LocalFrame)==false) return false;
     if (SurfaceCurve(r1Surface,x1Surface_LocalFrame)==false) return false;

     double a,b,c,alphaCone,x0Cone,t1,t2;
     double l[3],d,sqrt_d,EPS;
     int idim;


     for (idim=0;idim<3;idim++) l[idim]=x1LinearSegment_GlobalFrame[idim]-x0LinearSegment_GlobalFrame[idim];

     EPS=1.0E-10*sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
     alphaCone=(r1Surface-r0Surface)/(x1Surface_LocalFrame-x0Surface_LocalFrame);


       //transfer the initial point and pointing vector of the efge into the frame related to the surface
     double lLocalFrame[3]={0.0,0.0,0.0},x0LocalFrame[3]={0.0,0.0,0.0};


     lLocalFrame[0]=l[0]*e0[0]+l[1]*e0[1]+l[2]*e0[2];
     lLocalFrame[1]=l[0]*e1[0]+l[1]*e1[1]+l[2]*e1[2];
     lLocalFrame[2]=l[0]*e2[0]+l[1]*e2[1]+l[2]*e2[2];

     x0LocalFrame[0]=(x0LinearSegment_GlobalFrame[0]-OriginPosition[0])*e0[0]+(x0LinearSegment_GlobalFrame[1]-OriginPosition[1])*e0[1]+(x0LinearSegment_GlobalFrame[2]-OriginPosition[2])*e0[2];
     x0LocalFrame[1]=(x0LinearSegment_GlobalFrame[0]-OriginPosition[0])*e1[0]+(x0LinearSegment_GlobalFrame[1]-OriginPosition[1])*e1[1]+(x0LinearSegment_GlobalFrame[2]-OriginPosition[2])*e1[2];
     x0LocalFrame[2]=(x0LinearSegment_GlobalFrame[0]-OriginPosition[0])*e2[0]+(x0LinearSegment_GlobalFrame[1]-OriginPosition[1])*e2[1]+(x0LinearSegment_GlobalFrame[2]-OriginPosition[2])*e2[2];

     if (fabs(alphaCone)>1.0E-20) {
       x0Cone=x0Surface_LocalFrame-r0Surface/alphaCone;

       //!!!!! Maple worksheet used to generate the code is in utility/meshAMRinternalSurface_RotationBody.mw

       double t1 = pow(lLocalFrame[2], 0.2e1);
       double t2 = alphaCone * alphaCone;
       double t3 = pow(lLocalFrame[0], 0.2e1);
       double t5 = pow(lLocalFrame[1], 0.2e1);
       double t11 = t2 * x0LocalFrame[0];
       double t14 = pow(x0LocalFrame[1], 0.2e1);
       double t15 = pow(x0LocalFrame[0], 0.2e1);
       double t17 = pow(x0LocalFrame[2], 0.2e1);
       double t20 = x0Cone * x0Cone;
       a = t1 - t2 * t3 + t5;
       b = 0.2e1 * x0LocalFrame[1] * lLocalFrame[1] + 0.2e1 * t2 * lLocalFrame[0] * x0Cone + 0.2e1 * x0LocalFrame[2] * lLocalFrame[2] - 0.2e1 * t11 * lLocalFrame[0];
       c = t14 - t2 * t15 + t17 + 0.2e1 * t11 * x0Cone - t2 * t20;
     }
     else {
       double t1 = pow(lLocalFrame[2], 0.2e1);
       double t2 = pow(lLocalFrame[1], 0.2e1);
       double t9 = pow(x0LocalFrame[1], 0.2e1);
       double t10 = r0Surface * r0Surface;
       double t12 = pow(x0LocalFrame[2], 0.2e1);
       double t13 = r1Surface * r1Surface;
       a = t1 + t2;
       b = 0.2e1 * x0LocalFrame[1] * lLocalFrame[1] + 0.2e1 * x0LocalFrame[2] * lLocalFrame[2];
       c = t9 - 0.25e0 * t10 + t12 - 0.25e0 * t13 - 0.50e0 * r0Surface * r1Surface;
     }



     d=b*b-4.0*a*c;

     if (d<0.0) {
       if (4.0*a*pow(EPS,2)>-d) d=0.0; //equvalent to EPS/|l[:]| > sqrt(|d|)/(2a) -> the particle is within the distance of EPS from the surface's boundary
       else return false;
     }

     if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

     sqrt_d=sqrt(d);

     if ((b+sqrt_d>10.0*fabs(a)) && (fabs(b-sqrt_d)>10.0*fabs(a))) return false;

     t1=-(b+sqrt_d)/(2.0*a);
     t2=-2.0*c/(b+sqrt_d);

     double t=-1.0;

     if (t1>0.0) t=t1;
     if ( (t<0.0) || ((t2>0.0)&&(t>t2)) ) t=t2;


     if ((-EPS<t)&&(t<1.0+EPS)) {
       double xLocIntersection;

       for (idim=0;idim<3;idim++) IntersectionPoint_GlobalFrame[idim]=x0LinearSegment_GlobalFrame[idim]+t*l[idim];

       xLocIntersection=(IntersectionPoint_GlobalFrame[0]-OriginPosition[0])*e0[0]+
           (IntersectionPoint_GlobalFrame[1]-OriginPosition[1])*e0[1]+
           (IntersectionPoint_GlobalFrame[2]-OriginPosition[2])*e0[2];

       if ((xLocIntersection<x0SurfaceMin)||(x0SurfaceMax<xLocIntersection)) return false;

       return true;
     }




       return false;
   }


   bool LinearSegment_Surface_Intersection(double *x0LinearSegment_GlobalFrame,double *x1LinearSegment_GlobalFrame,double *IntersectionPoint_GlobalFrame) {
     int iAxis;

     for (iAxis=0;iAxis<nAxisSurfaceElements;iAxis++) {
       if (LinearSegment_Surface_Intersection(xAxisSurfaceElement[iAxis],xAxisSurfaceElement[iAxis+1],x0LinearSegment_GlobalFrame,x1LinearSegment_GlobalFrame,IntersectionPoint_GlobalFrame)==true) {
         return true;
       }
     }


     return false;
   }

   //=======================================================================
   //intersection of a block with the sphere
   int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
     int idim,iSurface,i,j,k,nCounter;
//     double x[3],locx[3],r0,rSurface;
//     double xBlockProjectionMax=0.0,xBlockProjectionMin=0.0;

     bool IntersectionPossibleFlag=false;
     int iSurfaceStart=-1,iSurfaceFinish=-1;


     //1 .check if the center point of the block is "far" from the surface
     //calculate positions of the block's nodes in the frame retlated to the 'body of rotation'
     double xBlockProjectionMin,xBlockProjectionMax,rBlockMin,rBlockMax;
     double x[3];

     static const int nNode[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},  {0,0,1},{1,0,1},{1,1,1},{0,1,1}};


     for (int pnode=0;pnode<8;pnode++) {
       double c0=0.0,c1=0.0;

       for (idim=0;idim<3;idim++) {
         x[idim]=((nNode[pnode][idim]==0) ? xBlockMin[idim] : xBlockMax[idim])-OriginPosition[idim];
         c0+=x[idim]*e0[idim];
       }

       for (idim=0;idim<3;idim++) c1+=pow(x[idim]-c0*e0[idim],2);

       c1=sqrt(c1);

       if (pnode==0) xBlockProjectionMin=c0,xBlockProjectionMax=c0,rBlockMin=c1,rBlockMax=c1;
       else {
         if (xBlockProjectionMin>c0) xBlockProjectionMin=c0;
         if (xBlockProjectionMax<c0) xBlockProjectionMax=c0;

         if (rBlockMin>c1) rBlockMin=c1;
         if (rBlockMax<c1) rBlockMax=c1;
       }
     }


     if (!((xBlockProjectionMax<xAxisMin)||(xAxisMax<xBlockProjectionMin))) {
       if (xBlockProjectionMin<xAxisMin) iSurfaceStart=0;
       if (xBlockProjectionMax>xAxisMax) iSurfaceFinish=nAxisSurfaceElements-1;


       for (iSurface=0;iSurface<nAxisSurfaceElements;iSurface++) {
         if ((xAxisSurfaceElement[iSurface]<=xBlockProjectionMin)&&(xBlockProjectionMin<xAxisSurfaceElement[iSurface+1])) iSurfaceStart=iSurface;
         if ((xAxisSurfaceElement[iSurface]<=xBlockProjectionMax)&&(xBlockProjectionMax<xAxisSurfaceElement[iSurface+1])) iSurfaceFinish=iSurface;

         if ((iSurfaceStart!=-1)&&(IntersectionPossibleFlag==false)) {
           //check if the intersection is possible
           double rmin,rmax,r0,r1;

           SurfaceCurve(r0,xAxisSurfaceElement[iSurface]);
           SurfaceCurve(r1,xAxisSurfaceElement[iSurface+1]);

           rmin=min(r0,r1);
           rmax=max(r0,r1);

           if ( ((rBlockMin<=rmin)&&(rmin<=rBlockMax)) || ((rBlockMin<=rmax)&&(rmax<=rBlockMax)) || ((rmin<=rBlockMin)&&(rBlockMin<=rmax)) || ((rmin<=rBlockMax)&&(rBlockMax<=rmax)) ) {
             IntersectionPossibleFlag=true;
           }

         }

         if ((iSurfaceStart!=-1)&&(iSurfaceFinish!=-1)) break;
       }
     }


     if (IntersectionPossibleFlag==true) {
       //the block is 'close' to  the surface ==> check its intersection with the surface
       //2. check if the block intersects the surface elements of the body of rotation
       double x[3],x0Projection,x1Projection,x2Projection,r1Surface,r2Surface,tgThetaCone,ezLocal[3],e1Local[3],e2Local[3],x0[3];
       int pnode;

       static const int nNode[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},  {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

       const int CylinderShape=0;
       const int ConeShape=1;

       int SurfaceElementShape=-1;


       for (iSurface=iSurfaceStart;iSurface<=iSurfaceFinish;iSurface++) {
         x1Projection=xAxisSurfaceElement[iSurface],x2Projection=xAxisSurfaceElement[iSurface+1];
         memcpy(ezLocal,e0,3*sizeof(double));

         if (SurfaceCurve(r1Surface,x1Projection)==false) continue;
         if (SurfaceCurve(r2Surface,x2Projection)==false) continue;

         tgThetaCone=(r2Surface-r1Surface)/(x2Projection-x1Projection);

         if (fabs(tgThetaCone)>1.0E-15) {
           //cone
           SurfaceElementShape=ConeShape;
           x0Projection=x1Projection-r1Surface/tgThetaCone;
           for (idim=0;idim<3;idim++) x0[idim]=OriginPosition[idim]+x0Projection*e0[idim];

           if (tgThetaCone>0.0) {
             x1Projection-=x0Projection,x2Projection-=x0Projection;
           }
           else {
             for (idim=0;idim<3;idim++) ezLocal[idim]*=-1.0;

             tgThetaCone*=-1.0;
             x1Projection=x0Projection-x1Projection;
             x2Projection=x0Projection-x2Projection;
           }
         }
         else {
           //cylinder
           SurfaceElementShape=CylinderShape;
           x0Projection=x1Projection;
           x1Projection-=x0Projection,x2Projection-=x0Projection;

           for (idim=0;idim<3;idim++) x0[idim]=OriginPosition[idim]+x0Projection*e0[idim];
         }

         for (pnode=0;pnode<8+1;pnode++) {
           if (pnode<=7) {
             for (idim=0;idim<3;idim++) {
               x[idim]=(nNode[pnode][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
             }
           }
           else for (idim=0;idim<3;idim++) x[idim]=0.5*(xBlockMin[idim]+xBlockMax[idim]);

           //reconstruct the frame of reference where e0Local -> along the axis of symmetry, e1Local -> taward x[], and e2Local = e0Local x e1Local (not used because not needed)
           double r=0.0,c0=0.0,c1=0.0,c2;

           for (idim=0;idim<3;idim++) c0+=ezLocal[idim]*(x[idim]-x0[idim]);

           for (idim=0;idim<3;idim++) {
             e1Local[idim]=x[idim]-x0[idim]-c0*ezLocal[idim];
             c1+=e1Local[idim]*e1Local[idim];
           }

           c1=sqrt(c1);
           if (c1<1.0E-15) continue; //the axis of symmetry is aligned with the vector 'x'-'x0Cone'

           for (idim=0;idim<3;idim++) e1Local[idim]/=c1;

           //get the polar and zenith angle of point 'x' in the spherical frame originated at 'x0Cone' with the orientation (e0Local,e1Local,e2Local);
           double theta;

           for (r=0.0,c1=0.0,idim=0;idim<3;idim++) {
             double t=x[idim]-x0[idim];

             r+=t*t;
             c1+=t*e1Local[idim];
           }

           r=sqrt(r);
           theta=acos(c0/r);

           //1. check the boundary points of the surface elements
           bool flag[2]={false,false};

           if (SurfaceElementShape==ConeShape) {
             for (idim=0;idim<3;idim++) {
               double t;

               t=x0[idim]+x1Projection*(ezLocal[idim]+tgThetaCone*e1Local[idim]);
               if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) flag[0]=true;

               t=x0[idim]+x2Projection*(ezLocal[idim]+tgThetaCone*e1Local[idim]);
               if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) flag[1]=true;
             }
           }
           else { //(SurfaceElementShape==CylinderShape)
             for (idim=0;idim<3;idim++) {
               double t;

               t=x0[idim]+r1Surface*e1Local[idim];
               if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) flag[0]=true;

               t=x0[idim]+x2Projection*ezLocal[idim]+r1Surface*e1Local[idim];
               if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) flag[1]=true;
             }
           }


           if ((flag[0]==false)||(flag[1]==false)) {
             //one of the points is within the block --> the block intersects the surface
             return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
           }

           //2. check the point that is 'closest' to 'x'
           if (SurfaceElementShape==ConeShape) {
             double ThetaCone,dTheta,xClosestProjection;

             ThetaCone=atan(tgThetaCone);
             dTheta=fabs(ThetaCone-theta);

             xClosestProjection=r*cos(dTheta)*cos(ThetaCone);

             if ((x1Projection<xClosestProjection) && (xClosestProjection<x2Projection)) {
               flag[0]=false;

               for (idim=0;idim<3;idim++) {
                 double t;

                 t=x0[idim]+xClosestProjection*(ezLocal[idim]+tgThetaCone*e1Local[idim]);
                 if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) {
                   flag[0]=true;
                   break;
                 }
               }

               if (flag[0]==false) {
                 //one of the points is within the block --> the block intersects the surface
                 return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
               }
             }
           }
           else { //(SurfaceElementShape==CylinderShape)
             if ((0.0<c0)&&(c0<x2Projection)) {
               flag[0]=false;

               for (idim=0;idim<3;idim++) {
                 double t;

                 t=x0[idim]+c0*ezLocal[idim]+r1Surface*e1Local[idim];
                 if ((t<xBlockMin[idim]) || (t>xBlockMax[idim])) {
                   flag[0]=true;
                   break;
                 }
               }

               if (flag[0]==false) {
                 //one of the points is within the block --> the block intersects the surface
                 return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
               }
             }

           }

         }
       }

     }

     //at this point: no intersection between the block and teh surface are determined ==> the block either entirely inside or outside of the body
     //check if point 'xBlockMin' is within the body
     long int nAxisElement,nAzimuthalElement;

     if (GetSurfaceElementProjectionIndex(xBlockMin,nAxisElement,nAzimuthalElement)==true) {
       double c=0.0,r2=0.0;
       double x1Projection,x2Projection,r1Surface,r2Surface,tgThetaCone;

       for (idim=0;idim<3;idim++) c+=(xBlockMin[idim]-OriginPosition[idim])*e0[idim];
       for (idim=0;idim<3;idim++) r2+=pow(xBlockMin[idim]-OriginPosition[idim]-c*e0[idim],2);

       x1Projection=xAxisSurfaceElement[nAxisElement],x2Projection=xAxisSurfaceElement[nAxisElement+1];
       SurfaceCurve(r1Surface,x1Projection);
       SurfaceCurve(r2Surface,x2Projection);

       tgThetaCone=(r2Surface-r1Surface)/(x2Projection-x1Projection);

       if (r1Surface+tgThetaCone*(c-x1Projection)>sqrt(r2)) {
         return _AMR_BLOCK_OUTSIDE_DOMAIN_;
       }
     }

     return _AMR_BLOCK_INSIDE_DOMAIN_;

/*     /////////////

     //2. check if a block's edge intersect the surface
     double x0[3],x1[3],xIntersection[3];

     static const int nX0Edge[12][3]={ {0,0,0},{0,1,0},{0,1,1},{0,0,1}, {0,0,0},{1,0,0},{1,0,1},{0,0,1}, {0,0,0},{1,0,0},{1,1,0},{0,1,0}};
     static const int nX1Edge[12][3]={ {1,0,0},{1,1,0},{1,1,1},{1,0,1}, {0,1,0},{1,1,0},{1,1,1},{0,1,1}, {0,0,1},{1,0,1},{1,1,1},{0,1,1}};


     for (int nedge=0;nedge<12;nedge++) {
       for (idim=0;idim<3;idim++) {
         x0[idim]=(nX0Edge[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
         x1[idim]=(nX1Edge[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
       }

       if (LinearSegment_Surface_Intersection(x0,x1,xIntersection)==true) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
     }


     //3. check intersection of the axis of symmetry with the faces of the block
     //the internal coordinated of the origin of the coordinate frame related to a face
     static const int nX0face[6][3]=  { {0,0,0},{1,0,0}, {0,0,0},{0,1,0}, {0,0,0},{0,0,1}};

     //the internal coordinate of the nodes that determine the coordinate vectors related to the frame
     static const int nE0[6][3]=  { {0,1,0},{1,1,0}, {1,0,0},{1,1,0}, {1,0,0},{1,0,1}};
     static const int nE1[6][3]=  { {0,0,1},{1,0,1}, {0,0,1},{0,1,1}, {0,1,0},{0,1,1}};

     //the direction to the local normal in the coordinate system related to the face
     static const int normFace[6][3]={ {-1,0,0},{1,0,0}, {0,-1,0},{0,1,0}, {0,0,-1},{0,0,1}};




     for (int nface=0;nface<6;nface++) {
       double r0[3],l[3],x0Face[3],e0Face[3],e1Face[3],normFace[3],rE0,rE1,lE0,lE1,lNorm,c,t,length=0.0,lProjection=0.0;

       double xFaceCenter[3],xLocalFaceCenterMin;

       for (idim=0,lNorm=0.0,lE0=0.0,lE1=0.0,xLocalFaceCenterMin=0.0;idim<3;idim++) {
         x0Face[idim]=(nX0face[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
         e0Face[idim]=((nE0[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0Face[idim];
         e1Face[idim]=((nE1[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0Face[idim];

         lE0+=pow(e0Face[idim],2);
         lE1+=pow(e1Face[idim],2);
       }


       //get the vector that is normal to 'e0' and minimize the distance from the axis of symmetry to the face
       double normSurfaceAxis[3];

       for (c=0.0,idim=0;idim<3;idim++) c+=normFace[nface][idim]*e0[idim];
       for (t=0.0,idim=0;idim<3;idim++) {
         normSurfaceAxis[idim]=normFace[nface][idim]-c*e0[idim];
         t+=pow(normSurfaceAxis[idim],2);
       }

       //determine the orientation of 'normSurfaceAxis' -> it should points toward the face
       for (c=0.0,idim=0;idim<3;idim) c+=(x0Face[idim]-OriginPosition[idim])*normSurfaceAxis[idim];

       t=sqrt(t);
       if (c<0.0) t*=-1.0;

       lE0=sqrt(lE0);
       lE1=sqrt(lE1);


       for (idim=0;idim<3;idim++) normSurfaceAxis[idim]/=t,e0Face[idim]/=lE0,e1Face[idim]/=lE1;


       for (i=0;i<nAxisSurfaceElements+1;i++) {
         //determine the point on the surface that corresponds to the surface element
         double rSurfaceGlobal[3],rSurfaceLocal[3]={0.0,0.0,0.0},rSurfaceLocal_last[3],rSurfacePoint_PositionIndicator;

         SurfaceCurve(rSurface,locx[0]);
         for (idim=0,rSurfacePoint_PositionIndicator=0.0;idim<3;idim++) rSurfaceGlobal[idim]=OriginPosition[idim]+xAxisSurfaceElement[i]*e0[idim] + normSurfaceAxis[idim]*rSurface;

         //conver the surface point into the frame related to the face
         for (idim=0;idim<3;idim++) {
           t=rSurfaceGlobal[idim]-OriginPosition[idim];

           rSurfaceLocal[0]+=t*e0Face[idim];
           rSurfaceLocal[1]+=t*e1Face[idim];
           rSurfaceLocal[2]+=t*normFace[nface][idim];
         }

         //check if the last two points are on the same side in respect to the face

       }



       //check if 'xLocalFaceCenterMin' within the surface
       if (SurfaceCurve(rSurface,xLocalFaceCenterMin)==false) continue;  //point 'xLocalFaceCenterMin' does not belong to the surface

       for (idim=0,locx[1]=0.0,locx[2]=0.0;idim<3;idim++) {
         double t;

         t=xFaceCenter[idim]-OriginPosition[idim];
         locx[1]+=t*e1[idim],locx[2]+=t*e2[idim];
       }

       if (rSurface*rSurface<locx[1]*locx[1]+locx[2]*locx[2]) continue; //there is no intersection between the surface and the face

       //intersection is found
       return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
     }

     //check is the surface axis's poins are inside the block
      for (i=0;i<nAxisSurfaceElements+1;i++) {
       for (idim=0;idim<3;idim++) {
         double t=OriginPosition[idim]+xAxisSurfaceElement[i]*e0[idim];

         if ((t<xBlockMin[idim])||(t>xBlockMax[idim])) break;
       }

       if (idim==3) {
         return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
       }
     }



      ////////


      //check if the block's points are inside the surface
      for (i=0,nCounter=0;i<3;i++) {
        x[0]=xBlockMin[0]+0.5*double(i)*(xBlockMax[0]-xBlockMin[0])-OriginPosition[0];

        for (j=0;j<3;j++) {
          x[1]=xBlockMin[1]+0.5*double(j)*(xBlockMax[1]-xBlockMin[1])-OriginPosition[1];

          for (k=0;k<3;k++) {
            x[2]=xBlockMin[2]+0.5*double(k)*(xBlockMax[2]-xBlockMin[2])-OriginPosition[2];

            locx[0]=x[0]*e0[0]+x[1]*e0[1]+x[2]*e0[2];
            locx[1]=x[0]*e1[0]+x[1]*e1[1]+x[2]*e1[2];
            locx[2]=x[0]*e2[0]+x[1]*e2[1]+x[2]*e2[2];

            if (SurfaceCurve(rSurface,locx[0])==false) continue;
            else {
              if (rSurface*rSurface>=locx[1]*locx[1]+locx[2]*locx[2]) {
                nCounter++;
              }
            }
          }
        }
      }

      if (nCounter==27) return _AMR_BLOCK_OUTSIDE_DOMAIN_;
     return _AMR_BLOCK_INSIDE_DOMAIN_;*/



   }


   double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
     int BlockIntersectionCode;
     double res=0.0,TotalResult=0.0;
     int idim;

     static const int nLevelMax=5;
     static const int nIntegrationIntervals=5;

     struct cLevelData {
       double xSubBlockMin[3],xSubBlockMax[3],dx[3];
       int i,j,k;
     };

     cLevelData LevelData[nLevelMax+1];
     cLevelData *levelDataPtr;
     int nLevel=0;


     //init the level 0 data
     for (levelDataPtr=LevelData,idim=0;idim<3;idim++) levelDataPtr->xSubBlockMax[idim]=xBlockMaxInit[idim],levelDataPtr->xSubBlockMin[idim]=xBlockMinInit[idim];


FunctionBegins:
     res=0.0;
     levelDataPtr=LevelData+nLevel;

     BlockIntersectionCode=BlockIntersection(levelDataPtr->xSubBlockMin,levelDataPtr->xSubBlockMax,EPS);
     *IntersectionStatus=BlockIntersectionCode;

     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) {
       if (nLevel==0) return 0.0;
       else {
         res=0.0;
         goto LevelProcessingDone;
       }
     }
     else if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
       for (res=1.0,idim=0;idim<3;idim++) res*=levelDataPtr->xSubBlockMax[idim]-levelDataPtr->xSubBlockMin[idim];

       #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
       if (res<0.0) exit(__LINE__,__FILE__,"Error: out of range");
       #endif

       if (nLevel==0) return res;
       else {
         goto LevelProcessingDone;
       }
     }


     //the block is intersected by the sphere
     if (nLevel<nLevelMax) {
       for (idim=0;idim<3;idim++) levelDataPtr->dx[idim]=(levelDataPtr->xSubBlockMax[idim]-levelDataPtr->xSubBlockMin[idim])/2.0;

       for (levelDataPtr->i=0;levelDataPtr->i<2;levelDataPtr->i++) {
         (levelDataPtr+1)->xSubBlockMin[0]=levelDataPtr->xSubBlockMin[0]+levelDataPtr->i*levelDataPtr->dx[0];
         (levelDataPtr+1)->xSubBlockMax[0]=levelDataPtr->xSubBlockMin[0]+(levelDataPtr->i+1)*levelDataPtr->dx[0];

         for (levelDataPtr->j=0;levelDataPtr->j<2;levelDataPtr->j++) {
           (levelDataPtr+1)->xSubBlockMin[1]=levelDataPtr->xSubBlockMin[1]+levelDataPtr->j*levelDataPtr->dx[1];
           (levelDataPtr+1)->xSubBlockMax[1]=levelDataPtr->xSubBlockMin[1]+(levelDataPtr->j+1)*levelDataPtr->dx[1];

           for (levelDataPtr->k=0;levelDataPtr->k<2;levelDataPtr->k++) {
             (levelDataPtr+1)->xSubBlockMin[2]=levelDataPtr->xSubBlockMin[2]+levelDataPtr->k*levelDataPtr->dx[2];
             (levelDataPtr+1)->xSubBlockMax[2]=levelDataPtr->xSubBlockMin[2]+(levelDataPtr->k+1)*levelDataPtr->dx[2];

             nLevel++;
             res=0.0;
             goto FunctionBegins;

LevelProcessingDone:
             TotalResult+=res;
             res=0.0;
             nLevel--;
             levelDataPtr=LevelData+nLevel;
           }
         }
       }
     }
     else {
       double dii=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])/nIntegrationIntervals;
       double djj=(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])/nIntegrationIntervals;
       double dkk=(levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2])/nIntegrationIntervals;
       int ii,jj,kk;
       bool IntersectionFound=false;
       double xLocal[3],x[3],rSurface,kkRes,jjRes;

       for (ii=0,res=0.0;ii<nIntegrationIntervals+1;ii++) {
         x[0]=levelDataPtr->xSubBlockMin[0]+ii*dii;

         for (jj=0,jjRes=0.0;jj<nIntegrationIntervals+1;jj++) {
           x[1]=levelDataPtr->xSubBlockMin[1]+jj*djj;

           for (kk=0,kkRes=0.0;kk<nIntegrationIntervals+1;kk++) {
             x[2]=levelDataPtr->xSubBlockMin[2]+kk*dkk;

             //transform the start point 'x0' to the frame related to the surface
             for (xLocal[0]=0.0,xLocal[1]=0.0,xLocal[2]=0.0,idim=0;idim<3;idim++) {
               double dx=x[idim]-OriginPosition[idim];

               xLocal[0]+=dx*e0[idim];
               xLocal[1]+=dx*e1[idim];
               xLocal[2]+=dx*e2[idim];
             }

             if (SurfaceCurve(rSurface,xLocal[0])==false) continue;
             else {
               if (rSurface*rSurface>=xLocal[1]*xLocal[1]+xLocal[2]*xLocal[2]) kkRes+=dkk*(((kk==0)||(kk==nIntegrationIntervals)) ? 0.5 : 1.0);
             }
           }

           jjRes+=kkRes*djj*(((jj==0)||(jj==nIntegrationIntervals)) ? 0.5 : 1.0);
         }

         res+=jjRes*dii*(((ii==0)||(ii==nIntegrationIntervals)) ? 0.5 : 1.0);
       }
     }

     if (nLevel!=0) goto LevelProcessingDone;

     return TotalResult;
   }

   double GetRemainedBlockVolumeMC(double *xBlockMin,double *xBlockMax,double EPS,double RelativeError,int *IntersectionStatus,int nTotalTests,int maxIntegrationLevel,int IntegrationLevel=0) {
     int BlockIntersectionCode,ntest;
     double VolumeL0=-1.0,VolumeL1=-1.0;

     BlockIntersectionCode=BlockIntersection(xBlockMin,xBlockMax,EPS);
     if (IntegrationLevel==0) *IntersectionStatus=BlockIntersectionCode;

     if ((BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_)||(BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_)) {


//TEST: BEGIN:
       {

         int cnt=0;

         for (int ii=0;ii<50000;ii++) {
         double x[3],r,locx[3]={0.0,0.0,0.0};

         for (int i=0;i<3;i++) x[i]=xBlockMin[i]+rnd()*(xBlockMax[i]-xBlockMin[i]);

         for (int i=0;i<3;i++) {
           locx[0]+=e0[i]*(x[i]-OriginPosition[i]);
           locx[1]+=e1[i]*(x[i]-OriginPosition[i]);
           locx[2]+=e2[i]*(x[i]-OriginPosition[i]);
         }

         if (SurfaceCurve(r,locx[0])==true) {
           if (r*r>=locx[1]*locx[1]+locx[2]*locx[2]) cnt++;
         }

         }

         if ((cnt!=0) && (cnt!=50000)) {
           BlockIntersectionCode=BlockIntersection(xBlockMin,xBlockMax,EPS);
         }

       }


//TEST: END

       return (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) ? 0.0 : (xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2]);
     }


     //perform integration
     if ((IntegrationLevel==0)||(IntegrationLevel==maxIntegrationLevel)) {
       for (ntest=0,VolumeL1=0.0;ntest<nTotalTests;ntest++) {
         double x[3],r,locx[3]={0.0,0.0,0.0};

         for (int i=0;i<3;i++) x[i]=xBlockMin[i]+rnd()*(xBlockMax[i]-xBlockMin[i]);

         for (int i=0;i<3;i++) {
           locx[0]+=e0[i]*(x[i]-OriginPosition[i]);
           locx[1]+=e1[i]*(x[i]-OriginPosition[i]);
           locx[2]+=e2[i]*(x[i]-OriginPosition[i]);
         }

         if (SurfaceCurve(r,locx[0])==true) {
           if (r*r>=locx[1]*locx[1]+locx[2]*locx[2]) VolumeL1++;
         }

       }

       VolumeL1=(1.0-VolumeL1/nTotalTests)*(xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2]);
     }

     if (IntegrationLevel!=maxIntegrationLevel) {
       int ii,jj,kk;

       if (IntegrationLevel==0) {
         if (maxIntegrationLevel==0) {
           return VolumeL1;
         }
         else {
           int upperIntegrationLevel=0;

           do {
             double x0[3],x1[3];

             VolumeL0=VolumeL1;
             VolumeL1=0.0;
             upperIntegrationLevel+=1;

             for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
               x0[0]=xBlockMin[0]+0.5*ii*(xBlockMax[0]-xBlockMin[0]);
               x0[1]=xBlockMin[1]+0.5*jj*(xBlockMax[1]-xBlockMin[1]);
               x0[2]=xBlockMin[2]+0.5*kk*(xBlockMax[2]-xBlockMin[2]);

               x1[0]=x0[0]+0.5*(xBlockMax[0]-xBlockMin[0]);
               x1[1]=x0[1]+0.5*(xBlockMax[1]-xBlockMin[1]);
               x1[2]=x0[2]+0.5*(xBlockMax[2]-xBlockMin[2]);

               VolumeL1+=GetRemainedBlockVolumeMC(x0,x1,EPS,RelativeError,NULL,nTotalTests,upperIntegrationLevel,1);
             }

             if (upperIntegrationLevel==maxIntegrationLevel) break;
             if (VolumeL1<1.0E-6*(xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2])) break;
             if (2.0*fabs(VolumeL1-VolumeL0)/(VolumeL1+VolumeL0)<RelativeError) break;
           } while (true);

           return VolumeL1;
         }
       }
       else {
         double x0[3],x1[3];
          VolumeL1=0.0;

          for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
            x0[0]=xBlockMin[0]+0.5*ii*(xBlockMax[0]-xBlockMin[0]);
            x0[1]=xBlockMin[1]+0.5*jj*(xBlockMax[1]-xBlockMin[1]);
            x0[2]=xBlockMin[2]+0.5*kk*(xBlockMax[2]-xBlockMin[2]);

            x1[0]=x0[0]+0.5*(xBlockMax[0]-xBlockMin[0]);
            x1[1]=x0[1]+0.5*(xBlockMax[1]-xBlockMin[1]);
            x1[2]=x0[2]+0.5*(xBlockMax[2]-xBlockMin[2]);

            VolumeL1+=GetRemainedBlockVolumeMC(x0,x1,EPS,RelativeError,NULL,nTotalTests,maxIntegrationLevel,IntegrationLevel+1);
          }
       }
     }

     return VolumeL1;
   }


   template <class cTetrahedron, class cQuadrilateral,class cCutBlock, class cTriangleCutFace,class cCutBlockNode,class cCutData,class cCutBlockSet>
   double GetRemainedBlockVolumeCutBlock(double *xBlockMin,double *xBlockMax,double EPS,double RelativeError,int *IntersectionStatus,cCutBlockSet* CutBlockSet,int maxIntegrationLevel,int IntegrationLevel=0) {
     int BlockIntersectionCode;
     double VolumeL0=-1.0,VolumeL1=-1.0;

     if (IntegrationLevel==0) CutBlockSet->clear();

     BlockIntersectionCode=BlockIntersection(xBlockMin,xBlockMax,EPS);
     if (IntegrationLevel==0) *IntersectionStatus=BlockIntersectionCode;

     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) {
       return 0.0;
     }
     else if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
       cQuadrilateral t;
       list<cCutBlockNode> NodeList;
       cCutBlockNode nd;

       static const int InternalNodeMap[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0}, {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

       for (int pnode=0;pnode<8;pnode++) {
         for (int i=0;i<3;i++) nd.x[i]=xBlockMin[i]+InternalNodeMap[pnode][i]*(xBlockMax[i]-xBlockMin[i]);

         NodeList.push_front(nd);
         t.node[pnode]=NodeList.begin();
       }

       CutBlockSet->AddQuadrilateral(t,EPS);

       return (xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2]);
     }


     //perform integration
     if ((IntegrationLevel==0)||(IntegrationLevel==maxIntegrationLevel)) {
       CutCell::cCutBlock bl;
       list<CutCell::cTetrahedron> indomainConnectivityList,outdomainConnectivityList;
       list<CutCell::cTriangleCutFace> TriangleCutConnectivity;
       typename list<CutCell::cTetrahedron>::iterator itr;

       VolumeL1=0.0;

       //populate the nodes of the block
       //nodes
       for (int i=0;i<2;i+=1) for (int j=0;j<2;j+=1) for (int k=0;k<2;k+=1) {
         double r,x[3],locx[3]={0.0,0.0,0.0};

         x[0]=xBlockMin[0]+i*(xBlockMax[0]-xBlockMin[0]);
         x[1]=xBlockMin[1]+j*(xBlockMax[1]-xBlockMin[1]);
         x[2]=xBlockMin[2]+k*(xBlockMax[2]-xBlockMin[2]);

         for (int ii=0;ii<3;ii++) {
           locx[0]+=e0[ii]*(x[ii]-OriginPosition[ii]);
           locx[1]+=e1[ii]*(x[ii]-OriginPosition[ii]);
           locx[2]+=e2[ii]*(x[ii]-OriginPosition[ii]);
         }

         if (SurfaceCurve(r,locx[0])) {
           if (r*r<=locx[1]*locx[1]+locx[2]*locx[2]) bl.AddNode(x,2*i,2*j,2*k);
         }
       }

       //add cut point on the edges
       //edges
       static const int EdgeMiddleNodeMap[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},  {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
       static const int x0EdgeMap[12][3]={{0,0,0},{0,2,0},{0,2,2},{0,0,2},  {0,0,0},{2,0,0},{2,0,2},{0,0,2},  {0,0,0},{2,0,0},{2,2,0},{0,2,0}};
       static const int x1EdgeMap[12][3]={{2,0,0},{2,2,0},{2,2,2},{2,0,2},  {0,2,0},{2,2,0},{2,2,2},{0,2,2},  {0,0,2},{2,0,2},{2,2,2},{0,2,2}};

       for (int nedge=0;nedge<12;nedge++) {
         double x0[3],x1[3],xIntersection[3];

         for (int iii=0;iii<3;iii++) {
           x0[iii]=xBlockMin[iii]+0.5*x0EdgeMap[nedge][iii]*(xBlockMax[iii]-xBlockMin[iii]);
           x1[iii]=xBlockMin[iii]+0.5*x1EdgeMap[nedge][iii]*(xBlockMax[iii]-xBlockMin[iii]);
         }

         if (LinearSegment_Surface_Intersection(x0,x1,xIntersection)==true) {
           int cnt=0;

           if (bl.node[x0EdgeMap[nedge][0]][x0EdgeMap[nedge][1]][x0EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
           if (bl.node[x1EdgeMap[nedge][0]][x1EdgeMap[nedge][1]][x1EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
           if ((cnt==0)||(cnt==2)) continue;

           bl.AddNode(xIntersection,EdgeMiddleNodeMap[nedge][0],EdgeMiddleNodeMap[nedge][1],EdgeMiddleNodeMap[nedge][2]);
         }
       }

       memcpy(bl.xBlockMin,xBlockMin,3*sizeof(double));
       memcpy(bl.xBlockMax,xBlockMax,3*sizeof(double));
       for (int i=0;i<3;i++) bl.dxBlock[i]=xBlockMax[i]-xBlockMin[i];

       CutCell::cutBlockTetrahedronConnectivity(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);
       CutBlockSet->AddTetrahedronList(indomainConnectivityList,EPS);

       for (itr=indomainConnectivityList.begin();itr!=indomainConnectivityList.end();itr++) {
         VolumeL1+=itr->Volume();
       }
     }

     if (IntegrationLevel!=maxIntegrationLevel) {
       int ii,jj,kk;

       if (IntegrationLevel==0) {
         if (maxIntegrationLevel==0) {
           return VolumeL1;
         }
         else {
           int upperIntegrationLevel=0;

           do {
             double x0[3],x1[3];

             VolumeL0=VolumeL1;
             VolumeL1=0.0;
             upperIntegrationLevel+=1;
             CutBlockSet->clear();

             for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
               x0[0]=xBlockMin[0]+0.5*ii*(xBlockMax[0]-xBlockMin[0]);
               x0[1]=xBlockMin[1]+0.5*jj*(xBlockMax[1]-xBlockMin[1]);
               x0[2]=xBlockMin[2]+0.5*kk*(xBlockMax[2]-xBlockMin[2]);

               x1[0]=x0[0]+0.5*(xBlockMax[0]-xBlockMin[0]);
               x1[1]=x0[1]+0.5*(xBlockMax[1]-xBlockMin[1]);
               x1[2]=x0[2]+0.5*(xBlockMax[2]-xBlockMin[2]);

               VolumeL1+=GetRemainedBlockVolumeCutBlock<cTetrahedron,cQuadrilateral,cCutBlock,cTriangleCutFace,cCutBlockNode,cCutData,cCutBlockSet>(x0,x1,EPS,RelativeError,IntersectionStatus,CutBlockSet,upperIntegrationLevel,1);
             }

             if (upperIntegrationLevel==maxIntegrationLevel) break;
             if (VolumeL1<1.0E-6*(xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2])) break;
             if (2.0*fabs(VolumeL1-VolumeL0)/(VolumeL1+VolumeL0)<RelativeError) break;
           } while (true);

           return VolumeL1;
         }
       }
       else {
         double x0[3],x1[3];
          VolumeL1=0.0;

          for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
            x0[0]=xBlockMin[0]+0.5*ii*(xBlockMax[0]-xBlockMin[0]);
            x0[1]=xBlockMin[1]+0.5*jj*(xBlockMax[1]-xBlockMin[1]);
            x0[2]=xBlockMin[2]+0.5*kk*(xBlockMax[2]-xBlockMin[2]);

            x1[0]=x0[0]+0.5*(xBlockMax[0]-xBlockMin[0]);
            x1[1]=x0[1]+0.5*(xBlockMax[1]-xBlockMin[1]);
            x1[2]=x0[2]+0.5*(xBlockMax[2]-xBlockMin[2]);

            VolumeL1+=GetRemainedBlockVolumeCutBlock<cTetrahedron,cQuadrilateral,cCutBlock,cTriangleCutFace,cCutBlockNode,cCutData,cCutBlockSet>(x0,x1,EPS,RelativeError,IntersectionStatus,CutBlockSet,maxIntegrationLevel,IntegrationLevel+1);
          }
       }
     }

     return VolumeL1;
   }


   void GetSurfaceTriangulation(CutCell::cTriangleFace* &SurfaceTriangulation,int &nSurfaceElements) {
     int iAxis,iAzimuthal,el,i;
     double x0[3],x1[3],x2[3],x3[3],FaceNorm[3],c;

     if (SurfaceTriangulation!=NULL) exit(__LINE__,__FILE__,"Error: redefinition of the surface triangulation");

     nSurfaceElements=2*nAxisSurfaceElements*nAzimuthalSurfaceElements;
     SurfaceTriangulation=new CutCell::cTriangleFace[nSurfaceElements];


     for (el=0,iAxis=0;iAxis<nAxisSurfaceElements;iAxis++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
       GetSurfaceCoordinate(x0,iAxis,iAzimuthal);
       GetSurfaceCoordinate(x1,iAxis,((iAzimuthal!=nAzimuthalSurfaceElements-1) ? iAzimuthal+1 : 0));
       GetSurfaceCoordinate(x2,iAxis+1,((iAzimuthal!=nAzimuthalSurfaceElements-1) ? iAzimuthal+1 : 0));
       GetSurfaceCoordinate(x3,iAxis+1,iAzimuthal);

       GetSurfaceElementNormal(FaceNorm,iAxis,iAzimuthal);

       (SurfaceTriangulation+el)->SetFaceNodes(x0,x1,x2);
       (SurfaceTriangulation+el)->attribute=0;
       for (c=0.0,i=0;i<3;i++) c+=FaceNorm[0]*(SurfaceTriangulation+el)->ExternalNormal[i];
       if (c>0) for (i=0;i<3;i++) (SurfaceTriangulation+el)->ExternalNormal[i]*=-1.0;
       ++el;

       (SurfaceTriangulation+el)->SetFaceNodes(x0,x3,x2);
       (SurfaceTriangulation+el)->attribute=0;
       for (c=0.0,i=0;i<3;i++) c+=FaceNorm[0]*(SurfaceTriangulation+el)->ExternalNormal[i];
       if (c>0) for (i=0;i<3;i++) (SurfaceTriangulation+el)->ExternalNormal[i]*=-1.0;
       ++el;
     }
   }


};





#endif
