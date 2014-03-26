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

       SurfaceCurve(r0,xAxisSurfaceElement[iAxis]);
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
      double x0,x1,x0Cone=xAxisSurfaceElement[nAxisElement]-alpha/r0;

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

  inline void GetSurfaceElementProjectionIndex(double *x,long int &nAxisElement,long int &nAzimuthalElement) {
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

    if (flag==false) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    r=sqrt(locx[1]*locx[1]+locx[2]*locx[2]);
    AzimuthAngle=acos(locx[1]/r);
    if (locx[2]<0.0) AzimuthAngle=2.0*Pi-AzimuthAngle;

    nAzimuthalElement=(long int)(AzimuthAngle/dAzimuthalAngle);
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

     double a,b,c,alphaCone,t1,t2;
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

       double t1 = alphaCone * alphaCone;
       double t2 = pow(lLocalFrame[0], 0.2e1);
       double t4 = pow(lLocalFrame[2], 0.2e1);
       double t5 = pow(lLocalFrame[1], 0.2e1);
       double t9 = t1 * lLocalFrame[0];
       double t10 = r1Surface - r0Surface;
       double t11 = 0.1e1 / t10;
       double t12 = r0Surface * t11;
       double t13 = t12 * x0Surface_LocalFrame;
       double t15 = t12 * x1Surface_LocalFrame;
       double t18 = t1 * x0LocalFrame[0];
       double t21 = pow(x0LocalFrame[1], 0.2e1);
       double t22 = r0Surface * r0Surface;
       double t23 = t1 * t22;
       double t24 = t10 * t10;
       double t25 = 0.1e1 / t24;
       double t30 = t1 * r0Surface;
       double t35 = pow(x0LocalFrame[2], 0.2e1);
       double t38 = x0Surface_LocalFrame * x0Surface_LocalFrame;
       double t43 = x0Surface_LocalFrame * x0Surface_LocalFrame;
       double t45 = pow(x0LocalFrame[0], 0.2e1);
       double t53 = x1Surface_LocalFrame * x1Surface_LocalFrame;
       double t56 = t21 + 0.2e1 * t23 * t25 * x1Surface_LocalFrame * x0Surface_LocalFrame + 0.2e1 * t30 * t11 * x1Surface_LocalFrame * x0Surface_LocalFrame + t35 - 0.2e1 * t18 * t15 - t23 * t25 * t38 + 0.2e1 * t18 * x0Surface_LocalFrame - t1 * t43 - t1 * t45 + 0.2e1 * t18 * t13 - 0.2e1 * t30 * t11 * x0Surface_LocalFrame * x0Surface_LocalFrame - t23 * t25 * t53;

       a = -t1 * t2 + t4 + t5;
       b = 0.2e1 * x0LocalFrame[1] * lLocalFrame[1] + 0.2e1 * x0LocalFrame[2] * lLocalFrame[2] + 0.2e1 * t9 * t13 - 0.2e1 * t9 * t15 + 0.2e1 * t9 * x0Surface_LocalFrame - 0.2e1 * t18 * lLocalFrame[0];
       c = t56;
     }
     else {
       double R=0.5*(r0Surface+r1Surface);

       double t2 = pow(lLocalFrame[2], 0.2e1);
       double t3 = pow(lLocalFrame[1], 0.2e1);
       double t10 = R * R;
       double t11 = pow(x0LocalFrame[2], 0.2e1);
       double t12 = pow(x0LocalFrame[1], 0.2e1);
       a = t2 + t3;
       b = 0.2e1 * x0LocalFrame[1] * lLocalFrame[1] + 0.2e1 * x0LocalFrame[2] * lLocalFrame[2];
       c = - t10 + t11 + t12;
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
     int idim,i,j,k,nCounter;
     double x[3],locx[3],r0,rSurface;
     double xBlockProjectionMax=0.0,xBlockProjectionMin=0.0;

     //check if the surface in entirely inside the block
     for (nCounter=0,i=0;i<3;i++) if ((xBlockMin[i]<=OriginPosition[i])&&(OriginPosition[i]<=xBlockMax[i])) nCounter++;
     if (nCounter==3) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;

     //if all corners of the block are within the body
     int nInsideCounter=0,nAboveMaxSurfaceRadiusCounter=0,nAboveMaxXCounter=0,nBelowMinXCounter=0;

     for (nCounter=0,i=0;i<3;i++) {
       x[0]=xBlockMin[0]+0.5*double(i)*(xBlockMax[0]-xBlockMin[0])-OriginPosition[0];

       for (j=0;j<3;j++) {
         x[1]=xBlockMin[1]+0.5*double(j)*(xBlockMax[1]-xBlockMin[1])-OriginPosition[1];

         for (k=0;k<3;k++) {
           x[2]=xBlockMin[2]+0.5*double(k)*(xBlockMax[2]-xBlockMin[2])-OriginPosition[2];

           locx[0]=x[0]*e0[0]+x[1]*e0[1]+x[2]*e0[2];

           for (r0=0.0,idim=0;idim<3;idim++) r0+=pow(x[idim]-locx[0]*e0[idim],2);
           r0=sqrt(r0);

           if (r0>rSurfaceMax) nAboveMaxSurfaceRadiusCounter++;

           if (locx[0]<xAxisMin) nBelowMinXCounter++;
           else if (locx[0]>xAxisMax) nAboveMaxXCounter++;
           else {
//             locx[1]=x[0]*e1[0]+x[1]*e1[1]+x[2]*e1[2];
//             locx[2]=x[0]*e2[0]+x[1]*e2[1]+x[2]*e2[2];

             SurfaceCurve(rSurface,locx[0]);
             if (r0<=rSurface) nInsideCounter++;
           }

           //calcualte the projection of the block's corners onto the axis of symmetry
           if ((i==0)&&(j==0)&&(k==0)) xBlockProjectionMax=locx[0],xBlockProjectionMin=locx[0];
           else {
             if (xBlockProjectionMax<locx[0]) xBlockProjectionMax=locx[0];
             if (xBlockProjectionMin>locx[0]) xBlockProjectionMin=locx[0];
           }
         }
       }
     }

     if (nInsideCounter==27) return _AMR_BLOCK_OUTSIDE_DOMAIN_;
     else if ((nAboveMaxSurfaceRadiusCounter==27)||(nAboveMaxXCounter==27)||(nBelowMinXCounter==27)) return _AMR_BLOCK_INSIDE_DOMAIN_;


     //check if the surface intersects the block
     double x0[3],x1[3],xIntersection[3];
     int iAxis,iAzimuth;
     double xSurface[3];

     for (iAxis=0;iAxis<nAxisSurfaceElements+1;iAxis++) {
        if ((xAxisSurfaceElement[iAxis]<xBlockProjectionMin)||(xBlockProjectionMax<xAxisSurfaceElement[iAxis])) continue;

        for (iAzimuth=0;iAzimuth<nAzimuthalSurfaceElements;iAzimuth++) {
          GetSurfaceCoordinate(xSurface,iAxis,iAzimuth);

          for (i=0;i<3;i++) if (! ((xBlockMin[i]<=xSurface[i])&&(xSurface[i]<=xBlockMax[i]))) continue;
          return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
        }
     }




     //internal nodes of the block that determine the edge
     static const int nX0Edge[12][3]={ {0,0,0},{0,1,0},{0,1,1},{0,0,1}, {0,0,0},{1,0,0},{1,0,1},{0,0,1}, {0,0,0},{1,0,0},{1,1,0},{0,1,0}};
     static const int nX1[12][3]={ {1,0,0},{1,1,0},{1,1,1},{1,0,1}, {0,1,0},{1,1,0},{1,1,1},{0,1,1}, {0,0,1},{1,0,1},{1,1,1},{0,1,1}};


     for (iAxis=0;iAxis<nAxisSurfaceElements;iAxis++) {
       if ( (xAxisSurfaceElement[iAxis+1]<xBlockProjectionMin) || (xAxisSurfaceElement[iAxis]>xBlockProjectionMax)) continue;

       for (int nedge=0;nedge<12;nedge++) {
         for (idim=0;idim<3;idim++) {
           x0[idim]=(nX0Edge[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
           x1[idim]=((nX1[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]);
         }

         if (LinearSegment_Surface_Intersection(xAxisSurfaceElement[iAxis],xAxisSurfaceElement[iAxis+1],x0,x1,xIntersection)==true) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
       }
     }


     return _AMR_BLOCK_INSIDE_DOMAIN_;
   }


   double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
     return fabs((xBlockMaxInit[0]-xBlockMinInit[0])*(xBlockMaxInit[1]-xBlockMinInit[1])*(xBlockMaxInit[2]-xBlockMinInit[2]));
     /*     int BlockIntersectionCode;
     double res=0.0,TotalResult=0.0;
     int idim;

     static const int nLevelMax=10;
     static const int nIntegrationIntervalsXYplane=10;
     static const int nIntegrationIterationsZplane=nIntegrationIntervalsXYplane;

     struct cLevelData {
       double xSubBlockMin[3],xSubBlockMax[3],dx,dy;
       int i,j;
     };

     //cLevelData LevelData[max(nLevelMax,nLevelEPS) +1];


     cLevelData LevelData[nLevelMax+1];
     cLevelData *levelDataPtr;
     int nLevel=0;


     //init the level 0 data
     for (levelDataPtr=LevelData,idim=0;idim<3;idim++) levelDataPtr->xSubBlockMax[idim]=xBlockMaxInit[idim],levelDataPtr->xSubBlockMin[idim]=xBlockMinInit[idim];



FunctionBegins:
     res=0.0;
     levelDataPtr=LevelData+nLevel;



     BlockIntersectionCode=BlockIntersection(levelDataPtr->xSubBlockMin,levelDataPtr->xSubBlockMax,EPS);

//     BlockIntersectionCode=-10;






     *IntersectionStatus=BlockIntersectionCode;




     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) {
       if (nLevel==0) return 0.0;
       else {
         res=0.0;
         goto LevelProcessingDone;
       }
     }


     if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
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



       levelDataPtr->dx=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])/2.0;
       levelDataPtr->dy=(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])/2.0;

       (levelDataPtr+1)->xSubBlockMin[2]=levelDataPtr->xSubBlockMin[2];
       (levelDataPtr+1)->xSubBlockMax[2]=levelDataPtr->xSubBlockMax[2];

       for (levelDataPtr->i=0;levelDataPtr->i<2;levelDataPtr->i++) {
         (levelDataPtr+1)->xSubBlockMin[0]=levelDataPtr->xSubBlockMin[0]+levelDataPtr->i*levelDataPtr->dx;
         (levelDataPtr+1)->xSubBlockMax[0]=levelDataPtr->xSubBlockMin[0]+(levelDataPtr->i+1)*levelDataPtr->dx;

         for (levelDataPtr->j=0;levelDataPtr->j<2;levelDataPtr->j++) {
           (levelDataPtr+1)->xSubBlockMin[1]=levelDataPtr->xSubBlockMin[1]+levelDataPtr->j*levelDataPtr->dy;
           (levelDataPtr+1)->xSubBlockMax[1]=levelDataPtr->xSubBlockMin[1]+(levelDataPtr->j+1)*levelDataPtr->dy;

           nLevel++;
           res=0.0;
           goto FunctionBegins;

//           res+=GetRemainedBlockVolume(xSubBlockMin,xSubBlockMax,EPS,&BlockIntersectionCode,nLevel+1);

LevelProcessingDone:
           TotalResult+=res;
           res=0.0;
           nLevel--;
           levelDataPtr=LevelData+nLevel;



         }
       }
     }
     else {
       double x0[3],l;
       int i;


       double dii=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])/nIntegrationIntervalsXYplane;
       double djj=(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])/nIntegrationIntervalsXYplane;
       int ii,jj;
       bool IntersectionFound=false;

       l=levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2];

       //transfor the pointing direction (z-axis) to the local frame of reference related to the body
       double locxStart[3],locxFinish[3],rStart,rFinish,locl[3]={e0[2],e1[2],e2[2]};
       double rStartSurface,rFinishSurface; //the parameters of the surface

       for (ii=0;ii<nIntegrationIntervalsXYplane+1;ii++) for (jj=0;jj<nIntegrationIntervalsXYplane+1;jj++) {
         x0[0]=levelDataPtr->xSubBlockMin[0]+ii*dii;
         x0[1]=levelDataPtr->xSubBlockMin[1]+jj*djj;
         x0[2]=levelDataPtr->xSubBlockMin[2];

         //transform the start point 'x0' to the frame related to the surface
         for (locxStart[0]=0.0,locxStart[1]=0.0,locxStart[2]=0.0,i=0;i<3;i++) {
           locxStart[0]+=(x0[i]-OriginPosition[i])*e0[i];
           locxStart[1]+=(x0[i]-OriginPosition[i])*e1[i];
           locxStart[2]+=(x0[i]-OriginPosition[i])*e2[i];
         }

         for (i=0;i<3;i++) locxFinish[i]=locxStart[i]+locl[i];

         rStart=sqrt(locxStart[1]*locxStart[1]+locxStart[2]*locxStart[2]);
         rFinish=sqrt(locxFinish[1]*locxFinish[1]+locxFinish[2]*locxFinish[2]);

         //check if the integration interval is outside of the body
         int cnt=0;
         bool StartPointFlag,FinishPointFlag;

         if ((StartPointFlag=SurfaceCurve(rStartSurface,locxStart[0]))==false) cnt++;
         else if (rStart>rStartSurface) {
           cnt++;
           StartPointFlag=false;
         }

         if ((FinishPointFlag=SurfaceCurve(rFinishSurface,locxFinish[0]))==false) cnt++;
         else if (rFinish>rFinishSurface) {
           cnt++;
           FinishPointFlag=false;
         }

         if (cnt==2) {
           //the interval is outside of the body
           res+=1.0;
           IntersectionFound=true;
           continue;
         }
         else if (cnt==0) {
           //the interval is invide the body
           continue;
         }
         else {
           //adjust the integration limits
           bool MiddlePointFlag;
           double MiddlePoint[3],rMiddlePoint,rMiddlePointSurface;
           double xDomainPoint[3],xBodyPoint[3],xInitDomainPoint[3];

           IntersectionFound=true;

           if (StartPointFlag==false) {
             memcpy(xDomainPoint,locxStart,3*sizeof(double));
             memcpy(xInitDomainPoint,locxStart,3*sizeof(double));
             memcpy(xBodyPoint,locxFinish,3*sizeof(double));
           }
           else {
             memcpy(xDomainPoint,locxFinish,3*sizeof(double));
             memcpy(xInitDomainPoint,locxFinish,3*sizeof(double));
             memcpy(xBodyPoint,locxStart,3*sizeof(double));
           }

           for (int niter=0;niter<nIntegrationIterationsZplane;niter++) {
             for (i=0;i<3;i++) MiddlePoint[i]=0.5*(xDomainPoint[i]+xBodyPoint[i]);
             rMiddlePoint=sqrt(MiddlePoint[1]*MiddlePoint[1]+MiddlePoint[2]*MiddlePoint[2]);

             if ((MiddlePointFlag=SurfaceCurve(rMiddlePointSurface,MiddlePoint[0]))==true) if (rMiddlePoint>rMiddlePointSurface) MiddlePointFlag=false;

             if (MiddlePointFlag==true) { //the point is inside the body
               memcpy(xBodyPoint,MiddlePoint,3*sizeof(double));
             }
             else {
               memcpy(xDomainPoint,MiddlePoint,3*sizeof(double));
             }
           }

           //increment the volume counter
           res+=fabs((xInitDomainPoint[0]-MiddlePoint[0])*e0[2]+(xInitDomainPoint[1]-MiddlePoint[1])*e1[2]+(xInitDomainPoint[2]-MiddlePoint[2])*e2[2])/l;
         }
       }

       res/=pow(nIntegrationIntervalsXYplane+1,2);
       if (IntersectionFound==false) res=pow(EPS,3.0);
       else res*=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])*(levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2]);

         #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
         if (res<0.0) {
           exit(__LINE__,__FILE__,"Error: out of range");
         }
         #endif


     }

     if (nLevel!=0) goto LevelProcessingDone;


     if (TotalResult<0.0) exit(__LINE__,__FILE__,"Error: out of range");

     return TotalResult;

   }



     */
   /*
   int cutBlockBinaryBlockConnectivity(FILE *nodesListFile,FILE *connectivityListFile,cBlock3d* bl,long int& nodenoCounter,long int& blocknoCounter) {
      int nedge,idim,i,j,k;
      int connectivityLength=0;
      double c,x[3];
      long int elementPosition;



  //-------


  #define _N_TEST_NODES_ 4


  bool PRINTEDFLAG=false;
  int BLOCKPRINTED=0;
  static int PRINTEDCELLS=0;

  //int nodeIDlist[_N_TEST_NODES_]={2082,23977,23676};

  //int nodeIDlist[_N_TEST_NODES_]={359410,359343,353757};
  //int nodeIDlist[_N_TEST_NODES_]={1791752,1790002,1791751};

  //int nodeIDlist[_N_TEST_NODES_]={44770,44773,44772,246963};

  int nodeIDlist[_N_TEST_NODES_]={73780,73780,73780,-1};

  int nnds=0;

  for (int i=0;i<_N_TEST_NODES_;i++) {
    for (int ii=0;ii<3;ii++) for (int jj=0;jj<3;jj++) for (int kk=0;kk<3;kk++) if (bl->node[ii][jj][kk]!=NULL) if (bl->node[ii][jj][kk]->TEMP_IDENTIFYER_NODE==nodeIDlist[i]) ++nnds;

    for (int nedge=0;nedge<12;nedge++) if (bl->edge[nedge]->cutPoint!=NULL) if (bl->edge[nedge]->cutPoint->TEMP_IDENTIFYER_NODE==nodeIDlist[i]) ++nnds;
  }





  ///--------

      int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};
      int nodeEdges[8][3]={{0,4,8},{0,5,9},{5,10,1},{4,11,1},{8,3,7},{3,9,6},{2,10,6},{7,11,2}};

      int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};

      int globalNodeMask[3][3][3];

      globalNodeMask[0][0][0]=0;
      globalNodeMask[2][0][0]=1;
      globalNodeMask[2][2][0]=2;
      globalNodeMask[0][2][0]=3;

      globalNodeMask[0][0][2]=4;
      globalNodeMask[2][0][2]=5;
      globalNodeMask[2][2][2]=6;
      globalNodeMask[0][2][2]=7;



      //the nodes od the block
      #define nConnectionsMax 20
      #define blNodeListLengthMax  20

      class cBlockNode {
      public:
        int nConnections;
        cBlockNode *Connection[nConnectionsMax];
        cAMRnode<cNode> *meshNode;

        bool ghostNode;

        cBlockNode() {nConnections=0,meshNode=NULL,ghostNode=false;}

        void assign(cAMRnode<cNode> *nd) {
          meshNode=nd;
          nConnections=0;
        }

        bool checkConnection(cBlockNode* nd) {
          for (int i=0;i<nConnections;i++) if (Connection[i]==nd) return true;

          return false;
        }

        void Connect(cBlockNode* nd) {
          bool errorflag=false;

          if (checkConnection(nd)==false) {
            if (nConnections==nConnectionsMax) errorflag=true;
            else Connection[nConnections++]=nd;
          }

          if (nd->checkConnection(this)==false) {
            if (nd->nConnections==nConnectionsMax) errorflag=true;
            else nd->Connection[nd->nConnections++]=this;
          }

          if (errorflag==true) {
            cout << "ERROR: too many connections (file=" << __FILE__ << ", line=" << __LINE__ << endl;
            ::exit(0);
          }
        }

        void removeConnection(cBlockNode* nd) {
          bool foundFirst=false,foundSecond=false;
          int i;

          for (i=0;i<nConnections;i++) if (Connection[i]==nd) {
            for (++i;i<nConnections;i++) Connection[i-1]=Connection[i];
            nConnections--;
            foundFirst=true;
          }

          for (i=0;i<nd->nConnections;i++) if (nd->Connection[i]==this) {
            for (++i;i<nd->nConnections;i++) nd->Connection[i-1]=nd->Connection[i];
            nd->nConnections--;
            foundSecond=true;
          }

          if ((foundFirst==false)||(foundSecond==false)) {
            cout << "Error: no connections found (file=" << __FILE__ << ", line=" << __LINE__ << endl;
            ::exit(0);
          }
        }
      };


      int blNodeListLength=0;
      cBlockNode blNodeList[blNodeListLengthMax];
      cBlockNode *nodeMesh[3][3][3];

      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) nodeMesh[i][j][k]=NULL;

      //get the cut information for the block
      cCellCut3d *cutData=NULL;

      //first set the cut information with the block data
      //cutData=bl->cutData;

      //if the block doe's not have a cut information, search the edjes of the block
      if (cutData==NULL) for (nedge=0;nedge<12;nedge++) if (bl->edge[nedge]->cutData!=NULL) {
        cutData=bl->edge[nedge]->cutData;
        break;
      }

      //if a cut is not found, return connectivity list based on four nodes of the block
      if (cutData==NULL) {

        exit(__LINE__,"Error: this is not a cutted block");

        #if _USE_CELLS_ == _USE_CELLS_ON_
        bl->cell->Volume=bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2];
        #endif


      for (i=0;i<3;i+=3) for (j=0;j<3;j++) for (k=0;k<3;k++) if (bl->node[i][j][k]!=NULL) {
       if (bl->node[i][j][k]->nodeno==-1) {
                bl->node[i][j][k]->nodeno=nodenoCounter++;
          elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
                fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

                elementPosition=blocks.GetEntryCountingNumber(bl);
                fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
       }

             }

            blocknoCounter++;

            //chack if all nodes are already counted
            for (int ii=0;ii<3;ii+=2) for (int jj=0;jj<3;jj+=2) for (int kk=0;kk<3;kk+=2) if (bl->node[ii][jj][kk]->nodeno<0) exit(__LINE__,"Un-numbered node is found");

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        return 1;
      }

      //init the mesh
      //corner nodes of the block
      int nnode,iedge;

      #define _IN_DOMAIN_      0
      #define _OUT_DOMAIN_     1
      #define _DEFAULT_DOMAIN_ 2

      int node_in_domain;

      bool cutBlockFlag=true;


      //a block is considered to be cut if at lease two face are cuted
      bool block_cut_flag=false;
      int nFaceIntersection=0;

      for (int nface=0;nface<6;nface++) {
        int nEdgeIntersection=0;

        for (int pedge=0;pedge<4;pedge++) if (bl->edge[faceEdges[nface][pedge]]->cutPoint!=NULL) nEdgeIntersection++;
        if (nEdgeIntersection>=2) nFaceIntersection++;
      }

      if (nFaceIntersection<2) { //the block cannot be considered to be cut: the regular connectivity list is printing
         std::cout << "ERROR: cutBlockConnectivity - found an uncut block" << std::endl;

        #if _USE_CELLS_ == _USE_CELLS_ON_
        bl->cell->Volume=bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2];
        #endif


        for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) if (bl->node[i][j][k]!=NULL) {
          if (bl->node[i][j][k]->nodeno==-1) {
             bl->node[i][j][k]->nodeno=nodenoCounter++;
             elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
             fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

             elementPosition=blocks.GetEntryCountingNumber(bl);
             fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
          }
        }

        blocknoCounter++;

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        return 1;

  cutBlockFlag=false;

      }




  //NEW PROCEDURE FOR POPULATING THE BLOCK'S NODE MESH
  //  1. the cut node on a edge determines weather the nearby corner nodes exists
  //  2. a) for a cut node on a edge : it is removed is two nearby corner nodes are not exists (the whole edge is removed)
  //  2. b) for a cut node on a edge : it is removed is two nearby corner nodes are  exists  (the edge sgould not be cut)

      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) globalNodeMask[i][j][k]=0;

      for (nedge=0;nedge<12;nedge++) {
        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        //check existance of the cornet nodes
        register int pnode;
        int ii,jj,kk;

        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

  //      if ((globalNodeMask[ii][jj][kk]==0)||(globalNodeMask[ii][jj][kk]==1)) {  //the node has not been removed from the block node mesh by this procedure
            if (bl->edge[nedge]->cutPoint!=NULL) {  // if the cut node doe not existes -> noth corner nodes should be on the mesh (if the edge itself on the mesh)
              //check if the corner nodes sould be on the mesh

              x[0]=bl->xmin[0]+bl->dxBlock[0]*ii/2.0;
              x[1]=bl->xmin[1]+bl->dxBlock[1]*jj/2.0;
              x[2]=bl->xmin[2]+bl->dxBlock[2]*kk/2.0;

              cutData=bl->edge[nedge]->cutData;

              for (c=0.0,idim=0;idim<3;idim++) c+=(x[idim]-cutData->x0[idim])*cutData->norm[idim];

                if (c==0.0) continue;

  //            if (cutBlockFlag==true) {



                if (globalNodeMask[ii][jj][kk]==0) globalNodeMask[ii][jj][kk]=(c>1.0E-6) ? -1 : 1;
                else {
                  int code=(c>1.0E-6) ? -1 : 1;

                  if (code!=globalNodeMask[ii][jj][kk]) {
                    code=1;

                    for (long int nSurfaceElement=0;nSurfaceElement<nSurfaceNASTRANelements;nSurfaceElement++) {
                      code=(surfaceNASTRANmesh[nSurfaceElement].checkPointInsideDomain(x,EPS)==true) ? 1 : -1;
                      if (code==-1) break;
                    }

                    globalNodeMask[ii][jj][kk]=code;
                  }

                }


  //            } else globalNodeMask[ii][jj][kk]=1;
            }
  //        }
        }
      }





  //#################  EXPLIVCITLY CHECK THE CORNER AND CUT NODES ################
  {
     //check the corner nodes of the block





  if (bl->TEMP_IDENTIFYER_BLOCK==22469) {
  cout << __LINE__ << endl;
  }

      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) { ///if (globalNodeMask[i][j][k]>=0) {

        globalNodeMask[i][j][k] = (checkPointInsideDomain(bl->node[i][j][k]->x)==true) ? 1 : -1;
      }






      //check the cut nodes

     int cnt,ii,jj,kk,pnode;

     for (nedge=0;nedge<12;nedge++) {
        cnt=0;

        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        //check existance of the cornet nodes
        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

          if (globalNodeMask[ii][jj][kk]==1) cnt++;
        }

        if (cnt==1) {

          //the cut should be here

          //check if the cut is steed up
          if (bl->edge[nedge]->cutData!=NULL) {
            double *n=bl->edge[nedge]->cutData->norm;

            if (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]>0.00001) continue;
          }


          //set the buffers
          if (bl->edge[nedge]->cutPoint==NULL) {
            bl->edge[nedge]->cutPoint=nodes.newElement();
          }

          if (bl->edge[nedge]->cutData==NULL) {
            bl->edge[nedge]->cutData=cutGeometryData.newElement();
          }


          //init the cut data
          double l;

          for (l=0.0,idim=0;idim<3;idim++) {
            bl->edge[nedge]->cutPoint->x[idim]=0.5*(bl->edge[nedge]->node[0]->x[idim]+bl->edge[nedge]->node[1]->x[idim]);
            bl->edge[nedge]->cutPoint->nodeat=bl->node[ii][jj][kk]->nodeat;

            bl->edge[nedge]->cutData->x0[idim]=bl->edge[nedge]->cutPoint->x[idim];
            bl->edge[nedge]->cutData->norm[idim]=bl->edge[nedge]->node[1]->x[idim]-bl->edge[nedge]->node[0]->x[idim];
            l+=pow(bl->edge[nedge]->cutData->norm[idim],2);
          }

          l=sqrt(l);

          if (globalNodeMask[ii][jj][kk]==1) {
            if (bl->node[ii][jj][kk]==bl->edge[nedge]->node[1]) l*=-1.0;
          }
          else {
            if (bl->node[ii][jj][kk]==bl->edge[nedge]->node[0]) l*=-1.0;
          }

          for (idim=0;idim<3;idim++) bl->edge[nedge]->cutData->norm[idim]/=l;
        }
        else {

          //the cut shpild not be here

          //check if the cut is setted up
          if (bl->edge[nedge]->cutData!=NULL) {
            double *n=bl->edge[nedge]->cutData->norm;

  //################  DEBUG #################

  if (bl->edge[nedge]->cutPoint->TEMP_IDENTIFYER_NODE==71215) {
  cout << __LINE__ << endl;
  }

  //################ END DEBUG ##############



            if (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]>0.00001) n[0]=0.0,n[1]=0.0,n[2]=0.0;
          }

        }
      }

  }
  //################ END: EXPLIVCITLY CHECK THE CORNER AND CUT NODES ################



      //set up the cut nodes on edges: a cut node can exists when only one corner node of the edge exists
      for (nedge=0;nedge<12;nedge++) {
        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        globalNodeMask[i][j][k]=-1;


        //the cut node exists if
        //1. bl->edge[nedge]->cutData!=NULL
        //2. the length of bl->edge[nedge]->cutData->norm > 0.0
        //3. only one of nodes that defines the edge if the block exists

        if (bl->edge[nedge]->cutData==NULL) continue;

        double *cutNorm=bl->edge[nedge]->cutData->norm;
        if (cutNorm[0]*cutNorm[0]+cutNorm[1]*cutNorm[1]+cutNorm[2]*cutNorm[2]<0.01) continue;



        //check existance of the cornet nodes
        register int pnode;
        int ii,jj,kk;

        int nCornerNodes=0;

        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

          if (globalNodeMask[ii][jj][kk]>=0) nCornerNodes++;
        }







        //add the cut node to the node mesh
        if (bl->edge[nedge]->cutPoint!=NULL) if ( (cutBlockFlag==false) || (nCornerNodes==1)) {   ///////((nCornerNodes==1)&&(bl->edge[nedge]->cutPoint!=NULL)) {
          if (blNodeListLength==blNodeListLengthMax) {
            exit(__LINE__,"Error: the nodeList is too short");
          }

          globalNodeMask[i][j][k]=1;

          blNodeList[blNodeListLength].assign(bl->edge[nedge]->cutPoint);
          nodeMesh[i][j][k]=blNodeList+blNodeListLength;
          nodeMesh[i][j][k]->ghostNode=false;
          blNodeListLength++;

          if (bl->edge[nedge]->cutPoint->nodeno==-1) {
            bl->edge[nedge]->cutPoint->nodeno=nodenoCounter++;
            elementPosition=nodes.GetEntryCountingNumber(bl->edge[nedge]->cutPoint);
            fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

            elementPosition=blocks.GetEntryCountingNumber(bl);
            fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
          }
        }
      }

      //add the corner nodes to the nodes mesh




      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) { ///if (globalNodeMask[i][j][k]>=0) {
        if (blNodeListLength==blNodeListLengthMax) {
          exit(__LINE__,"Error: the nodeList is too short");
        }




  if (globalNodeMask[i][j][k]==0) {
    globalNodeMask[i][j][k]=1;

    for (long int nSurfaceElement=0;nSurfaceElement<nSurfaceNASTRANelements;nSurfaceElement++) {
      if (surfaceNASTRANmesh[nSurfaceElement].checkPointInsideDomain(bl->node[i][j][k]->x,EPS)==false) {
        globalNodeMask[i][j][k]=-1;
        break;
      }
    }
  }



  //////globalNodeMask[i][j][k]=(surfaceNASTRANmesh.checkPointInsideDomain(bl->node[i][j][k]->x)==true) ? 1 : -1;







        blNodeList[blNodeListLength].assign(bl->node[i][j][k]);
        nodeMesh[i][j][k]=blNodeList+blNodeListLength;


  nodeMesh[i][j][k]->ghostNode=(globalNodeMask[i][j][k]>=0) ? false : true;

        blNodeListLength++;



        if ( (nodeMesh[i][j][k]->ghostNode==false) || (cutBlockFlag==false) ) if (bl->node[i][j][k]->nodeno==-1) {
          bl->node[i][j][k]->nodeno=nodenoCounter++;
          elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
          fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

          elementPosition=blocks.GetEntryCountingNumber(bl);
          fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
        }

      }






  // END OF NEW PROCEDURE FOR POPULATING THE BLOCK'S NODE MESH

      //connect nodes:nodes in corners of a block can be connected only along coordinate directions; nodes that cut edges sould be connected to another cutting node on the same face
      cBlockNode *connectNode,*nd;
      int n,ndir;
      bool foundflag;


      //connect the corner nodes
      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) if (nodeMesh[i][j][k]!=NULL) {
        nd=nodeMesh[i][j][k];

        for (ndir=0;ndir<3;ndir++) {
          connectNode=NULL;

          switch (ndir) {
          case 0:
            //check connection in the i-direction
            connectNode=(nodeMesh[1][j][k]==NULL) ? nodeMesh[(i==0) ? 2 : 0][j][k] : nodeMesh[1][j][k];
            break;
          case 1:
            //check connection in the j-direction
            connectNode=(nodeMesh[i][1][k]==NULL) ? nodeMesh[i][(j==0) ? 2 : 0][k] : nodeMesh[i][1][k];
            break;
          case 2:
            //check connection in the k-direction
            connectNode=(nodeMesh[i][j][1]==NULL) ? nodeMesh[i][j][(k==0) ? 2 : 0] : nodeMesh[i][j][1];
          }

          if (connectNode!=NULL) nd->Connect(connectNode);
        }
      }




      //connect cutting nodes: only those cutting nodes are connected that belongs to the same face
      cBlockNode *nd0,*nd1,*nd2,*nd3,*ndcnt;
      int nCuts,nface,pedge;







      //check connection of nodes that cuts the edges
      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) if ((i==1)||(j==1)||(k==1)) if (nodeMesh[i][j][k]!=NULL) if (nodeMesh[i][j][k]->nConnections==0) {
        exit(__LINE__,"Error: a cutting node has no connections");
      }


      //copy the cut nodes on the block's node map
      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) if ((i==1)||(j==1)||(k==1)) if (nodeMesh[i][j][k]!=NULL) bl->node[i][j][k]=nodeMesh[i][j][k]->meshNode;




      //combine the nodes into thetrahedrones
      int nConnections,ncnt;
      int nnodesNumber=blNodeListLength;

      double *xt0ptr,xt1[3],*xt1ptr,xt2[3],*xt2ptr,xt3[3],*xt3ptr,tVolume;
      cBlockNode *tn0;

      //reset the volume of the cell associated with the block
      #if _USE_CELLS_ == _USE_CELLS_ON_
      bl->cell->Volume=0.0;
      #endif

      while (nnodesNumber>3) {
        //choose the first point that has 3 connections or a point that has a maximum number of connections


  //choose the gost nodes

      for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if ((blNodeList[n].ghostNode==true)&&(blNodeList[n].nConnections==3)) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }


              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-25*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;

          }

        }


        if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].nConnections==3) {
          //check the volume of the prospective a tetrahedron
          tn0=blNodeList+n;

          xt0ptr=tn0->meshNode->x;
          xt1ptr=tn0->Connection[0]->meshNode->x;
          xt2ptr=tn0->Connection[1]->meshNode->x;
          xt3ptr=tn0->Connection[2]->meshNode->x;

          if ((tn0->Connection[0]->ghostNode==true)||(tn0->Connection[1]->ghostNode==true)||(tn0->Connection[2]->ghostNode==true)) continue;


          for (idim=0;idim<3;idim++) {
            xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
            xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
            xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
          }

          tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
          if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
            nd0=blNodeList+n;
            break;
          }
        }



      if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].ghostNode==true) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }


              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-25*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;

          }

        }




        if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].nConnections==3) {
          //check the volume of the prospective a tetrahedron
          tn0=blNodeList+n;

          xt0ptr=tn0->meshNode->x;
          xt1ptr=tn0->Connection[0]->meshNode->x;
          xt2ptr=tn0->Connection[1]->meshNode->x;
          xt3ptr=tn0->Connection[2]->meshNode->x;

          for (idim=0;idim<3;idim++) {
            xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
            xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
            xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
          }

          tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
          if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
            nd0=blNodeList+n;
            break;
          }
        }


         if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }

              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;
            }
          }


        if (nd0==NULL) {
          return connectivityLength;

          exit(__LINE__,"Error: all possible corner nodes has a zero volume tetrahedron among their connections");
        }

        if (nd0->nConnections==3) {
          //if the node nd0 has only 3 connections, combine them into a thetrahedral
          nd1=nd0->Connection[0];
          nd2=nd0->Connection[1];
          nd3=nd0->Connection[2];

          nd1->Connect(nd2);
          nd1->Connect(nd3);
          nd2->Connect(nd3);

          #if _USE_CELLS_ == _USE_CELLS_ON_
          bl->cell->Volume+=tVolume;
          #endif

          //print the connectivity list
        if (true) {   //////((FileOutputFlag==true)&&(printConnectivityList==true)) {

  /////
   //NBLOCKS++;
   BLOCKPRINTED++;


  if (PRINTEDFLAG==false) PRINTEDFLAG=true,PRINTEDCELLS++;
  ///





  if (  (cutBlockFlag==false) ||  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false))) {

    blocknoCounter++;

          if (nd1->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd1->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd2->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd2->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd3->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd3->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd0->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd0->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);






  ////////

  if ((1+nd0->meshNode->nodeno==0)||(1+nd1->meshNode->nodeno==0)||(1+nd2->meshNode->nodeno==0)||(1+nd3->meshNode->nodeno==0)) {
  cout << __LINE__ << endl;
  }
  //////


  }



            //the surface mesh should contain faces that:
            //1. are a part of a thetrahedrons that contain a ghost node (<-> the face is adjesent to the boundary)
            //2. the face itself should not contain any ghost nodes
            //3. all nodes that compose the cut face should have an attribute that is different from _DEFAULT_NODEAT_


            //prepare data for the surface mesh
            cAMRnode<cNode> *thetraNodes[4];
            int i,j,k;


            bool flag;
            cBlockNode *ndlist[4];



            //the thetrahedral face is a surface face only if each node of the face has the nodeat different from _DEFAULT_NODEAT_
            //the first thetrahedron
            thetraNodes[0]=nd1->meshNode;
            thetraNodes[1]=nd2->meshNode;
            thetraNodes[2]=nd3->meshNode;
            thetraNodes[3]=nd0->meshNode;

            ndlist[0]=nd1;
            ndlist[1]=nd2;
            ndlist[2]=nd3;
            ndlist[3]=nd0;




  if (cutBlockFlag==true) {
  //the block is cut


            for (i=0;i<4;i++)  {
              for (flag=false,j=0;j<4;j++) if (i!=j) if (ndlist[j]->ghostNode==true) {
                flag=true;
                break;
              }

              if (flag==true) continue;


              //if the node ndlist[i] is 'ghostNode' than the face 'ndlist[j]' belongs to the surface
              if (ndlist[i]->ghostNode==true) flag=true;



               //if 'flag==true' -> the surface face is found

               if (flag==true) {

                   cTriangularSurfaceFace newface;
                   cAMRnode<cNode> *upNode;


                   for (k=0,j=0;j<4;j++) {
                     if (i!=j) {
                       newface.node[k++]=thetraNodes[j];
                       if (thetraNodes[j]->nodeat!=_DEFAULT_NODEAT_) newface.faceat=thetraNodes[j]->nodeat;
                     }
                     else upNode=thetraNodes[j];
                   }

                   //get the external normal and the surface area of the cut face
                   double e0[3],e1[3],l,*x0,*x1,*x2,norm[3],c;

                   x0=newface.node[0]->x;
                   x1=newface.node[1]->x;
                   x2=newface.node[2]->x;

                   for (idim=0;idim<3;idim++) e0[idim]=x1[idim]-x0[idim],e1[idim]=x2[idim]-x0[idim];



                   norm[0]=e0[1]*e1[2]-e1[1]*e0[2];
                   norm[1]=e1[0]*e0[2]-e0[0]*e1[2];
                   norm[2]=e0[0]*e1[1]-e1[0]*e0[1];

                   l=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
                   newface.surfaceArea=0.5*l;

                   for (idim=0,c=0.0;idim<3;idim++) c+=norm[idim]*(upNode->x[idim]-x0[idim]);

                   if (c>0.0) l*=-1.0;

                   for (idim=0,c=0.0;idim<3;idim++) newface.externalNormal[idim]=norm[idim]/l;



                   //add the cut surface to the surface mesh
                   newface.nextSurfaceCutFace=bl->firstSurfaceCutFace;
                   newface.block=bl;
                   bl->firstSurfaceCutFace=SurfaceMeshData.nfaces;

                   SurfaceMeshData.faces.push_back(newface);
                   SurfaceMeshData.nfaces++;
                }
            }

  }
  else {
    //the block is not cut

  exit(__LINE__,"ERROR: the block is not cut - non finished section");


    //the surface face: 1. must contan a node that is an cutting node of an edge
    //2. must contain a node that is belongs to the computational domain

    int nface,pedge,nedge,nCorner,i[3],searchDirection,cnt,nd;
    cAMRnode<cNode> *firstNode,*secondNode,*thirdNode,*forthNode;

    static const int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};




  }



          }

          //inicrement the connectivity counter
  if ( (cutBlockFlag==false) || ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) )      ++connectivityLength;

          //disconnect the point nd0
          nd1->removeConnection(nd0);
          nd2->removeConnection(nd0);
          nd3->removeConnection(nd0);

          //udpate the number of left nodes
          if (nd0->nConnections==0) nnodesNumber--;
          if (nd1->nConnections==0) nnodesNumber--;
          if (nd2->nConnections==0) nnodesNumber--;
          if (nd3->nConnections==0) nnodesNumber--;
        }
        else if (nd0->nConnections>3) {
          //create a plane and sort the nodes in the plane
          //the point of the origin of the plane, normal to the plane aand coordinate vectors wrelated to the plane
          double x0[3]={0.0,0.0,0.0},norm[3],e0[3],e1[3];
          double *x,length;
          cBlockNode **nd0Connection=nd0->Connection;

          for (n=0;n<nConnections;n++) for (x=nd0Connection[n]->meshNode->x,idim=0;idim<3;idim++) x0[idim]+=x[idim];

          for (idim=0,length=0.0;idim<3;idim++) {
            x0[idim]/=nConnections;

            e0[idim]=nd0Connection[1]->meshNode->x[idim]-nd0Connection[0]->meshNode->x[idim];
            e1[idim]=nd0Connection[2]->meshNode->x[idim]-nd0Connection[0]->meshNode->x[idim];
          }


          norm[0]=e1[1]*e0[2]-e1[2]*e0[1];
          norm[1]=e1[2]*e0[0]-e1[0]*e0[2];
          norm[2]=e1[0]*e0[1]-e1[1]*e0[0];

          length=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);




          for (idim=0;idim<3;idim++) norm[idim]/=length;

          //construct a coordinate system in the plane: get the first coordinate vector
          if (fabs(norm[0])>1.0E-1) {
            register double t=sqrt(norm[0]*norm[0]+norm[1]*norm[1]);

            e0[0]=norm[1]/t,e0[1]=-norm[0]/t,e0[2]=0.0;
          }
          else {
            register double t=sqrt(norm[2]*norm[2]+norm[1]*norm[1]);

            e0[0]=0.0,e0[1]=-norm[2]/t,e0[2]=norm[1]/t;
          }

          //get the second coordinate vector
          e1[0]=norm[1]*e0[2]-norm[2]*e0[1];
          e1[1]=norm[2]*e0[0]-norm[0]*e0[2];
          e1[2]=norm[0]*e0[1]-norm[1]*e0[0];

          //get the angles of projections of the nodes in the plane (e0,e1,norm)
          double phi,xplane[2];
          double *phi_nd0Connection=new double[nConnections];

          for (n=0;n<nConnections;n++) {
            x=nd0Connection[n]->meshNode->x;

            for (idim=0,xplane[0]=0.0,xplane[1]=0.0;idim<3;idim++) xplane[0]+=(x[idim]-x0[idim])*e0[idim],xplane[1]+=(x[idim]-x0[idim])*e1[idim];

            phi_nd0Connection[n]=acos(xplane[0]/sqrt(xplane[0]*xplane[0]+xplane[1]*xplane[1]));
            if (xplane[1]<0.0) phi_nd0Connection[n]=2.0*Pi-phi_nd0Connection[n];
          }

          //sort the nodes with an increase of the angle
          int nmin,n1;
          double minphi;

          for (n=0;n<nConnections;n++) {
            for (n1=n+1,nmin=n,minphi=phi_nd0Connection[n];n1<nConnections;n1++) if (minphi>phi_nd0Connection[n1]) nmin=n1,minphi=phi_nd0Connection[n1];

            if (nmin!=n) {
              //swap the nodes in the list
              register cBlockNode *t;

              t=nd0Connection[n];
              nd0Connection[n]=nd0Connection[nmin];
              nd0Connection[nmin]=t;

              phi_nd0Connection[nmin]=phi_nd0Connection[n];
            }
          }

          delete [] phi_nd0Connection;

          //construct the set of thetrahedrons
          int nthetra;

          for (nd1=nd0Connection[0],nthetra=0;nthetra<nConnections-2;nthetra++) {
            nd2=nd0Connection[nthetra+1];
            nd3=nd0Connection[nthetra+2];

            nd1->Connect(nd2);
            nd1->Connect(nd3);
            nd2->Connect(nd3);


            //increment the volume of the cell assiciated with the block
            #if _USE_CELLS_ == _USE_CELLS_ON_
            xt0ptr=nd0->meshNode->x;
            xt1ptr=nd1->meshNode->x;
            xt2ptr=nd2->meshNode->x;
            xt3ptr=nd3->meshNode->x;

            for (idim=0;idim<3;idim++) {
              xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
              xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
              xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
            }

            tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
            bl->cell->Volume+=tVolume;
            #endif



            //print the connectivity list
        if (true) {  /////((FileOutputFlag==true)&&(printConnectivityList==true)) {

  if (nd1->meshNode->nodeno==-1) {
  cout << __LINE__ << endl;
  }


  if (  (cutBlockFlag==false) ||  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) ) {


    blocknoCounter++;

          if (nd1->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd1->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd2->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd2->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd3->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd3->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd0->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd0->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);






  ////////

  if ((1+nd0->meshNode->nodeno==0)||(1+nd1->meshNode->nodeno==0)||(1+nd2->meshNode->nodeno==0)||(1+nd3->meshNode->nodeno==0)) {
  cout << __LINE__ << endl;
  }
  //////



  }

              //prepare data for the surface mesh generation
              cAMRnode<cNode> *thetraNodes[4];
              int i,j,k;

            bool flag;
            cBlockNode *ndlist[4];

              //the thetrahedral face is a surface face only if each node of the face has the nodeat different from _DEFAULT_NODEAT_
              //the first thetrahedron
              thetraNodes[0]=nd1->meshNode;
              thetraNodes[1]=nd2->meshNode;
              thetraNodes[2]=nd3->meshNode;
              thetraNodes[3]=nd0->meshNode;

            ndlist[0]=nd1;
            ndlist[1]=nd2;
            ndlist[2]=nd3;
            ndlist[3]=nd0;


  if (cutBlockFlag==true) {
  //the block is cut


            for (i=0;i<4;i++)  {
              for (flag=false,j=0;j<4;j++) if (i!=j) if (ndlist[j]->ghostNode==true) {
                flag=true;
                break;
              }

              if (flag==true) continue;


              //if the node ndlist[i] is 'ghostNode' than the face 'ndlist[j]' belongs to the surface
              if (ndlist[i]->ghostNode==true) flag=true;


               if (flag==true) {

                   cTriangularSurfaceFace newface;
                   cAMRnode<cNode> *upNode;


                   for (k=0,j=0;j<4;j++) {
                     if (i!=j) {
                       newface.node[k++]=thetraNodes[j];
                       if (thetraNodes[j]->nodeat!=_DEFAULT_NODEAT_) newface.faceat=thetraNodes[j]->nodeat;
                     }
                     else upNode=thetraNodes[j];
                   }


                   //get the external normal and the surface area of the cut face
                   double e0[3],e1[3],l,*x0,*x1,*x2,norm[3],c;

                   x0=newface.node[0]->x;
                   x1=newface.node[1]->x;
                   x2=newface.node[2]->x;

                   for (idim=0;idim<3;idim++) e0[idim]=x1[idim]-x0[idim],e1[idim]=x2[idim]-x0[idim];

                   norm[0]=e0[1]*e1[2]-e1[1]*e0[2];
                   norm[1]=e1[0]*e0[2]-e0[0]*e1[2];
                   norm[2]=e0[0]*e1[1]-e1[0]*e0[1];


                   l=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
                   newface.surfaceArea=0.5*l;

                   for (idim=0,c=0.0;idim<3;idim++) c+=norm[idim]*(upNode->x[idim]-x0[idim]);

                   if (c>0.0) l*=-1.0;

                   for (idim=0,c=0.0;idim<3;idim++) newface.externalNormal[idim]=norm[idim]/l;


                    //add the cut surface to the surface mesh
                    newface.nextSurfaceCutFace=bl->firstSurfaceCutFace;
                    newface.block=bl;
                    bl->firstSurfaceCutFace=SurfaceMeshData.nfaces;

                    SurfaceMeshData.faces.push_back(newface);
                    SurfaceMeshData.nfaces++;
                }


              }

  }
  else {
    //the block is not cut


  exit(__LINE__,"ERROR: the block is not cut - non finished section");

    //the surface face: 1. must contan a node that is an cutting node of an edge
    //2. must contain a node that is belongs to the computational domain

    int nface,pedge,nedge,nCorner,i[3],searchDirection,cnt,nd;
    cAMRnode<cNode> *firstNode,*secondNode,*thirdNode,*forthNode;

    static const int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};



  }





            }

            //inicrement the connectivity counter
  if ( (cutBlockFlag==false) || ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) )        ++connectivityLength;



          }

          //disconnect the point nd0
          for (n=0;n<nConnections;n++) {
            nd1=nd0->Connection[0];
            nd0->removeConnection(nd1);
            if (nd1->nConnections==0) nnodesNumber--;
          }

          if (nd0->nConnections==0) nnodesNumber--;
        }
        else {
         //if the node nd0 has less than tree connections -> error
         exit(__LINE__,"Error: the maximum number of connections is less that three");
       }

      }








      return connectivityLength;
     */    }


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
       cCutBlock bl;
       list<cTetrahedron> indomainConnectivityList,outdomainConnectivityList;
       list<cTriangleCutFace> TriangleCutConnectivity;
       typename list<cTetrahedron>::iterator itr;

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

       CutCell::cutBlockTetrahedronConnectivity<cCutBlockNode,cCutBlock,cCutData,cTetrahedron>(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);
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
