//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal surfaces installed into the mesh


#ifndef _AMR_INTERNAL_SURFACE_SPHERE_
#define _AMR_INTERNAL_SURFACE_SPHERE_

#include "math.h"


#include "meshAMRdef.h"
#include "mpichannel.h"
#include "rnd.h"


//=======================================================================
//the class describes the data that defines the spherical internal boundary; the class may contains the user defined data
class cInternalSphericalData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalSphericalData_UserDefined
#endif
{
//protected:
public:
  double OriginPosition[3],Radius;
//  void *SurfaceData;

  static long int nZenithSurfaceElements,nAzimuthalSurfaceElements;
  static double dAzimuthalAngle;

  #if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
  static double dCosZenithAngle;
  #elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
  static double dZenithAngle;
  #endif

//public:

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
  fPrintDataStateVector PrintDataStateVector;

  typedef double (*fLocalResolution)(double *);
  fLocalResolution localResolution;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif


  void cleanDataBuffer() {
    OriginPosition[0]=0.0,OriginPosition[1]=0.0,OriginPosition[2]=0.0;
    Radius=0.0;

//    ,SurfaceData=NULL;

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL,localResolution=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalSphericalData ()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
  : cInternalSphericalData_UserDefined()
#endif
  {
    cleanDataBuffer();
    SetGeneralSurfaceMeshParameters(nZenithSurfaceElements,nAzimuthalSurfaceElements);
  }



  static void SetGeneralSurfaceMeshParameters(long int nZenithElements,long int nAzimuthalElements) {
    nZenithSurfaceElements=nZenithElements,nAzimuthalSurfaceElements=nAzimuthalElements;

    dAzimuthalAngle=2.0*Pi/nAzimuthalSurfaceElements;

    #if  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    dCosZenithAngle=2.0/nZenithSurfaceElements;
    #elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    dZenithAngle=Pi/nZenithSurfaceElements;
    #endif
  }


  /*
  void SetSurfaceDataPointer(void *DataPointer) {SurfaceData=DataPointer;}
  void *GetSurfaceDataPointer() {return SurfaceData;}
*/

  void SetSphereGeometricalParameters(double *x0,double r) {
     for (int idim=0;idim<3;idim++) OriginPosition[idim]=x0[idim];
     Radius=r;
  }

  void GetSphereGeometricalParameters(double* &x0,double &r) {
     x0=OriginPosition;
     r=Radius;
  }

  inline long int GetLocalSurfaceElementNumber(long int nZenithElement,long int nAzimuthalElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    #endif

    return nZenithElement+nZenithSurfaceElements*nAzimuthalElement;
  }

  inline void GetSurfaceElementIndex(int &nZenithElement,int &nAzimuthalElement,int nSurfaceElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nSurfaceElement<0)||(nSurfaceElement>=nZenithSurfaceElements*nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nSurfaceElement' is out of range");
    #endif

    nAzimuthalElement=(int)(nSurfaceElement/nZenithSurfaceElements); 
    nZenithElement=nSurfaceElement-nZenithSurfaceElements*nAzimuthalElement; 
  }

  long int GetTotalSurfaceElementsNumber() {return nZenithSurfaceElements*nAzimuthalSurfaceElements;}


  double GetSurfaceElementArea(int nZenithElement,int nAzimuthalElement) {
    double cosZenithAngleMin,cosZenithAngleMax;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    #endif

#if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    cosZenithAngleMin=1.0-dCosZenithAngle*nZenithElement;
    cosZenithAngleMax=1.0-dCosZenithAngle*(1+nZenithElement);
#elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    cosZenithAngleMin=cos(nZenithElement*dZenithAngle);
    cosZenithAngleMax=cos((1+nZenithElement)*dZenithAngle);
#else
    exit(__LINE__,__FILE__,"Error: wrong option");
#endif

    return (cosZenithAngleMin-cosZenithAngleMax)*dAzimuthalAngle*pow(Radius,2);
  }

  double GetSurfaceElementArea(int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    return GetSurfaceElementArea(nZenithElement,nAzimuthalElement);
  } 


  inline void GetSurfaceElementNormal(double *norm,int nZenithElement,int nAzimuthalElement) {
    double ZenithAngle,AzimuthalAngle;

    AzimuthalAngle=(nAzimuthalElement+0.5)*dAzimuthalAngle;

#if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    ZenithAngle=acos(1.0-dCosZenithAngle*(nZenithElement+0.5));
#elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    ZenithAngle=(0.5+nZenithElement)*dZenithAngle
#else
    exit(__LINE__,__FILE__,"Error: wrong option");
#endif

    norm[0]=sin(ZenithAngle)*cos(AzimuthalAngle);
    norm[1]=sin(ZenithAngle)*sin(AzimuthalAngle);
    norm[2]=cos(ZenithAngle);
  }

  inline void GetSurfaceElementNormal(double *norm,int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementNormal(norm,nZenithElement,nAzimuthalElement);
  }



  inline void GetSurfaceElementMiddlePoint(double *x,int nZenithElement,int nAzimuthalElement) {
    int idim;

    GetSurfaceElementNormal(x,nZenithElement,nAzimuthalElement);
    for (idim=0;idim<3;idim++) x[idim]=x[idim]*Radius+OriginPosition[idim];
  }


  inline void GetSurfaceElementMiddlePoint(double *x,int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementMiddlePoint(x,nZenithElement,nAzimuthalElement);
  }

  inline void GetSurfaceElementRandomDirection(double *x,int nZenithElement,int nAzimuthalElement) {
    double ZenithAngle,AzimuthalAngle;
    double c0,c1;


    AzimuthalAngle=(nAzimuthalElement+rnd())*dAzimuthalAngle;

#if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    c0=1.0-dCosZenithAngle*nZenithElement;
    c1=1.0-dCosZenithAngle*(nZenithElement+1);
#elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    c0=cos(nZenithElement*dZenithAngle);
    c1=cos((1+nZenithElement)*dZenithAngle);
#else
    exit(__LINE__,__FILE__,"Error: wrong option");
#endif

    ZenithAngle=acos(c0+rnd()*(c1-c0));

    x[0]=sin(ZenithAngle)*cos(AzimuthalAngle);
    x[1]=sin(ZenithAngle)*sin(AzimuthalAngle);
    x[2]=cos(ZenithAngle);
  }

  inline void GetSurfaceElementRandomPoint(double *x,int nZenithElement,int nAzimuthalElement) {
    GetSurfaceElementRandomDirection(x,nZenithElement,nAzimuthalElement);

    x[0]=Radius*x[0]+OriginPosition[0];
    x[1]=Radius*x[1]+OriginPosition[1];
    x[2]=Radius*x[2]+OriginPosition[2];
  }

  inline void GetSurfaceElementProjectionIndex(double *x,long int &nZenithElement,long int &nAzimuthalElement) {
    double r,r2,AzimuthalAngle,xNormalized[3];
    int idim;

    for (r2=0.0,idim=0;idim<3;idim++) r2+=pow(x[idim]-OriginPosition[idim],2);
    for (r=sqrt(r2),idim=0;idim<3;idim++) xNormalized[idim]=(x[idim]-OriginPosition[idim])/r;

    if ((r=pow(xNormalized[0],2)+pow(xNormalized[1],2))>0.0) {
      AzimuthalAngle=acos(xNormalized[0]/sqrt(r));
      if (xNormalized[1]<0.0) AzimuthalAngle=2.0*Pi-AzimuthalAngle;
    }
    else AzimuthalAngle=0.0;

    #if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    nZenithElement=(long int)((1.0-xNormalized[2])/dCosZenithAngle);
    #elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    double ZenithAngle=acos(xNormalized[2]);
    nZenithElement=(long int)(ZenithAngle/dZenithAngle);
    #else
    exit(__LINE__,__FILE__,"Error: wrong option");
    #endif


    nAzimuthalElement=(long int)(AzimuthalAngle/dAzimuthalAngle);

    if (nZenithElement==nZenithSurfaceElements) --nZenithElement;
    if (nAzimuthalElement==nAzimuthalSurfaceElements) --nAzimuthalElement;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    static unsigned long int FunctionCallCounter=0;

    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) {
      printf("$PREFIX:Error: out of range\n nZenithElement=%ld, nAzimuthalElement=%ld (%s@%i)\n",nZenithElement,nAzimuthalElement,__FILE__,__LINE__);
      printf("$PREFIX:x=%e, %e, %e\nCallCounter=%ld\n",x[0],x[1],x[2],FunctionCallCounter);
      exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    }

    FunctionCallCounter++;
    #endif
  }

  inline void GetSurfaceNormal(double *x,double iZenithPoint,double  iAzimutalPoint) {
    double cosZenithAngle,sinZenithAngle,AzimuthalAngle;

    #if _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
    cosZenithAngle=1.0-dCosZenithAngle*iZenithPoint;
    sinZenithAngle=sqrt(1.0-cosZenithAngle*cosZenithAngle);
    #elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
    double ZenithAngle=dZenithAngle*iZenithPoint;
    cosZenithAngle=cos(ZenithAngle);
    sinZenithAngle=sin(ZenithAngle);
    #else
    exit(__LINE__,__FILE__,"Error: wrong option");
    #endif

    AzimuthalAngle=dAzimuthalAngle*iAzimutalPoint;

    x[0]=sinZenithAngle*cos(AzimuthalAngle);
    x[1]=sinZenithAngle*sin(AzimuthalAngle);
    x[2]=cosZenithAngle;
  }

  inline void GetSurfaceCoordinate(double *x,double iZenithPoint,double  iAzimutalPoint) {
    GetSurfaceNormal(x,iZenithPoint,iAzimutalPoint);

    x[0]=Radius*x[0]+OriginPosition[0];
    x[1]=Radius*x[1]+OriginPosition[1];
    x[2]=Radius*x[2]+OriginPosition[2];
  }

  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
    long int iZenith,iAzimuthal;
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
      fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",(nZenithSurfaceElements+1)*nAzimuthalSurfaceElements,nZenithSurfaceElements*nAzimuthalSurfaceElements);
    }
    else pipe.openSend(0);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      GetSurfaceCoordinate(x,iZenith,iAzimuthal);
      if (ThisThread==0) fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iZenith==0) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else if (iZenith==nZenithSurfaceElements) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else {
          int iA,iZ,A[2],Z[2];

          Z[0]=iZenith-1,Z[1]=iZenith;

          A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1;
          A[1]=iAzimuthal;

          for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
        }


        PrintDataStateVector(fout,iZenith,iAzimuthal,InterpolationList,InterpolationListLength,this,nDataSet,&pipe,ThisThread,nTotalThreads);
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
      for (iZenith=0;iZenith<nZenithSurfaceElements;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
        iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
        iAzimuthalMin=iAzimuthal;

        nd0=1+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements;
        nd1=1+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements;
        nd2=1+iAzimuthalMax+(iZenith+1)*nAzimuthalSurfaceElements;
        nd3=1+iAzimuthalMin+(iZenith+1)*nAzimuthalSurfaceElements;


/*
        nd0=(iZenith!=0) ? 2+iAzimuthalMin+(iZenith-1)*nAzimuthalSurfaceElements : 1; //iZenithMin,iAzimuthalMin
        nd1=(iZenith!=nZenithSurfaceElements-1) ? 2+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements : 2+(nZenithSurfaceElements-1)*nAzimuthalSurfaceElements; //iZenithMax,iAzimuthalMin
        nd2=(iZenith!=nZenithSurfaceElements-1) ? 2+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements : 2+(nZenithSurfaceElements-1)*nAzimuthalSurfaceElements; //iZenithMax,iAzimuthalMax
        nd3=(iZenith!=0) ? 2+iAzimuthalMax+(iZenith-1)*nAzimuthalSurfaceElements : 1; //iZenithMin,iAzimuthalMax
*/

        fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
      }

      fclose(fout);
    }
  }


   void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

   //=======================================================================
   //intersection of a block with the sphere
   int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
     int idim,i,j,k,nCounter;
     double x[3];

     /*
     //if xBlockMax[]<OriginPosition-Radius or OriginPosition[]+Radius<xBlockMin -> the block does not cross the sphere
     for (idim=0;idim<3;idim++) if ((xBlockMax[idim]<OriginPosition[idim]-Radius)||(OriginPosition[idim]+Radius<xBlockMin[idim])) return _AMR_BLOCK_INSIDE_DOMAIN_;

*/

     //if all corners of the block are within the sphere -> the block is entirely within the sphere
     for (nCounter=0,i=0;i<2;i++) {
       x[0]=((i==0) ? xBlockMin[0] : xBlockMax[0])-OriginPosition[0];

       for (j=0;j<2;j++) {
         x[1]=((j==0) ? xBlockMin[1] : xBlockMax[1])-OriginPosition[1];

         for (k=0;k<2;k++) {
           x[2]=((k==0) ? xBlockMin[2] : xBlockMax[2])-OriginPosition[2];


//           cout << x[0]*x[0]+x[1]*x[1]+x[2]*x[2] << "   "  << pow(Radius+EPS,2) << endl;
           if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<=pow(Radius+EPS,2)) nCounter++;
           else goto Not_Inside_Sphere;

         }
       }
     }


Not_Inside_Sphere:
     if (nCounter==8) return _AMR_BLOCK_OUTSIDE_DOMAIN_;


     /*
     //check if the sphere is entirely within the block
     for (nCounter=0,idim=0;idim<3;idim++) {
       if ((OriginPosition[idim]+Radius<xBlockMax[idim]+EPS)&&(OriginPosition[idim]-Radius>xBlockMin[idim]-EPS)) nCounter++;
     }

     if (nCounter==3) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_; //the sphere is entirely inside the block


     //check if the block is intersected by the sphere
     double xSphereMin,xSphereMax;

     for (nCounter=0,idim=0;idim<3;idim++) {
       xSphereMin=OriginPosition[idim]-Radius,xSphereMax=OriginPosition[idim]+Radius;

       if ( ((xSphereMin<xBlockMin[idim]+EPS)&&(xSphereMax>xBlockMin[idim]-EPS)) || ((xSphereMin<xBlockMax[idim]+EPS)&&(xSphereMax>xBlockMax[idim]-EPS)) ) nCounter++;
     }

     if (nCounter==3) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
     */


     //check if the sphere intersects the block
     double r[3],x0[3],e0[3],e1[3],norm[3],rE0,rE1,lE0,lE1,lNorm,c;
     int nface;

     //the internal coordinated of the origin of the coordinate frame related to a face
     static const int nX0face[6][3]=  { {0,0,0},{1,0,0}, {0,0,0},{0,1,0}, {0,0,0},{0,0,1}};

     //the internal coordinate of the nodes that determine the coordinate vectors related to the frame
     static const int nE0[6][3]=  { {0,1,0},{1,1,0}, {1,0,0},{1,1,0}, {1,0,0},{1,0,1}};
     static const int nE1[6][3]=  { {0,0,1},{1,0,1}, {0,0,1},{0,1,1}, {0,1,0},{0,1,1}};

     //the direction to the local normal in the coordinate system related to the face
     static const int nNorm[6][3]={ {1,0,0},{0,0,0}, {0,1,0},{0,0,0}, {0,0,1},{0,0,0}};

     for (nface=0;nface<6;nface++) {
       for (idim=0,lNorm=0.0,lE0=0.0,lE1=0.0;idim<3;idim++) {
         x0[idim]=(nX0face[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];

         e0[idim]=((nE0[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0[idim];
         e1[idim]=((nE1[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0[idim];

         norm[idim]=((nNorm[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0[idim];

         lE0+=pow(e0[idim],2);
         lE1+=pow(e1[idim],2);
         lNorm+=pow(norm[idim],2);
       }

       for (c=0.0,idim=0,lNorm=sqrt(lNorm),lE1=sqrt(lE1),lE0=sqrt(lE0);idim<3;idim++) {
         r[idim]=OriginPosition[idim]-x0[idim];
         norm[idim]/=lNorm;

         c+=r[idim]*norm[idim];
       }

       //check the distance from the face to the point 'OriginPosition'
       if (fabs(c)<Radius+EPS) { //the distance betweenthe plane containing the face and the point 'OriginPosition' is less or equal 'Radius'
         for (rE0=0.0,rE1=0.0,idim=0;idim<3;idim++) {
           r[idim]-=c*norm[idim];

           e0[idim]/=lE0;
           e1[idim]/=lE1;

           rE0+=r[idim]*e0[idim];
           rE1+=r[idim]*e1[idim];
         }

         if ((-EPS<rE0)&&(rE0<lE0+EPS)&&(-EPS<rE1)&&(rE1<lE1+EPS)) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
       }

     }


     //check the intersection between the sphere and each edge of the block

     //internal nodes of the block that determine the edge
     static const int nX0Edge[12][3]={ {0,0,0},{0,1,0},{0,1,1},{0,0,1}, {0,0,0},{1,0,0},{1,0,1},{0,0,1}, {0,0,0},{1,0,0},{1,1,0},{0,1,0}};
     static const int nX1[12][3]={ {1,0,0},{1,1,0},{1,1,1},{1,0,1}, {0,1,0},{1,1,0},{1,1,1},{0,1,1}, {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

     int nedge;
     double l[3],dx[3],a=0.0,b=0.0,d,t1,t2,sqrt_d,tEPS;

     for (nedge=0;nedge<12;nedge++) {
       a=0.0,b=0.0,c=0.0;

       for (idim=0;idim<3;idim++) {
         x0[idim]=(nX0Edge[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
         l[idim]=((nX1[nedge][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0[idim];

         dx[idim]=x0[idim]-OriginPosition[idim];
         a+=pow(l[idim],2);
         b+=2*l[idim]*dx[idim];
         c+=pow(dx[idim],2);
       }

       c-=Radius*Radius;
       d=b*b-4.0*a*c;

       if (d<0.0) {
         if (4.0*a*pow(EPS,2)>-d) d=0.0; //equvalent to EPS/|l[:]| > sqrt(|d|)/(2a) -> the particle is within the distance of EPS from the surface's boundary
         else continue;
       }

       if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

       sqrt_d=sqrt(d);
       t1=-(b+sqrt_d)/(2.0*a);
       t2=-2.0*c/(b+sqrt_d);

       tEPS=EPS/sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);

       if (((-tEPS<t1)&&(t1<1.0+tEPS)) || ((-tEPS<t2)&&(t2<1.0+tEPS))) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
     }



     /*
     //check is the center of the sphere is outside of the block but a point of the sphere is inside the block and when the sphere is not intersecting the edges of the block
     for (i=0;i<3;i++) for (j=-1;j<2;j+=2) {
       x[0]=OriginPosition[0];
       x[1]=OriginPosition[1];
       x[2]=OriginPosition[2];

       x[i]+=j*Radius;

       //check if the point 'x' is in the block


     }
     */

     return _AMR_BLOCK_INSIDE_DOMAIN_;
   }


   double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
     int BlockIntersectionCode;
     double res=0.0,TotalResult=0.0;
     int idim;


//########  DEBUG ########
//return 0.0;
//########  END DEBUG #####




     static const int nLevelMax=6;

     /*
     int nLevelEPS=0;

     for (idim=0;idim<3;idim++) {
       register int t=1+(int)(log((xBlockMaxInit[idim]-xBlockMinInit[idim])/EPS)/log(2.0));
       if (t>nLevelEPS) nLevelEPS=t;
     }
*/


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
       double SegmentSplittingTime[4],t1,t2,x0[3],l,xy2,R2;
       double a,b,c,d,sqrt_d;
       int i,nSegments=1;

//########  DEBUG ########
//return 0.0;
//########  END DEBUG #####

       /*
       x0[0]=0.5*(levelDataPtr->xSubBlockMin[0]+levelDataPtr->xSubBlockMax[0]);
       x0[1]=0.5*(levelDataPtr->xSubBlockMin[1]+levelDataPtr->xSubBlockMax[1]);
       x0[2]=levelDataPtr->xSubBlockMin[2];
       */

       double dii=0.5*(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0]);
       double djj=0.5*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1]);
       int ii,jj;
       bool IntersectionFound=false;

       l=levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2];

       R2=Radius*Radius;

       for (ii=0;ii<3;ii++) for (jj=0;jj<3;jj++) {

         x0[0]=levelDataPtr->xSubBlockMin[0]+ii*dii;
         x0[1]=levelDataPtr->xSubBlockMin[1]+jj*djj;
         x0[2]=levelDataPtr->xSubBlockMin[2];


       a=l*l;
       b=2.0*l*(x0[2]-OriginPosition[2]);
       for (c=-R2,idim=0;idim<3;idim++) c+=pow(x0[idim]-OriginPosition[idim],2);

       d=b*b-4.0*a*c;

       if (d<=0.0) {
//         if (pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1],2)+pow(x0[2]-OriginPosition[2]+l*0.5,2)>R2) res+=1;

//         res=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])*(levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2]);
       }
       else {
         if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

         sqrt_d=sqrt(d);
         t1=-(b+sqrt_d)/(2.0*a);
         t2=-2.0*c/(b+sqrt_d);
         nSegments=0;

         SegmentSplittingTime[0]=0.0;
         if ((t1>0.0)&&(t1<1.0)) SegmentSplittingTime[++nSegments]=t1;

         if ((t2>0.0)&&(t2<1.0)) {
           if (t2>SegmentSplittingTime[nSegments]) SegmentSplittingTime[++nSegments]=t2;
           else {
             SegmentSplittingTime[2]=SegmentSplittingTime[1];
             SegmentSplittingTime[1]=t2;
             nSegments++;
           }
         }

         SegmentSplittingTime[++nSegments]=1.0;
         xy2=pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1],2);

         for (i=0;i<nSegments;i++) {
           t1=SegmentSplittingTime[i],t2=SegmentSplittingTime[i+1];

           if (xy2+pow(x0[2]-OriginPosition[2]+l*0.5*(t1+t2),2)>R2) {
             res+=t2-t1;
             IntersectionFound=true;
           }
         }

       }
       }

       res/=9.0;
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


/* reserve copy
   double GetRemainedBlockVolume(double *xBlockMin,double *xBlockMax,double EPS,int *IntersectionStatus,int nLevel=0) {
     int BlockIntersectionCode;
     double res=0.0;
     int idim;


//########  DEBUG ########
//return 0.0;
//########  END DEBUG #####




     static const int nLevelMax=7;


     BlockIntersectionCode=BlockIntersection(xBlockMin,xBlockMax,EPS);
     *IntersectionStatus=BlockIntersectionCode;

     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) return 0.0;
     if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
       for (res=1.0,idim=0;idim<3;idim++) res*=xBlockMax[idim]-xBlockMin[idim];
       return res;
     }


     //the block is intersected by the sphere
     if (nLevel<nLevelMax) {
       double xSubBlockMin[3],xSubBlockMax[3],dx,dy;
       int i,j;

       dx=(xBlockMax[0]-xBlockMin[0])/2.0;
       dy=(xBlockMax[1]-xBlockMin[1])/2.0;

       xSubBlockMin[2]=xBlockMin[2],xSubBlockMax[2]=xBlockMax[2];

       for (i=0;i<2;i++) {
         xSubBlockMin[0]=xBlockMin[0]+i*dx;
         xSubBlockMax[0]=xBlockMin[0]+(i+1)*dx;

         for (j=0;j<2;j++) {
           xSubBlockMin[1]=xBlockMin[1]+j*dy;
           xSubBlockMax[1]=xBlockMin[1]+(j+1)*dy;

           res+=GetRemainedBlockVolume(xSubBlockMin,xSubBlockMax,EPS,&BlockIntersectionCode,nLevel+1);
         }
       }
     }
     else {
       double SegmentSplittingTime[4],t1,t2,x0[3],l,xy2,R2;
       double a,b,c,d,sqrt_d;
       int i,nSegments=1;

//########  DEBUG ########
return 0.0;
//########  END DEBUG #####

       x0[0]=0.5*(xBlockMin[0]+xBlockMax[0]);
       x0[1]=0.5*(xBlockMin[1]+xBlockMax[1]);
       x0[2]=xBlockMin[2];

       l=xBlockMax[2]-xBlockMin[2];

       R2=Radius*Radius;

       a=l*l;
       b=2.0*l*(x0[2]-OriginPosition[2]);
       for (c=-R2,idim=0;idim<3;idim++) c+=pow(x0[idim]-OriginPosition[idim],2);

       d=b*b-4.0*a*c;

       if (d<=0.0) {
         res=(xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2]);
       }
       else {
         if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

         sqrt_d=sqrt(d);
         t1=-(b+sqrt_d)/(2.0*a);
         t2=-2.0*c/(b+sqrt_d);

         SegmentSplittingTime[0]=0.0;
         if ((t1>0.0)&&(t1<1.0)) SegmentSplittingTime[nSegments++]=t1;

         if ((t2>0.0)&&(t2<1.0)) {
           if (t2>SegmentSplittingTime[nSegments-1]) SegmentSplittingTime[nSegments++]=t2;
           else {
             SegmentSplittingTime[2]=SegmentSplittingTime[1];
             SegmentSplittingTime[1]=t2;
             nSegments++;
           }
         }

         SegmentSplittingTime[nSegments]=1.0;
         xy2=pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1],2);

         for (i=0;i<nSegments;i++) {
           t1=SegmentSplittingTime[i],t2=SegmentSplittingTime[i+1];

           if (xy2+pow(x0[2]+l*0.5*(t1+t2),2)>R2) res+=t2-t1;
         }

         res*=(xBlockMax[0]-xBlockMin[0])*(xBlockMax[1]-xBlockMin[1])*(xBlockMax[2]-xBlockMin[2]);
       }
     }


     return res;

   }
*/


};

//set up the surface parameters of the sphere
void inline SetGeneralSphereSurfaceMeshParameters(long int nZenithElements,long int nAzimuthalElements) {
  cInternalSphericalData::nAzimuthalSurfaceElements=nAzimuthalElements,cInternalSphericalData::nZenithSurfaceElements=nZenithElements;
  cInternalSphericalData::dAzimuthalAngle=2.0*Pi/cInternalSphericalData::nAzimuthalSurfaceElements;

#if  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
  cInternalSphericalData::dCosZenithAngle=2.0/cInternalSphericalData::nZenithSurfaceElements;
#elif _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ == _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_
  cInternalSphericalData::dZenithAngle=Pi/cInternalSphericalData::nZenithSurfaceElements;
#endif
}


#endif
