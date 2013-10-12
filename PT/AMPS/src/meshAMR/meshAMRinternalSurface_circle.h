//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal circle boundary


#ifndef _AMR_INTERNAL_SURFACE_CIRCLE_
#define _AMR_INTERNAL_SURFACE_CIRCLE_

#include "math.h"


#include "meshAMRdef.h"
#include "mpichannel.h"


class cInternalCircleData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalSphericalData_UserDefined
#endif
{
protected:
  double OriginPosition[3],Radius;

  static long int nPolarSurfaceElements;
  static double dPolarAngle;

public:

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
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

  cInternalCircleData ()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
  : cInternalSphericalData_UserDefined()
#endif
  {
    cleanDataBuffer();
  }

  void SetGeneralSurfaceMeshParameters(long int nPolarElements) {
    nPolarSurfaceElements=nPolarElements;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_
    dPolarAngle=Pi/nPolarSurfaceElements;
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    dPolarAngle=2.0*Pi/nPolarSurfaceElements;
#else
    exit(__LINE__,__FILE__,"Error: unknown option");
#endif
  }

  void SetCircleGeometricalParameters(double *x0,double r) {
     for (int idim=0;idim<2;idim++) OriginPosition[idim]=x0[idim];
     Radius=r;
  }

  void GetCircleGeometricalParameters(double* &x0,double &r) {
     x0=OriginPosition;
     r=Radius;
  }

  long int GetLocalSurfaceElementNumber(long int nPolarElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nPolarElement<0)||(nPolarElement>=nPolarSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    #endif

    return nPolarElement;
  }

  void GetSurfaceElementIndex(int &nPolarElement,int nSurfaceElement) {

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nSurfaceElement<0)||(nSurfaceElement>=nPolarSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nSurfaceElement' is out of range");
    #endif

    nPolarElement=nSurfaceElement;
  }

  long int GetTotalSurfaceElementsNumber() {return nPolarSurfaceElements;}


  double GetSurfaceElementArea(int nPolarElement) {
    double res=0.0;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_
    res=(cos(nPolarSurfaceElements*dPolarAngle)-cos((nPolarSurfaceElements+1)*dPolarAngle))*pow(Radius,2);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
    res=Radius*dPolarAngle;
#else
    exit(__LINE__,__FILE__,"Error: option is not found");
#endif

    return res;
  }

  void GetSurfaceElementProjectionIndex(double *x,long int &nPolarElement) {
    double r,r2,PolarAngle,xNormalized[2];
    int idim;

    for (r2=0.0,idim=0;idim<2;idim++) r2+=pow(x[idim]-OriginPosition[idim],2);
    for (r=sqrt(r2),idim=0;idim<2;idim++) xNormalized[idim]=(x[idim]-OriginPosition[idim])/r;

    PolarAngle=acos(xNormalized[0]);
    if (xNormalized[1]<0.0) PolarAngle=2.0*Pi-PolarAngle;

    nPolarElement=PolarAngle/dPolarAngle;
  }

  void GetSurfaceCoordinate(double *x,long int iPolarPoint) {
    double PolarAngle;

    PolarAngle=dPolarAngle*iPolarPoint;

    x[0]=Radius*cos(PolarAngle)+OriginPosition[0];
    x[1]=Radius*sin(PolarAngle)+OriginPosition[1];
  }



  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
    long int iPolar;
    FILE *fout=NULL;

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
      fprintf(fout,"VARIABLES=\"Polar Angle\"");
      if (PrintStateVectorFlag==true) {
        if (PrintVariableList==NULL) exit(__LINE__,__FILE__,"Error: PrintVariableList is not defined");
        PrintVariableList(fout);
      }

      fprintf(fout,"\n");
    }
    else pipe.openSend(0);

    //interpolate and print the state vector
    long int InterpolationList[nPolarSurfaceElements],InterpolationListLength=0;

    for (iPolar=0;iPolar<nPolarSurfaceElements+1;iPolar++) {
      if (ThisThread==0) fprintf(fout,"%e ",iPolar*dPolarAngle*180.0/Pi);

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iPolar==0) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(iPolar);

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nPolarSurfaceElements-1);
#endif
        }
        else if (iPolar==nPolarSurfaceElements) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(iPolar-1);

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0);
#endif
        }
        else {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(iPolar-1);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(iPolar);
        }


        if (PrintDataStateVector!=NULL) PrintDataStateVector(fout,iPolar,InterpolationList,InterpolationListLength,this,nDataSet,&pipe,ThisThread,nTotalThreads);
      }

      if (ThisThread==0) fprintf(fout,"\n");
    }

    //close the pipe
    if (ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();

    if (ThisThread==0) fclose(fout);
  }


  void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

  //=======================================================================
  //intersection of a block with the circle
  int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
    int idim,i,j,nCounter;
    double x[2];

    //if all corners of the block are within the sphere -> the block is entirely within the sphere
    for (nCounter=0,i=0;i<2;i++) {
      x[0]=((i==0) ? xBlockMin[0] : xBlockMax[0])-OriginPosition[0];

      for (j=0;j<2;j++) {
        x[1]=((j==0) ? xBlockMin[1] : xBlockMax[1])-OriginPosition[1];

        if (x[0]*x[0]+x[1]*x[1]<=pow(Radius+EPS,2)) nCounter++;
        else goto Not_Inside_Circle;
      }
    }


Not_Inside_Circle:
    if (nCounter==4) return _AMR_BLOCK_OUTSIDE_DOMAIN_;

    //check if the sphere intersects the block
    double dx[2],x0[2],e0[2],a,b,c,d,sqrt_d,t1,t2,tEPS;
    int nface;

    //the internal coordinated of the origin of the coordinate frame related to a face
    static const int nX0face[4][2]={{0,0},{1,0}, {0,0},{0,1}};

    //the internal coordinate of the nodes that determine the coordinate vectors related to the frame
    static const int nE0[4][2]={{0,1},{1,1}, {1,0},{1,1}};


    for (nface=0;nface<4;nface++) {
      a=0.0,b=0.0,c=0.0;

      for (idim=0;idim<2;idim++) {
        x0[idim]=(nX0face[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim];
        e0[idim]=((nE0[nface][idim]==0) ? xBlockMin[idim] : xBlockMax[idim]) - x0[idim];

        dx[idim]=x0[idim]-OriginPosition[idim];
        a+=pow(e0[idim],2);
        b+=2*e0[idim]*dx[idim];
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

      tEPS=EPS/sqrt(e0[0]*e0[0]+e0[1]*e0[1]);

      if (((-tEPS<t1)&&(t1<1.0+tEPS)) || ((-tEPS<t2)&&(t2<1.0+tEPS))) return _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
    }

    return _AMR_BLOCK_INSIDE_DOMAIN_;
  }


  double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
    int BlockIntersectionCode;
    double res=0.0,TotalResult=0.0;
    int idim;

    static const int nLevelMax=6;

    struct cLevelData {
      double xSubBlockMin[2],xSubBlockMax[2],dx;
      int i;
    };

    cLevelData LevelData[nLevelMax+1];
    cLevelData *levelDataPtr;
    int nLevel=0;


    exit(__LINE__,__FILE__,"error: the function is not finished");


    //init the level 0 data
    for (levelDataPtr=LevelData,idim=0;idim<2;idim++) levelDataPtr->xSubBlockMax[idim]=xBlockMaxInit[idim],levelDataPtr->xSubBlockMin[idim]=xBlockMinInit[idim];


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
      for (res=1.0,idim=0;idim<2;idim++) res*=levelDataPtr->xSubBlockMax[idim]-levelDataPtr->xSubBlockMin[idim];

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

      (levelDataPtr+1)->xSubBlockMin[1]=levelDataPtr->xSubBlockMin[1];
      (levelDataPtr+1)->xSubBlockMax[1]=levelDataPtr->xSubBlockMax[1];

      for (levelDataPtr->i=0;levelDataPtr->i<2;levelDataPtr->i++) {
        (levelDataPtr+1)->xSubBlockMin[0]=levelDataPtr->xSubBlockMin[0]+levelDataPtr->i*levelDataPtr->dx;
        (levelDataPtr+1)->xSubBlockMax[0]=levelDataPtr->xSubBlockMin[0]+(levelDataPtr->i+1)*levelDataPtr->dx;

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
    else {
      double SegmentSplittingTime[4],t1,t2,x0[2],l,x2,R2;
      double a,b,c,d,sqrt_d;
      int i,nSegments=1;

      double dii=0.5*(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0]);
      int ii;
      bool IntersectionFound=false;

      l=levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1];

      R2=Radius*Radius;

      for (ii=0;ii<2;ii++) {
        x0[0]=levelDataPtr->xSubBlockMin[0]+(ii+0.5)*dii;
        x0[1]=levelDataPtr->xSubBlockMin[1];

        a=l*l;
        b=2.0*l*(x0[1]-OriginPosition[1]);
        for (c=-R2,idim=0;idim<2;idim++) c+=pow(x0[idim]-OriginPosition[idim],2);

        d=b*b-4.0*a*c;

        if (d<=0.0) {
          if (pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1]+l*0.5,2)>R2) res+=1;
        }
        else {
          if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

          sqrt_d=sqrt(d);
          t1=-(b+sqrt_d)/(2.0*a);
          t2=-2.0*c/(b+sqrt_d);
          nSegments=0;

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
          x2=pow(x0[0]-OriginPosition[0],2);

          for (i=0;i<nSegments;i++) {
            t1=SegmentSplittingTime[i],t2=SegmentSplittingTime[i+1];

            if (x2+pow(x0[1]-OriginPosition[1]+l*0.5*(t1+t2),2)>R2) {
              IntersectionFound=true;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
              res+=dii*l*(t2-t1);
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_
              res+=Pi*dii*(pow(l*t2+levelDataPtr->xSubBlockMin[1],2)-pow(l*t1+levelDataPtr->xSubBlockMin[1],2));
#else
              exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif
            }
          }

          if (nSegments==0) if (x2+pow(x0[1]-OriginPosition[1]+l*0.5,2)>R2) {
            IntersectionFound=true;

#if _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_
            res+=dii*l;
#elif _AMR_SYMMETRY_MODE_ == _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_
            res+=Pi*dii*(pow(levelDataPtr->xSubBlockMax[1],2)-pow(levelDataPtr->xSubBlockMin[1],2));
#else
            exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif
          }
        }
      }

      if (IntersectionFound==false) res=1.0E-20;

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


};

#endif
