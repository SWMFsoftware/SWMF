/*
 * SphericalVolumeMesh.h
 *
 *  Created on: Oct 4, 2012
 *      Author: vtenishe
 */

//$Id$
//spherical volume mesh is used for sampling model data on a refined spherical mesh
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>
#include <dirent.h>

#include "constants.h"

#ifndef SPHERICALVOLUMEMESH_H_
#define SPHERICALVOLUMEMESH_H_


#define _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_       0
#define _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LOGARITHMIC_  1

#define _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_    0
#define _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__COSINE_    1


#define _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
#define _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_


class cSphericalVolumeMesh {
public:

  //radial distribution of points
  double minRadius,maxRadius,dR,dLogR;
  int nRadialPoints;

  //azimuthal and polar distributions of the mesh's pooints
  double dPolarAngle,dAzimuthAngle,dCosAzimuthAngle;
  int nPolarAnglePoints,nAzimuthAnglePoints;

  //sampling buffer
  int nSampledFieldVariables;
  double **SamplingBuffer;

  //definitions of the user-defined functions that print the variable list and output values of the field variables
  typedef int (*fPrintVariableList)(FILE*);
  typedef double (*fGetFieldVariable)(int,double*,int);

  fPrintVariableList funcPrintVariableList;
  fGetFieldVariable funcGetFieldVariable;

  cSphericalVolumeMesh() {
    minRadius=0.0,maxRadius=0.0,dR=0.0,dLogR=0.0;
    nRadialPoints=0;

    dPolarAngle=0.0,dAzimuthAngle=0.0,dCosAzimuthAngle=0.0;
    nPolarAnglePoints=0,nAzimuthAnglePoints=0;

    nSampledFieldVariables=0,SamplingBuffer=NULL;

    funcPrintVariableList=NULL,funcGetFieldVariable=NULL;
  }

  //flush the sampling buffer
  void FlushSamplingBuffer() {
    int var,el;
    const int eltotal=(nRadialPoints-1)*nPolarAnglePoints*(nAzimuthAnglePoints-1);

    for (var=0;var<nSampledFieldVariables;var++) for (el=0;el<eltotal;el++) SamplingBuffer[var][el]=0.0;
  }

  //init the mesh variables
  void init(double rmin,double rmax,int nRadialpoints,int nPolarpoints,int nAzimuthPoints,int nSampledVariables,fPrintVariableList UserDefinedPrintVariableList,fGetFieldVariable UserDefinedGetFieldVariable) {
    minRadius=rmin,maxRadius=rmax,nRadialPoints=nRadialpoints;
    nPolarAnglePoints=nPolarpoints,nAzimuthAnglePoints=nAzimuthPoints;

    funcPrintVariableList=UserDefinedPrintVariableList,funcGetFieldVariable=UserDefinedGetFieldVariable;

#if _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
    dR=(rmax-rmin)/(nRadialPoints-1);
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

    dPolarAngle=PiTimes2/nPolarAnglePoints;

#if _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_
    dAzimuthAngle=Pi/(nAzimuthAnglePoints-1);
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

    //init the sampling buffer;
    nSampledFieldVariables=nSampledVariables;
    SamplingBuffer=new double* [nSampledFieldVariables];

    for (int n=0;n<nSampledFieldVariables;n++) SamplingBuffer[n]=new double [(nRadialPoints-1)*nPolarAnglePoints*(nAzimuthAnglePoints-1)];
    FlushSamplingBuffer();
  }

  //calculate cartesian coordinate of a mesh point
  //the numbering order of the mesh points: radial point, poinar angle point, zenith angle point
  void getCartesianNodeCoordinate(int nradius,int npolar,int nzenith,double *x) {
    double r,phi,theta;

#if _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
    r=minRadius+dR*nradius;
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

#if _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_
   theta=dAzimuthAngle*nzenith;
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

   phi=dPolarAngle*npolar;

   x[0]=r*cos(phi)*sin(theta);
   x[1]=r*sin(phi)*sin(theta);
   x[2]=r*cos(theta);
  }

  void getCartesianCellCoordinate(int nradius,int npolar,int nzenith,double *x) {
    double r,phi,theta;

#if _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
    r=minRadius+dR*(nradius+0.5);
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

#if _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_
   theta=dAzimuthAngle*(nzenith+0.5);
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

   phi=dPolarAngle*(npolar+0.5);

   x[0]=r*cos(phi)*sin(theta);
   x[1]=r*sin(phi)*sin(theta);
   x[2]=r*cos(theta);
  }


  //convert indexces into the samplingelement number
  int GetSamplingElementNumber(int nRadial,int nPolar,int nAzimuth) {
    int el;

    el=nAzimuth+nPolar*(nAzimuthAnglePoints-1)+nRadial*(nPolarAnglePoints*(nAzimuthAnglePoints-1));
    if ((el<0)||(el>=(nRadialPoints-1)*nPolarAnglePoints*(nAzimuthAnglePoints-1))) exit(__LINE__,__FILE__,"Error: out of range");

    return el;
  }

  void GetSamplingElementIndexes(int& nRadial,int& nPolar,int& nAzimuth,int el) {
    nRadial=el/(nPolarAnglePoints*(nAzimuthAnglePoints-1));
    el-=(nPolarAnglePoints*(nAzimuthAnglePoints-1))*nRadial;

    nPolar=el/(nAzimuthAnglePoints-1);
    el-=(nAzimuthAnglePoints-1)*nPolar;

    nAzimuth=el;
  }

  int GetSamplingElementNumber(double* x) {
    double r,phi,theta;
    int nRadial,nPolar,nAzimuth;


    r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    if ((r<=minRadius)||(r>=maxRadius)) return -1;

    theta=acos(x[2]/r);

    phi=acos(x[0]/r/sin(theta));
    if (x[1]<0.0) phi=2.0*Pi-phi;

#if _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
    nRadial=(r-minRadius)/dR;
#else
    exit(__LINE__,__FILE__,"not implemented");
#endif

    nPolar=phi/dPolarAngle;

#if _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_
    nAzimuth=theta/dAzimuthAngle;
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif


    return GetSamplingElementNumber(nRadial,nPolar,nAzimuth);
  }

  double GetVolumeSamplingElement(int el) {
    int nRadial,nPolar,nAzimuth;
    double r,theta;
    double res=dPolarAngle;

    GetSamplingElementIndexes(nRadial,nPolar,nAzimuth,el);

#if _SPHERICAL_VOLUME_MESH_RADIAL_DISCTIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__RADIAL_DISTRIBUTION_MODE__LINEAR_
    r=minRadius+dR*nRadial;
    res*=(pow(r+dR,3)-pow(r,3))/3.0;
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

#if _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE_ == _SPHERICAL_VOLUME_MESH__AZIMUTHAL_DISTRIBUTION_MODE__LINEAR_
    theta=dAzimuthAngle*nAzimuth;
    res*=(cos(theta)-cos(theta+dAzimuthAngle));
#else
    exit(__LINE__,__FILE__,"Error: not implemented");
#endif

    return res;
  }

  //plot the mesh into a tecplot file
  void plot(const char *fname) {
    FILE *fout;
    int nRadial,nAzimuth,nPolar;


//    return;

    //move all sampled data to the root processor
    CMPI_channel pipe(1000000);
    int thread;
    int nvar,el,ntotalElements=(nRadialPoints-1)*nPolarAnglePoints*(nAzimuthAnglePoints-1);

    if (PIC::ThisThread!=0) {
      pipe.openSend(0);

      for (nvar=0;nvar<nSampledFieldVariables;nvar++) for (el=0;el<ntotalElements;el++) {
        double t;

        t=SamplingBuffer[nvar][el];

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
         if (!isfinite(t)) {
           exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
         }
#endif

        pipe.send(t);
      }

      pipe.closeSend();
      FlushSamplingBuffer();
      return;
    }
    else {
      pipe.openRecvAll();

      for (nvar=0;nvar<nSampledFieldVariables;nvar++) for (el=0;el<ntotalElements;el++) {
        double t;

        for (thread=1;thread<PIC::nTotalThreads;thread++) {
          pipe.recv(t,thread);

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
         if (!isfinite(t)) {
           exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
         }
#endif
          SamplingBuffer[nvar][el]+=t;
        }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
         if (!isfinite(SamplingBuffer[nvar][el])) {
           exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
         }
#endif
      }

      pipe.closeRecvAll();
    }

    //the number of field variables that will be printed (the varables are cell ceneterd)
    int nCellCenteredVariables=0;

    fout=fopen(fname,"w");

    //print the variable list
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");
    if (funcPrintVariableList!=NULL) {
      fprintf(fout,", ");
      nCellCenteredVariables=funcPrintVariableList(fout);
    }


    //print the definition of the zone
    fprintf(fout,"\nZone N=%i, E=%i\n",nRadialPoints*nPolarAnglePoints*nAzimuthAnglePoints,(nRadialPoints-1)*nPolarAnglePoints*(nAzimuthAnglePoints-1));

    if (nCellCenteredVariables==0) fprintf(fout,"ZONETYPE=FEBRICK DATAPACKING=BLOCK\n");
    else fprintf(fout,"ZONETYPE=FEBRICK DATAPACKING=BLOCK  VARLOCATION=([4-%i]=CELLCENTERED)\n",3+nCellCenteredVariables);

    //the numbering order of the mesh points: radial point, poinar angle point, zenith angle point
    double x[3];

    for (nRadial=0;nRadial<nRadialPoints;nRadial++) for (nPolar=0;nPolar<nPolarAnglePoints;nPolar++) {
      for (nAzimuth=0;nAzimuth<nAzimuthAnglePoints;nAzimuth++) {
        getCartesianNodeCoordinate(nRadial,nPolar,nAzimuth,x);
        fprintf(fout,"%e  ", x[0]);
      }

      fprintf(fout,"\n");
    }

    for (nRadial=0;nRadial<nRadialPoints;nRadial++) for (nPolar=0;nPolar<nPolarAnglePoints;nPolar++) {
      for (nAzimuth=0;nAzimuth<nAzimuthAnglePoints;nAzimuth++) {
        getCartesianNodeCoordinate(nRadial,nPolar,nAzimuth,x);
        fprintf(fout,"%e  ", x[1]);
      }

      fprintf(fout,"\n");
    }

    for (nRadial=0;nRadial<nRadialPoints;nRadial++) for (nPolar=0;nPolar<nPolarAnglePoints;nPolar++) {
      for (nAzimuth=0;nAzimuth<nAzimuthAnglePoints;nAzimuth++) {
        getCartesianNodeCoordinate(nRadial,nPolar,nAzimuth,x);
        fprintf(fout,"%e  ", x[2]);
      }

      fprintf(fout,"\n");
    }

    //print the field variables
    if ((nCellCenteredVariables!=0)&&(funcGetFieldVariable!=NULL)) {
      for (int nvar=0;nvar<nCellCenteredVariables;nvar++) {
        for (nRadial=0;nRadial<nRadialPoints-1;nRadial++) for (nPolar=0;nPolar<nPolarAnglePoints;nPolar++) {
          for (nAzimuth=0;nAzimuth<nAzimuthAnglePoints-1;nAzimuth++) {
            double x[3];

            getCartesianCellCoordinate(nRadial,nPolar,nAzimuth,x);
            fprintf(fout,"  %e", funcGetFieldVariable(nvar,x,GetSamplingElementNumber(nRadial,nPolar,nAzimuth)));
          }

          fprintf(fout,"\n");
        }
      }
    }

    //print the connectivity list
    //the numbering order of the mesh points: radial point, poinar angle point, zenith angle point
    for (nRadial=0;nRadial<nRadialPoints-1;nRadial++) for (nPolar=0;nPolar<nPolarAnglePoints;nPolar++) for (nAzimuth=0;nAzimuth<nAzimuthAnglePoints-1;nAzimuth++) {
      int n1,n2,n3,n4;
      int n5,n6,n7,n8;

      n1=1+nAzimuth+nPolar*nAzimuthAnglePoints+nRadial*(nPolarAnglePoints*nAzimuthAnglePoints);
      n2=1+nAzimuth+1+nPolar*nAzimuthAnglePoints+nRadial*(nPolarAnglePoints*nAzimuthAnglePoints);
      n3=1+nAzimuth+1+(nPolar==nPolarAnglePoints-1 ? 0 : nPolar+1)*nAzimuthAnglePoints+nRadial*(nPolarAnglePoints*nAzimuthAnglePoints);
      n4=1+nAzimuth+(nPolar==nPolarAnglePoints-1 ? 0 : nPolar+1)*nAzimuthAnglePoints+nRadial*(nPolarAnglePoints*nAzimuthAnglePoints);

      n5=1+nAzimuth+nPolar*nAzimuthAnglePoints+(nRadial+1)*(nPolarAnglePoints*nAzimuthAnglePoints);
      n6=1+nAzimuth+1+nPolar*nAzimuthAnglePoints+(nRadial+1)*(nPolarAnglePoints*nAzimuthAnglePoints);
      n7=1+nAzimuth+1+(nPolar==nPolarAnglePoints-1 ? 0 : nPolar+1)*nAzimuthAnglePoints+(nRadial+1)*(nPolarAnglePoints*nAzimuthAnglePoints);
      n8=1+nAzimuth+(nPolar==nPolarAnglePoints-1 ? 0 : nPolar+1)*nAzimuthAnglePoints+(nRadial+1)*(nPolarAnglePoints*nAzimuthAnglePoints);

      fprintf(fout,"%i %i %i %i   %i %i %i %i\n",n1,n2,n3,n4,n5,n6,n7,n8);
    }


    fclose(fout);
    FlushSamplingBuffer();
  }

};


#endif /* SPHERICALVOLUMEMESH_H_ */
