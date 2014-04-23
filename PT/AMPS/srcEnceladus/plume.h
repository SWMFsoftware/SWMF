/*
 *  plume.h
 *  plume-fitter
 *
 *  Created by Valeriy Tenishev on 2/24/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <time.h>
#include <iostream>
#include <fstream>

#include "rnd.h"
#include "constants.h"

//the limiting values for plume model parameters

#define PlumeSourceRateMin 1.0E13
#define PlumeSourceRateMax 4.5E16 

#define PlumeBulkVelocityMin 400.0
#define PlumeBulkVelocityMax 700.0

#define PlumeTemperatureMin  70.0
#define PlumeTemperatureMax  220.0

#define PlumeTiltAngleMax              30.0/180.0*Pi  
#define PlumeAzimuthAngleDeflectionMax 20.0/180.0*Pi 


#define twoPi 2.0*Pi


#ifndef _PLUME_
#define _PLUME_



#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))




inline double qgaus(double func(const double), const double a, const double b) {
  static const double x[]={0.1488743389816312,0.4333953941292472,0.6794095682990244,0.8650633666889845,0.9739065285171717}; 
  static const double w[]={0.2955242247147529,0.2692667193099963,0.2190863625159821,0.1494513491505806,0.0666713443086881}; 
  int j;
  double xr,xm,dx,s;

  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0;
  for (j=0;j<5;j++) {
    dx=xr*x[j];
    s+=w[j]*(func(xm+dx)+func(xm-dx));
  }

  return s*=xr;
}


namespace DistributionFunction {
  extern double cosTheta,sinTheta,plumeBulkVelocity,plumeTempetarure;
  extern double m_over_k;

  double func(double v);

  /*{
    return pow(v,2)* exp(-m_over_k/plumeTempetarure*(pow(v*cosTheta-plumeBulkVelocity,2)+pow(v*sinTheta,2)));
  }*/
}

class Cplume;

namespace ColumnDensityIntegtant {
  extern Cplume *plume;
  extern double func(double t);
  extern double *Ray,*rInit;
}


class Cplume {
public:
  double SourceRate,Temperature,BulkVelocity;
  double plumeDirection[3],plumePosition[3];
  double refTiltAngle,refAzimuthAngle;
  double TiltAngle,AzimuthAngle;
  

  #define PLUME_BUFFER_LENGTH 8
  double SourceRate_bestConfiguration,Temperature_bestConfiguration,BulkVelocity_bestConfiguration;
  double plumeDirection_bestConfiguration[3],TiltAngle_bestConfiguration,AzimuthAngle_bestConfiguration;

  void RefineModelParameters(double L) {
    double dVal,maxVal,minVal;
    int idim;

    //source rate
    double dLogVal,maxLogVal,minLogVal,logVal;
 
    maxLogVal=log(PlumeSourceRateMax);
    minLogVal=log(PlumeSourceRateMin);
    dLogVal=(maxLogVal-minLogVal)/pow(2.0,L);

    logVal=log(SourceRate_bestConfiguration);  

//    maxLogVal=min(maxLogVal,logVal+dLogVal);
//    minLogVal=max(minLogVal,logVal-dLogVal);


    maxLogVal=(maxLogVal<logVal+dLogVal) ? maxLogVal : logVal+dLogVal;
    minLogVal=(minLogVal>logVal-dLogVal) ? maxLogVal : logVal+dLogVal;


    logVal=minLogVal+rnd()*(maxLogVal-minLogVal);
    SourceRate=exp(logVal);

    //temperature
    dVal=(PlumeTemperatureMax-PlumeTemperatureMin)/pow(2.0,L);
    maxVal=min(PlumeTemperatureMax,Temperature_bestConfiguration+dVal);
    minVal=max(PlumeTemperatureMin,Temperature_bestConfiguration-dVal);
    Temperature=minVal+rnd()*(maxVal-minVal);

    //bulk velocity
    dVal=(PlumeBulkVelocityMax-PlumeBulkVelocityMin)/pow(2.0,L);
    maxVal=min(PlumeBulkVelocityMax,BulkVelocity_bestConfiguration+dVal);
    minVal=max(PlumeBulkVelocityMin,BulkVelocity_bestConfiguration-dVal);
    BulkVelocity=minVal+rnd()*(maxVal-minVal);

    //tilt angle: the cosine of the tilt angle is distributed uniformly  
    double cosTiltAngle,cosTiltAngleMax=max(0.0,cos(1.3*refTiltAngle));
  
    cosTiltAngle=cos(TiltAngle_bestConfiguration); 
    dVal=(1.0-cosTiltAngleMax)/pow(2.0,L);  

    maxVal=min(1.0,cosTiltAngle+dVal);
    minVal=max(cosTiltAngleMax,cosTiltAngle-dVal);
    cosTiltAngle=minVal+rnd()*(maxVal-minVal);
    TiltAngle=acos(cosTiltAngle);   

    //azimuth angle 
    double dAzimuthAngle=Pi/pow(2.0,L);
    AzimuthAngle=AzimuthAngle_bestConfiguration+(-1.0+2.0*rnd())*dAzimuthAngle;  

    if (AzimuthAngle<0.0) AzimuthAngle+=2.0*Pi;
    if (AzimuthAngle>2.0*Pi) AzimuthAngle-=2.0*Pi;

    //get the plume direction vector
    //construc the coordinate system for the plume

    //radial direction -> z-axys
    //direction to the local north -> x-axys
    //y-axys -> the "right" coordinate system

    double l,l1,coord[3][3];
    
    //get the direction to the local north
    const double north[]={0.0,0.0,1.0};
  
    for (idim=0,l=0.0;idim<3;idim++) l+=pow(plumePosition[idim],2);
    for (idim=0,l=sqrt(l);idim<3;idim++) coord[2][idim]=plumePosition[idim]/l;

    for (idim=0,l=0.0;idim<3;idim++) l+=north[idim]*coord[2][idim];

    for (idim=0,l1=0.0;idim<3;idim++) {
      coord[0][idim]=north[idim]-l*coord[2][idim];
      l1+=pow(coord[0][idim],2);
    }

    for (idim=0,l1=sqrt(l1);idim<3;idim++) coord[0][idim]/=l1;

    coord[1][0]=-(coord[2][1]*coord[0][2]-coord[0][1]*coord[2][2]);
    coord[1][1]=+(coord[2][0]*coord[0][2]-coord[0][0]*coord[2][2]);
    coord[1][2]=-(coord[2][0]*coord[0][1]-coord[0][0]*coord[2][1]);

    for (idim=0;idim<3;idim++) plumeDirection[idim]=cos(TiltAngle)*coord[2][idim];   //vertical direction
    for (idim=0;idim<3;idim++) plumeDirection[idim]+=sin(TiltAngle)*sin(AzimuthAngle)*coord[1][idim]; //direction to the east
    for (idim=0;idim<3;idim++) plumeDirection[idim]+=sin(TiltAngle)*cos(AzimuthAngle)*coord[0][idim]; //direction to the north
  }


  void InitModelParameters(double R,double Latitude,double WLongitude,double RefTiltAngle,double RefAzimuthAngle) { 
    SourceRate_bestConfiguration=(SourceRate=exp(0.5*(log(PlumeSourceRateMin)+log(PlumeSourceRateMax))));
    Temperature_bestConfiguration=(Temperature=0.5*(PlumeTemperatureMin+PlumeTemperatureMax));
    BulkVelocity_bestConfiguration=(BulkVelocity=0.5*(PlumeBulkVelocityMin+PlumeBulkVelocityMax));

    TiltAngle_bestConfiguration=(TiltAngle=RefTiltAngle),AzimuthAngle_bestConfiguration=(AzimuthAngle=RefAzimuthAngle);
    refTiltAngle=RefTiltAngle,refAzimuthAngle=RefAzimuthAngle;

    plumePosition[2]=R*sin(Latitude);
    plumePosition[1]=-R*cos(Latitude)*sin(WLongitude);
    plumePosition[0]=R*cos(Latitude)*cos(WLongitude);

    //get the direction of the plume 
    int idim;
    double l,l1,coord[3][3];
   
    //get the direction to the local north
    const double north[]={0.0,0.0,1.0};
 
    for (idim=0,l=0.0;idim<3;idim++) l+=pow(plumePosition[idim],2);
    for (idim=0,l=sqrt(l);idim<3;idim++) coord[2][idim]=plumePosition[idim]/l;

    for (idim=0,l=0.0;idim<3;idim++) l+=north[idim]*coord[2][idim];

    for (idim=0,l1=0.0;idim<3;idim++) {
      coord[0][idim]=north[idim]-l*coord[2][idim];
      l1+=pow(coord[0][idim],2);
    }

    for (idim=0,l1=sqrt(l1);idim<3;idim++) coord[0][idim]/=l1;

    coord[1][0]=-(coord[2][1]*coord[0][2]-coord[0][1]*coord[2][2]);
    coord[1][1]=+(coord[2][0]*coord[0][2]-coord[0][0]*coord[2][2]);
    coord[1][2]=-(coord[2][0]*coord[0][1]-coord[0][0]*coord[2][1]);

    for (idim=0;idim<3;idim++) plumeDirection[idim]=cos(TiltAngle)*coord[2][idim];   //vertical direction
    for (idim=0;idim<3;idim++) plumeDirection[idim]+=sin(TiltAngle)*sin(AzimuthAngle)*coord[1][idim]; //direction to the east
    for (idim=0;idim<3;idim++) plumeDirection[idim]+=sin(TiltAngle)*cos(AzimuthAngle)*coord[0][idim]; //direction to the north
  }
  
  void SaveBestConfiguration() {
    int i;
	
    SourceRate_bestConfiguration=SourceRate;
    Temperature_bestConfiguration=Temperature;
    BulkVelocity_bestConfiguration=BulkVelocity; 
    for (i=0;i<3;i++) plumeDirection_bestConfiguration[i]=plumeDirection[i];
    TiltAngle_bestConfiguration=TiltAngle;
    AzimuthAngle_bestConfiguration=AzimuthAngle;
  }
  
  void SetBestConfiguration(double *buffer) {
    int cnt=0,i;
	
    SourceRate_bestConfiguration=(SourceRate=buffer[cnt++]);
    Temperature_bestConfiguration=(Temperature=buffer[cnt++]);
    BulkVelocity_bestConfiguration=(BulkVelocity=buffer[cnt++]); 
    for (i=0;i<3;i++) plumeDirection_bestConfiguration[i]=(plumeDirection[i]=buffer[cnt++]);
    TiltAngle_bestConfiguration=(TiltAngle=buffer[cnt++]);
    AzimuthAngle_bestConfiguration=(AzimuthAngle=buffer[cnt++]);

    if (cnt!=PLUME_BUFFER_LENGTH) exit(0);
  } 
  
  
  int GetBestConfiguration(double *buffer=NULL) {
    int cnt=0,i;

    if (buffer!=0) {
      buffer[cnt++]=SourceRate_bestConfiguration;
      buffer[cnt++]=Temperature_bestConfiguration;
      buffer[cnt++]=BulkVelocity_bestConfiguration;
      for (i=0;i<3;i++) buffer[cnt++]=plumeDirection_bestConfiguration[i];
      buffer[cnt++]=TiltAngle_bestConfiguration;
      buffer[cnt++]=AzimuthAngle_bestConfiguration;

      if (cnt!=PLUME_BUFFER_LENGTH) exit(0);
    }

    return PLUME_BUFFER_LENGTH;
  }
  

  Cplume() {
    SourceRate=0.0,Temperature=0.0,BulkVelocity=0.0;
    refTiltAngle=-1.0,refAzimuthAngle=-1.0;
    TiltAngle=-1.0,AzimuthAngle=-1.0;
  }
  
  double GetDensity(double *x) {
    double rVectTimesPlumePosition,cosTheta,sinTheta,l,r,rVect[3];
    int idim;

    for (l=0.0,r=0.0,rVectTimesPlumePosition=0.0,idim=0;idim<3;idim++) {
      rVect[idim]=x[idim]-plumePosition[idim];
      l+=rVect[idim]*plumeDirection[idim],r+=pow(rVect[idim],2);
      rVectTimesPlumePosition+=rVect[idim]*plumePosition[idim];
    }
	
    r=sqrt(r);
    cosTheta=l/r;
    sinTheta=sqrt(1-pow(cosTheta,2));
	
    DistributionFunction::cosTheta=cosTheta;
    DistributionFunction::sinTheta=sinTheta;
    DistributionFunction::plumeBulkVelocity=BulkVelocity;  
    DistributionFunction::plumeTempetarure=Temperature; 

//    if (cosTheta<0.0) return 0.0;
    if (rVectTimesPlumePosition<0.0) return 0.0; 


	  
    const int nTotSteps=10;
    const double vmax=1.0E4;
	
    double res=0.0,dV=vmax/nTotSteps;
    int nstep;
	
    for (nstep=0;nstep<nTotSteps;nstep++) res+=qgaus(DistributionFunction::func,nstep*dV,(nstep+1)*dV);
	
    return SourceRate*res/pow(r,2); 
  }
  
  double GetColumnDesnity(double *rInit,double *Ray) {
    int nstep;
    double res=0.0;
	
    const int nTotalSteps=10;
    const double lengthMax=5.0E4*1.0E3;
    const double dLogLength=log(lengthMax)/nTotalSteps;
	
    ColumnDensityIntegtant::Ray=Ray;
    ColumnDensityIntegtant::rInit=rInit;
    ColumnDensityIntegtant::plume=this;
	
    for (nstep=0;nstep<nTotalSteps;nstep++) res+=qgaus(ColumnDensityIntegtant::func,exp(dLogLength*nstep)-1.0,exp(dLogLength*(nstep+1))-1.0);
  
    return res;
  }  

  double CalculateTotalProductionRate() {
    double theta,v,cosTheta,sinTheta,beta,res;
    long int ntest;

    const long int nTotTests=1000000;
    const double vmax=2.0E3;

    res=0.0;
    beta=DistributionFunction::m_over_k/Temperature_bestConfiguration; 

    for (ntest=0;ntest<nTotTests;ntest++) {
      v=rnd()*vmax;
      theta=rnd()*Pi/2.0;

      cosTheta=cos(theta);
      sinTheta=sin(theta);

      res+=pow(v,3)*exp(-beta*(pow(v*cosTheta-BulkVelocity_bestConfiguration,2)+pow(v*sinTheta,2)))*sinTheta;    
    } 

    return res/ntest*Pi/2.0*vmax*2.0*Pi*SourceRate_bestConfiguration;    
 }


  void StampPlumParameters(int nPlume,FILE* fout) {

    fprintf(fout,"\nPlume[%i].Temperature=%e;\n",nPlume,Temperature);
    fprintf(fout,"Plume[%i].BulkVelocity=%e;\n",nPlume,BulkVelocity);
    fprintf(fout,"Plume[%i].SourceRate=%e;\n",nPlume,CalculateTotalProductionRate());

    fprintf(fout,"Plume[%i].plumeDirection[0]=%e;\n",nPlume,plumeDirection[0]);
    fprintf(fout,"Plume[%i].plumeDirection[1]=%e;\n",nPlume,plumeDirection[1]);
    fprintf(fout,"Plume[%i].plumeDirection[2]=%e;\n",nPlume,plumeDirection[2]);

    fprintf(fout,"Plume[%i].plumePosition[0]=%e;\n",nPlume,plumePosition[0]);
    fprintf(fout,"Plume[%i].plumePosition[1]=%e;\n",nPlume,plumePosition[1]);
    fprintf(fout,"Plume[%i].plumePosition[2]=%e;\n\n",nPlume,plumePosition[2]);
  }
  
};


//===================================================================================================
/*
double ColumnDensityIntegtant::func(double t) {
  double x[3];
  int idim;
  
  for (idim=0;idim<3;idim++) x[idim]=rInit[idim]+t*Ray[idim];
  return plume->GetDensity(x);
}
*/
//double ColumnDensityIntegtant::func(double t);


#endif
