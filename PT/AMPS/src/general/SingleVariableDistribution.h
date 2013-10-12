//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * SingleVariableDistribution.h
 *
 *  Created on: Mar 24, 2012
 *      Author: vtenishe
 */

//$Id$
//the functions define a class that diftributes a variable according to a given distribution function

#ifndef SINGLEVARIABLEDISTRIBUTION_H_
#define SINGLEVARIABLEDISTRIBUTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <emmintrin.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>

#include "specfunc.h"


#define _SINGLE_VARIABLE_DISTRIBUTION_ON_    0
#define _SINGLE_VARIABLE_DISTRIBUTION_OFF_   1

#define _SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_  100000

#define _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_      0
#define _SIGNEL_VARIABLE_DISTRIBUTOIN__INTERPOLATION_MODE__EXPONENTIAL_ 1


#define _SINGLE_VARIABLE_DISTRIBUTION__MEMORY_PREFETCH_ _SINGLE_VARIABLE_DISTRIBUTION_OFF_



template <class T>
class cSingleVariableDistribution {
private:
  double xMin,xMax,DeltaX;
  int nCumulativeDistributionIntervals,ProbabilityDensityInterpolationMode;

  typedef double (*fProbabilityDistrvidution)(double,T*); //here T* is a pointer to constant parameter of the distribution function that is nessesaty for calcualtion (e.g the species number)
  fProbabilityDistrvidution ProbabilityDistributionFunction;

  struct cCumulativeDistributionTable {
    double x0,x1,p0,p1; //parameters of the probability density distribution
    double alpha,dF; //gradient of the probability density and increment of the cumiulative distribution
  } *CumulativeDistributionTable;

public:
  void Init(T* t) {
    double F,dF,dx,f1,f2;
    int i,imax,nInterval,nIntegrationStepsPerInterval;
    cCumulativeDistributionTable el;


    if (CumulativeDistributionTable!=NULL) exit(__LINE__,__FILE__,"Error: re-initialization of cSingleVariableDistribution object");

    DeltaX=xMax-xMin;
    CumulativeDistributionTable=new cCumulativeDistributionTable [nCumulativeDistributionIntervals];

    if (DeltaX==0.0) exit(__LINE__,__FILE__,"Error: parameters are not valid");

    //integrate the probability distribution function
    const int maxIntegrationStepMultiplier=500;
    int IntegrationStepMultiplier=50;

    begin:

    imax=IntegrationStepMultiplier*nCumulativeDistributionIntervals;
    dx=(xMax-xMin)/imax;

    f1=ProbabilityDistributionFunction(xMin,t);

    for (i=0,F=0.0;i<imax;i++) {
      f2=ProbabilityDistributionFunction(xMin+(i+1)*dx,t);

      if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
        if ((f1<0.0)||(f2<0.0)) exit(__LINE__,__FILE__,"Error: the probability density distribution must be positive!");

        F+=0.5*(f1+f2)*dx;
      }
      else if (ProbabilityDensityInterpolationMode==_SIGNEL_VARIABLE_DISTRIBUTOIN__INTERPOLATION_MODE__EXPONENTIAL_) {
        if ((f1<=0.0)||(f2<=0.0)) exit(__LINE__,__FILE__,"Error: the probability density distribution must be strictly positive!");

        exit(__LINE__,__FILE__,"Error: not implemented");
      }
      else exit(__LINE__,__FILE__,"Error: the option is not found");

      f1=f2;
    }

    //init the table of the cumilative distribution
    dF=F/nCumulativeDistributionIntervals;
    nIntegrationStepsPerInterval=0,nInterval=0;

    CumulativeDistributionTable[0].x0=0.0;
    CumulativeDistributionTable[0].p0=ProbabilityDistributionFunction(xMin,t);

    f1=ProbabilityDistributionFunction(xMin,t);

    for (i=0,F=0.0;i<imax;i++) {
      f2=ProbabilityDistributionFunction(xMin+(i+1)*dx,t);

      if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
        F+=0.5*(f1+f2)*dx;
      }
      else if (ProbabilityDensityInterpolationMode==_SIGNEL_VARIABLE_DISTRIBUTOIN__INTERPOLATION_MODE__EXPONENTIAL_) {
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
      else exit(__LINE__,__FILE__,"Error: the option is not found");

      nIntegrationStepsPerInterval++;
      f1=f2;

      if (F>(nInterval+1)*dF) {
        if (F>(nInterval+2)*dF) {
          if (IntegrationStepMultiplier<maxIntegrationStepMultiplier) {
            IntegrationStepMultiplier*=2;
            goto begin;
          }
          else exit(__LINE__,__FILE__,"Error: the probability distribution function is too steep");
        }

        if ((nIntegrationStepsPerInterval<10)&&(IntegrationStepMultiplier<maxIntegrationStepMultiplier)) {
          IntegrationStepMultiplier*=2;
          goto begin;
        }


        memcpy(&el,CumulativeDistributionTable+nInterval,sizeof(cCumulativeDistributionTable));

        el.x1=(i+1)*dx/DeltaX;
        el.p1=f2;

        if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
          el.alpha=(el.p1-el.p0)/(el.x1-el.x0)/DeltaX;
          el.dF=(el.p0+el.p1)*(el.x1-el.x0)*DeltaX/2.0;
        }
        else exit(__LINE__,__FILE__,"Error: not implemented");

        memcpy(CumulativeDistributionTable+nInterval,&el,sizeof(cCumulativeDistributionTable));

        nIntegrationStepsPerInterval=0;
        ++nInterval;

        if (nInterval<nCumulativeDistributionIntervals) {
          CumulativeDistributionTable[nInterval].x0=el.x1;
          CumulativeDistributionTable[nInterval].p0=el.p1;
        }

        if (nInterval>nCumulativeDistributionIntervals) exit(__LINE__,__FILE__,"Error: 'nInterval' index exeeds the number of elements in the table 'CumulativeDistributionTable'");
      }
    }

    if (nInterval<nCumulativeDistributionIntervals-1) exit(__LINE__,__FILE__,"Error: some of the elements in the table 'CumulativeDistributionTable' are not imitialied. Try to increase value of 'nCumulativeDistributionIntervals'");


    memcpy(&el,CumulativeDistributionTable+nCumulativeDistributionIntervals-1,sizeof(cCumulativeDistributionTable));

    el.x1=1.0;
    el.p1=ProbabilityDistributionFunction(xMax,t);

    if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
      el.alpha=(el.p1-el.p0)/(el.x1-el.x0)/DeltaX;
      el.dF=(el.p0+el.p1)*(el.x1-el.x0)*DeltaX/2.0;
    }
    else exit(__LINE__,__FILE__,"Error: not implemented");

    memcpy(CumulativeDistributionTable+nCumulativeDistributionIntervals-1,&el,sizeof(cCumulativeDistributionTable));

  }



  /*
  cSingleVariableDistribution(double xmin,double xmax,int nintervals,fProbabilityDistrvidution pfunc,int mode) {
    xMin=xmin,xMax=xmax,DeltaX=0.0,nCumulativeDistributionIntervals=nintervals;
    ProbabilityDistributionFunction=pfunc;
    ProbabilityDensityInterpolationMode=mode;
    CumulativeDistributionTable=NULL;
  }

  cSingleVariableDistribution(double xmin,double xmax,fProbabilityDistrvidution pfunc,int mode) {
    xMin=xmin,xMax=xmax,DeltaX=0.0,nCumulativeDistributionIntervals=_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_;
    ProbabilityDistributionFunction=pfunc;
    ProbabilityDensityInterpolationMode=mode;
    CumulativeDistributionTable=NULL;
  }
  */

  cSingleVariableDistribution() {
    xMin=0.0,xMax=0.0,DeltaX=0.0,nCumulativeDistributionIntervals=0;
    ProbabilityDistributionFunction=NULL;
    ProbabilityDensityInterpolationMode=0;
    CumulativeDistributionTable=NULL;
  }

  void Init(double xmin,double xmax,int nintervals,fProbabilityDistrvidution pfunc,int mode,T* t) {
    xMin=xmin,xMax=xmax,DeltaX=0.0,nCumulativeDistributionIntervals=nintervals;
    ProbabilityDistributionFunction=pfunc;
    ProbabilityDensityInterpolationMode=mode;
    CumulativeDistributionTable=NULL;

    Init(t);
  }

  void Init(double xmin,double xmax,fProbabilityDistrvidution pfunc,int mode,T* t) {
    xMin=xmin,xMax=xmax,DeltaX=0.0,nCumulativeDistributionIntervals=_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_;
    ProbabilityDistributionFunction=pfunc;
    ProbabilityDensityInterpolationMode=mode;
    CumulativeDistributionTable=NULL;

    Init(t);
  }

  inline double DistributeVariable() {
    int nInterval;
    double res=0.0;

#if _SINGLE_VARIABLE_DISTRIBUTION__MEMORY_PREFETCH_ == _SINGLE_VARIABLE_DISTRIBUTION_ON_
    _mm_prefetch(&xMin,_MM_HINT_T0);
#endif

    nInterval=(int)(rnd()*nCumulativeDistributionIntervals);

    if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
      cCumulativeDistributionTable t;

      memcpy(&t,CumulativeDistributionTable+nInterval,sizeof(cCumulativeDistributionTable));

      if (fabs(t.alpha)>0.0) {
        double dx,dp,A,b,c;

        dx=(t.x1-t.x0)*DeltaX;
        dp=t.p1-t.p0;
        b=dx*t.p0;

        A=2.0*dx*dp*rnd()*t.dF;
        c=A/(dx*dp)/(b*b+sqrt(b*b+A));
        res=xMin+DeltaX*t.x0+dx*c;
      }
      else res=xMin+DeltaX*t.x0+rnd()*(t.x1-t.x0)*DeltaX;
    }
    else exit(__LINE__,__FILE__,"Error: not implemented");

    return res;
  }

  void fPrintDistributionFunction(const char *fname,T* t) {
    FILE *fout;
    int i;
    double x,dx,norm=0.0,f1,f2;

    if (ProbabilityDistributionFunction==NULL) {
      std::cout << "$PREFIX:WARNING: ProbabilityDistributionFunction is NULL (file=" << __FILE__ << ", line=" << __LINE__ << ")" << std::endl;

      return;
    }

    fout=fopen(fname,"w");
    dx=(xMax-xMin)/nCumulativeDistributionIntervals;

    fprintf(fout,"TITLE=\"A single variable distribution function\"\n");
    fprintf(fout,"VARIABLES=\"x\", \"f(x)\"\n");

    //get the norm of the distribution function
    f1=ProbabilityDistributionFunction(xMin,t);

    for (i=0;i<nCumulativeDistributionIntervals;i++) {
      x=xMin+(i+1)*dx;
      f2=ProbabilityDistributionFunction(x,t);

      norm+=0.5*(f1+f2)*dx;
      f1=f2;
    }

    for (i=0;i<nCumulativeDistributionIntervals+1;i++) {
      x=xMin+i*dx;
      fprintf(fout,"%e  %e\n",x,ProbabilityDistributionFunction(x,t)/norm);
    }

    fclose(fout);
  }

  void fPrintCumulativeDistributionFunction(const char *fname) {
    FILE *fout;
    int i;
    double norm=0.0,summ=0.0;

    if (ProbabilityDistributionFunction==NULL) {
      std::cout << "$PREFIX:WARNING: ProbabilityDistributionFunction is NULL (file=" << __FILE__ << ", line=" << __LINE__ << ")" << std::endl;

      return;
    }

    fout=fopen(fname,"w");
    fprintf(fout,"TITLE=\"A single variable cumulative distribution function\"\n");
    fprintf(fout,"VARIABLES=\"x\", \"F(x)\"\n");

    //get the norm of the distribution function
    for (i=0;i<nCumulativeDistributionIntervals;i++) norm+=CumulativeDistributionTable[i].dF;

    for (i=0;i<nCumulativeDistributionIntervals;i++) {
      summ+=CumulativeDistributionTable[i].dF;
      fprintf(fout,"%e  %e\n",CumulativeDistributionTable[i].x0*DeltaX+xMin,summ/norm);
    }

    fprintf(fout,"%e  %e\n",xMax,summ/norm);
    fclose(fout);
  }


};



#endif /* SINGLEVARIABLEDISTRIBUTION_H_ */
