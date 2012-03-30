/*
 * SingleVariableDistribution.h
 *
 *  Created on: Mar 24, 2012
 *      Author: vtenishe
 */

//the functions define a class that diftributes a variable according to a given distribution function

#ifndef SINGLEVARIABLEDISTRIBUTION_H_
#define SINGLEVARIABLEDISTRIBUTION_H_

#define _SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_  100000

#define _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_      0
#define _SIGNEL_VARIABLE_DISTRIBUTOIN__INTERPOLATION_MODE__EXPONENTIAL_ 1

class cSingleVariableDistribution {
private:
  double xMin,xMax,DeltaX;
  int nCumulativeDistributionIntervals,ProbabilityDensityInterpolationMode;

  typedef double (*fProbabilityDistrvidution)(double);
  fProbabilityDistrvidution ProbabilityDistributionFunction;

  struct cCumulativeDistributionTable {
    float x0,x1,p0,p1; //parameters of the probability density distribution
    float alpha,dF; //gradient of the probability density and increment of the cumiulative distribution
  } *CumulativeDistributionTable;

  void Init() {
    double F,dF,dx,f1,f2;
    int i,imax,nInterval,nIntegrationStepsPerInterval;
    cCumulativeDistributionTable el;


    DeltaX=xMax-xMin;
    CumulativeDistributionTable=new cCumulativeDistributionTable [nCumulativeDistributionIntervals];

    //integrate the probability distribution function
    const int maxIntegrationStepMultiplier=500;
    int IntegrationStepMultiplier=50;

    begin:

    imax=IntegrationStepMultiplier*nCumulativeDistributionIntervals;
    dx=(xMax-xMin)/imax;

    f1=ProbabilityDistributionFunction(xMin);

    for (i=0,F=0.0;i<imax;i++) {
      f2=ProbabilityDistributionFunction(xMin+(i+1)*dx);

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
    CumulativeDistributionTable[0].p0=ProbabilityDistributionFunction(xMin);

    f1=ProbabilityDistributionFunction(xMin);

    for (i=0,F=0.0;i<imax;i++) {
      f2=ProbabilityDistributionFunction(xMin+(i+1)*dx);

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

          el.dF=(el.p0-el.alpha*(xMin+DeltaX*el.x0))*(el.x1-el.x0)*DeltaX+el.alpha/2.0*((el.x1-el.x0)*DeltaX)*((el.x1+el.x0)*DeltaX+2.0*xMin);
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
    el.p1=ProbabilityDistributionFunction(xMax);

    if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
      el.alpha=(el.p1-el.p0)/(el.x1-el.x0)/DeltaX;
      el.dF=(el.p0+el.p1)*dx/2.0;

      //el.dF=(el.p0-el.alpha*el.x0)*(el.x1-el.x0)*DeltaX+el.alpha/2.0*((el.x1-el.x0)*DeltaX)*((el.x1+el.x0)*DeltaX+2.0*xMin);
    }
    else exit(__LINE__,__FILE__,"Error: not implemented");

    memcpy(CumulativeDistributionTable+nCumulativeDistributionIntervals-1,&el,sizeof(cCumulativeDistributionTable));

  }

public:

  cSingleVariableDistribution(double xmin,double xmax,int nintervals,fProbabilityDistrvidution pfunc,int mode) {
    xMin=xmin,xMax=xmax,nCumulativeDistributionIntervals=nintervals;
    ProbabilityDistributionFunction=pfunc;
    ProbabilityDensityInterpolationMode=mode;
    CumulativeDistributionTable=NULL;

    Init();
  }

  cSingleVariableDistribution(double xmin,double xmax,fProbabilityDistrvidution pfunc,int mode) {
    cSingleVariableDistribution(xmin,xmax,_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_,pfunc,mode);
  }

  inline double DistributeVariable() {
    int nInterval;
    double res;

    nInterval=(int)(rnd()*nCumulativeDistributionIntervals);

    if (ProbabilityDensityInterpolationMode==_SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_) {
      cCumulativeDistributionTable t;

      memcpy(&t,CumulativeDistributionTable+nInterval,sizeof(cCumulativeDistributionTable));

      if (fabs(t.alpha)>0.0) res=t.x0+(-t.p0+sqrt(pow(t.p0,2)+2.0*t.alpha*t.dF*rnd()))/t.alpha;
      else res=t.x0+rnd()*(t.x1-t.x0);
    }
    else exit(__LINE__,__FILE__,"Error: not implemented");

    return res;
  }

  void fPrintDistributionFunction(const char *fname) {
    FILE *fout;
    int i;
    double x,dx,norm=0.0,f1,f2;

    fout=fopen(fname,"w");
    dx=(xMax-xMin)/nCumulativeDistributionIntervals;


    //get the norm of the distribution function
    f1=ProbabilityDistributionFunction(xMin);

    for (i=0;i<nCumulativeDistributionIntervals+1;i++) {
      x=xMin+i*dx;
      f2=ProbabilityDistributionFunction(x);

      norm+=0.5*(f1+f2)*dx;
      f1=f2;
    }

    for (i=0;i<nCumulativeDistributionIntervals+1;i++) {
      x=xMin+i*dx;
      fprintf(fout,"%e  %e\n",x,ProbabilityDistributionFunction(x)/norm);
    }

    fclose(fout);
  }

};



#endif /* SINGLEVARIABLEDISTRIBUTION_H_ */
