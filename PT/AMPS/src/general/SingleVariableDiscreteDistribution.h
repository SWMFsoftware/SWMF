//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

//$Id$
//the functions define a class that diftributes a variable according to a given distribution function (DISCRETE DISTRIBUTION)!!!!

#ifndef _SINGLEVARIABLE_DISCRETE_DISTRIBUTION_H_
#define _SINGLEVARIABLE_DISCRETE_DISTRIBUTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <emmintrin.h>
#include <errno.h>
#include <unistd.h>
#include <sys/time.h>





#define _SINGLE_VARIABLE_DISCRETE_DISTRIBUTION_ON_    0
#define _SINGLE_VARIABLE_DISCRETE_DISTRIBUTION_OFF_   1

#define _SINGLE_VARIABLE_DISCRETE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_  100000





template <class T>
class cSingleVariableDiscreteDistribution {
public:
  int xMin,xMax;
  int nCumulativeDistributionIntervals;

  typedef double (*fProbabilityDistrvidution)(int,T*);
  fProbabilityDistrvidution ProbabilityDistributionFunction;

  int *CumulativeDistributionTable;

public:
  void Init(T* t) {
    double F;
    int i,j,jStart=0,jFinish=0,iDeltaJmax=-1,DeltaJmax=-1,DeltaJ;

    if (CumulativeDistributionTable==NULL) exit(__LINE__,__FILE__,"Error: the cumulative distribution buffer is not allocated");


    //normalize the distribution function
    double norm=0.0,summ=0.0;
    for (i=xMin;i<=xMax;i++) norm+=ProbabilityDistributionFunction(i,t);

    if (norm<=0.0) {
      exit(__LINE__,__FILE__,"Error: the distribution norm must be positive");
    }

    //initialize the cumulative distribution function
    for (F=0.0,i=xMin;i<=xMax;i++) {
      summ+=ProbabilityDistributionFunction(i,t);
      F=summ/norm; //to minimize the effect of the round error, F is calculated the same way as 'norm' to garantee that the final value of F is 1
      jFinish=(int)(F*nCumulativeDistributionIntervals);

      for (DeltaJ=0;jStart<=std::min(jFinish,nCumulativeDistributionIntervals-1);jStart++) {
        CumulativeDistributionTable[jStart]=i;
        DeltaJ++;
      }

      if (DeltaJmax<DeltaJ) DeltaJmax=DeltaJ,iDeltaJmax=i;
      if (jStart>=nCumulativeDistributionIntervals) break;
    }

    //finish the rest of the cumulative distribution table
    if (jStart<nCumulativeDistributionIntervals) {
      if (jStart==0) {
        exit(__LINE__,__FILE__,"Error: the table is empty");
      }
      else {
        cout << "AMPS:: there is an empty element at the end of the cumulative distribution function table (file=" << __FILE__ << ", line=" << __LINE__ << ")" << endl;
      }

      for (;jStart<nCumulativeDistributionIntervals;jStart++) CumulativeDistributionTable[jStart]=iDeltaJmax;
    }

  }

  void InitArray(double *ProbabilityArray,int nProbabilityArray,int nintervals) {
    double F;
    int i,j,jStart=0,jFinish=0,iDeltaJmax=-1,DeltaJmax=-1,DeltaJ;

    //allocate the cumulative distribution table
    xMin=0;
    xMax=nProbabilityArray-1;
    nCumulativeDistributionIntervals=nintervals;

    if (CumulativeDistributionTable!=NULL) delete [] CumulativeDistributionTable;
    CumulativeDistributionTable=new int [nCumulativeDistributionIntervals];

    //normalize the distribution function
    double norm=0.0,summ=0.0;
    for (i=xMin;i<=xMax;i++) norm+=ProbabilityArray[i];

    if (norm<=0.0) {
      exit(__LINE__,__FILE__,"Error: the distribution norm must be positive");
    }

    //initialize the cumulative distribution function
    for (F=0.0,i=xMin;i<=xMax;i++) {
      summ+=ProbabilityArray[i];
      F=summ/norm; //to minimize the effect of the round error, F is calculated the same way as 'norm' to garantee that the final value of F is 1
      jFinish=(int)(F*nCumulativeDistributionIntervals);

      for (DeltaJ=0;jStart<=std::min(jFinish,nCumulativeDistributionIntervals-1);jStart++) {
        CumulativeDistributionTable[jStart]=i;
        DeltaJ++;
      }

      if (DeltaJmax<DeltaJ) DeltaJmax=DeltaJ,iDeltaJmax=i;
      if (jStart>=nCumulativeDistributionIntervals) break;
    }

    //finish the rest of the cumulative distribution table
    if (jStart<nCumulativeDistributionIntervals) {
      if (jStart==0) {
        exit(__LINE__,__FILE__,"Error: the table is empty");
      }
      else {
        cout << "AMPS:: there is an empty element at the end of the cumulative distribution function table (file=" << __FILE__ << ", line=" << __LINE__ << ")" << endl;
      }

      for (;jStart<nCumulativeDistributionIntervals;jStart++) CumulativeDistributionTable[jStart]=iDeltaJmax;
    }
  }


  cSingleVariableDiscreteDistribution(int xmin,int xmax,int nintervals,fProbabilityDistrvidution pfunc) {
    xMin=xmin,xMax=xmax,nCumulativeDistributionIntervals=nintervals;
    ProbabilityDistributionFunction=pfunc;
    CumulativeDistributionTable=new int [nCumulativeDistributionIntervals];
  }

  cSingleVariableDiscreteDistribution(int xmin,int xmax,fProbabilityDistrvidution pfunc) {
    xMin=xmin,xMax=xmax,nCumulativeDistributionIntervals=_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_;
    ProbabilityDistributionFunction=pfunc;
    CumulativeDistributionTable=new int [nCumulativeDistributionIntervals];
  }

  cSingleVariableDiscreteDistribution() {
    xMin=0,xMax=0,nCumulativeDistributionIntervals=_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_;
    ProbabilityDistributionFunction=NULL;
    CumulativeDistributionTable=NULL;
  }

  void SetLimits(int xmin,int xmax,fProbabilityDistrvidution pfunc) {
    xMin=xmin,xMax=xmax,nCumulativeDistributionIntervals=_SINGLE_VARIABLE_DISTRIBUTION__DEFAULT__CUMULATIVE_DISTRIBUTION_INTERVAL_;
    ProbabilityDistributionFunction=pfunc;
    CumulativeDistributionTable=new int [nCumulativeDistributionIntervals];
  }

  inline int DistributeVariable() {
    return CumulativeDistributionTable[(int)(rnd()*nCumulativeDistributionIntervals)];
  }

  void fPrintDistributionFunction(const char *fname, T* t) {
    FILE *fout;
    int i;
    double norm=0.0;

    fout=fopen(fname,"w");

    fprintf(fout,"TITLE=\"A single variable distribution function\"\n");
    fprintf(fout,"VARIABLES=\"x\", \"f(x)\"\n");

    //get the norm of the distribution function
    for (i=xMin;i<=xMax;i++) norm+=ProbabilityDistributionFunction(i,t);

    //output the normalized distribution function
    for (i=xMin;i<=xMax;i++) {
      fprintf(fout,"%i  %e\n",i,ProbabilityDistributionFunction(i,t)/norm);
    }

    fclose(fout);
  }

  void fPrintCumulativeDistributionFunction(const char *fname,T* t) {
    FILE *fout;
    int i;
    double norm=0.0,summ=0.0;

    fout=fopen(fname,"w");
    fprintf(fout,"TITLE=\"A single variable cumulative distribution function\"\n");
    fprintf(fout,"VARIABLES=\"x\", \"F(x)\"\n");

    //get the norm of the distribution function
    for (i=xMin;i<=xMax;i++) norm+=ProbabilityDistributionFunction(i,t);

    for (i=xMin;i<=xMax;i++) {
      summ+=ProbabilityDistributionFunction(i,t);
      fprintf(fout,"%i  %e\n",i,summ/norm);
    }

    fclose(fout);
  }


};



#endif
