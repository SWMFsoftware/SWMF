//$Id$
//procedure for sampling of the size distribution of the dust grains injected into the domain


/*
 * SampleSizeDistribution.cpp
 *
 *  Created on: Dec 28, 2015
 *      Author: vtenishe
 */

#include "pic.h"
#include "Dust.h"
#include "Comet.h"

double Comet::Sampling::InjectedDustSizeDistribution::dLogRadius=0.0;
double* Comet::Sampling::InjectedDustSizeDistribution::SamplingBuffer=NULL;

//init the sampling procedure
void Comet::Sampling::InjectedDustSizeDistribution::Init() {
  int i;

  //determine the value of dLogRadius and allocate the sampling buffer
  dLogRadius=(log(ElectricallyChargedDust::maxDustRadius)-log(ElectricallyChargedDust::minDustRadius))/SamplingBufferLength;

  SamplingBuffer=new double [SamplingBufferLength];
  for (i=0;i<SamplingBufferLength;i++) SamplingBuffer[i]=0.0;

  //set up the output procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(EmptyParticleDataSamplingFunction,OutputSampledData);
}

void Comet::Sampling::InjectedDustSizeDistribution::AddParticleData(double Radius,double Weight) {
  int i;

  i=(int)(log(Radius/ElectricallyChargedDust::minDustRadius)/dLogRadius);
  if (i>=SamplingBufferLength) return;

  SamplingBuffer[i]+=Weight;
}

void Comet::Sampling::InjectedDustSizeDistribution::OutputSampledData(int DataOutputFileNumber) {
  FILE* fout;
  int ierr,i;

  //combine sampled data on the root processor
  double buffer[SamplingBufferLength];

  ierr=MPI_Reduce(SamplingBuffer,buffer,SamplingBufferLength,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  //output sampled data
  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/Comet.InjectedDustSizeDistribution.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Dust Grain Radius\", \"f(a)\", \"a^{-%e}\"\n",ElectricallyChargedDust::SizeDistribution::PowerIndex);

    //normalize the distribution function
    double a,da,a0,a1=ElectricallyChargedDust::minDustRadius,norm=0.0;
    for (i=0;i<SamplingBufferLength;i++) norm+=buffer[i]*dLogRadius;

    if (norm==0) norm=1.0;

    //output of the distribution function
    for (i=0;i<SamplingBufferLength;i++) {
      a0=a1;
      a1=ElectricallyChargedDust::minDustRadius*pow(exp(dLogRadius),i+1);

      a=0.5*(a0+a1);
      da=(a1-a0);

      fprintf(fout,"%e %e %e\n",a,buffer[i]*dLogRadius/(norm*da),pow(a,-ElectricallyChargedDust::SizeDistribution::PowerIndex));
    }


    //close the output file
    fclose(fout);
  }

  //reset the sampling buffer
  for (i=0;i<SamplingBufferLength;i++) SamplingBuffer[i]=0.0;
}


void Comet::Sampling::InjectedDustSizeDistribution::EmptyParticleDataSamplingFunction() {}
