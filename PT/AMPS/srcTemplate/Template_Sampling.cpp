//$Id$

#include "pic.h"
#include "Template.h"

/*****************************************************************************
 * This file contains implementations for methods for sampling & printing
 * physical parameters to a separate file
 *****************************************************************************/

//============================================================================
// sampling of a generic parameter Gamma;
// if needed - use as a template, otherwise - remove / comment out

//----------------------------------------------------------------------------
// container for data
double Template::Sampling::Gamma::SampleBuffer[Template::Sampling::Gamma::nSampleInterval];

//----------------------------------------------------------------------------
// initilization function for sampling procedure
void   Template::Sampling::Gamma::Init() {
  // reset output files' counter
  nOutputFile = 0;
  // reset data container
  for(int iSampleInterval=0;iSampleInterval<nSampleInterval;iSampleInterval++){
    SampleBuffer[iSampleInterval] = 0.0;
  }
}

//----------------------------------------------------------------------------
// initilization function for sampling procedure
void   Template::Sampling::Gamma::SampleData() {
}

//----------------------------------------------------------------------------
// create output file
void   Template::Sampling::Gamma::PrintSampleData(int DataOutputFileNumber){
  //collect the distribution function from all processors on root
  double recvBuffer[nSampleInterval];
  MPI_Reduce(SampleBuffer,recvBuffer,nSampleInterval,
	     MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  // flush the data container once collected
  for(int iInterval=0;iInterval<nSampleInterval;iInterval++)
    SampleBuffer[iInterval] = 0.0;

  //..........................................................................
  // root processor prints data to file
  if(PIC::ThisThread==0) {
    // normalize the distribution function
    double norm = 0.0;
    for(int iInterval=0;iInterval<nSampleInterval;iInterval++)
      norm+=recvBuffer[iInterval];
    if(norm > 0.0) for(int iInterval=0;iInterval<nSampleInterval;iInterval++)
		     recvBuffer[iInterval]/=norm;

    //........................................................................
    char fname[_MAX_STRING_LENGTH_PIC_];
    FILE *fout;
    // open file
    sprintf(fname,"%s/amps.GammaDistribution.out=%i.dat",
	    PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fout=fopen(fname,"w");
    // create file's header in Teplot format
    fprintf(fout,"VARIABLES=\"Gamma\",\"f(Gamma)\"\n");
    // print data to file
    for(int iInterval=0;iInterval<nSampleInterval;iInterval++)
      fprintf(fout,"%e %e\n",iInterval,recvBuffer[iInterval]);
    //close file
    fclose(fout);
  }
}
