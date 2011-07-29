//====================================================
//$Id$
//====================================================
//the function for sampling of the particle's velocity/evergy distribution function

#include "pic.h"

const int PIC::DistributionFunctionSample::_LINEAR_SAMPLING_SCALE_=0,PIC::DistributionFunctionSample::_LOGARITHMIC_SAMPLING_SCALE_=1;
int PIC::DistributionFunctionSample::v2SamplingMode=_LINEAR_SAMPLING_SCALE_,PIC::DistributionFunctionSample::speedSamplingMode=_LINEAR_SAMPLING_SCALE_;
double PIC::DistributionFunctionSample::vMin=-1000.0,PIC::DistributionFunctionSample::vMax=1000.0;
long int PIC::DistributionFunctionSample::nSampledFunctionPoints=100;
double** PIC::DistributionFunctionSample::SamplingBuffer=NULL;
double** PIC::DistributionFunctionSample::SamplingLocations=NULL;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** PIC::DistributionFunctionSample::SampleNodes=NULL;
double PIC::DistributionFunctionSample::dV=0.0,PIC::DistributionFunctionSample::dV2=0.0,PIC::DistributionFunctionSample::dSpeed=0.0;
long int *PIC::DistributionFunctionSample::SampleLocalCellNumber=NULL;
int PIC::DistributionFunctionSample::nSamleLocations=0;
bool PIC::DistributionFunctionSample::SamplingInitializedFlag=false;

int PIC::DistributionFunctionSample::Sample_Velocity_Offset=0,PIC::DistributionFunctionSample::Sample_Speed_Offset=0;
int PIC::DistributionFunctionSample::Sample_V2_Offset=0,PIC::DistributionFunctionSample::SampleDataLength=0;

//====================================================
//init the sampling buffers
void PIC::DistributionFunctionSample::Init(double ProbeLocations[][DIM],int nProbeLocations) {
  int idim,nProbe,i,j,k;

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) printf("WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif


  nSamleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lenfths of the sampling intervals
  dV=1.05*(vMax-vMin)/(nSampledFunctionPoints-1);
  dSpeed=((speedSamplingMode==_LINEAR_SAMPLING_SCALE_) ? max(fabs(vMax),fabs(vMin))*sqrt(2.0) : log(max(fabs(vMax),fabs(vMin))*sqrt(2.0)))/(nSampledFunctionPoints-1);
  dV2=((v2SamplingMode==_LINEAR_SAMPLING_SCALE_) ? max(vMax*vMax,vMin*vMin)*sqrt(2.0) : log(max(vMax*vMax,vMin*vMin)*sqrt(2.0)))/(nSampledFunctionPoints-1);


  Sample_Velocity_Offset=0;
  SampleDataLength=DIM;

  Sample_Speed_Offset=SampleDataLength;
  SampleDataLength++;

  Sample_V2_Offset=SampleDataLength;
  SampleDataLength++;

  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nProbeLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nProbeLocations];
  SamplingLocations=new double* [nProbeLocations];
  SamplingLocations[0]=new double [DIM*nProbeLocations];

  SamplingBuffer=new double* [nProbeLocations];
  SamplingBuffer[0]=new double [nProbeLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1)];

  for (nProbe=1;nProbe<nProbeLocations;nProbe++) {
    SamplingLocations[nProbe]=SamplingLocations[nProbe-1]+DIM;
    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nProbeLocations;nProbe++) {
    for (idim=0;idim<DIM;idim++) SamplingLocations[nProbe][idim]=ProbeLocations[nProbe][idim];

    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void PIC::DistributionFunctionSample::flushSamplingBuffers() {
  long int i,TotalDataLength=nSamleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}
//====================================================
//return the offset where the sample data for the particular specie, sampling interval and the sampling point are located
long int PIC::DistributionFunctionSample::GetSampleDataOffset(int spec,int SampleVariableOffset) {
  long int offset;

  offset=spec*SampleDataLength*(nSampledFunctionPoints-1);
  offset+=SampleVariableOffset*(nSampledFunctionPoints-1);

  return offset;
}
//====================================================
//Sample the distribution function
void PIC::DistributionFunctionSample::SampleDistributionFnction(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* sampleNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,i,idim,offset;
  double LocalParticleWeight,v2;


    for (node=NULL,nProbe=0;nProbe<nSamleLocations;nProbe++) if (sampleNode==SampleNodes[nProbe]) {
      node=sampleNode;
      break;
    }

    if (node!=NULL) if (node->block!=NULL) {
      double *v;
      PIC::ParticleBuffer::byte *ParticleData;

      ptr=node->block->GetCenterNode(SampleLocalCellNumber[nProbe])->FirstCellParticle;

      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        v=PIC::ParticleBuffer::GetV(ParticleData);

        LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
        LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        for (v2=0.0,idim=0;idim<DIM;idim++) {
          i=(int)(v[idim]-vMin)/dV;
          v2+=pow(v[idim],2);
          offset=GetSampleDataOffset(spec,Sample_Velocity_Offset+idim);

          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight;
        }

        if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
          i=(int)(sqrt(v2)/dSpeed);
          offset=GetSampleDataOffset(spec,Sample_Speed_Offset);
          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight;

          i=(int)(v2/dV2);
          offset=GetSampleDataOffset(spec,Sample_V2_Offset);
          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight;
        }
        else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
          exit(__LINE__,__FILE__,"not implemented");
        }
        else exit(__LINE__,__FILE__,"Error: the option is unknown");

        ptr=PIC::ParticleBuffer::GetNext(ParticleData);
      }
    }

}


//====================================================
//print the distribution function into a file
void PIC::DistributionFunctionSample::printDistributionFunction(char *fname,int spec) {
  long int idim,nProbe,i,nVariable,thread,offset;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,dInterval=0.0;
  char str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);


  for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      sprintf(str,"%s.nSamplePoint=%ld.dat",fname,nProbe);
      fout=fopen(str,"w");

      printf("printing output file: %s.........         ",str);

      fprintf(fout,"\"TITLE=Distribution function at x=%e",SamplingLocations[nProbe][0]);
      for (idim=1;idim<DIM;idim++) fprintf(fout,", %e",SamplingLocations[nProbe][idim]);

      fprintf(fout,"\"\nVARIABLES=\"Sample Interval\",\"v\",\"f(vx)\"");
      if (DIM>1) fprintf(fout,",\"f(vy)\"");
      if (DIM>2) fprintf(fout,",\"f(vz)\"");

      fprintf(fout,",\"|v|\",\"f(|v|)\",\"v2\",\"f(v2)\"\n");

      //collect the sampled information from other processors
      for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]+=pipe.recv<double>(thread);
      }

      //normalize the distribution functions
      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        norm=0.0;
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) {
          if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
            if (nVariable<DIM) dInterval=dV;
            else if (nVariable==Sample_Speed_Offset) dInterval=dSpeed;
            else if (nVariable==Sample_V2_Offset) dInterval=dV2;
            else exit(__LINE__,__FILE__,"Error: unknown option");
          }
          else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
            exit(__LINE__,__FILE__,"not implemented");
          }
          else exit(__LINE__,__FILE__,"Error: the option is unknown");

          norm+=SamplingBuffer[nProbe][i+offset]*dInterval;
        }

        if (fabs(norm)>0.0) for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]/=norm;
      }

      //print the output file
      for (i=0;i<nSampledFunctionPoints-1;i++) {
        double v=0.0,v2=0.0,Speed=0.0;

        if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
          v=vMin+i*dV,v2=i*dV2,Speed=i*dSpeed;
        }
        else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
          exit(__LINE__,__FILE__,"not implemented");
        }
        else exit(__LINE__,__FILE__,"Error: the option is unknown");

        fprintf(fout,"%ld  %e ",i,v);

        for (idim=0;idim<DIM;idim++) {
          offset=GetSampleDataOffset(spec,idim);
          fprintf(fout,"  %e",SamplingBuffer[nProbe][i+offset]);
        }

        offset=GetSampleDataOffset(spec,Sample_Speed_Offset);
        fprintf(fout,"  %e  %e",Speed,SamplingBuffer[nProbe][i+offset]);

        offset=GetSampleDataOffset(spec,Sample_V2_Offset);
        fprintf(fout,"  %e  %e\n",v2,SamplingBuffer[nProbe][i+offset]);

      }

      //close the output file
      fclose(fout);
      printf("done.\n");
    }
    else {
      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        offset=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) {
          pipe.send(SamplingBuffer[nProbe][i+offset]);
          SamplingBuffer[nProbe][i+offset]=0.0; //this sampled information is stored by the root processor
        }
      }
    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_COMM_WORLD);
}

