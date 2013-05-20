//====================================================
//$Id$
//====================================================
//the function for sampling of the particle's velocity/evergy distribution function

#include "pic.h"

const int PIC::ParticleFluxDistributionSample::_LINEAR_SAMPLING_SCALE_=0,PIC::ParticleFluxDistributionSample::_LOGARITHMIC_SAMPLING_SCALE_=1;
int PIC::ParticleFluxDistributionSample::v2SamplingMode=_LINEAR_SAMPLING_SCALE_,PIC::ParticleFluxDistributionSample::speedSamplingMode=_LINEAR_SAMPLING_SCALE_;
double PIC::ParticleFluxDistributionSample::vMin=-1000.0,PIC::ParticleFluxDistributionSample::vMax=1000.0;
long int PIC::ParticleFluxDistributionSample::nSampledFunctionPoints=100;
double** PIC::ParticleFluxDistributionSample::SamplingBuffer=NULL;
double **PIC::ParticleFluxDistributionSample::SamplingLocations=NULL,**PIC::ParticleFluxDistributionSample::SamplingPointingDirections=NULL;
double** PIC::ParticleFluxDistributionSample::SamplingFlux=NULL;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** PIC::ParticleFluxDistributionSample::SampleNodes=NULL;
double PIC::ParticleFluxDistributionSample::dV2=0.0,PIC::ParticleFluxDistributionSample::dSpeed=0.0;
long int *PIC::ParticleFluxDistributionSample::SampleLocalCellNumber=NULL;
int PIC::ParticleFluxDistributionSample::nSamleLocations=0;
bool PIC::ParticleFluxDistributionSample::SamplingInitializedFlag=false;
double PIC::ParticleFluxDistributionSample::maxSamplingConeAngle=Pi,PIC::ParticleFluxDistributionSample::cosMaxSamplingConeAngle=-1.0;

int PIC::ParticleFluxDistributionSample::Sample_Speed_Offset=0;
int PIC::ParticleFluxDistributionSample::Sample_V2_Offset=0,PIC::ParticleFluxDistributionSample::SampleDataLength=0;

//====================================================
//init the sampling buffers
void PIC::ParticleFluxDistributionSample::Init(double ProbeLocations[][DIM],double ProbeDirections[][DIM],double maxConeAngles,int nProbeLocations) {
  int idim,nProbe,i,j,k;

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) printf("WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif

  maxSamplingConeAngle=maxConeAngles;
  cosMaxSamplingConeAngle=cos(maxSamplingConeAngle);

  nSamleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lenfths of the sampling intervals
  double t0,t1; //tempotary variables to satisfy intel c++ compiler

  dSpeed=((speedSamplingMode==_LINEAR_SAMPLING_SCALE_) ? max((t0=fabs(vMax)),(t1=fabs(vMin)))*sqrt(2.0) : log(max((t0=fabs(vMax)),(t1=fabs(vMin)))*sqrt(2.0)))/(nSampledFunctionPoints-1);
  dV2=((v2SamplingMode==_LINEAR_SAMPLING_SCALE_) ? max(t0=vMax*vMax,t1=vMin*vMin)*sqrt(2.0) : log(max(t0=vMax*vMax,t1=vMin*vMin)*sqrt(2.0)))/(nSampledFunctionPoints-1);

  Sample_Speed_Offset=0;
  SampleDataLength=1;

  Sample_V2_Offset=SampleDataLength;
  SampleDataLength++;

  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nProbeLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nProbeLocations];

  SamplingPointingDirections=new double* [nProbeLocations];
  SamplingPointingDirections[0]=new double [DIM*nProbeLocations];

  SamplingLocations=new double* [nProbeLocations];
  SamplingLocations[0]=new double [DIM*nProbeLocations];

  SamplingBuffer=new double* [nProbeLocations];
  SamplingBuffer[0]=new double [nProbeLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1)];

  SamplingFlux=new double *[nProbeLocations];
  SamplingFlux[0]=new double [nProbeLocations*PIC::nTotalSpecies];

  for (nProbe=1;nProbe<nProbeLocations;nProbe++) {
    SamplingLocations[nProbe]=SamplingLocations[nProbe-1]+DIM;
    SamplingPointingDirections[nProbe]=SamplingPointingDirections[nProbe-1]+DIM;

    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
    SamplingFlux[nProbe]=SamplingFlux[nProbe-1]+PIC::nTotalSpecies;
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nProbeLocations;nProbe++) {
    for (idim=0;idim<DIM;idim++) SamplingLocations[nProbe][idim]=ProbeLocations[nProbe][idim],SamplingPointingDirections[nProbe][idim]=ProbeDirections[nProbe][idim];
    for (int spec=0;spec<PIC::nTotalSpecies;spec++) SamplingFlux[nProbe][spec]=0.0;

    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void PIC::ParticleFluxDistributionSample::flushSamplingBuffers() {
  long int i,TotalDataLength=nSamleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
  for (i=0;i<nSamleLocations*PIC::nTotalSpecies;i++) SamplingFlux[0][i]=0.0;
}
//====================================================
//return the offset where the sample data for the particular specie, sampling interval and the sampling point are located
long int PIC::ParticleFluxDistributionSample::GetSampleDataOffset(int spec,int SampleVariableOffset) {
  long int offset;

  offset=spec*SampleDataLength*(nSampledFunctionPoints-1);
  offset+=SampleVariableOffset*(nSampledFunctionPoints-1);

  return offset;
}
//====================================================
//Sample the distribution function
void PIC::ParticleFluxDistributionSample::SampleDistributionFnction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,offset;

  for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) if (node->Thread==PIC::ThisThread) {
      double *v,speed,v2,cosPitchAngle,LocalParticleWeight;
      PIC::ParticleBuffer::byte *ParticleData;
      int i,j,k;
      double LocalTimeStep[PIC::nTotalSpecies],Measure;

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
      ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      if (ptr!=-1) {
        for (int s=0;s<PIC::nTotalSpecies;s++) LocalTimeStep[s]=node->block->GetLocalTimeStep(s);
        Measure=node->block->GetCenterNode(SampleLocalCellNumber[nProbe])->Measure;
      }

      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        v=PIC::ParticleBuffer::GetV(ParticleData);

        //get the angle between the particles's velocity and the direction of sampling
        double *l=SamplingPointingDirections[nProbe];

        v2=(v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]);
        speed=sqrt(v2);
        cosPitchAngle=-((v[0]*l[0])+(v[1]*l[1])+(v[2]*l[2]))/speed;

        if (cosPitchAngle>cosMaxSamplingConeAngle) {
          //the particle is within the cone of sample -> calculate the flux and sample the particle's properties

          LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
          LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

          SamplingFlux[nProbe][spec]+=speed*cosPitchAngle*LocalParticleWeight/(Measure*LocalTimeStep[spec]);

          if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
            i=(int)(speed/dSpeed);
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
        }

        ptr=PIC::ParticleBuffer::GetNext(ParticleData);
      }
    }

}

//====================================================
//print macroscopic parameters into a file
void PIC::ParticleFluxDistributionSample::printMacroscopicParameters(char *fname,int spec) {
  int nProbe,thread;
//  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  char str[_MAX_STRING_LENGTH_PIC_];


  //Print the distribution of fluxes
  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.openRecvAll();

    sprintf(str,"%s",fname);
    fout=fopen(str,"w");

    printf("printing output file: %s.........         ",str);
    fprintf(fout,"\"TITLE=Flux distribution");

    fprintf(fout,"\"\nVARIABLES=\"Sample Point\", \"Flux\" \n");
  }
  else pipe.openSend(0);

  double TotalFluxBuffer;


//  for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) {
  for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
    TotalFluxBuffer=SamplingFlux[nProbe][spec];

    if (PIC::ThisThread==0) {
      for (thread=1;thread<PIC::nTotalThreads;thread++) TotalFluxBuffer+=pipe.recv<double>(thread);

      fprintf(fout,"%i ",nProbe);
      fprintf(fout,"%e ",TotalFluxBuffer/PIC::LastSampleLength);

      fprintf(fout,"\n");
    }
    else {
      pipe.send(TotalFluxBuffer);
    }
  }


  //Print energy end velocity distribution of the sampled particles
  FILE *fout_v=NULL,*fout_v2=NULL;
  double vnorm[nSamleLocations],v2norm[nSamleLocations],fv,fv2;
  int voffset,v2offset,i;

  voffset=GetSampleDataOffset(spec,Sample_Speed_Offset);
  v2offset=GetSampleDataOffset(spec,Sample_V2_Offset);

  if (PIC::Mesh::mesh.ThisThread==0) {
    sprintf(str,"%s.f(v).dat",fname);
    fout_v=fopen(str,"w");
    fprintf(fout_v,"VARIABLES=\"v\"");
    for (nProbe=0;nProbe<nSamleLocations;nProbe++) fprintf(fout_v,", \"f(v,Sampling Set=%i)\"",nProbe);
    fprintf(fout_v,"\n");

    sprintf(str,"%s.f(v2).dat",fname);
    fout_v2=fopen(str,"w");
    fprintf(fout_v2,"VARIABLES=\"v^2\"");
    for (nProbe=0;nProbe<nSamleLocations;nProbe++) fprintf(fout_v2,", \"f(v^2,Sampling Set=%i)\"",nProbe);
    fprintf(fout_v2,"\n");



    //calculate the normal of distrinutions
//    for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) {
    for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
      vnorm[nProbe]=0.0,v2norm[nProbe]=0.0;

      for (i=0;i<nSampledFunctionPoints-1;i++) {
        vnorm[nProbe]=0.5*(SamplingBuffer[nProbe][voffset+i]+SamplingBuffer[nProbe][voffset+i+1])*dSpeed;
        v2norm[nProbe]=0.5*(SamplingBuffer[nProbe][v2offset+i]+SamplingBuffer[nProbe][v2offset+i+1])*dV2;
      }

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        vnorm[nProbe]+=pipe.recv<double>(thread);
        v2norm[nProbe]+=pipe.recv<double>(thread);
      }

      if (vnorm[nProbe]==0.0) vnorm[nProbe]=1.0;
      if (v2norm[nProbe]==0.0) v2norm[nProbe]=1.0;
    }

    //output the distribution function
    for (i=0;i<nSampledFunctionPoints;i++) {
      fprintf(fout_v,"%e ",i*dSpeed);
      fprintf(fout_v2,"%e ",i*dV2);

//      for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) {
      for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
        fv=SamplingBuffer[nProbe][voffset+i];
        fv2=SamplingBuffer[nProbe][v2offset+i];

        for (thread=1;thread<PIC::nTotalThreads;thread++) {
          fv+=pipe.recv<double>(thread);
          fv2+=pipe.recv<double>(thread);
        }

        fprintf(fout_v,"%e ",fv/vnorm[nProbe]);
        fprintf(fout_v2,"%e ",fv2/v2norm[nProbe]);
      }

      fprintf(fout_v,"\n");
      fprintf(fout_v2,"\n");
    }

    fclose(fout_v);
    fclose(fout_v2);
  }
  else {
    //calcualte the normal of the distribution function
//    for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) {
    for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
      vnorm[nProbe]=0.0,v2norm[nProbe]=0.0;

      for (i=0;i<nSampledFunctionPoints-1;i++) {
        vnorm[nProbe]=0.5*(SamplingBuffer[nProbe][voffset+i]+SamplingBuffer[nProbe][voffset+i+1])*dSpeed;
        v2norm[nProbe]=0.5*(SamplingBuffer[nProbe][v2offset+i]+SamplingBuffer[nProbe][v2offset+i+1])*dV2;
      }

      pipe.send(vnorm[nProbe]);
      pipe.send(v2norm[nProbe]);
    }

    //output the distribution function
    for (i=0;i<nSampledFunctionPoints;i++) {
//      for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) {
      for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
        fv=SamplingBuffer[nProbe][voffset+i];
        fv2=SamplingBuffer[nProbe][v2offset+i];

        pipe.send(fv);
        pipe.send(fv2);
      }
    }
  }


  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.closeRecvAll();
    fclose(fout);
    printf("done.\n");
  }
  else pipe.closeSend();

  MPI_Barrier(MPI_COMM_WORLD);
}

//====================================================
//print the distribution function into a file
void PIC::ParticleFluxDistributionSample::printDistributionFunction(char *fname,int spec) {

  /*

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
*/

  MPI_Barrier(MPI_COMM_WORLD);
}

