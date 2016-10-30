//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//sample energy distribution (relativistic)

#include "pic.h"

const int PIC::EnergyDistributionSampleRelativistic::_LINEAR_SAMPLING_SCALE_=0,PIC::EnergyDistributionSampleRelativistic::_LOGARITHMIC_SAMPLING_SCALE_=1;
double PIC::EnergyDistributionSampleRelativistic::eMin=0.0;
double PIC::EnergyDistributionSampleRelativistic::eMax=1000.0;
long int PIC::EnergyDistributionSampleRelativistic::nSampledFunctionPoints=100;
double** PIC::EnergyDistributionSampleRelativistic::SamplingBuffer=NULL;
double PIC::EnergyDistributionSampleRelativistic::SamplingLocations[][3]={{0.0,0.0,0.0}};
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** PIC::EnergyDistributionSampleRelativistic::SampleNodes=NULL;
double PIC::EnergyDistributionSampleRelativistic::dE=0.0;
long int *PIC::EnergyDistributionSampleRelativistic::SampleLocalCellNumber=NULL;
int PIC::EnergyDistributionSampleRelativistic::nSamleLocations=0;
bool PIC::EnergyDistributionSampleRelativistic::SamplingInitializedFlag=false;

int PIC::EnergyDistributionSampleRelativistic::EnergySamplingMode=PIC::EnergyDistributionSampleRelativistic::_LINEAR_SAMPLING_SCALE_;

int PIC::EnergyDistributionSampleRelativistic::Sample_Energy_Offset=0;
int PIC::EnergyDistributionSampleRelativistic::SampleDataLength=0;

vector<vector<int> > PIC::EnergyDistributionSampleRelativistic::CombinedSpeciesDistributionTable;

//====================================================
//init the sampling buffers
void PIC::EnergyDistributionSampleRelativistic::Init() {
  int idim,nProbe,i,j,k;

  if (SamplingInitializedFlag==true) exit(__LINE__,__FILE__,"Error: EnergyDistributionSampleRelativistic is already initialized");

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif


//  nSamleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lenfths of the sampling intervals
  double t0,t1; //tempotary variables to satisfy intel c++ compiler

  dE=1.05*(eMax-eMin)/(nSampledFunctionPoints-1);
  Sample_Energy_Offset=0;
  SampleDataLength=1;

  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nSamleLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nSamleLocations];

  SamplingBuffer=new double* [nSamleLocations];
  SamplingBuffer[0]=new double [nSamleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1)];

  for (nProbe=1;nProbe<nSamleLocations;nProbe++) {
    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void PIC::EnergyDistributionSampleRelativistic::flushSamplingBuffers() {
  long int i,TotalDataLength=nSamleLocations*PIC::nTotalSpecies*SampleDataLength*(nSampledFunctionPoints-1);
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}
//====================================================
//return the offset where the sample data for the particular specie, sampling interval and the sampling point are located
long int PIC::EnergyDistributionSampleRelativistic::GetSampleDataOffset(int spec,int SampleVariableOffset) {
  long int offset;

  offset=spec*SampleDataLength*(nSampledFunctionPoints-1);
  offset+=SampleVariableOffset*(nSampledFunctionPoints-1);

  return offset;
}
//====================================================
//Sample the distribution function
void PIC::EnergyDistributionSampleRelativistic::SampleDistributionFnction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,idim,offset;
  double LocalParticleWeight,e,mass,speed;

  for (node=SampleNodes[0],nProbe=0;nProbe<nSamleLocations;node=SampleNodes[++nProbe]) if (node->Thread==PIC::ThisThread) {
      double *v;
      PIC::ParticleBuffer::byte *ParticleData;
      int i,j,k;

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
      ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        v=PIC::ParticleBuffer::GetV(ParticleData);

        mass=PIC::MolecularData::GetMass(spec);

        LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
        LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        switch (EnergySamplingMode) {
        case _LINEAR_SAMPLING_SCALE_:
          speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
          e=Relativistic::Speed2E(speed,mass)*J2eV;

          i=(int)((e-eMin)/dE);
          offset=GetSampleDataOffset(spec,Sample_Energy_Offset);
          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight;

          break;
        case _LOGARITHMIC_SAMPLING_SCALE_:
          exit(__LINE__,__FILE__,"not implemented");

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is unknown");
        }

        ptr=PIC::ParticleBuffer::GetNext(ParticleData);
      }
    }

}


//====================================================
//print the distribution function into a file
void PIC::EnergyDistributionSampleRelativistic::printDistributionFunction(char *fname,int spec) {
  long int idim,nProbe,i,j,nVariable,thread,offset,s;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,e;
  char str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);

  //temporary sampling buffer
  double tempSamplingBuffer[SampleDataLength*nSampledFunctionPoints*PIC::nTotalSpecies];

  //find the list of species which energy distributions need to be paired
  vector<int> CombineDistributionSpecies;
  bool found=false;

  for (i=0;(i<CombinedSpeciesDistributionTable.size()) && (found==false);i++) {
    for (j=0;(j<CombinedSpeciesDistributionTable[i].size()) && (found==false);j++) if (CombinedSpeciesDistributionTable[i][j]==spec) {
      //a combination of species that need to be paired is found
      CombineDistributionSpecies=CombinedSpeciesDistributionTable[i];
      found=true;
    }
  }

  //if paired species are not found -> create a pair list that containes only one element
  if (found==false) CombineDistributionSpecies.push_back(spec);

  //output the distribution function
  for (nProbe=0;nProbe<nSamleLocations;nProbe++) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      sprintf(str,"%s.nSamplePoint=%ld.dat",fname,nProbe);
      fout=fopen(str,"w");

      fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);

      fprintf(fout,"\"TITLE=Distribution function at x=%e",SamplingLocations[nProbe][0]);
      for (idim=1;idim<DIM;idim++) fprintf(fout,", %e",SamplingLocations[nProbe][idim]);

      fprintf(fout,"\"\nVARIABLES=\"E [eV]\" \"f(E)\"\n");

      //init the tempSamplingBuffer
      for (nVariable=0;nVariable<SampleDataLength*nSampledFunctionPoints*PIC::nTotalSpecies;nVariable++) tempSamplingBuffer[nVariable]=0.0;

      //collect the sampled information from other processors
      int iCombinedSpecies,tempOffsetSpec;

      for (iCombinedSpecies=0;iCombinedSpecies<CombineDistributionSpecies.size();iCombinedSpecies++) {
        s=CombineDistributionSpecies[iCombinedSpecies];

        for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
          offset=GetSampleDataOffset(s,nVariable);
          tempOffsetSpec=GetSampleDataOffset(spec,nVariable);

          for (i=0;i<nSampledFunctionPoints-1;i++) if (thread==0) {
            tempSamplingBuffer[i+tempOffsetSpec]+=SamplingBuffer[nProbe][i+offset];
          }
          else {
            tempSamplingBuffer[i+tempOffsetSpec]+=pipe.recv<double>(thread);
          }
        }
      }

      //normalize the distribution functions
      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
        norm=0.0;
        tempOffsetSpec=GetSampleDataOffset(spec,nVariable);

        for (i=0;i<nSampledFunctionPoints-1;i++) {
          switch (EnergySamplingMode) {
          case _LINEAR_SAMPLING_SCALE_:
            norm+=tempSamplingBuffer[i+tempOffsetSpec]*dE;

            break;
          case _LOGARITHMIC_SAMPLING_SCALE_:
            exit(__LINE__,__FILE__,"not implemented");

            break;
          default:
            exit(__LINE__,__FILE__,"Error: the option is unknown");
          }
        }

        if (fabs(norm)>0.0) for (i=0;i<nSampledFunctionPoints-1;i++) tempSamplingBuffer[i+tempOffsetSpec]/=norm;
      }

      //print the output file
      for (i=0;i<nSampledFunctionPoints-1;i++) {
        double v=0.0,v2=0.0,Speed=0.0;

        switch (EnergySamplingMode) {
        case _LINEAR_SAMPLING_SCALE_:
          e=eMin+i*dE;

        break;
        case _LOGARITHMIC_SAMPLING_SCALE_:
          exit(__LINE__,__FILE__,"not implemented");

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is unknown");
         }

        tempOffsetSpec=GetSampleDataOffset(spec,Sample_Energy_Offset);
        fprintf(fout,"%e  %e\n",e,tempSamplingBuffer[i+tempOffsetSpec]);
      }

      //close the output file
      fclose(fout);
      fprintf(PIC::DiagnospticMessageStream,"done.\n");
    }
    else {
      //collect the sampled information from other processors
      int iCombinedSpecies;

      for (iCombinedSpecies=0;iCombinedSpecies<CombineDistributionSpecies.size();iCombinedSpecies++) {
        s=CombineDistributionSpecies[iCombinedSpecies];

        for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
          offset=GetSampleDataOffset(s,nVariable);

          for (i=0;i<nSampledFunctionPoints-1;i++) {
            pipe.send(SamplingBuffer[nProbe][i+offset]);
          }
        }
      }

    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}


void PIC::EnergyDistributionSampleRelativistic::AddCombinedCombinedParticleDistributionList(vector<int> CombinedSpecies) {
  CombinedSpeciesDistributionTable.push_back(CombinedSpecies);
}
