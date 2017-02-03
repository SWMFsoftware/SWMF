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
double PIC::EnergyDistributionSampleRelativistic::log10eMin=0.0;
double PIC::EnergyDistributionSampleRelativistic::log10eMax=0.0;
long int PIC::EnergyDistributionSampleRelativistic::nSampledFunctionPoints=100;
double*** PIC::EnergyDistributionSampleRelativistic::SamplingBuffer=NULL;
double PIC::EnergyDistributionSampleRelativistic::SamplingLocations[][3]={{0.0,0.0,0.0}};
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** PIC::EnergyDistributionSampleRelativistic::SampleNodes=NULL;
double PIC::EnergyDistributionSampleRelativistic::dE=0.0;
double PIC::EnergyDistributionSampleRelativistic::log10dE=0.0;
long int *PIC::EnergyDistributionSampleRelativistic::SampleLocalCellNumber=NULL;
int PIC::EnergyDistributionSampleRelativistic::nSamleLocations=0;
bool PIC::EnergyDistributionSampleRelativistic::SamplingInitializedFlag=false;

int PIC::EnergyDistributionSampleRelativistic::EnergySamplingMode=PIC::EnergyDistributionSampleRelativistic::_LINEAR_SAMPLING_SCALE_;

vector<vector<int> > PIC::EnergyDistributionSampleRelativistic::CombinedSpeciesDistributionTable;

//====================================================
//init the sampling buffers
void PIC::EnergyDistributionSampleRelativistic::Init() {
  int idim,iProbe,i,j,k,offset;

  if (SamplingInitializedFlag==true) exit(__LINE__,__FILE__,"Error: EnergyDistributionSampleRelativistic is already initialized");

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif


//  nSamleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lenfths of the sampling intervals
  double t0,t1; //tempotary variables to satisfy intel c++ compiler

  log10eMin=log10(eMin);
  log10eMax=log10(eMax);
  log10dE=(log10eMax-log10eMin)/nSampledFunctionPoints;

  dE=(eMax-eMin)/nSampledFunctionPoints;

  //de-allocate sampling buffers if they were allocated before
  if (SampleLocalCellNumber!=NULL) {
    delete [] SampleLocalCellNumber;
    SampleLocalCellNumber=NULL;
  }

  if (SampleNodes!=NULL) {
    delete [] SampleNodes;
    SampleNodes=NULL;
  }

  if (SamplingBuffer!=NULL) {
    delete [] SamplingBuffer[0][0];
    delete [] SamplingBuffer[0];
    delete [] SamplingBuffer;

    SamplingBuffer=NULL;
  }


  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nSamleLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nSamleLocations];


  //access to the sampling buffer : SamplingBuffer[iLocation][spec][iDistributionFunctionPoint]
  SamplingBuffer=new double** [nSamleLocations];
  SamplingBuffer[0]=new double* [nSamleLocations*PIC::nTotalSpecies];

  for (i=0,offset=0;i<nSamleLocations;i++) {
    SamplingBuffer[i]=SamplingBuffer[0]+offset;
    offset+=PIC::nTotalSpecies;
  }

  SamplingBuffer[0][0]=new double [nSamleLocations*PIC::nTotalSpecies*nSampledFunctionPoints];

  for (i=0,offset=0;i<nSamleLocations;i++) for (j=0;j<PIC::nTotalSpecies;j++) {
    SamplingBuffer[i][j]=SamplingBuffer[0][0]+offset;
    offset+=nSampledFunctionPoints;
  }


  //init the sampling informations
  for (iProbe=0;iProbe<nSamleLocations;iProbe++) {
    SampleNodes[iProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[iProbe]);
    if (SampleNodes[iProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[iProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[iProbe],i,j,k,SampleNodes[iProbe],false);
    if (SampleLocalCellNumber[iProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void PIC::EnergyDistributionSampleRelativistic::flushSamplingBuffers() {
  long int i,TotalDataLength=nSamleLocations*PIC::nTotalSpecies*nSampledFunctionPoints;
  double *ptr=SamplingBuffer[0][0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}

//====================================================
//Sample the distribution function
void PIC::EnergyDistributionSampleRelativistic::SampleDistributionFnction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,iProbe,spec,idim;
  double LocalParticleWeight,e,mass,speed,CellMeasure;

  for (node=SampleNodes[0],iProbe=0;iProbe<nSamleLocations;node=SampleNodes[++iProbe]) if (node->Thread==PIC::ThisThread) {
      double *v;
      PIC::ParticleBuffer::byte *ParticleData;
      int i,j,k;
      PIC::Mesh::cDataCenterNode *cell;

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[iProbe],i,j,k);
      cell=node->block->GetCenterNode(SampleLocalCellNumber[iProbe]);
      CellMeasure=cell->Measure;
      if (CellMeasure<=0.0) CellMeasure=1.0;

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
          if ((i>=0)&&(i<nSampledFunctionPoints)) SamplingBuffer[iProbe][spec][i]+=LocalParticleWeight/CellMeasure;

          break;
        case _LOGARITHMIC_SAMPLING_SCALE_:
          speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
          e=Relativistic::Speed2E(speed,mass)*J2eV;

          i=(int)((log10(e)-log10eMin)/log10dE);
          if ((i>=0)&&(i<nSampledFunctionPoints)) SamplingBuffer[iProbe][spec][i]+=LocalParticleWeight/CellMeasure;

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
  long int iProbe,idim,i,j,thread,s;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,e;
  char str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);

  //temporary sampling buffer
  double *tempSamplingBuffer=new double [nSampledFunctionPoints];

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
  for (iProbe=0;iProbe<nSamleLocations;iProbe++) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      sprintf(str,"%s.nSamplePoint=%ld.dat",fname,iProbe);
      fout=fopen(str,"w");

      fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);

      fprintf(fout,"\"TITLE=Distribution function at x=%e",SamplingLocations[iProbe][0]);
      for (idim=1;idim<DIM;idim++) fprintf(fout,", %e",SamplingLocations[iProbe][idim]);

      fprintf(fout,"\"\nVARIABLES=\"E [eV]\", \"v[m/s]\", \"f(E)\", \"Fluency(E)\"\n");

      //init the tempSamplingBuffer
      for (i=0;i<nSampledFunctionPoints;i++) tempSamplingBuffer[i]=0.0;

      //collect the sampled information from other processors
      int iCombinedSpecies,tempOffsetSpec;

      for (iCombinedSpecies=0;iCombinedSpecies<CombineDistributionSpecies.size();iCombinedSpecies++) {
        s=CombineDistributionSpecies[iCombinedSpecies];

        for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
          for (i=0;i<nSampledFunctionPoints;i++) {
            if (thread==0) {
              tempSamplingBuffer[i]+=SamplingBuffer[iProbe][s][i];
            }
            else {
              tempSamplingBuffer[i]+=pipe.recv<double>(thread);
            }
          }
        }

      }

      //normalize the energy distribution functions
      norm=0.0;

      for (i=0;i<nSampledFunctionPoints;i++) {
        switch (EnergySamplingMode) {
        case _LINEAR_SAMPLING_SCALE_:
          norm+=tempSamplingBuffer[i]*dE;

          break;
        case _LOGARITHMIC_SAMPLING_SCALE_:
          norm+=tempSamplingBuffer[i]*eMin*pow(10.0,i*log10dE)*(pow(10.0,log10dE)-1.0);

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is unknown");
        }
      }

      //check the normalization constant
      norm=(fabs(norm)>0.0) ? fabs(norm) : 1.0;

//      if (fabs(norm)>0.0) for (i=0;i<nSampledFunctionPoints;i++) tempSamplingBuffer[i+tempOffsetSpec]/=norm;

      //print the output file
      for (i=0;i<nSampledFunctionPoints;i++) {
        switch (EnergySamplingMode) {
        case _LINEAR_SAMPLING_SCALE_:
          e=eMin+(i+0.5)*dE;

          break;
        case _LOGARITHMIC_SAMPLING_SCALE_:
          e=eMin*pow(10.0,(i+0.5)*log10dE);

          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is unknown");
         }

        double Speed=Relativistic::E2Speed(e*eV2J,PIC::MolecularData::GetMass(spec));
        fprintf(fout,"%e  %e  %e  %e\n",e,Speed,tempSamplingBuffer[i]/norm,tempSamplingBuffer[i]/PIC::LastSampleLength*Speed);
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

        for (i=0;i<nSampledFunctionPoints;i++) {
          pipe.send(SamplingBuffer[iProbe][s][i]);
        }
      }

    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //deallocate the temporary buffer
  delete [] tempSamplingBuffer;
}


void PIC::EnergyDistributionSampleRelativistic::AddCombinedCombinedParticleDistributionList(vector<int> CombinedSpecies) {
  CombinedSpeciesDistributionTable.push_back(CombinedSpecies);
}
