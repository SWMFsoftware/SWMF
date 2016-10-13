//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions that model colecular collisions with the background atmosphere

#include "pic.h"

#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_

int PIC::MolecularCollisions::BackgroundAtmosphere::LocalTotalCollisionFreqSamplingOffset=-1;
int PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset=-1;

double **PIC::MolecularCollisions::BackgroundAtmosphere::TotalCollisionModelParticleSourceRate=NULL;
double **PIC::MolecularCollisions::BackgroundAtmosphere::TotalCollisionModelParticleLossRate=NULL;

//init the model
void PIC::MolecularCollisions::BackgroundAtmosphere::Init_BeforeParser() {
  //set up the output of the model
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
  PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);
}

void PIC::MolecularCollisions::BackgroundAtmosphere::Init_AfterParser() {
  int i,spec;

  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SampleModelData,OutputSampledModelData);

  //set up buffers for sampling the model particle spource and loss rates
  TotalCollisionModelParticleSourceRate=new double* [PIC::nTotalSpecies];
  TotalCollisionModelParticleSourceRate[0]=new double[PIC::nTotalSpecies*PIC::nTotalThreadsOpenMP];

  TotalCollisionModelParticleLossRate=new double* [PIC::nTotalSpecies];
  TotalCollisionModelParticleLossRate[0]=new double[PIC::nTotalSpecies*PIC::nTotalThreadsOpenMP];

  for (spec=1;spec<PIC::nTotalSpecies;spec++) {
    TotalCollisionModelParticleSourceRate[spec]=TotalCollisionModelParticleSourceRate[0]+spec*PIC::nTotalThreadsOpenMP;
    TotalCollisionModelParticleLossRate[spec]=TotalCollisionModelParticleLossRate[0]+spec*PIC::nTotalThreadsOpenMP;
  }

  for (i=0;i<PIC::nTotalSpecies*PIC::nTotalThreadsOpenMP;i++) {
    TotalCollisionModelParticleSourceRate[0][i]=0.0;
    TotalCollisionModelParticleLossRate[0][i]=0.0;
  }

}

//output sampled model paramters
void PIC::MolecularCollisions::BackgroundAtmosphere::SampleModelData() {
  //place holder for sampling rocedure for the background atmosphere model
}

void PIC::MolecularCollisions::BackgroundAtmosphere::OutputSampledModelData(int DataOutputFileNumber) {
  int spec,thread;

  //collect the thermalization and collision source rates
  double TotalCollsionLossRate[PIC::nTotalSpecies],TotalCollisonSourceRate[PIC::nTotalSpecies];
  double threadTotalCollsionLossRate[PIC::nTotalSpecies],threadTotalCollisonSourceRate[PIC::nTotalSpecies];

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    threadTotalCollsionLossRate[spec]=0.0;
    threadTotalCollisonSourceRate[spec]=0.0;

    for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      threadTotalCollsionLossRate[spec]+=TotalCollisionModelParticleLossRate[spec][thread];
      threadTotalCollisonSourceRate[spec]+=TotalCollisionModelParticleSourceRate[spec][thread];

      TotalCollisionModelParticleLossRate[spec][thread]=0.0;
      TotalCollisionModelParticleSourceRate[spec][thread]=0.0;
    }
  }

  MPI_Allreduce(threadTotalCollsionLossRate,TotalCollsionLossRate,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(threadTotalCollisonSourceRate,TotalCollisonSourceRate,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    cout << "$PREFIX:Background Atmosphere: \n$PREFIX: spec \t Collision Source Rate \t Collision Loss Rate \n";

    for (spec=0;spec<PIC::nTotalSpecies;spec++) cout << "$PREFIX:" << spec << "\t" << TotalCollisonSourceRate[spec]/PIC::LastSampleLength << "\t" << TotalCollsionLossRate[spec]/PIC::LastSampleLength << endl;
  }
}


//Request Sampling Buffers
int PIC::MolecularCollisions::BackgroundAtmosphere::RequestSamplingData(int offset) {
  int RequesterDataLength=0;

  LocalTotalCollisionFreqSamplingOffset=offset+RequesterDataLength;
  RequesterDataLength+=PIC::nTotalSpecies*sizeof(double);

  LocalEnergyTransferRateSamplingOffset=offset+RequesterDataLength;
  RequesterDataLength+=PIC::nTotalSpecies*sizeof(double);

  return RequesterDataLength;
}

void PIC::MolecularCollisions::BackgroundAtmosphere::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Total Background Atmosphere Collison Freq [per a real particle,per sec]\", \"Energy Exchange Rate [J/m^3/s]\"");
}

void PIC::MolecularCollisions::BackgroundAtmosphere::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

  struct cDataExchengeBuffer {
    double TotalCollisionFreq;
    double EnergyExchangeRate;
  };

  cDataExchengeBuffer buffer={0.0,0.0};

  if ((pipe->ThisThread==CenterNodeThread)&&(LocalTotalCollisionFreqSamplingOffset!=-1)) {
    buffer.TotalCollisionFreq=*(DataSetNumber+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
    buffer.EnergyExchangeRate=*(DataSetNumber+(double*)(LocalEnergyTransferRateSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv((char*)&buffer,sizeof(cDataExchengeBuffer),CenterNodeThread);
    }

    if (PIC::LastSampleLength!=0) {
      buffer.TotalCollisionFreq/=PIC::LastSampleLength;
      buffer.EnergyExchangeRate/=PIC::LastSampleLength;
    }

    fprintf(fout,"%e  %e ",buffer.TotalCollisionFreq,buffer.EnergyExchangeRate);
  }
  else {
    pipe->send((char*)&buffer,sizeof(cDataExchengeBuffer));
  }
}

void PIC::MolecularCollisions::BackgroundAtmosphere::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i,spec;
  double c;

  //exit in case sampling buffers are not allocated
  if (LocalTotalCollisionFreqSamplingOffset==-1) return;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //interpolate the local productions rate
    double InterpoaltedTotalCollisionFreq=0.0;
    double InterpolatedEnergyExchangeRate=0.0;

    //interpolate the sampled data
    for (i=0;i<nInterpolationCoeficients;i++) {
      c=InterpolationCoeficients[i];

      InterpoaltedTotalCollisionFreq+=c*(*(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
      InterpolatedEnergyExchangeRate+=c*(*(double*)(LocalEnergyTransferRateSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
    }

    //stored the interpolated data in the associated data buffer
    *(spec+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpoaltedTotalCollisionFreq;
    *(spec+(double*)(LocalEnergyTransferRateSamplingOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpolatedEnergyExchangeRate;
  }
}


void PIC::MolecularCollisions::BackgroundAtmosphere::CollisionProcessor() {
  int i,j,k,thread;

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  //sample collisino frequentcy
//  double TotalProspectiveCollisionParticleWeight[PIC::nTotalSpecies],TotalOccurringCollisionParticleWeight[PIC::nTotalSpecies];
//  double ModelParticleEnergyExchangeRate[PIC::nTotalSpecies];


  //init sampling buffers
  static double **TotalProspectiveCollisionParticleWeight=NULL,**TotalOccurringCollisionParticleWeight=NULL,**ModelParticleEnergyExchangeRate=NULL;

  if (TotalProspectiveCollisionParticleWeight==NULL) {
    TotalProspectiveCollisionParticleWeight=new double* [PIC::nTotalThreadsOpenMP];
    TotalOccurringCollisionParticleWeight=new double* [PIC::nTotalThreadsOpenMP];
    ModelParticleEnergyExchangeRate=new double* [PIC::nTotalThreadsOpenMP];

    TotalProspectiveCollisionParticleWeight[0]=new double [PIC::nTotalThreadsOpenMP*PIC::nTotalSpecies];
    TotalOccurringCollisionParticleWeight[0]=new double [PIC::nTotalThreadsOpenMP*PIC::nTotalSpecies];
    ModelParticleEnergyExchangeRate[0]=new double [PIC::nTotalThreadsOpenMP*PIC::nTotalSpecies];

    int s,offset=0;

    for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      TotalProspectiveCollisionParticleWeight[thread]=TotalProspectiveCollisionParticleWeight[0]+offset;
      TotalOccurringCollisionParticleWeight[thread]=TotalOccurringCollisionParticleWeight[0]+offset;
      ModelParticleEnergyExchangeRate[thread]=ModelParticleEnergyExchangeRate[0]+offset;

      offset+=PIC::nTotalSpecies;
    }
  }


  int LocalCellNumber;

  static int ParticleBufferLengthDefaut=50000;

  struct cCollidingParticleList {
    long int Particle;
    double CollisionTimeFraction;
  };



  if (initTableFlag==false) {
    //init thr center node table
    nTotalCenterNodes=0,initTableFlag=true;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataCenterNode *cell;

  //the buffer of particles that occuping the local cell
  long int modelParticle;
  int BackgroundSpecieNumber,spec,idim;
  PIC::ParticleBuffer::byte *modelParticleData;
  double vModelParticle[3],xModelParticle[3],vBackgroundParticle[3],particleCollisionTime,cr2;
  PIC::Mesh::cDataBlockAMR *block;

  //the temporary particle that represents the background atmosphere
  static long int *tempBackgroundAtmosphereParticle=NULL;
  static PIC::ParticleBuffer::byte **BackgroundAtmosphereParticleData=NULL;

  if (tempBackgroundAtmosphereParticle==NULL) {
    tempBackgroundAtmosphereParticle=new long int [PIC::nTotalThreadsOpenMP];
    BackgroundAtmosphereParticleData=new PIC::ParticleBuffer::byte *[PIC::nTotalThreadsOpenMP];
  }

  for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
    tempBackgroundAtmosphereParticle[thread]=PIC::ParticleBuffer::GetNewParticle();
    BackgroundAtmosphereParticleData[thread]=PIC::ParticleBuffer::GetParticleDataPointer(tempBackgroundAtmosphereParticle[thread]);
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)BackgroundAtmosphereParticleData[thread]);
  }

  static cCollidingParticleList **CollidingParticleList=NULL;
  static int *ParticleBufferLength=NULL;

  if (ParticleBufferLength==NULL) {
    ParticleBufferLength=new int [PIC::nTotalThreadsOpenMP];
    CollidingParticleList=new cCollidingParticleList* [PIC::nTotalThreadsOpenMP];

    for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      ParticleBufferLength[thread]=ParticleBufferLengthDefaut;
      CollidingParticleList[thread]=new cCollidingParticleList [ParticleBufferLength[thread]];
    }
  }

  //sample the processor load
//#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
//#endif

  //loop through all nodes and cells on the currect processor
  double MajorantCollisionFreq,SigmaCr,SigmaCrMax=0.0;
  double timeCounter,localTimeStep,TranslationalEnergy;
  double massModelParticle,massBackgroundParticle,Vrel[3]={0.0,0.0,0.0},Vcm[3]={0.0,0.0,0.0},am;


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    //reset the balancing counters
    for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      node=DomainBlockDecomposition::BlockTable[nLocalNode];
      if (node->block!=NULL) *(thread+(double*)(node->block->GetAssociatedDataBufferPointer()+PIC::Mesh::cDataBlockAMR::LoadBalancingMeasureOffset))=0.0;
    }
#endif //_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

#if _PIC__OPENMP_THREAD_SPLIT_MODE_ == _PIC__OPENMP_THREAD_SPLIT_MODE__BLOCKS_
#pragma omp parallel for schedule(dynamic,_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_) default (none) firstprivate (LocalCellNumber, \
      node,cell,modelParticle, BackgroundSpecieNumber,spec,idim,modelParticleData,vModelParticle,xModelParticle, \
      vBackgroundParticle,particleCollisionTime,cr2,block,tempBackgroundAtmosphereParticle,MajorantCollisionFreq,SigmaCr,SigmaCrMax,timeCounter,localTimeStep, \
      TranslationalEnergy,massModelParticle,massBackgroundParticle,Vrel,BackgroundAtmosphereParticleData,Vcm,am,EndTime,StartTime,thread) \
     \
    shared(centerNodeIndexTable_Glabal,nTotalCenterNodes,TotalProspectiveCollisionParticleWeight,TotalOccurringCollisionParticleWeight,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh,ModelParticleEnergyExchangeRate, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalTotalCollisionFreqSamplingOffset,PIC::Mesh::collectingCellSampleDataPointerOffset, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset,CollidingParticleList,ParticleBufferLength)
#else
#pragma omp parallel for schedule(dynamic,1) default (none) firstprivate (LocalCellNumber, \
      node,cell,modelParticle, BackgroundSpecieNumber,spec,idim,modelParticleData,vModelParticle,xModelParticle, \
      vBackgroundParticle,particleCollisionTime,cr2,block,tempBackgroundAtmosphereParticle,MajorantCollisionFreq,SigmaCr,SigmaCrMax,timeCounter,localTimeStep, \
      TranslationalEnergy,massModelParticle,massBackgroundParticle,Vrel,BackgroundAtmosphereParticleData,Vcm,am,EndTime,StartTime,thread) \
     \
    shared(centerNodeIndexTable_Glabal,nTotalCenterNodes,TotalProspectiveCollisionParticleWeight,TotalOccurringCollisionParticleWeight,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh,ModelParticleEnergyExchangeRate, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalTotalCollisionFreqSamplingOffset,PIC::Mesh::collectingCellSampleDataPointerOffset, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset,CollidingParticleList,ParticleBufferLength)
#endif  // _PIC__OPENMP_THREAD_SPLIT_MODE_
#endif  //_COMPILATION_MODE_
  for (int CellCounter=0;CellCounter<DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;CellCounter++) {
    int nLocalNode,ii=CellCounter;
    int kCell,jCell,iCell; //,i,j,k;
    long int nCollidingParticles=0;

    nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    kCell=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=kCell*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    jCell=ii/_BLOCK_CELLS_X_;
    ii-=jCell*_BLOCK_CELLS_X_;

    iCell=ii;

    StartTime=MPI_Wtime();
    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    block=node->block;
    if (block==NULL) continue;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    thread=omp_get_thread_num();
    #else
    thread=0;
    #endif



    {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
          cell=block->GetCenterNode(LocalCellNumber);
          if (cell==NULL) continue;

          modelParticle=block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];
          nCollidingParticles=0;



      for (spec=0;spec<PIC::nTotalSpecies;spec++) TotalProspectiveCollisionParticleWeight[thread][spec]=0.0,TotalOccurringCollisionParticleWeight[thread][spec]=0.0,ModelParticleEnergyExchangeRate[thread][spec]=0.0;


      while (modelParticle!=-1) {
        if (nCollidingParticles==ParticleBufferLength[thread]) {
          //exit(__LINE__,__FILE__,"Error: the value of 'ParticleBufferLength' is exeeded - too many particles in a cells. Increase the value of 'ParticleBufferLength'");

          //re-allocate the list of the colliding particles
          cCollidingParticleList *tmpCollidingParticleList=new cCollidingParticleList [2*ParticleBufferLength[thread]];
          memcpy(tmpCollidingParticleList,CollidingParticleList[thread],ParticleBufferLength[thread]*sizeof(cCollidingParticleList));
          delete [] CollidingParticleList[thread];

          CollidingParticleList[thread]=tmpCollidingParticleList;
          ParticleBufferLength[thread]*=2;
        }

        CollidingParticleList[thread][nCollidingParticles].Particle=modelParticle;
        CollidingParticleList[thread][nCollidingParticles].CollisionTimeFraction=1.0;

        ++nCollidingParticles;
        modelParticle=PIC::ParticleBuffer::GetNext(modelParticle);
      }

_StartParticleCollisionLoop_:
      while (nCollidingParticles!=0) {
        modelParticle=CollidingParticleList[thread][0].Particle;
        modelParticleData=PIC::ParticleBuffer::GetParticleDataPointer(modelParticle);

        spec=PIC::ParticleBuffer::GetI(modelParticleData);
        PIC::ParticleBuffer::GetV(vModelParticle,modelParticleData);
        PIC::ParticleBuffer::GetX(xModelParticle,modelParticleData);

        massModelParticle=PIC::MolecularData::GetMass(spec);
        localTimeStep=node->block->GetLocalTimeStep(spec);
        particleCollisionTime=CollidingParticleList[thread][0].CollisionTimeFraction*localTimeStep;

        TotalProspectiveCollisionParticleWeight[thread][spec]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

        //exclude the particle from the list of the model particle that still needs to be collided
        if (nCollidingParticles!=1) CollidingParticleList[thread][0]=CollidingParticleList[thread][nCollidingParticles-1];
        --nCollidingParticles;

        //simulate collisions with the background atmosphere
        for (BackgroundSpecieNumber=0;BackgroundSpecieNumber<nTotalBackgroundSpecies;BackgroundSpecieNumber++) {
          SigmaCrMax=GetSigmaCrMax(spec,BackgroundSpecieNumber,modelParticleData);
          MajorantCollisionFreq=GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,xModelParticle)*SigmaCrMax;


//TEST!!!!!!!!!!!!!!!
//MajorantCollisionFreq=(sqrt(xModelParticle[0]*xModelParticle[0]+xModelParticle[1]*xModelParticle[1]+xModelParticle[2]*xModelParticle[2])<700.0E3+_RADIUS_(_TARGET_)) ? 4E5*1.0E6*SigmaCrMax : 0.0;



#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_
          //if the majorant frequentcy large, reivalute it with the mean value of the background density calcualted along the trajectory of the particle with the sped xdMeanDensity;
          if (false) { //(MajorantCollisionFreq*particleCollisionTime>50.0) {
            double xProbe[3],meanValueOld,meanValueNew;
            int iDensityAveragingLevel;

            const int maxDensityAveragingLevels=6;


            for (idim=0;idim<DIM;idim++) xProbe[idim]=xModelParticle[idim]-vModelParticle[idim]*particleCollisionTime;
            meanValueOld=0.5*(GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,xModelParticle)+GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,xProbe));


            for (iDensityAveragingLevel=1;iDensityAveragingLevel<=maxDensityAveragingLevels;iDensityAveragingLevel++) {
              int idim,i,nTestPoints;

              nTestPoints=(1<<iDensityAveragingLevel);
              meanValueNew=meanValueOld*(1+(1<<(iDensityAveragingLevel-1)));

              for (i=1;i<=nTestPoints;i+=2) {
                for (idim=0;idim<DIM;idim++) xProbe[idim]=xModelParticle[idim]-vModelParticle[idim]*double(i)*particleCollisionTime/nTestPoints;     //'-' because the avaraging clong the 'past' trajectory

                meanValueNew+=GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,xProbe);
              }

              meanValueNew/=(1+nTestPoints);
              if (fabs(1.0-meanValueOld/meanValueNew)<0.2) break;
            }

            MajorantCollisionFreq=meanValueNew*SigmaCrMax;
          }
#endif


          timeCounter=0.0;

          massBackgroundParticle=BackgroundSpeciesMassTable[BackgroundSpecieNumber];
          am=massModelParticle+massBackgroundParticle;

          if (MajorantCollisionFreq>0.0) while ((timeCounter-=log(rnd())/MajorantCollisionFreq)<particleCollisionTime) {
            //generate particle that represents the background atmosphere and calcualte the cross section
            GenerateBackgoundAtmosphereParticle(BackgroundAtmosphereParticleData[thread],BackgroundSpecieNumber,cell,node);
            PIC::ParticleBuffer::GetV(vBackgroundParticle,BackgroundAtmosphereParticleData[thread]);

//TEST!!!!!!!!!!!!!!!
//for (idim=0;idim<3;idim++) vBackgroundParticle[idim]=0.0;

            for (idim=0,cr2=0.0;idim<3;idim++) {
              Vrel[idim]=vModelParticle[idim]-vBackgroundParticle[idim];
              Vcm[idim]=(massModelParticle*vModelParticle[idim]+massBackgroundParticle*vBackgroundParticle[idim])/am;

              cr2+=pow(Vrel[idim],2);
            }

            TranslationalEnergy=0.5*massModelParticle*massBackgroundParticle/(massModelParticle+massBackgroundParticle)*cr2;
            SigmaCr=GetCollisionCrossSectionBackgoundAtmosphereParticle(spec,BackgroundSpecieNumber,modelParticleData,BackgroundAtmosphereParticleData[thread],TranslationalEnergy,cr2)*sqrt(cr2);

            //check if the collision is possible
            if (rnd()*SigmaCrMax>SigmaCr) continue;

            //sample collision frecuentcy
            TotalOccurringCollisionParticleWeight[thread][spec]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

            //redistribute the relative velocity of the collided particles
            double Vrc,V[3];
            double CosKsi,SinKsi,CosEps,SinEps,D,c;

#if _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_
            CosKsi=2.0*rnd()-1.0;
            SinKsi=sqrt(1.0-CosKsi*CosKsi);
#elif _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ == _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__USER_DEFINED_
            double VelocityScatteringAngle;

            VelocityScatteringAngle=GetCollisionScatteringAngle(Vrel,TranslationalEnergy,spec,BackgroundSpecieNumber);
            CosKsi=cos(VelocityScatteringAngle);
            SinKsi=sin(VelocityScatteringAngle);
#else
            exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

            c=2*Pi*rnd();
            SinEps=sin(c);
            CosEps=cos(c);

            D=sqrt(Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

            if (D>1.0E-6) {
              Vrc=sqrt(Vrel[0]*Vrel[0]+Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);

              V[0]=CosKsi*Vrel[0]+SinKsi*SinEps*D;
              V[1]=CosKsi*Vrel[1]+SinKsi*(Vrc*Vrel[2]*CosEps-Vrel[0]*Vrel[1]*SinEps)/D;
              V[2]=CosKsi*Vrel[2]-SinKsi*(Vrc*Vrel[1]*CosEps+Vrel[0]*Vrel[2]*SinEps)/D;
            }
            else {
              V[0]=CosKsi*Vrel[0];
              V[1]=SinKsi*CosEps*Vrel[0];
              V[2]=SinKsi*SinEps*Vrel[0];
            }

            Vrel[0]=V[0];
            Vrel[1]=V[1];
            Vrel[2]=V[2];

            //the collision between the model particle and the particle from the background atmosphere has occured
            double InitialModelPArticleEnergy=0.0,FinalModelParticleEnergy=0.0;

            for (idim=0;idim<3;idim++) {
              InitialModelPArticleEnergy+=pow(vModelParticle[idim],2);

              vModelParticle[idim]=Vcm[idim]+massBackgroundParticle/am*Vrel[idim];
              vBackgroundParticle[idim]=Vcm[idim]-massModelParticle/am*Vrel[idim];

              FinalModelParticleEnergy+=pow(vModelParticle[idim],2);
            }

            //sample energy exchabge rate
            ModelParticleEnergyExchangeRate[thread][spec]+=(FinalModelParticleEnergy-InitialModelPArticleEnergy)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

            //update velocities of the particles
            PIC::ParticleBuffer::SetV(vBackgroundParticle,BackgroundAtmosphereParticleData[thread]);
            PIC::ParticleBuffer::SetV(vModelParticle,modelParticle);

            //check if the 'background' particle should be kept in the system
#if _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_PARTICLE_ACCEPTANCE_MODE_ == _PIC_MODE_ON_
            if (Background2ModelSpeciesConversionTable[BackgroundSpecieNumber]!=-1) {
              PIC::ParticleBuffer::SetI(Background2ModelSpeciesConversionTable[BackgroundSpecieNumber],BackgroundAtmosphereParticleData[thread]);

              if (KeepConditionModelParticle(BackgroundAtmosphereParticleData[thread])==true) {
                long int newParticle;
                PIC::ParticleBuffer::byte *newParticleData;

                newParticle=PIC::ParticleBuffer::GetNewParticle(block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);
                PIC::ParticleBuffer::CloneParticle(newParticle,tempBackgroundAtmosphereParticle[thread]);
                newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
                PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)newParticleData);

                //set the position of the new particle to be the position of the original model particle
                PIC::ParticleBuffer::SetX(xModelParticle,newParticleData);

                //add the new particle to the list of the particles that needes to be collided
                if (nCollidingParticles>ParticleBufferLength[thread]-1) {
                  //exit(__LINE__,__FILE__,"Error: the particle buffer is overflown");

                  //re-allocate the list of the colliding particles
                  cCollidingParticleList *tmpCollidingParticleList=new cCollidingParticleList [2*ParticleBufferLength[thread]];
                  memcpy(tmpCollidingParticleList,CollidingParticleList[thread],ParticleBufferLength[thread]*sizeof(cCollidingParticleList));
                  delete [] CollidingParticleList[thread];

                  CollidingParticleList[thread]=tmpCollidingParticleList;
                  ParticleBufferLength[thread]*=2;
                }

                CollidingParticleList[thread][nCollidingParticles].Particle=newParticle;
                CollidingParticleList[thread][nCollidingParticles].CollisionTimeFraction=1.0-timeCounter/localTimeStep;
                ++nCollidingParticles;

                //accout for the source of new exospheric particles
                TotalCollisionModelParticleSourceRate[Background2ModelSpeciesConversionTable[BackgroundSpecieNumber]][thread]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)*block->GetLocalParticleWeight(spec)/localTimeStep;

                //set particle weight
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
                double wcorrection;
                int bspec;

                bspec=PIC::ParticleBuffer::GetI(newParticleData);

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
                if (bspec>=PIC::nTotalSpecies) exit(__LINE__,__FILE__,"Error: the species number exeeds the total number of species");
#endif

                wcorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)*block->GetLocalParticleWeight(spec)/block->GetLocalParticleWeight(bspec);
                PIC::ParticleBuffer::SetIndividualStatWeightCorrection(wcorrection,newParticleData);
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
        exit(__LINE__,__FILE__,"Error: not implementd");
#else
               exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif
              }
            }
#endif

            //check if the model particle should be removed from the system
#if _PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_ == _PIC_MODE_ON_
            if (KeepConditionModelParticle(modelParticleData)==false) {
              //the particle should be removed
              long int next,prev;

              //accout for the termalization of the exospheric particles
              TotalCollisionModelParticleLossRate[spec][thread]+=node->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)/localTimeStep;

              next=PIC::ParticleBuffer::GetNext(modelParticleData);
              prev=PIC::ParticleBuffer::GetPrev(modelParticleData);

              //reconnect particles from the list
              if (prev==-1) block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=next;
              else PIC::ParticleBuffer::SetNext(next,prev);

              if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

              PIC::ParticleBuffer::DeleteParticle(modelParticle);
              goto _StartParticleCollisionLoop_;
            }
#endif
          }
        }

      }

      //sample total collision frequentcy
      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        if (TotalProspectiveCollisionParticleWeight[thread][spec]>0.0) {
          *(spec+(double*)(LocalTotalCollisionFreqSamplingOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=
            TotalOccurringCollisionParticleWeight[thread][spec]/TotalProspectiveCollisionParticleWeight[thread][spec]/node->block->GetLocalTimeStep(spec);
        }

        *(spec+(double*)(LocalEnergyTransferRateSamplingOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=
            ModelParticleEnergyExchangeRate[thread][spec]*node->block->GetLocalParticleWeight(spec)/node->block->GetLocalTimeStep(spec)/cell->Measure*PIC::MolecularData::GetMass(spec)/2.0;
      }

    }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

//    node=node->nextNodeThisThread;
  }

  //delete the temporary particle representing the background atmosphere
  for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) PIC::ParticleBuffer::DeleteParticle(tempBackgroundAtmosphereParticle[thread]);
}


void PIC::MolecularCollisions::BackgroundAtmosphere::RemoveThermalBackgroundParticles() {

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;
  long int next,prev,ptr,LocalCellNumber;
  PIC::ParticleBuffer::byte *pdata;

  if (initTableFlag==false) {
    nTotalCenterNodes=0,initTableFlag=true;
    int i,j,k;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];


  //filter particles
  while (node!=NULL) {
/*    for (centerNodeIndexCounter=0;centerNodeIndexCounter<nTotalCenterNodes;centerNodeIndexCounter++) {
      //sort the particle from the cell
      LocalCellNumber=centerNodeIndexTable[centerNodeIndexCounter];
      ptr=node->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;*/

    for (int kCell=0;kCell<_BLOCK_CELLS_Z_;kCell++)
       for (int jCell=0;jCell<_BLOCK_CELLS_Y_;jCell++)
         for (int iCell=0;iCell<_BLOCK_CELLS_X_;iCell++) {

           LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
           ptr=node->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];

      while (ptr!=-1) {
        pdata=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

        next=PIC::ParticleBuffer::GetNext(pdata);
        prev=PIC::ParticleBuffer::GetPrev(pdata);

        #if _PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_ == _PIC_MODE_ON_
        if (KeepConditionModelParticle(pdata)==false) {
          //the particle should be removed
          //reconnect particles from the list
          if (prev==-1) node->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=next;
          else PIC::ParticleBuffer::SetNext(next,prev);

          if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

          PIC::ParticleBuffer::DeleteParticle(ptr);
        }
        #endif

        ptr=next;
      }
    }

    node=node->nextNodeThisThread;
  }
}


//model decrease of energy of the simulated particles in the stopping power approxymation
void PIC::MolecularCollisions::BackgroundAtmosphere::StoppingPowerProcessor() {
  int i,j,k,thread;

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

  //init sampling buffers
  static double **ModelParticleEnergyExchangeRate=NULL;

  if (ModelParticleEnergyExchangeRate==NULL) {
    ModelParticleEnergyExchangeRate=new double* [PIC::nTotalThreadsOpenMP];
    ModelParticleEnergyExchangeRate[0]=new double [PIC::nTotalThreadsOpenMP*PIC::nTotalSpecies];

    int s,offset=0;

    for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      ModelParticleEnergyExchangeRate[thread]=ModelParticleEnergyExchangeRate[0]+offset;

      offset+=PIC::nTotalSpecies;
    }
  }


  int LocalCellNumber;


  if (initTableFlag==false) {
    //init thr center node table
    nTotalCenterNodes=0,initTableFlag=true;

#if DIM == 3
    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 2
    for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#elif DIM == 1
    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      centerNodeIndexTable_Glabal[nTotalCenterNodes++]=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    }
#endif
  }

  int centerNodeIndexTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

  memcpy(centerNodeIndexTable,centerNodeIndexTable_Glabal,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int));

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataCenterNode *cell;

  //the buffer of particles that occuping the local cell
  long int modelParticle;
  int BackgroundSpecieNumber,spec,idim;
  PIC::ParticleBuffer::byte *modelParticleData;
  double vModelParticle[3],xModelParticle[3],vBackgroundParticle[3],particleCollisionTime,cr2;
  PIC::Mesh::cDataBlockAMR *block;

  //sample the processor load
//#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
//#endif


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    //reset the balancing counters
    for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      node=DomainBlockDecomposition::BlockTable[nLocalNode];
      if (node->block!=NULL) *(thread+(double*)(node->block->GetAssociatedDataBufferPointer()+PIC::Mesh::cDataBlockAMR::LoadBalancingMeasureOffset))=0.0;
    }
#endif //_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

#if _PIC__OPENMP_THREAD_SPLIT_MODE_ == _PIC__OPENMP_THREAD_SPLIT_MODE__BLOCKS_
#pragma omp parallel for schedule(dynamic,_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_) default (none) firstprivate (LocalCellNumber, \
      node,cell,modelParticle, BackgroundSpecieNumber,spec,idim,modelParticleData,vModelParticle,xModelParticle, \
      block,localTimeStep,EndTime,StartTime,thread) \
     \
    shared(centerNodeIndexTable_Glabal,nTotalCenterNodes,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh,ModelParticleEnergyExchangeRate, \
      PIC::Mesh::collectingCellSampleDataPointerOffset, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset)
#else
#pragma omp parallel for schedule(dynamic,1) default (none) firstprivate (LocalCellNumber, \
      node,cell,modelParticle, BackgroundSpecieNumber,spec,idim,modelParticleData,vModelParticle,xModelParticle, block, \
      EndTime,StartTime,thread) \
     \
    shared(centerNodeIndexTable_Glabal,nTotalCenterNodes,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh,ModelParticleEnergyExchangeRate, \
      PIC::Mesh::collectingCellSampleDataPointerOffset, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset)
#endif  // _PIC__OPENMP_THREAD_SPLIT_MODE_
#endif  //_COMPILATION_MODE_
  for (int CellCounter=0;CellCounter<DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;CellCounter++) {
    int nLocalNode,ii=CellCounter;
    int kCell,jCell,iCell; //,i,j,k;
    long int nCollidingParticles=0;

    nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    kCell=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=kCell*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    jCell=ii/_BLOCK_CELLS_X_;
    ii-=jCell*_BLOCK_CELLS_X_;

    iCell=ii;

    StartTime=MPI_Wtime();
    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    block=node->block;
    if (block==NULL) continue;

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    thread=omp_get_thread_num();
    #else
    thread=0;
    #endif



    {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
          cell=block->GetCenterNode(LocalCellNumber);
          if (cell==NULL) continue;

          modelParticle=block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)];
          nCollidingParticles=0;



      for (spec=0;spec<PIC::nTotalSpecies;spec++) ModelParticleEnergyExchangeRate[thread][spec]=0.0;

      while (modelParticle!=-1) {
        double v[3],x[3],StoppingPower;
        int spec,idim;

        modelParticleData=PIC::ParticleBuffer::GetParticleDataPointer(modelParticle);
        modelParticle=PIC::ParticleBuffer::GetNext(modelParticleData);

        spec=PIC::ParticleBuffer::GetI(modelParticleData);
        PIC::ParticleBuffer::GetV(v,modelParticleData);
        PIC::ParticleBuffer::GetX(x,modelParticleData);

        #if _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_ == _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL__STOPPING_POWER_
        StoppingPower=GetStoppingPower(x,v,spec,cell,node);
        #else
        exit(__LINE__,__FILE__,"Error: something is wrong. This function should not be called when _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_ != _PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL__STOPPING_POWER_");
        #endif //_PIC_BACKGROUND_ATMOSPHERE__COLLISION_MODEL_

        //get the change of the particle energy and the new velocity vector
        double l,dE,speed2,SpeedOld,SpeedNew,c0;

        speed2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
        l=sqrt(speed2)*node->block->GetLocalTimeStep(spec);
        dE=StoppingPower*l;

        SpeedOld=sqrt(speed2);
        SpeedNew=speed2-dE*2.0/PIC::MolecularData::GetMass(spec);
        if (SpeedNew<0.0) SpeedNew=0.0;
        SpeedNew=sqrt(SpeedNew);

        for (idim=0,c0=SpeedNew/SpeedOld;idim<3;idim++) v[idim]*=c0;

        //save the velocity vector with in the particle data vector
        PIC::ParticleBuffer::SetV(v,modelParticleData);

        //sample the ebergy exchange rate
        //sample energy exchabge rate
        ModelParticleEnergyExchangeRate[thread][spec]+=(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]-speed2)*
            PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData);

        //check wether a particle with sucj energy needs to stay in the system
        //check if the model particle should be removed from the system
#if _PIC_BACKGROUND_ATMOSPHERE__MODEL_PARTICLE_REMOVAL_MODE_ == _PIC_MODE_ON_
        if (KeepConditionModelParticle(modelParticleData)==false) {
          //the particle should be removed
          long int next,prev;

          //accout for the termalization of the exospheric particles
          TotalCollisionModelParticleLossRate[spec][thread]+=node->block->GetLocalParticleWeight(spec)*
              PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)/node->block->GetLocalTimeStep(spec);

          next=PIC::ParticleBuffer::GetNext(modelParticleData);
          prev=PIC::ParticleBuffer::GetPrev(modelParticleData);

          //reconnect particles from the list
          if (prev==-1) block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]=next;
          else PIC::ParticleBuffer::SetNext(next,prev);

          if (next!=-1) PIC::ParticleBuffer::SetPrev(prev,next);

          PIC::ParticleBuffer::DeleteParticle(modelParticle);
          continue;
        }
#endif
      }

      //sample total energy exchange rate
      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        *(spec+(double*)(LocalEnergyTransferRateSamplingOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=
            ModelParticleEnergyExchangeRate[thread][spec]*node->block->GetLocalParticleWeight(spec)/node->block->GetLocalTimeStep(spec)/cell->Measure*PIC::MolecularData::GetMass(spec)/2.0;
      }

    }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif

  }
}



#endif
