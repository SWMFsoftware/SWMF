//$Id$
//stopping power model of the energetic particle interactins with a background atmosphere

/*
 * pic_stopping_power.cpp
 *
 *  Created on: May 8, 2017
 *      Author: vtenishe
 */


#include "pic.h"

int PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRateOffset=-1;
double **PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRate=NULL;


//init the model
void PIC::MolecularCollisions::StoppingPowerModel::Init_BeforeParser() {
  //set up the output of the model
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
  PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);
}

void PIC::MolecularCollisions::StoppingPowerModel::Init_AfterParser() {
  int i,spec;

  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SampleModelData,OutputSampledModelData);

  //set up buffers for sampling the model particle spource and loss rates
  TotalModelParticleEnergyLossRate=new double* [PIC::nTotalSpecies];
  TotalModelParticleEnergyLossRate[0]=new double[PIC::nTotalSpecies*PIC::nTotalThreadsOpenMP];

  for (spec=1;spec<PIC::nTotalSpecies;spec++) {
    TotalModelParticleEnergyLossRate[spec]=TotalModelParticleEnergyLossRate[0]+spec*PIC::nTotalThreadsOpenMP;
  }

  for (i=0;i<PIC::nTotalSpecies*PIC::nTotalThreadsOpenMP;i++) {
    TotalModelParticleEnergyLossRate[0][i]=0.0;
  }
}

//output sampled model paramters
void PIC::MolecularCollisions::StoppingPowerModel::SampleModelData() {
  //place holder for sampling rocedure for the background atmosphere model
}

void PIC::MolecularCollisions::StoppingPowerModel::OutputSampledModelData(int DataOutputFileNumber) {
  int spec,thread;

  //collect the thermalization and collision source rates
  double threadTotalModelParticleEnergyLossRate[PIC::nTotalSpecies],localTotalModelParticleEnergyLossRate[PIC::nTotalSpecies];

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    threadTotalModelParticleEnergyLossRate[spec]=0.0;

    for (thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      threadTotalModelParticleEnergyLossRate[spec]+=TotalModelParticleEnergyLossRate[spec][thread];
      TotalModelParticleEnergyLossRate[spec][thread]=0.0;
    }
  }

  MPI_Allreduce(threadTotalModelParticleEnergyLossRate,localTotalModelParticleEnergyLossRate,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    cout << "$PREFIX:Background Atmosphere: \n$PREFIX: spec  \t The Total Energy Loss Rate \n";

    for (spec=0;spec<PIC::nTotalSpecies;spec++) cout << "$PREFIX:" << spec << "\t" << localTotalModelParticleEnergyLossRate[spec]/PIC::LastSampleLength << endl;
  }
}


//Request Sampling Buffers
int PIC::MolecularCollisions::StoppingPowerModel::RequestSamplingData(int offset) {
  int RequesterDataLength=0;

  TotalModelParticleEnergyLossRateOffset=offset+RequesterDataLength;
  RequesterDataLength+=PIC::nTotalSpecies*sizeof(double);

  return RequesterDataLength;
}

void PIC::MolecularCollisions::StoppingPowerModel::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Total Energy Loss Rate [J/m^3/s]\"");
}

void PIC::MolecularCollisions::StoppingPowerModel::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

  struct cDataExchengeBuffer {
    double EnergyLossRate;
  };

  cDataExchengeBuffer buffer={0.0};

  if (pipe->ThisThread==CenterNodeThread) {
    buffer.EnergyLossRate=*(DataSetNumber+(double*)(TotalModelParticleEnergyLossRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv((char*)&buffer,sizeof(cDataExchengeBuffer),CenterNodeThread);
    }

    if (PIC::LastSampleLength!=0) {
      buffer.EnergyLossRate/=PIC::LastSampleLength;
    }

    fprintf(fout,"%e ",buffer.EnergyLossRate);
  }
  else {
    pipe->send((char*)&buffer,sizeof(cDataExchengeBuffer));
  }
}

void PIC::MolecularCollisions::StoppingPowerModel::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i,spec;
  double c;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //interpolate the local productions rate
    double InterpolatedEnergyLossRate=0.0;

    //interpolate the sampled data
    for (i=0;i<nInterpolationCoeficients;i++) {
      c=InterpolationCoeficients[i];

      InterpolatedEnergyLossRate+=c*(*(double*)(TotalModelParticleEnergyLossRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
    }

    //stored the interpolated data in the associated data buffer
    *(spec+(double*)(TotalModelParticleEnergyLossRateOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CenterNode->GetAssociatedDataBufferPointer()))=InterpolatedEnergyLossRate;
  }
}



//model decrease of energy of the simulated particles in the stopping power approxymation
void PIC::MolecularCollisions::StoppingPowerModel::ModelProcessor() {
  int i,j,k,thread,LocalCellNumber;

  //the table of increments for accessing the cells in the block
  static bool initTableFlag=false;
  static int centerNodeIndexTable_Glabal[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
  static int nTotalCenterNodes=-1;

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
      block,EndTime,StartTime,thread) \
     \
    shared(PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRate,PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRateOffset,centerNodeIndexTable_Glabal,nTotalCenterNodes,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh, \
      PIC::Mesh::collectingCellSampleDataPointerOffset, \
      PIC::MolecularCollisions::BackgroundAtmosphere::LocalEnergyTransferRateSamplingOffset)
#else
#pragma omp parallel for schedule(dynamic,1) default (none) firstprivate (LocalCellNumber, \
      node,cell,modelParticle, BackgroundSpecieNumber,spec,idim,modelParticleData,vModelParticle,xModelParticle, block, \
      EndTime,StartTime,thread) \
     \
    shared(PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRate,PIC::MolecularCollisions::StoppingPowerModel::TotalModelParticleEnergyLossRateOffset,centerNodeIndexTable_Glabal,nTotalCenterNodes,centerNodeIndexTable, \
      PIC::DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable,PIC::Mesh::mesh, \
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

      while (modelParticle!=-1) {
        double v[3],x[3],StoppingPower;
        int spec,idim;

        modelParticleData=PIC::ParticleBuffer::GetParticleDataPointer(modelParticle);
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

        if (SpeedNew<0.0) {
          //the particle has been stopped -> remove it
          long int nextModelParticle;

          nextModelParticle=PIC::ParticleBuffer::GetNext(modelParticle);
          PIC::ParticleBuffer::DeleteParticle(modelParticle,block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);

          modelParticle=nextModelParticle;
          continue;
        }

        SpeedNew=sqrt(SpeedNew);
        for (idim=0,c0=SpeedNew/SpeedOld;idim<3;idim++) v[idim]*=c0;

        //save the velocity vector with in the particle data vector
        PIC::ParticleBuffer::SetV(v,modelParticleData);

        //sample the ebergy exchange rate
        //sample energy exchabge rate
        double EnergyLoss=0.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]-speed2)*PIC::MolecularData::GetMass(spec)*
            node->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(modelParticleData)/
            node->block->GetLocalTimeStep(spec);

        TotalModelParticleEnergyLossRate[spec][thread]+=EnergyLoss;
        *(spec+(double*)(TotalModelParticleEnergyLossRateOffset+PIC::Mesh::collectingCellSampleDataPointerOffset+cell->GetAssociatedDataBufferPointer()))+=EnergyLoss/cell->Measure;

        modelParticle=PIC::ParticleBuffer::GetNext(modelParticle);
      }
    }


#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif
  }

}

