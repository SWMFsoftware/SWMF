//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the model of photolytic reactions

#include "pic.h"

double *PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime=NULL;
//PIC::ChemicalReactions::PhotolyticReactions::fReactionProcessor *PIC::ChemicalReactions::PhotolyticReactions::ReactionProcessorTable=NULL;
//PIC::ChemicalReactions::PhotolyticReactions::fTotalLifeTime *PIC::ChemicalReactions::PhotolyticReactions::TotalLifeTime=NULL;


void PIC::ChemicalReactions::PhotolyticReactions::Init() {

  //only one particle transformation model can be used
  if (_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_) exit(__LINE__,__FILE__,"Error: only one particle transformation model can be used");

  ConstantTotalLifeTime=new double [PIC::nTotalSpecies];

  for (int i=0;i<PIC::nTotalSpecies;i++) {
    ConstantTotalLifeTime[i]=-1.0;
  }
}


double PIC::ChemicalReactions::PhotolyticReactions::TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  ReactionAllowedFlag=(ConstantTotalLifeTime[spec]>0.0) ? true : false;
  return ConstantTotalLifeTime[spec];
}

/*
void PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(fReactionProcessor f,int spec) {
  if (TotalLifeTime==NULL) Init();
  if ((spec<0)||(spec>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: out of range");

  ReactionProcessorTable[spec]=f;
}
void PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(fTotalLifeTime f,int spec) {
  if (TotalLifeTime==NULL) Init();
  if ((spec<0)||(spec>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: out of range");

  TotalLifeTime[spec]=f;
}
*/

int PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int code=_PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  double p,lifetime,c;
  bool flag;

  lifetime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(x,spec,ptr,flag,node);
  if (flag==false) return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;


  c=exp(-TimeInterval/lifetime);
  p=1.0-c; //the probability for reaction to occur

  if (rnd()<p) {
    //determine the time offset when the reaction had occured
    TimeInterval=-lifetime*log(1.0+rnd()*(c-1.0));

    return _PHOTOLYTIC_REACTION_OCCURES_;
  }

  return code;
}

void PIC::ChemicalReactions::PhotolyticReactions::InitProductStatWeight() {
  bool NewWeightFound=false;
  int s,sSource,sProduct,nProducts,i,iReaction;
  bool UnknownWeightTableBefore[PIC::nTotalSpecies],UnknownWeightTableAfter[PIC::nTotalSpecies];
  int nUnknownWeightsBefore;

  double StatWeight[PIC::nTotalSpecies];

  for (s=0,nUnknownWeightsBefore=0;s<PIC::nTotalSpecies;s++) {
    StatWeight[s]=PIC::ParticleWeightTimeStep::GlobalParticleWeight[s];
    UnknownWeightTableBefore[s]=false;

    if (StatWeight[s]==0) UnknownWeightTableBefore[s]=true,nUnknownWeightsBefore++;
  }

  memcpy(UnknownWeightTableAfter,UnknownWeightTableBefore,PIC::nTotalSpecies*sizeof(bool));

  if (nUnknownWeightsBefore!=0) {
    do {
      NewWeightFound=false;

      for (sSource=0;sSource<PIC::nTotalSpecies;sSource++) if (UnknownWeightTableBefore[sSource]==false) {
        //check if weight of products of this reaction is unknown
        for (iReaction=0;iReaction<nTotalUnimolecularReactions;iReaction++) if (UnimoleculecularReactionDescriptor[iReaction].SourceSpecie==sSource) {
          nProducts=UnimoleculecularReactionDescriptor[iReaction].nProducts;

          for (i=0;i<nProducts;i++) {
            sProduct=UnimoleculecularReactionDescriptor[iReaction].ProductSpecies[i];

            if (UnknownWeightTableBefore[sProduct]==true) {
              //the weight of the ptoduct species is unknown yet
              StatWeight[sProduct]+=UnimoleculecularReactionDescriptor[iReaction].ReactionRate*StatWeight[sSource];
              UnknownWeightTableAfter[sProduct]=false;
              NewWeightFound=true;
            }
          }

        }
      }


      memcpy(UnknownWeightTableBefore,UnknownWeightTableAfter,PIC::nTotalSpecies*sizeof(bool));
    }
    while (NewWeightFound==true);
  }

  //update the stat weight
  for (s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=StatWeight[s];
}

//=======================================================================================================
//default photolytic reaction processor
//by default particle will stay in the system
void PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReactionProcessor_default(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
  PIC::ParticleBuffer::SetPrev(-1,ptr);

  if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
  FirstParticleCell=ptr;
}

//=======================================================================================================
//execute the photolytic reaction model
void PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k;
  long int oldFirstCellParticle,newFirstCellParticle,p,pnext;



  //simulate particle's collisions
  /*
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    //reset the balancing counters
    for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      node=DomainBlockDecomposition::BlockTable[nLocalNode];
      if (node->block!=NULL) *(thread+(double*)(node->block->GetAssociatedDataBufferPointer()+LoadBalancingMeasureOffset))=0.0;
    }
#endif //_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

#if _PIC__OPENMP_THREAD_SPLIT_MODE_ == _PIC__OPENMP_THREAD_SPLIT_MODE__BLOCKS_
#pragma omp parallel for schedule(dynamic,_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_) default (none)  \
    private (k,j,i,node,oldFirstCellParticle,newFirstCellParticle,p,pnext) \
    shared (DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable)
#else
#pragma omp parallel for schedule(dynamic,1) default (none) \
private (k,j,i,node,oldFirstCellParticle,newFirstCellParticle,p,pnext) \
shared (DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable)
#endif  // _PIC__OPENMP_THREAD_SPLIT_MODE_

#endif  //_COMPILATION_MODE_
    for (int CellCounter=0;CellCounter<DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;CellCounter++) {
    int nLocalNode,ii=CellCounter;

    nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    k=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=k*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    j=ii/_BLOCK_CELLS_X_;
    ii-=j*_BLOCK_CELLS_X_;

    i=ii;*/

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none)  \
    private (k,j,i,node,oldFirstCellParticle,newFirstCellParticle,p,pnext) \
    shared (DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable)
#endif

for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {


    double StartTime=MPI_Wtime();

    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block!=NULL) {

/*      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {*/

      {{{
              oldFirstCellParticle=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
              newFirstCellParticle=-1;

              if (oldFirstCellParticle!=-1) {
                p=oldFirstCellParticle;

                while (p!=-1) {
                  pnext=PIC::ParticleBuffer::GetNext(p);
                  _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(p,newFirstCellParticle,node);
                  p=pnext;
                }

                node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=newFirstCellParticle;
              }
           }
         }
      }

    if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
      node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
    }

    }
  }
}
