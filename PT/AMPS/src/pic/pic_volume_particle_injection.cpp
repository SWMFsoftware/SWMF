//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the functions describes process of particle injection from within the computational domain


#include "pic.h"

#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_

int PIC::VolumeParticleInjection::nRegistratedInjectionProcesses=0;
PIC::VolumeParticleInjection::cVolumeInjectionDescriptor PIC::VolumeParticleInjection::VolumeInjectionDescriptor[PIC::VolumeParticleInjection::nMaxInjectionProcessEntries];
double *PIC::VolumeParticleInjection::SourceRate=NULL;


void PIC::VolumeParticleInjection::GetRandomCellPosition(double *x,int iCell,int jCell,int kCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double xMinBlock[3],xMaxBlock[3];

  memcpy(xMinBlock,node->xmin,3*sizeof(double));
  memcpy(xMaxBlock,node->xmax,3*sizeof(double));

  #if DIM == 3
  const double dxCell[3]={(xMaxBlock[0]-xMinBlock[0])/_BLOCK_CELLS_X_,(xMaxBlock[1]-xMinBlock[1])/_BLOCK_CELLS_Y_,(xMaxBlock[2]-xMinBlock[2])/_BLOCK_CELLS_Z_};
  #elif DIM == 2
  const double dxCell[2]={(xMaxBlock[0]-xMinBlock[0])/_BLOCK_CELLS_X_,(xMaxBlock[1]-xMinBlock[1])/_BLOCK_CELLS_Y_};
  #elif DIM == 1
  const double dxCell[1]={(xMaxBlock[0]-xMinBlock[0])/_BLOCK_CELLS_X_};
  #else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
  #endif

  //generate the particle initial position
#if DIM == 3
  x[0]=xMinBlock[0]+dxCell[0]*(rnd()+iCell);
  x[1]=xMinBlock[1]+dxCell[1]*(rnd()+jCell);
  x[2]=xMinBlock[2]+dxCell[2]*(rnd()+kCell);
#elif DIM == 2
  x[0]=xMinBlock[0]+dxCell[0]*(rnd()+iCell);
  x[1]=xMinBlock[1]+dxCell[1]*(rnd()+jCell);
  x[2]=0.0;
#elif DIM == 1
  x[0]=xMinBlock[0]+dxCell[0]*(rnd()+iCell);
  x[1]=0.0,x[2]=0.0;
#else
exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif
}



void PIC::VolumeParticleInjection::InitTotalInjectionRate() {
  int iCell,jCell,kCell,nInjectionProcess,spec;
  bool InjectionFlag[PIC::nTotalSpecies];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  double res[PIC::nTotalSpecies],InjectionRate[PIC::nTotalSpecies];
  long int nd;

#if DIM == 3
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif


  if (nRegistratedInjectionProcesses==0) {
    printf("$PREFIX:WARNING: No Volume Injection Rate function is registreted (FILE=%s,LINE=%i)\n",__FILE__,__LINE__);
  }



  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
      cell=block->GetCenterNode(nd);

      for (spec=0;spec<PIC::nTotalSpecies;spec++) res[spec]=0.0;

      for (nInjectionProcess=0;nInjectionProcess<nRegistratedInjectionProcesses;nInjectionProcess++) {
        VolumeInjectionDescriptor[nInjectionProcess].SpeciesInjectionRate(InjectionFlag,InjectionRate,iCell,jCell,kCell,cell,node);

        for (spec=0;spec<PIC::nTotalSpecies;spec++) if (InjectionFlag[spec]==true) res[spec]+=InjectionRate[spec];
      }

      for (spec=0;spec<PIC::nTotalSpecies;spec++) *(spec+(double*)(PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset+cell->GetAssociatedDataBufferPointer()))=res[spec];
    }
  }
}



double PIC::VolumeParticleInjection::GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell) {
  if (nRegistratedInjectionProcesses==0) return 0.0;

  return *(spec+(double*)(PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset+cell->GetAssociatedDataBufferPointer()));
}

double PIC::VolumeParticleInjection::GetBlockInjectionRate(int spec,PIC::Mesh::cDataBlockAMR *block) {
  int iCell,jCell,kCell;
  PIC::Mesh::cDataCenterNode *cell;
  double res=0.0;
  long int nd;

#if DIM == 3
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif


  if (nRegistratedInjectionProcesses==0) return 0.0;

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=block->GetCenterNode(nd);

    res+=GetCellInjectionRate(spec,cell);
  }


  return res;
}

double PIC::VolumeParticleInjection::GetTotalInjectionRate(int spec) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  double res=0.0;

  if (nRegistratedInjectionProcesses==0) return 0.0;

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    res+=GetBlockInjectionRate(spec,node->block);
  }


  //collect the injection rate from all processors
  double buffer[PIC::nTotalThreads];

  MPI_Gather(&res,1,MPI_DOUBLE,buffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) for (int thread=1;thread<PIC::nTotalThreads;thread++) res+=buffer[thread];

  MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  return res;
}


double PIC::VolumeParticleInjection::GetTotalTimeStepInjection(int spec) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  double res=0.0,totalInjectionPerSecond=0.0,t;

  if (nRegistratedInjectionProcesses==0) return 0.0;

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    t=GetBlockInjectionRate(spec,node->block);

    res+=t*node->block->GetLocalTimeStep(spec);
    totalInjectionPerSecond+=t;
  }


  //collect the injection rate from all processors
  double buffer[PIC::nTotalThreads];

  MPI_Gather(&res,1,MPI_DOUBLE,buffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) for (int thread=1;thread<PIC::nTotalThreads;thread++) res+=buffer[thread];
  MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  MPI_Gather(&totalInjectionPerSecond,1,MPI_DOUBLE,buffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    for (int thread=1;thread<PIC::nTotalThreads;thread++) totalInjectionPerSecond+=buffer[thread];
    cout << "$PREFIX:The total prodiction rate for specie s=" << spec << " is " << totalInjectionPerSecond << " particles per second (" << __FILE__ << "@" << __LINE__ << ")" << endl;
  }

  return res;
}

//inject particles
long int PIC::VolumeParticleInjection::InjectParticle() {
  int iCell,jCell,kCell,nInjectionProcess;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nInjectedParticles=0;

  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];

  //the local copy of the volume injection processes
  int nRegistratedInjectionProcesses_local=nRegistratedInjectionProcesses;
  PIC::VolumeParticleInjection::cVolumeInjectionDescriptor localInjectionDescriptorList[nRegistratedInjectionProcesses_local];

  memcpy(localInjectionDescriptorList,VolumeInjectionDescriptor,nRegistratedInjectionProcesses*sizeof(PIC::VolumeParticleInjection::cVolumeInjectionDescriptor));



  //check if the namespace is initialized
  if (nRegistratedInjectionProcesses==0) exit(__LINE__,__FILE__,"Error: the model is not initialized");

#if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
#elif DIM == 2
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
#elif DIM == 1
  const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
#else
  exit(__LINE__,__FILE__,"Error: the option is not defined");
#endif

  //sample the processor load
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double EndTime,StartTime=MPI_Wtime();
#endif

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
//      cell=node->block->GetCenterNode(nd);
      cell=cellList[nd];

//      for (nInjectionProcess=0;nInjectionProcess<nRegistratedInjectionProcesses;nInjectionProcess++) nInjectedParticles+=VolumeInjectionDescriptor[nInjectionProcess].InjectionProcessor(iCell,jCell,kCell,cell,node);

      for (nInjectionProcess=0;nInjectionProcess<nRegistratedInjectionProcesses_local;nInjectionProcess++) nInjectedParticles+=localInjectionDescriptorList[nInjectionProcess].InjectionProcessor(iCell,jCell,kCell,cell,node);
    }

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    EndTime=MPI_Wtime();
    node->ParallelLoadMeasure+=EndTime-StartTime;
    StartTime=EndTime;
#endif
  }

  return nInjectedParticles;
}



void PIC::VolumeParticleInjection::Init() {

  //set up the velues of the local particle injection rate
  InitTotalInjectionRate();
}

#endif
