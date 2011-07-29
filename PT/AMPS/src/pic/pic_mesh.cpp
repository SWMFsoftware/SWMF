//====================================================
//$Id$
//====================================================
//the function for particle's data sampling

#include "pic.h"

//init the block's global data
int PIC::Mesh::cDataBlockAMR::LocalTimeStepOffset=0;
int PIC::Mesh::cDataBlockAMR::LocalParticleWeightOffset=0;
int PIC::Mesh::cDataBlockAMR::totalAssociatedDataLength=0;


//init the cells' sampling data
int PIC::Mesh::cDataCenterNode::totalAssociatedDataLength=0;


//the offsets to the sampled data stored in 'center nodes'
int PIC::Mesh::completedCellSampleDataPointerOffset=0,PIC::Mesh::collectingCellSampleDataPointerOffset=0;

int PIC::Mesh::sampledParticleWeghtRelativeOffset=0,PIC::Mesh::sampledParticleNumberRelativeOffset=0,PIC::Mesh::sampledParticleNumberDensityRelativeOffset=0;
int PIC::Mesh::sampledParticleVelocityRelativeOffset=0,PIC::Mesh::sampledParticleVelocity2RelativeOffset=0;
int PIC::Mesh::sampledExternalDataRelativeOffset=0;
int PIC::Mesh::sampleSetDataLength=0;

//the flag determines weather an external sampling procedure is used
bool PIC::Mesh::ExternalSamplingProcedureDefinedFlag=false;



void PIC::Mesh::SetCellSamplingDataRequest() {
  exit(__LINE__,__FILE__,"not implemented yet");
}


void PIC::Mesh::initCellSamplingDataBuffer() {

  if (cDataBlockAMR::totalAssociatedDataLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of the blocks associated data offsets");

  //local time step
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  cout << "Time step mode: specie dependent local time step" << endl;
  cDataBlockAMR::LocalTimeStepOffset=0;
  cDataBlockAMR::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  //do nothing for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //local particle weight
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  cout << "Particle weight mode: specie dependent local weight" << endl;
  cDataBlockAMR::LocalParticleWeightOffset=cDataBlockAMR::totalAssociatedDataLength;
  cDataBlockAMR::totalAssociatedDataLength+=sizeof(double)*PIC::nTotalSpecies;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  //do nothing for _SIMULATION_PARTICLE_WEIGHT_MODE_ ==_SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif

  //set up the offsets for 'center node' sampled data
  long int offset=0;

  PIC::Mesh::sampledParticleWeghtRelativeOffset=offset;
  offset+=sizeof(double)*PIC::nTotalSpecies;

  PIC::Mesh::sampledParticleNumberRelativeOffset=offset;
  offset+=sizeof(double)*PIC::nTotalSpecies;

  PIC::Mesh::sampledParticleNumberDensityRelativeOffset=offset;
  offset+=sizeof(double)*PIC::nTotalSpecies;

  PIC::Mesh::sampledParticleVelocityRelativeOffset=offset;
  offset+=3*sizeof(double)*PIC::nTotalSpecies;

  PIC::Mesh::sampledParticleVelocity2RelativeOffset=offset;
  offset+=3*sizeof(double)*PIC::nTotalSpecies;

  PIC::Mesh::sampledExternalDataRelativeOffset=offset;

  if (PIC::Mesh::ExternalSamplingProcedureDefinedFlag==true) {
    exit(__LINE__,__FILE__,"not implemented");
  }


  PIC::Mesh::sampleSetDataLength=offset;

  PIC::Mesh::completedCellSampleDataPointerOffset=0;
  PIC::Mesh::collectingCellSampleDataPointerOffset=PIC::Mesh::sampleSetDataLength;

  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength=2*PIC::Mesh::sampleSetDataLength;



}


//==============================================================
//get and set the local time step and particle weight
/*
double PIC::Mesh::GetLocalTimeStep(int spec,cDataBlockAMR* block) {
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  register double* ptr;
  ptr=(double*)(block->GetAssociatedDataBufferPointer()+cDataBlockAMR::LocalTimeStepOffset+spec*sizeof(double));
  return *ptr;
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

void PIC::Mesh::SetLocalTimeStep(double dt,int spec,cDataBlockAMR* block) {
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  register double* ptr;
  ptr=(double*)(block->GetAssociatedDataBufferPointer()+cDataBlockAMR::LocalTimeStepOffset+spec*sizeof(double));
  *ptr=dt;
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::GetLocalParticleWeight(int spec,cDataBlockAMR* block) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  register double* ptr;
  ptr=(double*)(block->GetAssociatedDataBufferPointer()+cDataBlockAMR::LocalParticleWeightOffset+spec*sizeof(double));
  return *ptr;
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

void PIC::Mesh::SetLocalParticleWeight(double weight,int spec,cDataBlockAMR* block) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  register double* ptr;
  ptr=(double*)(block->GetAssociatedDataBufferPointer()+cDataBlockAMR::LocalParticleWeightOffset+spec*sizeof(double));
  *ptr=weight;
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}
*/
//==============================================================
//flush and switch the sampling buffers in 'center' nodes
void PIC::Mesh::flushCompletedSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::flushCollectingSamplingBuffer(cDataCenterNode* node) {
  register int i,length=PIC::Mesh::sampleSetDataLength/sizeof(double);
  register double *ptr;

  for (i=0,ptr=(double*)(node->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset);i<length;i++,ptr++) *ptr=0.0;
}

void PIC::Mesh::switchSamplingBuffers() {
  int tempOffset;

  //exchange the CellSampleData offsets
  tempOffset=PIC::Mesh::completedCellSampleDataPointerOffset;
  PIC::Mesh::completedCellSampleDataPointerOffset=PIC::Mesh::collectingCellSampleDataPointerOffset;
  PIC::Mesh::collectingCellSampleDataPointerOffset=tempOffset;

  //switch the offsets for the internal spherical surfaces installed into the mesh
  PIC::BC::InternalBoundary::Sphere::switchSamplingBuffers();

}

//==============================================================
//init set and set up the computational mesh
void PIC::Mesh::Init(double* xMin,double* xMax,fLocalMeshResolution ResolutionFunction) {

  for (int idim=0;idim<DIM;idim++) xmin[idim]=xMin[idim],xmax[idim]=xmax[idim];

  LocalMeshResolution=ResolutionFunction;
  mesh.init(xMin,xMax,LocalMeshResolution);
}

void PIC::Mesh::buildMesh() {
  mesh.buildMesh();
}



void PIC::Mesh::cDataBlockAMR::sendBoundaryLayerBlockData(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->send(cell->associatedDataPointer,cell->totalAssociatedDataLength);
  }
}

void PIC::Mesh::cDataBlockAMR::sendMoveBlockAnotherProcessor(CMPI_channel *pipe) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  sendBoundaryLayerBlockData(pipe);

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  //send all blocks' data when the blocks is moved to another processor
  long int Particle,NextParticle;
  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);
    Particle=cell->FirstCellParticle;

    if  (Particle!=-1) {
      pipe->send(_CENTRAL_NODE_NUMBER_SIGNAL_);
      pipe->send(LocalCellNumber);

      while (Particle!=-1) {
        PIC::ParticleBuffer::PackParticleData(buffer,Particle);
        pipe->send(_NEW_PARTICLE_SIGNAL_);
         pipe->send(buffer,PIC::ParticleBuffer::ParticleDataLength);

        NextParticle=PIC::ParticleBuffer::GetNext(Particle);
        PIC::ParticleBuffer::DeleteParticle(Particle);
        Particle=NextParticle;
      }

      cell->FirstCellParticle=-1;
    }
  }

  pipe->send(_END_COMMUNICATION_SIGNAL_);
  delete [] buffer;
}

void PIC::Mesh::cDataBlockAMR::recvBoundaryLayerBlockData(CMPI_channel *pipe,int From) {
  int iCell,jCell,kCell;
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  #if DIM == 3
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
  #elif DIM == 2
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
  #elif DIM == 1
  static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
  #else
  exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
  #endif

  for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
    LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
    cell=GetCenterNode(LocalCellNumber);

    pipe->recv(cell->associatedDataPointer,cell->totalAssociatedDataLength,From);




//=========  DEBUG =========================

    double *p=(double*)0x11f8d4f40;

    p+=3;

    if (PIC::Mesh::mesh.ThisThread==3) if (*p==80) {
      cout << __FILE__ << __LINE__ << endl;
    }

//==========  END DEBUG =====================




  }
}

//recieve all blocks' data when the blocks is moved to another processo
void PIC::Mesh::cDataBlockAMR::recvMoveBlockAnotherProcessor(CMPI_channel *pipe,int From) {
  long int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell=NULL;

  recvBoundaryLayerBlockData(pipe,From);

  long int Particle;

  char *buffer=new char[PIC::ParticleBuffer::ParticleDataLength];

  int Signal;
  const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
  const int _NEW_PARTICLE_SIGNAL_=       2;
  const int _END_COMMUNICATION_SIGNAL_=  3;


  pipe->recv(Signal,From);
  LocalCellNumber=-1,cell=NULL;

  while (Signal!=_END_COMMUNICATION_SIGNAL_) {
    switch (Signal) {
    case _CENTRAL_NODE_NUMBER_SIGNAL_ :
      pipe->recv(LocalCellNumber,From);
      cell=GetCenterNode(LocalCellNumber);
      break;
    case _NEW_PARTICLE_SIGNAL_ :
      pipe->recv(buffer,PIC::ParticleBuffer::ParticleDataLength,From);

      Particle=PIC::ParticleBuffer::GetNewParticle(cell->FirstCellParticle);
      PIC::ParticleBuffer::UnPackParticleData(buffer,Particle);
      break;
    default :
      exit(__LINE__,__FILE__,"Error: unknown option");
    }

    pipe->recv(Signal,From);
  }

  delete [] buffer;
}

