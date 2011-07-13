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

