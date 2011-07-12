//====================================================
//$Id$
//====================================================
//the general functions for the pic solver

#include "pic.h"


//sampling variables

long int PIC::LastSampleLength=0,PIC::CollectingSampleCounter=0,PIC::RequiredSampleLength=100,PIC::DataOutputFileNumber=0;
int PIC::SamplingMode=_RESTART_SAMPLING_MODE_;

//====================================================
//perform one time step
void PIC::TimeStep() {

  //sampling of the particle data
  PIC::Sampling();

  //move existing particles
  PIC::Mover::MoveParticles();

  //injection boundary conditions
  PIC::BC::InjectionBoundaryConditions();

  //syncrinie processors and exchnge particle data
  PIC::Parallel::ExchangeParticleData();


}
//====================================================
//the general sampling procedure
void PIC::Sampling() {
  int s,i,j,k,idim;
  long int LocalCellNumber,ptr;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  char *SamplingData;
  double *v,LocalParticleWeight;

  //go through the 'local nodes'
  while (node!=NULL) {
    block=node->block;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            cell=block->GetCenterNode(LocalCellNumber);
            SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;

            ptr=cell->FirstCellParticle;

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              s=PIC::ParticleBuffer::GetI(ParticleData);
              v=PIC::ParticleBuffer::GetV(ParticleData);

              LocalParticleWeight=block->GetLocalParticleWeight(s);

              //make a correction in the weight if individual particle weight is used
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_  == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
              //do nothing
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ ==  _INDIVIDUAL_PARTICLE_WEIGHT_ON_
              exit(__LINE__,__FILE__,"not implemented");
#else
              exit(__LINE__,__FILE__,"Error: unknown option");
#endif

              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberRelativeOffset))+=1;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;

              for (idim=0;idim<3;idim++) {
                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocityRelativeOffset))+=v[idim]*LocalParticleWeight;
                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocity2RelativeOffset))+=v[idim]*v[idim]*LocalParticleWeight;
              }

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }
          }
       }
    }

    node=node->nextNodeThisThread;
  }

  //Increment the sample length
  CollectingSampleCounter++;

  //Output the data flow file if needed
  if (CollectingSampleCounter==RequiredSampleLength) {
    //exchnge the sampling data
    PIC::Mesh::mesh.ParallelBlockDataExchange(_PIC_SUBDOMAIN_BOUNDARY_LAYER_SAMPLING_DATA_EXCHANGE_TAG_);

    //check different sampling modes
    if (SamplingMode==_RESTART_SAMPLING_MODE_) {
      PIC::Mesh::switchSamplingBuffers();
      LastSampleLength=CollectingSampleCounter;
      CollectingSampleCounter=0;

      for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        block=node->block;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::flushCollectingSamplingBuffer(block->GetCenterNode(LocalCellNumber));
        }
      }

      //flush sampling buffers in the internal surfaces installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
  #elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();
//      PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere* Sphere;

      cInternalSphericalData *Sphere;

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
//        Sphere=(PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere*)PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->GetSurfaceDataPointer();

        Sphere=PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface);

        PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
      }


  #else
      exit(__LINE__,__FILE__,"Error: unknown option");
  #endif


    }
    else if (SamplingMode==_ACCUMULATE_SAMPLING_MODE_) {
      exit(__LINE__,__FILE__,"Error: the sampling mode '_ACCUMULATE_SAMPLING_MODE_' is not implemented");
    }
    else exit(__LINE__,__FILE__,"Error: the sampling mode is not defined");

    //print output file
    char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld.dat",ChemSymbol,s,DataOutputFileNumber);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("done.\n");
        fflush(stdout);
      }

      //print into the file the sampled data of the internal surfaces installed into the mesh
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
#elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
        sprintf(fname,"pic.Sphere=%ld.%s.s=%i.out=%ld.dat",iSphericalSurface,ChemSymbol,s,DataOutputFileNumber);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("printing output file: %s.........",fname);
          fflush(stdout);
        }

        PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->PrintSurfaceData(fname,s);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("done.\n");
          fflush(stdout);
        }
      }


#else
      exit(__LINE__,__FILE__,"Error: unknown option");
#endif
    }



    //increment the output file number
    DataOutputFileNumber++;
  }





}

//====================================================
//init the particle solver
void PIC::Init() {

  //init MPI variables
  MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
  MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);
}



//====================================================
//set up the particle weight and time step
void PIC::Mesh::cDataBlockAMR::SetLocalParticleWeight(double weight, int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  *(spec+(double *)(associatedDataPointer+LocalParticleWeightOffset))=weight;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  if (PIC::ParticleWeightTimeStep::GlobalParticleWeight==NULL) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]>weight)) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]=weight;
  }

  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalParticleWeight(int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalParticleWeightOffset);
  return *res;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  return PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

//set up the particle time step
void PIC::Mesh::cDataBlockAMR::SetLocalTimeStep(double dt, int spec) {
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  *(spec+(double *)(associatedDataPointer+LocalTimeStepOffset))=dt;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_

  if (PIC::ParticleWeightTimeStep::GlobalTimeStep==NULL) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]>dt)) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]=dt;
  }


  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalTimeStep(int spec) {

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalTimeStepOffset);
  return *res;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  return  PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif


}




