
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>


#include <sys/time.h>
#include <sys/resource.h>

//$Id$



//the particle class
#include "pic.h"
#include "constants.h"
#include "Earth.h"


void amps_init();
void amps_init_mesh();
void amps_time_step();

int main(int argc,char **argv) {
  static int LastDataOutputFileNumber=0;




  amps_init_mesh();

  Earth::CutoffRigidity::Init_BeforeParser();
  Earth::CutoffRigidity::AllocateCutoffRigidityTable();

  amps_init();




  PIC::RequiredSampleLength=60;


#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  int nTotalIterations = 0;
#else
  int nTotalIterations = 100000001;
#endif

  //estimate the total flux and rigidity in a set of the defined locations
  if (Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0) {
    int LacalParticleNumber,GlobalParticleNumber;
    int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
    int nTotalInjectedParticles=0;

    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

    do {
      //Inject new particles
      if (nTotalInjectedParticles<Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations) {
        nTotalInjectedParticles+=nIngectedParticlePerIteration;

        //inject the new portion of the particles
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
          double x[3],v[3];
          int idim,iCell,jCell,kCell;
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode;

          for (idim=0;idim<3;idim++) x[idim]=Earth::CutoffRigidity::IndividualLocations::xTestLocationTable[iLocation][idim];
          startNode=PIC::Mesh::mesh.findTreeNode(x);

          if (startNode->Thread==PIC::ThisThread) {
            //generate a new particle velocity
            double mass,speed,energy,rigidity,momentum;

            static double logMinEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
            static double logMaxEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);

            mass=PIC::MolecularData::GetMass(spec);

            if (PIC::Mesh::mesh.fingCellIndex(x,iCell,jCell,kCell,startNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

            for (int iNewParticle=0;iNewParticle<nIngectedParticlePerIteration;iNewParticle++) {
              energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

              speed=Relativistic::E2Speed(energy,mass);
              Vector3D::Distribution::Uniform(v,speed);

              //generate a new particle
              long int newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);
              PIC::ParticleBuffer::byte *newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

              PIC::ParticleBuffer::SetV(v,newParticleData);
              PIC::ParticleBuffer::SetX(x,newParticleData);
              PIC::ParticleBuffer::SetI(spec,newParticleData);

              *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex))=iLocation;
              *((double*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed))=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            }
          }
        }
      }


      //preform the next iteration
      amps_time_step();

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }


      //get the total number of particles in the system
      LacalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LacalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    }
    while (GlobalParticleNumber!=0);
  }

  //time step with the backward integration integration
  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::ActiveFlag==true) {
    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

    //particles will be injected only in the near Earth's region
    Earth::BoundingBoxInjection::BoundaryInjectionMode=false;
    Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=true;

    for (long int niter=0;(niter<nTotalIterations)&&(LastDataOutputFileNumber<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::LastActiveOutputCycleNumber);niter++) {
      amps_time_step();

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }
    }

    Earth::CutoffRigidity::DomainBoundaryParticleProperty::Gather();
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::SmoothSampleTable();
    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_OFF_;

    //partices will be injected from the boundary of the domain
    Earth::BoundingBoxInjection::BoundaryInjectionMode=true;
    Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=false;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::ApplyInjectionPhaseSpaceLimiting=true;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=false;
  }

  //time step with the forward integration
  for (long int niter=0;niter<nTotalIterations;niter++) {
    amps_time_step();
    
    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);
      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",
	     ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,
	     PIC::RequiredSampleLength,
	     PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
  }
  
  
  //output the particle statistics of the test run
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  char fname[300];
  sprintf(fname,"%s/test_Earth.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
#endif
  
  
  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  
  return EXIT_SUCCESS;
}
