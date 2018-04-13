
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


void SampleIndividualLocations(int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;

  //estimate the total flux and rigidity in a set of the defined locations
  if (Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0) {
    int LacalParticleNumber,GlobalParticleNumber;
    int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
    int nTotalInjectedParticles=0;

    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

    Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(std::max(1,Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations));

    do {
      //reset the partilce generation flag
      localParticleGenerationFlag=0;

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

              //set the particle generation flag
              localParticleGenerationFlag=1;

              //apply condition of tracking the particle
              if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
                PIC::ParticleTracker::InitParticleID(newParticleData);
                PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
              }

              *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex))=iLocation;
              *((double*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed))=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

              //set the initial value of the integrate path length
              if (Earth::CutoffRigidity::IntegratedPathLengthOffset!=-1) {
                *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=0.0;
              }

              //set up the particle rigidity
              if (Earth::CutoffRigidity::InitialRigidityOffset!=-1) {
                double momentum,charge,rigidity;

                charge=PIC::MolecularData::GetElectricCharge(spec);

                momentum=Relativistic::Speed2Momentum(speed,mass);
                rigidity=(charge>0.0) ? momentum/charge : 0.0;

                *((double*)(newParticleData+Earth::CutoffRigidity::InitialRigidityOffset))=rigidity;
              }

              //save the original location of the particle
              if (Earth::CutoffRigidity::InitialLocationOffset!=-1) {
                memcpy(newParticleData+Earth::CutoffRigidity::InitialLocationOffset,x,3*sizeof(double));
              }


            }
          }
        }
      }


      //limit the integration time (remove particle moving in trapped orbits)
      if ((true)&&(Earth::CutoffRigidity::InitialLocationOffset!=-1)) {
        //determine min altitude;
        static double rMax=-1.0;
        double r;
        int idim;

        if (rMax<0.0) {
          for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
            r=Vector3D::Length(Earth::CutoffRigidity::IndividualLocations::xTestLocationTable[iLocation]);

            if ((rMax<0.0)||(r>rMax)) rMax=r;
          }
        }

        //estimation of the domain length
        double IntegratioinLengthMax=sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
            pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
            pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

        //increase the limit of the total integrated oath length
        IntegratioinLengthMax/=4.0;

        IntegratioinLengthMax=2.0*Pi*rMax;

        //loop through all blocks
        for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
          int ptr,next;
          PIC::Mesh::cDataBlockAMR *block=node->block;
          double dt,l,v;
          PIC::ParticleBuffer::byte *ParticleData;

          if (sqrt(pow(node->xmax[0]+node->xmin[0],2)+pow(node->xmax[1]+node->xmin[1],2)+pow(node->xmax[2]+node->xmin[2],2))/2.0>rMax) continue;

          if (block!=NULL) for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
             ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

             while (ptr!=-1) {
               //loop through all particles and determin which are trapped and does not contribute to the incoming radiation
               ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

               v=Vector3D::Length(PIC::ParticleBuffer::GetV(ParticleData));
               next=PIC::ParticleBuffer::GetNext(ParticleData);

               dt=block->GetLocalTimeStep(PIC::ParticleBuffer::GetI(ParticleData));

               //update estimation of the total integration path length
               l=*((double*)(ParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset));
               l+=v*dt;
               *((double*)(ParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=l;


               if (l>IntegratioinLengthMax) {
                 //the particle is below the minimum altitude
                 PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
               }

               ptr=next;
             }

          }
        }
      }

      //preform the next iteration
      amps_time_step();

      static int LastDataOutputFileNumber=-1;

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }


      //get the total number of particles in the system
      LacalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LacalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

      //determine whether any particle has been generated during the current iteration
      MPI_Allreduce(&localParticleGenerationFlag,&globalParticleGenerationFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

      if (globalParticleGenerationFlag!=0) {
        //at least one particle has been generated -> reset the iteration counter
        IterationCounter=0;
      }
      else {
        //increment the iteration counter
        IterationCounter++;
      }


    }
    while ((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations));


    //determine the flux and eneregy spectra of the energetic particles in the poins of the observation
    if (true) {
      double ***EnergySpectrum,**TotalFlux,v[3],KineticEnergy,Speed,DiffFlux,dSurface,norm;
      int offset,iTestsLocation,spec,i,j,iface,iTable,jTable,Index,iBit,iByte,iE;

      const int nTotalEnergySpectrumIntervals=25;
      const double logMinEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
      const double logMaxEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);
      const double dE=(logMaxEnergyLimit-logMinEnergyLimit)/nTotalEnergySpectrumIntervals;

      //allocate the data buffers
      //TotalFlux[iLocation][spec]
      //EnergySpectrum[iLocation][spec][iEnergyInterval]

      TotalFlux=new double* [Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength];
      TotalFlux[0]=new double [Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength*PIC::nTotalSpecies];
      for (iTestsLocation=1;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) TotalFlux[iTestsLocation]=TotalFlux[iTestsLocation-1]+PIC::nTotalSpecies;

      EnergySpectrum=new double** [Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength];
      EnergySpectrum[0]=new double* [Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength*PIC::nTotalSpecies];
      for (iTestsLocation=1;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) EnergySpectrum[iTestsLocation]=EnergySpectrum[iTestsLocation-1]+PIC::nTotalSpecies;

      EnergySpectrum[0][0]=new double [Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength*PIC::nTotalSpecies*nTotalEnergySpectrumIntervals];

      for (offset=0,iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        EnergySpectrum[iTestsLocation][spec]=EnergySpectrum[0][0]+offset;
        offset+=nTotalEnergySpectrumIntervals;
      }

      //set the values of the buffers to zero
      for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        TotalFlux[iTestsLocation][spec]=0.0;

        for (i=0;i<nTotalEnergySpectrumIntervals;i++) EnergySpectrum[iTestsLocation][spec][i]=0.0;
      }

      //calculate the flux and energey spectrum
      for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) {
        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          for (iface=0;iface<6;iface++) {
            //surface area of the element of the surface mesh that covers the boundary of the computational domain
            dSurface=Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX[iface][0]*Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX[iface][1];

            //loop through the mesh that covers face 'iface' on the computational domain
            for (iTable=0;iTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;iTable++) {
              for (jTable=0;jTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;jTable++) {
                Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iTestsLocation][spec][iface][iTable][jTable].Gather();

                for (iByte=0;iByte<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iTestsLocation][spec][iface][iTable][jTable].FlagTableLength[0];iByte++) for (iBit=0;iBit<8;iBit++) {
                  Index=iBit+8*iByte;

                  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iTestsLocation][spec][iface][iTable][jTable].Test(Index)==true) {
                    //at least one particle that corrsponds to 'Index' has been detected. Add a contribution of such particles to the total energy spectrum and flux as observed at the point of the observation 'iTestsLocation'

                    Earth::CutoffRigidity::DomainBoundaryParticleProperty::ConvertVelocityVectorIndex2Velocity(spec,v,iface,Index);
                    Speed=Vector3D::Length(v);
                    KineticEnergy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

                    //probability density that particles has velocity 'v'
                    DiffFlux=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(KineticEnergy);

                    //determine the contributio of the particles into the 'observed' flux and energy spectrum
                    TotalFlux[iTestsLocation][spec]+=DiffFlux*dSurface;

                    //determine contribution of the particles to the energy flux
                    iE=(log10(KineticEnergy)-logMinEnergyLimit)/dE;
                    if (iE<0) iE=0;
                    if (iE>=nTotalEnergySpectrumIntervals) iE=nTotalEnergySpectrumIntervals-1;

                    EnergySpectrum[iTestsLocation][spec][iE]+=DiffFlux*dSurface;
                  }
                }
              }
            }
          }

          //normalize the energy spectrum
          for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) norm+=EnergySpectrum[iTestsLocation][spec][iE]*dE;
          if (norm>0) EnergySpectrum[iTestsLocation][spec][iE]/=norm;
        }
      }

      //output sampled particles flux and energy spectrum
      if (PIC::ThisThread==0) {
        //sampled energy spectrum
        FILE *fout;
        int spec,iTestsLocation,iE;
        char fname[400];

        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          sprintf(fname,"%s/EnergySpectrum[s=%i].dat",PIC::OutputDataFileDirectory,spec);
          fout=fopen(fname,"w");

          fprintf(fout,"VARIABLES=\"log10(Kinetic Energy[MeV]\"");
          for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) fprintf(fout,", \"Spectrum (iTestsLocation=%i)\"",iTestsLocation);
          fprintf(fout,"\n");

          for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) {
            double log10e=logMinEnergyLimit+iE*dE;
            double e=pow(10,log10e);

            e*=J2MeV;
            fprintf(fout,"%e  ",e);

            for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) fprintf(fout,"%e  ",EnergySpectrum[iTestsLocation][spec][iE]);
            fprintf(fout,"\n");
          }

          fclose(fout);
        }

        //The total energetic particle flux
        for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) {
          printf("spec=%i, iTestsLocation=%i: Flux=%e\n",spec,iTestsLocation,TotalFlux[iTestsLocation][spec]);
        }
      }



      //de-allocate the data buffers
      delete [] TotalFlux[0];
      delete [] TotalFlux;

      delete [] EnergySpectrum[0][0];
      delete [] EnergySpectrum[0];
      delete [] EnergySpectrum;

    }
    else {
      //non-uniform distribution of the injected particles
      exit(__LINE__,__FILE__,"Error: not implemented");
    }



  }

  //release sampling buffers
  Earth::CutoffRigidity::DomainBoundaryParticleProperty::Deallocate();
}



void SampleSphericalMaplLocations(double Radius,int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;
  int iLocation;
  double x[3]={0.0,0.0,0.0},v[3];


  int nZenithElements=50;
  int nAzimuthalElements=50;
  int nTotalInjectedParticlePerPoint=25;

  cInternalSphericalData Sphere;
  Sphere.SetGeneralSurfaceMeshParameters(nZenithElements,nAzimuthalElements);
  Sphere.SetSphereGeometricalParameters(x,Radius);

//  Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(nZenithElements*nAzimuthalElements);

  double ***EnergySpectrum,**TotalFlux,KineticEnergy,Speed,DiffFlux,dSurface,norm;
  int offset,iTestsLocation,spec,i,j,iface,iTable,jTable,Index,iBit,iByte,iE;

  const int nTotalEnergySpectrumIntervals=50;
  const double logMinEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
  const double logMaxEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);
  const double dE=(logMaxEnergyLimit-logMinEnergyLimit)/nTotalEnergySpectrumIntervals;

  TotalFlux=new double* [nZenithElements*nAzimuthalElements];
  TotalFlux[0]=new double [nZenithElements*nAzimuthalElements*PIC::nTotalSpecies];
  for (iTestsLocation=1;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) TotalFlux[iTestsLocation]=TotalFlux[iTestsLocation-1]+PIC::nTotalSpecies;

  EnergySpectrum=new double** [nZenithElements*nAzimuthalElements];
  EnergySpectrum[0]=new double* [nZenithElements*nAzimuthalElements*PIC::nTotalSpecies];
  for (iTestsLocation=1;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) EnergySpectrum[iTestsLocation]=EnergySpectrum[iTestsLocation-1]+PIC::nTotalSpecies;

  EnergySpectrum[0][0]=new double [nZenithElements*nAzimuthalElements*PIC::nTotalSpecies*nTotalEnergySpectrumIntervals];

  for (offset=0,iTestsLocation=0;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    EnergySpectrum[iTestsLocation][spec]=EnergySpectrum[0][0]+offset;
    offset+=nTotalEnergySpectrumIntervals;
  }

  //set the values of the buffers to zero
  for (iTestsLocation=0;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    TotalFlux[iTestsLocation][spec]=0.0;

    for (i=0;i<nTotalEnergySpectrumIntervals;i++) EnergySpectrum[iTestsLocation][spec][i]=0.0;
  }


  //estimate the total flux and rigidity in a set of the defined locations
  if (true) {
    int LacalParticleNumber,GlobalParticleNumber;
    int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
    int nTotalInjectedParticles=0;

    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;


    //loop through a part of the sphere
    int iZenithElementStart=0,iZenithElementFinish=0,iAzimutalElementStart=0,iAzimutalElementFinish=0;

    const int nSphereIndexIncrement=25;
    int iSphereIndex,iSphereIndexMin,iSphereIndexMax;

    iSphereIndexMin=0;
    iSphereIndexMax=std::min(nSphereIndexIncrement,nZenithElements*nAzimuthalElements);

    //allocate the buffers for output of the final resuts

    while ((iSphereIndexMin<nZenithElements*nAzimuthalElements)&&(iSphereIndexMax<nZenithElements*nAzimuthalElements)) {
      IterationCounter=0,nTotalInjectedParticles=0;

      //allocate the sampling buffer
      Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(iSphereIndexMax-iSphereIndexMin+1);


      do {
        //reset the partilce generation flag
        localParticleGenerationFlag=0;

        //Inject new particles
        int iZenithElement,iAzimutalElement;

        if (nTotalInjectedParticles<nTotalInjectedParticlePerPoint) {
          nTotalInjectedParticles+=nIngectedParticlePerIteration;

          //inject the new portion of the particles
          for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (iSphereIndex=iSphereIndexMin;iSphereIndex<=iSphereIndexMax;iSphereIndex++) {
            //          for (iZenithElement=0;iZenithElement<nZenithElements;iZenithElement++) for (iAzimutalElement=0;iAzimutalElement<nAzimuthalElements;iAzimutalElement++) {

            iZenithElement=iSphereIndex/nZenithElements;
            iAzimutalElement=iSphereIndex%nZenithElements;

            int idim,iCell,jCell,kCell;
            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode;

            iLocation=iSphereIndex-iSphereIndexMin; //Sphere.GetLocalSurfaceElementNumber(iZenithElement,iAzimutalElement);



//          if (startNode->Thread==PIC::ThisThread) {
            //generate a new particle velocity
            double mass,speed,energy,rigidity,momentum;
            int nd;

            static double logMinEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
            static double logMaxEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);

            mass=PIC::MolecularData::GetMass(spec);

            for (int iNewParticle=0;iNewParticle<nIngectedParticlePerIteration;iNewParticle++) {
              //location of the new particles
//               Vector3D::Distribution::Uniform(x,Radius);


              Sphere.GetSurfaceElementRandomPoint(x,iZenithElement,iAzimutalElement);
              startNode=PIC::Mesh::mesh.findTreeNode(x);

              if (startNode->Thread!=PIC::ThisThread) continue;

              if ((nd=PIC::Mesh::mesh.fingCellIndex(x,iCell,jCell,kCell,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

              #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
              PIC::Mesh::cDataCenterNode *cell;

              cell=startNode->block->GetCenterNode(nd);

              if (cell==NULL) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
              if (cell->Measure<=0.0) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
              #endif

              //generate energy of the new particle
              energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

              speed=Relativistic::E2Speed(energy,mass);
              Vector3D::Distribution::Uniform(v,speed);

              //generate a new particle
              long int newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);
              PIC::ParticleBuffer::byte *newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

              PIC::ParticleBuffer::SetV(v,newParticleData);
              PIC::ParticleBuffer::SetX(x,newParticleData);
              PIC::ParticleBuffer::SetI(spec,newParticleData);

              //set the particle generation flag
              localParticleGenerationFlag=1;

              //apply condition of tracking the particle
              if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
                PIC::ParticleTracker::InitParticleID(newParticleData);
                PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
              }

              *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex))=iLocation;
              *((double*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed))=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

              //set the initial value of the integrate path length
              if (Earth::CutoffRigidity::IntegratedPathLengthOffset!=-1) {
                *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=0.0;
              }

              //set up the particle rigidity
              if (Earth::CutoffRigidity::InitialRigidityOffset!=-1) {
                double momentum,charge,rigidity;

                charge=PIC::MolecularData::GetElectricCharge(spec);

                momentum=Relativistic::Speed2Momentum(speed,mass);
                rigidity=(charge>0.0) ? momentum/charge : 0.0;

                *((double*)(newParticleData+Earth::CutoffRigidity::InitialRigidityOffset))=rigidity;
              }

              //save the original location of the particle
              if (Earth::CutoffRigidity::InitialLocationOffset!=-1) {
                memcpy(newParticleData+Earth::CutoffRigidity::InitialLocationOffset,x,3*sizeof(double));
              }

            }
          }
        }

        //preform the next iteration
        amps_time_step();

        static int LastDataOutputFileNumber=-1;

        if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
          PIC::RequiredSampleLength*=2;
          if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


          LastDataOutputFileNumber=PIC::DataOutputFileNumber;
          if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
        }


        //get the total number of particles in the system
        LacalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
        MPI_Allreduce(&LacalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

        //determine whether any particle has been generated during the current iteration
        MPI_Allreduce(&localParticleGenerationFlag,&globalParticleGenerationFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

        if (globalParticleGenerationFlag!=0) {
          //at least one particle has been generated -> reset the iteration counter
          IterationCounter=0;
        }
        else {
          //increment the iteration counter
          IterationCounter++;
        }


      }
      while ((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations));


      //calculate the flux and energey spectrum
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none)  \
    private (iSphereIndex,spec,iface,dSurface,iTable,jTable,iByte,iBit,Index,Speed,KineticEnergy,DiffFlux,iE,norm,v) \
    shared (iSphereIndexMin,iSphereIndexMax,Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX,Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection, \
        Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable,TotalFlux,EnergySpectrum)
#endif
      for (iSphereIndex=iSphereIndexMin;iSphereIndex<=iSphereIndexMax;iSphereIndex++) {
        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          for (iface=0;iface<6;iface++) {
            //surface area of the element of the surface mesh that covers the boundary of the computational domain
            dSurface=Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX[iface][0]*Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX[iface][1];

            //loop through the mesh that covers face 'iface' on the computational domain
            for (iTable=0;iTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;iTable++) {
              for (jTable=0;jTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;jTable++) {
                Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iSphereIndex-iSphereIndexMin][spec][iface][iTable][jTable].Gather();

                for (iByte=0;iByte<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iSphereIndex-iSphereIndexMin][spec][iface][iTable][jTable].FlagTableLength[0];iByte++) for (iBit=0;iBit<8;iBit++) {
                  Index=iBit+8*iByte;

                  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[iSphereIndex-iSphereIndexMin][spec][iface][iTable][jTable].Test(Index)==true) {
                    //at least one particle that corrsponds to 'Index' has been detected. Add a contribution of such particles to the total energy spectrum and flux as observed at the point of the observation 'iTestsLocation'

                    Earth::CutoffRigidity::DomainBoundaryParticleProperty::ConvertVelocityVectorIndex2Velocity(spec,v,iface,Index);
                    Speed=Vector3D::Length(v);
                    KineticEnergy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

                    //probability density that particles has velocity 'v'
                    DiffFlux=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(KineticEnergy);

                    //determine the contributio of the particles into the 'observed' flux and energy spectrum
                    TotalFlux[iSphereIndex][spec]+=DiffFlux*dSurface;

                    //determine contribution of the particles to the energy flux
                    iE=(log10(KineticEnergy)-logMinEnergyLimit)/dE;
                    if (iE<0) iE=0;
                    if (iE>=nTotalEnergySpectrumIntervals) iE=nTotalEnergySpectrumIntervals-1;

                    EnergySpectrum[iSphereIndex][spec][iE]+=DiffFlux*dSurface;
                  }
                }
              }
            }
          }

          //normalize the energy spectrum
          for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) norm+=EnergySpectrum[iSphereIndex][spec][iE]*dE;
          if (norm>0) EnergySpectrum[iSphereIndex][spec][iE]/=norm;
        }
      }

      //Increment the sohere element index range
      iSphereIndexMin=iSphereIndexMax+1;
      iSphereIndexMax+=nSphereIndexIncrement;
      if (iSphereIndexMax>nZenithElements*nAzimuthalElements-1) iSphereIndexMax=nZenithElements*nAzimuthalElements-1;
    }



    //determine the flux and eneregy spectra of the energetic particles in the poins of the observation
    if (true) {
      //output sampled particles flux and energy spectrum
      if (PIC::ThisThread==0) {
        //sampled energy spectrum
        FILE *fout;
        int spec,iTestsLocation,iE;
        char fname[400];

        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          sprintf(fname,"%s/EnergySpectrum(R=%e)[s=%i].dat",PIC::OutputDataFileDirectory,Radius,spec);
          fout=fopen(fname,"w");

          fprintf(fout,"VARIABLES=\"log10(Kinetic Energy[MeV]\"");

          for (iTestsLocation=0;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) {
            double ZenithAngle,AzimuthalAngle;
            long int nZenithElement,nAzimuthalElement;

            Sphere.GetSurfaceElementMiddlePoint(x,iTestsLocation);
            Sphere.GetSurfaceElementProjectionIndex(x,ZenithAngle,nZenithElement,AzimuthalAngle,nAzimuthalElement);

            fprintf(fout,", \"Energy Spectrum (lon=%e, lat=%e)\"",ZenithAngle,AzimuthalAngle);
          }

          fprintf(fout,"\n");

          for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) {
            double log10e=logMinEnergyLimit+iE*dE;
            double e=pow(10,log10e);

            e*=J2MeV;
            fprintf(fout,"%e  ",e);

            for (iTestsLocation=0;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) fprintf(fout,"%e  ",EnergySpectrum[iTestsLocation][spec][iE]);
            fprintf(fout,"\n");
          }

          fclose(fout);
        }

        //The total energetic particle flux
        for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iTestsLocation=0;iTestsLocation<nZenithElements*nAzimuthalElements;iTestsLocation++) {
          printf("spec=%i, iTestsLocation=%i: Flux=%e\n",spec,iTestsLocation,TotalFlux[iTestsLocation][spec]);
        }
      }

      //output cutoff rigidity maps
      //IndividualLocations::CutoffRigidityTable
      CMPI_channel pipe(1000000);
      FILE *fout2d;

      if (PIC::ThisThread==0) {
        //sampled energy spectrum
        int spec,iTestsLocation,iE;
        char fname[400];

        pipe.openRecvAll();

        //open the output file
        sprintf(fname,"%s/CutoffRigidity[R=%e].dat",PIC::OutputDataFileDirectory,Radius);
        fout2d=fopen(fname,"w");

        fprintf(fout2d,"VARIABLES=\"Lon\", \"Lat\"");
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout2d,"\"Cutoff Rigidity (s=%i)\"",spec);

        fprintf(fout2d,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthalElements,nZenithElements+1);
      }
      else {
        pipe.openSend(0);
      }

      //interpolate and print the state vector
      long int InterpolationList[nZenithElements*nAzimuthalElements],InterpolationListLength=0;

      for (int iZenith=0;iZenith<nZenithElements+1;iZenith++) for (int iAzimuthal=0;iAzimuthal<nAzimuthalElements;iAzimuthal++) {
        Sphere.GetSurfaceCoordinate(x,iZenith,iAzimuthal);

        if (PIC::ThisThread==0) {
          double lon,lat;
          Sphere.GetSurfaceLonLatNormal(lon,lat,iZenith,iAzimuthal);
          fprintf(fout2d,"%e %e ",lon,lat);
        }


        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iZenith==0) {
          InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(0,iAzimuthal);
          InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalElements-1));
        }
        else if (iZenith==nZenithElements) {
          InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(nZenithElements-1,iAzimuthal);
          InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(nZenithElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalElements-1));
        }
        else {
          int iA,iZ,A[2],Z[2];

          Z[0]=iZenith-1,Z[1]=iZenith;

          A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalElements-1;
          A[1]=iAzimuthal;

          for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
        }

        //prepare and print the interpolated value of the cutoff rigidity
        //loop throught all species 
        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          //loop therough elements of the interpolation stencil
          double InterpolatedCutoffRigidity=0.0; 
          double norm=0.0;

          for (int el=0;el<InterpolationListLength;el++) {
            //loop through all MPI processes
            
            if (PIC::ThisThread==0) {
              //this process will output data
              double t,minElementRigidity=Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable[spec][InterpolationList[el]]; 

              //account for other MPI processes
              for (int thread=1;thread<PIC::nTotalThreads;thread++) {
                pipe.recv(t,thread);

                if ((t>0.0)&&(t<minElementRigidity)) minElementRigidity=t;  
              }

              //add the element cutoff rigidity to the interpiolated value
              norm+=1.0;
              InterpolatedCutoffRigidity+=minElementRigidity;
            }
            else {
              //send data to the "root" MPI process
              pipe.send(Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable[spec][InterpolationList[el]]); 
  
            }
          }

          //the interpolated value is calcualted by the root MPI process
          if (PIC::ThisThread==0) {
            InterpolatedCutoffRigidity/=norm;

            fprintf(fout2d,"  %e",InterpolatedCutoffRigidity);
          }
        }

        //the cutoff rigidity is computed for all speces at the given point in the map 
        if (PIC::ThisThread==0) fprintf(fout2d,"\n");
      }

      //close the pipe nad the file  
      if (ThisThread==0) {
        fclose(fout2d);
        pipe.closeRecvAll();
      }
      else pipe.closeSend();


      //de-allocate the data buffers
      delete [] TotalFlux[0];
      delete [] TotalFlux;

      delete [] EnergySpectrum[0][0];
      delete [] EnergySpectrum[0];
      delete [] EnergySpectrum;

    }

  }

  //release sampling buffers
  Earth::CutoffRigidity::DomainBoundaryParticleProperty::Deallocate();
}


int main(int argc,char **argv) {
  static int LastDataOutputFileNumber=0;


  Earth::CutoffRigidity::SampleRigidityMode=true;

  amps_init_mesh();

  Earth::CutoffRigidity::Init_BeforeParser();
  Earth::CutoffRigidity::AllocateCutoffRigidityTable();

  amps_init();




//  PIC::RequiredSampleLength=60;


#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  int nTotalIterations = 0;
#else
  int nTotalIterations = 100000001;
#endif


  int nMaxIterations=10; //000;

  //estimate the total flux and rigidity in a set of the defined locations

  SampleSphericalMaplLocations(_EARTH__RADIUS_+500.0E3,nMaxIterations);


  SampleIndividualLocations(nMaxIterations);



  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;

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
