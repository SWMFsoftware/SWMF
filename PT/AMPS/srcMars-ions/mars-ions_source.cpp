
//$Id$
//model of the ion source peocesses


/*
 * mars-ions_source.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: vtenishe
 */

#include "mars-ions.h"

double MarsIon::SourceProcesses::GetCellInjectionRate(int spec,double *xMiddle) {
  double res=0.0;
  double altitude=0.0;
  double rSeason=1.38758;
  double r2Season=pow(rSeason,2); 
  double Tnu_body = 134.0; //in K, neutral temperature
  double BodynDenNuSp_I_O= 8.0283e15; // per m^-3
  double BodynDenNuSp_I_Ox= 5.1736e14;
  double BodynDenNuSp_I_Oh= 6.3119e10;
  double BodynDenNuSp_I_Ohx= 3.9646e9;
  double HNuSpecies_I_O=13340; //scale height in m
  double HNuSpecies_I_Ox=50025;
  double HNuSpecies_I_Oh=290500;
  double HNuSpecies_I_Ohx=2436600;
  double Rate_I_O_Op= 6.346e-7/r2Season; //units s^-1 for O_hv__Op_em_
  double nDen_O;
  int i;    
  for (i=0;i<3;i++) {    
      altitude+=pow(xMiddle[i],2);
  }
  altitude=sqrt(altitude)-3396000.0; // in m where R_Mars=3396000.0 m

  if (altitude<0.0) return 0.0;

  nDen_O=BodynDenNuSp_I_O*exp(-altitude/HNuSpecies_I_O)+\
         BodynDenNuSp_I_Ox*exp(-altitude/HNuSpecies_I_Ox)+\
         BodynDenNuSp_I_Oh*exp(-altitude/HNuSpecies_I_Oh)+\
         BodynDenNuSp_I_Ohx*exp(-altitude/HNuSpecies_I_Ohx);
  //res=1.0;
  res=nDen_O*Rate_I_O_Op; //source rate per m^-3

  return res;
}

double MarsIon::SourceProcesses::GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell) {
  double res,*xMiddle;

  xMiddle=cell->GetX();
  res=GetCellInjectionRate(spec,xMiddle);

  return res*cell->Measure;
}


double MarsIon::SourceProcesses::GetBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double res=0.0;
  int i,j,k,LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;

  if (block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
    LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    cell=block->GetCenterNode(LocalCellNumber);

    if (cell!=NULL) res+=GetCellInjectionRate(spec,cell);
  }

  return res;
}

//inject model particles
long int MarsIon::SourceProcesses::InjectParticles() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;
  int i,j,k,LocalCellNumber,idim;
  char *data;
  double IonTemperature,*IonBulkVelocity,IonSourceRate;
  double *xCell,TimeCounter;
  int nInjectedParticles=0,newParticle;

  double ParticleWeight,LocalTimeStep,ParticleWeightCorrection;
  double v[3],x[3];

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_O_PLUS_SPEC_];
#else
  ParticleWeight=0.0;
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Exosphere::Planet->maxIntersectedNodeTimeStep[_O_PLUS_SPEC_];
#else
  LocalTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif



  while (node!=NULL) {
    double xMinBlock[3],xMaxBlock[3];

    block=node->block;
    memcpy(xMinBlock,node->xmin,3*sizeof(double));
    memcpy(xMaxBlock,node->xmax,3*sizeof(double));

    double dxCell[3]={(xMaxBlock[0]-xMinBlock[0])/_BLOCK_CELLS_X_,(xMaxBlock[1]-xMinBlock[1])/_BLOCK_CELLS_Y_,(xMaxBlock[2]-xMinBlock[2])/_BLOCK_CELLS_Z_};


    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
      LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      cell=block->GetCenterNode(LocalCellNumber);

      if (cell!=NULL) {
        data=cell->GetAssociatedDataBufferPointer();
        xCell=cell->GetX();

        IonSourceRate=(*(double*)(data+Output::OplusSource::RateOffset))*cell->Measure/ParticleWeight;
        IonBulkVelocity=(double*)(data+Output::OplusSource::BulkVelocityOffset);
        IonTemperature=*(double*)(data+Output::OplusSource::TemperatureOffset);

        //initial value of the time counter
        TimeCounter=rnd()*log(rnd())/IonSourceRate;

        //injection loop

        while ((TimeCounter+=-log(rnd())/IonSourceRate)<LocalTimeStep) {
          //get random position and velocity
          x[0]=xMinBlock[0]+dxCell[0]*(rnd()+i);
          x[1]=xMinBlock[1]+dxCell[1]*(rnd()+j);
          x[2]=xMinBlock[2]+dxCell[2]*(rnd()+k);

          //get random velocity vector
          PIC::Distribution::MaxwellianVelocityDistribution(v,IonBulkVelocity,IonTemperature,_O_PLUS_SPEC_);

          //generate new particle
          //the particle buffer used to set-up the new particle data
          PIC::ParticleBuffer::byte *newParticleData,tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
          PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

          //set the default value for the correction factor
          #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
          ParticleWeightCorrection=1.0;
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
          #endif

          #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
          if (block->GetLocalTimeStep(_O_PLUS_SPEC_)/LocalTimeStep<rnd()) continue;
          #endif

          //generate a particle
          PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
          PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
          PIC::ParticleBuffer::SetI(_O_PLUS_SPEC_,(PIC::ParticleBuffer::byte*)tempParticleData);

          #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
          #endif

          //apply condition of tracking the particle
          #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
          PIC::ParticleTracker::InitParticleID(tempParticleData);
          PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x_SO_OBJECT,v_SO_OBJECT,spec,tempParticleData);
          #endif

          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
          memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

          _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,block->GetLocalTimeStep(_O_PLUS_SPEC_)*rnd(),node,true);

          nInjectedParticles++;
        }



      }
    }

    node=node->nextNodeThisThread;
  }

  return nInjectedParticles;
}

