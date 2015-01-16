//$Id$

/*
 * Exosphere_BackgroundPlasmaSource.cpp
 *
 *  Created on: Dec 17, 2014
 *      Author: vtenishe
 */

//The model of the ion source that inject ions throught the boundary of the computational domain accoding to the distribution of the background plasma flow

#include "pic.h"
#include "Exosphere.h"
#include "constants.h"

long int Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::nTotalBoundaryInjectionFaces=-1;  //the number of the computationsl mesh blocks at the boundary of the domain
double **Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::BoundaryFaceProductionFraction=NULL; //the fraction of the total production rate that is due to a particular block (BoundaryBlockProductionFraction[spec][block]
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::maxLocalTimeStep=NULL; //the maximum value of the time step across the boundary of the computational domain (maxLocalTimeStep[spec])
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::minParticleWeight=NULL; //the minimum value of the particle weight across the computational domain (minParticleWeight[spec])
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::TotalInjectionRateTable=NULL;

//allocate the buffers and init the model
void Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::getMinMaxLimits() {

  //calcualte the number of the boundary injection faces
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  bool ExternalFaces[6];

  for (nTotalBoundaryInjectionFaces=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
      for (int nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
        nTotalBoundaryInjectionFaces++;
      }
    }
  }

  //allocate the buffers
  if (maxLocalTimeStep==NULL) {
    maxLocalTimeStep=new double [PIC::nTotalSpecies];
    minParticleWeight=new double [PIC::nTotalSpecies];
  }

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    maxLocalTimeStep[s]=-1.0;
    minParticleWeight[s]=-1.0;
  }

  //determine the max-times step and min-particle weight in the blocks that containes the boundary of the domain
  double res=0.0,LocalParticleWeight,LocalTimeStep;
  long int nBoundaryFace,nd;
  int spec;

  double PlasmaTemeprature,PlasmaBulkVelocity[3],PlasmaNumberDensity,ExternalNormal[3],BlockSurfaceArea,x[3],x0[3],e0[3],e1[3];
  int nface,idim,i,j,k;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    maxLocalTimeStep[spec]=-1.0,minParticleWeight[spec]=-1.0;

    for (nBoundaryFace=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
      node=*nodeptr;

      if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
        for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
          if (node->Thread==PIC::Mesh::mesh.ThisThread) {
            LocalTimeStep=node->block->GetLocalTimeStep(spec);
            LocalParticleWeight=node->block->GetLocalParticleWeight(spec);

            if (maxLocalTimeStep[spec]<LocalTimeStep) maxLocalTimeStep[spec]=LocalTimeStep;
            if ((minParticleWeight[spec]<0.0)||(minParticleWeight[spec]>LocalParticleWeight)) minParticleWeight[spec]=LocalParticleWeight;
          }
        }
      }
    }


    double tmp[PIC::nTotalThreads];
    int thread;

    MPI_Allgather(maxLocalTimeStep+spec,1,MPI_DOUBLE,tmp,1,MPI_DOUBLE,MPI_GLOBAL_COMMUNICATOR);
    for (thread=0;thread<PIC::nTotalThreads;thread++) if (maxLocalTimeStep[spec]<tmp[thread]) maxLocalTimeStep[spec]=tmp[thread];

    MPI_Allgather(minParticleWeight+spec,1,MPI_DOUBLE,tmp,1,MPI_DOUBLE,MPI_GLOBAL_COMMUNICATOR);
    for (thread=0;thread<PIC::nTotalThreads;thread++) if ((minParticleWeight[spec]<0.0)||(minParticleWeight[spec]>tmp[thread])) minParticleWeight[spec]=tmp[thread];
  }

}

//calculate of the source rate
double Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::GetTotalProductionRate(int spec) {
  double res=0.0,LocalParticleWeight,LocalTimeStep;
  long int nBoundaryFace,nd;
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  double PlasmaTemeprature,PlasmaBulkVelocity[3],PlasmaNumberDensity,ExternalNormal[3],BlockSurfaceArea,x[3],x0[3],e0[3],e1[3];
  int nface,idim,i,j,k;
  bool ExternalFaces[6];

  if (IonNumberDensityFraction[spec]==0.0) return 0.0;

#if _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE_ == _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE__STEADY_STATE_
  if (TotalInjectionRateTable!=NULL) {
    if (TotalInjectionRateTable[spec]>=0.0) return TotalInjectionRateTable[spec];
  }
  else {
    TotalInjectionRateTable=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) TotalInjectionRateTable[s]=-1.0;
  }
#endif

  //determine the total number of the boundary faces and allocate 'BoundaryFaceProductionFraction'
  if (BoundaryFaceProductionFraction==NULL) {
    //calculate the number of the faces
    for (nTotalBoundaryInjectionFaces=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
      node=*nodeptr;

      if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
        for (int nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
          nTotalBoundaryInjectionFaces++;
        }
      }
    }

    BoundaryFaceProductionFraction=new double* [PIC::nTotalSpecies];
    BoundaryFaceProductionFraction[0]=new double [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      if (s!=0) BoundaryFaceProductionFraction[s]=BoundaryFaceProductionFraction[s-1]+nTotalBoundaryInjectionFaces;
      for (long int i=0;i<nTotalBoundaryInjectionFaces;i++) BoundaryFaceProductionFraction[s][i]=-1.0;
    }
  }

  //calculate the source rate
  for (nBoundaryFace=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
      for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
        BoundaryFaceProductionFraction[spec][nBoundaryFace]=0.0;

        if (node->Thread==PIC::Mesh::mesh.ThisThread) {
          node->GetExternalNormal(ExternalNormal,nface);
          BlockSurfaceArea=node->GetBlockFaceSurfaceArea(nface);

          //get the position and the cell number for which the background lasma conditions will be evaluated
          PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,node);
          for (idim=0;idim<3;idim++) x[idim]=x0[idim]+0.5*(e0[idim]+e1[idim])-PIC::Mesh::mesh.EPS*ExternalNormal[idim];
          nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);

          //determine the bachground plasma conditions at the block's face
          PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity(x,nd,node);
          PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity,x,nd,node);
          PlasmaTemeprature=PIC::CPLR::GetBackgroundPlasmaTemperature(x,nd,node);

          if ( (isfinite(PlasmaNumberDensity)==false) || (isfinite(PlasmaTemeprature)==false) || (isfinite(PlasmaBulkVelocity[0])==false) || (isfinite(PlasmaBulkVelocity[1])==false) || (isfinite(PlasmaBulkVelocity[2])==false) ) {
            exit(__LINE__,__FILE__,"Error: a non-normalized number is found");
          }

          BoundaryFaceProductionFraction[spec][nBoundaryFace]=IonNumberDensityFraction[spec]*BlockSurfaceArea*PIC::BC::CalculateInjectionRate_MaxwellianDistribution(PlasmaNumberDensity,PlasmaTemeprature,PlasmaBulkVelocity,ExternalNormal,spec);
        }

        nBoundaryFace++;
      }
    }
  }

  //exchange the calculated block's source rates among all processors
  double buffer[nTotalBoundaryInjectionFaces];

  MPI_Allreduce(BoundaryFaceProductionFraction[spec],buffer,nTotalBoundaryInjectionFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  for (res=0.0,nBoundaryFace=0;nBoundaryFace<nTotalBoundaryInjectionFaces;nBoundaryFace++) res+=buffer[nBoundaryFace];

  #if _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE_ == _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE__STEADY_STATE_
  TotalInjectionRateTable[spec]=res;
  #endif

  //calculate the production rate fraction
  for (nBoundaryFace=0;nBoundaryFace<nTotalBoundaryInjectionFaces;nBoundaryFace++) BoundaryFaceProductionFraction[spec][nBoundaryFace]=buffer[nBoundaryFace]/res;

  return res;
}

//particle injection
long int Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::ParticleInjection(int spec) {
  long int nInjectedParticles=0;
  double ModelParticlesInjectionRate,TimeCounter=0.0,TimeIncrement,Ratio,p;
  double ProductionFractionSum=0.0,FaceProductionFraction=0.0;

  if (IonNumberDensityFraction[spec]==0.0) return 0;

  //init the injection model if neded
  if (maxLocalTimeStep==NULL) getMinMaxLimits();

  ModelParticlesInjectionRate=GetTotalProductionRate(spec);
  Ratio=minParticleWeight[spec]/maxLocalTimeStep[spec];
  ModelParticlesInjectionRate/=minParticleWeight[spec];

  long int nBoundaryFace,nd;
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  double PlasmaTemeprature,PlasmaBulkVelocity[3],PlasmaNumberDensity,ExternalNormal[3],BlockSurfaceArea,x[3],x0[3],e0[3],e1[3];
  int nface,idim,i,j,k;
  bool ExternalFaces[6];


  //inject particle
  nBoundaryFace=0,nface=0;
  nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin();
  end=PIC::BC::boundingBoxInjectionBlocksList.end();

  //determine the fist boundary face and FaceProductionFraction
  bool facefound=false;

  for (nBoundaryFace=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
      for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
        FaceProductionFraction=BoundaryFaceProductionFraction[spec][nBoundaryFace];

        facefound=true;
        break;
      }
    }

    if (facefound==true) break;
  }

  //generate the time interval till the next particle injection
  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<maxLocalTimeStep[spec]) {
    double TimeStepFraction=TimeCounter/maxLocalTimeStep[spec];

    //determine the face for the particle injection
    if (ProductionFractionSum+FaceProductionFraction<TimeStepFraction) {
      ProductionFractionSum+=FaceProductionFraction;
      ++nface;
      goto FindInjectionFace; //enter into the search loop at the place it was exited the last time
    }

    while (ProductionFractionSum+FaceProductionFraction<TimeStepFraction) { //loop through the face untill the condition in met
      for (;nodeptr!=end;nodeptr++) {
        node=*nodeptr;
        nface=0;

        if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {

FindInjectionFace:
          for (;nface<2*DIM;nface++) {
            if (ExternalFaces[nface]==true) {
              if (nBoundaryFace==nTotalBoundaryInjectionFaces) exit(__LINE__,__FILE__,"Error: nBoundaryFace exxeds the allowed value range");

              FaceProductionFraction=BoundaryFaceProductionFraction[spec][++nBoundaryFace];

              if (ProductionFractionSum+FaceProductionFraction>TimeStepFraction) {
                //exit out of the loop - the block and the face are found
                goto StartParticleInjection;
              }

              ProductionFractionSum+=FaceProductionFraction;
            }
          }


        }
      }
    }

StartParticleInjection:
    //inject new particle
    double v[3],c0,c1;
    long int newParticle,nd;
    PIC::ParticleBuffer::byte *newParticleData;
    double PlasmaNumberDensity,PlasmaBulkVelocity[3],PlasmaTemeprature,LocalTimeStep,LocalParticleWeight;
    const double ParticleWeightCorrection=1.0;

    if (node->Thread==PIC::Mesh::mesh.ThisThread) {
      //calculate the probability of injection of the new particle at this moment
      LocalTimeStep=node->block->GetLocalTimeStep(spec);
      LocalParticleWeight=node->block->GetLocalParticleWeight(spec);

      p=LocalTimeStep/LocalParticleWeight*Ratio;

      if (rnd()<p) {
        //generate the new particle position on the face
        PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,node);
        node->GetExternalNormal(ExternalNormal,nface);

        for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim]-PIC::Mesh::mesh.EPS*ExternalNormal[idim];

         //generate a particle
         newParticle=PIC::ParticleBuffer::GetNewParticle();
         newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
         nInjectedParticles++;

         Exosphere::Sampling::CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__BACKGROUND_PLASMA_ION_INJECTION_]+=ParticleWeightCorrection*LocalParticleWeight/LocalTimeStep;
         PIC::BC::nInjectedParticles[spec]+=1;
         PIC::BC::ParticleProductionRate[spec]+=ParticleWeightCorrection*LocalParticleWeight/LocalTimeStep;

         //get macrospcopic parameters of the plasma at the point of the injection
         nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node); //findCenterNodeIndex(x,i,j,k,node); //fingCellIndex(x,i,j,k,node);

         //determine the bachground plasma conditions at the block's face
         PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity(x,nd,node);
         PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity,x,nd,node);
         PlasmaTemeprature=PIC::CPLR::GetBackgroundPlasmaTemperature(x,nd,node);


         do {
           PIC::Distribution::InjectMaxwellianDistribution(v,PlasmaBulkVelocity,PlasmaTemeprature,ExternalNormal,spec,-1);
         }
         while (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>vmax*vmax);

         PIC::ParticleBuffer::SetX(x,newParticleData);
         PIC::ParticleBuffer::SetV(v,newParticleData);
         PIC::ParticleBuffer::SetI(spec,newParticleData);
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);

         Exosphere::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__BACKGROUND_PLASMA_ION_INJECTION_,(PIC::ParticleBuffer::byte*)newParticleData);

         //apply condition of tracking the particle
         #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
         PIC::ParticleTracker::InitParticleID(newParticleData);
         PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData);
         #endif

          //inject the particle into the system
         _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,node->block->GetLocalTimeStep(spec)*rnd(),node,true);
      }
    }

  }


  return nInjectedParticles;
}

long int Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::ParticleInjection() {
  long int nInjectedParticles=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) nInjectedParticles+=ParticleInjection(spec);

  return nInjectedParticles;
}
