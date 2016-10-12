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
double **Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::BoundaryFaceTotalInjectionRate=NULL; //the fraction of the total production rate that is due to a particular block (BoundaryBlockProductionFraction[spec][block]
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::maxLocalTimeStep=NULL; //the maximum value of the time step across the boundary of the computational domain (maxLocalTimeStep[spec])
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::minParticleWeight=NULL; //the minimum value of the particle weight across the computational domain (minParticleWeight[spec])
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::TotalInjectionRateTable=NULL;
double **Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::maxBoundaryFaceLocalInjectionRate=NULL;
double *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::maxBoundaryFaceTotalInjectionRate; //the maximum value of the boundary face injectino rate
Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::cBoundaryFaceDescriptor *Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::BoundaryFaceDescriptor=NULL;
double ***Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::BoundaryFaceLocalInjectionRate=NULL;

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
  if (BoundaryFaceTotalInjectionRate==NULL) {
    //calculate the number of the faces
    for (nTotalBoundaryInjectionFaces=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
      node=*nodeptr;

      if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
        for (int nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
          nTotalBoundaryInjectionFaces++;
        }
      }
    }

    //allocate the buffers for the total face injection rate
    BoundaryFaceTotalInjectionRate=new double* [PIC::nTotalSpecies];
    BoundaryFaceTotalInjectionRate[0]=new double [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies];

    maxBoundaryFaceTotalInjectionRate=new double [PIC::nTotalSpecies];


    //allocate the buffers for the local injectino rate
    BoundaryFaceLocalInjectionRate=new double** [PIC::nTotalSpecies];
    BoundaryFaceLocalInjectionRate[0]=new double* [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies];
    BoundaryFaceLocalInjectionRate[0][0]=new double [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies*nFaceInjectionIntervals*nFaceInjectionIntervals];

    maxBoundaryFaceLocalInjectionRate=new double* [PIC::nTotalSpecies];
    maxBoundaryFaceLocalInjectionRate[0]=new double [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies];

    //allocate the injection face descriptors
    BoundaryFaceDescriptor=new cBoundaryFaceDescriptor[nTotalBoundaryInjectionFaces];


    //init the arrays
    int s,i,j,offsetBoundaryFaceLocalInjectionRate=0;

    for (s=0;s<PIC::nTotalSpecies;s++) {
      maxBoundaryFaceTotalInjectionRate[s]=0.0;

      if (s!=0) {
        BoundaryFaceTotalInjectionRate[s]=BoundaryFaceTotalInjectionRate[s-1]+nTotalBoundaryInjectionFaces;
        maxBoundaryFaceLocalInjectionRate[s]=maxBoundaryFaceLocalInjectionRate[s-1]+nTotalBoundaryInjectionFaces;
        BoundaryFaceLocalInjectionRate[s]=BoundaryFaceLocalInjectionRate[s-1]+nTotalBoundaryInjectionFaces;
      }

      for (i=0;i<nTotalBoundaryInjectionFaces;i++) {
        BoundaryFaceTotalInjectionRate[s][i]=0.0;
        maxBoundaryFaceLocalInjectionRate[s][i]=0.0;

        BoundaryFaceLocalInjectionRate[s][i]=BoundaryFaceLocalInjectionRate[0][0]+offsetBoundaryFaceLocalInjectionRate;
        offsetBoundaryFaceLocalInjectionRate+=nFaceInjectionIntervals*nFaceInjectionIntervals;

        for (j=0;j<nFaceInjectionIntervals*nFaceInjectionIntervals;j++) BoundaryFaceLocalInjectionRate[s][i][j]=0.0;
      }
    }
  }

  //calculate the source rate
  for (nBoundaryFace=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
    node=*nodeptr;

    if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
      for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
        BoundaryFaceTotalInjectionRate[spec][nBoundaryFace]=0.0;
        maxBoundaryFaceLocalInjectionRate[spec][nBoundaryFace]=0.0;

        if (node->Thread==PIC::Mesh::mesh.ThisThread) {
          node->GetExternalNormal(ExternalNormal,nface);
          BlockSurfaceArea=node->GetBlockFaceSurfaceArea(nface);

          //get the position and the cell number for which the background lasma conditions will be evaluated
          //get a uniform mesh on the face with the number of the elements that corresponds to that of the face
          int ii,jj;
          double InjectionRate,FluxIntegrationIncrement=1.0/nFaceInjectionIntervals;

          PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,node);

          for (ii=0;ii<nFaceInjectionIntervals;ii++) for (jj=0;jj<nFaceInjectionIntervals;jj++) {
            for (idim=0;idim<3;idim++) x[idim]=x0[idim]+FluxIntegrationIncrement*((ii+0.5)*e0[idim]+(jj+0.5)*e1[idim])-PIC::Mesh::mesh.EPS*ExternalNormal[idim];

            //determine the bachground plasma conditions at the block's face
            PIC::CPLR::InitInterpolationStencil(x,node);
            PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
            PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);
            PlasmaTemeprature=PIC::CPLR::GetBackgroundPlasmaTemperature();

            if ( (isfinite(PlasmaNumberDensity)==false) || (isfinite(PlasmaTemeprature)==false) || (isfinite(PlasmaBulkVelocity[0])==false) || (isfinite(PlasmaBulkVelocity[1])==false) || (isfinite(PlasmaBulkVelocity[2])==false) ) {
             exit(__LINE__,__FILE__,"Error: a non-normalized number is found");
            }

            if (PlasmaTemeprature>0.0) {
              InjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(PlasmaNumberDensity,PlasmaTemeprature,PlasmaBulkVelocity,ExternalNormal,spec);

              //sample the local face injection rate
              BoundaryFaceLocalInjectionRate[spec][nBoundaryFace][ii+jj*nFaceInjectionIntervals]=InjectionRate;
              if (maxBoundaryFaceLocalInjectionRate[spec][nBoundaryFace]<InjectionRate) maxBoundaryFaceLocalInjectionRate[spec][nBoundaryFace]=InjectionRate;

              //sampel the global face injection rate
              BoundaryFaceTotalInjectionRate[spec][nBoundaryFace]+=InjectionRate;
            }

          }

          BoundaryFaceTotalInjectionRate[spec][nBoundaryFace]*=IonNumberDensityFraction[spec]*BlockSurfaceArea*pow(FluxIntegrationIncrement,2);
        }

        BoundaryFaceDescriptor[nBoundaryFace].nface=nface;
        BoundaryFaceDescriptor[nBoundaryFace].node=node;

        nBoundaryFace++;
      }
    }
  }

  //exchange the calculated block's source rates among all processors
  double buffer[nTotalBoundaryInjectionFaces*nFaceInjectionIntervals*nFaceInjectionIntervals];

  MPI_Allreduce(BoundaryFaceTotalInjectionRate[spec],buffer,nTotalBoundaryInjectionFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  for (res=0.0,nBoundaryFace=0;nBoundaryFace<nTotalBoundaryInjectionFaces;nBoundaryFace++) {
    res+=(BoundaryFaceTotalInjectionRate[spec][nBoundaryFace]=buffer[nBoundaryFace]);
    if (maxBoundaryFaceTotalInjectionRate[spec]<buffer[nBoundaryFace]) maxBoundaryFaceTotalInjectionRate[spec]=buffer[nBoundaryFace];
  }

  //collect maxBoundaryFaceLocalProduction table
  MPI_Allreduce(maxBoundaryFaceLocalInjectionRate[spec],buffer,nTotalBoundaryInjectionFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  for (nBoundaryFace=0;nBoundaryFace<nTotalBoundaryInjectionFaces;nBoundaryFace++) maxBoundaryFaceLocalInjectionRate[spec][nBoundaryFace]=buffer[nBoundaryFace];

  //collect the local face injection rate
  MPI_Allreduce(BoundaryFaceLocalInjectionRate[spec][0],buffer,nTotalBoundaryInjectionFaces*nFaceInjectionIntervals*nFaceInjectionIntervals,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  for (i=0;i<nTotalBoundaryInjectionFaces*nFaceInjectionIntervals*nFaceInjectionIntervals;i++) BoundaryFaceLocalInjectionRate[spec][0][i]=buffer[i];

  #if _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE_ == _EXOSPHERE_SOURCE__BACKGROUND_PLASMA_ION_INJECTION_MODE__STEADY_STATE_
  TotalInjectionRateTable[spec]=res;
  #endif

  return res;
}

//particle injection
long int Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::ParticleInjection(int spec) {
  long int nInjectedParticles=0;
  double ModelParticlesInjectionRate,TimeCounter=0.0,TimeIncrement,Ratio,p;
  double ProductionFractionSum=0.0,FaceProductionFraction=0.0;
  double FluxIntegrationIncrement=1.0/nFaceInjectionIntervals;

  if (IonNumberDensityFraction[spec]==0.0) return 0;

  //init the injection model if neded
  if (maxLocalTimeStep==NULL) getMinMaxLimits();

  ModelParticlesInjectionRate=GetTotalProductionRate(spec);
  Ratio=minParticleWeight[spec]/maxLocalTimeStep[spec];
  ModelParticlesInjectionRate/=minParticleWeight[spec];

  long int nBoundaryFace,nd;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  double PlasmaTemeprature,PlasmaBulkVelocity[3],PlasmaNumberDensity,ExternalNormal[3],BlockSurfaceArea,x[3],x0[3],e0[3],e1[3];
  int nface,idim,i,j,k;
  bool ExternalFaces[6];

  double v[3],c0,c1;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  double LocalTimeStep,LocalParticleWeight;
  const double ParticleWeightCorrection=1.0;


  //generate the time interval till the next particle injection
  TimeCounter-=-log(rnd())/ModelParticlesInjectionRate*rnd(); //shift back the beginig of the time counting to account for the stopping of the time counting on the previous iteration

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel
   {
#pragma omp single
     {
#endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<maxLocalTimeStep[spec]) {
    double TimeStepFraction=TimeCounter/maxLocalTimeStep[spec];

    //determine the block and face for the particle injection
    do {
      nBoundaryFace=(int)(rnd()*nTotalBoundaryInjectionFaces);
    }
    while (rnd()>BoundaryFaceTotalInjectionRate[spec][nBoundaryFace]/maxBoundaryFaceTotalInjectionRate[spec]);

    node=BoundaryFaceDescriptor[nBoundaryFace].node;
    nface=BoundaryFaceDescriptor[nBoundaryFace].nface;

    //consider injection of a particle
    if (node->Thread==PIC::ThisThread) {
      //calculate the probability of injection of the new particle at this moment
      LocalTimeStep=node->block->GetLocalTimeStep(spec);
      LocalParticleWeight=node->block->GetLocalParticleWeight(spec);

      p=LocalTimeStep/LocalParticleWeight*Ratio;

      if (rnd()<p) {

        //sample the source rate
        nInjectedParticles++;
        Exosphere::Sampling::CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__BACKGROUND_PLASMA_ION_INJECTION_]+=ParticleWeightCorrection*LocalParticleWeight/LocalTimeStep;
        PIC::BC::nInjectedParticles[spec]+=1;
        PIC::BC::ParticleProductionRate[spec]+=ParticleWeightCorrection*LocalParticleWeight/LocalTimeStep;

        //initiate a task for OpenMP threads
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp task default (none) shared (PIC::Mesh::mesh,BoundaryFaceLocalInjectionRate,maxBoundaryFaceLocalInjectionRate,FluxIntegrationIncrement) \
         firstprivate (node,nface,spec,nBoundaryFace) \
         private (x,x0,e0,e1,ExternalNormal,newParticle,newParticleData,PlasmaNumberDensity,PlasmaTemeprature,PlasmaBulkVelocity,idim,c0,c1,v)
         {

           int thisThreadOpenMP=omp_get_thread_num();
      #else
           int thisThreadOpenMP=0;
      #endif

          //generate the new particle position on the face
          PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,node);
          node->GetExternalNormal(ExternalNormal,nface);

          //determine the location of the injection
          double *InjectionFaceLocalRateTable=BoundaryFaceLocalInjectionRate[spec][nBoundaryFace];
          int ii,jj,nSurfaceElement;

          do {
            nSurfaceElement=(int)(rnd()*nFaceInjectionIntervals*nFaceInjectionIntervals);
          }
          while (rnd()>InjectionFaceLocalRateTable[nSurfaceElement]/maxBoundaryFaceLocalInjectionRate[spec][nBoundaryFace]);

          jj=nSurfaceElement/nFaceInjectionIntervals;
          ii=nSurfaceElement%nFaceInjectionIntervals;

          for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+FluxIntegrationIncrement*((ii+c0)*e0[idim]+(jj+c1)*e1[idim])-PIC::Mesh::mesh.EPS*ExternalNormal[idim];

          //generate a particle
          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);




          //get macrospcopic parameters of the plasma at the point of the injection
          //determine the bachground plasma conditions at the block's face
          PIC::CPLR::InitInterpolationStencil(x,node);

          PlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
          PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);
          PlasmaTemeprature=PIC::CPLR::GetBackgroundPlasmaTemperature();

          do {
            PIC::Distribution::InjectMaxwellianDistribution(v,PlasmaBulkVelocity,PlasmaTemeprature,ExternalNormal,spec);
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
          PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)node);
          #endif

           //inject the particle into the system
          _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,node->block->GetLocalTimeStep(spec)*rnd(),node,true);

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  }
#endif //the task section

      }
    }

  }

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
     }}  //patallel and single sections
#endif

  return nInjectedParticles;
}

long int Exosphere::SourceProcesses::BackgroundPlasmaBoundaryIonInjection::ParticleInjection() {
  long int nInjectedParticles=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) nInjectedParticles+=ParticleInjection(spec);

  return nInjectedParticles;
}
