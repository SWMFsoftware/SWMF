
//$Id$

/*
 * BoundaryInjection.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: vtenishe
 */


#include "pic.h"


bool Orbiter::UpstreamBC::BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
   bool ExternalFaces[6];
   double ExternalNormal[3],ModelParticlesInjectionRate;
   int nface;

   if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
     for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
       startNode->GetExternalNormal(ExternalNormal,nface);
       ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Orbiter::UpstreamBC::NumberDensity,Orbiter::UpstreamBC::Temperature,Orbiter::UpstreamBC::Velocity,ExternalNormal,0);

       if (ModelParticlesInjectionRate>0.0) return true;
     }
   }

   return false;
 }

//injection of model particles through the faces of the bounding box
long int Orbiter::UpstreamBC::BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
   bool ExternalFaces[6];
   double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
   int nface,idim;
   long int newParticle;
   PIC::ParticleBuffer::byte *newParticleData;
   long int nInjectedParticles=0;
   double v[3];
   double ModelParticlesInjectionRate;

   if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
     ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
     LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


     for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
       startNode->GetExternalNormal(ExternalNormal,nface);
       TimeCounter=0.0;

       //injection boundary conditions: -x -> inflow, other -> open boundary
       double c=0.0;

       for (idim=0;idim<3;idim++) c+=Orbiter::UpstreamBC::Velocity[idim]*ExternalNormal[idim];

       if (-c<0.9) {
         nInjectedParticles+=PIC::BC::ExternalBoundary::OpenFlow::InjectBlock(spec,startNode,nface);
         continue;
       }

       ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Orbiter::UpstreamBC::NumberDensity,Orbiter::UpstreamBC::Temperature,Orbiter::UpstreamBC::Velocity,ExternalNormal,spec);

       if (ModelParticlesInjectionRate>0.0) {
         ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

         PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

         #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
         #pragma omp parallel default(none)  shared(nInjectedParticles,ModelParticlesInjectionRate,LocalTimeStep,x0,e0,e1,Orbiter::UpstreamBC::Velocity,Orbiter::UpstreamBC::Temperature,ExternalNormal,spec,startNode) \
           private (idim,c0,c1,x,v,newParticle,newParticleData) firstprivate (TimeCounter)
           {
             #pragma omp single
             {
         #endif //_COMPILATION_MODE_

         while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
           //increment the particle counter
           nInjectedParticles++;

           PIC::BC::nInjectedParticles[spec]++;
           PIC::BC::ParticleProductionRate[spec]+=ParticleWeight/LocalTimeStep;
           PIC::BC::ParticleMassProductionRate[spec]+=ParticleWeight*PIC::MolecularData::GetMass(spec)/LocalTimeStep;

           #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
           #pragma omp task default (none) shared(ModelParticlesInjectionRate,LocalTimeStep,x0,e0,e1,Orbiter::UpstreamBC::Velocity,Orbiter::UpstreamBC::Temperature,ExternalNormal,spec,startNode) \
             private (idim,c0,c1,x,v,newParticle,newParticleData) firstprivate (TimeCounter)
             {
           #endif //_COMPILATION_MODE_

           //generate the new particle position on the face
           for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

           //generate a particle
           newParticle=PIC::ParticleBuffer::GetNewParticle();
           newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);


           //generate particles' velocity
           PIC::Distribution::InjectMaxwellianDistribution(v,Orbiter::UpstreamBC::Velocity,Orbiter::UpstreamBC::Temperature,ExternalNormal,spec,-1);

           PIC::ParticleBuffer::SetX(x,newParticleData);
           PIC::ParticleBuffer::SetV(v,newParticleData);
           PIC::ParticleBuffer::SetI(spec,newParticleData);
           PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

           //inject the particle into the system
           _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);

           //end task
           #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
             }
           #endif
         }

         //end OpenMP parallel section
         #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
         }}
         #endif

       }

     }
   }

   return nInjectedParticles;
 }

long int Orbiter::UpstreamBC::BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;

  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}

double Orbiter::UpstreamBC::BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  bool ExternalFaces[6];
  double ExternalNormal[3],BlockSurfaceArea;
  int nface;

  double ModelParticlesInjectionRate=0.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);
      ModelParticlesInjectionRate+=BlockSurfaceArea*PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Orbiter::UpstreamBC::NumberDensity,Orbiter::UpstreamBC::Temperature,Orbiter::UpstreamBC::Velocity,ExternalNormal,spec);
    }
  }

  return ModelParticlesInjectionRate;
}


