//The physical model of GCR at the boundary of the computational domain

/*
 * BoundaryInjection_GCR.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include "Earth.h"

double **Earth::BoundingBoxInjection::GCR::BoundaryFaceLocalInjectionRate=NULL;
double *Earth::BoundingBoxInjection::GCR::maxBoundaryFaceLocalInjectionRate=NULL;
double *Earth::BoundingBoxInjection::GCR::BoundaryFaceTotalInjectionRate=NULL;
cBoundaryFaceDescriptor *Earth::BoundingBoxInjection::BoundaryFaceDescriptor=NULL;
int Earth::BoundingBoxInjection::nTotalBoundaryInjectionFaces=0;

//get the total injection rate
double Earth::BoundingBoxInjection::GCR::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

//  exit(__LINE__,__FILE__,"Error: not implemented");

  return 1.0;
}

//init the boundary box injection tables
void Earth::BoundingBoxInjection::GCR::InitBoundingBoxInjectionTable() {
  int spec,iBoundaryFace,nface;
  bool ExternalFaces[6];
  double ExternalNormal[3],FaceSurfaceArea,FaceSourceRate;

  //allcate the tables if needed
  if (maxBoundaryFaceLocalInjectionRate==NULL) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
    list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
    bool ExternalFaces[6];

    //calculate the number of the faces
    for (nTotalBoundaryInjectionFaces=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
      node=*nodeptr;

      if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
       for (int nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
         nTotalBoundaryInjectionFaces++;
       }
      }
    }

    //allocate the tables
    //allocate the buffers for the total face injection rate
    BoundaryFaceTotalInjectionRate=new double [PIC::nTotalSpecies];
    maxBoundaryFaceLocalInjectionRate=new double [PIC::nTotalSpecies];

    //allocate the buffers for the local injectino rate
    BoundaryFaceLocalInjectionRate=new double* [PIC::nTotalSpecies];
    BoundaryFaceLocalInjectionRate[0]=new double [nTotalBoundaryInjectionFaces*PIC::nTotalSpecies];

    for (spec=0;spec<PIC::nTotalSpecies;spec++) BoundaryFaceLocalInjectionRate[spec]=BoundaryFaceLocalInjectionRate[0]+spec*nTotalBoundaryInjectionFaces;


    //allocate the injection face descriptors
    if (BoundaryFaceDescriptor==NULL) {
      BoundaryFaceDescriptor=new cBoundaryFaceDescriptor[nTotalBoundaryInjectionFaces];
    }
  }

  //init the tables with zeros
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    BoundaryFaceTotalInjectionRate[spec]=0.0;
    maxBoundaryFaceLocalInjectionRate[spec]=0.0;

    for (iBoundaryFace=0;iBoundaryFace<nTotalBoundaryInjectionFaces;iBoundaryFace++) {
      BoundaryFaceLocalInjectionRate[spec][iBoundaryFace]=0.0;
    }
  }

  //populate the source with the source rate values and calculate the total source rate
  list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* >::iterator end,nodeptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

   for (iBoundaryFace=0,nodeptr=PIC::BC::boundingBoxInjectionBlocksList.begin(),end=PIC::BC::boundingBoxInjectionBlocksList.end();nodeptr!=end;nodeptr++) {
     node=*nodeptr;

     if (PIC::Mesh::mesh.ExternalBoundaryBlock(node,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
       for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
         for (spec=0;spec<PIC::nTotalSpecies;spec++) {
           if (node->Thread==PIC::Mesh::mesh.ThisThread) {
             node->GetExternalNormal(ExternalNormal,nface);
             FaceSurfaceArea=node->GetBlockFaceSurfaceArea(nface);

             FaceSourceRate=InjectionRate(spec,nface,node)*FaceSurfaceArea;

             if (maxBoundaryFaceLocalInjectionRate[spec]<FaceSourceRate) maxBoundaryFaceLocalInjectionRate[spec]=FaceSourceRate;
             BoundaryFaceTotalInjectionRate[spec]+=FaceSourceRate;
             BoundaryFaceLocalInjectionRate[spec][iBoundaryFace]=FaceSourceRate;
           }
         }

         BoundaryFaceDescriptor[iBoundaryFace].nface=nface;
         BoundaryFaceDescriptor[iBoundaryFace].node=node;

         iBoundaryFace++;
       }
     }
   }

   //combine the maps from all processors
   double ExchangeBuffer[PIC::nTotalSpecies*nTotalBoundaryInjectionFaces];

   //the total production rate
   MPI_Allreduce(BoundaryFaceTotalInjectionRate,ExchangeBuffer,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
   memcpy(BoundaryFaceTotalInjectionRate,ExchangeBuffer,PIC::nTotalSpecies*sizeof(double));

   //source rate from each face
   MPI_Allreduce(BoundaryFaceLocalInjectionRate[0],ExchangeBuffer,PIC::nTotalSpecies*nTotalBoundaryInjectionFaces,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
   memcpy(BoundaryFaceLocalInjectionRate[0],ExchangeBuffer,PIC::nTotalSpecies*nTotalBoundaryInjectionFaces*sizeof(double));

   //maximum injection rate from the boundary faces
   if (PIC::ThisThread==0) {
     MPI_Status status;
     int thread;

     for (thread=1;thread<PIC::nTotalThreads;thread++) {
       MPI_Recv(ExchangeBuffer,PIC::nTotalSpecies,MPI_DOUBLE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

       for (spec=0;spec<PIC::nTotalSpecies;spec++) {
         if (maxBoundaryFaceLocalInjectionRate[spec]<ExchangeBuffer[spec]) maxBoundaryFaceLocalInjectionRate[spec]=ExchangeBuffer[spec];
       }
     }
   }
   else {
     MPI_Send(maxBoundaryFaceLocalInjectionRate,PIC::nTotalSpecies,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);
   }

   MPI_Bcast(maxBoundaryFaceLocalInjectionRate,PIC::nTotalSpecies,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
}

//generate properties of a new particle
void Earth::BoundingBoxInjection::GCR::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {
  double v[3],Energy,mass,WeightCorrectionFactor;

  //power law energy distribution will be used for testing
  //later more realisting models of GCR (Nymmic, etc) will be added to the library
  //of GCR boundary box injection models

  const double EnergyDistributionPowerIndex=3;

  //Distribute energy
  Energy=minEnergy*(1.0+rnd()*(maxEnergy-minEnergy));
  WeightCorrectionFactor=pow(Energy,-EnergyDistributionPowerIndex)/pow(minEnergy,-EnergyDistributionPowerIndex);


  //generate velocity of the new particle
  PIC::Distribution::InjectRingDistribution(v,Energy,ExternalNormal,spec);

  //save parameters of the new particle
  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,ParticleData);

}

//init the SEP injection model
void Earth::BoundingBoxInjection::GCR::Init() {
  exit(__LINE__,__FILE__,"Error: not implemented");
}
