//$Id$



/*
 * EnceladusMultiPlume_SourceTigerStripe.cpp
 *
 *  Created on: May 4, 2016
 *      Author: vtenishe
 */


#include "EnceladusMultiPlume.h"


double EnceladusMultiPlume::SourceModel::TigerStripes::TotalSourceRateTable[PIC::nTotalSpecies];
//double EnceladusMultiPlume::SourceModel::TigerStripes::TotalTigerStripeSourceRateTable[nTotalTigerStripes][PIC::nTotalSpecies];
double EnceladusMultiPlume::SourceModel::TigerStripes::maxTotalTigerStripeSourceRateTable[PIC::nTotalSpecies];

//the number of the Tiger Stripe's segments and their length
int *EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable;
double **EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable;
double *EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable;


//generate properties of a new model particle
bool EnceladusMultiPlume::SourceModel::TigerStripes::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
    double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {


  //1. Determine the Tiger Stripe through which the new particle is injected
  int iTigerStripe;

  do {
    do {
      iTigerStripe=(int)(rnd()*EnceladusMultiPlume::nTotalTigerStripes);
    }
    while (EnceladusMultiPlume::TigerStripeTable[iTigerStripe].ActiveFlag==false);
  }
  while (rnd()>EnceladusMultiPlume::TigerStripeTable[iTigerStripe].SourceRate[spec]/maxTotalTigerStripeSourceRateTable[spec]);


  //2. Determine segment for the new particle injection
  int iSegment;

  do {
    iSegment=(int)(rnd()*EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable[iTigerStripe]);
  }
  while(rnd()>EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable[iTigerStripe][iSegment]/EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable[iTigerStripe]);

  //3. Determine the coordinate of the injection
  double xLocal=rnd();

  for (int idim=0;idim<3;idim++) x_IAU_OBJECT[idim]=EnceladusMultiPlume::TigerStripeTable[iTigerStripe].Geomentry->cPoint3D[iSegment][idim]+
      xLocal*(EnceladusMultiPlume::TigerStripeTable[iTigerStripe].Geomentry->cPoint3D[iSegment+1][idim]-EnceladusMultiPlume::TigerStripeTable[iTigerStripe].Geomentry->cPoint3D[iSegment][idim]);

  //if the node location does not belongs to the current processor -> exit
  startNode=PIC::Mesh::mesh.findTreeNode(x_IAU_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

  //determine the new particle parameters
  //determine velocity of the new particle
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0},ExternalNorm[3];
  int idim;

  memcpy(ExternalNorm,x_IAU_OBJECT,3*sizeof(double));
  Vector3D::Normalize(ExternalNorm);

  for (int idim=0;idim<3;idim++) {
    vbulk[idim]=EnceladusMultiPlume::TigerStripeTable[iTigerStripe].BulkSpeed[spec]*ExternalNorm[idim];
    ExternalNorm[idim]*=-1.0;
  }

  SurfaceTemperature=EnceladusMultiPlume::TigerStripeTable[iTigerStripe].Temperature[spec];
  PIC::Distribution::InjectMaxwellianDistribution(v_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNorm,spec);


  //copy the location and velocity vectors into the provided buffers
  memcpy(v_SO_OBJECT,v_IAU_OBJECT,3*sizeof(double));
  memcpy(x_SO_OBJECT,x_IAU_OBJECT,3*sizeof(double));

  //return the sucess code
  return true;

}

