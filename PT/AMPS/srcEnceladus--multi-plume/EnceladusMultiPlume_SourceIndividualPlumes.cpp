
//$Id$


/*
 * EnceladusMultiPlume_SourceIndividualPlumes.cpp
 *
 *  Created on: May 4, 2016
 *      Author: vtenishe
 */

#include "EnceladusMultiPlume.h"

double EnceladusMultiPlume::SourceModel::IndividualPlumes::TotalSourceRateTable[PIC::nTotalSpecies];
double EnceladusMultiPlume::SourceModel::IndividualPlumes::maxTotalPlumeSourceRateTable[PIC::nTotalSpecies];


//generate properties of a new model particle
bool EnceladusMultiPlume::SourceModel::IndividualPlumes::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
    double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {


  //determine the plume for the next particle injection
  int nplume;

  do {
    do {
      nplume=rnd()*EnceladusMultiPlume::nTotalIndividualPlumes;
    }
    while (EnceladusMultiPlume::IndividualPlumeTable[nplume].ActiveFlag==false);
  }
  while (rnd()>EnceladusMultiPlume::IndividualPlumeTable[nplume].SourceRate[spec]/maxTotalPlumeSourceRateTable[spec]);


  //if the node location does not belongs to the current processor -> exit
  startNode=PIC::Mesh::mesh.findTreeNode(EnceladusMultiPlume::IndividualPlumeTable[nplume].xLoacation,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

  //determine the new particle parameters
  //determine velocity of the new particle
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0},ExternalNorm[3];
  int idim;

  memcpy(ExternalNorm,EnceladusMultiPlume::IndividualPlumeTable[nplume].lPointing,3*sizeof(double));
  Vector3D::Normalize(ExternalNorm);

  for (int idim=0;idim<3;idim++) {
    vbulk[idim]=EnceladusMultiPlume::IndividualPlumeTable[nplume].BulkSpeed[spec]*ExternalNorm[idim];
    ExternalNorm[idim]*=-1.0;
  }

  SurfaceTemperature=EnceladusMultiPlume::IndividualPlumeTable[nplume].Temperature[spec];
  PIC::Distribution::InjectMaxwellianDistribution(v_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNorm,spec);


  //copy the location and velocity vectors into the provided buffers
  memcpy(v_SO_OBJECT,v_IAU_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,EnceladusMultiPlume::IndividualPlumeTable[nplume].xLoacation,3*sizeof(double));
  memcpy(x_SO_OBJECT,x_IAU_OBJECT,3*sizeof(double));

  //return the sucess code
  return true;
}
