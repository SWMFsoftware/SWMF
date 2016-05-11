//$Id$

/*
 * EnceladusMultiPlume_SourceModel.cpp
 *
 *  Created on: May 4, 2016
 *      Author: vtenishe
 */


#include "EnceladusMultiPlume.h"

//initialize the source model
void EnceladusMultiPlume::SourceModel::Init() {
  int i,spec,n;

  //reset the source rate tables
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    IndividualPlumes::TotalSourceRateTable[spec]=0.0;
    IndividualPlumes::maxTotalPlumeSourceRateTable[spec]=0.0;

    TigerStripes::TotalSourceRateTable[spec]=0.0;
    TigerStripes::maxTotalTigerStripeSourceRateTable[spec]=0.0;
  }

  //calculate the total source rate due to all individual plumes
  for (n=0;n<EnceladusMultiPlume::nTotalIndividualPlumes;n++) if (EnceladusMultiPlume::IndividualPlumeTable[n].ActiveFlag==true) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    IndividualPlumes::TotalSourceRateTable[spec]+=EnceladusMultiPlume::IndividualPlumeTable[n].SourceRate[spec];
    if (IndividualPlumes::maxTotalPlumeSourceRateTable[spec]<EnceladusMultiPlume::IndividualPlumeTable[n].SourceRate[spec]) IndividualPlumes::maxTotalPlumeSourceRateTable[spec]=EnceladusMultiPlume::IndividualPlumeTable[n].SourceRate[spec];
  }


  //Initialize the source due to the Tiger Stripes
  for (n=0;n<EnceladusMultiPlume::nTotalTigerStripes;n++) if (EnceladusMultiPlume::TigerStripeTable[n].ActiveFlag==true) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    TigerStripes::TotalSourceRateTable[spec]+=EnceladusMultiPlume::TigerStripeTable[n].SourceRate[spec];
    if (TigerStripes::maxTotalTigerStripeSourceRateTable[spec]<EnceladusMultiPlume::TigerStripeTable[n].SourceRate[spec]) TigerStripes::maxTotalTigerStripeSourceRateTable[spec]=EnceladusMultiPlume::TigerStripeTable[n].SourceRate[spec];
  }

  //calculate the number of the TigerStripe's segments, their length, and the maximum length of individual Tiger Stripe segment
  EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable=new int [EnceladusMultiPlume::nTotalTigerStripes];
  EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable=new double *[EnceladusMultiPlume::nTotalTigerStripes];
  EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable=new double [EnceladusMultiPlume::nTotalTigerStripes];

  for (n=0;n<EnceladusMultiPlume::nTotalTigerStripes;n++) {
    EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable[n]=0;
    EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable[n]=0.0;
    EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable[n]=NULL;

    if (EnceladusMultiPlume::TigerStripeTable[n].ActiveFlag==true) {
      EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable[n]=EnceladusMultiPlume::TigerStripeTable[n].Geomentry->nPoints-1;
      EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable[n]=new double [EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable[n]];

      int iSegment,idim;
      double length;

      for (iSegment=0;iSegment<EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsNumberTable[n];iSegment++) {
        length=0.0;

        for (idim=0;idim<3;idim++) {
          length+=pow(EnceladusMultiPlume::TigerStripeTable[n].Geomentry->cPoint3D[iSegment][idim]-EnceladusMultiPlume::TigerStripeTable[n].Geomentry->cPoint3D[iSegment+1][idim],2);
        }

        length=sqrt(length);

        EnceladusMultiPlume::SourceModel::TigerStripes::TigerStripeSegmentsLengthTable[n][iSegment]=length;
        if (EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable[n]<length)  EnceladusMultiPlume::SourceModel::TigerStripes::maxTigerStripeSegmentsLengthTable[n]=length;
      }
    }
  }

}


//return the total source rate
double EnceladusMultiPlume::SourceModel::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
  return IndividualPlumes::TotalSourceRateTable[spec]+TigerStripes::TotalSourceRateTable[spec];
}

//generate properties of a new model particle
bool EnceladusMultiPlume::SourceModel::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
        double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
  bool res;

  if (rnd()<IndividualPlumes::TotalSourceRateTable[spec]/(IndividualPlumes::TotalSourceRateTable[spec]+TigerStripes::TotalSourceRateTable[spec])) {
    res=IndividualPlumes::GenerateParticleProperties(spec,tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);
  }
  else {
    res=TigerStripes::GenerateParticleProperties(spec,tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);
  }

  return res;
}
