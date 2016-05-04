/*
 * Enceladus.h
 *
 *  Created on: Jan 17, 2012
 *      Author: vtenishe
 */

//$Id$

#ifndef ENCELADUS_H_
#define ENCELADUS_H_

#include "pic.h"
#include "Exosphere.h"
#include "Dust.h"

#include "EnceladusMultiPlume.dfn"


namespace EnceladusMultiPlume {
  using namespace ElectricallyChargedDust;

  //definitions of the Tifer Stripes and Individual plume classes
  class cPoint3D {
  public:
    double x[3];
  };

  class cTigerStripeGeometry {
  public:
    int nPoints;
    double cPoint3D[25][3];
  };

  class cTigerStripe {
  public:
    cTigerStripeGeometry *Geomentry;
    double SourceRate[PIC::nTotalSpecies];
    double Temperature[PIC::nTotalSpecies];
    double BulkSpeed[PIC::nTotalSpecies];
    char ID[200];
    bool ActiveFlag;
  };

  class cIndividualPlume {
  public:
    char ID[200];
    double SourceRate[PIC::nTotalSpecies];
    double Temperature[PIC::nTotalSpecies];
    double BulkSpeed[PIC::nTotalSpecies];
    double Lat;
    double wLon;
    double TiltAngle;
    double AzimuthAngle;
    double xLoacation[3];
    double lPointing[3];
    bool ActiveFlag;
  };


  //the total number of the Tiger Stripes and Individual Plumes
  const int nTotalTigerStripes=1;
  const int nTotalIndividualPlumes=1;

  //definition of array of the TigerStriped and Individual plumes
  extern cTigerStripe TigerStripeTable[nTotalTigerStripes];
  extern cIndividualPlume IndividualPlumeTable[nTotalIndividualPlumes];

  //definition of the surface geomenty of the TigerStripes
  extern cTigerStripeGeometry TigerStripeGeometry__ALEXANDRIA;
  extern cTigerStripeGeometry TigerStripeGeometry__BAGDADHOT;
  extern cTigerStripeGeometry TigerStripeGeometry__BAGDAD;
  extern cTigerStripeGeometry TigerStripeGeometry__CAIROHOT;
  extern cTigerStripeGeometry TigerStripeGeometry__CAIRO;
  extern cTigerStripeGeometry TigerStripeGeometry__DAMASCUSHOT;
  extern cTigerStripeGeometry TigerStripeGeometry__DAMASCUS;





  //the boundary conditions on the sphere
  double sphereInjectionRate(int spec,int BoundaryElementType,void *SphereDataPointer);
  int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);


  //generate internal properties and inject particles
  bool GenerateInitialGrainParameters(double *x,double *v,double& GrainRadius, double& GrainWeightCorrection,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode);

  //functions that output the datafile
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);


  inline void Init_BeforeParser() {

    //init the dust model
    ElectricallyChargedDust::minDustRadius=0.01*_MICROMETER_;
    ElectricallyChargedDust::maxDustRadius=100.0*_MICROMETER_;

    ElectricallyChargedDust::Sampling::SetDustSamplingIntervals(10);

    ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity=100.0;
    ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity=600000.0;


    ElectricallyChargedDust::Init_BeforeParser();
  }

  inline void Init_AfterParser() {
    PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
    PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
    PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);

    //init the dust model
    ElectricallyChargedDust::TotalMassDustProductionRate=180.0;
    ElectricallyChargedDust::GenerateNewDustGrainInternalProperties=GenerateInitialGrainParameters;

    ElectricallyChargedDust::SizeDistribution::PowerIndex=4.0;



    ElectricallyChargedDust::Init_AfterParser();
  }

}

#endif /* ENCELADUS_H_ */
