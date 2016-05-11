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

  //definitions of the Tiger Stripes and Individual plume classes
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
  static const int nTotalTigerStripes=1;
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

  //convert the (Lat, wLon) location of the plumes into Cartesian coordinates
  void RecalculateSourceLocations();

  //source model of the Enceladus' exosphere
  namespace SourceModel {

    void Init();

    namespace IndividualPlumes {
      extern double TotalSourceRateTable[PIC::nTotalSpecies];
//      extern double TotalPlumeSourceRateTable[nTotalIndividualPlumes][PIC::nTotalSpecies];
      extern double maxTotalPlumeSourceRateTable[PIC::nTotalSpecies];

      bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
          double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
    }

    namespace TigerStripes {
      extern double TotalSourceRateTable[PIC::nTotalSpecies];
//      extern double TotalTigerStripeSourceRateTable[nTotalTigerStripes][PIC::nTotalSpecies];
      extern double maxTotalTigerStripeSourceRateTable[PIC::nTotalSpecies];

      //Information needed for distribution of the injeciton location on a Tiger Stripe
      extern int *TigerStripeSegmentsNumberTable;
      extern double **TigerStripeSegmentsLengthTable;
      extern double *maxTigerStripeSegmentsLengthTable;

      bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
          double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
    }

    //interface to AMPS
    double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer);
    bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
        double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
  }


  //Acceleration of the particles
  void TotalAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);

  //the boundary conditions on the sphere
  double sphereInjectionRate(int spec,int BoundaryElementType,void *SphereDataPointer);
  int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);


  //generate internal properties and inject particles
  bool GenerateInitialGrainParameters(double *x,double *v,double& GrainRadius, double& GrainWeightCorrection,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode);

  //functions that output the datafile
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);


  //Init MultiPlume Model
  void Init_BeforeParser();
  void Init_AfterParser();

}

#endif /* ENCELADUS_H_ */
