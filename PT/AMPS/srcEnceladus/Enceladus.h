/*
 * Enceladus.h
 *
 *  Created on: Jan 17, 2012
 *      Author: vtenishe
 */

//$Id$

#ifndef ENCELADUS_H_
#define ENCELADUS_H_

#include "plume.h"
#include "pic.h"
#include "pic__model__electrically_charged_dust.h"

#define nTotalPlumes 8


namespace Enceladus {
using namespace ElectricallyChargedDust;



 // static const double DustGrainRadius=1.0E-4;


  //parameters of the sources
  extern double sphericalSourceProductionRateConst;
  extern double sphericalSourceExpansionVelocityConst;
  extern double sphericalSource_Const0;
  extern double sphericalSource_Const1;

  extern Cplume Plume[nTotalPlumes];

  //the boundary conditions on the sphere
  double sphereInjectionRate(int spec,void *SphereDataPointer);
  int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);


  //generate internal properties and inject particles
  bool GenerateInitialGrainParameters(double *x,double *v,double& GrainRadius, double& GrainWeightCorrection,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode);

  //functions that output the datafile
  void PrintVariableList(FILE* fout,int DataSetNumber);
  void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);


  inline void Init_BeforeParser() {
    sphericalSourceProductionRateConst=4.683021e+25;
    sphericalSourceExpansionVelocityConst=6.000000e+02;
    sphericalSource_Const0=6.210543e+21;
    sphericalSource_Const1=4.705504e+09;


    Plume[0].Temperature=7.000000e+01;
    Plume[0].BulkVelocity=6.000000e+02;
    Plume[0].SourceRate=1.390928e+25;
    Plume[0].plumeDirection[0]=2.904524e-02;
    Plume[0].plumeDirection[1]=-1.103317e-01;
    Plume[0].plumeDirection[2]=-9.934703e-01;
    Plume[0].plumePosition[0]=3.160772e+04;
    Plume[0].plumePosition[1]=-1.914232e+04;
    Plume[0].plumePosition[2]=-2.472540e+05;


    Plume[1].Temperature=7.000000e+01;
    Plume[1].BulkVelocity=6.000000e+02;
    Plume[1].SourceRate=6.549067e+27;
    Plume[1].plumeDirection[0]=4.791380e-02;
    Plume[1].plumeDirection[1]=2.666430e-01;
    Plume[1].plumeDirection[2]=-9.626037e-01;
    Plume[1].plumePosition[0]=3.206783e+04;
    Plume[1].plumePosition[1]=3.414877e+04;
    Plume[1].plumePosition[2]=-2.455718e+05;


    Plume[2].Temperature=7.000000e+01;
    Plume[2].BulkVelocity=6.000000e+02;
    Plume[2].SourceRate=2.815701e+25;
    Plume[2].plumeDirection[0]=-3.825870e-01;
    Plume[2].plumeDirection[1]=4.092247e-01;
    Plume[2].plumeDirection[2]=-8.283492e-01;
    Plume[2].plumePosition[0]=1.567810e+04;
    Plume[2].plumePosition[1]=3.488536e+04;
    Plume[2].plumePosition[2]=-2.470571e+05;


    Plume[3].Temperature=7.000000e+01;
    Plume[3].BulkVelocity=6.000000e+02;
    Plume[3].SourceRate=3.594242e+25;
    Plume[3].plumeDirection[0]=-2.062738e-01;
    Plume[3].plumeDirection[1]=-1.826742e-01;
    Plume[3].plumeDirection[2]=-9.612915e-01;
    Plume[3].plumePosition[0]=-6.154404e+04;
    Plume[3].plumePosition[1]=-3.786215e+04;
    Plume[3].plumePosition[2]=-2.393299e+05;


    Plume[4].Temperature=7.000000e+01;
    Plume[4].BulkVelocity=6.000000e+02;
    Plume[4].SourceRate=1.846467e+25;
    Plume[4].plumeDirection[0]=-5.238592e-02;
    Plume[4].plumeDirection[1]=-1.319658e-01;
    Plume[4].plumeDirection[2]=-9.898690e-01;
    Plume[4].plumePosition[0]=1.464897e+04;
    Plume[4].plumePosition[1]=-4.674493e+04;
    Plume[4].plumePosition[2]=-2.451537e+05;


    Plume[5].Temperature=7.000000e+01;
    Plume[5].BulkVelocity=6.000000e+02;
    Plume[5].SourceRate=8.072171e+26;
    Plume[5].plumeDirection[0]=1.727698e-01;
    Plume[5].plumeDirection[1]=-7.120150e-02;
    Plume[5].plumeDirection[2]=-9.823853e-01;
    Plume[5].plumePosition[0]=-6.888722e+03;
    Plume[5].plumePosition[1]=1.060770e+04;
    Plume[5].plumePosition[2]=-2.496798e+05;


    Plume[6].Temperature=7.000000e+01;
    Plume[6].BulkVelocity=6.000000e+02;
    Plume[6].SourceRate=2.281627e+25;
    Plume[6].plumeDirection[0]=5.026898e-01;
    Plume[6].plumeDirection[1]=-2.431387e-01;
    Plume[6].plumeDirection[2]=-8.295701e-01;
    Plume[6].plumePosition[0]=5.775287e+04;
    Plume[6].plumePosition[1]=-3.188130e+04;
    Plume[6].plumePosition[2]=-2.411394e+05;


    Plume[7].Temperature=7.000000e+01;
    Plume[7].BulkVelocity=6.000000e+02;
    Plume[7].SourceRate=1.352574e+25;
    Plume[7].plumeDirection[0]=4.960674e-02;
    Plume[7].plumeDirection[1]=-1.221694e-01;
    Plume[7].plumeDirection[2]=-9.912688e-01;
    Plume[7].plumePosition[0]=-1.479285e+04;
    Plume[7].plumePosition[1]=-3.101386e+04;
    Plume[7].plumePosition[2]=-2.476274e+05;

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
