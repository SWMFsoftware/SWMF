//Physical model of SEP at the boundary of the computational domain
//$Id$

/*
 * BoundaryInjection_SEP.cpp
 *
 *  Created on: Oct 19, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include "Earth.h"

double Earth::BoundingBoxInjection::SEP::InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  if (nface!=1) {
    //all particles are injected in the anti-x direction
    return 0.0;
  }

  return 1.0;
}

void Earth::BoundingBoxInjection::SEP::GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double* x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal) {

  static const int nPoints=40;
  static const double EnergySpectrum[nPoints]={79689615981,67315760388,56807706047,47291021568,39943321713,32884523614,26817594563,21823180644,17338648107,14058924756,11102281687,
      8766389207,6838407502,5348068972,4150749740,3123281219,2409638480,1809604607,1392822418,1062211119,774325123,576246297,430250451.2,318204119.4,235098989.3,
      173236650,124923105.5,91664691.59,63573704.25,42179163.5,29175851.43,19139990.68,12127805.66,7982215.686,4707100.297,2087487.85,447284.6004,83203.81316};

  static double emin=100000.00*eV2J;
  static double emax=623550805.17*eV2J;
  static double de=(log(emax/emin)/(nPoints-1));

  double Speed,e,WeightCorrection=1.0,v[3];
  int iLevel;

  //distribute energy of the new particle
  e=emin+rnd()*(emax-emin);
  iLevel=(int)(log(e/emin)/de);

  if (iLevel<0) iLevel=0;
  if (iLevel>=nPoints) iLevel=nPoints-1;
  WeightCorrection=EnergySpectrum[iLevel]/EnergySpectrum[0];

  Speed=Relativistic::E2Speed(e,PIC::MolecularData::GetMass(spec));
  v[0]=-Speed,v[1]=0.0,v[2]=0.0;

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,ParticleData);
}

//init the SEP injection model
void Earth::BoundingBoxInjection::SEP::Init() {
  //nothing to do
}
