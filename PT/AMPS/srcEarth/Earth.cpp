/*
 * Europa.cpp
 *
 *  Created on: Feb 13, 2012
 *      Author: vtenishe
 */

//$Id$



#include "pic.h"
#include "Earth.h"



char Earth::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="";

//composition of the GCRs
cCompositionGroupTable *Earth::CompositionGroupTable=NULL;
int *Earth::CompositionGroupTableIndex=NULL;
int Earth::nCompositionGroups;


//manager of the data recovering
void Earth::DataRecoveryManager(list<pair<string,list<int> > >& SampledDataRecoveryTable ,int MinOutputFileNumber,int MaxOutputFilenumber) {
  int iOutputFile;
  char fname[_MAX_STRING_LENGTH_PIC_];
  pair<string,list<int> > DataRecoveryTableElement;

  DataRecoveryTableElement.second.push_back(0);

  for (iOutputFile=MinOutputFileNumber;iOutputFile<=MaxOutputFilenumber;iOutputFile++) {
    sprintf(fname,"pic.SamplingDataRestart.out=%i.dat",iOutputFile);

    DataRecoveryTableElement.first=fname;
    SampledDataRecoveryTable.push_back(DataRecoveryTableElement);
  }
}

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;
  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
//do nothing
}

double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {
  double res=1.0;
  ReemissionParticleFraction=0.0;

  if (spec==_O2_SPEC_) res=0.0;

  return res;
}

double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x) {
  return 300.0;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  //do nothing
}

//particle/sphere interactions
int Earth::BC::ParticleSphereInteraction(int spec,long int ptr,double *x,double *v, double &dtTotal, void *NodeDataPonter,void *SphereDataPointer)  {
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//the total injection rate from the Earth
double Earth::BC::sphereInjectionRate(int spec,void *SphereDataPointer) {
  double res=0.0;
  return res;
}

//init the Earth magnetosphere model
void Earth::Init() {
  //init the composition gourp tables
  //!!!!!!!!!!!!! For now only hydrogen is considered !!!!!!!!!!!!!!!!!

  //composition of the GCRs
  nCompositionGroups=1;
  CompositionGroupTable=new cCompositionGroupTable[nCompositionGroups];
  CompositionGroupTableIndex=new int[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) CompositionGroupTableIndex[spec]=0; //all simulated model species are hydrogen
  CompositionGroupTable[0].FistGroupSpeciesNumber=0;
  CompositionGroupTable[0].nModelSpeciesGroup=PIC::nTotalSpecies;

  CompositionGroupTable[0].minVelocity=Relativistic::E2Speed(Earth::BoundingBoxInjection::minEnergy,PIC::MolecularData::GetMass(0));
  CompositionGroupTable[0].maxVelocity=Relativistic::E2Speed(Earth::BoundingBoxInjection::maxEnergy,PIC::MolecularData::GetMass(0));

  CompositionGroupTable[0].GroupVelocityStep=(CompositionGroupTable[0].maxVelocity-CompositionGroupTable[0].minVelocity)/CompositionGroupTable[0].nModelSpeciesGroup;

  //init source models of SEP and GCR
  if (_PIC_EARTH_SEP__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::SEP::Init();
  if (_PIC_EARTH_GCR__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::GCR::Init();
}
