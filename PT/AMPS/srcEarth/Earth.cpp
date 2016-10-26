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
  //init source models of SEP and GCR
  if (_PIC_EARTH_SEP__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::SEP::Init();
  if (_PIC_EARTH_GCR__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::GCR::Init();
}
