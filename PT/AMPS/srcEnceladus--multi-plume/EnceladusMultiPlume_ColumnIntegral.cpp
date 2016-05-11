//$Id$

/*
 * EnceladusMultiPlume_ColumnIntegral.cpp
 *
 *  Created on: May 6, 2016
 *      Author: vtenishe
 */


#include "pic.h"
#include "EnceladusMultiPlume.h"

//double Exosphere::OrbitalMotion::GetTAA(double t) {return 0.0;}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;}

double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {return 100.0;}

int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int s,nvars=0;

  return nvars;
}


void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {

}
