
//$Id$

/*
 * OrbitCalculation.cpp
 *
 *  Created on: Jul 12, 2015
 *      Author: vtenishe
 */


#include "pic.h"
#include "Exosphere.h"

double Exosphere::OrbitalMotion::GetTAA(double t) {return 0.0;}
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;} 


