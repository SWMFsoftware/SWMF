/*
 * mars-ion_exosphere.cpp
 *
 *  Created on: May 26, 2015
 *      Author: vtenishe
 */

#include "mars-ions.h"

char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Mars";

//calculation of the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {return 0.0;}

//column integration
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}

//particle/surface interaction
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 1.0;}
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {return 100.0;}

//the name of the simulation frame of reference
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EOJ";




