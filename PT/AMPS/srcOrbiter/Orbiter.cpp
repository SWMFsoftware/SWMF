/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

//the object name and the names of the frames
/*char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Orbiter";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";*/

double Orbiter::UpstreamBC::Velocity[3]={-30.0E3,0.0,0.0};
double Orbiter::UpstreamBC::NumberDensity=1.5E21;
double Orbiter::UpstreamBC::Temperature=293.0;





/*double Exosphere::SurfaceInteraction::StickingProbability(int spec, double& ReemissionParticleFraction,double Temp) {
  double res=0.0;

   return res;
}


//surface temeprature of the planet
double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
  return 180;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int nVariables=0;
  return nVariables;
  }

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  double res=0.0;
  return res;
}*/


