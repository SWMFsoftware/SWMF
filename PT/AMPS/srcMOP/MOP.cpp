//$Id$
//model of the outer planet magnetospheres


/*
 * MOP.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "MOP.h"

char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Enceladus";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="SSO";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_SATURN";
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {return 0.0;}
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {return 0;}
void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {}
void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {}
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {return 0.0;}
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x) {return 0.0;}


//new particle injection functions
long int MOP::InjectParticles() {
  int res=0;

  switch (_PLANET_SYSTEM_) {
  case _SATURNIAN_SYSTEM_:
    res=MOP::SaturninanSystem::Enceladus::InjectParticles();
    break;
  default:
    exit(__LINE__,__FILE__,"Error: not implemented");
  }

  return res;
}

double MOP::SourceRate(int spec) {
  double res=0.0;

  switch (_PLANET_SYSTEM_) {
  case _SATURNIAN_SYSTEM_:
    res=MOP::SaturninanSystem::Enceladus::SourceRate(spec);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: not implemented");
  }

  return res;
}

//the requested local mesh resolution accounted for the topology of the planet system
double MOP::GetLocalMeshResolution(double *x) {
  double res;

  switch (_PLANET_SYSTEM_) {
  case _SATURNIAN_SYSTEM_:
    res=MOP::SaturninanSystem::GetLocalMeshResolution(x);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: not implemented");
  }

  return res;
}


