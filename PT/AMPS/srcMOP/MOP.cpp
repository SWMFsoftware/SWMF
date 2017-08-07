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

//init the model
void MOP::Init() {
  //set the printout procedures into the core
  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(Sampling::Output::PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(Sampling::Output::PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Sampling::Output::Interpolate);

  //get Saturn's axis of rotation
  SpiceDouble lz_IAU_SATURN[3]={0.0,0.0,1.0},rotate[3][3];

  utc2et_c(Exosphere::SimulationStartTimeString,&KMAG::et);
  pxform_c("IAU_SATURN","SSO",KMAG::et,rotate);
  mxv_c(rotate,lz_IAU_SATURN,MOP::SaturninanSystem::Saturn::RotationAxis);
}

//output of the sampled data into a file
void MOP::Sampling::Output::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Corotation Speed\"");
}

void MOP::Sampling::Output::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  //empty function. There is nothing to interpolate yet
  //only parameters derived from empirical modes are printed into a file
}

void MOP::Sampling::Output::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

  if (PIC::ThisThread==0) {
    //fprintf(fout," %e ",CenterNode->GetX()[0]);
    fprintf(fout," %e ",MOP::SaturninanSystem::Magnetosphere::GetCorotationSpeed(CenterNode->GetX()));
  }
}















