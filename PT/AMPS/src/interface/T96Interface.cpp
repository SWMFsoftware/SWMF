//$Id$
//Interface to the T96 model

/*
 * T96Interface.cpp
 *
 *  Created on: Jan 5, 2017
 *      Author: vtenishe
 */


#include "constants.h"
#include "constants.PlanetaryData.h"
#include "T96Interface.h"

/* 96::PARMOD:
 * c--------------------------------------------------------------------
C DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
C           (2) DST (NANOTESLA),  (3) BYIMF, AND (4) BZIMF (NANOTESLA).
c THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 4 ELEMENTS
c OF THE ARRAY PARMOD(10).
 *
 */

double T96::PS=0.170481;
double T96::PARMOD[11]={2.0,-50.0,1.0,-3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

extern "C"{
  void t96_01_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

void T96::GetMagneticField(double *B,double *x) {
  int idim,IOPT=0;

  //calculate the global magnetic field in the magnetosphere
  double xLocal[3],bT96[3];

  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;
  t96_01_(&IOPT,PARMOD,&PS,xLocal+0,xLocal+1,xLocal+2,bT96+0,bT96+1,bT96+2);

  //calcualte the Earth's internal magnetis field
  IGRF::GetMagneticField(B,x);

  //sum contributions of the internal Earth's magnetic field, and the global field in the magnetoshere
  for (idim=0;idim<3;idim++) B[idim]+=bT96[idim]*_NANO_;
}

void T96::SetSolarWindPressure(double SolarWindPressure) {PARMOD[0]=SolarWindPressure/_NANO_;}
void T96::SetDST(double DST) {PARMOD[1]=DST/_NANO_;}
void T96::SetBYIMF(double BYIMF) {PARMOD[2]=BYIMF/_NANO_;}
void T96::SetBZIMF(double BZIMF) {PARMOD[3]=BZIMF/_NANO_;}



