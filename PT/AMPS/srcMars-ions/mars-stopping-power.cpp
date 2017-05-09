//$Id$
//stopping power model for the solar wind protons, and the pickup ions in the Mars' atmosphere

/*
 * mars-stopping-power.cpp
 *
 *  Created on: May 8, 2017
 *      Author: vtenishe
 */

#include "mars-ions.h"

double PIC::MolecularCollisions::StoppingPowerModel::GetStoppingPower(double *x,double *v,int spec,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double E;

  E=0.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*PIC::MolecularData::GetMass(spec);

  static const int nDataPoints=79;
  static const double minE=10.0E3*KeV2J;
  static const double maxE=10.0*MeV2J;
  static const double dLogE=log10(maxE/minE)/nDataPoints;

  static const double BackgroundTargetReferenceDensity=3.0834E+22;

  static const double Table[nDataPoints] {
      8.40E-03,
      7.85E-03,
      7.37E-03,
      6.96E-03,
      6.60E-03,
      6.27E-03,
      5.98E-03,
      5.72E-03,
      5.48E-03,
      5.07E-03,
      4.64E-03,
      4.28E-03,
      3.98E-03,
      3.72E-03,
      3.50E-03,
      3.30E-03,
      3.13E-03,
      2.97E-03,
      2.71E-03,
      2.49E-03,
      2.31E-03,
      2.15E-03,
      2.02E-03,
      1.90E-03,
      1.71E-03,
      1.55E-03,
      1.42E-03,
      1.32E-03,
      1.22E-03,
      1.15E-03,
      1.08E-03,
      1.02E-03,
      9.64E-04,
      9.17E-04,
      8.74E-04,
      8.00E-04,
      7.24E-04,
      6.63E-04,
      6.11E-04,
      5.68E-04,
      5.30E-04,
      4.98E-04,
      4.69E-04,
      4.44E-04,
      4.01E-04,
      3.67E-04,
      3.38E-04,
      3.13E-04,
      2.92E-04,
      2.74E-04,
      2.44E-04,
      2.20E-04,
      2.01E-04,
      1.85E-04,
      1.71E-04,
      1.60E-04,
      1.50E-04,
      1.41E-04,
      1.33E-04,
      1.26E-04,
      1.20E-04,
      1.09E-04,
      9.84E-05,
      8.96E-05,
      8.24E-05,
      7.62E-05,
      7.10E-05,
      6.65E-05,
      6.25E-05,
      5.90E-05,
      5.31E-05,
      4.83E-05,
      4.44E-05,
      4.10E-05,
      3.82E-05,
      3.57E-05,
      3.17E-05,
      2.85E-05,
      2.59E-05};

  int iLevel;

  if (E<minE) iLevel=0;
  else if (E>=maxE) iLevel=nDataPoints-1;
  else {
    iLevel=(int)(log10(E/minE)/dLogE);
  }

  if (iLevel<0) iLevel=0;
  if (iLevel>=nDataPoints) iLevel=nDataPoints-1;

  //get local background atmosphere density
  double LocalAtmosphereDensity;


  double res=0.0;
  double altitude=0.0;
  double rSeason=1.38758;
  double r2Season=pow(rSeason,2);
  double Tnu_body = 134.0; //in K, neutral temperature
  double BodynDenNuSp_I_O= 8.0283e15; // per m^-3
  double BodynDenNuSp_I_Ox= 5.1736e14;
  double BodynDenNuSp_I_Oh= 6.3119e10;
  double BodynDenNuSp_I_Ohx= 3.9646e9;
  double HNuSpecies_I_O=13340; //scale height in m
  double HNuSpecies_I_Ox=50025;
  double HNuSpecies_I_Oh=290500;
  double HNuSpecies_I_Ohx=2436600;
  double Rate_I_O_Op= 6.346e-7/r2Season; //units s^-1 for O_hv__Op_em_
  double nDen_O;
  int i;

  for (i=0;i<3;i++) {
      altitude+=pow(x[i],2);
  }

  if (sqrt(altitude)<3.5E6) return 0.0;

  altitude=sqrt(altitude)-3396000.0; // in m where R_Mars=3396000.0 m

  if (altitude<0.0) return 0.0;

  nDen_O=BodynDenNuSp_I_O*exp(-altitude/HNuSpecies_I_O)+\
         BodynDenNuSp_I_Ox*exp(-altitude/HNuSpecies_I_Ox)+\
         BodynDenNuSp_I_Oh*exp(-altitude/HNuSpecies_I_Oh)+\
         BodynDenNuSp_I_Ohx*exp(-altitude/HNuSpecies_I_Ohx);


  return Table[iLevel]*(MeV2J*100)*LocalAtmosphereDensity/BackgroundTargetReferenceDensity;
}




