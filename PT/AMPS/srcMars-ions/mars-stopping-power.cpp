//$Id$
//stopping power model for the solar wind protons, and the pickup ions in the Mars' atmosphere

/*
 * mars-stopping-power.cpp
 *
 *  Created on: May 8, 2017
 *      Author: vtenishe
 */

#include "mars-ions.h"

//tables containing the stopping power data
namespace StoppingPowerData_H_in_CO2 {
  const double emin=1.1*eV2J;
  const double emax=50.0*eV2J;
  const double ReferenceTargetNumberDensity=7.5615E+19*1.0E6;
  const int nDataPoints=123;

  double dlog10e=log10(emax/emin)/nDataPoints;

  double Table[nDataPoints]={4.79E-03,5.00E-03,5.21E-03,5.40E-03,5.59E-03,5.77E-03,5.94E-03,6.10E-03,6.41E-03,6.77E-03,7.11E-03,7.42E-03,7.72E-03,8.00E-03,8.26E-03,
    8.51E-03,8.74E-03,9.18E-03,9.59E-03,9.97E-03,1.03E-02,1.06E-02,1.10E-02,1.15E-02,1.20E-02,1.25E-02,1.29E-02,1.33E-02,1.36E-02,1.40E-02,1.43E-02,1.46E-02,
    1.48E-02,1.51E-02,1.56E-02,1.61E-02,1.66E-02,1.70E-02,1.73E-02,1.77E-02,1.80E-02,1.83E-02,1.85E-02,1.90E-02,1.94E-02,1.97E-02,2.00E-02,2.03E-02,2.06E-02,
    2.09E-02,2.13E-02,2.15E-02,2.17E-02,2.19E-02,2.20E-02,2.21E-02,2.22E-02,2.23E-02,2.23E-02,2.23E-02,2.23E-02,2.23E-02,2.23E-02,2.22E-02,2.20E-02,2.19E-02,
    2.18E-02,2.16E-02,2.14E-02,2.11E-02,2.07E-02,2.04E-02,2.01E-02,1.97E-02,1.94E-02,1.88E-02,1.82E-02,1.76E-02,1.71E-02,1.66E-02,1.61E-02,1.57E-02,1.53E-02,
    1.49E-02,1.46E-02,1.42E-02,1.36E-02,1.29E-02,1.23E-02,1.17E-02,1.12E-02,1.08E-02,1.04E-02,1.00E-02,9.66E-03,9.05E-03,8.52E-03,8.06E-03,7.65E-03,7.29E-03,
    6.97E-03,6.41E-03,5.95E-03,5.55E-03,5.22E-03,4.92E-03,4.66E-03,4.43E-03,4.23E-03,4.04E-03,3.87E-03,3.72E-03,3.45E-03,3.17E-03,2.94E-03,2.74E-03,2.57E-03,
    2.42E-03,2.29E-03,2.18E-03,2.07E-03,1.89E-03,1.75E-03};
}

namespace StoppingPowerData_O_in_CO2 {
  const double emin=1.1*eV2J;
  const double emax=50.0*eV2J;
  const double ReferenceTargetNumberDensity=7.5615E+19*1.0E6;
  const int nDataPoints=123;

  double dlog10e=log10(emax/emin)/nDataPoints;

  double Table[nDataPoints]={4.54E-02,4.80E-02,5.04E-02,5.28E-02,5.51E-02,5.73E-02,5.95E-02,6.16E-02,6.57E-02,7.06E-02,7.52E-02,7.96E-02,8.39E-02,8.79E-02,
    9.19E-02,9.57E-02,9.93E-02,1.06E-01,1.13E-01,1.19E-01,1.25E-01,1.31E-01,1.37E-01,1.47E-01,1.57E-01,1.67E-01,1.76E-01,1.84E-01,1.92E-01,2.00E-01,
    2.07E-01,2.14E-01,2.21E-01,2.28E-01,2.40E-01,2.55E-01,2.69E-01,2.82E-01,2.94E-01,3.06E-01,3.17E-01,3.28E-01,3.38E-01,3.57E-01,3.74E-01,3.91E-01,
    4.06E-01,4.21E-01,4.35E-01,4.60E-01,4.83E-01,5.05E-01,5.24E-01,5.43E-01,5.60E-01,5.76E-01,5.91E-01,6.05E-01,6.19E-01,6.32E-01,6.56E-01,6.83E-01,
    7.07E-01,7.29E-01,7.50E-01,7.69E-01,7.86E-01,8.02E-01,8.17E-01,8.45E-01,8.69E-01,8.91E-01,9.10E-01,9.28E-01,9.44E-01,9.72E-01,9.97E-01,1.02E+00,
    1.04E+00,1.05E+00,1.07E+00,1.08E+00,1.09E+00,1.10E+00,1.11E+00,1.11E+00,1.13E+00,1.14E+00,1.15E+00,1.16E+00,1.16E+00,1.16E+00,1.17E+00,1.17E+00,
    1.17E+00,1.17E+00,1.16E+00,1.16E+00,1.15E+00,1.14E+00,1.14E+00,1.12E+00,1.10E+00,1.08E+00,1.06E+00,1.04E+00,1.03E+00,1.01E+00,9.91E-01,9.75E-01,
    9.59E-01,9.43E-01,9.14E-01,8.80E-01,8.48E-01,8.19E-01,7.93E-01,7.68E-01,7.45E-01,7.23E-01,7.03E-01,6.67E-01,6.35E-01};
}


double PIC::MolecularCollisions::StoppingPowerModel::GetStoppingPower(double *x,double *v,int spec,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double res;

  class cInterpolator {
  public:
    double emin,emax,dLog10E,ReferenceTargetNumberDensity;
    double *Table;
    int nPoints;

    cInterpolator(double eminIn,double emaxIn,int nPointsIn,double ReferenceTargetNumberDensityIn,double *TableIn) {
      emin=eminIn,emax=emaxIn,nPoints=nPointsIn,Table=TableIn,ReferenceTargetNumberDensity=ReferenceTargetNumberDensityIn;
      dLog10E=log10(emax/emin)/nPoints;
    }

    double EvaluateLocation(double *x,double *v,int spec) {
      double LocalAtmosphereDensity,E;
      int iLevel;

      E=0.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*PIC::MolecularData::GetMass(spec);

      if (E<emin) return 0.0;
      else if (E>=emax) iLevel=nPoints-1;
      else {
        iLevel=(int)(log10(E/emin)/dLog10E);
      }

      if (iLevel<0) iLevel=0;
      if (iLevel>=nPoints) iLevel=nPoints-1;

      LocalAtmosphereDensity=MarsIon::GetBackgroundAtmosphereDensity(x,_CO2_SPEC_);

      return Table[iLevel]*(MeV2J*100)*LocalAtmosphereDensity/ReferenceTargetNumberDensity;
    }
  };

  static cInterpolator H_in_CO2(StoppingPowerData_H_in_CO2::emin,StoppingPowerData_H_in_CO2::emax,StoppingPowerData_H_in_CO2::nDataPoints,StoppingPowerData_H_in_CO2::ReferenceTargetNumberDensity,StoppingPowerData_H_in_CO2::Table);
  static cInterpolator O_in_CO2(StoppingPowerData_O_in_CO2::emin,StoppingPowerData_O_in_CO2::emax,StoppingPowerData_O_in_CO2::nDataPoints,StoppingPowerData_O_in_CO2::ReferenceTargetNumberDensity,StoppingPowerData_O_in_CO2::Table);

  switch (spec) {
  case _O_PLUS_SPEC_:
    res=O_in_CO2.EvaluateLocation(x,v,spec);
    break;

  case _H_PLUS_SPEC_:
    res=H_in_CO2.EvaluateLocation(x,v,spec);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the species is unknown");
  }

  return res;
}




