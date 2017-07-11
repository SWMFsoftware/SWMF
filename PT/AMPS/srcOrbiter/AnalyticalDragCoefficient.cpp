//$Id$
//closed form expressions for the drag coefficient


/*
 * AnalyticalDragCoefficient.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: vtenishe
 */

#include "Orbiter.h"

//analytical drag coefficient of a sphere (Metha-2013-ASR,Mehta-2014-JSR)
double Orbiter::DragCoefficent::ClosedFormExpressions::AccommodationCoefficient::Sphere(double Mass,double Speed,double AccommodationCoefficient,double Tw,double Tinf) {
  double s,s2,s3,s4,T_kr,T_ki,v_mp;

  T_ki=Mass*pow(Speed,2)/(3.0*Kbol);
  T_kr=(1.0-AccommodationCoefficient)*T_ki+AccommodationCoefficient*Tw;

  v_mp=sqrt(2.0*Kbol*Tinf/Mass);
  s=Speed/v_mp;

  s2=s*s;
  s3=s2*s;
  s2=s2*s2;

  return (4.0*s4+4.0*s2-1.0)/(2.0*s4)*erf(s)+(2.0*s2+1.0)/(sqrtPi*s3)*exp(-s2)+(2.0*sqrtPi)/(3.0*s)*sqrt(T_kr/Tinf);
}

double Orbiter::DragCoefficent::ClosedFormExpressions::AccommodationCoefficient::Plate(double Mass,double Speed,double AccommodationCoefficient,double Tw,double Tinf) {
  double s,s2,s3,s4,T_kr,T_ki,v_mp;

  T_ki=Mass*pow(Speed,2)/(3.0*Kbol);
  T_kr=(1.0-AccommodationCoefficient)*T_ki+AccommodationCoefficient*Tw;

  v_mp=sqrt(2.0*Kbol*Tinf/Mass);
  s=Speed/v_mp;

  s2=s*s;
  s3=s2*s;
  s2=s2*s2;

  return (2.0+1.0/s2)*erf(s)+2.0/(sqrtPi*s)*exp(-s2)+sqrtPi/s*sqrt(T_kr/Tinf);
}

double Orbiter::DragCoefficent::ClosedFormExpressions::AccommodationCoefficient::Cylinder(double Mass,double Speed,double AccommodationCoefficient,double Tw,double Tinf,double Length,double Diameter) {
  double res,s,s2,s3,s4,T_kr,T_ki,v_mp;

  T_ki=Mass*pow(Speed,2)/(3.0*Kbol);
  T_kr=(1.0-AccommodationCoefficient)*T_ki+AccommodationCoefficient*Tw;

  v_mp=sqrt(2.0*Kbol*Tinf/Mass);
  s=Speed/v_mp;

  s2=s*s;
  s3=s2*s;
  s2=s2*s2;

  return (2.0+1.0/s2)*erf(s)+2.0/(sqrtPi*s)*exp(-s2)+sqrtPi/s*sqrt(T_kr/Tinf) +4.0/(sqrtPi*s)*Length/Diameter;
}




