/*
 * plume.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: vtenishe
 */

#include "plume.h"

double DistributionFunction::cosTheta,DistributionFunction::sinTheta,DistributionFunction::plumeBulkVelocity,DistributionFunction::plumeTempetarure;
double DistributionFunction::m_over_k=2.99E-26/1.380658E-23/2.0;

double DistributionFunction::func(double v) {
  return pow(v,2)* exp(-m_over_k/plumeTempetarure*(pow(v*cosTheta-plumeBulkVelocity,2)+pow(v*sinTheta,2)));
}



Cplume *ColumnDensityIntegtant::plume;
double ColumnDensityIntegtant::func(double t);
double *ColumnDensityIntegtant::Ray,*ColumnDensityIntegtant::rInit;
