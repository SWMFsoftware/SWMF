//$Id$
//interface to T96 model

/*
 * T96Interface.h
 *
 *  Created on: Jan 5, 2017
 *      Author: vtenishe
 */

#ifndef _INTERFACE_T96INTERFACE_H_
#define _INTERFACE_T96INTERFACE_H_

#include "GeopackInterface.h"

namespace T96 {
  using namespace Geopack;

  //dipole tilt angle
  extern double PS;

  //the set of themodel paramters
  extern double PARMOD[11];
  void SetSolarWindPressure(double SolarWindPressure);
  void SetDST(double DST);
  void SetBYIMF(double BYIMF);
  void SetBZIMF(double BZIMF);


  void GetMagneticField(double *B,double *x);


}



#endif /* SRC_INTERFACE_T96INTERFACE_H_ */
