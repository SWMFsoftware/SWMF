//Interface to the Geopack package
//$Id$

/*
 * GeopackInterface.h
 *
 *  Created on: Jan 4, 2017
 *      Author: vtenishe
 */

#ifndef _SRC_INTERFACE_GEOPACKINTERFACE_H_
#define _SRC_INTERFACE_GEOPACKINTERFACE_H_

namespace Geopack {

  //IGRF mgnetoc field model
  namespace IGRF {
    void GetMagneticField(double *B,double *x);
  }

  void GetMagneticField(double *B,double *x);
  void Init(const char* Epoch,double *SolarWindVelocity);
}



#endif /* SRC_INTERFACE_GEOPACKINTERFACE_H_ */
