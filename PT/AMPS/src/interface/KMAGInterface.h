//$Id$
//interface to Khurana model of magnetis field in Saturn's magnetosphre

/*
 * KMAGInterface.h
 *
 *  Created on: Aug 1, 2017
 *      Author: vtenishe
 */

#ifndef _INTERFACE_KMAGINTERFACE_H_
#define _INTERFACE_KMAGINTERFACE_H_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "constants.h"

namespace KMAG {

  //time and epoch
  extern double et;
  extern char Epoch[100];
  extern char FRAME[100];

  //solar wind magnetic field and dynamic pressure
  namespace SolarWind {
    extern double IMF[3],DynamicPressure;
  }

  void SetEpochTime(double etIn,const char* EpochIn);
  void GetMagneticField(double *B,double *x);
}



#endif /* _INTERFACE_KMAGINTERFACE_H_ */
