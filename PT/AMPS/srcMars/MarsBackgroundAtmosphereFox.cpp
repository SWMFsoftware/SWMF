/*
 * MarsBackgroundAtmosphereFox.cpp
 *
 *  Created on: Dec 16, 2011
 *      Author: vtenishe
 */

#include "rnd.h"
#include "constants.h"

#include "pic.h"
#include "MarsBackgroundAtmosphereFox.h"


  double MARS_BACKGROUND_ATMOSPHERE_J_FOX_::GetTotalOInjectionRate(double *x) {
    return 2.0*GetTotalO2PlusDissociativeRecombinationRate(x);
  }
