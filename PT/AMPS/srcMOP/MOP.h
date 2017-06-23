//$Id$
//heander for modelingof the Jovian system radiation environemt



/*
 * MOP.h
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#ifndef _SRC_MOP_H_
#define _SRC_MOP_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "pic.h"
#include "Dust.h"
#include "Exosphere.h"
#include "constants.h"

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif

#include "MOP.dfn"


namespace MOP {

  //Saturnian system
  namespace SaturninanSystem {

    namespace Saturn {
    }

    //Enceladus is the main source of water in the Jovian system
    namespace Enceladus {
      //location and velocity of Encelsdus
      extern double xEnceladus[3],vEnceladus[3];

      //the total souce rate
      extern double TotalSourceRate,SourceTemperature;

      //model of the volatiles injection from Enceladus
      //PIC::ParticleWeightTimeStep::fUserDefinedExtraSourceRate PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate
      //PIC::BC::UserDefinedParticleInjectionFunction=
      long int InjectParticles();
      double SourceRate(int spec);
    }

    //the requested local mesh resolution accounted for the topology of the planet system
    double GetLocalMeshResolution(double *x);
  }

  //Jupiter system
  namespace JovianSystem {

    namespace Jupiter {

    }

    namespace Europa {

    }
  }

  //interface for MOP
  double SourceRate(int spec);
  long int InjectParticles();

  //the requested local mesh resolution accounted for the topology of the planet system
  double GetLocalMeshResolution(double *x);
}



#endif /* _SRC_MOP_H_ */
