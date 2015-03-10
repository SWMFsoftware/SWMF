#ifndef _OH_H_
#define _OH_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

// SPICE is not used for this application
// use empty functions to fill in
#include "SpiceEmptyDefinitions.h"

namespace OH {
  using namespace Exosphere;

  namespace Sampling{
    using namespace Exosphere::Sampling;
  }

  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {

    accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

  }
  

}

#endif
