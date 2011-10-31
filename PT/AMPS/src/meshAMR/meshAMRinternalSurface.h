//=======================================================================
//$Id$
//=======================================================================
//the file combines the definitions for the internal surfaces

#ifndef _AMR_INTERNAL_SURFACE_
#define _AMR_INTERNAL_SURFACE_

#include "math.h"


#include "meshAMRdef.h"
#include "mpichannel.h"

//include the user defined data for the internal boundaries
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
#define _LOAD_USER_DEFINITIONS_
#endif

#ifdef _LOAD_USER_DEFINITIONS_
#include "UserDefinition.meshAMR.h"
#endif

#include "meshAMRinternalSurface_sphere.h"
#include "meshAMRinternalSurface_circle.h"
#include "meshAMRinternalSurface_1D_sphere.h"


//=======================================================================
//the descriptor of the internal boundary conditions
class cInternalBoundaryConditionsDescriptor {
public:
  unsigned char BondaryType;
  void *BoundaryElement;
  cInternalBoundaryConditionsDescriptor *nextInternalBCelement;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  void cleanDataBuffer() {
    BondaryType=_INTERNAL_BOUNDARY_TYPE_UNDEFINED_;
    BoundaryElement=NULL,nextInternalBCelement=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalBoundaryConditionsDescriptor() {
    cleanDataBuffer();
  }
};

#endif
