//======================================================================
//$Id$
//======================================================================
//the file contains the user defined setting for compilation of the model

#include "UserDefinition.PIC.dfn"

#include "MarsBackgroundAtmosphereFox.h"

#include "DifferentialCrossSection.h"
extern cDiffCrossSection ForwardCollisionCrossSection;

extern unsigned int _O_SPEC_;
extern int maxLocalBackdroundDensityOffset;

//the plenat
#define _TARGET_ _MARS_

//define the macro for the user-defined mode's header
#define _PIC__USER_DEFINED__USER_PHYSICAL_MODEL_LIST_ "newMars.h"

//BEGIN: the global parameters of the forward scattering cross section
extern double TotalIntegratedForwardScatteringCrossSection;
extern int CumulativeDistributionMaskList;
extern int *CumulativeDistributionMask;
extern int nForwardScatteringCrossSectionLines;

#ifndef _cForwardScatteringCrossSection_
#define _cForwardScatteringCrossSection_
struct cForwardScatteringCrossSection {
  double Angle,DifferentialCrossSection,CumulativeDistributionFunction,deltaCumulativeDistributionFunction;
};
#endif

extern cForwardScatteringCrossSection *ForwardScatteringCrossSectionData;
//END: the global parameters of the forward scattering cross section





