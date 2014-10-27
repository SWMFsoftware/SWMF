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

#undef _PIC_BACKGROUND_ATMOSPHERE_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_MODE_ _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_

#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_


//the method for removing the thermalized particles from the system
#define _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__LOCAL_BACKGROUND_THERMAL_SPEED_
//#define _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__ESCAPE_SPEED_

//collision cross section type
#define _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__HARD_SPHERE_

//the volume injection of the particles and the type of injection of model particles into a cell
#undef _PIC_VOLUME_PARTICLE_INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION_MODE_ _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_

//the mode of volume injection
#undef _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_  _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE__RATE_DEPENDENT_

//the distriburion of the collision frequency in a cell
#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_


//the model for the baclground atmosphere
//define _MARS_BACKGROUND_ATMOSPHERE_MODEL_ _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_
#define _MARS_BACKGROUND_ATMOSPHERE_MODEL_ _MARS_BACKGROUND_ATMOSPHERE_MODEL__MTGCM_

//the mode of counting the escaping particles:
#define _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ _MARS_ESCAPE_PARTICLES_COUNTING_MODE__ESCAPE_SPEED_


#undef _PIC_DYNAMIC_LOAD_BALANCING_MODE_
//#define _PIC_DYNAMIC_LOAD_BALANCING_MODE_  _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
#define _PIC_DYNAMIC_LOAD_BALANCING_MODE_ _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

#undef _PIC_DEBUGGER_MODE_
#define _PIC_DEBUGGER_MODE_ _PIC_DEBUGGER_MODE_ON_

#undef _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_
#define _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode) \
    newMars::TotalParticleAcceleration(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);


//include the header of the Mars model
#undef _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_
#define _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__ON_


#undef _PIC_MEMORY_PREFETCH_MODE_
#define _PIC_MEMORY_PREFETCH_MODE_  _PIC_MEMORY_PREFETCH_MODE__OFF_


#define _PIC_USER_DEFING_PARTICLE_SAMPLING_(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s) \
    newMars::SampleParticleData(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s);




