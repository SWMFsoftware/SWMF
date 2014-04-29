//======================================================================
//$Id$
//======================================================================
//the file contains the user defined setting for compilation of the model

//#include "UserDefinition.PIC.dfn"

//the species table
//#define _NA_SPEC_     0
//#define _NA_PLUS_SPEC_ 1

// Europa addition

//possible species used in the Europa model
#ifndef _H2O_SPEC_
#define _H2O_SPEC_             -1
#endif

#ifndef _O2PLUS_SPEC_
#define _O2PLUS_SPEC_         -2
#endif

#ifndef _OPLUS_HIGH_SPEC_
#define _OPLUS_HIGH_SPEC_     -3
#endif

#ifndef _OPLUS_THERMAL_SPEC_
#define _OPLUS_THERMAL_SPEC_  -4
#endif


//extern unsigned int _NA_SPEC_;
extern int maxLocalBackdroundDensityOffset;

//the plenat
#define _TARGET_ _EUROPA_


//define the macro for the user-defined mode's header

//ICES
#undef _PIC_ICES_SWMF_MODE_
#define _PIC_ICES_SWMF_MODE_ _PIC_ICES_MODE_ON_

#undef _PIC_BACKGROUND_ATMOSPHERE_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_MODE_ _PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_

#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_


//uniform time step in the computational domain
#undef _SIMULATION_TIME_STEP_MODE_
#define _SIMULATION_TIME_STEP_MODE_ _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_

//the method for removing the thermalized particles from the system
//#de fine _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__LOCAL_BACKGROUND_THERMAL_SPEED_
#define _THERMALIZED_PARTICLE_REMOVE_CRITERION_ _THERMALIZED_PARTICLE_REMOVE_CRITERION__ESCAPE_SPEED_

//collision cross section type
#define _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION_ _BACKGROUND_ATMOSPHERE_COLLISION_CROSS_SECTION__HARD_SPHERE_

//the volume injection of the particles and the type of injection of model particles into a cell
#undef _PIC_VOLUME_PARTICLE_INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION_MODE_ _PIC_VOLUME_PARTICLE_INJECTION_MODE__OFF_

//the mode of volume injection
#undef _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_
#define _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE_  _PIC_VOLUME_PARTICLE_INJECTION__INJECTION_MODE__RATE_DEPENDENT_

//the distriburion of the collision frequency in a cell
#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_FREQ_MODE__LOCAL_BACKGROUND_DENSITY_


//photolytic reactions
#undef _PIC_PHOTOLYTIC_REACTIONS_MODE_
#define _PIC_PHOTOLYTIC_REACTIONS_MODE_ _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_

#undef _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_
#define _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode) \
    Comet::TotalParticleAcceleration(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);

#undef _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_
#define _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_(xMiddle,vMiddle,spec,ptr,ParticleData,dtMin,TransformationTimeStepLimitFlag,startNode) \
  ElectricallyChargedDust::DustChargingProcessorIndicator(xMiddle,vMiddle,spec,ptr,ParticleData,dtMin,TransformationTimeStepLimitFlag,startNode);

#undef _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_
#define _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xinit,x,v,spec,ptr,ParticleData,dt,startNode) \
  ElectricallyChargedDust::DustChargingProcessor_Implicit_SecondOrder(xinit,x,v,spec,ptr,ParticleData,dt,startNode);

#define _PIC_USER_DEFING_PARTICLE_SAMPLING_(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s) \
  ElectricallyChargedDust::Sampling::SampleParticleData(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s);

#undef _PIC_MODEL__DUST__MODE_
#define _PIC_MODEL__DUST__MODE_ _PIC_MODEL__DUST__MODE__ON_

#undef _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_
#define _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_


//particle sampling


//execute user defined MPI data exchange function after the syncronizatino point at the particle exchange routine
















