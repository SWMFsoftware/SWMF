//======================================================================
//$Id$
//======================================================================
//the file contains the user defined setting for compilation of the model

#include "UserDefinition.PIC.dfn"


//extern unsigned int _DUST_SPEC_;
extern int maxLocalBackdroundDensityOffset;

//the planet
#define _TARGET_ _ENCELADUS_

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
#define _PIC_BACKGROUND_ATMOSPHERE_MODE_ _PIC_BACKGROUND_ATMOSPHERE_MODE__OFF_

#undef _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_
#define _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION_ _PIC_BACKGROUND_ATMOSPHERE_COLLISION_VELOCITY_REDISTRIBUTION__ISOTROPIC_


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


//the model for the baclground atmosphere
#define _MARS_BACKGROUND_ATMOSPHERE_MODEL_ _MARS_BACKGROUND_ATMOSPHERE_MODEL__FOX_

//the mode of counting the escaping particles:
#define _MARS_ESCAPE_PARTICLES_COUNTING_MODE_ _MARS_ESCAPE_PARTICLES_COUNTING_MODE__ESCAPE_SPEED_

//read the plasma data from ICES
#undef _PIC_ICES_SWMF_MODE_
#define _PIC_ICES_SWMF_MODE_ _PIC_ICES_MODE_ON_

//allow dust grain's charging
#undef _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_
#define _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_  _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_

//re-define functions that contrals the particle motion
#undef _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_
#define _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_INDICATOR_(xMiddle,vMiddle,spec,ptr,ParticleData,dtMin,TransformationTimeStepLimitFlag,startNode) \
    ElectricallyChargedDust::DustChargingProcessorIndicator(xMiddle,vMiddle,spec,ptr,ParticleData,dtMin,TransformationTimeStepLimitFlag,startNode);

#undef _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_
#define _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode) \
    ElectricallyChargedDust::TotalGrainAcceleration(acclMiddle,spec,ptr,xMiddle,vMiddle,middleNode);

#undef _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_
#define _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xinit,x,v,spec,ptr,ParticleData,dt,startNode) \
    ElectricallyChargedDust::DustChargingProcessor_Implicit_SecondOrder(xinit,x,v,spec,ptr,ParticleData,dt,startNode);

#define _PIC_USER_DEFING_PARTICLE_SAMPLING_(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s) \
  ElectricallyChargedDust::Sampling::SampleParticleData(tempParticleData,LocalParticleWeight,tempSamplingBuffer,s);

//the model of the electrically charged dust
#undef _PIC_MODEL__DUST__MODE_
#define _PIC_MODEL__DUST__MODE_ _PIC_MODEL__DUST__MODE__ON_

#undef _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_
#define _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_

