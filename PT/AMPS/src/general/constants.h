//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================
//The file containes constants that can be used in the model

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "constants.AtomicData.h"
#include "constants.PlanetaryData.h"

//geometrical constants
#define Pi       3.14159265358979323846264338327950288419716939937510582
#define sqrtPi   1.7724538509055160272981674833411
#define PiTimes2 6.28318530717958647692528676655900576839433879875021

//length conversion
#define _MICROMETER_ 1.0E-6

//metric conversion coefficients
#define _NANO_  1.0E-9
#define _MICRO_ 1.0E-6
#define _A_     1.0E-10


//physical constants
#define Kbol                 1.3806503E-23
#define ElectronCharge       1.602176565E-19
#define ElectronMass         9.10938291E-31
#define ProtonMass           1.67262158E-27
#define eV2J                 ElectronCharge
#define EV2J                 ElectronCharge
#define GravityConstant      6.67300E-11
#define VacuumPermittivity   8.854187817620E-12
#define PlanckConstant       6.62606957E-34
#define SpeedOfLight         299792458.0
#define cm2m                 1.0E-2
#define m2cm                 100.0
#define J2eV                 (1.0/eV2J)
#define KeV2J                (1.0E3*eV2J)
#define MeV2J                (1.0E6*eV2J)
#define J2KeV                (1.0/KeV2J)
#define J2MeV                (1.0/MeV2J)

//convert time units
#define Min2Sec 60.0
#define Hour2Sec 60.0*Min2Sec


//macro funcrions for generating planet spacifica calls
#define _MASS1_(x) x##_MASS_
#define _MASS_(x) _MASS1_(x)

#define _RADIUS1_(x) x##_RADIUS_
#define _RADIUS_(x) _RADIUS1_(x)

#define _TARGET_ID1_(x) x##_ID_
#define _TARGET_ID_(x) _TARGET_ID1_(x)

//get the string value of a macro
#define _MACRO_STR_VALUE1_(arg) #arg
#define _MACRO_STR_VALUE_(arg) _MACRO_STR_VALUE1_(arg)

#endif /* CONSTANTS_H_ */
