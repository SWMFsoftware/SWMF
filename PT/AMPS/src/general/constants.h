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


//physical constants
#define Kbol                 1.3806503E-23
#define ElectronCharge       1.602176565E-19
#define ElectronMass         9.10938291E-31
#define ProtonMass           1.67262158E-27
#define eV2J                 ElectronCharge
#define GravityConstant      6.67300E-11
#define VacuumPermittivity   8.854187817620E-12



//macro funcrions for generating planet spacifica calls
#define _MASS1_(x) x##_MASS_
#define _MASS_(x) _MASS1_(x)

#define _RADIUS1_(x) x##_RADIUS_
#define _RADIUS_(x) _RADIUS1_(x)

//get the string value of a macro
#define _MACRO_STR_VALUE1_(arg) #arg
#define _MACRO_STR_VALUE_(arg) _MACRO_STR_VALUE1_(arg)

#endif /* CONSTANTS_H_ */
