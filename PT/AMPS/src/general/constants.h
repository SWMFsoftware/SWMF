//===================================================
//$Id$
//===================================================
//The file containes constants that can be used in the model

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

//geometrical constants
#define Pi       3.14159265358979323846264338327950288419716939937510582
#define sqrtPi   1.7724538509055160272981674833411
#define PiTimes2 6.28318530717958647692528676655900576839433879875021

//physical constants
#define Kbol              1.3806503E-23
#define ElectronCharge    1.602176565E-19
#define ElectronMass      9.10938291E-31
#define ProtonMass        1.67262158E-27
#define eV2J              ElectronCharge
#define GravityConstant   6.67300E-11

//macro funcrions for generating planet spacifica calls
#define _MASS1_(x) x##_MASS_
#define _MASS_(x) _MASS1_(x)

#define _RADIUS1_(x) x##_RADIUS_
#define _RADIUS_(x) _RADIUS1_(x)

//Planetary objects
//MARS
//#define MarsMass   6.4185E23
#define _MARS__RADIUS_ 3.396E6

#define _MARS__MASS_ 6.4185E23

#endif /* CONSTANTS_H_ */
