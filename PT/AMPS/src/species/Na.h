//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the file contains functions descpibing specific phhysical functions for sodium

#ifndef _Na_PHYSICAL_PARAMETERS_
#define _Na_PHYSICAL_PARAMETERS_

//the solar dariation pressure as a function of heliocentric distrance and velocty
//taken from Combi-1997-icarus.pdf

//input parameters: heliocentric velociy (m/s) and heliocentric distance (m)
//return the radiation pressure (m/s^2)
//the positive velocity is in the direction out from the sun
//the acceleration is in teh direction out of the sun

double SodiumRadiationPressureAcceleration__Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
double SodiumGfactor__D1D2__Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);

double SodiumGfactor__5891_58A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance);
double SodiumGfactor__5897_56A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance);

#endif

