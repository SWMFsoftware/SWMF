//=======================================================================
//$Id$
//=======================================================================
//the file contains functions descpibing specific phhysical functions for sodium

#ifndef _Na_PAHYSICAL_PARAMETERS_
#define _Na_PAHYSICAL_PARAMETERS_

//the solar dariation pressure as a function of heliocentric distrance and velocty
//taken from Combi-1997-icarus.pdf

//input parameters: heliocentric velociy (m/s) and heliocentric distance (m)
//return the radiation pressure (m/s^2)
//the positive velocity is in the direction out from the sun
//the acceleration is in teh direction out of the sun

double SodiumRadiationPressureAcceleration__Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance); 

#endif

