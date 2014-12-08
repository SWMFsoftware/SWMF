//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the file contains functions descpibing specific phhysical functions for sodium

#ifndef _He_PHYSICAL_PARAMETERS_
#define _He_PHYSICAL_PARAMETERS_

//the solar dariation pressure as a function of heliocentric distrance and velocty
//taken from Combi-1997-icarus.pdf

//input parameters: heliocentric velociy (m/s) and heliocentric distance (m)
//return the radiation pressure (m/s^2)
//the positive velocity is in the direction out from the sun
//the acceleration is in teh direction out of the sun

double RadiationPressureAcceleration__Chamberlian_TOPA_Eq7_2_9(double HeliocentricVelocity, double HeliocentricDistance);
double HeliumGfactor__584_35A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance);

#endif
