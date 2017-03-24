//$Id$
#include "ChargeExchange.h"

double ChargeExchange::LindsayStebbings::LifeTime(int      spec,
						  double* vParticle,
						  double* vPlasma,
						  double   PlasmaTemperature,
						  double   PlasmaNumberDensity) {
  // for notation see Lindsay & Stebbings, JGR, vol. 110, A12213, 2005
  double a1 = 4.15, a2 = 0.531, a3 = 67.3, v, v2, energy, omega=0.0, sigma;

  // this model is for atomic hydrogen and three charge exchange species only
  if(spec != _H_SPEC_ && spec != _H_ENA_V1_SPEC_ && spec != _H_ENA_V2_SPEC_ && spec != _H_ENA_V3_SPEC_) return 1E+100;
 
  for(int idim = 0; idim < 3; idim++ )
    v2 += pow(vParticle[idim]-vPlasma[idim], 2.0);
  if(v2 < 1E-8) return 1E+20; // to avoid division by zero
  energy = 0.5 * _MASS_(_H_) * v2 /eV2J / 1000; // energy in keV
  double c;
  if (energy > 4.0)
    c = pow(1 - exp( -a3/energy ), 4.5);
  else
    c = 1.0; // error < 1E-6
  if(energy > 1E-10)
    sigma = pow(a1 - a2*log(energy), 2.0) * c * 1E-20; // m^2
  else
    // to avoid computing log() in vicinity of energy = 0: 
    // substitute log(energy) with (-10*log(10))
    sigma = pow(a1 + a2*23.026, 2.0) * 1E-20; // m^2
  return 1.0 / (PlasmaNumberDensity * pow(v2, 0.5) * sigma);
} 

