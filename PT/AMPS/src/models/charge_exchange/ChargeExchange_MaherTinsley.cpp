#include "ChargeExchange.h"

double ChargeExchange::MaherTinsley::LifeTime(int      spec,
					      double* vParticle,
					      double* vPlasma,
					      double   PlasmaTemperature,
					      double   PlasmaNumberDensity) {
  // for notation see Heerikhuisen et al., JGR, Vol.111, A06110
  double v_th, v_rel, omega=0.0, sigma;

  // this model is for atomic hydrogen only
  if(spec != _H_SPEC_) return 1E+100;
 
  v_th = sqrt(2.0 * Kbol * PlasmaTemperature / _MASS_(_H_));
  for(int idim = 0; idim < 3; idim++ )
    omega += pow(vParticle[idim]-vPlasma[idim], 2.0);
  omega = pow(omega, 0.5) / v_th;
  if(omega > 1E-2)
    v_rel = v_th* ( exp(-omega*omega)/sqrtPi + (omega + 0.5/omega)*erf(omega));
  else
    // Taylor expansion in vicinity of omega = 0 to avoid divison by zero;
    // relative error < 1E-4
    v_rel = v_th* 2.0 * (1.0 + omega*omega/3.0) / sqrtPi;

  sigma = pow(1.6 - 0.0695 * log(v_rel), 2.0)*1E-14*1E-4; // m^2
  return 1.0 / (PlasmaNumberDensity * v_rel * sigma);
} 

