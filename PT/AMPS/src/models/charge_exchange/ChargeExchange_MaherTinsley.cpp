#include "ChargeExchange.h"

double ChargeExchange::MaherTinsley::CrossSection(int      spec,
						  double* vParticle,
						  double* vPlasma,
						  double   PlasmaTemperature) {
  // for notation see Heerikhuisen et al., JGR, Vol.111, A06110
  double v_th, v_rel, omega=0.0;

  // this mdel is for atomic hydrogen only
  if(spec != _H_SPEC_) return 0.0;
 
  v_th = sqrt(2.0 * PlasmaTemperature / _MASS_(_H_));
  for(int idim = 0; idim < 3; idim++ )
    omega += pow(vParticle[idim]-vPlasma[idim], 2.0);
  omega = pow(omega, 0.5) / v_th;
  v_rel = v_th * ( exp(-omega*omega)/sqrtPi + (omega + 0.5/omega)*erf(omega));

  return (1.6 - 0.0695 * pow(log(v_rel), 2.0))*1E-14*1E-4; // meters 
} 

#endif//_CHARGEEXCHANGE_H_
