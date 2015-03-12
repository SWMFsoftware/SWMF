#include "ChargeExchange.h"


double ChargeExchange::CrossSection(int      spec,
				    double* vParticle,
				    double* vPlasma,
				    double   PlasmaTemperature) {
#if _CHARGE_EXCHANGE__MODEL_ == _CHARGE_EXCHANGE__MAHER_TINSLEY_
  return ChargeExchange::MaherTinsley::CrossSection(spec, vParticle, vPlasma, PlasmaTemperature);
#else
  exit(__LINE__,__FILE__,"Error: charge exchange model is not recognized");
#endif
}
