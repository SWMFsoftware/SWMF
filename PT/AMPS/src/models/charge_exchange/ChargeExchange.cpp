//$Id$
#include "ChargeExchange.h"


double ChargeExchange::LifeTime(int      spec,
				double* vParticle,
				double* vPlasma,
				double   PlasmaTemperature,
				double   PlasmaNumberDensity) {
#if _CHARGE_EXCHANGE__MODEL_ == _CHARGE_EXCHANGE__MAHER_TINSLEY_
  return ChargeExchange::MaherTinsley::LifeTime(spec, vParticle, vPlasma, PlasmaTemperature,PlasmaNumberDensity);
#elif _CHARGE_EXCHANGE__MODEL_ == _CHARGE_EXCHANGE__LINDSAY_STEBBINGS_
  return ChargeExchange::LindsayStebbings::LifeTime(spec, vParticle, vPlasma, PlasmaTemperature,PlasmaNumberDensity);
#else
  exit(__LINE__,__FILE__,"Error: charge exchange model is not recognized");
#endif
}
