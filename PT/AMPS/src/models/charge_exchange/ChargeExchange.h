//$Id$
// Physical model for charge exchange
// HOW TO ADD A NEW MODEL
// Two additions are required:
// 1) add a new namespace(s):
//    - namespace ModelName  
//    source file Sputtering_<ModelName>.cpp contains yield table / functions
// 2) add a corresponding branch in:
//    - ChargeExchange::Init,
//    - ChargeExchange::GetCrossSection,
//------------------------------------------------------------------------

#ifndef _CHARGE_EXCHANGE_H_ 
#define _CHARGE_EXCHANGE_H_

#include "pic.h"
#include "ChargeExchange.dfn"


namespace ChargeExchange {

  // Maher & Tinsley model (JGR, vol.82, 689-695)==============================
  namespace MaherTinsley {

    double LifeTime(int      spec,
		    double* vParticle,
		    double* vPlasma,
		    double   PlasmaTemperature,
		    double   PlasmaNumberDensity);
  } 
  // namespace MaherTinsley====================================================

  // Lindsay & Stebbings model (JGR, vol.110, A12213, 2005)====================
  namespace LindsayStebbings {

    double LifeTime(int      spec,
		    double* vParticle,
		    double* vPlasma,
		    double   PlasmaTemperature,
		    double   PlasmaNumberDensity);
  } 
  // namespace LindsayStabbings================================================


  // wrapper function: 
  double LifeTime(int      spec,
		  double* vParticle,
		  double* vPlasma,
		  double   PlasmaTemperature,
		  double   PlasmaNumberDensity);
}

#endif//_CHARGEEXCHANGE_H_
