//$Id$

/*
 * DustScatteringEfficientcy.h
 *
 *  Created on: Jan 12, 2016
 *      Author: vtenishe
 */

#ifndef _DUSTSCATTERINGEFFICIENTCY_H_
#define _DUSTSCATTERINGEFFICIENTCY_H_


namespace LK {
  namespace Ice2Dust0_899999976__Porosity0_649122834 {
    const int nDataPoints=11;
    const double Data[nDataPoints][2]={
        {0.100000001, 0.0922335312},
        {0.199526235, 0.666842222},
        {0.398107201, 2.35393715},
        {0.794328272, 2.27319002},
        {1.58489335, 1.52125061},
        {3.16227818, 1.14077652},
        {6.3095746, 1.11582184 },
        {12.5892582, 1.10241592},
        {25.1188698, 1.09419072 },
        {50.1187286, 1.08739305 },
        {100., 1.08209419}};
  }

  namespace   Ice2Dust0_0500000007__Porosity0_828326166 {
    const int nDataPoints=11;
    const double Data[nDataPoints][2]={
      {0.100000001,  0.001184531721},
      {0.199526235,  0.00813579094},
      {0.398107201,  0.0401802771},
      {0.794328272,  0.161976323},
      {1.58489335,  0.546364427},
      {3.16227818,  1.35141349},
      {6.3095746,  1.55201817},
      {12.5892582,  1.08713472},
      {25.1188698,  1.02689314},
      {50.1187286,  1.01919293},
      {100.,  1.01795852}};
  }

  //the wave length of the scattaring efficientcy model
  const double ModelWavalength=0.65; //microns

  double GetScatteringEfficeintcy (double GrainRadius, const double ScatteringData[][2],int nScatteringDataPoints,double WaveLength);
}



#endif /* SRCCG_UTILITY_POSTPROCESSOR_DUSTSCATTERINGEFFICIENTCY_H_ */
