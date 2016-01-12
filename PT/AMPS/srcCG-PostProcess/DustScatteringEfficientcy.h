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

  double GetScatteringEfficeintcy (double GrainRadius, const double ScatteringData[][2],int nScatteringDataPoints);
}



#endif /* SRCCG_UTILITY_POSTPROCESSOR_DUSTSCATTERINGEFFICIENTCY_H_ */
