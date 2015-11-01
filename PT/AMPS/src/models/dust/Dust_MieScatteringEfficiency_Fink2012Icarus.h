//$Id$

/*
 * Dust_MieScatteringEfficiency_Fink2012Icarus.h
 *
 *  Created on: Oct 31, 2015
 *      Author: vtenishe
 */

#ifndef _SRC_MODELS_DUST_DUST_SCATTERINGEFFICIENCY_FINK2012ICARUS_H_
#define _SRC_MODELS_DUST_DUST_SCATTERINGEFFICIENCY_FINK2012ICARUS_H_

//dust grain Mie scattering effcientcy (Extracted from Fink-2012-icarus)
namespace Fink2012Icarus {
  //available cases
  const int n20i0001=0;   //the case for n=2.0 and i=0.001
  const int n20i01=1;     //the case for n=2.0 and i=0.10
  const int n20i04=2;     //the case for n=2.0 and i=0.40

  const int nCases=3;     //the total number of the digitized curves

  const int nPoints=51;
  const double aMinLog10=-1.0;
  const double aMaxLog10=3.0;
  const double dLog10=(aMaxLog10-aMinLog10)/(nPoints-1);

  const double DataArray[51][3]={
      {0.013, 0.007, 0.007},
      {0.013, 0.007, 0.007},
      {0.013, 0.007, 0.007},
      {0.013, 0.007, 0.013},
      {0.006, 0.007, 0.013},
      {0.019, 0.019, 0.013},
      {0.026, 0.026, 0.019},
      {0.026, 0.026, 0.026},
      {0.045, 0.052, 0.052},
      {0.096, 0.097, 0.116},
      {0.180, 0.187, 0.199},
      {0.495, 0.450, 0.380},
      {0.990, 0.836, 0.617},
      {1.543, 1.235, 0.900},
      {2.347, 1.762, 1.177},
      {3.190, 2.296, 1.402},
      {3.614, 2.534, 1.524},
      {3.543, 2.450, 1.518},
      {3.370, 2.225, 1.460},
      {3.164, 1.987, 1.415},
      {2.907, 1.724, 1.351},
      {2.695, 1.505, 1.299},
      {2.553, 1.402, 1.286},
      {2.457, 1.383, 1.280},
      {2.392, 1.325, 1.273},
      {2.367, 1.306, 1.267},
      {2.341, 1.267, 1.261},
      {2.289, 1.209, 1.254},
      {2.244, 1.196, 1.254},
      {2.180, 1.190, 1.248},
      {2.096, 1.183, 1.241},
      {2.032, 1.183, 1.241},
      {1.994, 1.171, 1.235},
      {1.981, 1.171, 1.228},
      {1.955, 1.171, 1.222},
      {1.910, 1.164, 1.222},
      {1.871, 1.164, 1.216},
      {1.814, 1.164, 1.222},
      {1.756, 1.164, 1.209},
      {1.704, 1.151, 1.209},
      {1.646, 1.151, 1.209},
      {1.595, 1.138, 1.209},
      {1.543, 1.151, 1.203},
      {1.492, 1.145, 1.203},
      {1.428, 1.145, 1.196},
      {1.383, 1.145, 1.196},
      {1.331, 1.138, 1.196},
      {1.293, 1.138, 1.196},
      {1.254, 1.132, 1.196},
      {1.235, 1.126, 1.196},
      {1.215, 1.126, 1.196}
  };

  inline double Get(double GrainRadius,double WaweLength,int Case) {
    double aLog;

    aLog=log10(2.0*Pi*GrainRadius/WaweLength);

    if (GrainRadius<aMinLog10) return DataArray[0][Case];
    else if (GrainRadius>aMaxLog10) return DataArray[nPoints-1][Case];
    else {
      double w,x;
      int i;

      x=(aLog-aMinLog10)/dLog10;
      i=(int)x;

      w=x-i;
      return w*DataArray[i+1][Case]+(1.0-w)*DataArray[i][Case];
    }
  }


}



#endif /* _SRC_MODELS_DUST_DUST_SCATTERINGEFFICIENCY_FINK2012ICARUS_H_ */
