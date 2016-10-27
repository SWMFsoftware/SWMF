//$Id$
//Energy speactra of Hydrogen GCR extracted from Badavi-2011-ASR (Fig. 5)

#include "specfunc.h"
#include "constants.h"

namespace GCR_BADAVI2011ASR {

  namespace Hydrogen {
    const double Emin=1.0; //MeV
    const double Emax=1.0E5; //Mev
    const int nPoints=51; //uniformly spaced in the logarithmic energy space
    extern double dLogE;

    const double DiffFlux[nPoints]={1.414,1.802,2.295,2.875,3.642,4.491,5.690,7.130,9.159,11.156,14.310,18.146,22.814,28.319,35.434,44.186,55.525,70.345,87.384,108.134,
        129.273,156.858,184.130,211.829,238.265,251.075,256.797,249.766,228.098,206.029,172.182,136.890,103.534,74.702,51.625,34.966,22.509,
        14.173,8.807,5.173,3.145,1.782,1.056,0.592,0.325,0.174,0.100,0.055,0.030,0.016,0.009};

    double GetDiffFlux(double E); 
  }
}



