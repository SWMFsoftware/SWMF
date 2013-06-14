//=======================================================================
//$Id$
//=======================================================================
//the file contains functions descpibing specific phhysical functions for sodium


#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "constants.h"
#include "He.h"


double RadiationPressureAcceleration__Chamberlian_TOPA_Eq7_2_9(double HeliocentricVelocity, double HeliocentricDistance) {
//TOPA Theory of Planetary Atmospheres

	double res=0.0, gfactor=0.0;
	static const double frequency__584_35A = 5.1303578E15;

	gfactor = HeliumGfactor__584_35A__Killen_2009_AJSS(HeliocentricVelocity, HeliocentricDistance);


	  //convert the acceleration from [cm s^{-2}] to [m s^{-2}] and scale it to the required heliocentric distance
	  res = PlanckConstant*frequency__584_35A*gfactor/1.00794/_AMU_/SpeedOfLight;
      return res;

}

double HeliumGfactor__584_35A__Killen_2009_AJSS(double HeliocentricVelocity,double HeliocentricDistance) {

  struct cGFactor {
    double RadiaVelocity,gFactor;
  };

  static const int nDataElements=49;
  static const double DataHeliocentricDistance=0.352*_AU_;

  static const cGFactor gFactor584_35A__Killen_2009_AJSS__Fig_3[nDataElements] = {
		  {-10.062411, 3.68739E-05},
		  {-9.605661,	3.73882E-05},
		  {-9.148911,	3.78862E-05},
		  {-8.692161,	3.83962E-05},
		  {-8.235411,	3.88327E-05},
		  {-7.778661,	3.93587E-05},
		  {-7.321911,	3.98157E-05},
		  {-6.865161,	4.02665E-05},
		  {-6.408411,	4.07376E-05},
		  {-5.951661,	4.11913E-05},
		  {-5.494911,	4.16437E-05},
		  {-5.038161,	4.20787E-05},
		  {-4.581411,	4.2548E-05},
		  {-4.124661,	4.30009E-05},
		  {-3.667911,	4.33886E-05},
		  {-3.211161,	4.38068E-05},
		  {-2.754411,	4.42214E-05},
		  {-2.297661,	4.46334E-05},
		  {-1.840911,	4.49946E-05},
		  {-1.384161,	4.5378E-05},
		  {-0.927411,	4.57672E-05},
		  {-0.470661,	4.60842E-05},
		  {-0.013911,	4.64032E-05},
		  {0.442839,	4.67652E-05},
		  {0.899589,	4.70868E-05},
		  {1.356339,	4.73678E-05},
		  {1.813089,	4.76521E-05},
		  {2.269839,	4.79448E-05},
		  {2.726589,	4.82298E-05},
		  {3.183339,	4.84553E-05},
		  {3.640089,	4.8691E-05},
		  {4.096839,	4.89079E-05},
		  {4.553589,	4.91402E-05},
		  {5.010339,	4.93365E-05},
		  {5.467089,	4.94684E-05},
		  {5.923839,	4.96457E-05},
		  {6.380589,	4.98187E-05},
		  {6.837339,	4.99137E-05},
		  {7.294089,	5.00432E-05},
		  {7.750839,	5.01827E-05},
		  {8.207589,	5.02665E-05},
		  {8.664339,	5.03102E-05},
		  {9.121089,	5.03847E-05},
		  {9.577839,	5.04431E-05},
		  {10.034589,	5.04304E-05},
		  {10.491339,	5.04292E-05},
		  {10.948089,	5.0428E-05},
		  {11.404839,	5.03947E-05},
		  {11.861589,	5.03423E-05}};

  static const double dv=(gFactor584_35A__Killen_2009_AJSS__Fig_3[1].RadiaVelocity-gFactor584_35A__Killen_2009_AJSS__Fig_3[0].RadiaVelocity)*1.0E3;
  static const double vmin=gFactor584_35A__Killen_2009_AJSS__Fig_3[0].RadiaVelocity*1.0E3;

  int Interval=(int)((HeliocentricVelocity-vmin)/dv);
  double res;

  if (Interval<0) res=gFactor584_35A__Killen_2009_AJSS__Fig_3[0].gFactor;
  else if (Interval>=nDataElements) res=gFactor584_35A__Killen_2009_AJSS__Fig_3[nDataElements-1].gFactor;
  else res=gFactor584_35A__Killen_2009_AJSS__Fig_3[Interval].gFactor;

  res*=pow(DataHeliocentricDistance/HeliocentricDistance,2);

  return res;
}


