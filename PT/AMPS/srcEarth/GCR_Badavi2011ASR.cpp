//$Id$
//Energy speactra of Hydrogen GCR extracted from Badavi-2011-ASR (Fig. 5)

#include "specfunc.h"
#include "constants.h"
#include "GCR_Badavi2011ASR.h"

double GCR_BADAVI2011ASR::Hydrogen::dLogE=
  (log10(GCR_BADAVI2011ASR::Hydrogen::Emax)-log10(GCR_BADAVI2011ASR::Hydrogen::Emin))/(GCR_BADAVI2011ASR::Hydrogen::nPoints-1);

double GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(double E) {
  double logE,a,eMeV=E*J2MeV;
  int n;

  const int EPS=1.0E-4;

  if ((eMeV<Emin)||(eMeV>Emax)) {
    char msg[200];

    if (eMeV>Emin*(1.0-EPS)) eMeV=Emin;
    else if (eMeV<Emax*(1.0+EPS)) eMeV=Emax; 
    else {
      sprintf(msg,"Error: out of the energy limit (E=%e MeV)",eMeV);
      exit(__LINE__,__FILE__,msg);
     }
  }

  logE=log10(eMeV/Emin);
  n=(int)(logE/dLogE);

  a=(log10(eMeV)-log10(Emin))/dLogE;
  n=(int)(a);
  a-=n;

  return (n!=nPoints-1) ? (1.0-a)*DiffFlux[n]+a*DiffFlux[n+1] : DiffFlux[n];
}


