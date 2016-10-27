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

  if ((eMeV<Emin)||(eMeV>Emax)) exit(__LINE__,__FILE__,"Error: out of the energy limit");

  logE=log10(eMeV/Emin);
  n=(int)(logE/dLogE);

  a=(log10(eMeV)-log10(Emin))/dLogE;
  n=(int)(a);
  a-=n;

  return (n!=nPoints-1) ? (1.0-a)*DiffFlux[n]+a*DiffFlux[n+1] : DiffFlux[n];
}


