//===================================================
//$Id$
//===================================================

#ifdef MPI_ON
#include "mpi.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "vss.h"
#include "const.dfn"
#include "data.h"
#include "dsmc.dfn"
#include "specfunc.h"

//===================================================
void Cvss::init(unsigned char total_species_number) {
  alfa.resize(total_species_number);
  for (int s1=0;s1<NS;s1++) {
    alfa[s1].resize(total_species_number);
    for (int s2=0;s2<NS;s2++) alfa[s1][s2]=0.0;
  }

  TempIndex.resize(total_species_number);
  for (int s1=0;s1<NS;s1++) {
    TempIndex[s1].resize(total_species_number);
    for (int s2=0;s2<NS;s2++) TempIndex[s1][s2]=0.0;
  }

  Tref=new double[total_species_number];
  mass=new double[total_species_number]; 
  diam=new double[total_species_number];
} 

//===================================================
void  Cvss::SetAlfa(double x, unsigned char s1,unsigned char s2) {
  alfa[s1][s2]=x;
  alfa[s2][s1]=x;
}

//===================================================
double Cvss::GetAlfa(unsigned char s1,unsigned char s2) {
  return alfa[s1][s2];
}

//===================================================
void Cvss::GetPCVel(double* Vrel,unsigned char s1,unsigned char s2) {
  double Vrc,V[3];
  double CosKsi,SinKsi,CosEps,SinEps,D,c;

  CosKsi=2.0*pow(rnd(),1.0/alfa[s1][s2])-1.0;
  SinKsi=sqrt(1.0-CosKsi*CosKsi);

  c=2*Pi*rnd();
  SinEps=sin(c);
  CosEps=cos(c);

  D=sqrt(Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);
  if (D>1.0E-6) {
    Vrc=sqrt(Vrel[0]*Vrel[0]+Vrel[1]*Vrel[1]+Vrel[2]*Vrel[2]);
    V[0]=CosKsi*Vrel[0]+SinKsi*SinEps*D;
    V[1]=CosKsi*Vrel[1]+
      SinKsi*(Vrc*Vrel[2]*CosEps-Vrel[0]*Vrel[1]*SinEps)/D;
    V[2]=CosKsi*Vrel[2]-
      SinKsi*(Vrc*Vrel[1]*CosEps+Vrel[0]*Vrel[2]*SinEps)/D;}
  else {
    V[0]=CosKsi*Vrel[0];
    V[1]=SinKsi*CosEps*Vrel[0];
    V[2]=SinKsi*SinEps*Vrel[0];
  }

  Vrel[0]=V[0];
  Vrel[1]=V[1];
  Vrel[2]=V[2];
}
