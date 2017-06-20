//===================================================
//$Id$
//===================================================

#ifdef MPI_ON
#include "mpi.h"
#endif


#include <stdio.h>
#include <stdlib.h>

#include "hs.h"
#include "const.dfn"
#include "data.h"
#include "dsmc.dfn"
#include "specfunc.h"


//===================================================
void Chs::init(unsigned char total_species_number) {
  mass=new double[total_species_number];
  diam=new double[total_species_number];
  charge=new int[total_species_number];

  for (int s=0;s<total_species_number;s++) charge[s]=0;
} 

//===================================================
void Chs::SetCharge(int q,unsigned char s) {
  charge[s]=q;
}

//===================================================
int Chs::GetCharge(unsigned char s) {
  return charge[s]; 
}

//===================================================
void  Chs::SetRefDiam(double d, unsigned char s) {
  diam[s]=d;
}

//===================================================
void Chs::SetMass(double m, unsigned char s) {
  mass[s]=m;
}

//===================================================
double Chs::GetRefDiam(unsigned char s) {
  return diam[s];
}

//===================================================
double Chs::GetDiam(double Vrel,unsigned char s1,unsigned char s2) {
  return 0.5*(diam[s1]+diam[s2]);
}

//===================================================
double Chs::GetMass(unsigned char s) {
  return mass[s];
}

//===================================================
double Chs::GetTotalCrossSect(double Vrel,ParticlePtr ptr1,ParticlePtr ptr2) {
  double d;
  unsigned char s1,s2;

  s1=ParticleSolver.pbuffer.GetI(ptr1);
  s2=ParticleSolver.pbuffer.GetI(ptr2);
  d=0.5*(diam[s1]+diam[s2]);
  return Pi*d*d;
}

//===================================================
void Chs::GetPCVel(double* Vrel,unsigned char s1,unsigned char s2) {
  double Vrc,V[3];
  double CosKsi,SinKsi,CosEps,SinEps,D,c;

  CosKsi=2.0*rnd()-1.0;
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

//===================================================
double Chs::GetTotalCrossSect_Cr_MAX(unsigned char s1, 
  unsigned char s2, long int nsubcl,long int ncell,double CrossSection) {
  double res;
  
  ParticlePtr ptr;
  float v[3],vmax[3],vmin[3],gmax;
  int idim;
  unsigned char spec;
  bool initflag=false;

  ptr=ParticleSolver.grid.GetFirstPtr(nsubcl,ncell);
  while (ptr>=0) {
    spec=ParticleSolver.pbuffer.GetI(ptr); 

    if ((spec==s1)||(spec==s2)) { 
      ParticleSolver.pbuffer.GetV(v,ptr);

      if (initflag==false) for (initflag=true,idim=0;idim<3;idim++) vmax[idim]=v[idim],vmin[idim]=v[idim];  
      else for (idim=0;idim<3;idim++) { 
        if (v[idim]>vmax[idim]) vmax[idim]=v[idim];
        if (v[idim]<vmin[idim]) vmin[idim]=v[idim];
      }
    }

    ptr=ParticleSolver.pbuffer.GetNext(ptr);
  }

  gmax=(initflag==true) ? sqrt(pow(vmax[0]-vmin[0],2)+pow(vmax[1]-vmin[1],2)+pow(vmax[2]-vmin[2],2)) : 0.0; 

  res=(CrossSection>=0.0) ? gmax*CrossSection : gmax*pow(0.5*(diam[s1]+diam[s2]),2);  
  return res; 
 } 
