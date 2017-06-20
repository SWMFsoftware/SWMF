//===================================================
//$Id$
//===================================================


#ifdef MPI_ON
#include "mpi.h"
#endif


#include <stdlib.h>
#include <vector>

#include "vhs.h"
#include "const.dfn"
#include "data.h"
#include "dsmc.dfn"
#include "specfunc.h"

//===================================================
void Cvhs::init(unsigned char total_species_number) {
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
void  Cvhs::SetTempIndex(double x, unsigned char s1,unsigned char s2) {
  TempIndex[s1][s2]=x;
  TempIndex[s2][s1]=x;
}

//===================================================
void Cvhs::SetRefTemp(double t, unsigned char s) {
  Tref[s]=t;
}

//===================================================
//===================================================
double Cvhs::GetDiam(double rel_speed,unsigned char s1,unsigned char s2) {
  double d,T,w,m;

  if (rel_speed<1.0E-5) return 0.0;

  d=0.5*(diam[s1]+diam[s2]);
  T=0.5*(Tref[s1]+Tref[s2]);
  w=TempIndex[s1][s2];
  m=mass[s1]*mass[s2]/(mass[s1]+mass[s2]);

//Bird. Eq 4.63
  return d*sqrt(pow(2*Kbol*T/(m*rel_speed*rel_speed),w-0.5)/gamma(2.5-w));
}

//===================================================
double Cvhs::GetTempIndex(unsigned char s1,unsigned char s2) {
  return TempIndex[s1][s2];
}

//===================================================
double Cvhs::GetRefTemp(unsigned char s) {
  return Tref[s];
}
//===================================================
//===================================================
double Cvhs::GetTotalCrossSect(double Vrel,ParticlePtr ptr1,ParticlePtr ptr2) {
  double d;
  unsigned char s1,s2;

  s1=ParticleSolver.pbuffer.GetI(ptr1);
  s2=ParticleSolver.pbuffer.GetI(ptr2);

  d=GetDiam(Vrel,s1,s2);
  return Pi*d*d;
}

//===================================================
double Cvhs::GetTotalCrossSect_Cr_MAX(unsigned char s1, 
  unsigned char s2, long int nsubcl,long int ncell,ECSFunc ExternalCrossSectFunc) {
  
  unsigned char s;
  ParticlePtr ptr;
  long int npart_s1,npart_s2;
  ParticlePtr *part_buff_s1,*part_buff_s2; 
  ParticlePtr part_pos_s1,part_pos_s2;

  npart_s1=0;
  npart_s2=0;
  ptr=ParticleSolver.grid.GetFirstPtr(nsubcl,ncell);
  while (ptr>=0) {
    s=ParticleSolver.pbuffer.GetI(ptr);
    if (s==s1) npart_s1++;
    if (s==s2) npart_s2++; 
    ptr=ParticleSolver.pbuffer.GetNext(ptr);
  }

  if ((npart_s1==0)||(npart_s2==0)) return 0.0; 
  part_buff_s1=new ParticlePtr[npart_s1];
  part_buff_s2=new ParticlePtr[npart_s2];

  part_pos_s1=0;
  part_pos_s2=0;

  ptr=ParticleSolver.grid.GetFirstPtr(nsubcl,ncell);
  while (ptr>=0) {
    s=ParticleSolver.pbuffer.GetI(ptr);
    if (s==s1) {
      part_buff_s1[part_pos_s1]=ptr;
      part_pos_s1++;
    }
    if (s==s2) {
      part_buff_s2[part_pos_s2]=ptr;
      part_pos_s2++;
    }
    ptr=ParticleSolver.pbuffer.GetNext(ptr);
  }

  double ncombinations;
  ParticlePtr i_s1,i_s2,ptr1,ptr2;
  int idim;
  float v1[3],v2[3];
  double Vrel[3],Vrc,dot_pr,res,c;
  long int n_try;

  res=0.0;
  ncombinations=npart_s1*npart_s2;
  if (ncombinations<5000.0) {
    for (i_s1=0;i_s1<npart_s1;i_s1++)
      for (i_s2=0;i_s2<npart_s2;i_s2++) {
        if ((s1==s2)&&(i_s1==i_s2)) continue;
        ParticleSolver.pbuffer.GetV(v1,part_buff_s1[i_s1]);
        ParticleSolver.pbuffer.GetV(v2,part_buff_s2[i_s2]); 

        dot_pr=0.0;
        for (idim=0;idim<3;idim++) {
          Vrel[idim]=v1[idim]-v2[idim];
          dot_pr+=Vrel[idim]*Vrel[idim];
        }
        Vrc=sqrt(dot_pr);
        c=Vrc*((ExternalCrossSectFunc==NULL) ? GetTotalCrossSect(Vrc,part_buff_s1[i_s1],part_buff_s2[i_s2]) : (*ExternalCrossSectFunc)(Vrc,part_buff_s1[i_s1],part_buff_s2[i_s2],ncell));
        if (c>res) res=c; 
      } 
  }
  else {
    n_try=5000;
    while (n_try-->0) {
      do {
        ptr1=part_buff_s1[(ParticlePtr)(rnd()*npart_s1)]; 
        ptr2=part_buff_s2[(ParticlePtr)(rnd()*npart_s2)]; 
      }
      while ((s1==s2)&&(ptr1==ptr2));

      ParticleSolver.pbuffer.GetV(v1,ptr1);
      ParticleSolver.pbuffer.GetV(v2,ptr2);

      dot_pr=0.0;
      for (idim=0;idim<3;idim++) {
        Vrel[idim]=v1[idim]-v2[idim];
        dot_pr+=Vrel[idim]*Vrel[idim];
      }
      Vrc=sqrt(dot_pr);
      c=Vrc*((ExternalCrossSectFunc==NULL) ? GetTotalCrossSect(Vrc,ptr1,ptr2) : (*ExternalCrossSectFunc)(Vrc,ptr1,ptr2,ncell));
      if (c>res) res=c;
    }
  }

  delete [] part_buff_s1;
  delete [] part_buff_s2;
  return res;
}
