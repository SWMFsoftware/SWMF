//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===================================================
//$Id$
//===================================================

#include <math.h>
#include "cell.h"
#include "specfunc.h"
#include "array_1d.h"
#include "data.h"

#include "const.dfn"

using namespace std;

extern int DIM;

double Ccell::measure_dim_0=1.0E-6;

//===================================================
void Ccell::SetVolume_dim0(float volume) {
  measure_dim_0=volume;
}

//===================================================
void Ccell::InitCellNormal() {
  int nface,nnode;
  float m;
  array_1d<float> r1(DIM),r2(DIM),n(DIM);

  for (nface=0;nface<DIM+1;nface++) {
    switch (DIM) {
    case 3 :
      r1=face[nface]->node[1]->X()-face[nface]->node[0]->X(); 
      r2=face[nface]->node[2]->X()-face[nface]->node[0]->X();
      n=cross_product(r1,r2);
      n.normalize();
      m=-1.0;
      for (nnode=0;nnode<DIM;nnode++) {
        r1=face[nface]->node[nnode]->X()-node[nface]->X();
        if (dot_product(r1,n)>0.0) m=1.0;
      }
      n*=m; 
      face[nface]->SetNormal(n);
      break;
    case 2 :
      r1=face[nface]->node[1]->X()-face[nface]->node[0]->X(); 
      n(0)=r1(1);n(1)=-r1(0);
      n.normalize();
      m=-1.0;
      for (nnode=0;nnode<DIM;nnode++) {
        r1=face[nface]->node[nnode]->X()-node[nface]->X();
        if (dot_product(r1,n)>0.0) m=1.0; 
      }
      n*=m;
      face[nface]->SetNormal(n);
      break;
    case 1 :
      n=node[1]->X()-node[0]->X();
      if (nface==1) n*=-1;
      n.normalize();
      face[nface]->SetNormal(n);
      break;
    default :
      printf("proc. Ccell::InitCellNormal()\n");
      printf("wrong DIM value: DIM=%i\n",DIM);
      exit(__LINE__,__FILE__);
    } 
  } 
}

//===================================================
double Ccell::Measure() {
  static array_1d<float> x_1d(1);
  static array_1d<float> x1_2d(2),x2_2d(2);
  static array_1d<float> x1_3d(3),x2_3d(3),x3_3d(3);

  if (measure_value>0.0) return measure_value;

  switch (DIM) {
  case 0 :
    measure_value=measure_dim_0;
    break;
  case 1 :
    if (SymmetryMode==no_symmetry) {
      x_1d=node[1]->X()-node[0]->X();
      measure_value=x_1d.abs();
    }
    else {
      double nd0[3],nd1[3];
      node[0]->GetX(nd0);node[1]->GetX(nd1);
      measure_value=4.0/3.0*Pi*fabs(pow(nd0[0],3)-pow(nd1[0],3));
    } 

    break;
  case 2 :
    x1_2d=node[1]->X()-node[0]->X();
    x2_2d=node[2]->X()-node[0]->X();
    measure_value=0.5*fabs(x1_2d(0)*x2_2d(1)-x1_2d(1)*x2_2d(0));

    if (SymmetryMode==cylindrical_symmetry) {
      double r=(node[0]->X()(1)+node[1]->X()(1)+node[2]->X()(1))/3.0;  
      measure_value*=2.0*Pi*r;
    } 
    break;
  case 3 :
    x1_3d=node[1]->X()-node[0]->X();
    x2_3d=node[2]->X()-node[0]->X();
    x3_3d=node[3]->X()-node[0]->X();
    measure_value=fabs(mix_product(x1_3d,x2_3d,x3_3d))/6.0;
    break;
  }

  return measure_value;
}

//===================================================
double Ccell::CharacteristicSize() {
  double nfaces,l,nd0[3],nd1[3],res=0.0;
  int i,j,idim;

  for (i=0;i<DIM+1;i++) {
    node[i]->GetX(nd0);  

    for (j=i+1;j<DIM+1;j++) { 
      node[j]->GetX(nd1);
      for (l=0.0,idim=0;idim<DIM;idim++) l+=pow(nd0[idim]-nd1[idim],2);

      res+=sqrt(l);
    }
  }

  nfaces=DIM*(1+DIM)/2.0;
  if (DIM!=0) res/=nfaces;

  return res;
}

//===================================================
void Ccell::RandomPosition(float* x) {
  float f1,f2,f3; 

  switch (DIM) {
  case 3 :
    f1=rnd(); if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999; f1=pow(f1,(float)0.3333333333);
    f2=rnd(); if (f2<0.00001) f2=0.00001; if (f2>0.99999) f2=0.99999; f2=sqrt(f2);
    f3=rnd(); if (f3<0.00001) f3=0.00001; if (f3>0.99999) f3=0.99999;
    x[0]=f1*(1.0-f2);
    x[1]=f1*f2*(1.0-f3);
    x[2]=f1*f2*f3;
    return;
  case 2 :
    f1=rnd(); if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999; f1=sqrt(f1);
    f2=rnd(); if (f2<0.00001) f2=0.00001; if (f2>0.99999) f2=0.99999;
    x[0]=f1*(1.0-f2);
    x[1]=f1*f2;
    return;
  case 1 : 
    f1=rnd(); if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999;
    x[0]=f1;
    return;
  case 0 :
    return;
  default :
    printf("proc. Cface::RandomPosition()\n");
    printf("wrong DIM value: DIM=%i\n",DIM);
    exit(__LINE__,__FILE__);
  }
}


//===================================================
bool CellsIntersection(Ccell& cell1,Ccell& cell2) {
  bool res; 
  int idim,pnode;
  float cells1_xmin[3],cells1_xmax[3];
  float cells2_xmin[3],cells2_xmax[3];
  array_1d<float> x(DIM);

//==============
  x=cell1.node[0]->X();
  for (idim=0;idim<DIM;idim++) {
    cells1_xmin[idim]=x(idim);
    cells1_xmax[idim]=x(idim);
  }
  
  for (pnode=1;pnode<DIM+1;pnode++) {
    x=cell1.node[pnode]->X();
    for (idim=0;idim<DIM;idim++) {
      if (cells1_xmin[idim]>x(idim)) cells1_xmin[idim]=x(idim);
      if (cells1_xmax[idim]<x(idim)) cells1_xmax[idim]=x(idim);
    }
  }

//==============
  x=cell2.node[0]->X();
  for (idim=0;idim<DIM;idim++) {
    cells2_xmin[idim]=x(idim);
    cells2_xmax[idim]=x(idim);
  }
  
  for (pnode=1;pnode<DIM+1;pnode++) {
    x=cell2.node[pnode]->X();
    for (idim=0;idim<DIM;idim++) {
      if (cells2_xmin[idim]>x(idim)) cells2_xmin[idim]=x(idim);
      if (cells2_xmax[idim]<x(idim)) cells2_xmax[idim]=x(idim);
    }
  }

//==============
  res=true;
  for (idim=0;idim<DIM;idim++)
    if ((cells1_xmax[idim]<cells2_xmin[idim])||(cells1_xmin[idim]>cells2_xmax[idim])) res=false;

  return res;
}
