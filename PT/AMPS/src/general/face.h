#ifndef FACE
#define FACE

#include "node.h"
#include <iostream.h>
#include "array_1d.h"
#include "specfunc.h"

class Cface{
  float nrml[3];
public:
  unsigned char faceat,surface_group;
  long int nodeno[3];
  Cnode* node[3];

//===================================================
  void SetNormal(array_1d<float>& n) {
    int idim;

    for (idim=0;idim<DIM;idim++) nrml[idim]=n(idim);
  };
//===================================================
  array_1d<float> Normal() {
    array_1d<float> normal(DIM);
    int idim;

    for (idim=0;idim<DIM;idim++) normal(idim)=nrml[idim];
    return normal;
  };

//===================================================
  void GetNormal(float* n) {
    for (int idim=0;idim<DIM;idim++) n[idim]=nrml[idim];
  }; 
//===================================================
  double Measure() {
    static double measure;
    static array_1d<float> x_2d(2);
    static array_1d<float> x1_3d(3),x2_3d(3);

    switch(DIM) {
    case 1:  
      measure=1.0;
      break;
    case 2:
      x_2d=(*node[1]).X()-(*node[0]).X();
      measure=x_2d.abs();
      break;
    case 3:
      x1_3d=(*node[1]).X()-(*node[0]).X();
      x2_3d=(*node[2]).X()-(*node[0]).X();
      measure=0.5*fabs(cross_product(x1_3d,x2_3d).abs());
      break;
    default : 
      cout << "proc. Cface::Measure()" << endl; 
      cout << "wrong DIM value: DIM=" << DIM << endl;
      exit(0);
    }

    return measure;
  };
//===================================================
  void RandomPosition(float* x) {
    float f1,f2;

    switch (DIM) {
    case 3 :
      f1=rnd(); if (f1<0.001) f1=0.001; if (f1>0.999) f1=0.999; f1=sqrt(f1);
      f2=rnd(); if (f2<0.001) f2=0.001; if (f2>0.999) f2=0.999;
      x[0]=f1*(1.0-f2);
      x[1]=f1*f2;
      return;
    case 2 :
      f1=rnd(); if (f1<0.001) f1=0.001; if (f1>0.999) f1=0.999;
      x[0]=f1;
      return;
    default :
      cout << "proc. Cface::RandomPosition()" << endl;
      cout << "wrong DIM value: DIM=" << DIM << endl;
      exit(0);
    }

  };
//===================================================
};
#endif
