//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//==================================================
//$Id$
//==================================================

#ifndef NODE 
#define NODE

#include "array_1d.h"

extern int DIM;

template <class T=float>
class Cnode{
  T node_coordinate[3];
public:
  T* InterpolationWeight;
  long int* InterpolationMask;
  long int nodeno;

  Cnode() {
    InterpolationWeight=NULL;
    InterpolationMask=NULL;
  };

/*
  ~Cnode () {
    if (InterpolationWeight!=NULL) delete [] InterpolationWeight;
    if (InterpolationMask!=NULL) delete [] InterpolationMask;
  };
*/

  array_1d<T> X() {
    int idim;
    array_1d<T> x(DIM);

    for(idim=0;idim<DIM;idim++) x(idim)=node_coordinate[idim];
    return x;
  };

  void SetX(T* x) {
    int idim;
    for(idim=0;idim<DIM;idim++) node_coordinate[idim]=x[idim];
  };

  void GetX(float* x) {
    int idim;
    for(idim=0;idim<DIM;idim++) x[idim]=node_coordinate[idim];
  };

  void GetX(double* x) {
    int idim;
    for(idim=0;idim<DIM;idim++) x[idim]=node_coordinate[idim];
  };

  T* GetX() {
    return node_coordinate;
  };
};

#endif
