#ifndef NODE 
#define NODE

#include "array_1d.h"

extern int DIM;

class Cnode{
  float node_coordinate[3];
public:
  float* InterpolationWeight;
  long int* InterpolationMask;
  long int nodeno;

  Cnode() {
    InterpolationWeight=NULL;
    InterpolationMask=NULL;
  };

  ~Cnode () {
    if (InterpolationWeight!=NULL) delete [] InterpolationWeight;
    if (InterpolationMask!=NULL) delete [] InterpolationMask;
  };

  array_1d<float> X() {
    int idim;
    array_1d<float> x(DIM);

    for(idim=0;idim<DIM;idim++) x(idim)=node_coordinate[idim];
    return x;
  };

  void SetX(float* x) {
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
};

#endif
