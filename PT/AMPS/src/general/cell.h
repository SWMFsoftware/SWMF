//===================================================
//$Id$
//===================================================

#ifndef CELL
#define CELL

#include <iostream.h>
#include "face.h"
#include "dsmc.dfn"

class Ccell{
private:
  static double measure_dim_0;

public:
  ParticlePtr first_ptr;
  int thread,sbdm;
  long int nodeno[4],faceno[4],neighbour_cellno[4];
  Cnode* node[4];
  Cface* face[4];

  void SetVolume_dim0(float);
  void InitCellNormal();
  double Measure();
  void RandomPosition(float*); 
  friend bool CellsIntersection(Ccell&,Ccell&);
};

#endif
