//===================================================
//$Id$
//===================================================

#ifndef CELL
#define CELL

#include "face.h"
#include "dsmc.dfn"

class Ccell{
private:
  static double measure_dim_0;

public:
  long int nodeno[4],faceno[4],neighbour_cellno[4];
  Cnode* node[4];
  Cface* face[4];

  void SetVolume_dim0(float);
  void InitCellNormal();
  double Measure();
  void RandomPosition(float*); 
  friend bool CellsIntersection(Ccell&,Ccell&);

  Ccell() {
    for (int i=0;i<4;i++) nodeno[i]=-1,faceno[i]=-1,neighbour_cellno[i]=-1,node[i]=NULL,face[i]=NULL;
  }; 
};

#endif
