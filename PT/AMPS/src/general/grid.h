//==================================================
//$Id$
//==================================================

#ifndef GRID 
#define GRID

#include <vector>

#include "face.h"
#include "cell.h"
#include "node.h"
#include "array_4d.h"
#include "array_3d.h"

#include "dsmc.h"

struct cells_connection_data_type {
  int ne[4][4];
  int ne1[2][4][4],ne2[2][4][4];
  int  nbasis[4][4];
  int pface[4][4];
};

struct surface_interpolation_element_data_type {
  long int nnode;
  vector <long int> face;
  vector <double> weight;
};

struct surface_interpolation_data_type {
  long int nfaces,nnodes;
  vector<surface_interpolation_element_data_type> node;
};

class Cgrid{
private:

  void InitCellConnectionData();
public:
  long int nnodes,nfaces,ncells;

  Cnode* node;
  Cface* face;
  Ccell* cell;
  array_4d<float> TMatrix;
  cells_connection_data_type* cells_connection_data;
  array_4d<double> bvector;
  vector<surface_interpolation_data_type> surface_interpolation_data;

  Cgrid () {
    nnodes=0;node=NULL;
    nfaces=0;face=NULL;
    ncells=0;cell=NULL;
    cells_connection_data=NULL; 
  };

  ~Cgrid() {

     if (nnodes>0) {
       while (--nnodes>=0) node[nnodes].~Cnode();
       while (--nfaces>=0) face[nfaces].~Cface();
       while (--ncells>=0) cell[ncells].~Ccell();

       delete [] node;
       delete [] face;
       delete [] cell;

       delete [] cells_connection_data;
     }
  };

  double Measure();
  void LoadFromFile(char*);
  void InitGridData();
  void InitInterpolationData();
  void InitSurfaceInterpolationData(vector< vector<long int> >&);
  long int GetNCell(float*);
  void OutputInFile(char*);
  void SaveImageFile(int);
  void LoadImageFile(int);

  void GetTMatrix3D();
  void ChangeLocalVector3D(float*,char,char);
  void ChangeLocalPositionVector3D(float*,char,char);

  void GetTMatrix2D();
  void ChangeLocalVector2D(float*,char,char);
  void ChangeLocalPositionVector2D(float*,char,char); 

  void GetTMatrix1D();
  void ChangeLocalVector1D(float*,char,char);
  void ChangeLocalPositionVector1D(float*,char,char);
  }; 

#endif
