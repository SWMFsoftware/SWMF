//===================================================
//$Id$
//===================================================

#ifndef CELL
#define CELL

#include <math.h>
#include "specfunc.h"
#include "array_1d.h"
#include "node.h"
#include "face.h"

#include "const.dfn"

#define UNDEFINEDMODEL -1
#define UNDERDEFINING  -2 

#define DSMCGASMODEL  0
#define EULERGASMODEL 1 

using namespace std;
extern int DIM;

template <class DataType=float,class NodeType=Cnode<DataType>,class FaceType=Cface<DataType,NodeType> >
class Ccell{
protected:
  double measure_value;

public:
  long int nodeno[4],faceno[4],neighbour_cellno[4];
  NodeType* node[4];
  FaceType* face[4];

  int LocalGasModel;

  double midpoint[3];
  bool midpointinitflag;

  Ccell() {
    midpointinitflag=false;
    measure_value=-1.0;
    LocalGasModel=-1;
    for (int i=0;i<4;i++) nodeno[i]=-1,faceno[i]=-1,neighbour_cellno[i]=-1,node[i]=NULL,face[i]=NULL;
  };

  double *GetMidPoint() {
    if (midpointinitflag==false) {
      midpointinitflag=true;
      GetCellCenter(midpoint);
    }

    return midpoint;
  };

  void GetCellCenter(double* x) {
    int i,pnode;
    double xnode[3],c=1.0/(DIM+1);

    for (i=0;i<DIM;i++) x[i]=0.0;

    for (pnode=0;pnode<DIM+1;pnode++) {
      node[pnode]->GetX(xnode);
      for (i=0;i<DIM;i++) x[i]+=c*xnode[i];
    } 
  }; 

  void GetCellCenter(float* x) {
    int i;
    float center[3];

   GetCellCenter(center);
   for (i=0;i<DIM;i++) x[i]=center[i]; 
  }; 

  void SetLocalGasModel(int model) {LocalGasModel=model;}; 
  int GetLocalGasModel() {return LocalGasModel;}; 

  void SetVolume_dim0(double volume) {
    measure_value=volume;
  }; 

  void InitCellNormal() {
    int nface,nnode;
    DataType m;
    array_1d<DataType> r1(DIM),r2(DIM),n(DIM);

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
        exit(0);
      } 
    } 
  }; 

  double Measure() {
    static array_1d<DataType> x_1d(1);
    static array_1d<DataType> x1_2d(2),x2_2d(2);
    static array_1d<DataType> x1_3d(3),x2_3d(3),x3_3d(3);

    if (measure_value>0.0) return measure_value;

    switch (DIM) {
    case 0 :
      printf("Error: Ccell::Measure()\n");
      printf("measure_value dor DIM=0 is not defined\n");
      exit(0); 

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
  }; 

  double CharacteristicSize() {
    double res;
    static array_1d<DataType> x_1d(1);
    static array_1d<DataType> x1_2d(2),x2_2d(2);
    static array_1d<DataType> x1_3d(3),x2_3d(3),x3_3d(3);


    switch (DIM) {
    case 0 :
      res=1.0;
      break;
    case 1 :
      x_1d=node[1]->X()-node[0]->X();
      res=x_1d.abs();  

      break;
    case 2 :
      x1_2d=node[1]->X()-node[0]->X();
      x2_2d=node[2]->X()-node[0]->X();
      res=sqrt(0.5*fabs(x1_2d(0)*x2_2d(1)-x1_2d(1)*x2_2d(0))/Pi);

      break;
    case 3 :
      x1_3d=node[1]->X()-node[0]->X();
      x2_3d=node[2]->X()-node[0]->X();
      x3_3d=node[3]->X()-node[0]->X();
      res=pow(fabs(mix_product(x1_3d,x2_3d,x3_3d))/8.0/Pi,0.3333333);

      break;
    }

    return res;
  }; 

  void RandomPosition(float* x) {
    double xtmp[3];
    int i;

    RandomPosition(xtmp);
    for (i=0;i<DIM;i++) x[i]=xtmp[i];
  };

  void RandomPosition(double* x) {
    double f1,f2,f3; 

    switch (DIM) {
    case 3 :
      f1=rnd(); if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999; f1=pow(f1,(double)0.3333333333);
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
      exit(0);
    }
  }; 

  friend bool CellsIntersection(Ccell<DataType,NodeType,FaceType>& cell1,Ccell<DataType,NodeType,FaceType>& cell2) {
    bool res; 
    int idim,pnode;
    DataType cells1_xmin[3],cells1_xmax[3];
    DataType cells2_xmin[3],cells2_xmax[3];
    array_1d<DataType> x(DIM);

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

    res=true;
    for (idim=0;idim<DIM;idim++)
      if ((cells1_xmax[idim]<cells2_xmin[idim])||(cells1_xmin[idim]>cells2_xmax[idim])) res=false;

    return res;
  }; 
};

#endif
