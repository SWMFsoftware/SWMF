#ifndef FACE
#define FACE

#include "node.h"
#include "array_1d.h"
#include "specfunc.h"

#include "const.dfn"
#include "global.dfn"

extern int SymmetryMode;

using namespace std;

template<class DataType=float,class NodeType=Cnode<DataType> >
class Cface{
  DataType nrml[3];
  double measure;

  double midpoint[3]; 
  bool midpointinitflag;
   
public:
  unsigned char faceat,surface_group;
  long int faceno,nodeno[3];
  NodeType* node[3];

  Cface() {
    midpointinitflag=false;
    measure=-1.0;
  };

  double* GetMidPoint() {
    if (midpointinitflag==false) {
      DataType nodex[3];
      int pnode,idim;

      for (idim=0;idim<DIM;idim++) midpoint[idim]=0.0;

      for (pnode=0;pnode<DIM;pnode++) {
        node[pnode]->GetX(nodex);
        for (idim=0;idim<DIM;idim++) midpoint[idim]+=nodex[idim]/DIM;       
      }

      midpointinitflag=true;
    }

    return midpoint;
  };

//===================================================
  void SetNormal(array_1d<DataType>& n) {
    int idim;

    for (idim=0;idim<DIM;idim++) nrml[idim]=n(idim);
  };
//===================================================
  array_1d<DataType> Normal() {
    array_1d<DataType> normal(DIM);
    int idim;

    for (idim=0;idim<DIM;idim++) normal(idim)=nrml[idim];
    return normal;
  };

//===================================================
  void GetNormal(DataType* n) {
    for (int idim=0;idim<DIM;idim++) n[idim]=nrml[idim];
  }; 
//===================================================
  double Measure() {
    if (measure>0.0) return measure;

    array_1d<DataType> x_2d(2);
    array_1d<DataType> x1_3d(3),x2_3d(3);

    switch(DIM) {
    case 1:  
      if (SymmetryMode==no_symmetry) measure=1.0;
      else {
        double x[3];
        node[0]->GetX(x);
        measure=4.0*Pi*pow(x[0],2);
      }

      break;
    case 2:
      if (SymmetryMode==no_symmetry) {
        x_2d=(*node[1]).X()-(*node[0]).X();
        measure=x_2d.abs(); 
      }
      else if (SymmetryMode==cylindrical_symmetry) {
        double l,alfa,e[2],nd0[2];
  
        node[0]->GetX(nd0);node[1]->GetX(e);
        e[0]-=nd0[0],e[1]-=nd0[1];
        l=sqrt(e[0]*e[0]+e[1]*e[1]); 

        if (fabs(e[0]/l)<1.0E-8) measure=Pi*fabs((2.0*nd0[1]+e[1])*e[1]);
        else {
          alfa=e[1]/e[0];
          measure=2.0*Pi*fabs(e[0]*(nd0[1]+alfa*e[0]/2.0))*sqrt(1.0+alfa*alfa);
        } 
      }

      break;
    case 3:
      x1_3d=(*node[1]).X()-(*node[0]).X();
      x2_3d=(*node[2]).X()-(*node[0]).X();
      measure=0.5*fabs(cross_product(x1_3d,x2_3d).abs());
      break;
    default : 
      printf("proc. Cface::Measure()\n"); 
      printf("wrong DIM value: DIM=%i\n",DIM);
      exit(__LINE__,__FILE__);
    }

    return measure;
  };
//===================================================
  template <class T> void RandomPosition(T* x) {

    switch (DIM) {
    case 3 :
      x[0]=1.0-sqrt(rnd());
      x[1]=rnd()*(1.0-x[0]); 

      return;
    case 2 :
      if (SymmetryMode==no_symmetry) x[0]=rnd();     
      else if (SymmetryMode==cylindrical_symmetry) {
        register double dr,l,nd0[2],nd1[2],e[2];

        node[0]->GetX(nd0);node[1]->GetX(nd1);
        e[0]=nd1[0]-nd0[0],e[1]=nd1[1]-nd0[1];
        l=sqrt(e[0]*e[0]+e[1]*e[1]);
        dr=e[1]; 

        x[0]=(fabs(dr/l)<1.0E-8) ? rnd() : (-nd0[1]+sqrt(pow(nd0[1],2)+2.0*dr*rnd()*(nd0[1]+dr/2.0)))/dr;    
      }

      if (x[0]<0.00001E0) x[0]=0.00001E0; 
      if (x[0]>0.99999E0) x[0]=0.99999E0;
      return;
    default :
      exit(__LINE__,__FILE__,"wrong DIM value");
    }

  };
//===================================================
};
#endif
