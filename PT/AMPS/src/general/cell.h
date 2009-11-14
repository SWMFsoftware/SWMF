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
  long int cellno,nodeno[4],faceno[4],neighbour_cellno[4];
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
    double center[3];

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
//    array_1d<DataType> r1(DIM),r2(DIM),n(DIM);

    for (nface=0;nface<DIM+1;nface++) {
      switch (DIM) {
      case 3 :

exit(__LINE__,__FILE__);

/*
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
*/
        break;
      case 2 :

/*
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
*/


double r1vect[3],r2vect[3],nvect[3],nnorm;

face[nface]->node[1]->GetX(r1vect);
face[nface]->node[0]->GetX(r2vect);
r1vect[0]-=r2vect[0],r1vect[1]-=r2vect[1],r1vect[2]-=r2vect[2]; 
nvect[0]=r1vect[1],nvect[1]=-r1vect[0];
nnorm=sqrt(pow(nvect[0],2)+pow(nvect[1],2));
nvect[0]/=nnorm,nvect[1]/=nnorm;

m=-1.0;
for (nnode=0;nnode<DIM;nnode++) {
  face[nface]->node[nnode]->GetX(r1vect);
  node[nface]->GetX(r2vect);
  r1vect[0]-=r2vect[0],r1vect[1]-=r2vect[1],r1vect[2]-=r2vect[2];
  if (r1vect[0]*nvect[0]+r1vect[1]*nvect[1]>0.0) m=1.0; 
}
nvect[0]*=m,nvect[1]*=m;
face[nface]->SetNormal(nvect);





        break;
      case 1 :

//exit(__LINE__,__FILE__);



///////////
 {
 double r1vect[3]={0.0,0.0,0.0};
 double r0vect[3]={0.0,0.0,0.0};
 
 node[1]->GetX(r1vect);
 node[0]->GetX(r0vect);
 r1vect[0]-=r0vect[0];
 if (nface==1) r1vect[0]*=-1.0;
 r1vect[0]/=fabs(r1vect[0]);
 face[nface]->SetNormal(r1vect); 
 }
 
 
 
 ////////////

/*

        n=node[1]->X()-node[0]->X();
        if (nface==1) n*=-1;
        n.normalize();
        face[nface]->SetNormal(n);
*/

        break;
      default :
        printf("proc. Ccell::InitCellNormal()\n");
        printf("wrong DIM value: DIM=%i\n",DIM);
        exit(__LINE__,__FILE__);
      } 
    } 
  }; 

  double Measure() {
    if (measure_value>0.0) return measure_value;

/*
    array_1d<DataType> x_1d(1);
    array_1d<DataType> x1_2d(2),x2_2d(2);
    array_1d<DataType> x1_3d(3),x2_3d(3),x3_3d(3);
*/

    switch (DIM) {
    case 0 :
      measure_value=1.0;

      break;
    case 1 :

//exit(__LINE__,__FILE__);

/*

      if (SymmetryMode==no_symmetry) {
        x_1d=node[1]->X()-node[0]->X();
        measure_value=x_1d.abs();
      }
      else {
        double nd0[3],nd1[3];
        node[0]->GetX(nd0);node[1]->GetX(nd1);
        measure_value=4.0/3.0*Pi*fabs(pow(nd0[0],3)-pow(nd1[0],3));
      } 
*/


      if (SymmetryMode==no_symmetry) {
        double nd0[3],nd1[3];
        node[0]->GetX(nd0);node[1]->GetX(nd1);
        measure_value=fabs(nd0[0]-nd1[0]);
      }
      else {
        double nd0[3],nd1[3];
        node[0]->GetX(nd0);node[1]->GetX(nd1);
        measure_value=4.0/3.0*Pi*fabs(pow(nd0[0],3)-pow(nd1[0],3));
      }


      break;
    case 2 :

/* 
      x1_2d=node[1]->X()-node[0]->X();
      x2_2d=node[2]->X()-node[0]->X();
      measure_value=0.5*fabs(x1_2d(0)*x2_2d(1)-x1_2d(1)*x2_2d(0));

      if (SymmetryMode==cylindrical_symmetry) {
        double r=(node[0]->X()(1)+node[1]->X()(1)+node[2]->X()(1))/3.0;  
        measure_value*=2.0*Pi*r;
      } 
*/

      double x1_2d[3],x2_2d[3],t[3];

      node[0]->GetX(t);
      node[1]->GetX(x1_2d);
      node[2]->GetX(x2_2d);

      x1_2d[0]-=t[0],x1_2d[1]-=t[1];
      x2_2d[0]-=t[0],x2_2d[1]-=t[1];

      measure_value=0.5*fabs(x1_2d[0]*x2_2d[1]-x1_2d[1]*x2_2d[0]);

      if (SymmetryMode==cylindrical_symmetry) {
//        double r=(node[0]->X()(1)+node[1]->X()(1)+node[2]->X()(1))/3.0;

        measure_value*=2.0*Pi*(3.0*t[1]+x1_2d[1]+x2_2d[1])/3.0;
      }


      break;
    case 3 :

exit(__LINE__,__FILE__);

/*
      x1_3d=node[1]->X()-node[0]->X();
      x2_3d=node[2]->X()-node[0]->X();
      x3_3d=node[3]->X()-node[0]->X();
      measure_value=fabs(mix_product(x1_3d,x2_3d,x3_3d))/6.0;

*/
      break;
    }

    return measure_value;
  }; 

  double CharacteristicSize() {
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
  }; 

  template <class T> void RandomPosition(T* x) {
    register double f1,f2; 

    switch (DIM) {
    case 3 :
      x[0]=1.0-pow(rnd(),(double)0.33333333);
      x[1]=(1-x[0])*(1.0-sqrt(rnd()));
      x[2]=rnd()*(1.0-x[0]-x[1]);   

      return;
    case 2 :
      if (SymmetryMode==no_symmetry) {
        x[0]=1.0-sqrt(rnd());
        x[1]=(1.0-x[0])*rnd();
      }
      else if (SymmetryMode==cylindrical_symmetry) {
        register double r0,ltot,p,pmax,nd0[2],nd1[2],nd2[2],e0[2],e1[2];

        node[0]->GetX(nd0);node[1]->GetX(nd1);node[2]->GetX(nd2); 
        e0[0]=nd1[0]-nd0[0],e0[1]=nd1[1]-nd0[1];
        e1[0]=nd2[0]-nd0[0],e1[1]=nd2[1]-nd0[1];

        r0=nd0[1];

        //distribute the first coordinate
        ltot=sqrt(pow(e1[0],2)+pow(e1[1],2))+sqrt(pow(e0[0],2)+pow(e0[1],2));

        f2=(fabs(e0[1]-2.0*e1[1])>1.0E-8*ltot) ? f2=(r0+e0[1]-e1[1])/(e0[1]-2.0*e1[1]) : 10.0;
        pmax=((f2>0.0)&&(f2<1.0)) ? max(r0+e0[1]/2.0,(1.0-f2)*(r0+f2*e1[1]+e0[1]*(1.0-f2)/2.0)) : r0+e0[1]/2.0; 

        do {
          f2=rnd();
          p=(1.0-f2)*(r0+f2*e1[1]+e0[1]*(1.0-f2)/2.0);
        }
        while (p/pmax<rnd()); 

        //distribute the second coordinate
        pmax=max(r0+e0[1]*(1.0-f2)+e1[1]*f2,r0+e1[1]*f2);

        do {
          f1=(1.0-f2)*rnd();
          p=r0+e0[1]*f1+e1[1]*f2;
        }
        while (p/pmax<rnd());

        if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999;
        if (f2<0.00001) f2=0.00001; if (f2>0.99999) f2=0.99999;

        x[0]=f1,x[1]=f2;
      }

      return;
    case 1 : 
      f1=rnd(); if (f1<0.00001) f1=0.00001; if (f1>0.99999) f1=0.99999;
      x[0]=f1;
      return;
    default :
      exit(__LINE__,__FILE__,"wrong DIM value");
    }
  }  

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
