//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//the cut-cell's functions
#include <iostream>
#include <list>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "meshAMRdef.h"


#ifndef _CUT_CELL_MESHAMR_
#define _CUT_CELL_MESHAMR_

namespace CutCell {

  struct cNodeCoordinates {
    double *x;
    int id,pic__shadow_attribute;

    int nface;
  };

  struct cFaceNodeConnection {
    list<cNodeCoordinates>::iterator node[3];
  };

  struct cNASTRANnode {
    double x[3];
    double BallAveragedExternalNormal[3];
    int id;
  };

  struct cNASTRANface {
    int node[3],faceat;
    double externalNormal[3];
  };

  class cCutBlockNode {
  public:
    double x[3];
    int id;

    cCutBlockNode() {
      for (int idim=0;idim<3;idim++) x[idim]=0.0;
      id=-1;
    }
  };

  class cCutData {
  public:
    double x0[3],norm[3];

    cCutData() {
      for (int idim=0;idim<3;idim++) x0[idim]=0.0,norm[idim]=0.0;
    }
  };



  class cCutEdge {
  public:
    cCutBlockNode* cutPoint;
    cCutData* cutData;

    cCutEdge() {
      cutPoint=NULL,cutData=NULL;
    }
  };

  class cTetrahedron {
  public:
    list<cCutBlockNode>::iterator node[4];

    cTetrahedron() {
  //    for (int i=0;i<3;i++) node[i]=NULL;
    }

    double Volume() {
      double e1[3],e2[3],e3[3],*x0,*x1,*x2,*x3;

      x0=node[0]->x,x1=node[1]->x,x2=node[2]->x,x3=node[3]->x;

      for (int i=0;i<3;i++) {
        e1[i]=x1[i]-x0[i];
        e2[i]=x2[i]-x0[i];
        e3[i]=x3[i]-x0[i];
      }

      return fabs(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6.0;
    }
  };


  class cCutBlock {
  public:
    list<cCutBlockNode>::iterator node[3][3][3];
    cCutEdge* edge[12];
    double dxBlock[3],xBlockMin[3],xBlockMax[3];

    list<cCutBlockNode> NodeBuffer;

    cCutBlock() {
      for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) node[i][j][k]=NodeBuffer.end();
      for (int i=0;i<12;i++) edge[i]=NULL;
    }

    list<cCutBlockNode>::iterator AddNode(double *x,int i,int j,int k) {
       if (node[i][j][k]!=NodeBuffer.end()) exit(__LINE__,__FILE__,"Error: redefinition if the node");

       cCutBlockNode nd;

       memcpy(nd.x,x,3*sizeof(double));
       NodeBuffer.push_front(nd);

       node[i][j][k]=NodeBuffer.begin();
       return node[i][j][k];
    }

    list<cCutBlockNode>::iterator AddNode(double x0,double x1,double x2,int i,int j,int k) {
      double x[3];

      x[0]=x0,x[1]=x1,x[2]=x2;
      return AddNode(x,i,j,k);
    }

    void Reset() {
      for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) node[i][j][k]=NodeBuffer.end();
    }

    void ClearNode(int i,int j, int k) {
      if (node[i][j][k]!=NodeBuffer.end()) {
        NodeBuffer.erase(node[i][j][k]);
        node[i][j][k]=NodeBuffer.end();
      }
    }
  };


  class cTriangleFace;

  class cTriangleEdge {
    cNASTRANnode *node[2];
    cTriangleFace *face[2];

    cTriangleEdge() {
      for (int i=0;i<2;i++) node[i]=NULL,face[i]=NULL;
    }
  };

  class cTriangleFace {
  public:
  //  list<cCutBlockNode>::iterator node[3];
    cNASTRANnode *node[3];
    list<cTriangleEdge>::iterator edge[3];


    double ExternalNormal[3],SurfaceArea;
    int attribute;

    double x0Face[3],x1Face[3],x2Face[3];
    double e0[3],e1[3],c00,c01,c11,c,e0Length,e1Length;

    double CharacteristicFaceSize;

    cTriangleFace *next,*prev;

    //the variables used by AMPS to determine the surface elements that are in the shadow. The values are modified by pic__ray_tracing.cpp
    unsigned int pic__shadow_attribute,pic__RayTracing_TestDirectAccessCounterValue;

    void GetCenterPosition(double *x) {
      for (int idim=0;idim<3;idim++) x[idim]=x0Face[idim]+x1Face[idim]/2.0+x2Face[idim]/4.0;
    }

    void GetRandomPosition(double *x,double EPS=0.0) {
      double xLocal[2];

      xLocal[0]=1.0-sqrt(rnd());
      xLocal[1]=rnd()*(1.0-xLocal[0]);

      for (int idim=0;idim<3;idim++) x[idim]=x0Face[idim]+xLocal[0]*e0[idim]+xLocal[1]*e1[idim]  +   EPS*ExternalNormal[idim];
    }

    void GetRandomPosition(double *x,double *LocalNorm,double EPS=0.0) {
      double xLocal[2];

      //get the local position of the point
      xLocal[0]=1.0-sqrt(rnd());
      xLocal[1]=rnd()*(1.0-xLocal[0]);

      //get the interpolated value of the BollAvaragedExternalNormal
      double f01,f12,l=0.0;
      int idim;

      for (idim=0;idim<3;idim++) {
        f01=(1.0-xLocal[0])*node[0]->BallAveragedExternalNormal[idim]+xLocal[0]*node[1]->BallAveragedExternalNormal[idim];
        f12=(1.0-xLocal[0])*node[2]->BallAveragedExternalNormal[idim]+xLocal[0]*node[1]->BallAveragedExternalNormal[idim];

        LocalNorm[idim]=(1.0-xLocal[1])*f01+xLocal[1]*f12;
        l+=pow(LocalNorm[idim],2);
      }

      for (l=sqrt(l),idim=0;idim<3;idim++) LocalNorm[idim]/=l;

      for (int idim=0;idim<3;idim++) x[idim]=x0Face[idim]+xLocal[0]*e0[idim]+xLocal[1]*e1[idim]  +   EPS*ExternalNormal[idim];
    }


    void SetFaceNodes(double *x0,double *x1,double *x2) {
      int i;
      double l;


      memcpy(x0Face,x0,3*sizeof(double));
      memcpy(x1Face,x1,3*sizeof(double));
      memcpy(x2Face,x2,3*sizeof(double));


      for (i=0;i<3;i++) {
        e0[i]=x1Face[i]-x0Face[i];
        e1[i]=x2Face[i]-x0Face[i];
      }

      ExternalNormal[0]=+(e0[1]*e1[2]-e1[1]*e0[2]);
      ExternalNormal[1]=-(e0[0]*e1[2]-e1[0]*e0[2]);
      ExternalNormal[2]=+(e0[0]*e1[1]-e1[0]*e0[1]);

      l=sqrt(ExternalNormal[0]*ExternalNormal[0]+ExternalNormal[1]*ExternalNormal[1]+ExternalNormal[2]*ExternalNormal[2]);
      SurfaceArea=l/2.0;
      ExternalNormal[0]/=l,ExternalNormal[1]/=l,ExternalNormal[2]/=l;

      c00=e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2];
      c11=e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2];
      c01=e0[0]*e1[0]+e0[1]*e1[1]+e0[2]*e1[2];
      c=c00*c11-c01*c01;
      e0Length=sqrt(c00),e1Length=sqrt(c11);

      //calculate the characteristic size of the face
      double l0=0.0,l1=0.0,l2=0.0;
      int idim;

      for (idim=0;idim<3;idim++) {
        l0+=pow(x0Face[idim]-x1Face[idim],2);
        l1+=pow(x0Face[idim]-x2Face[idim],2);
        l2+=pow(x2Face[idim]-x1Face[idim],2);
      }

      CharacteristicFaceSize=0.3*(sqrt(l0)+sqrt(l1)+sqrt(l2));

    }

    inline void GetProjectedLocalCoordinates(double *xLocal,double *x) {
      double c0=0.0,c1=0.0,t;

      for (int i=0;i<3;i++) {
        t=x[i]-x0Face[i];
        c0+=t*e0[i],c1+=t*e1[i];
      }

      c=c11*c00-c01*c01;
      xLocal[0]=(c0*c11-c1*c01)/c;
      xLocal[1]=(c1*c00-c01*c0)/c;
    }

    bool RayIntersection(double *x0,double *l,double &t,double EPS) {
      double length,lNorm;

      lNorm=l[0]*ExternalNormal[0]+l[1]*ExternalNormal[1]+l[2]*ExternalNormal[2];
      length=(x0[0]-x0Face[0])*ExternalNormal[0]+(x0[1]-x0Face[1])*ExternalNormal[1]+(x0[2]-x0Face[2])*ExternalNormal[2];

/*      if (fabs(length)<EPS) {
        t=0.0;
      }
      else {
        if (lNorm==0.0) {
          t=0.0;
          return false;
        }

        t=-length/lNorm;
        if (t<0.0) return false;
      }*/

/*      if (lNorm==0.0) {
        t=0.0;
        return true;
      }*/


      if (lNorm==0.0) {
        t=0.0;
        return false;
      }
      else {
        t=-length/lNorm;
      }

      if (t<0.0) return false;

      //find position of the intersection in the internal frame related to the face
/*      double c0=0.0,c1=0.0,xLocal[2];

      for (int i=0;i<3;i++) {
        double x=x0[i]+l[i]*t-x0Face[i];

        c0+=x*e0[i],c1+=x*e1[i];
      }

      c=c11*c00-c01*c01;
      xLocal[0]=(c0*c11-c1*c01)/c;
      xLocal[1]=(c1*c00-c01*c0)/c;*/

      double x[3],xLocal[2];

      for (int i=0;i<3;i++) x[i]=x0[i]+l[i]*t;
      GetProjectedLocalCoordinates(xLocal,x);

      //determine weather the node in outside of the face
/*      double a0=e0Length*xLocal[0],a1=e1Length*xLocal[1];

      if ((a0<-EPS)||(a0>EPS+e0Length)) return false;
      else if ((a1<-EPS)||(a1>EPS+e1Length)) return false;
      else if (xLocal[0]+xLocal[1]>1.0) return false;*/

      if ((xLocal[0]<0.0)||(xLocal[0]>1.0) || (xLocal[1]<0.0)||(xLocal[1]>1.0) || (xLocal[0]+xLocal[1]>1.0)) return false;

      return true;
    }

    inline bool RayIntersection(double *x0,double *l,double EPS) {
      double t;

      return RayIntersection(x0,l,t,EPS);
    }

    inline bool IntervalIntersection(double *x0,double *x1,double& t,double EPS) {
      double l[3];
      bool res=false;

      l[0]=x1[0]-x0[0],l[1]=x1[1]-x0[1],l[2]=x1[2]-x0[2];
      if (RayIntersection(x0,l,t,EPS)==true) if ((0.0<=t)&&(t<=1.0)) res=true;

      return res;
    }

    inline bool IntervalIntersection(double *x0,double *x1,double *xIntersection,double EPS) {
      double t;
      bool res;

      res=IntervalIntersection(x0,x1,t,EPS);
      if (res==true) for (int i=0;i<3;i++) xIntersection[i]=x0[i]+t*(x1[i]-x0[i]);

      return res;
    }

    inline bool IntervalIntersection(double *x0,double *x1,double EPS) {
      double t;

      return IntervalIntersection(x0,x1,t,EPS);
    }


    bool BlockIntersection(double *xmin,double *xmax,double EPS) {
      int idim,pedge,pface;
      bool flag;

      //check if any of the face's node is within the block
      bool x0flag=true,x1flag=true,x2flag=true;

      for (idim=0;idim<3;idim++) {
        if ((x0Face[idim]<xmin[idim])||(xmax[idim]<x0Face[idim])) x0flag=false;
        if ((x1Face[idim]<xmin[idim])||(xmax[idim]<x1Face[idim])) x1flag=false;
        if ((x2Face[idim]<xmin[idim])||(xmax[idim]<x2Face[idim])) x2flag=false;
      }

      if ((x0flag==true)||(x1flag==true)||(x2flag==true)) return true;

      //check intersection of the face with the edges of the block
      static const double x0EdgeMap[12][3]={{0,0,0},{0,1,0},{0,1,1},{0,0,1},  {0,0,0},{1,0,0},{1,0,1},{0,0,1},  {0,0,0},{1,0,0},{1,1,0},{0,1,0}};
      static const double x1EdgeMap[12][3]={{1,0,0},{1,1,0},{1,1,1},{1,0,1},  {0,1,0},{1,1,0},{1,1,1},{0,1,1},  {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

      double x0[3],x1[3],t;

      for (pedge=0;pedge<12;pedge++) {
        for (idim=0;idim<3;idim++) {
          x0[idim]=xmin[idim]+x0EdgeMap[pedge][idim]*(xmax[idim]-xmin[idim]);
          x1[idim]=xmin[idim]+x1EdgeMap[pedge][idim]*(xmax[idim]-xmin[idim]);
        }

        if (IntervalIntersection(x0,x1,EPS)==true) return true;
      }

      //check intersection of the edges of the face with the faces of the block
      static const double x0Block[6][3]={{0,0,0},{1,0,0},   {0,0,0},{0,1,0},   {0,0,0},{0,0,1}};
      static const double normBlock[6][3]={{1,0,0},{1,0,0},   {0,1,0},{0,1,0},    {0,0,1},{0,0,1}};
  //    static const double e0Block[6][3]={{0,1,0},{1,1,0},    {1,0,0}, {1,1,0},   {1,0,0},{1,0,1}};
  //    static const double e1Block[6][3]={{0,0,1},{1,0,1},     {0,0,1},{0,1,1},   {0,1,0},{0,1,1}};
      static const int skipBlockFaceDirection[6]={0,0,1,1,2,2};

      double x0Edge[3],lEdge[3],lNorm,xNorm;

      for (pface=0;pface<6;pface++) for (int pFaceEdge=0;pFaceEdge<3;pFaceEdge++) {
        lNorm=0.0,xNorm=0.0;

        switch (pFaceEdge) {
        case 0:
          memcpy(x0Edge,x1Face,3*sizeof(double));
          for (int i=0;i<3;i++) lEdge[i]=x2Face[i]-x1Face[i];
          break;
        case 1:
          memcpy(x0Edge,x0Face,3*sizeof(double));
          for (int i=0;i<3;i++) lEdge[i]=x2Face[i]-x0Face[i];
          break;
        case 2:
          memcpy(x0Edge,x0Face,3*sizeof(double));
          for (int i=0;i<3;i++) lEdge[i]=x1Face[i]-x0Face[i];
        }

        for (int i=0;i<3;i++) {
          lNorm+=lEdge[i]*normBlock[pface][i];
          xNorm+=(x0Edge[i]- (xmin[i]+x0Block[pface][i]*(xmax[i]-xmin[i])) )*normBlock[pface][i];
        }

        if (fabs(lNorm)<=fabs(xNorm)) continue;

        if (fabs(xNorm)<EPS) t=0.0;
        else {
          t=-xNorm/lNorm;
        }

        if (t<0.0) continue;



        for (flag=true,idim=0;idim<3;idim++) if (idim!=skipBlockFaceDirection[pface]) {
          double x=x0Edge[idim]+t*lEdge[idim];

          if ((x<xmin[idim])||(xmax[idim]<x)) {
            flag=false;
            break;
          }
        }

        if (flag==true) return true;
      }

      return false;
    }



    cTriangleFace() {
      SurfaceArea=0.0,CharacteristicFaceSize=0.0;
      for (int i=0;i<3;i++) ExternalNormal[i]=0.0;

      pic__shadow_attribute=0,pic__RayTracing_TestDirectAccessCounterValue=0;
    }

    inline double CharacteristicSize() {
      return CharacteristicFaceSize;
    }
  };

  class cQuadrilateral {
  public:
    list<cCutBlockNode>::iterator node[8];
  };

  class cTriangleCutFace : public cTriangleFace {
  public:
    list<cCutBlockNode>::iterator node[3];
  };

  extern cTriangleFace *BoundaryTriangleFaces;
  extern int nBoundaryTriangleFaces;

  extern cNASTRANnode *BoundaryTriangleNodes;
  extern int nBoundaryTriangleNodes;

  extern list<cTriangleEdge> BoundaryTriangleEdges;

  void PrintSurfaceTriangulationMesh(const char *fname,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS);
  void PrintSurfaceTriangulationMesh(const char *fname);

  void ReadNastranSurfaceMeshLongFormat(const char *fname,double *xSurfaceMin,double *xSurfaceMax,double EPS=0.0);
  void ReadNastranSurfaceMeshLongFormat(const char *fname);

  bool CheckPointInsideDomain(double *x,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS);
  bool GetClosestSurfaceIntersectionPoint(double *x0,double *lSearch,double *xIntersection,double &tIntersection,cTriangleFace* &FaceIntersection,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS=0.0);


  //cderive the connectivity list
  class cConnectivityListTriangleFace;
  class cConnectivityListTriangleEdge;
  class cConnectivityListTriangleNode;
  class cConnectivityListTriangleEdgeDescriptor;

  class cConnectivityListTriangleNode {
  public:
    cNASTRANnode *node;
    list<cConnectivityListTriangleEdgeDescriptor>::iterator ballEdgeList;
  };

  class cConnectivityListTriangleEdgeDescriptor {
  public:
    list<cConnectivityListTriangleEdge>::iterator edge;
    list<cConnectivityListTriangleEdgeDescriptor>::iterator next;
  };

  class cConnectivityListTriangleEdge {
  public:
    list<cConnectivityListTriangleFace>::iterator face[2];
    list<cConnectivityListTriangleNode>::iterator node[2];

    void GetMiddleNode(double *x) {
      for (int i=0;i<3;i++) x[i]=0.5*(node[0]->node->x[i]+node[1]->node->x[i]);
    }

    double GetLength() {
      double s=0.0;

      for (int i=0;i<3;i++) s+=pow(node[0]->node->x[i]+node[1]->node->x[i],2);
      return sqrt(s);
    }
  };

  class cConnectivityListTriangleFace {
  public:
    list<cConnectivityListTriangleFace>::iterator neib[3];
    list<cConnectivityListTriangleEdge>::iterator edge[3];
    list<cConnectivityListTriangleNode>::iterator node[3];

    double GetLength() {
      double s=0.0;

      for (int i=0;i<3;i++) s+=edge[i]->GetLength();
      return s/3.0;
    }

  };

  void ReconstructConnectivityList(list<cConnectivityListTriangleNode>& nodes,list<cConnectivityListTriangleEdge>& edges,list<cConnectivityListTriangleFace>& faces,list<cConnectivityListTriangleEdgeDescriptor>& RecoveredEdgeDescriptorList);

  //refine the smooth the surface mesh
  class cLocalTriangle;


  class cLocalNode {
  public:
    cNASTRANnode *OriginalNode;
    double x[3],xSmooth[3];
    vector<vector<cLocalTriangle>::iterator> ball;
    int OriginalNodeID,nodeno;

    cLocalNode() {
      OriginalNode=NULL;
      for (int idim=0;idim<3;idim++) x[idim]=0.0,xSmooth[idim]=0.0;
      OriginalNodeID=-1,nodeno=-1;
    }
  };

  class cLocalEdge {
  public:
    vector<cLocalNode>::iterator CornerNode[2],MiddleNode;
    vector<cLocalEdge>::iterator downEdge[2];
    bool processed;

    cLocalEdge() {
      processed=false;
    }
  };

  class cLocalTriangle{
  public:
    vector<cLocalNode>::iterator node[3];
    vector<cLocalEdge>::iterator edge[3];
    vector<cLocalTriangle>::iterator upTriangle;
    cTriangleFace *TriangleFace;

    cLocalTriangle() {
      TriangleFace=NULL;
    }
  };

  void SmoothRefine(double SmoothingCoefficient);
  void SmoothMeshResolution(double MaxNeibSizeRatio);


  class cTriangleFaceDescriptor {
  public:
    cTriangleFace *TriangleFace;
    cTriangleFaceDescriptor *next,*prev;

    int Temp_ID;

    void cleanDataBuffer() {
      TriangleFace=NULL,Temp_ID=-1,next=NULL,prev=NULL;
    }

    cTriangleFaceDescriptor() {
      cleanDataBuffer();
    }
  };

  extern cAMRstack<cTriangleFaceDescriptor> BoundaryTriangleFaceDescriptor;

  double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,cTriangleFaceDescriptor* TriangleCutFaceDescriptorList);
  //double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,cTriangleFaceDescriptor* TriangleCutFaceDescriptorList,int maxIntegrationLevel,int IntegrationLevel);
  double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, list<cTriangleFace*>& BlockTriangulationSet,int maxIntegrationLevel,int IntegrationLevel);

  int cutBlockTetrahedronConnectivity(CutCell::cCutBlock* bl,list<CutCell::cTetrahedron>& indomainConnectivityList,list<CutCell::cTetrahedron>& outdomainConnectivityList,list<CutCell::cTriangleCutFace> &TriangleCutFaceConnectivity);

  //===============================================================================================================
  //print the volume mesh


  template <class T>
  struct cNodeData {
    typename list<T>::iterator node;
    int id;
  };


  template <class T>
  struct cCellData {

    typename list<T>::iterator node[4];

  };

  template <class cCutBlockNode,class cTetrahedron>
  void cutBlockPrintVolumeMesh(const char *fname,list<cTetrahedron>& ConnectivityList) {
    FILE *fout;
    int i;
    typename list<cTetrahedron>::iterator ConnectivityIterator;
    bool found;




    list<cNodeData<cCutBlockNode> > nodes;
    list<cCellData<cNodeData<cCutBlockNode> > > cells;

    typename list<cNodeData<cCutBlockNode> >::iterator nodeitr;
    typename list<cCellData<cNodeData<cCutBlockNode> > >::iterator cellitr;
    int idMax=0;



    for (ConnectivityIterator=ConnectivityList.begin();ConnectivityIterator!=ConnectivityList.end();ConnectivityIterator++) {
      cCellData<cNodeData<cCutBlockNode> > cl;

      for (i=0;i<4;i++) {
        for (found=false,nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) if (ConnectivityIterator->node[i]==nodeitr->node) {
          found=true;
          cl.node[i]=nodeitr;
          break;
        }

        if (found==false) {
          //add the node to the list
          cNodeData<cCutBlockNode> nd;

          nd.node=ConnectivityIterator->node[i];

          nodes.push_front(nd);
          cl.node[i]=nodes.begin();
        }
      }

      cells.push_back(cl);
    }


    //numerate the nodes
    for (idMax=0,nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) nodeitr->id=++idMax;

    //output the mesh
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"\nZONE N=%ld, E=%ld,DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON\n",idMax,cells.size());

    for (nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) fprintf(fout,"%e %e %e\n",nodeitr->node->x[0],nodeitr->node->x[1],nodeitr->node->x[2]);

    for (cellitr=cells.begin();cellitr!=cells.end();cellitr++) {
      fprintf(fout,"%i %i %i %i\n",cellitr->node[0]->id,cellitr->node[1]->id,cellitr->node[2]->id,cellitr->node[3]->id);
    }

    fclose(fout);
  }

  template <class cCutBlockNode,class cTetrahedron>
  bool CheckConnectivityListIntersection(double *xmin,double *xmax,list<cTetrahedron>** ConnectivityList,int nConnectivityList,double &ConnectivityListVolume,double &VolumeFraction) {
    double x[3],e1[3],e2[3],e3[3],rhs[3],*x0,*x1,*x2,*x3,LocalCoordinates[3],A;
    int i,ntest,nList,cntIntersection=0,cnt;
    typename list<cTetrahedron>::iterator ConnectivityIterator;
    list<cTetrahedron> IntersectionConnectivity;

    static const int nTotalChecks=100000;

    //check cTetrahedron intersection
    for (ntest=0;ntest<nTotalChecks;ntest++) {
      for (i=0;i<3;i++) x[i]=xmin[i]+rnd()*(xmax[i]-xmin[i]);

      for (cnt=0,nList=0;nList<nConnectivityList;nList++) for (ConnectivityIterator=ConnectivityList[nList]->begin();ConnectivityIterator!=ConnectivityList[nList]->end();ConnectivityIterator++) {
        x0=ConnectivityIterator->node[0]->x;
        x1=ConnectivityIterator->node[1]->x;
        x2=ConnectivityIterator->node[2]->x;
        x3=ConnectivityIterator->node[3]->x;

        for (i=0;i<3;i++) {
          e1[i]=x1[i]-x0[i];
          e2[i]=x2[i]-x0[i];
          e3[i]=x3[i]-x0[i];

          rhs[i]=x[i]-x0[i];
        }

        A=e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]);
        if (fabs(A)<1.0E-15) exit(__LINE__,__FILE__,"Error: found Tetrahedron with zero volume");

        LocalCoordinates[0]=rhs[0]*(e2[1]*e3[2]-e2[2]*e3[1])-rhs[1]*(e2[0]*e3[2]-e2[2]*e3[0])+rhs[2]*(e2[0]*e3[1]-e2[1]*e3[0]);
        LocalCoordinates[1]=e1[0]*(rhs[1]*e3[2]-rhs[2]*e3[1])-e1[1]*(rhs[0]*e3[2]-rhs[2]*e3[0])+e1[2]*(rhs[0]*e3[1]-rhs[1]*e3[0]);
        LocalCoordinates[2]=e1[0]*(e2[1]*rhs[2]-e2[2]*rhs[1])-e1[1]*(e2[0]*rhs[2]-e2[2]*rhs[0])+e1[2]*(e2[0]*rhs[1]-e2[1]*rhs[0]);

        if ((LocalCoordinates[0]>=0.0)&&(LocalCoordinates[1]>=0.0)&&(LocalCoordinates[2]>=0.0)&&(LocalCoordinates[0]+LocalCoordinates[1]+LocalCoordinates[2]<=1.0)) {
          cnt++;
        }
      }

      if (cnt>=2) {
        cntIntersection++;
        IntersectionConnectivity.push_back(*ConnectivityIterator);
      }
    }

    //calculate the total connectivity list volume
    for (ConnectivityListVolume=0.0,nList=0;nList<nConnectivityList;nList++) for (ConnectivityIterator=ConnectivityList[nList]->begin();ConnectivityIterator!=ConnectivityList[nList]->end();ConnectivityIterator++) {
      x0=ConnectivityIterator->node[0]->x;
      x1=ConnectivityIterator->node[1]->x;
      x2=ConnectivityIterator->node[2]->x;
      x3=ConnectivityIterator->node[3]->x;

      for (i=0;i<3;i++) {
        e1[i]=x1[i]-x0[i];
        e2[i]=x2[i]-x0[i];
        e3[i]=x3[i]-x0[i];

        rhs[i]=x[i]-x0[i];
      }

      ConnectivityListVolume+=fabs(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6.0;
    }

    VolumeFraction=(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2])/ConnectivityListVolume;

    if (cntIntersection!=0) {
      cutBlockPrintVolumeMesh<cCutBlockNode,cTetrahedron>("CutBlockInteresection.dat",IntersectionConnectivity);
      return true;
    }

    return false;
  }
}


#endif
