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
    int id;
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




  class cTriangleFace {
  public:
  //  list<cCutBlockNode>::iterator node[3];
    cNASTRANnode *node[3];

    double ExternalNormal[3],SurfaceArea;
    int attribute;

    double x0Face[3],x1Face[3],x2Face[3];
    double e0[3],e1[3],c00,c01,c11,c,e0Length,e1Length;

    double CharacteristicFaceSize;

    cTriangleFace *next,*prev;


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


    bool RayIntersection(double *x0,double *l,double &t,double EPS) {
      double length,lNorm;

      lNorm=l[0]*ExternalNormal[0]+l[1]*ExternalNormal[1]+l[2]*ExternalNormal[2];
      length=(x0[0]-x0Face[0])*ExternalNormal[0]+(x0[1]-x0Face[1])*ExternalNormal[1]+(x0[2]-x0Face[2])*ExternalNormal[2];

      if (fabs(length)<EPS) {
        t=0.0;
      }
      else {
        if (lNorm==0.0) {
          t=0.0;
          return false;
        }

        t=-length/lNorm;
        if (t<0.0) return false;
      }

      //find position of the intersection in the internal frame related to the face
      double c0=0.0,c1=0.0,xLocal[2];

      for (int i=0;i<3;i++) {
        double x=x0[i]+l[i]*t-x0Face[i];

        c0+=x*e0[i],c1+=x*e1[i];
      }

      c=c11*c00-c01*c01;
      xLocal[0]=(c0*c11-c1*c01)/c;
      xLocal[1]=(c1*c00-c01*c0)/c;

      //determine weather the node in outside of the face
      double a0=e0Length*xLocal[0],a1=e1Length*xLocal[1];

      if ((a0<-EPS)||(a0>EPS+e0Length)) return false;
      else if ((a1<-EPS)||(a1>EPS+e1Length)) return false;
      else if (xLocal[0]+xLocal[1]>1.0) return false;

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

  void PrintSurfaceTriangulationMesh(const char *fname,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS);
  void ReadNastranSurfaceMeshLongFormat(const char *fname,double *xSurfaceMin,double *xSurfaceMax,double EPS=0.0);
  bool CheckPointInsideDomain(double *x,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS=0.0);
  bool GetClosestSurfaceIntersectionPoint(double *x0,double *lSearch,double *xIntersection,double &tIntersection,cTriangleFace* &FaceIntersection,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS=0.0);



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
  double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,cTriangleFaceDescriptor* TriangleCutFaceDescriptorList,int maxIntegrationLevel,int IntegrationLevel);

  template <class cCutBlockNode,class cCutBlock,class cCutData,class cTetrahedron,class cTriangleCutFace>
  int cutBlockTetrahedronConnectivity(cCutBlock* bl,list<cTetrahedron>& indomainConnectivityList,list<cTetrahedron>& outdomainConnectivityList,list<cTriangleCutFace> TriangleCutFaceConnectivity) {
      int nedge,idim,i,j,k;
      int connectivityLength=0;
      double c,x[3];
      long int elementPosition;



  //-------
  /*


  #define _N_TEST_NODES_ 4


  bool PRINTEDFLAG=false;
  int BLOCKPRINTED=0;
  static int PRINTEDCELLS=0;

  //int nodeIDlist[_N_TEST_NODES_]={2082,23977,23676};

  //int nodeIDlist[_N_TEST_NODES_]={359410,359343,353757};
  //int nodeIDlist[_N_TEST_NODES_]={1791752,1790002,1791751};

  //int nodeIDlist[_N_TEST_NODES_]={44770,44773,44772,246963};

  int nodeIDlist[_N_TEST_NODES_]={73780,73780,73780,-1};

  int nnds=0;

  for (int i=0;i<_N_TEST_NODES_;i++) {
    for (int ii=0;ii<3;ii++) for (int jj=0;jj<3;jj++) for (int kk=0;kk<3;kk++) if (bl->node[ii][jj][kk]!=NULL) if (bl->node[ii][jj][kk]->TEMP_IDENTIFYER_NODE==nodeIDlist[i]) ++nnds;

    for (int nedge=0;nedge<12;nedge++) if (bl->edge[nedge]->cutPoint!=NULL) if (bl->edge[nedge]->cutPoint->TEMP_IDENTIFYER_NODE==nodeIDlist[i]) ++nnds;
  }

  if (nnds==3) {






  //if (bl->TEMP_IDENTIFYER_BLOCK==116267) {

  //if ((9.76-0.01<bl->xmin[0])&&(bl->xmin[0]<9.76+0.01)&& (8.12-0.01<bl->xmin[1])&&(bl->xmin[1]<8.12+0.01)&& (2.1-0.01<bl->xmin[2])&&(bl->xmin[2]<2.1+0.01)) {

  //if (bl->node[2][2][2]->nodeno==10985) {

  //for (int i=0;i<3;i++) for (int j=0;j<3;j++) for (int k=0;k<3;k++) if (bl->node[i][j][k]!=NULL) if (bl->node[i][j][k]->nodeno==3116) {

  //for (int ii=0;ii<3;ii++) for (int jj=0;jj<3;jj++) for (int kk=0;kk<3;kk++) if (bl->node[ii][jj][kk]!=NULL) if (bl->node[ii][jj][kk]->TEMP_IDENTIFYER_NODE==13635) {

  //for (int ii=0;ii<3;ii++) for (int jj=0;jj<3;jj++) for (int kk=0;kk<3;kk++) if (bl->node[ii][jj][kk]!=NULL) if (bl->node[ii][jj][kk]->TEMP_IDENTIFYER_NODE==10985) {

  //if (bl->node[2][2][2]->TEMP_IDENTIFYER_NODE==10356) {

  cout <<"line=" << __LINE__ << endl;
  }
  */



  ///--------

      static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};
      static const int nodeEdges[8][3]={{0,4,8},{0,5,9},{5,10,1},{4,11,1},{8,3,7},{3,9,6},{2,10,6},{7,11,2}};

      static const int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};

      int globalNodeMask[3][3][3];

      globalNodeMask[0][0][0]=0;
      globalNodeMask[2][0][0]=1;
      globalNodeMask[2][2][0]=2;
      globalNodeMask[0][2][0]=3;

      globalNodeMask[0][0][2]=4;
      globalNodeMask[2][0][2]=5;
      globalNodeMask[2][2][2]=6;
      globalNodeMask[0][2][2]=7;



      //the nodes od the block
      #define nConnectionsMax 20
      #define blNodeListLengthMax  20

      class cBlockNode {
      public:
        int nConnections;
        cBlockNode *Connection[nConnectionsMax];
        typename list<cCutBlockNode>::iterator meshNode;
        int ExtendedNodeMapID;

        bool ghostNode;

        cBlockNode() {nConnections=0,ghostNode=false;}

        void assign(typename list<cCutBlockNode>::iterator nd) {
          meshNode=nd;
          nConnections=0;
        }

        bool checkConnection(cBlockNode* nd) {
          for (int i=0;i<nConnections;i++) if (Connection[i]==nd) return true;

          return false;
        }

        void Connect(cBlockNode* nd) {
          bool errorflag=false;

          if (checkConnection(nd)==false) {
            if (nConnections==nConnectionsMax) errorflag=true;
            else Connection[nConnections++]=nd;
          }

          if (nd->checkConnection(this)==false) {
            if (nd->nConnections==nConnectionsMax) errorflag=true;
            else nd->Connection[nd->nConnections++]=this;
          }

          if (errorflag==true) {
            cout << "ERROR: too many connections (file=" << __FILE__ << ", line=" << __LINE__ << endl;
            ::exit(0);
          }
        }

        void removeConnection(cBlockNode* nd) {
          bool foundFirst=false,foundSecond=false;
          int i;

          for (i=0;i<nConnections;i++) if (Connection[i]==nd) {
            for (++i;i<nConnections;i++) Connection[i-1]=Connection[i];
            nConnections--;
            foundFirst=true;
          }

          for (i=0;i<nd->nConnections;i++) if (nd->Connection[i]==this) {
            for (++i;i<nd->nConnections;i++) nd->Connection[i-1]=nd->Connection[i];
            nd->nConnections--;
            foundSecond=true;
          }

          if ((foundFirst==false)||(foundSecond==false)) {
            cout << "Error: no connections found (file=" << __FILE__ << ", line=" << __LINE__ << endl;
            ::exit(0);
          }
        }
      };


      int blNodeListLength=0;
      cBlockNode blNodeList[blNodeListLengthMax];
      cBlockNode *nodeMesh[3][3][3];

      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) nodeMesh[i][j][k]=NULL;

      //get the cut information for the block
  /*    cCutData *cutData=NULL;

      //first set the cut information with the block data
      //cutData=bl->cutData;

      //if the block doe's not have a cut information, search the edjes of the block
      if (cutData==NULL) for (nedge=0;nedge<12;nedge++) if (bl->edge[nedge]->cutData!=NULL) {
        cutData=bl->edge[nedge]->cutData;
        break;
      }*/


      bool CutBlockFlag=false;

      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) if ((i==1)||(j==1)||(k==1)) if (bl->node[i][j][k]!=bl->NodeBuffer.end()) CutBlockFlag=true;

      //if a cut is not found, return connectivity list based on four nodes of the block
      if (CutBlockFlag==false) {   ////(cutData==NULL) {

  /*
        exit(__LINE__,"Error: this is not a cutted block");

        #if _USE_CELLS_ == _USE_CELLS_ON_
        bl->cell->Volume=bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2];
        #endif
  */

  /*
      for (i=0;i<3;i+=3) for (j=0;j<3;j++) for (k=0;k<3;k++) if (bl->node[i][j][k]!=NULL) {
         if (bl->node[i][j][k]->nodeno==-1) {
                  bl->node[i][j][k]->nodeno=nodenoCounter++;
            elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
                  fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

                  elementPosition=blocks.GetEntryCountingNumber(bl);
                  fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
         }

             }

            blocknoCounter++;

            //check if all nodes are already counted
            for (int ii=0;ii<3;ii+=2) for (int jj=0;jj<3;jj+=2) for (int kk=0;kk<3;kk+=2) if (bl->node[ii][jj][kk]->nodeno<0) exit(__LINE__,"Un-numbered node is found");

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][0]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

            elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][2]);
            fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        return 1;*/

        //the block is not cut -> return the triangulation of an un-cut cube
        static const int nonCutBlockThetrahedronConnectivity[6][4]={{0,1,2,5},{0,2,4,5},{6,2,4,5},     {4,6,7,3},{2,6,4,3},{2,0,4,3}};
        static const int InternalNodeMap[8][3]={{0,0,0},{2,0,0},{2,2,0},{0,2,0},    {0,0,2},{2,0,2},{2,2,2},{0,2,2}};



        for (int nTetrahedron=0;nTetrahedron<6;nTetrahedron++) {
          cTetrahedron tt;

          for (int iNode=0;iNode<4;iNode++) {
            int i=InternalNodeMap[nonCutBlockThetrahedronConnectivity[nTetrahedron][iNode]][0];
            int j=InternalNodeMap[nonCutBlockThetrahedronConnectivity[nTetrahedron][iNode]][1];
            int k=InternalNodeMap[nonCutBlockThetrahedronConnectivity[nTetrahedron][iNode]][2];

            tt.node[iNode]=bl->node[i][j][k];
          }

          indomainConnectivityList.push_back(tt);
        }

        return 6;
      }

      //init the mesh
      //corner nodes of the block
      int nnode,iedge;

      #define _IN_DOMAIN_      0
      #define _OUT_DOMAIN_     1
      #define _DEFAULT_DOMAIN_ 2

  /*    int node_in_domain;

      bool cutBlockFlag=true;


      //a block is considered to be cut if at lease two face are cuted
      bool block_cut_flag=false;
      int nFaceIntersection=0;

      for (int nface=0;nface<6;nface++) {
        int nEdgeIntersection=0;

        for (int pedge=0;pedge<4;pedge++) if (bl->edge[faceEdges[nface][pedge]]->cutPoint!=NULL) nEdgeIntersection++;
        if (nEdgeIntersection>=2) nFaceIntersection++;
      }

      if (nFaceIntersection<2) { //the block cannot be considered to be cut: the regular connectivity list is printing
         std::cout << "ERROR: cutBlockConnectivity - found an uncut block" << std::endl;

        #if _USE_CELLS_ == _USE_CELLS_ON_
        bl->cell->Volume=bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2];
        #endif


        for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) if (bl->node[i][j][k]!=NULL) {
          if (bl->node[i][j][k]->nodeno==-1) {
             bl->node[i][j][k]->nodeno=nodenoCounter++;
             elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
             fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

             elementPosition=blocks.GetEntryCountingNumber(bl);
             fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
          }
        }

        blocknoCounter++;

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][0]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][0][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][0][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[2][2][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        elementPosition=nodes.GetEntryCountingNumber(bl->node[0][2][2]);
        fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

        return 1;

  cutBlockFlag=false;

      }*/




  //NEW PROCEDURE FOR POPULATING THE BLOCK'S NODE MESH
  //  1. the cut node on a edge determines weather the nearby corner nodes exists
  //  2. a) for a cut node on a edge : it is removed is two nearby corner nodes are not exists (the whole edge is removed)
  //  2. b) for a cut node on a edge : it is removed is two nearby corner nodes are  exists  (the edge sgould not be cut)

      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) globalNodeMask[i][j][k]=0;



      //test consistency of the edge cut
      for (nedge=0;nedge<12;nedge++) {
         i=edgeCutNodeCoordinates[nedge][0];
         j=edgeCutNodeCoordinates[nedge][1];
         k=edgeCutNodeCoordinates[nedge][2];

         //check existance of the cornet nodes
         register int pnode;
         int ii,jj,kk;

         //if the cut node does not exists ->
         int cnt=0;

         for (pnode=-1;pnode<2;pnode+=2) {
           if (i==1) ii=i+pnode,jj=j,kk=k;
           if (j==1) ii=i,jj=j+pnode,kk=k;
           if (k==1) ii=i,jj=j,kk=k+pnode;

           if (bl->node[ii][jj][kk]!=bl->NodeBuffer.end()) cnt++;
         }

         if (bl->node[i][j][k]!=bl->NodeBuffer.end()) {
           if ((cnt==0)||(cnt==2)) exit(__LINE__,__FILE__,"Error: a cut-point is found on a non-cut edge");

           globalNodeMask[i][j][k]=1;

           if (blNodeListLength==blNodeListLengthMax) exit(__LINE__,__FILE__,"Error: blNodeListLength reached its maximum value");
           blNodeList[blNodeListLength].assign(bl->node[i][j][k]);
           nodeMesh[i][j][k]=blNodeList+blNodeListLength;
           nodeMesh[i][j][k]->ghostNode=false;
           blNodeListLength++;
         }
         else {
           if (cnt==1) exit(__LINE__,__FILE__,"Error: a cut-point is not found on a cut-edge");

           globalNodeMask[i][j][k]=-1;
         }
      }

  /*

      for (nedge=0;nedge<12;nedge++) {
        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        //check existance of the cornet nodes
        register int pnode;
        int ii,jj,kk;

        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

  //      if ((globalNodeMask[ii][jj][kk]==0)||(globalNodeMask[ii][jj][kk]==1)) {  //the node has not been removed from the block node mesh by this procedure
            if (bl->edge[nedge]->cutPoint!=NULL) {  // if the cut node doe not existes -> noth corner nodes should be on the mesh (if the edge itself on the mesh)
              //check if the corner nodes sould be on the mesh

              x[0]=bl->xmin[0]+bl->dxBlock[0]*ii/2.0;
              x[1]=bl->xmin[1]+bl->dxBlock[1]*jj/2.0;
              x[2]=bl->xmin[2]+bl->dxBlock[2]*kk/2.0;

              cutData=bl->edge[nedge]->cutData;

              for (c=0.0,idim=0;idim<3;idim++) c+=(x[idim]-cutData->x0[idim])*cutData->norm[idim];

                if (c==0.0) continue;

  //            if (cutBlockFlag==true) {



                if (globalNodeMask[ii][jj][kk]==0) globalNodeMask[ii][jj][kk]=(c>1.0E-6) ? -1 : 1;
                else {
                  int code=(c>1.0E-6) ? -1 : 1;

                  if (code!=globalNodeMask[ii][jj][kk]) {
                    code=1;

                    for (long int nSurfaceElement=0;nSurfaceElement<nSurfaceNASTRANelements;nSurfaceElement++) {
                      code=(surfaceNASTRANmesh[nSurfaceElement].checkPointInsideDomain(x,EPS)==true) ? 1 : -1;
                      if (code==-1) break;
                    }

                    globalNodeMask[ii][jj][kk]=code;
                  }

                }


  //            } else globalNodeMask[ii][jj][kk]=1;
            }
  //        }
        }
      }





  //#################  EXPLIVCITLY CHECK THE CORNER AND CUT NODES ################
  {
     //check the corner nodes of the block





  if (bl->TEMP_IDENTIFYER_BLOCK==22469) {
  cout << __LINE__ << endl;
  }

      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) { ///if (globalNodeMask[i][j][k]>=0) {

        globalNodeMask[i][j][k] = (checkPointInsideDomain(bl->node[i][j][k]->x)==true) ? 1 : -1;
      }






      //check the cut nodes

     int cnt,ii,jj,kk,pnode;

     for (nedge=0;nedge<12;nedge++) {
        cnt=0;

        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        //check existance of the cornet nodes
        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

          if (globalNodeMask[ii][jj][kk]==1) cnt++;
        }

        if (cnt==1) {

          //the cut should be here

          //check if the cut is steed up
          if (bl->edge[nedge]->cutData!=NULL) {
            double *n=bl->edge[nedge]->cutData->norm;

            if (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]>0.00001) continue;
          }


          //set the buffers
          if (bl->edge[nedge]->cutPoint==NULL) {
            bl->edge[nedge]->cutPoint=nodes.newElement();
          }

          if (bl->edge[nedge]->cutData==NULL) {
            bl->edge[nedge]->cutData=cutGeometryData.newElement();
          }


          //init the cut data
          double l;

          for (l=0.0,idim=0;idim<3;idim++) {
            bl->edge[nedge]->cutPoint->x[idim]=0.5*(bl->edge[nedge]->node[0]->x[idim]+bl->edge[nedge]->node[1]->x[idim]);
            bl->edge[nedge]->cutPoint->nodeat=bl->node[ii][jj][kk]->nodeat;

            bl->edge[nedge]->cutData->x0[idim]=bl->edge[nedge]->cutPoint->x[idim];
            bl->edge[nedge]->cutData->norm[idim]=bl->edge[nedge]->node[1]->x[idim]-bl->edge[nedge]->node[0]->x[idim];
            l+=pow(bl->edge[nedge]->cutData->norm[idim],2);
          }

          l=sqrt(l);

          if (globalNodeMask[ii][jj][kk]==1) {
            if (bl->node[ii][jj][kk]==bl->edge[nedge]->node[1]) l*=-1.0;
          }
          else {
            if (bl->node[ii][jj][kk]==bl->edge[nedge]->node[0]) l*=-1.0;
          }

          for (idim=0;idim<3;idim++) bl->edge[nedge]->cutData->norm[idim]/=l;
        }
        else {

          //the cut shpild not be here

          //check if the cut is setted up
          if (bl->edge[nedge]->cutData!=NULL) {
            double *n=bl->edge[nedge]->cutData->norm;

  //################  DEBUG #################

  if (bl->edge[nedge]->cutPoint->TEMP_IDENTIFYER_NODE==71215) {
  cout << __LINE__ << endl;
  }

  //################ END DEBUG ##############



            if (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]>0.00001) n[0]=0.0,n[1]=0.0,n[2]=0.0;
          }

        }
      }

  }
  //################ END: EXPLIVCITLY CHECK THE CORNER AND CUT NODES ################



      //set up the cut nodes on edges: a cut node can exists when only one corner node of the edge exists
      for (nedge=0;nedge<12;nedge++) {
        i=edgeCutNodeCoordinates[nedge][0];
        j=edgeCutNodeCoordinates[nedge][1];
        k=edgeCutNodeCoordinates[nedge][2];

        globalNodeMask[i][j][k]=-1;


        //the cut node exists if
        //1. bl->edge[nedge]->cutData!=NULL
        //2. the length of bl->edge[nedge]->cutData->norm > 0.0
        //3. only one of nodes that defines the edge if the block exists

        if (bl->edge[nedge]->cutData==NULL) continue;

        double *cutNorm=bl->edge[nedge]->cutData->norm;
        if (cutNorm[0]*cutNorm[0]+cutNorm[1]*cutNorm[1]+cutNorm[2]*cutNorm[2]<0.01) continue;



        //check existance of the cornet nodes
        register int pnode;
        int ii,jj,kk;

        int nCornerNodes=0;

        //if the cut node does not exists ->
        for (pnode=-1;pnode<2;pnode+=2) {
          if (i==1) ii=i+pnode,jj=j,kk=k;
          if (j==1) ii=i,jj=j+pnode,kk=k;
          if (k==1) ii=i,jj=j,kk=k+pnode;

          if (globalNodeMask[ii][jj][kk]>=0) nCornerNodes++;
        }







        //add the cut node to the node mesh
        if (bl->edge[nedge]->cutPoint!=NULL) if ( (cutBlockFlag==false) || (nCornerNodes==1)) {   ///////((nCornerNodes==1)&&(bl->edge[nedge]->cutPoint!=NULL)) {
          if (blNodeListLength==blNodeListLengthMax) {
            exit(__LINE__,"Error: the nodeList is too short");
          }

          globalNodeMask[i][j][k]=1;

          blNodeList[blNodeListLength].assign(bl->edge[nedge]->cutPoint);
          nodeMesh[i][j][k]=blNodeList+blNodeListLength;
          nodeMesh[i][j][k]->ghostNode=false;
          blNodeListLength++;

          if (bl->edge[nedge]->cutPoint->nodeno==-1) {
            bl->edge[nedge]->cutPoint->nodeno=nodenoCounter++;
            elementPosition=nodes.GetEntryCountingNumber(bl->edge[nedge]->cutPoint);
            fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

            elementPosition=blocks.GetEntryCountingNumber(bl);
            fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
          }
        }
      }*/

      //add the corner nodes to the nodes mesh
      list<cCutBlockNode> GhostNodeList;


      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) {
        globalNodeMask[i][j][k]=(bl->node[i][j][k]!=bl->NodeBuffer.end()) ? 1 : -1;

        if (globalNodeMask[i][j][k]==-1) {
          cCutBlockNode GhostNode;

           GhostNode.x[0]=bl->xBlockMin[0]+0.5*i*bl->dxBlock[0];
           GhostNode.x[1]=bl->xBlockMin[1]+0.5*j*bl->dxBlock[1];
           GhostNode.x[2]=bl->xBlockMin[2]+0.5*k*bl->dxBlock[2];

           GhostNodeList.push_front(GhostNode);
        }

        if (blNodeListLength==blNodeListLengthMax) exit(__LINE__,__FILE__,"Error: blNodeListLength reached its maximum value");

        blNodeList[blNodeListLength].assign(((globalNodeMask[i][j][k]!=-1) ? bl->node[i][j][k] : GhostNodeList.begin()));

        nodeMesh[i][j][k]=blNodeList+blNodeListLength;
        nodeMesh[i][j][k]->ghostNode=(globalNodeMask[i][j][k]>=0) ? false : true;
        blNodeListLength++;
      }




  /*

      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) { ///if (globalNodeMask[i][j][k]>=0) {
        if (blNodeListLength==blNodeListLengthMax) {
          exit(__LINE__,__FILE__,"Error: the nodeList is too short");
        }




  if (globalNodeMask[i][j][k]==0) {
    globalNodeMask[i][j][k]=1;

    for (long int nSurfaceElement=0;nSurfaceElement<nSurfaceNASTRANelements;nSurfaceElement++) {
      if (surfaceNASTRANmesh[nSurfaceElement].checkPointInsideDomain(bl->node[i][j][k]->x,EPS)==false) {
        globalNodeMask[i][j][k]=-1;
        break;
      }
    }
  }



  //////globalNodeMask[i][j][k]=(surfaceNASTRANmesh.checkPointInsideDomain(bl->node[i][j][k]->x)==true) ? 1 : -1;







        blNodeList[blNodeListLength].assign(bl->node[i][j][k]);
        nodeMesh[i][j][k]=blNodeList+blNodeListLength;


  nodeMesh[i][j][k]->ghostNode=(globalNodeMask[i][j][k]>=0) ? false : true;

        blNodeListLength++;



        if ( (nodeMesh[i][j][k]->ghostNode==false) || (cutBlockFlag==false) ) if (bl->node[i][j][k]->nodeno==-1) {
          bl->node[i][j][k]->nodeno=nodenoCounter++;
          elementPosition=nodes.GetEntryCountingNumber(bl->node[i][j][k]);
          fwrite(&elementPosition,sizeof(long int),1,nodesListFile);

          elementPosition=blocks.GetEntryCountingNumber(bl);
          fwrite(&elementPosition,sizeof(long int),1,nodesListFile);
        }

      }

  */






  // END OF NEW PROCEDURE FOR POPULATING THE BLOCK'S NODE MESH

      //connect nodes:nodes in corners of a block can be connected only along coordinate directions; nodes that cut edges sould be connected to another cutting node on the same face
      cBlockNode *connectNode,*nd;
      int n,ndir;
      bool foundflag;


      //connect the corner nodes
      for (i=0;i<3;i+=2) for (j=0;j<3;j+=2) for (k=0;k<3;k+=2) if (nodeMesh[i][j][k]!=NULL) {
        nd=nodeMesh[i][j][k];

        for (ndir=0;ndir<3;ndir++) {
          connectNode=NULL;

          switch (ndir) {
          case 0:
            //check connection in the i-direction
            connectNode=(nodeMesh[1][j][k]==NULL) ? nodeMesh[(i==0) ? 2 : 0][j][k] : nodeMesh[1][j][k];
            break;
          case 1:
            //check connection in the j-direction
            connectNode=(nodeMesh[i][1][k]==NULL) ? nodeMesh[i][(j==0) ? 2 : 0][k] : nodeMesh[i][1][k];
            break;
          case 2:
            //check connection in the k-direction
            connectNode=(nodeMesh[i][j][1]==NULL) ? nodeMesh[i][j][(k==0) ? 2 : 0] : nodeMesh[i][j][1];
          }

          if (connectNode!=NULL) nd->Connect(connectNode);
        }
      }




      //connect cutting nodes: only those cutting nodes are connected that belongs to the same face
      cBlockNode *nd0,*nd1,*nd2,*nd3,*ndcnt;
      int nCuts,nface,pedge;







      //check connection of nodes that cuts the edges
      for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) if ((i==1)||(j==1)||(k==1)) if (nodeMesh[i][j][k]!=NULL) if (nodeMesh[i][j][k]->nConnections==0) {
        exit(__LINE__,"Error: a cutting node has no connections");
      }


  /*
      //copy the cut nodes on the block's node map
    for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) if ((i==1)||(j==1)||(k==1)) if (nodeMesh[i][j][k]!=NULL) bl->node[i][j][k]=nodeMesh[i][j][k]->meshNode;
  */




      //combine the nodes into thetrahedrones
      int nConnections,ncnt;
      int nnodesNumber=blNodeListLength;

      double *xt0ptr,xt1[3],*xt1ptr,xt2[3],*xt2ptr,xt3[3],*xt3ptr,tVolume;
      cBlockNode *tn0;

      //reset the volume of the cell associated with the block
  /*
      #if _USE_CELLS_ == _USE_CELLS_ON_
      bl->cell->Volume=0.0;
      #endif
  */

      while (nnodesNumber>3) {
        //choose the first point that has 3 connections or a point that has a maximum number of connections


  //choose the gost nodes

      for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if ((blNodeList[n].ghostNode==true)&&(blNodeList[n].nConnections==3)) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }


              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-25*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;

          }

        }


        if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].nConnections==3) {
          //check the volume of the prospective a tetrahedron
          tn0=blNodeList+n;

          xt0ptr=tn0->meshNode->x;
          xt1ptr=tn0->Connection[0]->meshNode->x;
          xt2ptr=tn0->Connection[1]->meshNode->x;
          xt3ptr=tn0->Connection[2]->meshNode->x;

          if ((tn0->Connection[0]->ghostNode==true)||(tn0->Connection[1]->ghostNode==true)||(tn0->Connection[2]->ghostNode==true)) continue;


          for (idim=0;idim<3;idim++) {
            xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
            xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
            xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
          }

          tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
          if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
            nd0=blNodeList+n;
            break;
          }
        }



      if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].ghostNode==true) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }


              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-25*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;

          }

        }




        if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) if (blNodeList[n].nConnections==3) {
          //check the volume of the prospective a tetrahedron
          tn0=blNodeList+n;

          xt0ptr=tn0->meshNode->x;
          xt1ptr=tn0->Connection[0]->meshNode->x;
          xt2ptr=tn0->Connection[1]->meshNode->x;
          xt3ptr=tn0->Connection[2]->meshNode->x;

          for (idim=0;idim<3;idim++) {
            xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
            xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
            xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
          }

          tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
          if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
            nd0=blNodeList+n;
            break;
          }
        }


         if (nd0==NULL) for (n=0,nd0=NULL,nConnections=-1;n<blNodeListLength;n++) {
            int itn1,itn2,itn3,tn0Connections;
            bool zeroVolume=false;

            //check the volume of possible tetraedrons
            tn0=blNodeList+n;
            tn0Connections=tn0->nConnections;

            for (itn1=0;(itn1<tn0Connections)&&(zeroVolume==false);itn1++)  for (itn2=itn1+1;(itn2<tn0Connections)&&(zeroVolume==false);itn2++) for (itn3=itn2+1;(itn3<tn0Connections)&&(zeroVolume==false);itn3++) {
              xt0ptr=tn0->meshNode->x;
              xt1ptr=tn0->Connection[itn1]->meshNode->x;
              xt2ptr=tn0->Connection[itn2]->meshNode->x;
              xt3ptr=tn0->Connection[itn3]->meshNode->x;

              for (idim=0;idim<3;idim++) {
                xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
                xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
                xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
              }

              tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
              if (tVolume>1.0E-15*bl->dxBlock[0]*bl->dxBlock[1]*bl->dxBlock[2]) {
                nd0=blNodeList+n,nConnections=blNodeList[n].nConnections;
              }
              else zeroVolume=true;
            }
          }


        if (nd0==NULL) {
          return connectivityLength;

          exit(__LINE__,"Error: all possible corner nodes has a zero volume tetrahedron among their connections");
        }

        if (nd0->nConnections==3) {
          //if the node nd0 has only 3 connections, combine them into a thetrahedral
          nd1=nd0->Connection[0];
          nd2=nd0->Connection[1];
          nd3=nd0->Connection[2];

          //combine the nodes into a thetrahedron
          cTetrahedron tt;

          tt.node[0]=nd0->meshNode;
          tt.node[1]=nd1->meshNode;
          tt.node[2]=nd2->meshNode;
          tt.node[3]=nd3->meshNode;

          if ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) {
            indomainConnectivityList.push_back(tt);
          }
          else {
            outdomainConnectivityList.push_back(tt);
          }

          //connect the nodes
          nd1->Connect(nd2);
          nd1->Connect(nd3);
          nd2->Connect(nd3);

  /*
          #if _USE_CELLS_ == _USE_CELLS_ON_
          bl->cell->Volume+=tVolume;
          #endif
  */

          //print the connectivity list
        if (true) {   //////((FileOutputFlag==true)&&(printConnectivityList==true)) {

  /*
  /////
   //NBLOCKS++;
   BLOCKPRINTED++;


  if (PRINTEDFLAG==false) PRINTEDFLAG=true,PRINTEDCELLS++;
  ///





  if (  (cutBlockFlag==false) ||  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false))) {

    blocknoCounter++;

          if (nd1->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd1->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd2->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd2->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd3->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd3->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd0->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd0->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);






  ////////

  if ((1+nd0->meshNode->nodeno==0)||(1+nd1->meshNode->nodeno==0)||(1+nd2->meshNode->nodeno==0)||(1+nd3->meshNode->nodeno==0)) {
  cout << __LINE__ << endl;
  }
  //////


  }
  */



            //the surface mesh should contain faces that:
            //1. are a part of a thetrahedrons that contain a ghost node (<-> the face is adjesent to the boundary)
            //2. the face itself should not contain any ghost nodes
            //3. all nodes that compose the cut face should have an attribute that is different from _DEFAULT_NODEAT_


            //prepare data for the surface mesh
  //          cAMRnode<cNode> *thetraNodes[4];

          typename list<cCutBlockNode>::iterator thetraNodes[4];
          int i,j,k;


            bool flag;
            cBlockNode *ndlist[4];



            //the thetrahedral face is a surface face only if each node of the face has the nodeat different from _DEFAULT_NODEAT_
            //the first thetrahedron
            thetraNodes[0]=nd1->meshNode;
            thetraNodes[1]=nd2->meshNode;
            thetraNodes[2]=nd3->meshNode;
            thetraNodes[3]=nd0->meshNode;

            ndlist[0]=nd1;
            ndlist[1]=nd2;
            ndlist[2]=nd3;
            ndlist[3]=nd0;




  if (true) { ///(cutBlockFlag==true) {
  //the block is cut


            for (i=0;i<4;i++)  {
              for (flag=false,j=0;j<4;j++) if (i!=j) if (ndlist[j]->ghostNode==true) {
                flag=true;
                break;
              }

              if (flag==true) continue;


              //if the node ndlist[i] is 'ghostNode' than the face 'ndlist[j]' belongs to the surface
              if (ndlist[i]->ghostNode==true) flag=true;



               //if 'flag==true' -> the surface face is found

               if (flag==true) {

  /*                 cTriangularSurfaceFace newface;
                   cAMRnode<cNode> *upNode;*/

                   cTriangleCutFace newface;
                   typename list<cCutBlockNode>::iterator upNode;


                   for (k=0,j=0;j<4;j++) {
                     if (i!=j) {
                       newface.node[k++]=thetraNodes[j];

  /*                     if (thetraNodes[j]->nodeat!=_DEFAULT_NODEAT_) newface.faceat=thetraNodes[j]->nodeat;*/

                     }
                     else upNode=thetraNodes[j];
                   }

                   //get the external normal and the surface area of the cut face
                   double e0[3],e1[3],l,*x0,*x1,*x2,norm[3],c;

                   x0=newface.node[0]->x;
                   x1=newface.node[1]->x;
                   x2=newface.node[2]->x;

                   for (idim=0;idim<3;idim++) e0[idim]=x1[idim]-x0[idim],e1[idim]=x2[idim]-x0[idim];



                   norm[0]=e0[1]*e1[2]-e1[1]*e0[2];
                   norm[1]=e1[0]*e0[2]-e0[0]*e1[2];
                   norm[2]=e0[0]*e1[1]-e1[0]*e0[1];

                   l=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
                   newface.SurfaceArea=0.5*l;

                   for (idim=0,c=0.0;idim<3;idim++) c+=norm[idim]*(upNode->x[idim]-x0[idim]);

                   if (c>0.0) l*=-1.0;

                   for (idim=0,c=0.0;idim<3;idim++) newface.ExternalNormal[idim]=norm[idim]/l;



                   //add the cut surface to the surface mesh
  /*                 newface.nextSurfaceCutFace=bl->firstSurfaceCutFace;
                   newface.block=bl;
                   bl->firstSurfaceCutFace=SurfaceMeshData.nfaces;

                   SurfaceMeshData.faces.push_back(newface);
                   SurfaceMeshData.nfaces++;*/

                   TriangleCutFaceConnectivity.push_back(newface);
                }
            }

  }
  else {
    //the block is not cut

  exit(__LINE__,"ERROR: the block is not cut - non finished section");


    //the surface face: 1. must contan a node that is an cutting node of an edge
    //2. must contain a node that is belongs to the computational domain

  /*  int nface,pedge,nedge,nCorner,i[3],searchDirection,cnt,nd;
    cCutBlockNode *firstNode,*secondNode,*thirdNode,*forthNode;

    static const int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};*/




  }



          }

          //inicrement the connectivity counter
  if (  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) )      ++connectivityLength;

          //disconnect the point nd0
          nd1->removeConnection(nd0);
          nd2->removeConnection(nd0);
          nd3->removeConnection(nd0);

          //udpate the number of left nodes
          if (nd0->nConnections==0) nnodesNumber--;
          if (nd1->nConnections==0) nnodesNumber--;
          if (nd2->nConnections==0) nnodesNumber--;
          if (nd3->nConnections==0) nnodesNumber--;
        }
        else if (nd0->nConnections>3) {
          //create a plane and sort the nodes in the plane
          //the point of the origin of the plane, normal to the plane aand coordinate vectors wrelated to the plane
          double x0[3]={0.0,0.0,0.0},norm[3],e0[3],e1[3];
          double *x,length;
          cBlockNode **nd0Connection=nd0->Connection;

          for (n=0;n<nConnections;n++) for (x=nd0Connection[n]->meshNode->x,idim=0;idim<3;idim++) x0[idim]+=x[idim];

          for (idim=0,length=0.0;idim<3;idim++) {
            x0[idim]/=nConnections;

            e0[idim]=nd0Connection[1]->meshNode->x[idim]-nd0Connection[0]->meshNode->x[idim];
            e1[idim]=nd0Connection[2]->meshNode->x[idim]-nd0Connection[0]->meshNode->x[idim];
          }


          norm[0]=e1[1]*e0[2]-e1[2]*e0[1];
          norm[1]=e1[2]*e0[0]-e1[0]*e0[2];
          norm[2]=e1[0]*e0[1]-e1[1]*e0[0];

          length=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);




          for (idim=0;idim<3;idim++) norm[idim]/=length;

          //construct a coordinate system in the plane: get the first coordinate vector
          if (fabs(norm[0])>1.0E-1) {
            register double t=sqrt(norm[0]*norm[0]+norm[1]*norm[1]);

            e0[0]=norm[1]/t,e0[1]=-norm[0]/t,e0[2]=0.0;
          }
          else {
            register double t=sqrt(norm[2]*norm[2]+norm[1]*norm[1]);

            e0[0]=0.0,e0[1]=-norm[2]/t,e0[2]=norm[1]/t;
          }

          //get the second coordinate vector
          e1[0]=norm[1]*e0[2]-norm[2]*e0[1];
          e1[1]=norm[2]*e0[0]-norm[0]*e0[2];
          e1[2]=norm[0]*e0[1]-norm[1]*e0[0];

          //get the angles of projections of the nodes in the plane (e0,e1,norm)
          double phi,xplane[2];
          double *phi_nd0Connection=new double[nConnections];

          for (n=0;n<nConnections;n++) {
            x=nd0Connection[n]->meshNode->x;

            for (idim=0,xplane[0]=0.0,xplane[1]=0.0;idim<3;idim++) xplane[0]+=(x[idim]-x0[idim])*e0[idim],xplane[1]+=(x[idim]-x0[idim])*e1[idim];

            phi_nd0Connection[n]=acos(xplane[0]/sqrt(xplane[0]*xplane[0]+xplane[1]*xplane[1]));
            if (xplane[1]<0.0) phi_nd0Connection[n]=2.0*Pi-phi_nd0Connection[n];
          }

          //sort the nodes with an increase of the angle
          int nmin,n1;
          double minphi;

          for (n=0;n<nConnections;n++) {
            for (n1=n+1,nmin=n,minphi=phi_nd0Connection[n];n1<nConnections;n1++) if (minphi>phi_nd0Connection[n1]) nmin=n1,minphi=phi_nd0Connection[n1];

            if (nmin!=n) {
              //swap the nodes in the list
              register cBlockNode *t;

              t=nd0Connection[n];
              nd0Connection[n]=nd0Connection[nmin];
              nd0Connection[nmin]=t;

              phi_nd0Connection[nmin]=phi_nd0Connection[n];
            }
          }

          delete [] phi_nd0Connection;

          //construct the set of thetrahedrons
          int nthetra;

          for (nd1=nd0Connection[0],nthetra=0;nthetra<nConnections-2;nthetra++) {
            nd2=nd0Connection[nthetra+1];
            nd3=nd0Connection[nthetra+2];

            nd1->Connect(nd2);
            nd1->Connect(nd3);
            nd2->Connect(nd3);

            //combine the nodes into a thetrahedron
            cTetrahedron tt;

            tt.node[0]=nd0->meshNode;
            tt.node[1]=nd1->meshNode;
            tt.node[2]=nd2->meshNode;
            tt.node[3]=nd3->meshNode;

            if ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) {
              indomainConnectivityList.push_back(tt);
            }
            else {
              outdomainConnectivityList.push_back(tt);
            }


            //increment the volume of the cell assiciated with the block
            #if _USE_CELLS_ == _USE_CELLS_ON_
            xt0ptr=nd0->meshNode->x;
            xt1ptr=nd1->meshNode->x;
            xt2ptr=nd2->meshNode->x;
            xt3ptr=nd3->meshNode->x;

            for (idim=0;idim<3;idim++) {
              xt1[idim]=xt1ptr[idim]-xt0ptr[idim];
              xt2[idim]=xt2ptr[idim]-xt0ptr[idim];
              xt3[idim]=xt3ptr[idim]-xt0ptr[idim];
            }

            tVolume=fabs(xt1[0]*(xt2[1]*xt3[2]-xt2[2]*xt3[1])-xt1[1]*(xt2[0]*xt3[2]-xt2[2]*xt3[0])+xt1[2]*(xt2[0]*xt3[1]-xt2[1]*xt3[0]))/6.0;
  //          bl->cell->Volume+=tVolume;
            #endif



            //print the connectivity list
        if (true) {  /////((FileOutputFlag==true)&&(printConnectivityList==true)) {



  /*
  if (  (cutBlockFlag==false) ||  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) ) {


    blocknoCounter++;

          if (nd1->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd1->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd2->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd2->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd3->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd3->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);

          if (nd0->meshNode->nodeno<0) exit(__LINE__,"Un-numbered node is found");
          elementPosition=nodes.GetEntryCountingNumber(nd0->meshNode);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);
          fwrite(&elementPosition,sizeof(long int),1,connectivityListFile);






  ////////

  if ((1+nd0->meshNode->nodeno==0)||(1+nd1->meshNode->nodeno==0)||(1+nd2->meshNode->nodeno==0)||(1+nd3->meshNode->nodeno==0)) {
  cout << __LINE__ << endl;
  }
  //////



  }
  */

              //prepare data for the surface mesh generation
          typename list<cCutBlockNode>::iterator thetraNodes[4];
              int i,j,k;

            bool flag;
            cBlockNode *ndlist[4];

              //the thetrahedral face is a surface face only if each node of the face has the nodeat different from _DEFAULT_NODEAT_
              //the first thetrahedron
              thetraNodes[0]=nd1->meshNode;
              thetraNodes[1]=nd2->meshNode;
              thetraNodes[2]=nd3->meshNode;
              thetraNodes[3]=nd0->meshNode;

            ndlist[0]=nd1;
            ndlist[1]=nd2;
            ndlist[2]=nd3;
            ndlist[3]=nd0;


  if (true) { //(cutBlockFlag==true) {
  //the block is cut


            for (i=0;i<4;i++)  {
              for (flag=false,j=0;j<4;j++) if (i!=j) if (ndlist[j]->ghostNode==true) {
                flag=true;
                break;
              }

              if (flag==true) continue;


              //if the node ndlist[i] is 'ghostNode' than the face 'ndlist[j]' belongs to the surface
              if (ndlist[i]->ghostNode==true) flag=true;


               if (flag==true) {

                 cTriangleCutFace newface;
                 typename list<cCutBlockNode>::iterator upNode;


                   for (k=0,j=0;j<4;j++) {
                     if (i!=j) {
                       newface.node[k++]=thetraNodes[j];
  //                     if (thetraNodes[j]->nodeat!=_DEFAULT_NODEAT_) newface.faceat=thetraNodes[j]->nodeat;
                     }
                     else upNode=thetraNodes[j];
                   }


                   //get the external normal and the surface area of the cut face
                   double e0[3],e1[3],l,*x0,*x1,*x2,norm[3],c;

                   x0=newface.node[0]->x;
                   x1=newface.node[1]->x;
                   x2=newface.node[2]->x;

                   for (idim=0;idim<3;idim++) e0[idim]=x1[idim]-x0[idim],e1[idim]=x2[idim]-x0[idim];

                   norm[0]=e0[1]*e1[2]-e1[1]*e0[2];
                   norm[1]=e1[0]*e0[2]-e0[0]*e1[2];
                   norm[2]=e0[0]*e1[1]-e1[0]*e0[1];


                   l=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
                   newface.SurfaceArea=0.5*l;

                   for (idim=0,c=0.0;idim<3;idim++) c+=norm[idim]*(upNode->x[idim]-x0[idim]);

                   if (c>0.0) l*=-1.0;

                   for (idim=0,c=0.0;idim<3;idim++) newface.ExternalNormal[idim]=norm[idim]/l;


                    //add the cut surface to the surface mesh
  /*                  newface.nextSurfaceCutFace=bl->firstSurfaceCutFace;
                    newface.block=bl;
                    bl->firstSurfaceCutFace=SurfaceMeshData.nfaces;

                    SurfaceMeshData.faces.push_back(newface);
                    SurfaceMeshData.nfaces++;*/

                   TriangleCutFaceConnectivity.push_back(newface);
                }


              }

  }
  else {
    //the block is not cut


  exit(__LINE__,"ERROR: the block is not cut - non finished section");

    //the surface face: 1. must contan a node that is an cutting node of an edge
    //2. must contain a node that is belongs to the computational domain

  /*  int nface,pedge,nedge,nCorner,i[3],searchDirection,cnt,nd;
    cCutBlockNode *firstNode,*secondNode,*thirdNode,*forthNode;

    static const int edgeCutNodeCoordinates[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},   {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};*/



  }





            }

            //inicrement the connectivity counter
  if (  ((nd0->ghostNode==false)&&(nd1->ghostNode==false)&&(nd2->ghostNode==false)&&(nd3->ghostNode==false)) )        ++connectivityLength;



          }

          //disconnect the point nd0
          for (n=0;n<nConnections;n++) {
            nd1=nd0->Connection[0];
            nd0->removeConnection(nd1);
            if (nd1->nConnections==0) nnodesNumber--;
          }

          if (nd0->nConnections==0) nnodesNumber--;
        }
        else {
         //if the node nd0 has less than tree connections -> error
         exit(__LINE__,"Error: the maximum number of connections is less that three");
       }

      }






  //////
  /*
  if (PRINTEDCELLS==252) {
  cout << "printed bloks: "<< NBLOCKS << endl;

  //fclose(fout);
  //exit(0);
  }
  */
  /////

      return connectivityLength;
    }



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
