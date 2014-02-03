


//the particle class
#include "pic.h"
#include "constants.h"

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


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"


//$Id$


void amps_init();
void amps_time_step();





#include "meshAMRinternalSurface_RotationBody.h"


bool Radius(double &r,double x);


/*bool Radius(double &r,double x) {


 // r=1.0;
 // return true;


  if ((x>=0.0-1.0E-15)&&(x<=Pi+1.0E-15)) {
    r=sin(x);
    return true;
  }

  return false;
}*/


//the cut-block classes
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

class cQuadrilateral {
public:
  list<cCutBlockNode>::iterator node[8];
};


/*class cTriangleFace {
public:
//  list<cCutBlockNode>::iterator node[3];
  double ExternalNormal[3],SurfaceArea;

  double x0Face[3],x1Face[3],x2Face[3];
  double e0[3],e1[3],c00,c01,c11,c,e0Length,e1Length;

  cTriangleFace *next,*prev;

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
  }


  bool RayIntersection(double *x0,double *l,double &t,double EPS) {
    double length,lNorm;

    lNorm=l[0]*ExternalNormal[0]+l[1]*ExternalNormal[1]+l[2]*ExternalNormal[2];
    length=(x0[0]-x0Face[0])*ExternalNormal[0]+(x0[1]-x0Face[1])*ExternalNormal[1]+(x0[2]-x0Face[2])*ExternalNormal[2];

    if (fabs(length)<EPS) {
      t=0.0;
    }
    else {
      t=-length/lNorm;
      if (t<0.0) return false;
    }

    //find position of the intersection in the internal frame related to the face
    double c0=0.0,c1=0.0,xLocal[2];

    for (int i=0;i<3;i++) {
      double x=x0[i]+l[i]*t;

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

    l[0]=x1[0]-x0[0],l[1]=x1[1]-x0[1],l[2]=x1[2]-x0[2];
    RayIntersection(x0,l,t,EPS);
    return ((0<=t)&&(t<=1.0)) ? true : false;
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
    static const double x0Block[8][3]={{0,0,0},{1,0,0},   {0,0,0},{0,1,0},   {0,0,0},{0,0,1}};
    static const double normBlock[8][3]={{1,0,0},{1,0,0},   {0,1,0},{0,1,0},    {0,0,1},{0,0,1}};
//    static const double e0Block[8][3]={{0,1,0},{1,1,0},    {1,0,0}, {1,1,0},   {1,0,0},{1,0,1}};
//    static const double e1Block[8][3]={{0,0,1},{1,0,1},     {0,0,1},{0,1,1},   {0,1,0},{0,1,1}};
    static const int skipBlockFaceDirection[8]={0,0,1,1,2,2};

    double x0Edge[3],lEdge[3],lNorm,xNorm;

    for (pface=0;pface<8;pface++) for (int pFaceEdge=0;pFaceEdge<3;pFaceEdge++) {
      lNorm=0.0,xNorm=0.0;

      switch (pFaceEdge) {
      case 0:
        memcpy(x0Edge,x1Face,3*sizeof(double));
        for (int i=0;i<3;i++) lEdge[i]=x2Face[i]-x1Face[i];
        break;
      case 1:
        memcpy(x0Edge,x0Face,3*sizeof(double));
        for (int i=0;i<3;i++) lEdge[i]=x0Face[i]-x2Face[i];
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
    SurfaceArea=0.0;
    for (int i=0;i<3;i++) ExternalNormal[i]=0.0;
  }
};

*/

class cTriangleCutFace : public cTriangleFace {
public:
  list<cCutBlockNode>::iterator node[3];
};



double BulletLocalResolution(double *x) {
  int idim;


  //check weather the point is inside the surface
  if (CheckPointInsideDomain(x,BoundaryTriangleFaces,nBoundaryTriangleFaces,1.0E-8)==false) {
    return 100.0;
  }

  //the point is within the domain
  //get the distrance of the point from the point to the surface
  double lOrigin[3],l,xSurfaceIntersection[3],tSurfaceIntersection;
  cTriangleFace* FaceIntersection;

  for (idim=0;idim<3;idim++) lOrigin[idim]=-x[idim];

  GetClosestSurfaceIntersectionPoint(x,lOrigin,xSurfaceIntersection,tSurfaceIntersection,FaceIntersection,BoundaryTriangleFaces,nBoundaryTriangleFaces,1.0E-8);
  for (l=0.0,idim=0;idim<3;idim++) l+=pow(x[idim]-xSurfaceIntersection[idim],2);

  l=sqrt(l);


  return (l<10) ? max(1.0,l) : 100;
}


int main(int argc,char **argv) {


  PIC::InitMPI();

/*
  double x0[3]={0.0,0.0,0.0};
  double l[3]={1.0,0.0,0.0};
  cInternalRotationBodyData Nucleus;

  Nucleus.SetGeometricalParameters(x0,l,0.0,Pi,Radius,50,50);
  Nucleus.PrintSurfaceMesh("surface.dat");

//  Nucleus.GetSurfaceTriangulation(BoundaryTriangleFaces,nBoundaryTriangleFaces);
//  PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",BoundaryTriangleFaces,nBoundaryTriangleFaces,1.0E-8);

  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  ReadNastranSurfaceMeshLongFormat("BulletSurface.nas",BoundaryTriangleFaces,nBoundaryTriangleFaces,xmin,xmax,1.0E-8);
  if (PIC::ThisThread==0) PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",BoundaryTriangleFaces,nBoundaryTriangleFaces,1.0E-8);


  for (int i=0;i<3;i++) xmin[i]*=4.0,xmax[i]*=4.0;

  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.buildMesh();
  if (PIC::ThisThread==0) PIC::Mesh::mesh.outputMeshTECPLOT("VolumeMesh.dat");



  int nTetrahedron=0;


  cCutBlockNode nodes[27];
  cCutData CutData[12];
  cCutEdge edges[12];

  cCutBlock bl;
  list<cTetrahedron> indomainConnectivityList,outdomainConnectivityList;
  list<cTriangleCutFace> TriangleCutConnectivity;

  cCutBlockSet<cCutBlockNode,cTetrahedron,cQuadrilateral> CutBlockSet;

  bl.Reset();
  for (int i=0;i<2;i+=1) for (int j=0;j<2;j+=1) for (int k=0;k<2;k+=1) {
    double x[3];

    x[0]=xmin[0]+i*(xmax[0]-xmin[0]);
    x[1]=xmin[1]+j*(xmax[1]-xmin[1]);
    x[2]=xmin[2]+k*(xmax[2]-xmin[2]);


    bl.AddNode(x,2*i,2*j,2*k);
  }


  nTetrahedron=cutBlockTetrahedronConnectivity<cCutBlockNode,cCutBlock,cCutData,cTetrahedron>(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);


  CutBlockSet.AddTetrahedronList(indomainConnectivityList,1.0E-5);

  double TetrahedronVolume=0.0;

  for (list<cTetrahedron>::iterator itr=indomainConnectivityList.begin();itr!=indomainConnectivityList.end();itr++) {
    TetrahedronVolume+=itr->Volume();
  }





  for (int ii=0;ii<10;ii++) {
    int ntest,cnt;
    int IntersectionStatus,i,nAllTest=5000,nMaxAllTest=200000;
    double r,x[3],locx,volume=0,v0=-1.0;
    int IntersectionFlag;

    do {
      for (ntest=0,cnt=0,IntersectionFlag=false;ntest<nAllTest;ntest++) {
        for (i=0,locx=0.0;i<3;i++) {
          x[i]=xmin[i]+rnd()*(xmax[i]-xmin[i]);
          locx+=Nucleus.e0[i]*(x[i]-Nucleus.OriginPosition[i]);
        }

        if (Nucleus.SurfaceCurve(r,locx)) {
          if (r*r>=x[1]*x[1]+x[2]*x[2]) ++cnt,IntersectionFlag=true;
        }
      }

      volume=double(cnt)/double(nAllTest);
      if (cnt!=0) if (fabs((volume-v0)/volume)<0.01) break;
      if (ntest==nMaxAllTest) break;

      v0=volume;

      nAllTest*=4;
      if (nAllTest>nMaxAllTest) nAllTest=nMaxAllTest;
    }
    while (true);


    //triandulation: begin
    //combine the cut-block
    bl.Reset();
    indomainConnectivityList.clear();
    outdomainConnectivityList.clear();
    TriangleCutConnectivity.clear();

    //nodes
    for (int i=0;i<2;i+=1) for (int j=0;j<2;j+=1) for (int k=0;k<2;k+=1) {
      double x[3],locx;
      int ii;

      x[0]=xmin[0]+i*(xmax[0]-xmin[0]);
      x[1]=xmin[1]+j*(xmax[1]-xmin[1]);
      x[2]=xmin[2]+k*(xmax[2]-xmin[2]);

      for (ii=0,locx=0.0;ii<3;ii++) locx+=Nucleus.e0[ii]*(x[ii]-Nucleus.OriginPosition[ii]);

      if (Nucleus.SurfaceCurve(r,locx)) {
        if (r*r<=x[1]*x[1]+x[2]*x[2]) bl.AddNode(x,2*i,2*j,2*k);
      }
    }

    memcpy(bl.xBlockMin,xmin,3*sizeof(double));
    memcpy(bl.xBlockMax,xmax,3*sizeof(double));
    for (i=0;i<3;i++) bl.dxBlock[i]=xmax[i]-xmin[i];

    //edges
    static const int EdgeMiddleNodeMap[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},  {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int x0EdgeMap[12][3]={{0,0,0},{0,2,0},{0,2,2},{0,0,2},  {0,0,0},{2,0,0},{2,0,2},{0,0,2},  {0,0,0},{2,0,0},{2,2,0},{0,2,0}};
    static const int x1EdgeMap[12][3]={{2,0,0},{2,2,0},{2,2,2},{2,0,2},  {0,2,0},{2,2,0},{2,2,2},{0,2,2},  {0,0,2},{2,0,2},{2,2,2},{0,2,2}};

    for (int nedge=0;nedge<12;nedge++) {
      double x0[3],x1[3],xIntersection[3];


      for (int iii=0;iii<3;iii++) {
        x0[iii]=xmin[iii]+0.5*x0EdgeMap[nedge][iii]*(xmax[iii]-xmin[iii]);
        x1[iii]=xmin[iii]+0.5*x1EdgeMap[nedge][iii]*(xmax[iii]-xmin[iii]);
      }


      if (Nucleus.LinearSegment_Surface_Intersection(x0,x1,xIntersection)==true) {
        int cnt=0;

        if (bl.node[x0EdgeMap[nedge][0]][x0EdgeMap[nedge][1]][x0EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
        if (bl.node[x1EdgeMap[nedge][0]][x1EdgeMap[nedge][1]][x1EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
        if ((cnt==0)||(cnt==2)) continue;


        bl.AddNode(xIntersection,EdgeMiddleNodeMap[nedge][0],EdgeMiddleNodeMap[nedge][1],EdgeMiddleNodeMap[nedge][2]);
      }
    }

    nTetrahedron=cutBlockTetrahedronConnectivity<cCutBlockNode,cCutBlock,cCutData,cTetrahedron>(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);



    CutBlockSet.clear();
    CutBlockSet.AddTetrahedronList(indomainConnectivityList,1.0E-5);
    CutBlockSet.PrintTECPLOT("mesh.dat");


    double ConnectivityListVolume,VolumeFraction;
    list<cTetrahedron>* ConnectivityListArray[]={&indomainConnectivityList};


    if (CheckConnectivityListIntersection<cCutBlockNode,cTetrahedron>(bl.xBlockMin,bl.xBlockMax,ConnectivityListArray,1,ConnectivityListVolume,VolumeFraction)==true) {
      exit(__LINE__,__FILE__);
    }

    cutBlockPrintVolumeMesh<cCutBlockNode,cTetrahedron>("inCutBlock.dat",indomainConnectivityList);
//    cutBlockPrintVolumeMesh<cCutBlockNode,cTetrahedron>("outCutBlock.dat",outdomainConnectivityList);


    double TetrahedronVolume=0.0;

    for (list<cTetrahedron>::iterator itr=indomainConnectivityList.begin();itr!=indomainConnectivityList.end();itr++) {
      TetrahedronVolume+=itr->Volume();
    }
//triangulation: end


    double meshvolume=0.0; //Nucleus.GetRemainedBlockVolume(xmin,xmax,1e-10,&IntersectionStatus);
    double newMCvolume=Nucleus.GetRemainedBlockVolumeMC(xmin,xmax,1e-10,1.0E-4,&IntersectionStatus,50000,5);
    double newCutBlockVolume=Nucleus.GetRemainedBlockVolumeCutBlock<cTetrahedron,cQuadrilateral,cCutBlock,cTriangleCutFace,cCutBlockNode,cCutData,cCutBlockSet<cCutBlockNode,cTetrahedron,cQuadrilateral> >(xmin,xmax,1e-10,1.0E-4,&IntersectionStatus,&CutBlockSet,5);

    CutBlockSet.PrintTECPLOT("mesh.dat");

    IntersectionStatus=Nucleus.BlockIntersection(xmin,xmax,1.0E-10);

    std::cout << "Volume: MC=" << (1.0-volume)*(xmax[0]-xmin[0])*(xmax[1]-xmin[1])*(xmax[2]-xmin[2]) << ", new MCvolume="  << newMCvolume << ", newCutCells=" << newCutBlockVolume << ", mesh=" << meshvolume << std::endl;
    std::cout << "Intersection:  MC=" << IntersectionFlag <<", mesh=" << ((IntersectionStatus==_AMR_BLOCK_INSIDE_DOMAIN_) ? 0 : 1) << std::endl;
    std::cout << "Triangulation Volume = " << TetrahedronVolume << std::endl;


    //theoretical value
    double theta,S=0.0;

    if (xmin[2]<=1.0) {
      theta=acos(xmin[2]);
      S=Pi*theta/Pi-xmin[2]*xmin[2]*tan(theta);
    }

    std::cout << "Theoretical value=" << 2.0-S << std::endl;

    double dz=0.1*(xmax[2]-xmin[2]);
    xmax[2]-=dz;
    xmin[2]-=dz;
  }

*/




  //================================

  amps_init();


  for (long int niter=0;niter<100000001 /*000001*/;niter++) {

    amps_time_step();
  }


  char fname[400];

  sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
