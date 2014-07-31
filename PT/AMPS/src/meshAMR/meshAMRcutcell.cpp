//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

#include "constants.h"
#include "specfunc.h"
#include "ifileopr.h"
#include "meshAMRcutcell.h"

CutCell::cTriangleFace *CutCell::BoundaryTriangleFaces=NULL;
int CutCell::nBoundaryTriangleFaces=0;

CutCell::cNASTRANnode *CutCell::BoundaryTriangleNodes=NULL;
int CutCell::nBoundaryTriangleNodes=0;

cAMRstack<CutCell::cTriangleFaceDescriptor> CutCell::BoundaryTriangleFaceDescriptor;

/*struct cNodeCoordinates {
  double *x;
  int id;
};

struct cFaceNodeConnection {
  list<cNodeCoordinates>::iterator node[3];
};*/

void CutCell::PrintSurfaceTriangulationMesh(const char *fname,CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
  int nface,pnode,cnt;
  bool flag;
  double *xNode,*xFace;

  list<cNodeCoordinates> nodeCoordinates;
  list<cNodeCoordinates>::iterator nodeitr;
  cFaceNodeConnection *FaceNodeConnection=new cFaceNodeConnection[nSurfaceTriangulation];

  //reconstruct the node list
  for (nface=0;nface<nSurfaceTriangulation;nface++) for (pnode=0;pnode<3;pnode++) {
    flag=false;

    switch (pnode) {
    case 0:
      xFace=(SurfaceTriangulation+nface)->x0Face;
      break;
    case 1:
      xFace=(SurfaceTriangulation+nface)->x1Face;
      break;
    case 2:
      xFace=(SurfaceTriangulation+nface)->x2Face;
      break;
    }

    for (nodeitr=nodeCoordinates.begin();nodeitr!=nodeCoordinates.end();nodeitr++) {
      xNode=nodeitr->x;

      if (pow(xFace[0]-xNode[0],2)+pow(xFace[1]-xNode[1],2)+pow(xFace[2]-xNode[2],2)<EPS*EPS) {
        flag=true;
        FaceNodeConnection[nface].node[pnode]=nodeitr;
        break;
      }
    }

    if (flag==false) {
      //the node is not found -> create new node
      cNodeCoordinates nd;

      nd.id=0;
      nd.x=xFace;
      nd.pic__shadow_attribute=SurfaceTriangulation[nface].pic__shadow_attribute;

      nd.nface=nface;

      nodeCoordinates.push_front(nd);

      FaceNodeConnection[nface].node[pnode]=nodeCoordinates.begin();
    }
  }


  //print the mesh
  FILE *fout=fopen(fname,"w");
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\",\"Surface shadow attribute\", \"nface\"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",(int)nodeCoordinates.size(),nSurfaceTriangulation);


  for (cnt=1,nodeitr=nodeCoordinates.begin();nodeitr!=nodeCoordinates.end();nodeitr++) {
    nodeitr->id=cnt++;
    fprintf(fout,"%e %e %e %i %i\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2],nodeitr->pic__shadow_attribute,nodeitr->nface);
  }

  for (nface=0;nface<nSurfaceTriangulation;nface++) {
    fprintf(fout,"%i %i %i\n",FaceNodeConnection[nface].node[0]->id,FaceNodeConnection[nface].node[1]->id,FaceNodeConnection[nface].node[2]->id);
  }

  fclose(fout);
  delete [] FaceNodeConnection;
}

/*struct cNASTRANnode {
  double x[3];
  int id;
};

struct cNASTRANface {
  int node[3],faceat;
  double externalNormal[3];
};*/

void CutCell::ReadNastranSurfaceMeshLongFormat(const char *fname,double *xSurfaceMin,double *xSurfaceMax,double EPS) {
  CiFileOperations ifile;
  char str[10000],dat[10000],*endptr;
  long int i,j,idim,nnodes=0,nfaces=0;

  if (BoundaryTriangleFaces!=NULL) exit(__LINE__,__FILE__,"Error: redifinition of the surface triangulation array");


  cNASTRANnode node;
  vector<cNASTRANnode> nodes;

  cNASTRANface face;
  vector<cNASTRANface> faces;

  ifile.openfile(fname);
  ifile.GetInputStr(str,sizeof(str));

  //read nodes from the file

  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(dat,str);

  //load the nodes
  while (strcmp("GRID*",dat)==0) {
    ifile.CutInputStr(dat,str);
    node.id=strtol(dat,&endptr,10);

    ifile.CutInputStr(dat,str);
    ifile.CutInputStr(dat,str);
    node.x[0]=strtod(dat,&endptr);

    ifile.CutInputStr(dat,str);
    node.x[1]=strtod(dat,&endptr);

    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);

    ifile.CutInputStr(dat,str);
    node.x[2]=strtod(dat,&endptr);


    if (nnodes==0) for (idim=0;idim<3;idim++) xSurfaceMin[idim]=node.x[idim],xSurfaceMax[idim]=node.x[idim];
    else {
      for (idim=0;idim<3;idim++) {
        if (xSurfaceMin[idim]>node.x[idim]) xSurfaceMin[idim]=node.x[idim];
        if (xSurfaceMax[idim]<node.x[idim]) xSurfaceMax[idim]=node.x[idim];
      }
    }

    nodes.push_back(node);
    nnodes++;

    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);
  }

  //find the beginig for the face information
  while (strcmp("CBAR",dat)==0) {
    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);
  }

  //read the face information
  while (strcmp("CTRIA3",dat)==0) {
    ifile.CutInputStr(dat,str);

    ifile.CutInputStr(dat,str);
    face.faceat=strtol(dat,&endptr,10);

    for (idim=0;idim<3;idim++) {
      ifile.CutInputStr(dat,str);
      face.node[idim]=strtol(dat,&endptr,10);
    }

    nfaces++;
    faces.push_back(face);

    ifile.GetInputStr(str,sizeof(str));
    ifile.CutInputStr(dat,str);
  }


  //check the distances between nodes
  long int nNode0,nNode1;
  long int nd,nfc,id;
  double d,*xNode0,*xNode1;

  for (nNode0=0;nNode0<nnodes;nNode0++) {
    xNode0=nodes[nNode0].x;

    for (nNode1=nNode0+1;nNode1<nnodes;nNode1++) {
      xNode1=nodes[nNode1].x;
      d=sqrt(pow(xNode1[0]-xNode0[0],2)+pow(xNode1[1]-xNode0[1],2)+pow(xNode1[2]-xNode0[2],2));

      if (d<EPS) {
        id=nodes[nNode1].id;

        //change the node's of the faces
        for (nfc=0;nfc<nfaces;nfc++) for (idim=0;idim<3;idim++) if (faces[nfc].node[idim]==id) faces[nfc].node[idim]=nodes[nNode0].id;

        //remove the node
        nodes.erase(nodes.begin()+nNode1);
        --nNode1;
        --nnodes;
      }
    }
  }

  //remove faces that have at least two identical node
  int *nlist;
  bool removeface;

  for (nfc=0;nfc<nfaces;nfc++) {
    nlist=faces[nfc].node;
    removeface=false;

    for (i=0;i<3;i++) for (j=i+1;j<3;j++) if (nlist[i]==nlist[j]) removeface=true; //remove the face

    if (removeface==true) {
      faces.erase(faces.begin()+nfc);
      --nfc;
      --nfaces;
    }
  }


  //renumerate nodes
  for (nfc=0;nfc<nfaces;nfc++) for (idim=0;idim<3;idim++) {
    id=faces[nfc].node[idim];

    for (nd=0;nd<nnodes;nd++) if (nodes[nd].id==id) {
      faces[nfc].node[idim]=nd;
      break;
    }
  }

  for (nd=0;nd<nnodes;nd++) nodes[nd].id=nd;


  //create the surface triangulation array
  nBoundaryTriangleFaces=nfaces;
  BoundaryTriangleFaces=new cTriangleFace[nfaces];

  nBoundaryTriangleNodes=nnodes;
  BoundaryTriangleNodes=new cNASTRANnode[nnodes];

  //copy the local nodes' vector into the array
  for (nd=0;nd<nnodes;nd++) {
    BoundaryTriangleNodes[nd]=nodes[nd];
    for (int idim=0;idim<3;idim++) BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim]=0.0;
  }

  for (nfc=0;nfc<nfaces;nfc++) {
    BoundaryTriangleFaces[nfc].SetFaceNodes(nodes[faces[nfc].node[0]].x,nodes[faces[nfc].node[1]].x,nodes[faces[nfc].node[2]].x);

    for (int idim=0;idim<3;idim++) BoundaryTriangleFaces[nfc].node[idim]=BoundaryTriangleNodes+faces[nfc].node[idim];
  }

  //calculate external normals to the faces
  double x0[3],e0[3],e1[3],SearchDirection[3],l,l0;
  long int nface;
  cNASTRANface *fcptr;
  cNASTRANnode *nd0,*nd1,*nd2;
  double xLocal[2];

  const double angleCosMin=cos(85.0/180.0*3.141592654);

  int nStartFace,nFinishFace,nTotalThreads,ThisThread,nFaceThread;

  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

  nFaceThread=nfaces/nTotalThreads;
  nStartFace=nFaceThread*ThisThread;
  nFinishFace=nStartFace+nFaceThread;
  if (ThisThread==nTotalThreads-1) nFinishFace=nfaces;



  for (nface=nStartFace;nface<nFinishFace;nface++) {
    fcptr=&faces[nface];

    nd0=&nodes[fcptr->node[0]];
    nd1=&nodes[fcptr->node[1]];
    nd2=&nodes[fcptr->node[2]];

    do {
      xLocal[0]=sqrt(rnd());
      xLocal[1]=(1.0-xLocal[0])*rnd();
    }
    while ((xLocal[0]<1.0E-4)||(xLocal[1]<1.0E-4)||(1.0-xLocal[0]-xLocal[1]<1.0E-4));

    for (idim=0,l=0.0,l0=0.0;idim<3;idim++) {
      e0[idim]=nd1->x[idim]-nd0->x[idim],e1[idim]=nd2->x[idim]-nd0->x[idim];
      x0[idim]=nd0->x[idim]+xLocal[0]*e0[idim]+xLocal[1]*e1[idim];

      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
      l0+=SearchDirection[idim]*BoundaryTriangleFaces[nface].ExternalNormal[idim];
    }

    l=sqrt(l);
    if (l0<0.0) l*=-1.0;

    for (idim=0;idim<3;idim++) SearchDirection[idim]/=l;

    //count face intersections
    int nIntersections=0;

    for (nfc=0;nfc<nfaces;nfc++) if (nfc!=nface) {
      if (BoundaryTriangleFaces[nfc].RayIntersection(x0,SearchDirection,EPS)==true) nIntersections++;
    }

    if (nIntersections%2!=0) {
      //the norm has to be reversed
      for (idim=0;idim<3;idim++) BoundaryTriangleFaces[nface].ExternalNormal[idim]*=-1.0;
    }
  }

  //collect the surface normals
  double sendBuffer[3*2*nFaceThread];
  int thread,cnt;

  for (thread=0;thread<nTotalThreads;thread++) {
    nStartFace=nFaceThread*thread;
    nFinishFace=nStartFace+nFaceThread;
    if (thread==nTotalThreads-1) nFinishFace=nfaces;

    if (thread==ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(sendBuffer+3*cnt,BoundaryTriangleFaces[nface].ExternalNormal,3*sizeof(double));
    }

    MPI_Bcast(sendBuffer,3*(nFinishFace-nStartFace),MPI_DOUBLE,thread,MPI_GLOBAL_COMMUNICATOR);

    if (thread!=ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(BoundaryTriangleFaces[nface].ExternalNormal,sendBuffer+3*cnt,3*sizeof(double));
    }
  }

  //set up the face attribute
  for (nfc=0;nfc<nfaces;nfc++) {
    memcpy(faces[nfc].externalNormal,BoundaryTriangleFaces[nfc].ExternalNormal,3*sizeof(double));
    BoundaryTriangleFaces[nfc].attribute=faces[nfc].faceat;
  }

  //calculate the ball averaged normal for the surface nodes
  for (nfc=0;nfc<nfaces;nfc++) for (int nd=0;nd<3;nd++) for (int idim=0;idim<3;idim++) {
    BoundaryTriangleFaces[nfc].node[nd]->BallAveragedExternalNormal[idim]+=BoundaryTriangleFaces[nfc].ExternalNormal[idim];
  }

  for (int nd=0;nd<nBoundaryTriangleNodes;nd++) {
    double l;
    int idim;

    for (idim=0,l=0.0;idim<3;idim++) l+=pow(BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim],2);

    if (l>1.0E-50) {
      for (idim=0,l=sqrt(l);idim<3;idim++) BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim]/=l;
    }
  }

}

//check weather a point (x0) in insed the domain:
//if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
//l -> is a random ray (intersection search) direction
bool CutCell::CheckPointInsideDomain(double *x,CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS) {
  int nface,nfaceStart,nfaceFinish,iIntersections;
  double SearchDirection[3],l;
  int idim;
  bool flag=true;

  if (SurfaceTriangulation==NULL) return true;

  static bool initflag=false;
  static int ThisThread,nTotalThreads;

  if (initflag==false) {
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
    MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);
    initflag=true;
  }

  if (ParallelCheck==true) {
    nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
    nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
    if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;
  }
  else nfaceStart=0,nfaceFinish=nSurfaceTriangulation;

  bool flagbuffer[nTotalThreads];

  do {
    //distribute ditrction of the search
    for (l=0.0,idim=0;idim<3;idim++) {
      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
    }

    if (ParallelCheck==true) MPI_Bcast(SearchDirection,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for (l=sqrt(l),idim=0;idim<3;idim++) SearchDirection[idim]/=l;
    iIntersections=0;
    flag=true;

    //find intersections with the faces on the mesh
    for (nface=nfaceStart;nface<nfaceFinish;nface++) {
      if (SurfaceTriangulation[nface].RayIntersection(x,SearchDirection,EPS)==true) iIntersections++;

      for (l=0.0,idim=0;idim<3;idim++) l+=pow(SurfaceTriangulation[nface].ExternalNormal[idim]*SearchDirection[idim],2);
      if (l<1.0E-10) {
        flag=false;
        break;
      }

    }

    if (ParallelCheck==true) {
      MPI_Gather(&flag,sizeof(bool),MPI_CHAR,flagbuffer,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
      if (ThisThread==0) for (int thread=1;thread<nTotalThreads;thread++) if (flagbuffer[thread]==false) flag=false;
      MPI_Bcast(&flag,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
    }
  }
  while (flag==false);

  if (ParallelCheck==true) {
    int t;

    MPI_Allreduce(&iIntersections,&t,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    iIntersections=t;
  }


  return (2*(iIntersections/2)==iIntersections) ? true : false;
}

bool CutCell::GetClosestSurfaceIntersectionPoint(double *x0,double *lSearch,double *xIntersection,double &tIntersection,CutCell::cTriangleFace* &FaceIntersection,CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
  int nface,nfaceStart,nfaceFinish,nFaceIntersection=-1;
  double t;

  FaceIntersection=NULL,tIntersection=0.0;

  static bool initflag=false;
  static int ThisThread,nTotalThreads;

  if (initflag==false) {
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
    MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);
    initflag=true;
  }

  MPI_Bcast(lSearch,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

  nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
  nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
  if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;

  for (nface=nfaceStart;nface<nfaceFinish;nface++) {
    if (SurfaceTriangulation[nface].RayIntersection(x0,lSearch,t,EPS)==true) {
      if (FaceIntersection==NULL) {
        FaceIntersection=SurfaceTriangulation+nface,tIntersection=t,nFaceIntersection=nface;
        for (int idim=0;idim<3;idim++) xIntersection[idim]=x0[idim]+t*lSearch[idim];
      }
      else {
        if (t<tIntersection) {
          FaceIntersection=SurfaceTriangulation+nface,tIntersection=t,nFaceIntersection=nface;
          for (int idim=0;idim<3;idim++) xIntersection[idim]=x0[idim]+t*lSearch[idim];
        }
      }
    }
  }

  struct cBuffer {
    int nFaceIntersection;
    double tIntersection;
  };

  cBuffer buffer[nTotalThreads],t1;


  buffer[0].nFaceIntersection=nFaceIntersection;
  buffer[0].tIntersection=tIntersection;

  cBuffer bufferRecv[nTotalThreads];

  MPI_Gather(buffer,sizeof(cBuffer),MPI_CHAR,bufferRecv,sizeof(cBuffer),MPI_CHAR,0,MPI_COMM_WORLD);
  memcpy(buffer,bufferRecv,nTotalThreads*sizeof(cBuffer));

  if (ThisThread==0) {
    for (int thread=1;thread<nTotalThreads;thread++) if (buffer[thread].nFaceIntersection!=-1) if ((tIntersection>buffer[thread].tIntersection)||(nFaceIntersection==-1)) {
      tIntersection=buffer[thread].tIntersection;
      nFaceIntersection=buffer[thread].nFaceIntersection;
    }
  }

  t1.nFaceIntersection=nFaceIntersection;
  t1.tIntersection=tIntersection;
  MPI_Bcast(&t1,sizeof(cBuffer),MPI_CHAR,0,MPI_COMM_WORLD);

  if (t1.nFaceIntersection!=-1) {
    FaceIntersection=SurfaceTriangulation+t1.nFaceIntersection,tIntersection=t1.tIntersection,nFaceIntersection=t1.nFaceIntersection;
    for (int idim=0;idim<3;idim++) xIntersection[idim]=x0[idim]+tIntersection*lSearch[idim];
  }

  return (FaceIntersection==NULL) ? false : true;
}


//calculate the part of the block that is within the computational domain

/*class cCutBlockNode {
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

class cQuadrilateral {
public:
  list<cCutBlockNode>::iterator node[8];
};

class cTriangleCutFace : public cTriangleFace {
public:
  list<cCutBlockNode>::iterator node[3];
};*/


double CutCell::GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, list<cTriangleFace*>& BlockTriangulationSet,int maxIntegrationLevel,int IntegrationLevel) {
  double VolumeL1=0.0;

/*
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
 */

  //determine is the block is within/outside of domain
  list<cTriangleFace*>::iterator BlockTriangulationItr;
  bool IntersectionFound=false;

  for (BlockTriangulationItr=BlockTriangulationSet.begin();BlockTriangulationItr!=BlockTriangulationSet.end();BlockTriangulationItr++) {
    if ((*BlockTriangulationItr)->BlockIntersection(xCellMin,xCellMax,EPS)==true) {
      IntersectionFound=true;
      break;
    }
  }

  //if no intersection is found -> determine weather the block is inside or outside of the domain
  if (IntersectionFound==false) {
  if (CheckPointInsideDomain(xCellMin,BoundaryTriangleFaces,nBoundaryTriangleFaces,false,EPS)==true) {
    //the block is within the domain
    return (xCellMax[0]-xCellMin[0])*(xCellMax[1]-xCellMin[1])*(xCellMax[2]-xCellMin[2]);
  }
  else return 0.0;
  }



  //perform integration
  if (IntegrationLevel==maxIntegrationLevel) {
    cCutBlock bl;
    list<cTetrahedron> indomainConnectivityList,outdomainConnectivityList;
    list<cTriangleCutFace> TriangleCutConnectivity;
    list<cTetrahedron>::iterator itr;

    VolumeL1=0.0;

    //populate the nodes of the block
    //nodes

    //check if the node is within the domain
    double xCellCenter[3];
    int IntersectionCounter;
    bool CellCenterInsideDomain;

    for (int ii=0;ii<3;ii++) xCellCenter[ii]=(xCellMax[ii]+xCellMin[ii])/2.0;

    CellCenterInsideDomain=CheckPointInsideDomain(xCellCenter,BoundaryTriangleFaces,nBoundaryTriangleFaces,false,EPS);


    for (int i=0;i<2;i+=1) for (int j=0;j<2;j+=1) for (int k=0;k<2;k+=1) {
      double x[3];
      list<cTriangleFace*>::iterator BlockTriangulationItr;

      x[0]=xCellMin[0]+i*(xCellMax[0]-xCellMin[0]);
      x[1]=xCellMin[1]+j*(xCellMax[1]-xCellMin[1]);
      x[2]=xCellMin[2]+k*(xCellMax[2]-xCellMin[2]);

      for (IntersectionCounter=0,BlockTriangulationItr=BlockTriangulationSet.begin();BlockTriangulationItr!=BlockTriangulationSet.end();BlockTriangulationItr++) {
        if ((*BlockTriangulationItr)->IntervalIntersection(x,xCellCenter,EPS)==true) {
          IntersectionCounter++;
        }
      }

      if ( ((CellCenterInsideDomain==true)&&(IntersectionCounter%2==0)) || ((CellCenterInsideDomain==false)&&(IntersectionCounter%2!=0)) ) {
        bl.AddNode(x,2*i,2*j,2*k);
      }
    }

    //add cut point on the edges
    //edges
    static const int EdgeMiddleNodeMap[12][3]={{1,0,0},{1,2,0},{1,2,2},{1,0,2},  {0,1,0},{2,1,0},{2,1,2},{0,1,2},  {0,0,1},{2,0,1},{2,2,1},{0,2,1}};
    static const int x0EdgeMap[12][3]={{0,0,0},{0,2,0},{0,2,2},{0,0,2},  {0,0,0},{2,0,0},{2,0,2},{0,0,2},  {0,0,0},{2,0,0},{2,2,0},{0,2,0}};
    static const int x1EdgeMap[12][3]={{2,0,0},{2,2,0},{2,2,2},{2,0,2},  {0,2,0},{2,2,0},{2,2,2},{0,2,2},  {0,0,2},{2,0,2},{2,2,2},{0,2,2}};

    for (int nedge=0;nedge<12;nedge++) {
      double x0[3],x1[3],xIntersection[3];
      list<cTriangleFace*>::iterator BlockTriangulationItr;

      for (int iii=0;iii<3;iii++) {
        x0[iii]=xCellMin[iii]+0.5*x0EdgeMap[nedge][iii]*(xCellMax[iii]-xCellMin[iii]);
        x1[iii]=xCellMin[iii]+0.5*x1EdgeMap[nedge][iii]*(xCellMax[iii]-xCellMin[iii]);
      }


      for (IntersectionCounter=0,BlockTriangulationItr=BlockTriangulationSet.begin();BlockTriangulationItr!=BlockTriangulationSet.end();BlockTriangulationItr++) {
        if ((*BlockTriangulationItr)->IntervalIntersection(x0,x1,xIntersection,EPS)==true) {
          int cnt=0;

          if (bl.node[x0EdgeMap[nedge][0]][x0EdgeMap[nedge][1]][x0EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
          if (bl.node[x1EdgeMap[nedge][0]][x1EdgeMap[nedge][1]][x1EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
          if ((cnt==0)||(cnt==2)) continue;

          bl.AddNode(xIntersection,EdgeMiddleNodeMap[nedge][0],EdgeMiddleNodeMap[nedge][1],EdgeMiddleNodeMap[nedge][2]);

          break;
        }
      }

    }

    //check the consistency of the cut edges


    static long int nCall=0;

    nCall++;

    if (nCall==239) {
      double a;

      a=10.0;
    }


    for (int nedge=0;nedge<12;nedge++) {
        int cnt=0;

        if (bl.node[x0EdgeMap[nedge][0]][x0EdgeMap[nedge][1]][x0EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
        if (bl.node[x1EdgeMap[nedge][0]][x1EdgeMap[nedge][1]][x1EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;

        if ((cnt==0)||(cnt==2)) {
          if (bl.node[EdgeMiddleNodeMap[nedge][0]][EdgeMiddleNodeMap[nedge][1]][EdgeMiddleNodeMap[nedge][2]]!=bl.NodeBuffer.end()) {
            //remove the cut-point
            bl.node[EdgeMiddleNodeMap[nedge][0]][EdgeMiddleNodeMap[nedge][1]][EdgeMiddleNodeMap[nedge][2]]=bl.NodeBuffer.end();
          }
        }
        else {
          if (bl.node[EdgeMiddleNodeMap[nedge][0]][EdgeMiddleNodeMap[nedge][1]][EdgeMiddleNodeMap[nedge][2]]==bl.NodeBuffer.end()) {
            double xIntersection[3];

            for (int iii=0;iii<3;iii++) {
              xIntersection[iii]=xCellMin[iii]+0.5*EdgeMiddleNodeMap[nedge][iii]*(xCellMax[iii]-xCellMin[iii]);
            }

            bl.AddNode(xIntersection,EdgeMiddleNodeMap[nedge][0],EdgeMiddleNodeMap[nedge][1],EdgeMiddleNodeMap[nedge][2]);
          }
        }
    }



    memcpy(bl.xBlockMin,xCellMin,3*sizeof(double));
    memcpy(bl.xBlockMax,xCellMax,3*sizeof(double));
    for (int i=0;i<3;i++) bl.dxBlock[i]=xCellMax[i]-xCellMin[i];

    cutBlockTetrahedronConnectivity(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);

 //   cutBlockTetrahedronConnectivity<cCutBlockNode,cCutBlock,cCutData,cTetrahedron,cTriangleCutFace>(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);
//  int cutBlockTetrahedronConnectivity(cCutBlock* bl,list<cTetrahedron>& indomainConnectivityList,list<cTetrahedron>& outdomainConnectivityList,list<cTriangleCutFace> TriangleCutFaceConnectivity) {

    for (itr=indomainConnectivityList.begin();itr!=indomainConnectivityList.end();itr++) {
      VolumeL1+=itr->Volume();
    }
  }
  else {
    //the level of the resolution does not reached the maximum allowed value
    //increase the refinment level
    list<cTriangleFace*> BlockTriangulationSubSet;
    list<cTriangleFace*>::iterator BlockTriangulationSetItr;
    int ii,jj,kk;
    double x0[3],x1[3];

      for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
        x0[0]=xCellMin[0]+0.5*ii*(xCellMax[0]-xCellMin[0]);
        x0[1]=xCellMin[1]+0.5*jj*(xCellMax[1]-xCellMin[1]);
        x0[2]=xCellMin[2]+0.5*kk*(xCellMax[2]-xCellMin[2]);

        x1[0]=x0[0]+0.5*(xCellMax[0]-xCellMin[0]);
        x1[1]=x0[1]+0.5*(xCellMax[1]-xCellMin[1]);
        x1[2]=x0[2]+0.5*(xCellMax[2]-xCellMin[2]);

        //determine the sub-set of the cut faces that intersects the refined cell
        BlockTriangulationSubSet.clear();

        for (BlockTriangulationItr=BlockTriangulationSet.begin();BlockTriangulationItr!=BlockTriangulationSet.end();BlockTriangulationItr++) {
          if ((*BlockTriangulationItr)->BlockIntersection(x0,x1,EPS)==true) {
          BlockTriangulationSubSet.push_back(*BlockTriangulationItr);
          }
        }

        //calcualte the volume
        VolumeL1+=GetRemainedBlockVolume(x0,x1,EPS,RelativeError,BlockTriangulationSubSet,maxIntegrationLevel,IntegrationLevel+1);
      }
  }



  return VolumeL1;
}


double CutCell::GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,CutCell::cTriangleFaceDescriptor* TriangleCutFaceDescriptorList) {
  int maxIntegrationLevel;

  maxIntegrationLevel=1+(int)(log(1.0/(pow(RelativeError,0.333))));

  if (maxIntegrationLevel<0) maxIntegrationLevel=0;
  if (maxIntegrationLevel>_AMR__CUT_CELL_VOLUME_CALCULATION__MAX_REFINMENT_LEVEL_) maxIntegrationLevel=_AMR__CUT_CELL_VOLUME_CALCULATION__MAX_REFINMENT_LEVEL_;

  //create the list of the intersecting triangular elements
  list<cTriangleFace*> BlockTriangulationSet;
  cTriangleFaceDescriptor *FaceDescriptor;

  for (FaceDescriptor=TriangleCutFaceDescriptorList;FaceDescriptor!=NULL;FaceDescriptor=FaceDescriptor->next) BlockTriangulationSet.push_back(FaceDescriptor->TriangleFace)  ;

  return GetRemainedBlockVolume(xCellMin,xCellMax,EPS,RelativeError,BlockTriangulationSet,maxIntegrationLevel,0);
}


void CutCell::SmoothRefine(double SmoothingCoefficient) {
  int nd,cl,idim;

  vector<cLocalNode> LocalNode;
  vector<cLocalEdge> LocalEdge;
  vector<cLocalTriangle> LocalTriangle,refinedLocalTriangle;

  //estimate the needed capacity of the vectors
  int EdgeNumber=(int)(1.5*3.0*nBoundaryTriangleFaces/2.0); //reserve 50% more just in case.....

  LocalNode.reserve(nBoundaryTriangleNodes+EdgeNumber);
  LocalEdge.reserve(EdgeNumber);
  LocalTriangle.reserve(nBoundaryTriangleFaces);
  refinedLocalTriangle.reserve(4*nBoundaryTriangleFaces);

  //recover the original triangulation
  for (nd=0;nd<nBoundaryTriangleNodes;nd++) {
    cLocalNode t;

    t.OriginalNode=BoundaryTriangleNodes+nd;
    t.OriginalNodeID=BoundaryTriangleNodes[nd].id;
    memcpy(t.x,BoundaryTriangleNodes[nd].x,3*sizeof(double));

    LocalNode.push_back(t);
  }

  for (cl=0;cl<nBoundaryTriangleFaces;cl++) {
    cLocalTriangle t;
    cLocalEdge edge;

    //get the nodes
    for (idim=0;idim<3;idim++) {
      nd=BoundaryTriangleFaces[cl].node[idim]-BoundaryTriangleNodes;
      t.node[idim]=LocalNode.begin()+nd;
    }

    //get the edges
    for (idim=0;idim<3;idim++) {
      nd=idim+1;
      if (nd==3) nd=0;
      edge.CornerNode[0]=t.node[nd];

      nd+=1;
      if (nd==3) nd=0;
      edge.CornerNode[1]=t.node[nd];

      //check if such edge already exists
      bool foundflag=false;

      for (vector<cLocalEdge>::iterator e=LocalEdge.begin();e!=LocalEdge.end();e++) {
        if ( ((e->CornerNode[0]==edge.CornerNode[0])&&(e->CornerNode[1]==edge.CornerNode[1])) || ((e->CornerNode[1]==edge.CornerNode[0])&&(e->CornerNode[0]==edge.CornerNode[1])) ) {
          //the edge exists:
          t.edge[idim]=e;
          foundflag=true;
          break;
        }
      }

      if (foundflag==false) {
        //create new edge
        cLocalNode newNode;

        //get the middle point
        for (int i=0;i<3;i++) newNode.x[i]=0.5*(edge.CornerNode[0]->x[i]+edge.CornerNode[1]->x[i]);

        newNode.OriginalNodeID=edge.CornerNode[0]->OriginalNodeID;
        LocalNode.push_back(newNode);
        edge.MiddleNode=LocalNode.begin()+LocalNode.size()-1;

        LocalEdge.push_back(edge);
        t.edge[idim]=LocalEdge.begin()+LocalEdge.size()-1;
      }
    }

    //add triangle to the list
    t.TriangleFace=BoundaryTriangleFaces+cl;
    LocalTriangle.push_back(t);
  }

  //refine cells
  for (vector<cLocalTriangle>::iterator tr=LocalTriangle.begin();tr!=LocalTriangle.end();tr++) {
    cLocalTriangle newTriangle;

    newTriangle.upTriangle=tr;

    //triangle: nd0,e2,e1
    newTriangle.node[0]=tr->node[0];
    newTriangle.node[1]=tr->edge[2]->MiddleNode;
    newTriangle.node[2]=tr->edge[1]->MiddleNode;
    refinedLocalTriangle.push_back(newTriangle);

    //triangle:e2,nd1,e0
    newTriangle.node[0]=tr->edge[2]->MiddleNode;
    newTriangle.node[1]=tr->node[1];
    newTriangle.node[2]=tr->edge[0]->MiddleNode;
    refinedLocalTriangle.push_back(newTriangle);

    //triangle: e1,e0,dn2;
    newTriangle.node[0]=tr->edge[1]->MiddleNode;
    newTriangle.node[1]=tr->edge[0]->MiddleNode;
    newTriangle.node[2]=tr->node[2];
    refinedLocalTriangle.push_back(newTriangle);


    //triangle:e0,e1,e2
    newTriangle.node[0]=tr->edge[0]->MiddleNode;
    newTriangle.node[1]=tr->edge[1]->MiddleNode;
    newTriangle.node[2]=tr->edge[2]->MiddleNode;
    refinedLocalTriangle.push_back(newTriangle);
  }

  //determine the list of triangles that share the same point
  for (vector<cLocalTriangle>::iterator tr=refinedLocalTriangle.begin();tr!=refinedLocalTriangle.end();tr++) {
    for (idim=0;idim<3;idim++) tr->node[idim]->ball.push_back(tr);
  }

  //determine smoother positions of the surface nodes
  for (vector<cLocalNode>::iterator node=LocalNode.begin();node!=LocalNode.end();node++) {
    int i,cnt=0;
    double x[3]={0.0,0.0,0.0};


    for (vector<vector<cLocalTriangle>::iterator>::iterator p=node->ball.begin();p!=node->ball.end();p++) {
      ++cnt;

      for (idim=0;idim<3;idim++) if ((*p)->node[idim]!=node) for (i=0;i<3;i++) x[i]+=0.5*(*p)->node[idim]->x[i];
    }

    for (idim=0;idim<3;idim++) node->xSmooth[idim]=(1.0-SmoothingCoefficient)*node->x[idim]+SmoothingCoefficient*x[idim]/(double)cnt;
  }

  //init the new array of the surface triangualtion
  cTriangleFace *newBoundaryTriangleFaces=new cTriangleFace[refinedLocalTriangle.size()];
  cNASTRANnode  *newBoundaryTriangleNodes=new cNASTRANnode[LocalNode.size()];
  int nnode,nface;
  vector<cLocalNode>::iterator node;
  vector<cLocalTriangle>::iterator tr;

  for (nnode=0,node=LocalNode.begin();node!=LocalNode.end();node++,nnode++) {
    newBoundaryTriangleNodes[nnode].id=node->OriginalNodeID;
    node->nodeno=nnode;
    memcpy(newBoundaryTriangleNodes[nnode].x,node->xSmooth,3*sizeof(double));
  }


  for (nface=0,tr=refinedLocalTriangle.begin();tr!=refinedLocalTriangle.end();tr++,nface++) {
    double c=0.0;

    newBoundaryTriangleFaces[nface].SetFaceNodes(tr->node[0]->xSmooth,tr->node[1]->xSmooth,tr->node[2]->xSmooth);
    newBoundaryTriangleFaces[nface].attribute=tr->upTriangle->TriangleFace->attribute;

    for (idim=0;idim<3;idim++) {
      newBoundaryTriangleFaces[nface].node[idim]=newBoundaryTriangleNodes+tr->node[idim]->nodeno;

      c+=newBoundaryTriangleFaces[nface].ExternalNormal[idim]*tr->upTriangle->TriangleFace->ExternalNormal[idim];
    }

    //the direction of the external notmal must coinside with that of the priginal face
    if (c<0.0) for (idim=0;idim<3;idim++) newBoundaryTriangleFaces[nface].ExternalNormal[idim]*=-1.0;
  }


  //update the arrays
  delete [] BoundaryTriangleFaces;
  delete [] BoundaryTriangleNodes;

  BoundaryTriangleFaces=newBoundaryTriangleFaces;
  nBoundaryTriangleFaces=refinedLocalTriangle.size();

  BoundaryTriangleNodes=newBoundaryTriangleNodes;
  nBoundaryTriangleNodes=LocalNode.size();


  //calculate the ball averaged normal for the surface nodes
  for (int nd=0;nd<nBoundaryTriangleNodes;nd++) for (int idim=0;idim<3;idim++) BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim]=0.0;

  for (int nfc=0;nfc<nBoundaryTriangleFaces;nfc++) for (int nd=0;nd<3;nd++) for (int idim=0;idim<3;idim++) {
    BoundaryTriangleFaces[nfc].node[nd]->BallAveragedExternalNormal[idim]+=BoundaryTriangleFaces[nfc].ExternalNormal[idim];
  }

  for (int nd=0;nd<nBoundaryTriangleNodes;nd++) {
    double l;
    int idim;

    for (idim=0,l=0.0;idim<3;idim++) l+=pow(BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim],2);

    if (l>1.0E-50) {
      for (idim=0,l=sqrt(l);idim<3;idim++) BoundaryTriangleNodes[nd].BallAveragedExternalNormal[idim]/=l;
    }
  }



}

int CutCell::cutBlockTetrahedronConnectivity(CutCell::cCutBlock* bl,list<CutCell::cTetrahedron>& indomainConnectivityList,list<CutCell::cTetrahedron>& outdomainConnectivityList,list<CutCell::cTriangleCutFace> &TriangleCutFaceConnectivity) {
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
      list<cCutBlockNode>::iterator meshNode;
      int ExtendedNodeMapID;

      bool ghostNode;

      cBlockNode() {nConnections=0,ghostNode=false;}

      void assign(list<cCutBlockNode>::iterator nd) {
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

        list<cCutBlockNode>::iterator thetraNodes[4];
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
                 list<cCutBlockNode>::iterator upNode;


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
        list<cCutBlockNode>::iterator thetraNodes[4];
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
               list<cCutBlockNode>::iterator upNode;


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



