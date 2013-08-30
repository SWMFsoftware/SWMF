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

cTriangleFace *BoundaryTriangleFaces=NULL;
int nBoundaryTriangleFaces=0;
cAMRstack<cTriangleFaceDescriptor> BoundaryTriangleFaceDescriptor;

struct cNodeCoordinates {
  double *x;
  int id;
};

struct cFaceNodeConnection {
  list<cNodeCoordinates>::iterator node[3];
};

void PrintSurfaceTriangulationMesh(const char *fname,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
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
      nodeCoordinates.push_front(nd);

      FaceNodeConnection[nface].node[pnode]=nodeCoordinates.begin();
    }
  }


  //print the mesh
  FILE *fout=fopen(fname,"w");
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",(int)nodeCoordinates.size(),nSurfaceTriangulation);


  for (cnt=1,nodeitr=nodeCoordinates.begin();nodeitr!=nodeCoordinates.end();nodeitr++) {
    nodeitr->id=cnt++;
    fprintf(fout,"%e %e %e\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2]);
  }

  for (nface=0;nface<nSurfaceTriangulation;nface++) {
    fprintf(fout,"%i %i %i\n",FaceNodeConnection[nface].node[0]->id,FaceNodeConnection[nface].node[1]->id,FaceNodeConnection[nface].node[2]->id);
  }

  fclose(fout);
  delete [] FaceNodeConnection;
}

struct cNASTRANnode {
  double x[3];
  int id;
};

struct cNASTRANface {
  int node[3],faceat;
  double externalNormal[3];
};

void ReadNastranSurfaceMeshLongFormat(const char *fname,cTriangleFace* &SurfaceTriangulation,int &nSurfaceTriangulation,double *xSurfaceMin,double *xSurfaceMax,double EPS) {
  CiFileOperations ifile;
  char str[10000],dat[10000],*endptr;
  long int i,j,idim,nnodes=0,nfaces=0;

  if (SurfaceTriangulation!=NULL) exit(__LINE__,__FILE__,"Error: redifinition of the surface triangulation array");


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


  //calculate external normals to the faces
  bool faceIntersectionFlag;
  double x[3],x0[3],angle,e0[3],e1[3],norm[3],searchDirection[3],l;
  long int nface,nIntersectionsForward,nIntersectionsBackward;
  cNASTRANface *fcptr,*fc;
  cNASTRANnode *nd0,*nd1,*nd2;
  double t,c,c0,c1,c00,c11,c01,xLocal[2];



  const double angleCosMin=cos(85.0/180.0*3.141592654);

  for (nface=0;nface<nfaces;nface++) {
    fcptr=&faces[nface];

    nd0=&nodes[fcptr->node[0]];
    nd1=&nodes[fcptr->node[1]];
    nd2=&nodes[fcptr->node[2]];

    for (idim=0;idim<3;idim++) {
      e0[idim]=nd1->x[idim]-nd0->x[idim],e1[idim]=nd2->x[idim]-nd0->x[idim];
      x0[idim]=nd0->x[idim]+0.25*(e0[idim]+e1[idim]);
    }

    norm[0]=e0[1]*e1[2]-e0[2]*e1[1];
    norm[1]=-(e0[0]*e1[2]-e0[2]*e1[0]);
    norm[2]=e0[0]*e1[1]-e0[1]*e1[0];

    l=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
    for (idim=0;idim<3;idim++) fcptr->externalNormal[idim]=norm[idim]/l;

    //scan trought the faces to determine the external normal
    do {
      nIntersectionsForward=0,nIntersectionsBackward=0;
      faceIntersectionFlag=true;

      //generate the search direction
      do {
        for (l=0.0,idim=0;idim<3;idim++) {
          searchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
          l+=pow(searchDirection[idim],2);
        }

        for (angle=0.0,l=sqrt(l),idim=0;idim<3;idim++) {
          searchDirection[idim]/=l;
          angle+=searchDirection[idim]*fcptr->externalNormal[idim];
        }
      }
      while (fabs(angle)<angleCosMin);

      //get the time and local coordinates of an intersection between the cut face and the direction of the search
      for (nfc=0;nfc<nfaces;nfc++) if (nfc!=nface) {
        fc=&faces[nfc];

        nd0=&nodes[fc->node[0]];
        nd1=&nodes[fc->node[1]];
        nd2=&nodes[fc->node[2]];

        for (idim=0;idim<3;idim++) e0[idim]=nd1->x[idim]-nd0->x[idim],e1[idim]=nd2->x[idim]-nd0->x[idim];

        //get a normal to the face
        norm[0]=e0[1]*e1[2]-e0[2]*e1[1];
        norm[1]=-(e0[0]*e1[2]-e0[2]*e1[0]);
        norm[2]=e0[0]*e1[1]-e0[1]*e1[0];

        for (c0=0.0,c1=0.0,idim=0;idim<3;idim++) c0+=(x0[idim]-nd0->x[idim])*norm[idim],c1+=searchDirection[idim]*norm[idim];

        if (fabs(c1)<1.0E-15*fabs(c0)) {
          faceIntersectionFlag=false;
          break;
        }

        t=-c0/c1;

        //get local coordanated of the point of the intersection
        for (c0=0.0,c1=0.0,c00=0.0,c01=0.0,c11=0.0,idim=0;idim<3;idim++) {
          x[idim]=x0[idim]+searchDirection[idim]*t-nd0->x[idim];

          c1+=x[idim]*e1[idim],c0+=x[idim]*e0[idim];
          c00+=e0[idim]*e0[idim],c11+=e1[idim]*e1[idim],c01+=e0[idim]*e1[idim];
        }

        c=c11*c00-c01*c01;
        xLocal[0]=(c0*c11-c1*c01)/c;
        xLocal[1]=(c1*c00-c01*c0)/c;


        //determine weather the node in outside of the face
        if (((fabs(xLocal[0])<1.0E-5)||(fabs(xLocal[0]-1.0)<1.0E-5))&&(1.0E-5<xLocal[1])&&(xLocal[0]+xLocal[1]<1.0-1.0E-5)) {
          faceIntersectionFlag=false;
          break;
        }

        if (((fabs(xLocal[1])<1.0E-5)||(fabs(xLocal[1]-1.0)<1.0E-5))&&(1.0E-5<xLocal[0])&&(xLocal[0]+xLocal[1]<1.0-1.0E-5)) {
          faceIntersectionFlag=false;
          break;
        }

        if ((xLocal[0]>-1.0E-5)&&(xLocal[1]>-1.0E-5)&&(fabs(xLocal[0]+xLocal[1]-1.0)<1.0E-5)) {
          faceIntersectionFlag=false;
          break;
        }

        if ((xLocal[0]>0.0)&&(xLocal[1]>0.0)&&(xLocal[0]+xLocal[1]<1.0)) {
          if (t>0.0) nIntersectionsForward++; else nIntersectionsBackward++;
        }
      }
    }
    while (faceIntersectionFlag==false);

    //deternime the direction of the external normal vector

    if (angle<0.0) nIntersectionsForward=nIntersectionsBackward;
    if (faceIntersectionFlag==true) if (2*(nIntersectionsForward/2)==nIntersectionsForward) {
      //the forward direction along the chosen line of the search does not coinsides with the direction of the external normal
      for (idim=0;idim<3;idim++) fcptr->externalNormal[idim]*=-1.0;
    }
  }

  //create the surface triangulation array
  nSurfaceTriangulation=nfaces;
  SurfaceTriangulation=new cTriangleFace[nfaces];

  for (nfc=0;nfc<nfaces;nfc++) {
    SurfaceTriangulation[nfc].SetFaceNodes(nodes[faces[nfc].node[0]].x,nodes[faces[nfc].node[1]].x,nodes[faces[nfc].node[2]].x);

    memcpy(SurfaceTriangulation[nfc].ExternalNormal,faces[nfc].externalNormal,3*sizeof(double));
    SurfaceTriangulation[nfc].attribute=faces[nfc].faceat;
  }
}

//check weather a point (x0) in insed the domain:
//if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
//l -> is a random ray (intersection search) direction
bool CheckPointInsideDomain(double *x,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
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

  nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
  nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
  if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;

  bool flagbuffer[nTotalThreads];

  do {
    //distribute ditrction of the search
    for (l=0.0,idim=0;idim<3;idim++) {
      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
    }

    MPI_Bcast(SearchDirection,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

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

    MPI_Gather(&flag,sizeof(bool),MPI_CHAR,flagbuffer,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
    if (ThisThread==0) for (int thread=1;thread<nTotalThreads;thread++) if (flagbuffer[thread]==false) flag=false;
    MPI_Bcast(&flag,sizeof(bool),MPI_CHAR,0,MPI_COMM_WORLD);
  }
  while (flag==false);

  MPI_Bcast(&iIntersections,1,MPI_INT,0,MPI_COMM_WORLD);
  return (2*(iIntersections/2)==iIntersections) ? true : false;
}

bool GetClosestSurfaceIntersectionPoint(double *x0,double *lSearch,double *xIntersection,double &tIntersection,cTriangleFace* &FaceIntersection,cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
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
  MPI_Gather(buffer,sizeof(cBuffer),MPI_CHAR,buffer,sizeof(cBuffer),MPI_CHAR,0,MPI_COMM_WORLD);

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

class cQuadrilateral {
public:
  list<cCutBlockNode>::iterator node[8];
};

class cTriangleCutFace : public cTriangleFace {
public:
  list<cCutBlockNode>::iterator node[3];
};


double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,cTriangleFaceDescriptor* TriangleCutFaceDescriptorList,int maxIntegrationLevel,int IntegrationLevel) {
  int BlockIntersectionCode;
  double VolumeL0=-1.0,VolumeL1=-1.0;




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


  //determine is the block is within/outside of domain
  cTriangleFaceDescriptor *FaceDescriptor;
  bool IntersectionFound=false;

  for (FaceDescriptor=TriangleCutFaceDescriptorList;FaceDescriptor!=NULL;FaceDescriptor=FaceDescriptor->next) {
    if (FaceDescriptor->TriangleFace->BlockIntersection(xCellMin,xCellMax,EPS)==true) {
      IntersectionFound=true;
      break;
    }
  }

  //if no intersection is found -> determin weather the block is insode or outside of the domain
  if (IntersectionFound==false) {
    int idim;
    double c=0.0,*x0=TriangleCutFaceDescriptorList->TriangleFace->x0Face,*norm=TriangleCutFaceDescriptorList->TriangleFace->ExternalNormal;

    for (idim=0;idim<3;idim++) {
      c+=(xCellMin[idim]-x0[idim])*norm[idim];
    }

    if (c>0.0) {
      //the block is outsed of the domain
      return 0.0;
    }
    else {
      //the block is within the domain
      return (xCellMax[0]-xCellMin[0])*(xCellMax[1]-xCellMin[1])*(xCellMax[2]-xCellMin[2]);
    }
  }



  //perform integration
  if ((IntegrationLevel==0)||(IntegrationLevel==maxIntegrationLevel)) {
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

    CellCenterInsideDomain=CheckPointInsideDomain(xCellCenter,SurfaceTriangulation,nSurfaceTriangulation,EPS);


    for (int i=0;i<2;i+=1) for (int j=0;j<2;j+=1) for (int k=0;k<2;k+=1) {
      double r,x[3],locx[3]={0.0,0.0,0.0};

      x[0]=xCellMin[0]+i*(xCellMax[0]-xCellMin[0]);
      x[1]=xCellMin[1]+j*(xCellMax[1]-xCellMin[1]);
      x[2]=xCellMin[2]+k*(xCellMax[2]-xCellMin[2]);

      for (IntersectionCounter=0,FaceDescriptor=TriangleCutFaceDescriptorList;FaceDescriptor!=NULL;FaceDescriptor=FaceDescriptor->next) {
        if (FaceDescriptor->TriangleFace->IntervalIntersection(x,xCellCenter,EPS)==true) {
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

      for (int iii=0;iii<3;iii++) {
        x0[iii]=xCellMin[iii]+0.5*x0EdgeMap[nedge][iii]*(xCellMax[iii]-xCellMin[iii]);
        x1[iii]=xCellMin[iii]+0.5*x1EdgeMap[nedge][iii]*(xCellMax[iii]-xCellMin[iii]);
      }


      for (IntersectionCounter=0,FaceDescriptor=TriangleCutFaceDescriptorList;FaceDescriptor!=NULL;FaceDescriptor=FaceDescriptor->next) {
        if (FaceDescriptor->TriangleFace->IntervalIntersection(x0,x1,xIntersection,EPS)==true) {
          int cnt=0;

          if (bl.node[x0EdgeMap[nedge][0]][x0EdgeMap[nedge][1]][x0EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
          if (bl.node[x1EdgeMap[nedge][0]][x1EdgeMap[nedge][1]][x1EdgeMap[nedge][2]]!=bl.NodeBuffer.end()) cnt++;
          if ((cnt==0)||(cnt==2)) continue;

          bl.AddNode(xIntersection,EdgeMiddleNodeMap[nedge][0],EdgeMiddleNodeMap[nedge][1],EdgeMiddleNodeMap[nedge][2]);

          break;
        }
      }

    }

    memcpy(bl.xBlockMin,xCellMin,3*sizeof(double));
    memcpy(bl.xBlockMax,xCellMax,3*sizeof(double));
    for (int i=0;i<3;i++) bl.dxBlock[i]=xCellMax[i]-xCellMin[i];

//    cutBlockTetrahedronConnectivity<cCutBlockNode,cCutBlock,cCutData,cTetrahedron>(&bl,indomainConnectivityList,outdomainConnectivityList,TriangleCutConnectivity);

    for (itr=indomainConnectivityList.begin();itr!=indomainConnectivityList.end();itr++) {
      VolumeL1+=itr->Volume();
    }
  }

  if (IntegrationLevel!=maxIntegrationLevel) {
    int ii,jj,kk;

    if (IntegrationLevel==0) {
      if (maxIntegrationLevel==0) {
        return VolumeL1;
      }
      else {
        int upperIntegrationLevel=0;

        do {
          double x0[3],x1[3];

          VolumeL0=VolumeL1;
          VolumeL1=0.0;
          upperIntegrationLevel+=1;

          for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
            x0[0]=xCellMin[0]+0.5*ii*(xCellMax[0]-xCellMin[0]);
            x0[1]=xCellMin[1]+0.5*jj*(xCellMax[1]-xCellMin[1]);
            x0[2]=xCellMin[2]+0.5*kk*(xCellMax[2]-xCellMin[2]);

            x1[0]=x0[0]+0.5*(xCellMax[0]-xCellMin[0]);
            x1[1]=x0[1]+0.5*(xCellMax[1]-xCellMin[1]);
            x1[2]=x0[2]+0.5*(xCellMax[2]-xCellMin[2]);

            VolumeL1+=GetRemainedBlockVolume(x0,x1,EPS,RelativeError,SurfaceTriangulation,nSurfaceTriangulation,TriangleCutFaceDescriptorList,upperIntegrationLevel,1);
          }

          if (upperIntegrationLevel==maxIntegrationLevel) break;
          if (VolumeL1<1.0E-6*(xCellMax[0]-xCellMin[0])*(xCellMax[1]-xCellMin[1])*(xCellMax[2]-xCellMin[2])) break;
          if (2.0*fabs(VolumeL1-VolumeL0)/(VolumeL1+VolumeL0)<RelativeError) break;
        } while (true);

        return VolumeL1;
      }
    }
    else {
      double x0[3],x1[3];
       VolumeL1=0.0;

       for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
         x0[0]=xCellMin[0]+0.5*ii*(xCellMax[0]-xCellMin[0]);
         x0[1]=xCellMin[1]+0.5*jj*(xCellMax[1]-xCellMin[1]);
         x0[2]=xCellMin[2]+0.5*kk*(xCellMax[2]-xCellMin[2]);

         x1[0]=x0[0]+0.5*(xCellMax[0]-xCellMin[0]);
         x1[1]=x0[1]+0.5*(xCellMax[1]-xCellMin[1]);
         x1[2]=x0[2]+0.5*(xCellMax[2]-xCellMin[2]);

         VolumeL1+=GetRemainedBlockVolume(x0,x1,EPS,RelativeError,SurfaceTriangulation,nSurfaceTriangulation,TriangleCutFaceDescriptorList,maxIntegrationLevel,IntegrationLevel+1);
       }
    }
  }

  return VolumeL1;
}


double GetRemainedBlockVolume(double* xCellMin,double* xCellMax,double EPS,double RelativeError, cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,cTriangleFaceDescriptor* TriangleCutFaceDescriptorList) {
  int maxIntegrationLevel;

  maxIntegrationLevel=1+(int)(log(1.0/(pow(RelativeError,0.333))));

  if (maxIntegrationLevel<0) maxIntegrationLevel=0;
  if (maxIntegrationLevel>10) maxIntegrationLevel=10;

  return GetRemainedBlockVolume(xCellMin,xCellMax,EPS,RelativeError,SurfaceTriangulation,nSurfaceTriangulation,TriangleCutFaceDescriptorList,maxIntegrationLevel,0);
}

