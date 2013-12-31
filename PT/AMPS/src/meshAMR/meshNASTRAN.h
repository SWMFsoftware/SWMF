//$Id$
//the class that reads the NASTRAN mesh

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

#include <list>
#include <vector>


#include "ifileopr.h"

#define MAXSTR 1000

using namespace std;



class cNASTRANmesh {
public:

  struct cNASTRANnode {
    double x[3];
    long int id;
  };

  struct cNASTRANface {
    long int node[3];
    long int faceat;
    double externalNormal[3];
  };

  long int nnodes,nfaces; ///,ncells;
  vector<cNASTRANnode> nodes;
  vector<cNASTRANface> faces;


  cNASTRANmesh() {
    nnodes=0,nfaces=0;   ////,ncells=0; 
    nodes.clear();
    faces.clear();
  }

  double rnd() {
    return ((double)(random())+1.0)/2147483649.0;
  }


//--------------- read NASTRAN mesh (long format) 
  void readMesh_longFormat(const char *fname,double *xGlobalMin,double *dxGlobal,int maxLevel,double EPS,double* rootBlock_Xmin,double* rootBlock_dX) {
    CiFileOperations ifile;
    char str[MAXSTR],dat[MAXSTR],*endptr;
    long int i,j,idim;
 
    cNASTRANnode node;
    cNASTRANface face;

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


      //check if coordinates of the node are too close to the nodes of the mesh 
      double dX,dXmin;

      for (idim=0;idim<3;idim++) {
        dXmin=rootBlock_dX[idim]/(1<<maxLevel);
        dX=(node.x[idim]-rootBlock_Xmin[idim])/dXmin; 

/*
        dX=(dX-(long int)dX)*dXmin;

        if (dX<EPS) node.x[idim]+=2.0*EPS; 
        else if (dXmin-dX<EPS) node.x[idim]-=2.0*EPS; 
*/

        dX-=(long int)dX;
        if (dX<0.1) node.x[idim]+=0.1*dXmin;
        else if (dX>0.9) node.x[idim]-=0.1*dXmin;


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
    register long int nNode0,nNode1;
    register long int nd,nfc,id;
    double d,*xNode0,*xNode1;

    for (nNode0=0;nNode0<nnodes;nNode0++) {
      xNode0=nodes[nNode0].x;
 
      for (nNode1=nNode0+1;nNode1<nnodes;nNode1++) { 
        xNode1=nodes[nNode1].x; 
        d=sqrt(pow(xNode1[0]-xNode0[0],2)+pow(xNode1[1]-xNode0[1],2)+pow(xNode1[2]-xNode0[2],2)); 

        if (d<3.0*EPS) {
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
    long int *nlist;
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

          if (fabs(t)<EPS) continue; //not intersection  

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
  }


//--------------- read NASTAN mesh (sort format) 
  void readMesh(const char *fname,double *xGlobalMin,double *dxGlobal,int maxLevel) {
    FILE *fd;
    long int initPosition,position,nd,i,idim,nfc;
    char str[MAXSTR],dat[MAXSTR],*endptr;

    //open the file 
    errno=0;
    fd=fopen(fname,"r");

    if (errno!=0) {
      printf("cannot open the mesh file %s\n",fname);
      exit(0);
    }

    //cout the number of nodes
    bool foundflag=false;

    while (!feof(fd)) {
      position=ftell(fd);
      fgets(str,MAXSTR,fd);

      if ((str[0]=='G')&&(str[3]=='D')) {
        if (foundflag==false) foundflag=true,initPosition=position;

        ++nnodes; 
      }
      else {
        if (foundflag==true) break;
      }
    } 

    if (foundflag==false) {
      cout << "ERROR: something wrong with the NASTRAN file -> too short" << endl;
      exit(0);
    }

    //load the nodes
    fseek(fd,initPosition,SEEK_SET);
    nodes.resize(nnodes);

    for (nd=0;nd<nnodes;nd++) {
      fgets(str,MAXSTR,fd);

      for (i=0;i<12;i++) dat[i]=str[i+4];
      dat[i]='\n';
      nodes[nd].id=strtol(dat,&endptr,10);

      for (idim=0;idim<3;idim++) {
        for (i=0;i<8;i++) dat[i]=str[i+24+8*idim];
        dat[i]='\n';

        //search for an exponential part in the node coordinte 
        for (i=1;i<8;i++) if ((dat[i]=='+')||(dat[i]=='-')) {
          bool foundflag=false;

          //search for symbols 'E' or 'e' in the line 
          for (i=1;i<8;i++) if ((dat[i]=='E')||(dat[i]=='e')) foundflag=true; 
          if (foundflag==true) break;

          //if 'E' or 'e'  is not found in the line -> insert them
          for (i=1;i<8;i++) if ((dat[i]=='+')||(dat[i]=='-')) {
            for (int j=9;j>i;j--) dat[j]=dat[j-1];
            dat[i]='E';
            break;
          }

          break;
        }

        for (i=1;i<8;i++) if (dat[i]==' ') {
          dat[i]='\n';
          break;
        }

        nodes[nd].x[idim]=strtod(dat,&endptr);
      }
    }


    //count the number of faces 
    foundflag=false;

    while (!feof(fd)) {
      position=ftell(fd);
      fgets(str,MAXSTR,fd);

      if ((str[0]=='C')&&(str[3]=='I')) {
        if (foundflag==false) foundflag=true,initPosition=position;

        ++nfaces;
      }
      else {
        if (foundflag==true) break;
      }
    }

    if (foundflag==false) {
      cout << "ERROR: something wrong with the NASTRAN file -> too short" << endl;
      exit(0);
    }

    //load the faces 
    fseek(fd,initPosition,SEEK_SET);
    faces.resize(nfaces);

    for (nfc=0;nfc<nfaces;nfc++) {
      fgets(str,MAXSTR,fd);

      for (i=0;i<8;i++) dat[i]=str[i+8];
      dat[i]='\n';
      /// cell.id=strtol(dat,&endptr,10);

      for (i=0;i<8;i++) dat[i]=str[i+16];
      dat[i]='\n';
      faces[nfc].faceat=strtol(dat,&endptr,10);

      for (idim=0;idim<3;idim++) {
        for (i=0;i<8;i++) dat[i]=str[i+24+8*idim];
        dat[i]='\n';
        faces[nfc].node[idim]=strtol(dat,&endptr,10);
      }
    }



    //renumerate nodes
    register long int id;


    for (nfc=0;nfc<nfaces;nfc++) for (idim=0;idim<3;idim++) {
      id=faces[nfc].node[idim];

      for (nd=0;nd<nnodes;nd++) if (nodes[nd].id==id) faces[nfc].node[idim]=nd;
    } 




    //calculate external normals to the faces 
    bool faceIntersectionFlag;
    double x[3],x0[3],angle,e0[3],e1[3],norm[3],searchDirection[3],l;
    long int nface,nIntersectionsForward,nIntersectionsBackward;
    cNASTRANface *fcptr,*fc;
    cNASTRANmesh::cNASTRANnode *nd0,*nd1,*nd2;
    double t,c,c0,c1,c00,c11,c01,xLocal[2];
    
  

    const double angleCosMin=cos(85.0/180.0*3.141592654);

    for (nface=0;nface<nfaces;nface++) {
      fcptr=&faces[nface];

      nd0=&nodes[fcptr->node[0]];
      nd1=&nodes[fcptr->node[1]];
      nd2=&nodes[fcptr->node[2]];


////////

if ((-970<nd0->x[0])&&(nd0->x[0]<-636)&&(0<nd0->x[1])&&(nd0->x[1]<300)&&(20<nd0->x[2])&&(nd0->x[2]<150)) {
cout << __LINE__ << endl;
}  
 

/////



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
  }


  //intersection of a ray (x=x0+l*t, t in [0,\infty] ) with a face

  #define _NODE_INTERSECTION_NOTDETERMINED_NASTRAN_  -1 
  #define _NODE_INTERSECTION_FALSE_NASTRAN_           0 
  #define _NODE_INTERSECTION_TRUE_NASTRAN_            1
  #define _NODE_INTERSECTION_BOUNDARY_NASTRAN_        2

  int faceIntersection(double *x0,double *l,cNASTRANface *face,double EPS) {
    int idim;
    double length,e0[3],e1[3],norm[3],xLocal[2];  
    double t,c0,c1,c00,c11,c01,c,x[3];
    cNASTRANnode *nd0,*nd1,*nd2;

    nd0=&nodes[face->node[0]];
    nd1=&nodes[face->node[1]];
    nd2=&nodes[face->node[2]];

    for (idim=0;idim<3;idim++) e0[idim]=nd1->x[idim]-nd0->x[idim],e1[idim]=nd2->x[idim]-nd0->x[idim];  

    //get a normal to the face 
    norm[0]=e0[1]*e1[2]-e0[2]*e1[1];
    norm[1]=-(e0[0]*e1[2]-e0[2]*e1[0]);
    norm[2]=e0[0]*e1[1]-e0[1]*e1[0];

    length=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
    
    norm[0]/=length;
    norm[1]/=length;
    norm[2]/=length;

    //compare the found locations woth EPS
    double eEPS[2];

    eEPS[0]=EPS/sqrt(e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2]);
    eEPS[1]=EPS/sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);

   //check is the node is a boundary node ar the face: 1. get the distance between the node and the face and 2. calculate the prohection of the node on the face
   //get the distance between the node and the face 
   for (idim=0,length=0.0;idim<3;idim++) length+=pow((x0[idim]-nd0->x[idim])*norm[idim],2);
  
   if (length<EPS*EPS) {
     //the point "x" is within the distance 'EPS' from the plane that contains the cut face 'face'
     //get the local coordinated of the point projection of the face 

     for (c0=0.0,c1=0.0,c00=0.0,c01=0.0,c11=0.0,idim=0;idim<3;idim++) {
       x[idim]=x0[idim]-nd0->x[idim];

       c1+=x[idim]*e1[idim],c0+=x[idim]*e0[idim];
       c00+=e0[idim]*e0[idim],c11+=e1[idim]*e1[idim],c01+=e0[idim]*e1[idim];
     }

     c=c11*c00-c01*c01;
     xLocal[0]=(c0*c11-c1*c01)/c;
     xLocal[1]=(c1*c00-c01*c0)/c;

     //determine weather the node in outside of the face
     for (idim=0;idim<2;idim++) if ((xLocal[idim]<-eEPS[idim])||(xLocal[idim]>1.0+eEPS[idim])) return _NODE_INTERSECTION_FALSE_NASTRAN_;
     if (xLocal[0]+xLocal[1]>1.0+max(eEPS[0],eEPS[1])) return _NODE_INTERSECTION_FALSE_NASTRAN_; 

     if ((fabs(xLocal[0])<eEPS[0])||(fabs(xLocal[1])<eEPS[1])||(fabs(xLocal[0]+xLocal[1]-1.0)<max(eEPS[0],eEPS[1]))) return _NODE_INTERSECTION_NOTDETERMINED_NASTRAN_;  

     return _NODE_INTERSECTION_BOUNDARY_NASTRAN_;
   }


    //the initial point cannot be a "boundary node"
    //get the time of the intersection of the ray with the  
    for (c0=0.0,c1=0.0,idim=0;idim<3;idim++) c0+=(x0[idim]-nd0->x[idim])*norm[idim],c1+=l[idim]*norm[idim];

    if (fabs(c1)<1.0E-15*fabs(c0)) return _NODE_INTERSECTION_FALSE_NASTRAN_;
    t=-c0/c1;

    if (t<0.0) return _NODE_INTERSECTION_FALSE_NASTRAN_; 

    //get local coordanated of the point of the intersection
    for (c0=0.0,c1=0.0,c00=0.0,c01=0.0,c11=0.0,idim=0;idim<3;idim++) {
      x[idim]=x0[idim]+l[idim]*t-nd0->x[idim];
      
      c1+=x[idim]*e1[idim],c0+=x[idim]*e0[idim];
      c00+=e0[idim]*e0[idim],c11+=e1[idim]*e1[idim],c01+=e0[idim]*e1[idim];
    }

    c=c11*c00-c01*c01;
    xLocal[0]=(c0*c11-c1*c01)/c;
    xLocal[1]=(c1*c00-c01*c0)/c;  


    //determine weather the node in outside of the face
    for (idim=0;idim<2;idim++) if ((xLocal[idim]<-eEPS[idim])||(xLocal[idim]>1.0+eEPS[idim])) return _NODE_INTERSECTION_FALSE_NASTRAN_;
    if (xLocal[0]+xLocal[1]>1.0+max(eEPS[0],eEPS[1])) return _NODE_INTERSECTION_FALSE_NASTRAN_;

    //determine weather the node in on the boundary of the face 
    for (idim=0;idim<2;idim++) if ((fabs(xLocal[idim])<eEPS[idim])||(fabs(xLocal[idim]-1.0)<eEPS[idim])) {
      return _NODE_INTERSECTION_NOTDETERMINED_NASTRAN_;  
    }
  
    if (fabs(xLocal[0]+xLocal[1]-1.0)<max(eEPS[0],eEPS[1])) {
      return _NODE_INTERSECTION_NOTDETERMINED_NASTRAN_;   
    }

    return _NODE_INTERSECTION_TRUE_NASTRAN_;
  }


  //check weather a point (x0) in insed the domain:
  //if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
  //l -> is a random ray (intersection search) direction
  bool checkPointInsideDomain(double *x,double EPS) {
    long int nface,iIntersections=0;
    double searchDirection[3],l;
    int idim,returnCode;

    do {
      //distribute ditrction of the search  
      for (l=0.0,idim=0;idim<3;idim++) {
        searchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
        l+=pow(searchDirection[idim],2);
      }

      for (l=sqrt(l),idim=0;idim<3;idim++) searchDirection[idim]/=l; 

      iIntersections=0; 
    
      //find intersections with the faces on the mesh 
      for (nface=0;nface<nfaces;nface++) { 
        returnCode=faceIntersection(x,searchDirection,&faces[nface],EPS);

        if (returnCode==_NODE_INTERSECTION_BOUNDARY_NASTRAN_) {
          return true; 
        }


        if (returnCode==_NODE_INTERSECTION_NOTDETERMINED_NASTRAN_) break; 
        if (returnCode==_NODE_INTERSECTION_TRUE_NASTRAN_) iIntersections++; 
      }
    }
    while (returnCode==_NODE_INTERSECTION_NOTDETERMINED_NASTRAN_);

    return (2*(iIntersections/2)==iIntersections) ? true : false;
  }  

 


 


};
 



