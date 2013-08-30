//$Id$

#ifndef _CCUTBLOCKSET_
#define _CCUTBLOCKSET_


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


template <class cCutBlockNode, class cTetrahedron, class cQuadrilateral>
class cCutBlockSet {
public:
  int nTotalNodes,nTotalTetrahedronCells,nTotalQuadrilateralCells;
  list <cCutBlockNode> nodes;
  list <cTetrahedron> TetrahedronCells;
  list <cQuadrilateral> QuadrilateralCells;


  void clear() {
    nTotalNodes=0,nTotalTetrahedronCells=0,nTotalQuadrilateralCells=0;

    nodes.clear();
    TetrahedronCells.clear();
    QuadrilateralCells.clear();
  }

  cCutBlockSet() {
    clear();
  }

  void AddTetrahedron(cTetrahedron t,double EPS) {
    bool flag;
    int pnode;
    double *x,*xt;
    typename list<cCutBlockNode>::iterator nodeitr;

    //determine weather the nodes are in the node-list already
    for (pnode=0;pnode<4;pnode++) {
      flag=false;
      xt=t.node[pnode]->x;

      for (nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) {
        x=nodeitr->x;

        if ( ((x[0]-xt[0])*(x[0]-xt[0])) + ((x[1]-xt[1])*(x[1]-xt[1])) + ((x[2]-xt[2])*(x[2]-xt[2])) < EPS*EPS) {
          flag=true;
          t.node[pnode]=nodeitr;
          break;
        }
      }

      if (flag==false) {
        //the node is not in the list -> add it to the list
        nodes.push_front(*t.node[pnode]);
        nTotalNodes++;

        t.node[pnode]=nodes.begin();
      }
    }

    TetrahedronCells.push_back(t);
    nTotalTetrahedronCells++;
  }



  void AddTetrahedronList(list<cTetrahedron>& ConnectivityList,double EPS) {
    typename list<cTetrahedron>::iterator itr;

    for (itr=ConnectivityList.begin();itr!=ConnectivityList.end();itr++) AddTetrahedron(*itr,EPS);
  }

  void AddQuadrilateral(cQuadrilateral t,double EPS) {
    bool flag;
    int pnode;
    double *x,*xt;
    typename list<cCutBlockNode>::iterator nodeitr;

    //determine weather the nodes are in the node-list already
    for (pnode=0;pnode<8;pnode++) {
      flag=false;
      xt=t.node[pnode]->x;

      for (nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) {
        x=nodeitr->x;

        if ( ((x[0]-xt[0])*(x[0]-xt[0])) + ((x[1]-xt[1])*(x[1]-xt[1])) + ((x[2]-xt[2])*(x[2]-xt[2])) < EPS*EPS) {
          flag=true;
          t.node[pnode]=nodeitr;
          break;
        }
      }

      if (flag==false) {
        //the node is not in the list -> add it to the list
        nodes.push_front(*t.node[pnode]);
        nTotalNodes++;

        t.node[pnode]=nodes.begin();
      }
    }

    QuadrilateralCells.push_back(t);
    nTotalQuadrilateralCells++;
  }

  void AddQuadrilateralList(list<cQuadrilateral>& ConnectivityList,double EPS) {
    typename list<cQuadrilateral>::iterator itr;

    for (itr=ConnectivityList.begin();itr!=ConnectivityList.end();itr++) AddQuadrilateral(*itr,EPS);
  }


  void PrintTECPLOT(const char* fname) {
    FILE *fout;

    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"\n");

    //numerate and output the nodes
    typename list<cCutBlockNode>::iterator nodeitr;
    typename list<cTetrahedron>::iterator TetrahedronItr;
    typename list<cQuadrilateral>::iterator QuadrilateralItr;
    int cnt,pnode;

    //plot TetrahedronCells
    if (nTotalTetrahedronCells!=0) {
      for (nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) nodeitr->id=-1;

      for (cnt=0,TetrahedronItr=TetrahedronCells.begin();TetrahedronItr!=TetrahedronCells.end();TetrahedronItr++) {
        for (pnode=0;pnode<4;pnode++) if (TetrahedronItr->node[pnode]->id==-1) TetrahedronItr->node[pnode]->id=cnt++;
      }

      fprintf(fout,"ZONE N=%i, E=%i,DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON\n",cnt,nTotalTetrahedronCells);

      for (cnt=1,nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) if (nodeitr->id!=-1) {
        nodeitr->id=cnt++;
        fprintf(fout,"%e %e %e\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2]);
      }

      for (TetrahedronItr=TetrahedronCells.begin();TetrahedronItr!=TetrahedronCells.end();TetrahedronItr++) {
        fprintf(fout,"%i %i %i %i\n",TetrahedronItr->node[0]->id,TetrahedronItr->node[1]->id,TetrahedronItr->node[2]->id,TetrahedronItr->node[3]->id);
      }
    }


    //plot QuadrilateralCells
    if (nTotalQuadrilateralCells!=0) {
      for (nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) nodeitr->id=-1;

      for (cnt=0,QuadrilateralItr=QuadrilateralCells.begin();QuadrilateralItr!=QuadrilateralCells.end();QuadrilateralItr++) {
        for (pnode=0;pnode<8;pnode++) if (QuadrilateralItr->node[pnode]->id==-1) QuadrilateralItr->node[pnode]->id=cnt++;
      }

      fprintf(fout,"ZONE N=%i, E=%i,DATAPACKING=POINT, ZONETYPE=FEBRICK\n",cnt,nTotalQuadrilateralCells);

      for (cnt=1,nodeitr=nodes.begin();nodeitr!=nodes.end();nodeitr++) if (nodeitr->id!=-1) {
        nodeitr->id=cnt++;
        fprintf(fout,"%e %e %e\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2]);
      }

      for (QuadrilateralItr=QuadrilateralCells.begin();QuadrilateralItr!=QuadrilateralCells.end();QuadrilateralItr++) {
        fprintf(fout,"%i %i %i %i    %i %i %i %i\n",QuadrilateralItr->node[0]->id,QuadrilateralItr->node[1]->id,QuadrilateralItr->node[2]->id,QuadrilateralItr->node[3]->id,QuadrilateralItr->node[4]->id,QuadrilateralItr->node[5]->id,QuadrilateralItr->node[6]->id,QuadrilateralItr->node[7]->id);
      }
    }

    fclose(fout);
  }


};











#endif
