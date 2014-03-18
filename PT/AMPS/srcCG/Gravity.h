#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "ifileopr.h"

#define MAXSTR 1000

using namespace std;

#ifndef nucleusGravityCalculator
#define nucleusGravityCalculator

namespace nucleusGravity {

  struct cNASTRANnode {
    double x[3];
    long int id;
  };

  struct cNASTRANtetra {
    long int node[4];
  };

  extern cNASTRANnode * nodes;
  extern cNASTRANtetra * tetras;

  extern long int nnodes,ntetras; ///,ncells;

  extern double density;
    
  void setDensity(double d); 

  void readMesh_longformat(char *fname);

  void readMesh(const char *fname);

  
  void  gravity(double * gravitypermass, double * position); 
  
  void printMeshFile();
  
  void clearArrays();

};
#endif
