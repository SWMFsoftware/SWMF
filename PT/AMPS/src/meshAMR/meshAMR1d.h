
//$Id$ 
//1D version of the AMR mesh

#ifndef _MESH_AMR_
#define _MESH_AMR_

//define the dimension of the mesh 
#define _MESH_DIMENSION_ 1

#include "meshAMRdef.h"




#include "meshAMRgeneric.h"

template <class cCornerNode,class cBlockAMR> 
class cMeshAMR1d : public cMeshAMRgeneric<cCornerNode,cBlockAMR> {
public:

  cMeshAMR1d(double *xmin,double *xmax,double (*localResolutionFunction)(double*)) : cMeshAMRgeneric<cCornerNode,cBlockAMR> (xmin,xmax,localResolutionFunction) {
  } 

  

};


#endif
