
//$Id$
//the c-side interface routine for the cell centered liniar interpolation library

#include "pic.h"

// procedure to find block containing the point of interest and its parameters
extern "C"{

  void interface__cell_centered_linear_interpolation__find_cpp_(int* nDim, double* Xyz_D, int* iProc, int* iBlock, double* XyzCorner_D, double* Dxyz_D, int* IsOut) {

    //check correctness
    if(DIM != *nDim)
      exit(__LINE__,__FILE__,"Error: inconsistent number of dimensions called by AMR interpolation procedure");

    // find the node containg the point
    static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    if(PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent==0) node = NULL;

    node=PIC::Mesh::mesh.findTreeNode(Xyz_D,node);

    // if node is not found => point is outside of the domain, break
    if(node==NULL) {
      //check if it is the original point => ERROR
      if(PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent==0)
	exit(__LINE__, __FILE__, "ERROR: point is outside of the domain");
      //otherwise it is simply close to the boundary (handled by the algorithm)
      *IsOut=1;
      return;
    }

    //node has been found, proceed
    *iProc=node->Thread;
    *iBlock = PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent;
    PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent++] = node;

    memcpy(XyzCorner_D, node->xmin, DIM*sizeof(double));

    Dxyz_D[0] = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
    Dxyz_D[1] = (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
    Dxyz_D[2] = (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;

    *IsOut=0;

    //change coordinates of the point by subtracting the corner coordinates
    for(int iDim = 0; iDim < *nDim; iDim++)
      Xyz_D[iDim] -= XyzCorner_D[iDim];
  }
}



