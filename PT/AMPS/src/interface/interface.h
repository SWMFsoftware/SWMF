//$Id$  
// header for interfacing to external modules written in different languages
#ifndef _INTERFACE_
#define _INTERFACE_

#include "pic.h"
#include "interface.dfn"

namespace INTERFACE {
  extern void Init();
}

#if _INTERFACE__INTERPOLATION_AMR__MODE_ == _INTERFACE_MODE_ON_
// INTERPOLATION_AMR procedure --------------------------------------------  

namespace INTERFACE{

  namespace INTERPOLATION_AMR {

    //list of pointers to nodes to identify them as integers
    //less than 2*2^3 integers might be needed to get an interpolation stencils
    const  int nBlockFoundMax=16;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* BlockFound[nBlockFoundMax];
    //index of the current position in the list
    extern int iBlockFoundCurrent;

    //list of pointers to cells in interpolation stencil
    const  int    nCellStencilMax=8;
    extern int    nCellStencil; 
    extern double WeightStencil[nCellStencilMax];
    extern PIC::Mesh::cDataCenterNode* CellStencil[nCellStencilMax];
    //the actual function to be called
    void interpolate_amr(double* XyzIn_D, int& nCell,
			 PIC::Mesh::cDataCenterNode** Cell, double* Weight);
  }

}

extern "C"{
  // interface function implemented in Fortran
  void INTERFACE__INTERPOLATION_AMR__interpolate_amr_(int* nDim, double* XyzIn_D, int* nIndexes, int* nCell_D,  int* nGridOut, double* Weight_I, int* iIndexes_II, bool* IsSecondOrder, bool* UseGhostCell);
  // interface function implemented in C called by interpolate_amr_
  void INTERFACE__INTERPOLATION_AMR__find_(int* nDim, double* Xyz_D, int* iProc, int* iBlock, double* XyzCorner_D, double* Dxyz_D, bool* IsOut);
}

#endif//_INTERFACE__INTERPOLATION_AMR_MODE_ == _INTERFACE_MODE_ON_

#endif//_INTERFACE_
