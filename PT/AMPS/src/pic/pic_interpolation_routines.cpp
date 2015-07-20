
//$Id$
//interpolation routines

/*
 * pic_interpolation_routines.cpp
 *
 *  Created on: Jul 18, 2015
 *      Author: vtenishe
 */

#include "pic.h"

PIC::InterpolationRoutines::CellCentered::cStencil PIC::InterpolationRoutines::CellCentered::Stencil;
int PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent=0;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::nBlockFoundMax];
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::last=NULL;

//definition of the interface functions that are used to access the FORTRAN written linear cell cenetered interpolation library
extern "C"{
  void interface__cell_centered_linear_interpolation__init_stencil_(int* nDim, double* XyzIn_D, int* nIndexes, int* nCell_D,  int* nGridOut, double* Weight_I, int* iIndexes_II, int* IsSecondOrder, int* UseGhostCell);
}


//determine stencil for the cell centered piecewise constant interpolation
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int i,j,k;
  PIC::Mesh::cDataCenterNode *cell;

  //flush the stencil
  PIC::InterpolationRoutines::CellCentered::Stencil.flush();

  //find the block
  if (node==NULL) {
    node=PIC::Mesh::mesh.findTreeNode(x);

    if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
    if (node->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not allocated");
  }

  //find cell
  cell=node->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

  //add the cell to the stencil
  if (cell!=NULL) PIC::InterpolationRoutines::CellCentered::Stencil.AddCell(1.0,cell);

  return &PIC::InterpolationRoutines::CellCentered::Stencil;
}

//determine the stencil for the cell centered linear interpolation using interpolation library from ../share/Library/src/
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(double *XyzIn_D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  //re-init variables in the INTERFACE and flush the Stencil
  INTERFACE::iBlockFoundCurrent=0;
  PIC::InterpolationRoutines::CellCentered::Stencil.flush();
  INTERFACE::last=node;

  //prepare the call of FORTRAN interpolate_amr subroutine
  int nDim       = DIM;
  // number of indices to identify cell ON A GIVEN PROCESSOR:
  // block id + DIM indices of cell in a block
  int nIndexes   = DIM + 1;
  int nCell_D[3] = {_BLOCK_CELLS_X_, _BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  // IMPORTANT NOTE: iIndexes_II has an extra index per cell which is
  // a processor ID => need one more index in addition to nIndexes
  int iIndexes_II[(DIM+1+1)*PIC::InterpolationRoutines::CellCentered::nMaxStencilLength] = {-1};
  double WeightStencil[PIC::InterpolationRoutines::CellCentered::nMaxStencilLength];
  int IsSecondOrder;
  int UseGhostCell=true;
  int nCellStencil;

  // call the interpolation subroutine

  interface__cell_centered_linear_interpolation__init_stencil_(&nDim,XyzIn_D,&nIndexes, nCell_D, &nCellStencil, WeightStencil, iIndexes_II, &IsSecondOrder, &UseGhostCell);


  // size of the stencil and weights are known
  // need to identify cells in the stencil
  // NOTE: FORTRAN is column-major for 2D arrays like iIndexes_II
  for (int iCellStencil = 0; iCellStencil < nCellStencil; iCellStencil++) {
    int ind[3]={0,0,0};
    PIC::Mesh::cDataBlockAMR  *block;
    PIC::Mesh::cDataCenterNode *cell;
    int iThread = iIndexes_II[0    +iCellStencil*(nIndexes+1)];

    for(int i = 0; i < DIM; i++) ind[i]=iIndexes_II[1+i  +iCellStencil*(nIndexes+1)];

    int iBlock  = iIndexes_II[1+DIM+iCellStencil*(nIndexes+1)];

    //retrieve the pointer to the current cell
    PIC::InterpolationRoutines::CellCentered::Stencil.AddCell(WeightStencil[iCellStencil],INTERFACE::BlockFound[iBlock]->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(ind[0],ind[1],ind[2])));
  }


  return &PIC::InterpolationRoutines::CellCentered::Stencil;
}
