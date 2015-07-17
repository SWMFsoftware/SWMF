//$Id$  

//implementation of interface functions for AMR interpolation procedure

#include "interface.h"

//block found while constructing an interpolation stencil
int INTERFACE::INTERPOLATION_AMR::iBlockFoundCurrent=-1;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* INTERFACE::INTERPOLATION_AMR::BlockFound[nBlockFoundMax]={NULL};

//parameters of the final interpolation stencil
int    INTERFACE::INTERPOLATION_AMR::nCellStencil=0;
double INTERFACE::INTERPOLATION_AMR::WeightStencil[nCellStencilMax]={0.0};
PIC::Mesh::cDataCenterNode* INTERFACE::INTERPOLATION_AMR::CellStencil[nCellStencilMax];

void INTERFACE::INTERPOLATION_AMR::interpolate_amr(double* XyzIn_D, int& nCell, PIC::Mesh::cDataCenterNode** Cell, double* Weight){
  //reset the list of found blocks
  iBlockFoundCurrent = 0;
  for(int iBlock = 0; iBlock < nBlockFoundMax; iBlock++)
    BlockFound[iBlock] = NULL;

  //prepare the call of FORTRAN interpolate_amr subroutine
  int nDim       = DIM;
  // number of indices to identify cell ON A GIVEN PROCESSOR:
  // block id + DIM indices of cell in a block
  int nIndexes   = DIM + 1; 
  int nCell_D[3] = {_BLOCK_CELLS_X_,
		    _BLOCK_CELLS_Y_,
		    _BLOCK_CELLS_Z_};
  // IMPORTANT NOTE: iIndexes_II has an extra index per cell which is 
  // a processor ID => need one more index in addition to nIndexes
  int iIndexes_II[(DIM+1+1)*nCellStencilMax] = {-1};
  bool IsSecondOrder;
  bool UseGhostCell=true;

  // call the interpolation subroutine
  INTERFACE__INTERPOLATION_AMR__interpolate_amr_(&nDim, XyzIn_D, &nIndexes, 
						 nCell_D, &nCellStencil, WeightStencil, iIndexes_II, &IsSecondOrder, &UseGhostCell);

  // size of the stencil and weights are known
  // need to identify cells in the stencil
  // NOTE: FORTRAN is column-major for 2D arrays like iIndexes_II
  for(int iCellStencil = 0; iCellStencil < nCellStencil; iCellStencil++){
    int ind[3]  = {0};
    PIC::Mesh::cDataBlockAMR  *block;
    PIC::Mesh::cDataCenterNode *cell;
    int iThread = iIndexes_II[0    +iCellStencil*(nIndexes+1)];
    for(int i = 0; i < DIM; i++)
      ind[i]    = iIndexes_II[1+i  +iCellStencil*(nIndexes+1)];
    int iBlock  = iIndexes_II[1+DIM+iCellStencil*(nIndexes+1)];
    //retrieve the pointer to the current cell
    CellStencil[iCellStencil] = BlockFound[iBlock]->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(ind[0],ind[1],ind[2]));
  }

  //write the output parameters
  nCell  = nCellStencil;
  Weight = WeightStencil;
  Cell   = CellStencil;
}


// procedure to find block containing the point of interest and its parameters
void INTERFACE__INTERPOLATION_AMR__find_(int* nDim, double* Xyz_D, int* iProc, int* iBlock, double* XyzCorner_D, double* Dxyz_D, bool* IsOut){

  //check correctness
  if(DIM != *nDim)
    exit(__LINE__,__FILE__,"Error: inconsistent number of dimensions called by AMR interpolation procedure");

  // find the node containg the point
  static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node = NULL;
  node=PIC::Mesh::mesh.findTreeNode(Xyz_D, node);

  // if node is not found => point is outside of the domain, break
  if(node==NULL) {
    *IsOut=true;
    return;
  }

  //node has been found, proceed
  *iProc=node->Thread;
  INTERFACE::INTERPOLATION_AMR::BlockFound[INTERFACE::INTERPOLATION_AMR::iBlockFoundCurrent++] = node;
  memcpy(XyzCorner_D, node->xmin, DIM*sizeof(double));
  Dxyz_D[0] = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
  Dxyz_D[1] = (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
  Dxyz_D[2] = (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;
  *IsOut=false;
  
  //change coordinates of the point by subtracting the corner coordinates
  for(int iDim = 0; iDim < *nDim; iDim++)
    Xyz_D[iDim] -= XyzCorner_D[iDim];

}
