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


// macro switch is needed in the case some other interpolation is used
// and interface function is not compiled
#if _PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_
//definition of the interface functions that are used to access the FORTRAN written linear cell cenetered interpolation library
extern "C"{
  void interface__cell_centered_linear_interpolation__init_stencil_(int* nDim, double* XyzIn_D, int* nIndexes, int* nCell_D,  int* nGridOut, double* Weight_I, int* iIndexes_II, int* IsSecondOrder, int* UseGhostCell);
}
#endif//_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_


//determine stencil for the cell centered piecewise constant interpolation
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int i,j,k;
  long int nd;
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
  nd = PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false);
  cell=node->block->GetCenterNode(nd);//PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

  //add the cell to the stencil
  if (cell!=NULL) PIC::InterpolationRoutines::CellCentered::Stencil.AddCell(1.0,cell);

  return &PIC::InterpolationRoutines::CellCentered::Stencil;
}

//determine the stencil for the cell centered linear interpolation using interpolation library from ../share/Library/src/
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(double *XyzIn_D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // macro switch is needed in the case some other interpolation is used
  // and interface function is not compiled
#if _PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_

  //find the block if needed
  if (node==NULL) node=PIC::Mesh::mesh.findTreeNode(XyzIn_D);


  //check whether the point is located deep in the block -> use three linear interpolation
  double iLoc,jLoc,kLoc;
  double xmin[3],xmax[3];

  memcpy(xmin,node->xmin,3*sizeof(double));
  memcpy(xmax,node->xmax,3*sizeof(double));

  iLoc=(XyzIn_D[0]-xmin[0])/(xmax[0]-xmin[0])*_BLOCK_CELLS_X_;
  jLoc=(XyzIn_D[1]-xmin[1])/(xmax[1]-xmin[1])*_BLOCK_CELLS_Y_;
  kLoc=(XyzIn_D[2]-xmin[2])/(xmax[2]-xmin[2])*_BLOCK_CELLS_Z_;

  //if the point of interest is deep inside the block or all neighbors of the block has the same resolution level -> use simple trilinear interpolation
  if ( ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) ||
      (0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5) ) {
    return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);
  }

  //re-init variables in the INTERFACE and flush the Stencil
  INTERFACE::iBlockFoundCurrent=0;
  for(int iBlock = 0; iBlock < INTERFACE::nBlockFoundMax; iBlock++ ) INTERFACE::BlockFound[iBlock] = NULL;

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
  int iIndexes_II[(DIM+1+1)*PIC::InterpolationRoutines::CellCentered::nMaxStencilLength];
  double WeightStencil[PIC::InterpolationRoutines::CellCentered::nMaxStencilLength];
  int IsSecondOrder;
  int UseGhostCell=1;
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
    for(int i = 0; i < DIM; i++)
      //cell indices are 1-based in FORTRAN 
      ind[i]    = iIndexes_II[1+i  +iCellStencil*(nIndexes+1)]-1;
    int iBlock  = iIndexes_II[1+DIM+iCellStencil*(nIndexes+1)];

    //retrieve the pointer to the current cell
    if ((block=INTERFACE::BlockFound[iBlock]->block)==NULL) {
      return PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node);
    }

    cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(ind[0],ind[1],ind[2]));

    if (cell==NULL) {
      //there is no enough information to build up linear interpolation -> use constant interpolation
      return PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node);
    }

    PIC::InterpolationRoutines::CellCentered::Stencil.AddCell(WeightStencil[iCellStencil],INTERFACE::BlockFound[iBlock]->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(ind[0],ind[1],ind[2])));
  }
#else
  exit(__LINE__,__FILE__,"ERROR: cell centered linear interpolation is currently available only through interface, add corresponding block to the input file!");
#endif//_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_

  return &PIC::InterpolationRoutines::CellCentered::Stencil;
}

//triliniar interpolation used inside blocks
PIC::InterpolationRoutines::CellCentered::cStencil *PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationStencil(double iLoc,double jLoc,double kLoc,double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;
  if (block==NULL) exit(__LINE__,__FILE__,"Error: the block is node allocated");

  //flush the stencil
  PIC::InterpolationRoutines::CellCentered::Stencil.flush();

  //determine the aray of the cell's pointer that will be used in the interpolation stencil
  double w[3],InterpolationWeight,totalInterpolationWeight=0.0;
  int i,j,k,i0,j0,k0,nd;

  i0=(iLoc<0.5) ? -1 : (int)(iLoc-0.50);
  j0=(jLoc<0.5) ? -1 : (int)(jLoc-0.50);
  k0=(kLoc<0.5) ? -1 : (int)(kLoc-0.50);

  //get coefficients used in determening of the interpolation weights
  w[0]=iLoc-(i0+0.5);
  w[1]=jLoc-(j0+0.5);
  w[2]=kLoc-(k0+0.5);

  for (i=0;i<2;i++) for (j=0;j<2;j++) for (k=0;k<2;k++) {
    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i0+i,j0+j,k0+k);
    cell=block->GetCenterNode(nd);

    if (cell!=NULL) {
      switch (i+2*j+4*k) {
      case 0+0*2+0*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*(1.0-w[2]);
        break;
      case 1+0*2+0*4:
        InterpolationWeight=w[0]*(1.0-w[1])*(1.0-w[2]);
        break;
      case 0+1*2+0*4:
        InterpolationWeight=(1.0-w[0])*w[1]*(1.0-w[2]);
        break;
      case 1+1*2+0*4:
        InterpolationWeight=w[0]*w[1]*(1.0-w[2]);
        break;

      case 0+0*2+1*4:
        InterpolationWeight=(1.0-w[0])*(1.0-w[1])*w[2];
        break;
      case 1+0*2+1*4:
        InterpolationWeight=w[0]*(1.0-w[1])*w[2];
        break;
      case 0+1*2+1*4:
        InterpolationWeight=(1.0-w[0])*w[1]*w[2];
        break;
      case 1+1*2+1*4:
        InterpolationWeight=w[0]*w[1]*w[2];
        break;

      default:
        exit(__LINE__,__FILE__,"Error: the option is not defined");
      }

      PIC::InterpolationRoutines::CellCentered::Stencil.AddCell(InterpolationWeight,cell);
      totalInterpolationWeight+=InterpolationWeight;
    }
  }

  if (PIC::InterpolationRoutines::CellCentered::Stencil.Length==0) {
    //no cell have been found -> use a canstarnt interpolation stencil
    return PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node);
  }
  else if (PIC::InterpolationRoutines::CellCentered::Stencil.Length!=8) {
    //the interpolated stencil containes less that 8 elements -> the interpolation weights have to be renormalized
    for (int i=0;i<PIC::InterpolationRoutines::CellCentered::Stencil.Length;i++) {
      PIC::InterpolationRoutines::CellCentered::Stencil.Weight[i]/=totalInterpolationWeight;
    }
  }

  return &PIC::InterpolationRoutines::CellCentered::Stencil;
}
