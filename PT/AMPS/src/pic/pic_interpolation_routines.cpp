//$Id$
//interpolation routines

/*
 * pic_interpolation_routines.cpp
 *
 *  Created on: Jul 18, 2015
 *      Author: vtenishe
 */

#include "pic.h"

PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::StencilTable=NULL;
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


//initialize the interpolation module
void PIC::InterpolationRoutines::Init() {

  //init the stencil table
  CellCentered::StencilTable=new CellCentered::cStencil[PIC::nTotalThreadsOpenMP];
}

//determine stencil for the cell centered piecewise constant interpolation
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int i,j,k;
  long int nd;
  PIC::Mesh::cDataCenterNode *cell;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].flush();

  //find the block
  if (node==NULL) {
    node=PIC::Mesh::mesh.findTreeNode(x);

    if (node==NULL) exit(__LINE__,__FILE__,"Error: the location is outside of the computational domain");
    if (node->block==NULL) {
      char msg[200];

      sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
      exit(__LINE__,__FILE__,msg);
    }
  }
  else if (node->block==NULL) {
    char msg[200];

    sprintf(msg,"Error: the block is not allocated\nCurrent MPI Process=%i\nnode->Thread=%i\n",PIC::ThisThread,node->Thread);
    exit(__LINE__,__FILE__,msg);
  }

  //find cell
  nd = PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false);
  cell=node->block->GetCenterNode(nd);//PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

  //add the cell to the stencil
  if (cell!=NULL) {
    PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].AddCell(1.0,cell);
  }
  else exit(__LINE__,__FILE__,"Error: cell is not initialized");

  return PIC::InterpolationRoutines::CellCentered::StencilTable+ThreadOpenMP;
}

//determine the stencil for the cell centered linear interpolation using interpolation library from ../share/Library/src/
PIC::InterpolationRoutines::CellCentered::cStencil* PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(double *XyzIn_D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();

  if  (_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_) {
    exit(__LINE__,__FILE__,"Error: the procedure is not adapted fro using with OpenMP: INTERFACE::BlockFound... need to be converted into an array as PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP] ");
  }
  #else
  int ThreadOpenMP=0;
  #endif



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

  #if _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
  //if the point of interest is deep inside the block or all neighbors of the block has the same resolution level -> use simple trilinear interpolation
  if ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) {
    //all blocks around has the same level of the resolution
    return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);
  }
  else if ((1.0<iLoc)&&(iLoc<_BLOCK_CELLS_X_-1) && (1.0<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-1) && (1.0<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-1)) {
    //the point is deep inside the block
    return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);
  }
  else if (node->RefinmentLevel==node->minNeibRefinmentLevel) {
    if  ((0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5)) {
      //1. the point is deeper than a half cell size insode the block
      //2. the block is geometrically largest among the neibours
      //3. -> it is only the intermal points of the block that will be used in the stencil
      return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);
    }
    else  {
      //the block is largest between the neibours
      //getermine the size limit for the interpolation stencil
      double xStencilMin[3],xStencilMax[3];
      double dxCell[3];

      dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

      for (int idim=0;idim<3;idim++) {
        int iInterval;

        iInterval=(int)((XyzIn_D[idim]-(xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
        xStencilMin[idim]=(xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
        xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
      }

      return GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,node);
    }
  }
  else {
    //1. find a coarser block that is close to the point, and can be used for interpolation
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *CoarserBlock=NULL;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *NeibNode;
    int idim,iFace;
    double dxCell[3];

    dxCell[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
    dxCell[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
    dxCell[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

    for (idim=0;idim<3;idim++) if (CoarserBlock==NULL) {
      switch (idim) {
      case 0:
        if (iLoc<=1.0) iFace=0;
        else if (iLoc>=_BLOCK_CELLS_X_-1.0) iFace=1;
        else continue;

        break;
      case 1:
        if (jLoc<=1.0) iFace=2;
        else if (jLoc>=_BLOCK_CELLS_Y_-1.0) iFace=3;
        else continue;

        break;
      case 2:
        if (kLoc<=1.0) iFace=4;
        else if (kLoc>=_BLOCK_CELLS_Z_-1.0) iFace=5;
        else continue;
      }

      //check blocks connected through the face
      NeibNode=node->GetNeibFace(iFace,0,0);

      if (NeibNode!=NULL) if (NeibNode->RefinmentLevel<node->RefinmentLevel) {
        //found a coarser block
        int cnt=0;

        for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

        if (cnt==3) {
          //the block can be used for interpolation
           CoarserBlock=NeibNode;
        }
      }

      //check blocks connected though the edges
      static const int faceEdges[6][4]={{4,11,7,8},{5,10,6,9},{0,9,3,8},{1,10,2,11},{0,5,1,4},{3,6,2,7}};

      if (CoarserBlock==NULL) for (int iEdge=0;iEdge<4;iEdge++) {
        NeibNode=node->GetNeibEdge(faceEdges[iFace][iEdge],0);

        if (NeibNode!=NULL) if (NeibNode->RefinmentLevel<node->RefinmentLevel) {
          //found a coarser block
          int cnt=0;

          for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

          if (cnt==3) {
            //the block can be used for interpolation
            CoarserBlock=NeibNode;
            break;
          }
        }
      }

      //check blocks connected through the corners
      static const int FaceNodeMap[6][4]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};

      if (CoarserBlock==NULL) for (int iCorner=0;iCorner<4;iCorner++) {
        NeibNode=node->GetNeibCorner(FaceNodeMap[iFace][iCorner]);

        if (NeibNode!=NULL) if (NeibNode->RefinmentLevel<node->RefinmentLevel) {
          //found a coarser block
          int cnt=0;

          for (int ii=0;ii<3;ii++) if ((NeibNode->xmin[ii]-dxCell[ii]<=XyzIn_D[ii]) && (NeibNode->xmax[ii]+dxCell[ii]>=XyzIn_D[ii]) ) cnt++;

          if (cnt==3) {
            //the block can be used for interpolation
            CoarserBlock=NeibNode;
            break;
          }
        }
      }
    }

    if (CoarserBlock!=NULL) {
      //getermine the size limit for the interpolation stencil
      double xStencilMin[3],xStencilMax[3];
      double dxCell[3];

      dxCell[0]=(CoarserBlock->xmax[0]-CoarserBlock->xmin[0])/_BLOCK_CELLS_X_;
      dxCell[1]=(CoarserBlock->xmax[1]-CoarserBlock->xmin[1])/_BLOCK_CELLS_Y_;
      dxCell[2]=(CoarserBlock->xmax[2]-CoarserBlock->xmin[2])/_BLOCK_CELLS_Z_;

      for (idim=0;idim<3;idim++) {
        int iInterval;

        iInterval=(int)((XyzIn_D[idim]-(xmin[idim]-dxCell[idim]/2.0))/dxCell[idim]);
        xStencilMin[idim]=(xmin[idim]-dxCell[idim]/2.0)+iInterval*dxCell[idim];
        xStencilMax[idim]=xStencilMin[idim]+dxCell[idim];
      }

      return GetTriliniarInterpolationMutiBlockStencil(XyzIn_D,xStencilMin,xStencilMax,CoarserBlock);
    }

    //2.threre is no a coarser block that can be used for interpolation
    // check if all cells are available for tri-liniar interpolation using the data from the block
    int i,j,k,i0,j0,k0,nd;
    bool found=true;

    i0=(iLoc<0.5) ? -1 : (int)(iLoc-0.50);
    j0=(jLoc<0.5) ? -1 : (int)(jLoc-0.50);
    k0=(kLoc<0.5) ? -1 : (int)(kLoc-0.50);

    for (i=0;(i<2)&&(found==true);i++) for (j=0;(j<2)&&(found==true);j++) for (k=0;(k<2)&&(found==true);k++) {
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i0+i,j0+j,k0+k);
      if (node->block->GetCenterNode(nd)==NULL) found=false;
    }

    if (found==true) return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);

    //3. if not: use a constant interpolation spencil
    return PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(XyzIn_D,node);
  }

  #elif _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ == _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__SWMF_
  if ( ((node->RefinmentLevel==node->minNeibRefinmentLevel)&&(node->RefinmentLevel==node->maxNeibRefinmentLevel)) ||
      (0.5<iLoc)&&(iLoc<_BLOCK_CELLS_X_-0.5) && (0.5<jLoc)&&(jLoc<_BLOCK_CELLS_Y_-0.5) && (0.5<kLoc)&&(kLoc<_BLOCK_CELLS_Z_-0.5) ) {
    return GetTriliniarInterpolationStencil(iLoc,jLoc,kLoc,XyzIn_D,node);
  }

  // if the point of interest is very close to the cell center
  // then truncated stencil to the single point
  if( fabs(iLoc-0.5-(long)iLoc) < PrecisionCellCenter &&
      fabs(jLoc-0.5-(long)jLoc) < PrecisionCellCenter &&
      fabs(kLoc-0.5-(long)kLoc) < PrecisionCellCenter   ){
    return PIC::InterpolationRoutines::CellCentered::
      Constant::InitStencil(XyzIn_D,node);
  }

  //re-init variables in the INTERFACE and flush the Stencil
  INTERFACE::iBlockFoundCurrent=0;
  for(int iBlock = 0; iBlock < INTERFACE::nBlockFoundMax; iBlock++ ) INTERFACE::BlockFound[iBlock] = NULL;

  PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].flush();
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

    PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].AddCell(WeightStencil[iCellStencil],INTERFACE::BlockFound[iBlock]->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(ind[0],ind[1],ind[2])));
  }
  #else  //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
  exit(__LINE__,__FILE__,"Error: the option is unknown");
  #endif //_PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE_ _PIC_CELL_CENTERED_LINEAR_INTERPOLATION_ROUTINE__AMPS_
#else
  exit(__LINE__,__FILE__,"ERROR: cell centered linear interpolation is currently available only through interface, add corresponding block to the input file!");
#endif//_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_


  return PIC::InterpolationRoutines::CellCentered::StencilTable+ThreadOpenMP;
}

//triliniar interpolation used inside blocks
PIC::InterpolationRoutines::CellCentered::cStencil *PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationStencil(double iLoc,double jLoc,double kLoc,double *x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;
  if (block==NULL) exit(__LINE__,__FILE__,"Error: the block is node allocated");

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].flush();

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

      PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].AddCell(InterpolationWeight,cell);
      totalInterpolationWeight+=InterpolationWeight;
    }
  }

  if (PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Length==0) {
    //no cell have been found -> use a canstarnt interpolation stencil
    return PIC::InterpolationRoutines::CellCentered::Constant::InitStencil(x,node);
  }
  else if (PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Length!=8) {
    //the interpolated stencil containes less that 8 elements -> the interpolation weights have to be renormalized
    for (int i=0;i<PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Length;i++) {
      PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].Weight[i]/=totalInterpolationWeight;
    }
  }

  return PIC::InterpolationRoutines::CellCentered::StencilTable+ThreadOpenMP;
}



PIC::InterpolationRoutines::CellCentered::cStencil *PIC::InterpolationRoutines::CellCentered::Linear::GetTriliniarInterpolationMutiBlockStencil(double *x,double *xStencilMin,double *xStencilMax,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int iStencil,jStencil,kStencil,i,j,k,nd,idim;
  double xLoc[3],xStencil[3],dx[3];

  const int nStencilElementsMax=64;
  double StencilWeight[nStencilElementsMax],StencilElementWeight,summStencilElementWeight=0.0;
  PIC::Mesh::cDataCenterNode *cell,*StencilCellTable[nStencilElementsMax];
  int StencilElementCounter=0;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *StencilNode=node;

  //calculate the local coordinates of the point where the stancil is constructed
  for (idim=0;idim<3;idim++) {
    xLoc[idim]=(x[idim]-xStencilMin[idim])/(xStencilMax[idim]-xStencilMin[idim]);

    if ((xLoc[idim]<0.0)||(xLoc[idim]>1.0)) exit(__LINE__,__FILE__,"Error: the local coordinate is out or range");
  }

  //cell sizes
  dx[0]=(xStencilMax[0]-xStencilMin[0]);
  dx[1]=(xStencilMax[1]-xStencilMin[1]);
  dx[2]=(xStencilMax[2]-xStencilMin[2]);


  //build interpolation stencil
  for (iStencil=0;iStencil<2;iStencil++) {
    xStencil[0]=xStencilMin[0]+iStencil*(xStencilMax[0]-xStencilMin[0]);

    for (jStencil=0;jStencil<2;jStencil++) {
      xStencil[1]=xStencilMin[1]+jStencil*(xStencilMax[1]-xStencilMin[1]);

      for (kStencil=0;kStencil<2;kStencil++) {
        xStencil[2]=xStencilMin[2]+kStencil*(xStencilMax[2]-xStencilMin[2]);

        StencilNode=PIC::Mesh::mesh.findTreeNode(xStencil,StencilNode);

        switch (iStencil+2*jStencil+4*kStencil) {
        case 0+0*2+0*4:
          StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*(1.0-xLoc[2]);
          break;
        case 1+0*2+0*4:
          StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*(1.0-xLoc[2]);
          break;
        case 0+1*2+0*4:
          StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*(1.0-xLoc[2]);
          break;
        case 1+1*2+0*4:
          StencilElementWeight=xLoc[0]*xLoc[1]*(1.0-xLoc[2]);
          break;

        case 0+0*2+1*4:
          StencilElementWeight=(1.0-xLoc[0])*(1.0-xLoc[1])*xLoc[2];
          break;
        case 1+0*2+1*4:
          StencilElementWeight=xLoc[0]*(1.0-xLoc[1])*xLoc[2];
          break;
        case 0+1*2+1*4:
          StencilElementWeight=(1.0-xLoc[0])*xLoc[1]*xLoc[2];
          break;
        case 1+1*2+1*4:
          StencilElementWeight=xLoc[0]*xLoc[1]*xLoc[2];
          break;

        default:
          exit(__LINE__,__FILE__,"Error: the option is not defined");
        }

        //determine contribution of the node to the stencil
        if (StencilNode->RefinmentLevel==node->RefinmentLevel) {
          //both blocks have the save refinment levels

          nd=PIC::Mesh::mesh.fingCellIndex(xStencil,i,j,k,StencilNode,false);
          cell=StencilNode->block->GetCenterNode(nd);//PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) {
            //add the cell to the stencil
            if (StencilElementCounter==nStencilElementsMax) exit(__LINE__,__FILE__,"Error: StencilElementCounter==nStencilElementsMax, try to increase nStencilElementsMax");
            StencilCellTable[StencilElementCounter]=cell;
            StencilWeight[StencilElementCounter]=StencilElementWeight;
            StencilElementCounter++;
            summStencilElementWeight+=StencilElementWeight;
          }
        }
        else if (StencilNode->RefinmentLevel>node->RefinmentLevel) {
          //get left coordinates of the cell block that will be used as a part of the stencil
          int iCellNeib,jCellNeib,kCellNeib,ii,jj,kk;

          iCellNeib=2*((int)((xStencil[0]-StencilNode->xmin[0])/dx[0]));
          jCellNeib=2*((int)((xStencil[1]-StencilNode->xmin[1])/dx[1]));
          kCellNeib=2*((int)((xStencil[2]-StencilNode->xmin[2])/dx[2]));

          for (ii=0;ii<2;ii++) for (jj=0;jj<2;jj++) for (kk=0;kk<2;kk++) {
            nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(ii+iCellNeib,jj+jCellNeib,kk+kCellNeib);
            cell=StencilNode->block->GetCenterNode(nd);//PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

            if (cell!=NULL) {
              //add the cell to the stencil
              if (StencilElementCounter==nStencilElementsMax) exit(__LINE__,__FILE__,"Error: StencilElementCounter==nStencilElementsMax, try to increase nStencilElementsMax");
              StencilCellTable[StencilElementCounter]=cell;
              StencilWeight[StencilElementCounter]=StencilElementWeight/8.0;
              summStencilElementWeight+=StencilElementWeight/8.0;
              StencilElementCounter++;
            }
          }

        }
        else {
          exit(__LINE__,__FILE__,"Error: something is wrong. The conditions StencilNode->RefinmentLevel>=node->RefinmentLevel must hold");
        }

      }

    }
  }

  //construct the stencil
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int ThreadOpenMP=omp_get_thread_num();
  #else
  int ThreadOpenMP=0;
  #endif

  //flush the stencil
  PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].flush();

  for (int iCell=0;iCell<StencilElementCounter;iCell++) {
    PIC::InterpolationRoutines::CellCentered::StencilTable[ThreadOpenMP].AddCell(StencilWeight[iCell]/summStencilElementWeight,StencilCellTable[iCell]);
  }


  return PIC::InterpolationRoutines::CellCentered::StencilTable+ThreadOpenMP;
}


