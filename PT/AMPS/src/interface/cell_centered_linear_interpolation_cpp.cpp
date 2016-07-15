
//$Id$
//the c-side interface routine for the cell centered liniar interpolation library

#include "pic.h"

// procedure to find block containing the point of interest and its parameters
extern "C"{

  void interface__cell_centered_linear_interpolation__find_cpp_(int* nDim, double* Xyz_D, int* iProc, int* iBlock, double* XyzCorner_D, double* Dxyz_D, int* IsOut) {
    //check correctness
    if(DIM != *nDim)
      exit(__LINE__,__FILE__,"Error: inconsistent number of dimensions called by AMR interpolation procedure");

    //corner and middle coordinates of the first found block
    static double xMinOriginal[3], xMidOriginal[3],xMaxOriginal[3];

    static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    if(PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::iBlockFoundCurrent==0) {
      // find the node containg the point
      node = NULL;
      node=PIC::Mesh::Search::FindBlock(Xyz_D);
      for(int idim=0; idim < DIM; idim++){
	xMinOriginal[idim] = node->xmin[idim];
	xMaxOriginal[idim] = node->xmax[idim];
	xMidOriginal[idim] = 0.5*(xMinOriginal[idim]+xMaxOriginal[idim]);
      }
    }
    else{
      //the block containing the point has enough info to find all subsequent blocks
      // index of displacement with respect to the first found block and type of displacement
      int disp[3]={0}, type=0;
      //compute displacements
      for(int idim=0; idim < DIM; idim++){
	if     (xMaxOriginal[idim] <=Xyz_D[idim]){ disp[idim] = 1; type+= 1;}
	else if(xMinOriginal[idim] > Xyz_D[idim]){ disp[idim] =-1; type+= 1;}
      }
      switch(type)
	{
	case 0:
	  node = PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[0];
	  break;
	case 1:
	  //find which face
	  int nface;
	  if      (disp[0]==-1) nface=0; else if (disp[0]==1) nface=1;
	  else if (disp[1]==-1) nface=2; else if (disp[1]==1) nface=3;
	  else if (disp[2]==-1) nface=4; else if (disp[2]==1) nface=5;
	  //local coordinates: matters if there are 4 neighbors, all the same if just 1
	  int i,j;
	  if     (disp[0]!=0){i=(xMidOriginal[1]>Xyz_D[1])?0:1; j=(xMidOriginal[2]>Xyz_D[2])?0:1;}
	  else if(disp[1]!=0){i=(xMidOriginal[0]>Xyz_D[0])?0:1; j=(xMidOriginal[2]>Xyz_D[2])?0:1;}
	  else if(disp[2]!=0){i=(xMidOriginal[0]>Xyz_D[0])?0:1; j=(xMidOriginal[1]>Xyz_D[1])?0:1;}
	  else exit(__LINE__,__FILE__,"Error: incorrect neighboring block indices");
	  node = PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[0]->GetNeibFace(nface,i,j);
	  break;
	case 2:
	  //find which edge
	  int nedge;
	  if      (disp[0]==-1) {if (disp[1]!=0) nedge=(disp[1]==-1) ? 8 : 11;
 	                         else            nedge=(disp[2]==-1) ? 4 : 7;
	  }
	  else if (disp[0]==1) {if  (disp[1]!=0) nedge=(disp[1]==-1) ? 9 : 10;
	                        else             nedge=(disp[2]==-1) ? 5 : 6;
	  }
	  else if (disp[1]==-1) nedge=(disp[2]==-1) ? 0 : 3;
	  else if (disp[1]== 1) nedge=(disp[2]==-1) ? 1 : 2;

	  //local coordinates: matters if there are 2 neighbors, all the same if just 1
	  if     (disp[0]==0){i=(xMidOriginal[0]>Xyz_D[0])?0:1;}
	  else if(disp[1]==0){i=(xMidOriginal[1]>Xyz_D[1])?0:1;}
	  else if(disp[2]==0){i=(xMidOriginal[2]>Xyz_D[2])?0:1;}
	  else exit(__LINE__,__FILE__,"Error: incorrect neighboring block indices");
	  node = PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[0]->GetNeibEdge(nedge,i);
	  break;
	case 3:
	  //find which corner
	  int ncorner;
	  if(disp[0]>0) ncorner = 1; else ncorner = 0;
	  if(disp[1]>0) ncorner+= 2;
	  if(disp[2]>0) ncorner+= 4;
	  node = PIC::InterpolationRoutines::CellCentered::Linear::INTERFACE::BlockFound[0]->GetNeibCorner(ncorner);
	  break;
	}

    }


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



