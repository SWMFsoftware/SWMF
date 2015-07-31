
//$Id$
//model of the ion source peocesses


/*
 * mars-ions_source.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: vtenishe
 */

#include "mars-ions.h"

double MarsIon::SourceProcesses::GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell) {
  double res=0.0;

  res=1.0;

  return res*cell->Measure;
}


double MarsIon::SourceProcesses::GetBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double res=0.0;
  int i,j,k,LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;

  if (block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
    LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
    cell=block->GetCenterNode(LocalCellNumber);

    if (cell!=NULL) res+=GetCellInjectionRate(spec,cell);
  }

  return res;
}



