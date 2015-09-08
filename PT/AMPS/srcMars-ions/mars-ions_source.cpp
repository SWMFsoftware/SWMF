
//$Id$
//model of the ion source peocesses


/*
 * mars-ions_source.cpp
 *
 *  Created on: Jul 31, 2015
 *      Author: vtenishe
 */

#include "mars-ions.h"

double MarsIon::SourceProcesses::GetCellInjectionRate(int spec,double *xMiddle) {
  double res=0.0;
  double altitude=0.0;
  double rSeason=1.38758;
  double r2Season=pow(rSeason,2); 
  double Tnu_body = 134.0; //in K, neutral temperature
  double BodynDenNuSp_I_O= 8.0283e15; // per m^-3
  double BodynDenNuSp_I_Ox= 5.1736e14;
  double BodynDenNuSp_I_Oh= 6.3119e10;
  double BodynDenNuSp_I_Ohx= 3.9646e9;
  double HNuSpecies_I_O=13340; //scale height in m
  double HNuSpecies_I_Ox=50025;
  double HNuSpecies_I_Oh=290500;
  double HNuSpecies_I_Ohx=2436600;
  double Rate_I_O_Op= 6.346e-7/r2Season; //units s^-1 for O_hv__Op_em_
  double nDen_O;
  int i;    
  for (i=0;i<3;i++) {    
      altitude+=pow(xMiddle[i],2);
  }
  altitude=sqrt(altitude)-3396000.0; // in m where R_Mars=3396000.0 m

  if (altitude<0.0) return 0.0;

  nDen_O=BodynDenNuSp_I_O*exp(-altitude/HNuSpecies_I_O)+\
         BodynDenNuSp_I_Ox*exp(-altitude/HNuSpecies_I_Ox)+\
         BodynDenNuSp_I_Oh*exp(-altitude/HNuSpecies_I_Oh)+\
         BodynDenNuSp_I_Ohx*exp(-altitude/HNuSpecies_I_Ohx);
  //res=1.0;
  res=nDen_O*Rate_I_O_Op; //source rate per m^-3

  return res;
}

double MarsIon::SourceProcesses::GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell) {
  double res,*xMiddle;

  xMiddle=cell->GetX();
  res=GetCellInjectionRate(spec,xMiddle);

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



