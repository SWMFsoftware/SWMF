/*
 * OH_Coupler.cpp
 *
 *  Created on: Mar 20, 2015
 *      Author: vtenishe
 */
//$Id$
//Coupler functions for OH

#include "OH.h"
double OH::Coupling::TimeAfterCoupling[PIC::nTotalSpecies] = {0.0};

void OH::Coupling::Send(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI) {
  int i0=0,i1=0;
  char vname[200];

  //get the offsets
  int Density_AMPS2OH  =-1;
  int MomentumX_AMPS2OH=-1;
  int MomentumY_AMPS2OH=-1;
  int MomentumZ_AMPS2OH=-1;
  int Energy_AMPS2OH   =-1;

  while ((NameVar[i0]!=0)&&(NameVar[i0]==' ')) i0++;

  for (int n=0;n<(*nVarIn);n++) {
    i1=i0;
    while (NameVar[i1]!=' ') {
      vname[i1-i0]=tolower(NameVar[i1]);
      i1++;
    }

    vname[i1-i0]=0;

    if (strcmp(vname,"srho")==0) Density_AMPS2OH  =n;
    if (strcmp(vname,"smx" )==0) MomentumX_AMPS2OH=n;
    if (strcmp(vname,"smy" )==0) MomentumY_AMPS2OH=n;
    if (strcmp(vname,"smz" )==0) MomentumZ_AMPS2OH=n;
    if (strcmp(vname,"se" ) ==0) Energy_AMPS2OH   =n;
    
    i0=i1;

    while ((NameVar[i0]!=0)&&(NameVar[i0]==' ')) i0++;
  }

  //save the data
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;
  PIC::Mesh::cDataCenterNode *cell;
  char *AssociatedDataPointer;
  int i,j,k,LocalCellNumber;

  for (int pt=0;pt<(*nPoint);pt++) {
    double x[3]={0.0,0.0,0.0}; //the location of the point; to be sure that it has all 3 components

    //find the block
    memcpy(x,Xyz_DI+pt*(*nDimIn),(unsigned int)(*nDimIn)*sizeof(double));
    if ((node=PIC::Mesh::mesh.findTreeNode(x,node))==NULL) exit(__LINE__,__FILE__,"Error: can not find the block");
    if (node->Thread!=PIC::ThisThread) exit(__LINE__,__FILE__,"Error: the point data is located on another processor");

    //find the cell
    if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");

    cell=node->block->GetCenterNode(LocalCellNumber);
    AssociatedDataPointer=cell->GetAssociatedDataBufferPointer();

    if (PIC::LastSampleLength!=0) {
      // copy density source
      Data_VI[Density_AMPS2OH+pt*(*nVarIn)]=
          (*(0+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceDensityOffset)))/PIC::LastSampleLength;
      // copy momentum source
      Data_VI[MomentumX_AMPS2OH+pt*(*nVarIn)]=
          (*(0+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceMomentumOffset)))/PIC::LastSampleLength;
      Data_VI[MomentumY_AMPS2OH+pt*(*nVarIn)]=
          (*(1+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceMomentumOffset)))/PIC::LastSampleLength;
      Data_VI[MomentumZ_AMPS2OH+pt*(*nVarIn)]=
          (*(2+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceMomentumOffset)))/PIC::LastSampleLength;
      // copy energy source
      Data_VI[Energy_AMPS2OH+pt*(*nVarIn)]=
          (*(0+(double*)(AssociatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceEnergyOffset)))/PIC::LastSampleLength;
    }
    else {
      Data_VI[Density_AMPS2OH  +pt*(*nVarIn)]=0.0;
      Data_VI[MomentumX_AMPS2OH+pt*(*nVarIn)]=0.0;
      Data_VI[MomentumY_AMPS2OH+pt*(*nVarIn)]=0.0;
      Data_VI[MomentumZ_AMPS2OH+pt*(*nVarIn)]=0.0;
      Data_VI[Energy_AMPS2OH   +pt*(*nVarIn)]=0.0;
    }

  }

  // reset time after last coupling session
  for(int spec=0; spec < PIC::nTotalSpecies; spec++)
    OH::Coupling::TimeAfterCoupling[spec] = 0.0;
}



