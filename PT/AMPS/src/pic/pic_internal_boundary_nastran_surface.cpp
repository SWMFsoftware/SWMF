
//$Id$
//the general functions for the "nastran surface" internal boundary installed into the mesh


#include "pic.h"

cAMRheap<cInternalNastranSurfaceData> PIC::BC::InternalBoundary::NastranSurface::InternalNastranSurface;

//====================================================
//register a new nastran body on the mesh
cInternalBoundaryConditionsDescriptor PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface() {
  cInternalNastranSurfaceData *newNastranSurface;
  cInternalBoundaryConditionsDescriptor NastranSurfaceDescriptor;
  int i;


//  if (TotalSampleSetLength==0) exit(__LINE__,__FILE__,"Error: need to initialize PIC::BC::InternalBoundary::Sphere first");

  //set the new sphere
  newNastranSurface=InternalNastranSurface.newElement();

  //set up the output functions
  newNastranSurface->PrintVariableList=PrintDefaultVariableList;
  newNastranSurface->PrintTitle=PrintDefaultTitle;
  newNastranSurface->PrintDataStateVector=PrintDefaultDataStateVector;

  newNastranSurface->localResolution=NULL;

  //init the internal parameters of the new sphere
  newNastranSurface->faceat=0;
  newNastranSurface->maxIntersectedNodeTimeStep=new double [PIC::nTotalSpecies];

  for (i=0;i<PIC::nTotalSpecies;i++) newNastranSurface->maxIntersectedNodeTimeStep[i]=-1.0;

  //set the descriptor
  NastranSurfaceDescriptor.BoundaryElement=(void*)newNastranSurface;
  NastranSurfaceDescriptor.BondaryType=_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_;

  PIC::Mesh::mesh.RegisterInternalBoundary(NastranSurfaceDescriptor);

  return NastranSurfaceDescriptor;
}

void PIC::BC::InternalBoundary::NastranSurface::PrintDefaultDataStateVector(FILE* fout,long int nElement,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalNastranSurfaceData *NastranSurface,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {

}

