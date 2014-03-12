
//$Id$
//the general functions for the "rotation body" internal boundary installed into the mesh


#include "pic.h"

cAMRheap<cInternalRotationBodyData> PIC::BC::InternalBoundary::RotationBody::InternalRotationBody;


//====================================================
//register a new spherical body on the mesh
cInternalBoundaryConditionsDescriptor PIC::BC::InternalBoundary::RotationBody::RegisterInternalRotationBody() {
  cInternalRotationBodyData *newRotationBody;
  cInternalBoundaryConditionsDescriptor RotationBodyDescriptor;


  if (TotalSampleSetLength==0) exit(__LINE__,__FILE__,"Error: need to initialize PIC::BC::InternalBoundary::Sphere first");

  //set the new sphere
  newRotationBody=InternalRotationBody.newElement();

  //set up the output functions
  newRotationBody->PrintVariableList=PrintDefaultVariableList;
  newRotationBody->PrintTitle=PrintDefaultTitle;
  newRotationBody->PrintDataStateVector=PrintDefaultDataStateVector;

  //init the internal parameters of the new sphere
  double *sBuffer;
  long int i,sBufferTotalLength;

  TotalSurfaceElementNumber=newRotationBody->GetTotalSurfaceElementsNumber();
  sBufferTotalLength=2*TotalSurfaceElementNumber*TotalSampleSetLength;


  newRotationBody->faceat=0;
  newRotationBody->SamplingBuffer=new double [sBufferTotalLength];
  newRotationBody->maxIntersectedNodeTimeStep=new double [PIC::nTotalSpecies];

  for (i=0,sBuffer=newRotationBody->SamplingBuffer;i<sBufferTotalLength;i++) sBuffer[i]=0.0;
  for (i=0;i<PIC::nTotalSpecies;i++) newRotationBody->maxIntersectedNodeTimeStep[i]=-1.0;

  //set the descriptor
  RotationBodyDescriptor.BoundaryElement=(void*)newRotationBody;
  RotationBodyDescriptor.BondaryType=_INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_;

  PIC::Mesh::mesh.RegisterInternalBoundary(RotationBodyDescriptor);

  return RotationBodyDescriptor;
}

void PIC::BC::InternalBoundary::RotationBody::PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalRotationBodyData *RotationBody,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {

}

