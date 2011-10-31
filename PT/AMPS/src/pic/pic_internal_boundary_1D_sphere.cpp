//====================================================
//$Id$
//====================================================
//the general functions for the spherical (1D) internal boundary installed into the mesh

#include "pic.h"




int PIC::BC::InternalBoundary::Sphere_1D::completedCellSampleDataPointerOffset=0,PIC::BC::InternalBoundary::Sphere_1D::collectingCellSampleDataPointerOffset=0;

int PIC::BC::InternalBoundary::Sphere_1D::sampledFluxDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere_1D::sampledFluxUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere_1D::sampledMeanVelocityDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere_1D::sampledMeanVelocityUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere_1D::sampledMeanEnergyDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere_1D::sampledMeanEnergyUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere_1D::sampledSurfaceNumberDensityRelativeOffset=0;

long int PIC::BC::InternalBoundary::Sphere_1D::TotalSampleSetLength=0;
long int *PIC::BC::InternalBoundary::Sphere_1D::SpeciesSampleDataOffset=NULL;
long int *PIC::BC::InternalBoundary::Sphere_1D::SpeciesSampleUserDefinedDataOffset=NULL;
long int PIC::BC::InternalBoundary::Sphere_1D::TotalSurfaceElementNumber=0;

bool PIC::BC::InternalBoundary::Sphere_1D::UserDefinedSamplingProcedureFlag=false;
cAMRheap<cInternalSphere1DData> PIC::BC::InternalBoundary::Sphere_1D::InternalSpheres;

PIC::BC::InternalBoundary::Sphere_1D::fPrintVariableList PIC::BC::InternalBoundary::Sphere_1D::PrintUserDefinedVariableList=NULL;
PIC::BC::InternalBoundary::Sphere_1D::fPrintTitle PIC::BC::InternalBoundary::Sphere_1D::PrintUserDefinedTitle=NULL;
PIC::BC::InternalBoundary::Sphere_1D::fPrintDataStateVector PIC::BC::InternalBoundary::Sphere_1D::PrintUserDefinedDataStateVector=NULL;
PIC::BC::InternalBoundary::Sphere_1D::fSampleParticleData PIC::BC::InternalBoundary::Sphere_1D::SampleParticleData=NULL;


//====================================================
//init the data for the spherical body installed into the mesh:
void PIC::BC::InternalBoundary::Sphere_1D::Init(long int* RequestedSamplingSetDataLength,long int* UserDefinedSampleDataRelativeOffset) {

  if (TotalSampleSetLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of PIC::BC::InternalBoundary::Sphere_1D:");

  //init the sampling offsets
  int DefaultDataLength=0;

  SpeciesSampleDataOffset=new long int[PIC::nTotalSpecies];
  SpeciesSampleUserDefinedDataOffset=new long int[PIC::nTotalSpecies];


  //set the offsets to the default sampled data
  sampledFluxDownRelativeOffset=DefaultDataLength;
  DefaultDataLength+=1;

  sampledFluxUpRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  sampledMeanVelocityDownRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  sampledMeanVelocityUpRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  sampledMeanEnergyDownRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  sampledMeanEnergyUpRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  sampledSurfaceNumberDensityRelativeOffset+=DefaultDataLength;
  DefaultDataLength+=1;

  //set offsets to the user defined data
  for (int s=0;s<PIC::nTotalSpecies;s++) {
    SpeciesSampleDataOffset[s]=TotalSampleSetLength;
    TotalSampleSetLength+=DefaultDataLength;

    SpeciesSampleUserDefinedDataOffset[s]=TotalSampleSetLength;
    if (UserDefinedSampleDataRelativeOffset!=NULL) UserDefinedSampleDataRelativeOffset[s]=TotalSampleSetLength;
    if (RequestedSamplingSetDataLength!=NULL) TotalSampleSetLength+=RequestedSamplingSetDataLength[s];
  }


  //init the offsets for the completed and collecting samples
  completedCellSampleDataPointerOffset=0,collectingCellSampleDataPointerOffset=TotalSampleSetLength;
}

void PIC::BC::InternalBoundary::Sphere_1D::Init() {
  Init(NULL,NULL);
}

//====================================================
//register a new spherical body on the mesh
cInternalBoundaryConditionsDescriptor PIC::BC::InternalBoundary::Sphere_1D::RegisterInternalSphere() {
  cInternalSphere1DData *newSphere;
  cInternalBoundaryConditionsDescriptor SphereDescriptor;
//  cSurfaceDataSphere_1D* SurfaceData;



  if (TotalSampleSetLength==0) exit(__LINE__,__FILE__,"Error: need to initialize PIC::BC::InternalBoundary::Sphere_1D first");

  //set the new Sphere_1D
  newSphere=InternalSpheres.newElement();
//  SurfaceData=(cSurfaceDataSphere_1D*)malloc(sizeof(cSurfaceDataSphere_1D));
//  newSphere_1D->SetSurfaceDataPointer((void*)SurfaceData);

  //set up the output functions
  newSphere->PrintVariableList=PrintDefaultVariableList;
  newSphere->PrintTitle=PrintDefaultTitle;
  newSphere->PrintDataStateVector=PrintDefaultDataStateVector;

  //init the internal parameters of the new Sphere_1D
  double *sBuffer;
  long int i,sBufferTotalLength;

  TotalSurfaceElementNumber=newSphere->GetTotalSurfaceElementsNumber();
  sBufferTotalLength=2*TotalSurfaceElementNumber*TotalSampleSetLength;


  newSphere->faceat=0;
  newSphere->SamplingBuffer=new double [sBufferTotalLength];
  newSphere->maxIntersectedNodeTimeStep=new double [PIC::nTotalSpecies];

  for (i=0,sBuffer=newSphere->SamplingBuffer;i<sBufferTotalLength;i++) sBuffer[i]=0.0;
  for (i=0;i<PIC::nTotalSpecies;i++) newSphere->maxIntersectedNodeTimeStep[i]=-1.0;

  //set the descriptor
  SphereDescriptor.BoundaryElement=(void*)newSphere;
  SphereDescriptor.BondaryType=_INTERNAL_BOUNDARY_TYPE_1D_SPHERE_;

  PIC::Mesh::mesh.RegisterInternalBoundary(SphereDescriptor);

  return SphereDescriptor;
}


//====================================================
//clear sampling buffer
void PIC::BC::InternalBoundary::Sphere_1D::flushCollectingSamplingBuffer(cInternalSphericalData* Sphere) {
  long int offset,i,nSurfaceElement,nSurfaceElementsMax;
  double *SamplingBuffer;

//  SamplingBuffer=((cSurfaceDataSphere_1D*)Sphere_1D->GetSurfaceDataPointer())->SamplingBuffer;
  SamplingBuffer=Sphere->SamplingBuffer;
  nSurfaceElementsMax=Sphere->GetTotalSurfaceElementsNumber();

  for (nSurfaceElement=0;nSurfaceElement<nSurfaceElementsMax;nSurfaceElement++) {
    offset=collectingSpecieSamplingDataOffset(0,nSurfaceElement);

    for (i=0;i<TotalSampleSetLength;i++) *(SamplingBuffer+offset+i)=0.0;
  }
}

//====================================================
//get access to the internal data of the Sphere_1D
/*
PIC::BC::InternalBoundary::Sphere_1D::cSurfaceDataSphere_1D* PIC::BC::InternalBoundary::Sphere_1D::GetSphere_1DSurfaceData(cInternalBoundaryConditionsDescriptor Descriptor) {
  if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_Sphere_1D_) exit(__LINE__,__FILE__,"Error: wrong boundary type");

  return (cSurfaceDataSphere_1D*)(((cInternalSphericalData*)Descriptor.BoundaryElement)->GetSurfaceDataPointer());

}
*/

double* PIC::BC::InternalBoundary::Sphere_1D::GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
//  return GetSphere_1DSurfaceData(Descriptor)->SamplingBuffer+completedCellSampleDataPointerOffset;


  return ((cInternalSphericalData*)Descriptor.BoundaryElement)->SamplingBuffer+completedCellSampleDataPointerOffset;
}


double* PIC::BC::InternalBoundary::Sphere_1D::GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
//  return GetSphere_1DSurfaceData(Descriptor)->SamplingBuffer+collectingCellSampleDataPointerOffset;

  return ((cInternalSphericalData*)Descriptor.BoundaryElement)->SamplingBuffer+collectingCellSampleDataPointerOffset;
}

//====================================================
//the offset of sampling data for a particular specie
int PIC::BC::InternalBoundary::Sphere_1D::completeSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  if (SurfaceElement!=0) exit(__LINE__,__FILE__,"Error: out of range");

  return 2*SurfaceElement*TotalSampleSetLength+completedCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Sphere_1D::collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  if (SurfaceElement!=0) exit(__LINE__,__FILE__,"Error: out of range");

  return 2*SurfaceElement*TotalSampleSetLength+collectingCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Sphere_1D::SurfaceElementSamplingSetLength() {
  return TotalSampleSetLength;
}

void PIC::BC::InternalBoundary::Sphere_1D::switchSamplingBuffers() {
  long int tempSamplingOffset;

  tempSamplingOffset=completedCellSampleDataPointerOffset;
  completedCellSampleDataPointerOffset=collectingCellSampleDataPointerOffset;
  collectingCellSampleDataPointerOffset=tempSamplingOffset;
}

//====================================================
//particle-spherical surface interaction
//specular reflection of the particle velocity vector on the surface of the Sphere_1D
int PIC::BC::InternalBoundary::Sphere_1D::ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPointer,void* SphereDataPointer)  {
   double radiusSphere,*x0Sphere,l[3],r,vNorm,c;
   cInternalSphericalData *Sphere;
   cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*startNode;
//   cSurfaceDataSphere_1D *SurfaceData;
   int idim;


   Sphere=(cInternalSphericalData*)SphereDataPointer;
   startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPointer;


   #if DIM != 3
   exit(__LINE__,__FILE__,"Error: wrong dimension");
   #endif


   Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

   for (r=0.0,idim=0;idim<DIM;idim++) {
     l[idim]=x[idim]-x0Sphere[idim];
     r+=pow(l[idim],2);
   }

   for (r=sqrt(r),vNorm=0.0,idim=0;idim<DIM;idim++) vNorm+=v[idim]*l[idim]/r;
   for (c=2.0*vNorm/r,idim=0;idim<DIM;idim++) v[idim]-=c*l[idim];

   //sample the particle data
   double *SampleData;
   long int nSurfaceElement,nZenithElement,nAzimuthalElement;

   Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
   nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
//   SurfaceData=Sphere_1D->SamplingBuffer;
   SampleData=Sphere->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


///###########  DEBUG #############

//   SampleData=SurfaceData->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,0);

/*
   if (x[2]>0.20) {
     cout << __LINE__ << endl;
     Sphere_1D->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
     nSurfaceElement=Sphere_1D->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

   }
*/

//########### END DEBUG #######





   SampleData[sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);
   return _PARTICLE_REJECTED_ON_THE_FACE_;
}


//====================================================
//print default surface data
void PIC::BC::InternalBoundary::Sphere_1D::PrintDefaultVariableList(FILE* fout) {
  fprintf(fout,", \"Flux Down\"");
}
void PIC::BC::InternalBoundary::Sphere_1D::PrintDefaultTitle(FILE* fout) {
  fprintf(fout,"TITPLE=\"SurfaceData\"");
}

void PIC::BC::InternalBoundary::Sphere_1D::PrintDefaultDataStateVector(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  double *SampleData;
  double FluxDown=0.0;

  SampleData=Sphere->SamplingBuffer+completeSpecieSamplingDataOffset(spec,0);

  FluxDown+=SampleData[sampledFluxDownRelativeOffset];





  if (ThisThread==0)  {

    //collect sampled data from all processors
    for (int thread=1;thread<nTotalThreads;thread++) {
      FluxDown+=pipe->recv<double>(thread);
    }



    fprintf(fout," %e ",FluxDown);
  }
  else pipe->send(FluxDown);
}


/*
//Sampling of the particles data

void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* Sphere_1DDescriptor);

*/
