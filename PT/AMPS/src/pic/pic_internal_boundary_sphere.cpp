//====================================================
//$Id$
//====================================================
//the general functions for the spherical internal boundary installed into the mesh

#include "pic.h"




int PIC::BC::InternalBoundary::Sphere::completedCellSampleDataPointerOffset=0,PIC::BC::InternalBoundary::Sphere::collectingCellSampleDataPointerOffset=0;

int PIC::BC::InternalBoundary::Sphere::sampledFluxDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere::sampledFluxUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere::sampledMeanVelocityDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere::sampledMeanVelocityUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere::sampledMeanEnergyDownRelativeOffset=0,PIC::BC::InternalBoundary::Sphere::sampledMeanEnergyUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Sphere::sampledSurfaceNumberDensityRelativeOffset=0;

long int PIC::BC::InternalBoundary::Sphere::TotalSampleSetLength=0;
long int *PIC::BC::InternalBoundary::Sphere::SpeciesSampleDataOffset=NULL;
long int *PIC::BC::InternalBoundary::Sphere::SpeciesSampleUserDefinedDataOffset=NULL;
long int PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber=0;

bool PIC::BC::InternalBoundary::Sphere::UserDefinedSamplingProcedureFlag=false;
cAMRheap<cInternalSphericalData> PIC::BC::InternalBoundary::Sphere::InternalSpheres;

PIC::BC::InternalBoundary::Sphere::fPrintVariableList PIC::BC::InternalBoundary::Sphere::PrintUserDefinedVariableList=NULL;
PIC::BC::InternalBoundary::Sphere::fPrintTitle PIC::BC::InternalBoundary::Sphere::PrintUserDefinedTitle=NULL;
PIC::BC::InternalBoundary::Sphere::fPrintDataStateVector PIC::BC::InternalBoundary::Sphere::PrintUserDefinedDataStateVector=NULL;
PIC::BC::InternalBoundary::Sphere::fParticleSphereInteraction PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction=PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction_SpecularReflection;
PIC::BC::InternalBoundary::Sphere::fSampleParticleData PIC::BC::InternalBoundary::Sphere::SampleParticleData=NULL;


//====================================================
//init the data for the spherical body installed into the mesh:
void PIC::BC::InternalBoundary::Sphere::Init(long int* RequestedSamplingSetDataLength,long int* UserDefinedSampleDataRelativeOffset) {

  if (TotalSampleSetLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of PIC::BC::InternalBoundary::Sphere:");

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

void PIC::BC::InternalBoundary::Sphere::Init() {
  Init(NULL,NULL);
}

//====================================================
//register a new spherical body on the mesh
cInternalBoundaryConditionsDescriptor PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere() {
  cInternalSphericalData *newSphere;
  cInternalBoundaryConditionsDescriptor SphereDescriptor;
  cSurfaceDataSphere* SurfaceData;



  if (TotalSampleSetLength==0) exit(__LINE__,__FILE__,"Error: need to initialize PIC::BC::InternalBoundary::Sphere first");

  //set the new sphere
  newSphere=InternalSpheres.newElement();
  SurfaceData=(cSurfaceDataSphere*)malloc(sizeof(cSurfaceDataSphere));
  newSphere->SetSurfaceDataPointer((void*)SurfaceData);

  //set up the output functions
  newSphere->PrintVariableList=PrintDefaultVariableList;
  newSphere->PrintTitle=PrintDefaultTitle;
  newSphere->PrintDataStateVector=PrintDefaultDataStateVector;

  //init the internal parameters of the new sphere
  double *sBuffer;
  long int i,sBufferTotalLength;

  TotalSurfaceElementNumber=newSphere->GetTotalSurfaceElementsNumber();
  sBufferTotalLength=2*TotalSurfaceElementNumber*TotalSampleSetLength;


  SurfaceData->faceat=0;
  SurfaceData->SamplingBuffer=new double [sBufferTotalLength];

  for (i=0,sBuffer=SurfaceData->SamplingBuffer;i<sBufferTotalLength;i++) sBuffer[i]=0.0;

  //set the descriptor
  SphereDescriptor.BoundaryElement=(void*)newSphere;
  SphereDescriptor.BondaryType=_INTERNAL_BOUNDARY_TYPE_SPHERE_;

  PIC::Mesh::mesh.RegisterInternalBoundary(SphereDescriptor);

  return SphereDescriptor;
}


//====================================================
//clear sampling buffer
void PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(cInternalSphericalData* Sphere) {
  long int offset,i,nSurfaceElement,nSurfaceElementsMax;
  double *SamplingBuffer;

  SamplingBuffer=((cSurfaceDataSphere*)Sphere->GetSurfaceDataPointer())->SamplingBuffer;
  nSurfaceElementsMax=Sphere->GetTotalSurfaceElementsNumber();

  for (nSurfaceElement=0;nSurfaceElement<nSurfaceElementsMax;nSurfaceElement++) {
    offset=collectingSpecieSamplingDataOffset(0,nSurfaceElement);

    for (i=0;i<TotalSampleSetLength;i++) *(SamplingBuffer+offset+i)=0.0;
  }
}

//====================================================
//get access to the internal data of the sphere
PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere* PIC::BC::InternalBoundary::Sphere::GetSphereSurfaceData(cInternalBoundaryConditionsDescriptor Descriptor) {
  if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: wrong boundary type");

  return (cSurfaceDataSphere*)(((cInternalSphericalData*)Descriptor.BoundaryElement)->GetSurfaceDataPointer());

}

double* PIC::BC::InternalBoundary::Sphere::GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
  return GetSphereSurfaceData(Descriptor)->SamplingBuffer+completedCellSampleDataPointerOffset;
}
double* PIC::BC::InternalBoundary::Sphere::GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
  return GetSphereSurfaceData(Descriptor)->SamplingBuffer+collectingCellSampleDataPointerOffset;
}

//====================================================
//the offset of sampling data for a particular specie
int PIC::BC::InternalBoundary::Sphere::completeSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  return 2*SurfaceElement*TotalSampleSetLength+completedCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  return 2*SurfaceElement*TotalSampleSetLength+collectingCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Sphere::SurfaceElementSamplingSetLength() {
  return TotalSampleSetLength;
}

void PIC::BC::InternalBoundary::Sphere::switchSamplingBuffers() {
  long int tempSamplingOffset;

  tempSamplingOffset=completedCellSampleDataPointerOffset;
  completedCellSampleDataPointerOffset=collectingCellSampleDataPointerOffset;
  collectingCellSampleDataPointerOffset=tempSamplingOffset;
}

//====================================================
//particle-spherical surface interaction
//specular reflection of the particle velocity vector on the surface of the sphere
void PIC::BC::InternalBoundary::Sphere::ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor)  {
   double radiusSphere,*x0Sphere,l[3],r,vNorm,c;
   cInternalSphericalData *Sphere;
   cSurfaceDataSphere *SurfaceData;
   int idim;

   #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
   if (sphereDescriptor->BondaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: wrong boundary type");
   #endif

   #if DIM != 3
   exit(__LINE__,__FILE__,"Error: wrong dimension");
   #endif

   Sphere=(cInternalSphericalData*)sphereDescriptor->BoundaryElement;
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
   SurfaceData=(cSurfaceDataSphere*)Sphere->GetSurfaceDataPointer();
   SampleData=SurfaceData->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


///###########  DEBUG #############

//   SampleData=SurfaceData->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,0);

/*
   if (x[2]>0.20) {
     cout << __LINE__ << endl;
     Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
     nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

   }
*/

//########### END DEBUG #######





   SampleData[sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);

}


//====================================================
//print default surface data
void PIC::BC::InternalBoundary::Sphere::PrintDefaultVariableList(FILE* fout) {
  fprintf(fout,", \"Flux Down\"");
}
void PIC::BC::InternalBoundary::Sphere::PrintDefaultTitle(FILE* fout) {
  fprintf(fout,"SurfaceData");
}

void PIC::BC::InternalBoundary::Sphere::PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement;
  double InterpolationNormalization=0.0,InterpolationCoefficient;
  cSurfaceDataSphere *SurfaceData;
  double *SampleData;

  double FluxDown=0.0;

  SurfaceData=(cSurfaceDataSphere*)Sphere->GetSurfaceDataPointer();


  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];



//##### DEBUG #####

    /*
if (nSurfaceElement==311) {
  cout << __LINE__ << endl;
}

*/
//##### END DEBUG ######



    SampleData=SurfaceData->SamplingBuffer+completeSpecieSamplingDataOffset(spec,nSurfaceElement);
    InterpolationCoefficient=Sphere->GetSurfaceElementArea(nSurfaceElement);

    FluxDown+=SampleData[sampledFluxDownRelativeOffset]*InterpolationCoefficient;

    InterpolationNormalization+=InterpolationCoefficient;
  }


  FluxDown/=InterpolationNormalization;



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

void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);

*/
