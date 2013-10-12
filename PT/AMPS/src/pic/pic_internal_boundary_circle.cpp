//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the general functions for the circle internal boundary installed into the mesh

#include "pic.h"




int PIC::BC::InternalBoundary::Circle::completedCellSampleDataPointerOffset=0,PIC::BC::InternalBoundary::Circle::collectingCellSampleDataPointerOffset=0;

int PIC::BC::InternalBoundary::Circle::sampledFluxDownRelativeOffset=0,PIC::BC::InternalBoundary::Circle::sampledFluxUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Circle::sampledMeanVelocityDownRelativeOffset=0,PIC::BC::InternalBoundary::Circle::sampledMeanVelocityUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Circle::sampledMeanEnergyDownRelativeOffset=0,PIC::BC::InternalBoundary::Circle::sampledMeanEnergyUpRelativeOffset=0;
int PIC::BC::InternalBoundary::Circle::sampledSurfaceNumberDensityRelativeOffset=0;

long int PIC::BC::InternalBoundary::Circle::TotalSampleSetLength=0;
long int *PIC::BC::InternalBoundary::Circle::SpeciesSampleDataOffset=NULL;
long int *PIC::BC::InternalBoundary::Circle::SpeciesSampleUserDefinedDataOffset=NULL;
long int PIC::BC::InternalBoundary::Circle::TotalSurfaceElementNumber=0;

bool PIC::BC::InternalBoundary::Circle::UserDefinedSamplingProcedureFlag=false;
cAMRheap<cInternalCircleData> PIC::BC::InternalBoundary::Circle::InternalCircles;

PIC::BC::InternalBoundary::Circle::fPrintVariableList PIC::BC::InternalBoundary::Circle::PrintUserDefinedVariableList=NULL;
PIC::BC::InternalBoundary::Circle::fPrintTitle PIC::BC::InternalBoundary::Circle::PrintUserDefinedTitle=NULL;
PIC::BC::InternalBoundary::Circle::fPrintDataStateVector PIC::BC::InternalBoundary::Circle::PrintUserDefinedDataStateVector=NULL;
PIC::BC::InternalBoundary::Circle::fSampleParticleData PIC::BC::InternalBoundary::Circle::SampleParticleData=NULL;


//====================================================
//init the data for the circular body installed into the mesh:
void PIC::BC::InternalBoundary::Circle::Init(long int* RequestedSamplingSetDataLength,long int* UserDefinedSampleDataRelativeOffset) {

  if (TotalSampleSetLength!=0) exit(__LINE__,__FILE__,"Error: reinitialization of PIC::BC::InternalBoundary::Circle:");

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

void PIC::BC::InternalBoundary::Circle::Init() {
  Init(NULL,NULL);
}

//====================================================
//register a new circular body on the mesh
cInternalBoundaryConditionsDescriptor PIC::BC::InternalBoundary::Circle::RegisterInternalCircle() {
  cInternalCircleData *newCircle;
  cInternalBoundaryConditionsDescriptor CircleDescriptor;
//  cSurfaceDataCircle* SurfaceData;



  if (TotalSampleSetLength==0) exit(__LINE__,__FILE__,"Error: need to initialize PIC::BC::InternalBoundary::Circle first");

  //set the new Circle
  newCircle=InternalCircles.newElement();
//  SurfaceData=(cSurfaceDataCircle*)malloc(sizeof(cSurfaceDataCircle));
//  newCircle->SetSurfaceDataPointer((void*)SurfaceData);

  //set up the output functions
  newCircle->PrintVariableList=PrintDefaultVariableList;
  newCircle->PrintTitle=PrintDefaultTitle;
  newCircle->PrintDataStateVector=PrintDefaultDataStateVector;

  //init the internal parameters of the new Circle
  double *sBuffer;
  long int i,sBufferTotalLength;

  TotalSurfaceElementNumber=newCircle->GetTotalSurfaceElementsNumber();
  sBufferTotalLength=2*TotalSurfaceElementNumber*TotalSampleSetLength;


  newCircle->faceat=0;
  newCircle->SamplingBuffer=new double [sBufferTotalLength];
  newCircle->maxIntersectedNodeTimeStep=new double [PIC::nTotalSpecies];

  for (i=0,sBuffer=newCircle->SamplingBuffer;i<sBufferTotalLength;i++) sBuffer[i]=0.0;
  for (i=0;i<PIC::nTotalSpecies;i++) newCircle->maxIntersectedNodeTimeStep[i]=-1.0;

  //set the descriptor
  CircleDescriptor.BoundaryElement=(void*)newCircle;
  CircleDescriptor.BondaryType=_INTERNAL_BOUNDARY_TYPE_CIRCLE_;

  PIC::Mesh::mesh.RegisterInternalBoundary(CircleDescriptor);

  return CircleDescriptor;
}


//====================================================
//clear sampling buffer
void PIC::BC::InternalBoundary::Circle::flushCollectingSamplingBuffer(cInternalCircleData* Circle) {
  long int offset,i,nSurfaceElement,nSurfaceElementsMax;
  double *SamplingBuffer;

//  SamplingBuffer=((cSurfaceDataCircle*)Circle->GetSurfaceDataPointer())->SamplingBuffer;
  SamplingBuffer=Circle->SamplingBuffer;
  nSurfaceElementsMax=Circle->GetTotalSurfaceElementsNumber();

  for (nSurfaceElement=0;nSurfaceElement<nSurfaceElementsMax;nSurfaceElement++) {
    offset=collectingSpecieSamplingDataOffset(0,nSurfaceElement);

    for (i=0;i<TotalSampleSetLength;i++) *(SamplingBuffer+offset+i)=0.0;
  }
}

//====================================================
//get access to the internal data of the Circle
/*
PIC::BC::InternalBoundary::Circle::cSurfaceDataCircle* PIC::BC::InternalBoundary::Circle::GetCircleSurfaceData(cInternalBoundaryConditionsDescriptor Descriptor) {
  if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_Circle_) exit(__LINE__,__FILE__,"Error: wrong boundary type");

  return (cSurfaceDataCircle*)(((cInternalCircleData*)Descriptor.BoundaryElement)->GetSurfaceDataPointer());

}
*/

double* PIC::BC::InternalBoundary::Circle::GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
//  return GetCircleSurfaceData(Descriptor)->SamplingBuffer+completedCellSampleDataPointerOffset;


  return ((cInternalCircleData*)Descriptor.BoundaryElement)->SamplingBuffer+completedCellSampleDataPointerOffset;
}


double* PIC::BC::InternalBoundary::Circle::GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {
//  return GetCircleSurfaceData(Descriptor)->SamplingBuffer+collectingCellSampleDataPointerOffset;

  return ((cInternalCircleData*)Descriptor.BoundaryElement)->SamplingBuffer+collectingCellSampleDataPointerOffset;
}

//====================================================
//the offset of sampling data for a particular specie
int PIC::BC::InternalBoundary::Circle::completeSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  return 2*SurfaceElement*TotalSampleSetLength+completedCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Circle::collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
  return 2*SurfaceElement*TotalSampleSetLength+collectingCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
}

int PIC::BC::InternalBoundary::Circle::SurfaceElementSamplingSetLength() {
  return TotalSampleSetLength;
}

void PIC::BC::InternalBoundary::Circle::switchSamplingBuffers() {
  long int tempSamplingOffset;

  tempSamplingOffset=completedCellSampleDataPointerOffset;
  completedCellSampleDataPointerOffset=collectingCellSampleDataPointerOffset;
  collectingCellSampleDataPointerOffset=tempSamplingOffset;
}

//====================================================
//particle-circle surface interaction
//specular reflection of the particle velocity vector on the surface of the Circle
int PIC::BC::InternalBoundary::Circle::ParticleCircleInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPointer,void* CircleDataPointer)  {
   double radiusCircle,*x0Circle,l[3],r,vNorm,c;
   cInternalCircleData *Circle;
   cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*startNode;
//   cSurfaceDataCircle *SurfaceData;
   int idim;


   Circle=(cInternalCircleData*)CircleDataPointer;
   startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPointer;


   #if DIM != 3
   exit(__LINE__,__FILE__,"Error: wrong dimension");
   #endif


   Circle->GetCircleGeometricalParameters(x0Circle,radiusCircle);

   for (r=0.0,idim=0;idim<DIM;idim++) {
     l[idim]=x[idim]-x0Circle[idim];
     r+=pow(l[idim],2);
   }

   for (r=sqrt(r),vNorm=0.0,idim=0;idim<DIM;idim++) vNorm+=v[idim]*l[idim]/r;
   for (c=2.0*vNorm/r,idim=0;idim<DIM;idim++) v[idim]-=c*l[idim];

   //sample the particle data
   double *SampleData;
   long int nSurfaceElement,nPolarElement;

   Circle->GetSurfaceElementProjectionIndex(x,nPolarElement);
   nSurfaceElement=Circle->GetLocalSurfaceElementNumber(nPolarElement);
//   SurfaceData=Circle->SamplingBuffer;
   SampleData=Circle->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


///###########  DEBUG #############

//   SampleData=SurfaceData->SamplingBuffer+collectingSpecieSamplingDataOffset(spec,0);

/*
   if (x[2]>0.20) {
     cout << __LINE__ << endl;
     Circle->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
     nSurfaceElement=Circle->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

   }
*/

//########### END DEBUG #######





   SampleData[sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Circle->GetSurfaceElementArea(nPolarElement);
   return _PARTICLE_REJECTED_ON_THE_FACE_;
}


//====================================================
//print default surface data
void PIC::BC::InternalBoundary::Circle::PrintDefaultVariableList(FILE* fout) {
  fprintf(fout,", \"Flux Down\"");
}
void PIC::BC::InternalBoundary::Circle::PrintDefaultTitle(FILE* fout) {
  fprintf(fout,"TITLE=\"SurfaceData\"");
}

void PIC::BC::InternalBoundary::Circle::PrintDefaultDataStateVector(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement;
  double InterpolationNormalization=0.0,InterpolationCoefficient;
//  cSurfaceDataCircle *SurfaceData;
  double *SampleData;

  double FluxDown=0.0;

//  SurfaceData=(cSurfaceDataCircle*)Circle->GetSurfaceDataPointer();


  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];



//##### DEBUG #####

    /*
if (nSurfaceElement==311) {
  cout << __LINE__ << endl;
}

*/
//##### END DEBUG ######



    SampleData=Circle->SamplingBuffer+completeSpecieSamplingDataOffset(spec,nSurfaceElement);
    InterpolationCoefficient=Circle->GetSurfaceElementArea(nSurfaceElement);

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

void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* CircleDescriptor);

*/
