//$Id$

#include "OH.h"
#include "pic.h"

// user defined global time step
double OH::UserGlobalTimeStep = -1.0;

//  injection boundary condition
double OH::InjectionVelocity[3] = {26.3E3, 0.0, -2.3E3};
double OH::InjectionNDensity    = 0.18E6;
double OH::InjectionTemperature = 6519;

// computational domain size
double OH::DomainXMin[3] = {-2.25E14,-2.25E14,-2.25E14};
double OH::DomainXMax[3] = { 2.25E14, 2.25E14, 2.25E14};
double OH::DomainDXMin   = 1.8E13;
double OH::DomainDXMax   = 1.8E13;

// Declaring origin offset variable
long int OH::OffsetOriginTag = -1;

// OUTPUT ---------------------------------------------------------------------
int OH::Output::TotalDataLength = 0; 
int OH::Output::ohSourceDensityOffset =-1; 
int OH::Output::ohSourceMomentumOffset=-1;
int OH::Output::ohSourceEnergyOffset  =-1;


void OH::Output::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"ohSourceDensity\",\"ohSourceMomentumX\",\"ohSourceMomentumY\",\"ohSourceMomentumZ\",\"ohSourceEnergy\"");

  //if DataSetNumber is 0, and  OH::Sampling::OriginLocation::nSampledOriginLocations!=-1 -> output density of particles produced in each source region
  if (DataSetNumber==0) for (int i=0;i<OH::Sampling::OriginLocation::nSampledOriginLocations;i++) fprintf(fout,", \"ENA Density (Source Region ID=%i)\"",i);
}

void OH::Output::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode){

  double S1=0.0, S2[3]={0.0,0.0,0.0}, S3=0.0;
  int i,idim;
  char *offset;

  for (i=0;i<nInterpolationCoeficients;i++) {
    
    offset = 
      InterpolationList[i]->GetAssociatedDataBufferPointer() + 
      PIC::Mesh::completedCellSampleDataPointerOffset;

    S1+=(*((double*)(offset+OH::Output::ohSourceDensityOffset)))*InterpolationCoeficients[i];
    
    for(idim=0 ; idim<3; idim++)
      S2[idim]+=(*(idim+(double*)(offset+OH::Output::ohSourceMomentumOffset)))*InterpolationCoeficients[i];

    S3+=(*((double*)(offset+OH::Output::ohSourceEnergyOffset)))*InterpolationCoeficients[i];
  }
  
  offset = 
    CenterNode->GetAssociatedDataBufferPointer() + 
    PIC::Mesh::completedCellSampleDataPointerOffset;

  memcpy(offset+OH::Output::ohSourceDensityOffset, &S1,  sizeof(double));
  memcpy(offset+OH::Output::ohSourceMomentumOffset,&S2,3*sizeof(double));
  memcpy(offset+OH::Output::ohSourceEnergyOffset,  &S3,  sizeof(double));

  //evaluate density of ENAs produced in each source region
  for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations;iSource++) {
    double t=0.0;

    for (i=0;i<nInterpolationCoeficients;i++) {
      offset=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
      t+=(*(iSource+(double*)(offset+OH::Sampling::OriginLocation::OffsetDensitySample)))*InterpolationCoeficients[i];
    }

    offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
    *(iSource+(double*)(offset+OH::Sampling::OriginLocation::OffsetDensitySample))=t;
  }

}

void OH::Output::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode){
  double t;

  //SourceDensity
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceDensityOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //SourceMomentum
  for(int idim=0; idim < 3; idim++){
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(idim+(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceMomentumOffset));
    }
    
    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);
      
      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }

  //SourceEnergy
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Output::ohSourceEnergyOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //density of the ENAs produced in specific origin regions
  if (DataSetNumber==0) for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations;iSource++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(iSource+(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+OH::Sampling::OriginLocation::OffsetDensitySample));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t/PIC::LastSampleLength);
    }
    else pipe->send(t);
  }

}

int OH::Output::RequestDataBuffer(int offset){
  OH::Output::ohSourceDensityOffset=offset;
  OH::Output::TotalDataLength = 1;
  offset+=sizeof(double);

  OH::Output::ohSourceMomentumOffset=offset;
  OH::Output::TotalDataLength+=3;
  offset+=3*sizeof(double);

  OH::Output::ohSourceEnergyOffset=offset;
  OH::Output::TotalDataLength++;
  offset+=sizeof(double);

  return OH::Output::TotalDataLength*sizeof(double);
}


void OH::Output::Init() {
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestSamplingData.push_back(OH::Output::RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(OH::Output::PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(OH::Output::PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(OH::Output::Interpolate);
}


// Loss -----------------------------------------------------------------------

double OH::Loss::LifeTime(double *x, int spec, long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node){

  double PlasmaNumberDensity, PlasmaPressure, PlasmaTemperature;
  double PlasmaBulkVelocity[3];

  double lifetime=0.0;

  //in case of running in a stand-along mode
  #ifdef _OH_STAND_ALONG_MODE_
  #if _OH_STAND_ALONG_MODE_ == _PIC_MODE_ON_
  PhotolyticReactionAllowedFlag=false;
  return lifetime;
  #endif
  #endif



  PIC::CPLR::InitInterpolationStencil(x,node);

  PlasmaNumberDensity = PIC::CPLR::GetBackgroundPlasmaNumberDensity();
  PlasmaPressure      = PIC::CPLR::GetBackgroundPlasmaPressure();
  PlasmaTemperature   = PlasmaPressure / (2*Kbol * PlasmaNumberDensity);
  PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);

  // change this to false to turn off charge-exchange
  PhotolyticReactionAllowedFlag=true;

  // this model has hydrogen and three charge-exchange species only
  if (spec!=_H_SPEC_ && spec != _H_ENA_V1_SPEC_ && spec != _H_ENA_V2_SPEC_ && spec != _H_ENA_V3_SPEC_) {
    PhotolyticReactionAllowedFlag=false;
    return -1.0;
  }

  // velocity of a particle
  double v[3];
  PIC::ParticleBuffer::GetV(v,ptr);  
  //-------------------
  switch (spec) {
  case _H_SPEC_:
    lifetime= ChargeExchange::LifeTime(_H_SPEC_, v, PlasmaBulkVelocity, PlasmaTemperature, PlasmaNumberDensity);
    break;
  case _H_ENA_V1_SPEC_:
    lifetime= ChargeExchange::LifeTime(_H_ENA_V1_SPEC_, v, PlasmaBulkVelocity, PlasmaTemperature, PlasmaNumberDensity);
    break;
  case _H_ENA_V2_SPEC_:
    lifetime= ChargeExchange::LifeTime(_H_ENA_V2_SPEC_, v, PlasmaBulkVelocity, PlasmaTemperature, PlasmaNumberDensity);
    break;
  case _H_ENA_V3_SPEC_:
    lifetime= ChargeExchange::LifeTime(_H_ENA_V3_SPEC_, v, PlasmaBulkVelocity, PlasmaTemperature, PlasmaNumberDensity);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown species");
  }

  return lifetime;

}

void OH::Loss::ReactionProcessor(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node){
  //as a result of the reaction only velocity of a particle is changed
  //----------------------------------------------------------------------

  //the life time of the original particle
  int spec;
  PIC::ParticleBuffer::byte *ParticleData;
  double xParent[3],vParent[3],ParentLifeTime,ParentTimeStep;
  bool ReactionOccurredFlag;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  PIC::ParticleBuffer::GetX(xParent,ParticleData);
  PIC::ParticleBuffer::GetV(vParent,ParticleData);

  ParentLifeTime=LifeTime(xParent,spec,ptr,ReactionOccurredFlag,node);

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
  ParentTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  if (ReactionOccurredFlag==true) {
    //reaction is theoretically possible;
    //determine whether the reaction has occurred by comparing the local time step with the partent particle life time

    if (rnd()<exp(-ParentTimeStep/ParentLifeTime)) {
      //the particle is remain in the system
      ReactionOccurredFlag=false;
    }
  }

  if (ReactionOccurredFlag==true) {
    //inject the products of the reaction
    double ParentTimeStep,ParentParticleWeight;

  #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
  #else
    ParentParticleWeight=0.0;
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
  #endif

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
  #else
    ParentTimeStep=0.0;
    exit(__LINE__,__FILE__,"Error: the time step node is not defined");
  #endif

    //account for the parent particle correction factor
    ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

    // new particle comes from solar wind and has random velocity from proton Maxwellian distribution
    double PlasmaBulkVelocity[3];
    double PlasmaNumberDensity, PlasmaPressure, PlasmaTemperature;
    double vp[3];
    
    PIC::CPLR::InitInterpolationStencil(xParent,node);
    PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity);
    PlasmaNumberDensity = PIC::CPLR::GetBackgroundPlasmaNumberDensity();
    PlasmaPressure      = PIC::CPLR::GetBackgroundPlasmaPressure();
    PlasmaTemperature   = PlasmaPressure / (2*Kbol * PlasmaNumberDensity);
    
    // calculating the random velocity of the proton from the maxwellian velocity of the local plasma
    PIC::Distribution::MaxwellianVelocityDistribution(vp,PlasmaBulkVelocity,PlasmaTemperature,spec);

    // charge exchange process transfers momentum and energy to plasma
    PIC::Mesh::cDataCenterNode *CenterNode;
    char *offset;

    CenterNode=PIC::Mesh::Search::FindCell(xParent); ///node->block->GetCenterNode(nd);
    offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
    double vh2 = 0.0, vp2 = 0.0;
    double c = ParentParticleWeight
      / PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]
      / CenterNode->Measure;
    *((double*)(offset+OH::Output::ohSourceDensityOffset)) += 0.0;

    for(int idim=0; idim<3; idim++){
      *(idim + (double*)(offset+OH::Output::ohSourceMomentumOffset)) +=
	c*_MASS_(_H_)*(vParent[idim]-vp[idim]);
      vh2      +=vParent[idim]*vParent[idim];
      vp2+=vp[idim]*vp[idim];
    }
    *((double*)(offset+OH::Output::ohSourceEnergyOffset)) +=
      c*0.5*_MASS_(_H_)*(vh2-vp2);
   
    // creating new neutral particle with the velocity of the selected proton
    PIC::ParticleBuffer::SetV(vp,ParticleData);
    // adding neutral to correct species depending on its velocity
    if (_H_ENA_V3_SPEC_ >= 0 && sqrt(vp2) >= 500.0E3) {
      PIC::ParticleBuffer::SetI(_H_ENA_V3_SPEC_,ParticleData);
    }
    else if (_H_ENA_V2_SPEC_ >=0 && sqrt(vp2)>=150.0E3 && sqrt(vp2)<500.0E3) {
      PIC::ParticleBuffer::SetI(_H_ENA_V2_SPEC_,ParticleData);
    }
    else if (_H_ENA_V1_SPEC_ >=0 && sqrt(vp2)>=50.0E3 && sqrt(vp2)<150.0E3) {
      PIC::ParticleBuffer::SetI(_H_ENA_V1_SPEC_,ParticleData);
    }
    else {
      PIC::ParticleBuffer::SetI(_H_SPEC_,ParticleData);
    }
    // tagging the particle to the right population that it was created
    OH::SetOriginTag(OH::GetEnaOrigin(PlasmaNumberDensity,PlasmaPressure,PlasmaBulkVelocity), ParticleData);
  }

  //add the particle to the list of the particles existing in the system after reaction
  PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
  PIC::ParticleBuffer::SetPrev(-1,ptr);

  if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
  FirstParticleCell=ptr;
}


void OH::Init_BeforeParser(){
  Exosphere::Init_BeforeParser();
  PIC::ParticleBuffer::RequestDataStorage(OH::OffsetOriginTag, sizeof(int));
  OH::Output::Init();

  //set the coupling procedure
  PIC::CPLR::SWMF::SendCenterPointData.push_back(Coupling::Send);

  //request sampling data
  PIC::IndividualModelSampling::RequestSamplingData.push_back(Sampling::OriginLocation::RequestSamplingData);
}

// User defined functions -----------------------------------------------------
int OH::user_set_face_boundary(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  //setting User defined function to process particles leaving domain at certain faces
  //useful for 1D runs if just want flows in one direction
  
  // removing particles if hit faces perpandiulat to x axis
  if (nIntersectionFace == 0 || nIntersectionFace == 1) return _PARTICLE_DELETED_ON_THE_FACE_; 

  // keep and reflect particles if hit face perpandiculat to y or z axes
  // y axis reflection
  if (nIntersectionFace == 2 || nIntersectionFace == 3) {
    vInit[1]=-vInit[1];
    return _PARTICLE_REJECTED_ON_THE_FACE_; // particles are not deleted but remain in domain
  }

  // z axis reflection
  if (nIntersectionFace == 4 || nIntersectionFace == 5) {
    vInit[2]=-vInit[2];
    return _PARTICLE_REJECTED_ON_THE_FACE_;
  }
}

//-----------------------------------------------------------------------------
//substitutes for Exosphere functions
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_],Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_],Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_];

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=1;
  }

  return nVariables;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;


  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    //    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  //do nothing
}

double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {
  ReemissionParticleFraction=0.0;

  return 1.0;
}


double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {


  return 100.0;
}

int OH::GetEnaOrigin(double PlasmaNumberDensity, double PlasmaPressure, double *PlasmaBulkVelocity){
  // determines which population the ENA should be added to based on where in the heliosphere it was created
  // 0 = SSSW
  // 1 = inner heliosheath
  // 2 = disturbed ISM / Outer heliosheath
  // 3 = Prisine ISM

  double mach=0.0;
  double PlasmaBulkSpeed=0.0,PlasmaBulkSpeed2=0.0;
  double PlasmaTemperature=0.0;

  PlasmaTemperature   = PlasmaPressure / (2*Kbol * PlasmaNumberDensity);

  for(int idim=0; idim<3; idim++) PlasmaBulkSpeed2 += PlasmaBulkVelocity[idim]*PlasmaBulkVelocity[idim];
  PlasmaBulkSpeed = sqrt(PlasmaBulkSpeed2);

  mach=PlasmaBulkSpeed*sqrt(3.0*PlasmaNumberDensity*_MASS_(_H_)/(5.0*PlasmaPressure));

  // adding neutral to correct population where it underwent charge exchange
  if (mach >= 1.0 && PlasmaBulkSpeed < 100.0E3) {
    // particle is in the pristine ISM - population 4
    return 3;
  }
  else if (mach < 1.0 && PlasmaTemperature < 3.28E5 && PlasmaBulkSpeed < 100.0E3) {
    // particle is in the disturbed ISM / Outer heliosheath - population 3
    return 2;
  }
  else if (mach < 1.0 && PlasmaTemperature >= 3.28E5) {
    // particle is in the HS / inner heliosheath - population 2
    return 1;
  }
  else {
    // particle is in the super sonic SW - population 1
    return 0;
  }
}

//=====================================================================================================
//sampling of the ENAs density individually for each origin region
int OH::Sampling::OriginLocation::nSampledOriginLocations=-1;
int OH::Sampling::OriginLocation::OffsetDensitySample=-1;

int OH::Sampling::OriginLocation::RequestSamplingData(int offset) {
  int res=0;

  if (nSampledOriginLocations!=-1) {
    OffsetDensitySample=offset;
    res=nSampledOriginLocations*sizeof(double);
  }

  return res;
}

void OH::Sampling::OriginLocation::SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec) {
  int OriginID;

  if (nSampledOriginLocations!=-1) {
    OriginID=OH::GetOriginTag((PIC::ParticleBuffer::byte*)ParticleData);

    if (OriginID>=nSampledOriginLocations) exit(__LINE__,__FILE__,"Error: OriginID is out of range. Update the input file with the corrected value of the number of the source regions");
    *(OriginID+(double*)(SamplingBuffer+OffsetDensitySample))+=LocalParticleWeight;
  }
}
