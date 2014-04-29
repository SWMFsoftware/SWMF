//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * ElectricallyChargedDust.cpp
 *
 *  Created on: Jan 20, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_



#include "pic__model__electrically_charged_dust.h"


//the offsets to the internal properties of the dust grains
double ElectricallyChargedDust::SampledDustMassInjectionRate=0.0;
double ElectricallyChargedDust::TotalMassDustProductionRate=0.0;
ElectricallyChargedDust::fGenerateNewDustGrainInternalProperties ElectricallyChargedDust::GenerateNewDustGrainInternalProperties=NULL;
ElectricallyChargedDust::fGrainDragCoeffiient ElectricallyChargedDust::GrainDragCoeffiient_UserDefinedFunction=NULL;
ElectricallyChargedDust::fGrainMassDensity ElectricallyChargedDust::GrainMassDensity_UserDefinedFunction=NULL;


//int ElectricallyChargedDust::_DUST_SPEC_=-1;
double ElectricallyChargedDust::minDustRadius=-1.0,ElectricallyChargedDust::maxDustRadius=-1.0;

//velocity groups of the grains

double ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity=0.0,ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity=0.0;
double ElectricallyChargedDust::GrainVelocityGroup::logMinGrainVelocity=0.0,ElectricallyChargedDust::GrainVelocityGroup::logMaxGrainVelocity=0.0;
double ElectricallyChargedDust::GrainVelocityGroup::dLogGrainVelocity=0.0;
int ElectricallyChargedDust::GrainVelocityGroup::nGroups=0;


//size distribution of the dust grains
double ElectricallyChargedDust::SizeDistribution::PowerIndex=0.0,ElectricallyChargedDust::SizeDistribution::NormalizationFactor=0.0;
double ElectricallyChargedDust::SizeDistribution::LogMinDustGrainRadius=0.0,ElectricallyChargedDust::SizeDistribution::LogMaxDustGrainRadius=0.0;
double ElectricallyChargedDust::SizeDistribution::maxValue__non_normalized=0.0;

int ElectricallyChargedDust::Sampling::nDustSizeSamplingIntervals=0;
double ElectricallyChargedDust::Sampling::dLogDustSamplingIntervals=0.0;
int ElectricallyChargedDust::Sampling::CellSamplingDataOffset=-1;
int ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset=-1,ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Velocity_Offset=-1;
int ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Speed_Offset=-1;
int ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCherge_Offset=-1,ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCurrent_Offset=-1;
int ElectricallyChargedDust::Sampling::TotalDustElectricChargeDensitySamplingOffset=-1,ElectricallyChargedDust::Sampling::TotalDustElectricCurrentDensitySamplingOffset=-1;

//sampling size distribution function
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::nSamplingIntervals=0;
double ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::dLog10DustRadius=0.0;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::nSamplingLocations=0;
double **ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SamplingLocations=NULL;
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SampleNodes=NULL;
double **ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SamplingBuffer=NULL;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SamplingIntervalDataLength=0;
long int *ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SampleLocalCellNumber=NULL;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::NumberDensitySamplingOffset=0;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::VelocitySamplingOffset=0;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SpeedSamplingOffset=0;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::ElectricChargeSamplingOffset=0;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::ElectricCurentSamplingOffset=0;
int ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::AbsoluteElectricCurrentValueSamplingOffset=0;
#endif


//set up the sampling procedures of the model
void ElectricallyChargedDust::Sampling::SetDustSamplingIntervals(int nIntervals) {
  if ((minDustRadius<0.0)||(maxDustRadius<0.0)) exit(__LINE__,__FILE__,"Error: The function can be used only after the grain's size limit is defeiend");
  if (CellSamplingDataOffset!=-1) exit(__LINE__,__FILE__,"Error: the function must be used before the sampling buffers for the dust are allocated");

  nDustSizeSamplingIntervals=nIntervals;
  dLogDustSamplingIntervals=log(maxDustRadius/minDustRadius)/nIntervals;
}

int ElectricallyChargedDust::Sampling::RequestSamplingData(int offset) {
  int SamplingLength=0;

  if (CellSamplingDataOffset!=-1) exit(__LINE__,__FILE__,"Error: second request for the sampling data");

  if (nDustSizeSamplingIntervals>0) {
    CellSamplingDataOffset=offset;

    DustSizeSamplingInterval_NumberDensity_Offset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=sizeof(double)*nDustSizeSamplingIntervals;

    DustSizeSamplingInterval_Velocity_Offset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=3*sizeof(double)*nDustSizeSamplingIntervals;

    DustSizeSamplingInterval_Speed_Offset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=sizeof(double)*nDustSizeSamplingIntervals;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    DustSizeSamplingInterval_ElectricCherge_Offset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=sizeof(double)*nDustSizeSamplingIntervals;

    DustSizeSamplingInterval_ElectricCurrent_Offset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=3*sizeof(double)*nDustSizeSamplingIntervals;

    TotalDustElectricChargeDensitySamplingOffset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=sizeof(double);

    TotalDustElectricCurrentDensitySamplingOffset=CellSamplingDataOffset+SamplingLength;
    SamplingLength+=3*sizeof(double);
#endif
  }
  else CellSamplingDataOffset=0;

  return SamplingLength;
}

//init the data for the model
void ElectricallyChargedDust::Init_BeforeParser() {
  PIC::IndividualModelSampling::RequestSamplingData.push_back(Sampling::RequestSamplingData);
}



void ElectricallyChargedDust::Init_AfterParser() {

  //calcualte the number of the dust time-step groups
  bool DustSpeciesBegin=false,DustSpacesEnd=false;
  int spec;
  char ChemSymbol[_MAX_STRING_LENGTH_PIC_];

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);

    if ((DustSpeciesBegin==false)&&(DustSpacesEnd==false)&&(strcmp(ChemSymbol,"DUST")==0)) {
      //the first dust group
      DustSpeciesBegin=true;
      GrainVelocityGroup::nGroups=1; //,_DUST_SPEC_=spec;
    }
    else if ((DustSpeciesBegin==true)&&(DustSpacesEnd==false)&&(strcmp(ChemSymbol,"DUST")==0)) {
      //the list of the dust time-step groups continues
      GrainVelocityGroup::nGroups++;
    }
    else if ((DustSpeciesBegin==true)&&(DustSpacesEnd==false)&&(strcmp(ChemSymbol,"DUST")!=0)) {
      //the list of the dust groups stopped
      DustSpacesEnd=true;
    }
    else if ((DustSpeciesBegin==true)&&(DustSpacesEnd==true)&&(strcmp(ChemSymbol,"DUST")==0)) {
     exit(__LINE__,__FILE__,"Error: a discontinuous list of the DUST groups in the imput file is found");
    }
    else if (strcmp(ChemSymbol,"DUST")==0) exit(__LINE__,__FILE__,"Error: the conditions list does not account for the occered event");
  }

  //init dust sub-models
  SizeDistribution::Init();
  GrainVelocityGroup::Init();

  //request the particle data storage
  // NO ADDITIONAL DATA IS REQUESTED

  //init the sampling procedues
  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SampleModelData,OutputSampledModelData);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);

  //suppress output of the data files for the additional groups of the dust specie
  if (PIC::Sampling::SaveOutputDataFile==NULL) exit(__LINE__,__FILE__,"Error: PIC::Sampling::SaveOutputDataFile is not allocated");
  for (int spec=_DUST_SPEC_+1;spec<_DUST_SPEC_+GrainVelocityGroup::nGroups;spec++) PIC::Sampling::SaveOutputDataFile[spec]=false;
}


/*-------------------------------------------------  Sampling Procedires ----------------------------------------------*/

void ElectricallyChargedDust::OutputSampledModelData(int DataOutputFileNumber) {
  double buffer[PIC::nTotalThreads];

  MPI_Gather(&SampledDustMassInjectionRate,1,MPI_DOUBLE,buffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    for (int thread=1;thread<PIC::nTotalThreads;thread++) buffer[0]+=buffer[thread];

    fprintf(PIC::DiagnospticMessageStream,"The total dust production rate: %e [kg/sec]\n",buffer[0]/PIC::LastSampleLength);
  }

  if (PIC::SamplingMode==_RESTART_SAMPLING_MODE_) SampledDustMassInjectionRate=0.0;
}

void ElectricallyChargedDust::SampleModelData() {
//do nothing
}



/*------------------------------------------------  Evaluate the local times step --------------------------------------*/
//evaluate the local time step: the function returns 'false' when the evaluation of the time step is requested for a NON-DUST specie
bool ElectricallyChargedDust::EvaluateLocalTimeStep(int spec,double &dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CharacteristicCellSize,minGroupSpeed,maxGroupSpeed;
  int DustGroup=spec-_DUST_SPEC_;

  if ((DustGroup<0)||(DustGroup>=GrainVelocityGroup::nGroups)) return false;

  GrainVelocityGroup::GetGroupVelocityRange(minGroupSpeed,maxGroupSpeed,DustGroup);
  CharacteristicCellSize=startNode->GetCharacteristicCellSize();

  dt=0.1*CharacteristicCellSize/maxGroupSpeed;
  return true;
}

/*----------------------------------------------  Total Acceleration of Dust Grains ------------------------------------*/
/*
void ElectricallyChargedDust::TotalGrainAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  char ParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double GrainCharge,GrainMass,GrainRadius;

  memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

  GrainCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);
  GrainMass=GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);
  GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);


  char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
  PIC::Mesh::cDataCenterNode *CenterNode;
  double *E,*B;
  int i,j,k;

  //copy to local variables
  double accl_LOCAL[3]={0.0,0.0,0.0},x_LOCAL[3],v_LOCAL[3];

  memcpy(x_LOCAL,x,3*sizeof(double));
  memcpy(v_LOCAL,v,3*sizeof(double));


#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  long int nd;

  if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
    exit(__LINE__,__FILE__,"Error: the cell is not found");
  }

  if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif

  CenterNode=startNode->block->GetCenterNode(nd);
  memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+CenterNode->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);

  //Lorentz force
  E=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::ElectricFieldOffset);
  B=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::MagneticFieldOffset);

  accl_LOCAL[0]+=GrainCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/GrainMass;
  accl_LOCAL[1]+=GrainCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/GrainMass;
  accl_LOCAL[2]+=GrainCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/GrainMass;

  //the gravity force
  double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
  double c=GravityConstant*_MASS_(_TARGET_)/pow(r2,3.0/2.0);

  accl_LOCAL[0]-=c*x_LOCAL[0];
  accl_LOCAL[1]-=c*x_LOCAL[1];
  accl_LOCAL[2]-=c*x_LOCAL[2];

  //Drag force
  / *
  double A,cr2;
  double *BackgroundAtmosphereBulkVelocity=(double*)(ICES_AssociatedData+PIC::CPLR::ICES::NeutralBullVelocityOffset);
  double BackgroundAtmosphereNumberDensity=*((double*)(ICES_AssociatedData+PIC::CPLR::ICES::NeutralNumberDensityOffset));

  cr2=(v_LOCAL[0]-BackgroundAtmosphereBulkVelocity[0])*(v_LOCAL[0]-BackgroundAtmosphereBulkVelocity[0])+
    (v_LOCAL[1]-BackgroundAtmosphereBulkVelocity[1])*(v_LOCAL[1]-BackgroundAtmosphereBulkVelocity[1])+
    (v_LOCAL[2]-BackgroundAtmosphereBulkVelocity[2])*(v_LOCAL[2]-BackgroundAtmosphereBulkVelocity[2]);


  A=Pi*pow(GrainRadius,2)/2.0*GrainDragCoefficient*sqrt(cr2)/GrainMass*BackgroundAtmosphereNumberDensity*_MASS_(_H2O_);

  accl_LOCAL[0]+=A*(BackgroundAtmosphereBulkVelocity[0]-v_LOCAL[0]);
  accl_LOCAL[1]+=A*(BackgroundAtmosphereBulkVelocity[1]-v_LOCAL[1]);
  accl_LOCAL[2]+=A*(BackgroundAtmosphereBulkVelocity[2]-v_LOCAL[2]);
   * /

  //copy the accaleration vector from the internal buffer
  memcpy(accl,accl_LOCAL,3*sizeof(double));
}
*/
/*----------------------------------------------  Dust Charging ---------------------------------------*/

//charging of the dust grains
/*int ElectricallyChargedDust::DustChargingProcessorIndicator(double *x,double *v,int spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double &dt,bool &TransformationTimeStepLimitFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
}
*/

/*
int ElectricallyChargedDust::DustChargingProcessor(double *xInit,double *xFinal,double *v,int &spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *initNode) {
  double electronCurrent,ionCurrent,plasmaTemperature,plasmaNumberDensity;
  double GrainElectricCharge;

return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;

  PIC::Mesh::cDataCenterNode* cell;
  int i,j,k;
  long int LocalCellNumber;

  if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(xInit,i,j,k,initNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
  cell=initNode->block->GetCenterNode(LocalCellNumber);

  //ICES plasma data
  char ICES_AssociatedData[PIC::CPLR::ICES::TotalAssociatedDataLength];
  double *swVel;

  memcpy(ICES_AssociatedData,PIC::CPLR::ICES::AssociatedDataOffset+cell->GetAssociatedDataBufferPointer(),PIC::CPLR::ICES::TotalAssociatedDataLength);
  plasmaTemperature=*((double*)(PIC::CPLR::ICES::PlasmaTemperatureOffset+ICES_AssociatedData));
  swVel=(double*)(PIC::CPLR::ICES::PlasmaBulkVelocityOffset+ICES_AssociatedData);
  plasmaNumberDensity=*((double*)(PIC::CPLR::ICES::PlasmaNumberDensityOffset+ICES_AssociatedData));

  if (plasmaTemperature<=0.0) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;

  //get the grain electric potential
  char localParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double M,GrainRadius,dustPotential;

  memcpy((void*)localParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
  GrainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)localParticleData);
  GrainElectricCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)localParticleData);
  dustPotential=GrainElectricCharge/(4.0*Pi*VacuumPermittivity*GrainRadius);

  //get the electron current (Horanyi-1996-ARAA, Eq. 4)
  double xi=-ElectronCharge*dustPotential/(Kbol*plasmaTemperature);

  electronCurrent=-ElectronCharge*  4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ElectronMass));
  electronCurrent*=(xi>=0.0) ? exp(-xi) : 1.0-xi;


  //get the ion current (Horanyi-1996-ARAA, Eq. 4,14)

  xi=ElectronCharge*dustPotential/(Kbol*plasmaTemperature);
  M=sqrt((pow(swVel[0]-v[0],2)+pow(swVel[1]-v[1],2)+pow(swVel[2]-v[2],2))/(2.0*Kbol*plasmaTemperature/ProtonMass));

  ionCurrent=ElectronCharge* 4.0*Pi*pow(GrainRadius,2)*plasmaNumberDensity*sqrt(Kbol*plasmaTemperature/(PiTimes2*ProtonMass));

  if (dustPotential<0.0) {
    ionCurrent*=((pow(M,2)+0.5-xi)*sqrt(Pi)/M*erf(M)+exp(-pow(M,2)))/2.0;
  }
  else {
    double sqrt_xi=sqrt(xi);

    ionCurrent*=(
        (pow(M,2)+0.5-xi)*sqrt(Pi)/M*(erf(M+sqrt_xi)+erf(M-sqrt_xi))  +
        (sqrt(xi/M)+1.0)*exp(-pow(M-sqrt_xi,2.0)) -
        (sqrt(xi/M)-1.0)*exp(-pow(M+sqrt_xi,2.0))
        )/2.0;
  }

  //update the grain charge (Horanyi-1996-ARAA, Eq.1)
  GrainElectricCharge+=(ionCurrent+electronCurrent)*dt;

  SetGrainCharge(GrainElectricCharge,ParticleData);

  //move the particle into diferent velocity group if needed
  int oldVelocityGroup,newVelocityGroup;

  oldVelocityGroup=spec-_DUST_SPEC_;
  newVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);

  if (oldVelocityGroup!=newVelocityGroup) {
    //move the particle into different velocity group
    double GrainWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

    GrainWeightCorrection*=initNode->block->GetLocalTimeStep(_DUST_SPEC_+newVelocityGroup)/initNode->block->GetLocalTimeStep(_DUST_SPEC_+oldVelocityGroup);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+oldVelocityGroup]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+newVelocityGroup];
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,ParticleData);
    spec=_DUST_SPEC_+newVelocityGroup;
    PIC::ParticleBuffer::SetI(spec,ParticleData);
  }



  return _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_;
}
*/

/*-----------------------------------------   Inject dust grains from a sphere -----------------------------------*/
long int ElectricallyChargedDust::DustInjection__Sphere(void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ParticleWeight,LocalTimeStep,x[3],v[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  PIC::Mesh::cDataBlockAMR *block;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double GrainRadius,GrainMass,GrainWeightCorrection;
  bool InjectionFlag;
  ElectricallyChargedDust::fGenerateNewDustGrainInternalProperties GenerateNewDustGrainInternalProperties_LOCAL=GenerateNewDustGrainInternalProperties;
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  int GrainVelocityGroup;

  const int ParticleDataLength=PIC::ParticleBuffer::ParticleDataLength;
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  //============  DEBUG =========================
  static double GrainInjectedMass=0.0;


//  if (nDustGroups!=1) exit(__LINE__,__FILE__,"Error: the injection procedure doesn't work correctly for nDustGroups!=1");

  //============ END DEBUG =========================

  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups-1];
#else
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[_DUST_SPEC_];
#endif


  GrainInjectedMass+=TotalMassDustProductionRate*LocalTimeStep;


  while (GrainInjectedMass>0.0) {
    startNode=NULL;
    InjectionFlag=GenerateNewDustGrainInternalProperties_LOCAL(x,v,GrainRadius,GrainWeightCorrection,startNode);

    GrainMass=4.0/3.0*Pi*MeanDustDensity*pow(GrainRadius,3);
    GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;


    if (InjectionFlag==false) continue;
    if ((block=startNode->block)->GetLocalTimeStep(_DUST_SPEC_)/LocalTimeStep<rnd()) continue;

    //determine the velocity group of the injected grain;
    //calculate additional particle weight correction because the particle will be placed in a different weight group
    GrainVelocityGroup=GrainVelocityGroup::GetGroupNumber(v);
    GrainWeightCorrection*=block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/block->GetLocalTimeStep(_DUST_SPEC_);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


    //generate a particle:init the temporary buffer with the particle's properties
    PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::SetI(_DUST_SPEC_+GrainVelocityGroup,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);

    SetGrainCharge(0.0,(PIC::ParticleBuffer::byte*)tempParticleData);
    SetGrainMass(GrainMass,(PIC::ParticleBuffer::byte*)tempParticleData);
    SetGrainRadius(GrainRadius,(PIC::ParticleBuffer::byte*)tempParticleData);


    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

    memcpy((void*)newParticleData,(void*)tempParticleData,ParticleDataLength);
    nInjectedParticles++;



    //sample the injection flux
    //sample the particle data
    double *SampleData;
    long int nSurfaceElement,nZenithElement,nAzimuthalElement;

    Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
    nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
    SampleData=Sphere->SamplingBuffer+PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(_DUST_SPEC_,nSurfaceElement);


    SampleData[PIC::BC::InternalBoundary::Sphere::sampledFluxUpRelativeOffset]+=GrainMass*ParticleWeight*GrainWeightCorrection/block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);



    //inject the particle into the system
    SampledDustMassInjectionRate+=GrainMass*ParticleWeight*GrainWeightCorrection/startNode->block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup);
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,rnd()*startNode->block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup),startNode,true);
  }

  return nInjectedParticles;
}


/*----------------------------------------------- Print Output File -------------------------------------------*/
//functions that output the data file
void ElectricallyChargedDust::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Total Dust Number Density [m^{-3}]\", \"Mean Grains Bulk Velocity[0]\", \"Mean Grains Bulk Velocity[1]\", \"Mean Grains Bulk Velocity[2]\", \"Mean Grains Speed\"");

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  fprintf(fout,", \"Total Dust Electric Charge Desnity\", \"Dust Electric Current Density[0]\", \"Dust Electric Current Density[1]\", \"Dust Electric Current Density[2]\"");
#endif


  if (Sampling::nDustSizeSamplingIntervals!=0) {
    //print variable list for each dust size interval

    for (int nInterval=0;nInterval<Sampling::nDustSizeSamplingIntervals;nInterval++) {
      fprintf(fout,", \"Dust Number Density (%e<a<%e)\"",minDustRadius*exp(nInterval*Sampling::dLogDustSamplingIntervals),minDustRadius*exp((1+nInterval)*Sampling::dLogDustSamplingIntervals));

      for (int idim=0;idim<DIM;idim++) {
        fprintf(fout,", \"Dust Velocity[%i] (%e<a<%e)\"",idim,minDustRadius*exp(nInterval*Sampling::dLogDustSamplingIntervals),minDustRadius*exp((1+nInterval)*Sampling::dLogDustSamplingIntervals));
      }

      fprintf(fout,", \"Dust Grains Speed (%e<a<%e)\"",minDustRadius*exp(nInterval*Sampling::dLogDustSamplingIntervals),minDustRadius*exp((1+nInterval)*Sampling::dLogDustSamplingIntervals));

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
      fprintf(fout,", \"Dust Charge Density (%e<a<%e)\"",minDustRadius*exp(nInterval*Sampling::dLogDustSamplingIntervals),minDustRadius*exp((1+nInterval)*Sampling::dLogDustSamplingIntervals));

      for (int idim=0;idim<DIM;idim++) {
        fprintf(fout,", \"Dust Current Density[%i] (%e<a<%e)\"",idim,minDustRadius*exp(nInterval*Sampling::dLogDustSamplingIntervals),minDustRadius*exp((1+nInterval)*Sampling::dLogDustSamplingIntervals));
      }
#endif
    }
  }
}


void ElectricallyChargedDust::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int i,idim,nInterval;
  double c,c0,c1,Measure;
  char *SamplingBuffer,*CellNodeSamplingBuffer;

  CellNodeSamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

  //interpolate the total electric charged and current densities
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  double chargeTotal=0.0,currentTotal[3]={0.0,0.0,0.0},TotalMeasure=0.0;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

    Measure=InterpolationList[i]->Measure;
    if (Measure<=0.0) exit(__LINE__,__FILE__,"Error: non-positive cell volume is found");

    TotalMeasure+=Measure;

    chargeTotal+=*((double*)(SamplingBuffer+Sampling::TotalDustElectricChargeDensitySamplingOffset));
    for (idim=0;idim<3;idim++) currentTotal[idim]+=*(idim+(double*)(SamplingBuffer+Sampling::TotalDustElectricCurrentDensitySamplingOffset));
  }



  if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) chargeTotal/=PIC::LastSampleLength*TotalMeasure;
  *((double*)(CellNodeSamplingBuffer+Sampling::TotalDustElectricChargeDensitySamplingOffset))=chargeTotal;


  for (idim=0;idim<3;idim++) {
    if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) currentTotal[idim]/=PIC::LastSampleLength*TotalMeasure;
    *(idim+(double*)(CellNodeSamplingBuffer+Sampling::TotalDustElectricCurrentDensitySamplingOffset))=currentTotal[idim];
  }
#endif

  for (nInterval=0;nInterval<Sampling::nDustSizeSamplingIntervals;nInterval++) {
    //interpolate dust number density
    for (c0=0.0,i=0;i<nInterpolationCoeficients;i++) {
      c=InterpolationCoeficients[i];
      SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
      Measure=InterpolationList[i]->Measure;

      if (Measure<=0.0) exit(__LINE__,__FILE__,"Error: non-positive cell volume is found");
      if (PIC::LastSampleLength!=0) c0+=c*(*(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset))/PIC::LastSampleLength/Measure);
    }

    *(nInterval+(double*)(CellNodeSamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset))=c0;

    //grains' bulk velocity
    for (idim=0;idim<DIM;idim++) {
      for (c0=0.0,c1=0.0,i=0;i<nInterpolationCoeficients;i++) {
        SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

        c0+=*(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset));
        c1+=*(idim+3*nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Velocity_Offset));
      }

      if (c0>0.0) c1/=c0;
      *(idim+3*nInterval+(double*)(CellNodeSamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Velocity_Offset))=c1;
    }

    //grain's speed
    for (c0=0.0,c1=0.0,i=0;i<nInterpolationCoeficients;i++) {
      SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

      c0+=*(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset));
      c1+=*(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Speed_Offset));
    }

    if (c0>0.0) c1/=c0;
    *(nInterval+(double*)(CellNodeSamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Speed_Offset))=c1;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    //dust grains' charge density
    for (c0=0.0,i=0;i<nInterpolationCoeficients;i++) {
      c=InterpolationCoeficients[i];
      SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
      Measure=InterpolationList[i]->Measure;

      if (PIC::LastSampleLength!=0) c0+=c*(*(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCherge_Offset))/PIC::LastSampleLength/Measure);
    }

    *(nInterval+(double*)(CellNodeSamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCherge_Offset))=c0;

    //dust grain's electric current density
    for (idim=0;idim<DIM;idim++) {
      for (c0=0.0,i=0;i<nInterpolationCoeficients;i++) {
        c=InterpolationCoeficients[i];
        SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
        Measure=InterpolationList[i]->Measure;

        if (PIC::LastSampleLength!=0) c0+=c*(*(idim+3*nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCurrent_Offset))/PIC::LastSampleLength/Measure);
      }

      *(idim+3*nInterval+(double*)(CellNodeSamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCurrent_Offset))=c0;
    }
#endif

  }



}


void ElectricallyChargedDust::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double TotalDustNumberDensity=0.0;
  double t,numberDensity;
  char *SamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;


//=======  DEBUG BEGIN  ==============
static long int nCallCounter=0;

nCallCounter++;

/*
if (nCallCounter==61723) {
  cout << __FILE__ << "%" << __LINE__ << endl;
}
*/
//=======  DEBUG END  ==============


  //total density
  if (pipe->ThisThread==CenterNodeThread) for (int s=_DUST_SPEC_;s<_DUST_SPEC_+GrainVelocityGroup::nGroups;s++) TotalDustNumberDensity+=CenterNode->GetNumberDensity(s);

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(TotalDustNumberDensity,CenterNodeThread);

    fprintf(fout,"%e ",TotalDustNumberDensity);
  }
  else pipe->send(TotalDustNumberDensity);

  //the total mean dust velocity vector and speed
  double MeanBulkVelocity[3]={0.0,0.0,0.0},Speed=0.0;

  if (pipe->ThisThread==CenterNodeThread) {
    double vGroup[3],w,wtot=0.0;

    for (int s=_DUST_SPEC_;s<_DUST_SPEC_+GrainVelocityGroup::nGroups;s++) {
      w=CenterNode->GetNumberDensity(s);
      wtot+=w;
      Speed+=CenterNode->GetMeanParticleSpeed(s)*w;

      CenterNode->GetBulkVelocity(vGroup,s);
      for (int idim=0;idim<DIM;idim++) MeanBulkVelocity[idim]+=vGroup[idim]*w;
    }

    if (wtot>0.0) {
      Speed/=wtot;
      for (int idim=0;idim<DIM;idim++) MeanBulkVelocity[idim]/=wtot;
    }
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv(Speed,CenterNodeThread);
      pipe->recv(MeanBulkVelocity,DIM,CenterNodeThread);
    }

    for (int idim=0;idim<3;idim++) fprintf(fout,"%e ",MeanBulkVelocity[idim]);
    fprintf(fout,"%e ",Speed);
  }
  else {
    pipe->send(Speed);
    pipe->send(MeanBulkVelocity,DIM);
  }




  //total dust elelctric charge and current densities
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  //total electric charge
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::TotalDustElectricChargeDensitySamplingOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //electric current density
  for (int idim=0;idim<3;idim++) {
    //grain's bulk velocity
    if (pipe->ThisThread==CenterNodeThread) {
      t=*(idim+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::TotalDustElectricCurrentDensitySamplingOffset));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }


#endif

  //macroscopic parameters for particular intervals
  for (int nInterval=0;nInterval<Sampling::nDustSizeSamplingIntervals;nInterval++) {
    //number density
    if (pipe->ThisThread==CenterNodeThread) {
      numberDensity= *(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_NumberDensity_Offset));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(numberDensity,CenterNodeThread);

      fprintf(fout,"%e ",numberDensity);
    }
    else pipe->send(numberDensity);

    for (int idim=0;idim<DIM;idim++) {
      //grain's bulk velocity
      if (pipe->ThisThread==CenterNodeThread) {
        t=*(idim+3*nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Velocity_Offset));
      }

      if (pipe->ThisThread==0) {
        if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

        fprintf(fout,"%e ",t);
      }
      else pipe->send(t);
    }

    //grain's speed
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_Speed_Offset));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    //electric charge density
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCherge_Offset));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);

    //electric current density
    for (int idim=0;idim<DIM;idim++) {
      //grain's bulk velocity
      if (pipe->ThisThread==CenterNodeThread) {
        t=*(idim+3*nInterval+(double*)(SamplingBuffer+ElectricallyChargedDust::Sampling::DustSizeSamplingInterval_ElectricCurrent_Offset));
      }

      if (pipe->ThisThread==0) {
        if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

        fprintf(fout,"%e ",t);
      }
      else pipe->send(t);
    }


#endif
  }
}




/*----------------------------------------------- Sampling Size Distribution Function : BEGIN -------------------------------------------*/
void ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::Init(double ProbeLocations[][DIM],int nProbeLocations,int nIntervals) {
  int idim,nProbe,i,j,k;

  if (ElectricallyChargedDust::minDustRadius<=0.0) exit(__LINE__,__FILE__,"Error: the dust model sould be initialied first");
  if (nSamplingIntervals!=0) exit(__LINE__,__FILE__,"Error: re-initialization of the sampling module");

  nSamplingLocations=nProbeLocations;
  nSamplingIntervals=nIntervals;
  dLog10DustRadius=log10(ElectricallyChargedDust::maxDustRadius/ElectricallyChargedDust::minDustRadius)/nSamplingIntervals;

  //calculate the sampling offsets
  NumberDensitySamplingOffset=0;
  SamplingIntervalDataLength+=1;

  VelocitySamplingOffset=SamplingIntervalDataLength;
  SamplingIntervalDataLength+=3;

  SpeedSamplingOffset=SamplingIntervalDataLength;
  SamplingIntervalDataLength+=1;

  #if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  ElectricChargeSamplingOffset=SamplingIntervalDataLength;
  SamplingIntervalDataLength+=1;

  ElectricCurentSamplingOffset=SamplingIntervalDataLength;
  SamplingIntervalDataLength+=3;

  AbsoluteElectricCurrentValueSamplingOffset=SamplingIntervalDataLength;
  SamplingIntervalDataLength+=1;
  #endif


  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nProbeLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nProbeLocations];
  SamplingLocations=new double* [nProbeLocations];
  SamplingLocations[0]=new double [DIM*nProbeLocations];

  SamplingBuffer=new double* [nProbeLocations];
  SamplingBuffer[0]=new double [nProbeLocations*SamplingIntervalDataLength*nSamplingIntervals];

  for (nProbe=1;nProbe<nProbeLocations;nProbe++) {
    SamplingLocations[nProbe]=SamplingLocations[nProbe-1]+DIM;
    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+SamplingIntervalDataLength*nSamplingIntervals;
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nProbeLocations;nProbe++) {
    for (idim=0;idim<DIM;idim++) SamplingLocations[nProbe][idim]=ProbeLocations[nProbe][idim];

    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

void ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::flushSamplingBuffers() {
  long int i,TotalDataLength=nSamplingLocations*SamplingIntervalDataLength*nSamplingIntervals;
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}

void ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::SampleDistributionFnction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,idim;
  double LocalParticleWeight,grainRadius,v[3],speed;
  int ParticleDataLength=PIC::ParticleBuffer::ParticleDataLength,grainRadiusInterval;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::ParticleBuffer::byte *ParticleData;
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  double *SamplingDataPointer=NULL;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  double grainElecticCharge;
#endif

  for (node=SampleNodes[0],nProbe=0;nProbe<nSamplingLocations;node=SampleNodes[++nProbe]) if ((block=node->block)!=NULL) {
    //    ptr=block->GetCenterNode(SampleLocalCellNumber[nProbe])->FirstCellParticle;
    int i,j,k;

    PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
    ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (ptr!=-1) {
      ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      spec=PIC::ParticleBuffer::GetI(ParticleData);

      //check if the particle is the dust grain
      if ((spec<_DUST_SPEC_)||(spec>=_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) {
        ptr=PIC::ParticleBuffer::GetNext(ParticleData);
        continue;
      }

      //sample internal parameters of the dust grain
      memcpy((void*)tempParticleData,(void*)ParticleData,ParticleDataLength);

      LocalParticleWeight=block->GetLocalParticleWeight(spec);
      LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)tempParticleData);

      grainRadius=GetGrainRadius((PIC::ParticleBuffer::byte*)tempParticleData);

      grainRadiusInterval=(int)(log10(grainRadius/minDustRadius)/dLog10DustRadius);
      if (grainRadiusInterval<0) grainRadiusInterval=0;
      if (grainRadiusInterval>=nSamplingIntervals) grainRadiusInterval=nSamplingIntervals-1;

      //the pointer to the sampling data
      SamplingDataPointer=SamplingBuffer[nProbe]+grainRadiusInterval*SamplingIntervalDataLength;

      //sample the particles' number
      SamplingDataPointer[NumberDensitySamplingOffset]+=LocalParticleWeight;

      //sample the particle velocity
      PIC::ParticleBuffer::GetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
      for (idim=0,speed=0.0;idim<3;idim++) {
        v[idim]*=LocalParticleWeight;
        speed+=pow(v[idim],2);

        SamplingDataPointer[VelocitySamplingOffset+idim]+=v[idim];
      }

      speed=sqrt(speed);
      SamplingDataPointer[SpeedSamplingOffset]+=speed;

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
      //electric charge
      grainElecticCharge=GetGrainCharge((PIC::ParticleBuffer::byte*)tempParticleData);
      SamplingDataPointer[ElectricChargeSamplingOffset]+=grainElecticCharge*LocalParticleWeight;

      for (idim=0;idim<3;idim++) {
        SamplingDataPointer[ElectricCurentSamplingOffset+idim]+=v[idim]*grainElecticCharge;
      }

      SamplingDataPointer[AbsoluteElectricCurrentValueSamplingOffset]+=speed*grainElecticCharge;
#endif

      ptr=PIC::ParticleBuffer::GetNext(ParticleData);
    }
  }
}



void ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::printDistributionFunction(int DataOutputFileNumber) {
  int idim,nProbe,nVariable,thread,offset,iInterval;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,c;
  char str[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);


  for (nProbe=0;nProbe<nSamplingLocations;nProbe++) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      sprintf(str,"%s/pic.DUST.SizeDistribution.out=%i.nSamplePoint=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber,nProbe);
      fout=fopen(str,"w");

      fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);

      fprintf(fout,"\"TITLE=Dust size distribution function at x=%e",SamplingLocations[nProbe][0]);
      for (idim=1;idim<DIM;idim++) fprintf(fout,", %e",SamplingLocations[nProbe][idim]);

      fprintf(fout,"\"\nVARIABLES=\"Dust GrainRadius\", \"f_size(a)\",  \"Vx(a)\",\"Vy(a)\",\"Vz(a)\",\"Speed(a)\"");

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
      fprintf(fout,"    \"Electric Charge(a)\",  \"Vx times Electric Charge (a)\",\"Vy times Electric Charge (a)\",\"Vz times Electric Charge (a)\",\"Speed times Electric Charge (a)\"\n");
#else
      fprintf(fout,"\n");
#endif


      //collect the sampled information from other processors
      for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (iInterval=0;iInterval<nSamplingIntervals;iInterval++)  {
        offset=iInterval*SamplingIntervalDataLength;

        for (nVariable=0;nVariable<SamplingIntervalDataLength;nVariable++) SamplingBuffer[nProbe][nVariable+offset]+=pipe.recv<double>(thread);
      }



      //determine the mean values of the size-dependent dust parameters per dust grain
      //size dependent velocity
      for (iInterval=0;iInterval<nSamplingIntervals;iInterval++) {
        double t=SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+0];

        if (t>0.0) for (nVariable=1;nVariable<SamplingIntervalDataLength;nVariable++) SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+nVariable]/=t;
      }

      //diffirentiate the sample of grains's sizes to get the distribution function
      for (iInterval=0;iInterval<nSamplingIntervals;iInterval++) {
        double dl=minDustRadius*(pow(10,(iInterval+1)*dLog10DustRadius)-pow(10,iInterval*dLog10DustRadius));

        SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+0]/=dl;
      }


      //normalize the distribution fucntion
      for (norm=0.0,iInterval=0;iInterval<nSamplingIntervals;iInterval++)  {
        double dl=minDustRadius*(pow(10,(iInterval+1)*dLog10DustRadius)-pow(10,iInterval*dLog10DustRadius));

        norm+=SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+0]*dl;
      }

      if (norm>0.0) for (iInterval=0;iInterval<nSamplingIntervals;iInterval++) SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+0]/=norm;


      //print the output file
      for (iInterval=0;iInterval<nSamplingIntervals+1;iInterval++) {
        fprintf(fout,"%e  ",minDustRadius*pow(10,iInterval*dLog10DustRadius));

        for (nVariable=0;nVariable<SamplingIntervalDataLength;nVariable++) {
          if (iInterval==0) c=SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+nVariable];
          else if (iInterval<nSamplingIntervals) {
            c=0.5*(SamplingBuffer[nProbe][(iInterval-1)*SamplingIntervalDataLength+nVariable]+SamplingBuffer[nProbe][iInterval*SamplingIntervalDataLength+nVariable]);
          }
          else c=SamplingBuffer[nProbe][(iInterval-1)*SamplingIntervalDataLength+nVariable];

          fprintf(fout,"%e  ",c);
        }

        fprintf(fout,"\n");
      }

      //close the output file
      fclose(fout);
      fprintf(PIC::DiagnospticMessageStream,"done.\n");
    }
    else {
      for (iInterval=0;iInterval<nSamplingIntervals;iInterval++)  {
        offset=iInterval*SamplingIntervalDataLength;

        for (nVariable=0;nVariable<SamplingIntervalDataLength;nVariable++) pipe.send(SamplingBuffer[nProbe][nVariable+offset]);
      }
    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}


/*----------------------------------------------- Sampling Size Distribution Function : END -------------------------------------------*/


#endif




