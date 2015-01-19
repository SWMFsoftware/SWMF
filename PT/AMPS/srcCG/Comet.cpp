/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Comet";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";

int Comet::GravityFieldOffset=-1;


/*
int Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=-1;
int Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=-1;
int Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=-1;
*/

static bool probabilityFunctionDefinedJet=false,probabilityFunctionDefined=false,probabilityFunctionDefinedWaist=false,probabilityFunctionDefinedHartley2=false,probabilityFunctionDefinedUniformNASTRAN=false;
//static double productionDistributionJet[360][180],cumulativeProductionDistributionJet[360][180];
static double productionDistributionJet[6000],cumulativeProductionDistributionJet[6000];
static double productionDistributionWaist[6000],cumulativeProductionDistributionWaist[6000];
static double productionDistributionHartley2[6000],cumulativeProductionDistributionHartley2[6000];
static double productionDistributionUniformNASTRAN[200000],cumulativeProductionDistributionUniformNASTRAN[200000];
#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_
static double productionDistributionJetNASTRAN[3][200000],cumulativeProductionDistributionJetNASTRAN[3][200000],fluxDFMS[3][200000];
static bool definedFluxDFMS[3],probabilityFunctionDefinedJetNASTRAN[3];
static double DFMSproduction[3];
#else
static double productionDistributionJetNASTRAN[200000],cumulativeProductionDistributionJetNASTRAN[200000];
static bool probabilityFunctionDefinedJetNASTRAN=false;
#endif
static double productionDistribution[180],cumulativeProductionDistribution[180];
static double angle;
static double azimuthCenter;
static double zenithCenter;
static cInternalRotationBodyData* Nucleus;

double subSolarPointAzimuth=0.0; //53.0*Pi/180; //0.0;
double subSolarPointZenith=0.0;

double DustSizeMin=1.0e-7;
double DustSizeMax=1.0e-2;
double DustTotalMassProductionRate=0.0;
int DustSampleIntervals=10;
double DustSizeDistribution=0.0;

double fluxBjorn[90];
double nightSideFlux;

#if _DFMS_RATIO_MODE_ == _DFMS_RATIO_MODE_ON_
static double ratioBjornSpec[3][200000];
static double productionDistributionNASTRAN[3][200000],cumulativeProductionDistributionNASTRAN[3][200000],fluxBjornDFMS[3][200000];
static bool definedFluxBjornDFMS[3],probabilityFunctionDefinedNASTRAN[3];
static double DFMSBjornProduction[3];
#else
static double productionDistributionNASTRAN[200000],cumulativeProductionDistributionNASTRAN[200000];
static bool probabilityFunctionDefinedNASTRAN=false;
#endif

static double gravityAngle[200000];
static bool gravityAngleInitialized=false;

long int offsetSurfaceElement;

void Comet::Init_BeforeParser() {
#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  // Keep track of original surface element the particle was created from
  PIC::ParticleBuffer::RequestDataStorage(offsetSurfaceElement,sizeof(int));
#endif


#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
#endif

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model                                                                                                                                                                           
  ElectricallyChargedDust::minDustRadius=DustSizeMin; //0.1*_MICROMETER_;
  ElectricallyChargedDust::maxDustRadius=DustSizeMax; //1.0e4*_MICROMETER_;
  ElectricallyChargedDust::Sampling::SetDustSamplingIntervals(DustSampleIntervals);
  ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity=0.01;
  ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity=4.0;
  ElectricallyChargedDust::TotalMassDustProductionRate=DustTotalMassProductionRate;
  ElectricallyChargedDust::SizeDistribution::PowerIndex=DustSizeDistribution;
  ElectricallyChargedDust::Init_BeforeParser();
#endif

}

void Comet::Init_AfterParser() {

  /*  //set up the Chamberlen model
  double ExosphereEscapeRate[PIC::nTotalSpecies],ExospehreTemsprature[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) { 
    ExosphereEscapeRate[spec]=1.0e25,ExospehreTemsprature[spec]=1000.0;
    //ExosphereEscapeRate[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceRate[spec];
    //ExospehreTemsprature[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceTemeprature[spec];
    }*/

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model                                                                                                                                                                           
  ElectricallyChargedDust::Init_AfterParser();
#endif
  
  //init Gravity
#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  InitGravityData();
#endif

  
  //  Exosphere::ChamberlainExosphere::Init(ExosphereEscapeRate,ExospehreTemsprature);

}

int Comet::RequestDataBuffer(int offset) {
  int TotalDataLength;

  GravityFieldOffset=offset;
  TotalDataLength=3;

  return TotalDataLength*sizeof(double);
}

void Comet::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"Gx\",\"Gy\",\"Gz\"");
}

void Comet::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  //Gravity Field
  for (idim=0;idim<3;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}

void Comet::InitGravityData(){
  int thread,i,j,k,idim,offset,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  double gravityAccl[3],*position;

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
	  cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
	  if (cell!=NULL) {  
	    position=cell->GetX();
	    nucleusGravity::gravity(gravityAccl,position);
	    for (idim=0;idim<3;idim++) {
	      *((double*)(cell->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)))=gravityAccl[idim];
	      gravityAccl[idim]=0.0;
	    }
	  }
	}
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
      block=node->block;

      for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
	for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
	  for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
	    cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
	    if (cell!=NULL) {
	      position=cell->GetX();
	      nucleusGravity::gravity(gravityAccl,position);
	      for (idim=0;idim<3;idim++) {
		*((double*)(cell->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)))=gravityAccl[idim];
		gravityAccl[idim]=0.0;
	      }
	    }
	  }
      }
    }


  return ;
}

void Comet::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double G[3]={0.0,0.0,0.0};
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {

    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+GravityFieldOffset;idim<3;idim++) G[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+GravityFieldOffset,G,3*sizeof(double));
}

void Comet::GetGravityAcceleration(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  register int idim;
  register double *offset=(double*)(GravityFieldOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

  for (idim=0;idim<3;idim++) x[idim]=offset[idim];

}				   

double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
#if _COMET_TEMPERATURE_MODE_ == _COMET_TEMPERATURE_MODE__BJORN_
  const double DistanceFromTheSun[6]={1.3,2.0,2.7,3.0,3.25,3.5};
  const double minTemp[6]={172.0,163.0,150.0,145.0,139.0,133.0};
  double res,r,zenith,azimuth;
  int angle;
  double angletemp;

  //determine if the point on the night side of the Comet 
  if (CosSubSolarAngle<0.0) return minTemp[Comet::ndist];

  if (Comet::ndist<5) {
    //return the day-side temeprature
    angle=(int) (acos(CosSubSolarAngle)*180.0/Pi);
    if(angle>89) angle=89;
    res=(SurfaceTemp[angle][Comet::ndist+1]>minTemp[Comet::ndist]) ?  SurfaceTemp[angle][Comet::ndist+1] : minTemp[Comet::ndist];
  }else{
    angletemp=acos(CosSubSolarAngle);
    angle=(int) (acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[Comet::ndist],2.0)*cos(angletemp))*180/Pi);
    if(angle>89) angle=89;
    res=(SurfaceTemp[angle][4+1]>minTemp[Comet::ndist]) ?  SurfaceTemp[angle][4+1] : minTemp[Comet::ndist];
  }
  return res;

#elif _COMET_TEMPERATURE_MODE_ == _COMET_TEMPERATURE_MODE__CONSTANT_
  double surfaceTemperatureConstant=180.0;
  return surfaceTemperatureConstant;

#elif _COMET_TEMPERATURE_MODE_ == _COMET_TEMPERATURE_MODE__ANALYTICAL_
  double Tmax=182.09922,T75=136.11149,Tmin=133;
  double res;

  double angle,a,b;

  a=(Tmax-T75)/(1-1/cos(75*Pi/180.0));
  b=Tmax-a;
  
  angle=acos(CosSubSolarAngle);
  
  res=(angle<=75*Pi/180.0) ? max(Tmin,a/cos(angle)+b) : Tmin; 
  return res;
#else
  exit(__LINE__,__FILE__,"Temperature mode does not exist");
#endif
}

long int Comet::InjectionBoundaryModel_Limited() {
  int spec;
  long int res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=InjectionBoundaryModel_Limited(spec);

  return res;
}

long int Comet::InjectionBoundaryModel_Limited(int spec) {
  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_SO_OBJECT[3],x_IAU_OBJECT[3],v_SO_OBJECT[3],v_IAU_OBJECT[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag=false;
  int SourceProcessID;

  double totalProductionRate=Comet::GetTotalProductionRateBjornNASTRAN(spec)+Comet::GetTotalProductionRateUniformNASTRAN(spec)+Comet::GetTotalProductionRateJetNASTRAN(spec);

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    LocalTimeStep=Comet::CG->maxIntersectedNodeTimeStep[spec];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  ModelParticlesInjectionRate=totalProductionRate/ParticleWeight;

  //  if (ModelParticlesInjectionRate*ParticleWeight*LocalTimeStep<1.0E-10) return 0;

  if (ModelParticlesInjectionRate*LocalTimeStep>nMaxInjectedParticles) {
    ParticleWeightCorrection=ModelParticlesInjectionRate*LocalTimeStep/nMaxInjectedParticles;
    ModelParticlesInjectionRate/=ParticleWeightCorrection;
  }

  //definition of indexes TEMPORARY!!!!!
  //int _EXOSPHERE__SOURCE_MAX_ID_VALUE_=2;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_=0;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_=1;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_=2;
  int _exosphere__SOURCE_MAX_ID_VALUE_=2;


  //calcualte probabilities of each source processes                                                                 
  double TotalFlux,FluxSourceProcess[1+_exosphere__SOURCE_MAX_ID_VALUE_]; //,ProbabilitySourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];                                                                                               
int iSource;

for (iSource=0;iSource<1+_exosphere__SOURCE_MAX_ID_VALUE_;iSource++) FluxSourceProcess[iSource]=0.0; //,ProbabilitySourceProcess[iSource]=0.0;                                                                                         

TotalFlux=totalProductionRate;

//only Used defined source here since we only want the Bjorn model so far
//calculate the source rate due to user defined source functions                                                   
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=Comet::GetTotalProductionRateBjornNASTRAN(spec);
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]=Comet::GetTotalProductionRateUniformNASTRAN(spec);
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_]=Comet::GetTotalProductionRateJetNASTRAN(spec);

//Distribution of dust injection correlated with water
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
 if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=Comet::GetTotalProductionRateBjornNASTRAN(_H2O_SPEC_);
   FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]=Comet::GetTotalProductionRateUniformNASTRAN(_H2O_SPEC_);
   FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_]=Comet::GetTotalProductionRateJetNASTRAN(_H2O_SPEC_);
 }
 TotalFlux=Comet::GetTotalProductionRateBjornNASTRAN(_H2O_SPEC_)+Comet::GetTotalProductionRateUniformNASTRAN(_H2O_SPEC_)+Comet::GetTotalProductionRateJetNASTRAN(_H2O_SPEC_);;
#endif

 double CalculatedSourceRate[PIC::nTotalSpecies][1+_exosphere__SOURCE_MAX_ID_VALUE_];
 CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=0.0;
 CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]=0.0;
 CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_]=0.0;


#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
 if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   static double GrainInjectedMass=0.0;
   PIC::Mesh::cDataBlockAMR *block;
   double GrainRadius,GrainMass,GrainWeightCorrection;
   int GrainVelocityGroup;

   GrainInjectedMass+=ElectricallyChargedDust::TotalMassDustProductionRate*LocalTimeStep;

   while (GrainInjectedMass>0.0) {
     startNode=NULL;

     do {
       SourceProcessID=(int)(rnd()*(1+_exosphere__SOURCE_MAX_ID_VALUE_));
     }
     while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());

     //generate a particle                                                                                             
     char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
     PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

     if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_) {
       flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
       GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
       GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
       SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
       if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
     }
     else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_) {
       flag=Comet::GenerateParticlePropertiesUniformNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
       GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
       GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
       SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_;
       if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
     }
     else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_) {
       flag=Comet::GenerateParticlePropertiesJetNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
       GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
       GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
       SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_;
       if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
     }
     else {
       continue;
     }
     if (flag==false) continue;
     if ((block=startNode->block)->GetLocalTimeStep(_DUST_SPEC_)/LocalTimeStep<rnd()) continue;

     //determine the velocity group of the injected grain;                                                                                                                                                                                                                     
     //calculate additional particle weight correction because the particle will be placed in a different weight group                                                                                                                                                         
     GrainVelocityGroup=ElectricallyChargedDust::GrainVelocityGroup::GetGroupNumber(v_SO_OBJECT);
     GrainWeightCorrection*=block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/block->GetLocalTimeStep(_DUST_SPEC_);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
     GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup];
#else
     exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

    //determine the surface element of the particle origin                                                            
     PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetI(_DUST_SPEC_+GrainVelocityGroup,(PIC::ParticleBuffer::byte*)tempParticleData);

  ElectricallyChargedDust::SetGrainCharge(0.0,(PIC::ParticleBuffer::byte*)tempParticleData);
  ElectricallyChargedDust::SetGrainMass(GrainMass,(PIC::ParticleBuffer::byte*)tempParticleData);
  ElectricallyChargedDust::SetGrainRadius(GrainRadius,(PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);


  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  nInjectedParticles++;

  //inject the particle into the system                                                                             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
  if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is n\
ot defined");
#endif

  _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
 }
 
 }else{
#endif
while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
  //determine the source process to generate a particle's properties                                               
  do {
    SourceProcessID=(int)(rnd()*(1+_exosphere__SOURCE_MAX_ID_VALUE_));
  }
  while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());

  //generate a particle                                                                                             
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

  //to satisfy the compiler and fit the while structure                                                             
  if (false) {}

  //Add the user defined particle gineration                                                                        
  else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_) {
    flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
    SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
    if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  }

  else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_) {
    flag=Comet::GenerateParticlePropertiesUniformNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
    SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_;
    if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  }  

  else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_) {
    flag=Comet::GenerateParticlePropertiesJetNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
    SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_;
    if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  }  

  else {
    continue;
  }

  if (flag==false) continue;

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  if (startNode->block->GetLocalTimeStep(spec)/LocalTimeStep<rnd()) continue;
#endif
 
  //determine the surface element of the particle origin                                                            
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);


  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  nInjectedParticles++;

  //inject the particle into the system                                                                             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
  if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
  if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is n\
ot defined");
#endif

  _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
 }

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
}
#endif

 return nInjectedParticles;
}

double Comet::GetTotalProductionRateBjornNASTRAN(int spec){
  double rSphere=1980,c=0.0,X=0.0,totalProductionRate=0.0;
  double positionSun[3],x[3],norm[3];
  long int totalSurfaceElementsNumber,i;
  double percentageActive=0.05;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  const double DistanceFromTheSun[6]={1.3,2.0,2.7,3.0,3.25,3.5};
  double HeliocentricDistance=3.3*_AU_;
  int idim;
  double ProductionRateScaleFactor,TableTotalProductionRate=0.0,BjornTotalProductionRate;
  int angle;
  double angletemp;

#if _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ ==  _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ANALYTICAL_
  double Qmin=4.5e17,Qmax=6.3e18;

  for (i=0;i<90;i++) {
    angle=(double) i;
    fluxBjorn[i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
  }

  nightSideFlux=Qmin;  

#else

  if (Comet::ndist<5) {
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      TableTotalProductionRate+=ProductionRate[i][2+ndist];
    }
  }else{
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      angletemp=(double) i;
      angletemp*=Pi/180.0;
      angle=(int) (acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[ndist],2.0)*cos(angletemp))*180/Pi);
      TableTotalProductionRate+=ProductionRate[angle][2+4]/ProductionRate[angle][1]*ProductionRate[i][1];
    }
  }

  ProductionRateScaleFactor=(1.0-2.0*NightSideProduction[Comet::ndist]);

  BjornTotalProductionRate=TableTotalProductionRate/ProductionRateScaleFactor;

  if (Comet::ndist<5) {
    for (i=0;i<90;i++) {
      fluxBjorn[i]=(ProductionRateScaleFactor/TableTotalProductionRate*BjornTotalProductionRate*ProductionRate[i][2+Comet::ndist]/ProductionRate[i][1]+2.0*NightSideProduction[Comet::ndist]*BjornTotalProductionRate/(4*Pi*rSphere*rSphere))*percentageActive;
    }
  }else{
    for (i=0;i<90;i++) {
      angletemp=(double) i;
      angletemp*=Pi/180.0;
      angle=(int) (acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[ndist],2.0)*cos(angletemp))*180/Pi);
      fluxBjorn[i]=(ProductionRateScaleFactor/TableTotalProductionRate*BjornTotalProductionRate*ProductionRate[angle][2+4]/ProductionRate[angle][1]+2.0*NightSideProduction[Comet::ndist]*BjornTotalProductionRate/(4*Pi*rSphere*rSphere))*percentageActive;
    }
  }

  nightSideFlux=BjornTotalProductionRate*NightSideProduction[Comet::ndist]/(2*Pi*rSphere*rSphere)*percentageActive;
#endif

#if _DFMS_RATIO_MODE_ == _DFMS_RATIO_MODE_ON_
  //if (spec==_CO_SPEC_ || spec==_CO2_SPEC_) {
  if (definedFluxBjornDFMS[spec]==false) {

    if (_H2O_SPEC_!=0) exit(__LINE__,__FILE__,"Error: H2O needs to have the index 0 for this mode to work");

    long int File_Header = 0;
    static double **Data=NULL;
    long int Data_length;
    
    FILE *fH;
    fH = fopen("majorSpeciesFromDFMS.txt","r");
    
    //readDATAlength
    char str[10000];
    Data_length=0;
    rewind(fH);
    while (!feof(fH)){
      fgets(str,10000,fH);
      Data_length++;
    }
    Data_length -= (1+File_Header);
    
    if (PIC::ThisThread==0) printf("Data_length: %li \n",Data_length);
    
    Data = new double* [Data_length];
    Data[0] = new double [Data_length*9];
    for (int i=0;i<Data_length;i++) { 
      Data[i]=Data[0]+i*9;  
      for (int j=0;j<9;j++) Data[i][j]=0.0;
    }
    
    //read DATA
    long int nline,i=0;
    double f1,f2,f3,f4,f5,f6,f7,f8,f9;
    char strH[10000];
    
    rewind(fH);
    
    // Header
    for (nline=0;nline<File_Header;nline++) {
      fgets(strH,10000,fH);
    }
    // Data
    for(i=0;i<Data_length;i++,nline++) {
      fscanf(fH,"%le%le%le%le%le%le%le%le%le\n",&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8,&f9);
      Data[i][0]=f1;
      Data[i][1]=f2;
      Data[i][2]=f3;
      Data[i][3]=f4;
      Data[i][4]=f5;
      Data[i][5]=f6;
      Data[i][6]=f7;
      Data[i][7]=f8;
      Data[i][8]=f9;
    }
    
    int ct=0,idim;
    double scalar=0.0,normal=0.0,xSpacecraft[3],angle,r=0.0,x[3];
    
    //initialization to zero
    for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) fluxBjornDFMS[spec][i]=0.0;
    
    for (ct=0;ct<Data_length;ct++) {    
      xSpacecraft[0]=cos(Data[ct][3]*Pi/180.0)*cos(Data[ct][4]*Pi/180.0);
      xSpacecraft[1]=sin(Data[ct][3]*Pi/180.0)*cos(Data[ct][4]*Pi/180.0);
      xSpacecraft[2]=sin(Data[ct][4]*Pi/180.0);
      
      for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
	CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);
	scalar=0.0,normal=0.0,r=0.0;
	
	for (idim=0;idim<3;idim++) {
	  scalar+=x[idim]*xSpacecraft[idim];
	  normal+=pow(x[idim],2.0);
	  r+=pow(x[idim],2.0);
	}
	normal=sqrt(normal);
	angle=acos(scalar/normal)*180.0/Pi;
	r=sqrt(r);
	
	if (spec==_H2O_SPEC_) {
	  if(angle<5.0 && Data[ct][6]>0.0 && Data[ct][6]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<7.0e26) fluxBjornDFMS[spec][i]=Data[ct][6]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
	}
	else if (spec==_CO_SPEC_) {
	  if(angle<5.0 && Data[ct][7]>0.0 && Data[ct][7]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<1.0e26) fluxBjornDFMS[spec][i]=Data[ct][7]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
	}
	else if (spec==_CO2_SPEC_) {
	  if(angle<5.0 && Data[ct][8]>0.0 && Data[ct][8]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<1.0e26) fluxBjornDFMS[spec][i]=Data[ct][8]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
	}
      }
    }
    
    for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
      if (fluxBjornDFMS[spec][i]==0.0) {
	if (spec==_H2O_SPEC_) fluxBjornDFMS[spec][i]=1.0e14;
	if (spec==_CO_SPEC_) fluxBjornDFMS[spec][i]=1.0e13;
	if (spec==_CO2_SPEC_) fluxBjornDFMS[spec][i]=1.0e13;
      }
    }
    
    for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) ratioBjornSpec[spec][i]=(fluxBjornDFMS[spec][i]/fluxBjornDFMS[_H2O_SPEC_][i]<10) ? fluxBjornDFMS[spec][i]/fluxBjornDFMS[_H2O_SPEC_][i]:1.0;

    DFMSBjornProduction[spec]=0.0;
    
    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);
    
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS); //I had middle element on body rotation...
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
	c+=norm[idim]*(positionSun[idim]-x[idim]);
	X+=pow(positionSun[idim]-x[idim],2.0);
      }
      
      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow
	DFMSBjornProduction[spec]+=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea*ratioBjornSpec[spec][i];
      }else{
	double angleProd;
	int angleProdInt;
	angleProd=acos(c/sqrt(X))*180/Pi;
	angleProdInt=(int) angleProd;
	DFMSBjornProduction[spec]+=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*ratioBjornSpec[spec][i];
      }
    }
    if (PIC::ThisThread==0) printf("totalProductionRate[%i]=%e flux(0)=%e nightsideFlux=%e \n",spec,DFMSBjornProduction[spec],fluxBjorn[0],nightSideFlux);
    
    definedFluxBjornDFMS[spec]=true;
  }
  return DFMSBjornProduction[spec];

  
#endif
  

#if _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ ==  _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ON_
  return Comet::Bjorn_SourceRate[spec];

#else

  positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

  totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;

  for (i=0;i<totalSurfaceElementsNumber;i++) {
    for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS); //I had middle element on body rotation...
    for (c=0.0,X=0.0,idim=0;idim<3;idim++){
      c+=norm[idim]*(positionSun[idim]-x[idim]);
      X+=pow(positionSun[idim]-x[idim],2.0);
    }

    if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow
      totalProductionRate+=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
    }else{
      double angleProd;
      int angleProdInt;
      angleProd=acos(c/sqrt(X))*180/Pi;
      angleProdInt=(int) angleProd;
      totalProductionRate+=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
    }
  }
  if (probabilityFunctionDefinedNASTRAN==false)  printf("totalProductionRate=%e flux(0)=%e nightsideFlux=%e \n",totalProductionRate,fluxBjorn[0],nightSideFlux);

  return totalProductionRate;

#endif

}

#if _DFMS_RATIO_MODE_ == _DFMS_RATIO_MODE_ON_
bool Comet::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.5*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  double area;
  double totalProdNightSide=0.0,totalProdDaySide=0.0,scalingFactor,scalingFactorDay,totalSurfaceInShadow=0.0,totalSurfaceInDayLight=0.0;


  if (probabilityFunctionDefinedNASTRAN[spec]==false) {

    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;

    total=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS);
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }
      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
        productionDistributionNASTRAN[spec][i]=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea*ratioBjornSpec[spec][i];
        total+=productionDistributionNASTRAN[spec][i];
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        productionDistributionNASTRAN[spec][i]=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*ratioBjornSpec[spec][i];
        total+=productionDistributionNASTRAN[spec][i];
      }
    }
   
    cumulativeProductionDistributionNASTRAN[spec][0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionNASTRAN[spec][i]+=productionDistributionNASTRAN[spec][i]/total;
      }else{
	cumulativeProductionDistributionNASTRAN[spec][i]=cumulativeProductionDistributionNASTRAN[spec][i-1]+productionDistributionNASTRAN[spec][i]/total;
      }
      
    }
    
    /*    FILE *out;
    out = fopen("GasFlux_%i.dat","w",spec);
    fprintf(out,"Angle Flux(m-2.s-1) \n");
    for (i=0;i<90;i++) {
      fprintf(out,"%e %e \n",ProductionRate[i][0],fluxBjorn[i]);
    }
    fclose(out);

    FILE *fout;
    fout = fopen("SurfTemp_%i.dat","w",spec);
    fprintf(fout,"Angle Temp (K) \n");
    double angletemp;
    for (i=0;i<90;i++) {
      angletemp=(double) i*Pi/180.0;
      fprintf(fout,"%e %e \n",ProductionRate[i][0],Exosphere::GetSurfaceTemeprature(cos(angletemp),x));
    }
    fclose(fout);
    */
    
    probabilityFunctionDefinedNASTRAN[spec]=true;
  }
  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
  while (gamma>cumulativeProductionDistributionNASTRAN[spec][i]){
    i++;
  }

#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  Comet::SetParticleSurfaceElement(i,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif
    
  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh.EPS);
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,ExternalNormal,1.0); //1e-4

  for (idim=0;idim<3;idim++) ExternalNormal[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};
  if(CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  double r2Tang=0.0;;
  double xFace[3];
  double vDustInit=0.01;
  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  #else
  double vInit=500;
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vInit*ExternalNormal[idim];
  #endif

#endif
  
  //init the internal degrees of freedom if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::InitRotTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
  PIC::IDF::InitVibTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}
#else
bool Comet::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.5*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  double area;
  double totalProdNightSide=0.0,totalProdDaySide=0.0,scalingFactor,scalingFactorDay,totalSurfaceInShadow=0.0,totalSurfaceInDayLight=0.0;


  if (probabilityFunctionDefinedNASTRAN==false) {

    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;

    total=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS);
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }
      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
        productionDistributionNASTRAN[i]=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
        total+=productionDistributionNASTRAN[i];
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        productionDistributionNASTRAN[i]=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
        total+=productionDistributionNASTRAN[i];
      }
    }
   
    cumulativeProductionDistributionNASTRAN[0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionNASTRAN[i]+=productionDistributionNASTRAN[i]/total;
      }else{
	cumulativeProductionDistributionNASTRAN[i]=cumulativeProductionDistributionNASTRAN[i-1]+productionDistributionNASTRAN[i]/total;
      }
      
    }
    /*
    FILE *out;
    out = fopen("GasFlux.dat","w");
    fprintf(out,"Angle Flux(m-2.s-1) \n");
    for (i=0;i<90;i++) {
      fprintf(out,"%e %e \n",ProductionRate[i][0],fluxBjorn[i]);
    }
    fclose(out);

    FILE *fout;
    fout = fopen("SurfTemp.dat","w");
    fprintf(fout,"Angle Temp (K) \n");
    double angletemp;
    for (i=0;i<90;i++) {
      angletemp=(double) i*Pi/180.0;
      fprintf(fout,"%e %e \n",ProductionRate[i][0],Exosphere::GetSurfaceTemeprature(cos(angletemp),x));
    }
    fclose(fout);
    
    */
    probabilityFunctionDefinedNASTRAN=true;
  }
  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
  while (gamma>cumulativeProductionDistributionNASTRAN[i]){
    i++;
  }

#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  Comet::SetParticleSurfaceElement(i,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif
    
  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh.EPS);
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,ExternalNormal,1.0); //1e-4

  for (idim=0;idim<3;idim++) ExternalNormal[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};
  if(CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  double r2Tang=0.0;;
  double xFace[3];
  double vDustInit=0.01;
  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  #else
  double vInit=500;
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vInit*ExternalNormal[idim];
  #endif

#endif
  
  //init the internal degrees of freedom if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::InitRotTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
  PIC::IDF::InitVibTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}
#endif

double Comet::GetTotalProductionRateUniformNASTRAN(int spec){
  return Comet::Uniform_SourceRate[spec];
}

bool Comet::GenerateParticlePropertiesUniformNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.3*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  double area;

  if (probabilityFunctionDefinedUniformNASTRAN==false) {       
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      TableTotalProductionRate+=ProductionRate[i][2+Comet::ndist];
    }
    
    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);
    
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    
    total=0.0;      
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS); //I had middle element on body rotation...
      
      productionDistributionUniformNASTRAN[i]=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
      total+=productionDistributionUniformNASTRAN[i];
    }
      
    cumulativeProductionDistributionUniformNASTRAN[0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionUniformNASTRAN[i]+=productionDistributionUniformNASTRAN[i]/total;
      }else{
	cumulativeProductionDistributionUniformNASTRAN[i]=cumulativeProductionDistributionUniformNASTRAN[i-1]+productionDistributionUniformNASTRAN[i]/total;
      }
      
    }
    probabilityFunctionDefinedUniformNASTRAN=true;
  }
  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
  while (gamma>cumulativeProductionDistributionUniformNASTRAN[i]){
    i++;
  }

#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  Comet::SetParticleSurfaceElement(i,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif
    
  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh.EPS);
  //    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,ExternalNormal,1.0);

  for (idim=0;idim<3;idim++) ExternalNormal[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};
  if(CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  double r2Tang=0.0;
  double xFace[3];
  double vDustInit=0.01;
  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  #else
  double vInit=500;
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vInit*ExternalNormal[idim];
  #endif

#endif
  
  //init the internal degrees of freedom if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::InitRotTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
  PIC::IDF::InitVibTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}


double Comet::GetTotalProductionRateJetNASTRAN(int spec){
  /*  double rSphere=1980,c=0.0,X=0.0,totalProductionRate=0.0;
  double positionSun[3],x[3],norm[3];
  long int totalSurfaceElementsNumber,i;
  double percentageActive=1.0;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  const double DistanceFromTheSun[6]={1.3,2.0,2.7,3.0,3.25,3.5};
  double HeliocentricDistance=3.3*_AU_;
  int idim;
  double ProductionRateScaleFactor,TableTotalProductionRate=0.0,BjornTotalProductionRate;
  int angle;
  double angletemp;

  if (Comet::ndist<5) {
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      TableTotalProductionRate+=ProductionRate[i][2+ndist];
    }
  }else{
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      angletemp=(double) i;
      angletemp*=Pi/180.0;
      angle=(int) (acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[ndist],2.0)*cos(angletemp))*180/Pi);
      TableTotalProductionRate+=ProductionRate[angle][2+4]/ProductionRate[angle][1]*ProductionRate[i][1];
    }
  }

  ProductionRateScaleFactor=(1.0-2.0*NightSideProduction[Comet::ndist]);

  BjornTotalProductionRate=TableTotalProductionRate/ProductionRateScaleFactor;

  if (Comet::ndist<5) {
    for (i=0;i<90;i++) {
      fluxBjorn[i]=(ProductionRateScaleFactor/TableTotalProductionRate*BjornTotalProductionRate*ProductionRate[i][2+Comet::ndist]/ProductionRate[i][1]+2.0*NightSideProduction[Comet::ndist]*BjornTotalProductionRate/(4*Pi*rSphere*rSphere))*percentageActive;
    }
  }else{
    for (i=0;i<90;i++) {
      angletemp=(double) i;
      angletemp*=Pi/180.0;
      angle=(int) (acos(pow(DistanceFromTheSun[4]/DistanceFromTheSun[ndist],2.0)*cos(angletemp))*180/Pi);
      fluxBjorn[i]=(ProductionRateScaleFactor/TableTotalProductionRate*BjornTotalProductionRate*ProductionRate[angle][2+4]/ProductionRate[angle][1]+2.0*NightSideProduction[Comet::ndist]*BjornTotalProductionRate/(4*Pi*rSphere*rSphere))*percentageActive;
    }
  }

  nightSideFlux=BjornTotalProductionRate*NightSideProduction[Comet::ndist]/(2*Pi*rSphere*rSphere)*percentageActive;


  positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
  positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

  totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;

  for (i=0;i<totalSurfaceElementsNumber;i++) {
    for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
    //    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS); //I had middle element on body rotation...                                                                    
    CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);
    if (gravityAngleInitialized==false) {
      int nd,u,v,w;
      double accl[3]={0.0,0.0,0.0};
      double scalar=0.0,angle=0.0,normAcclVector=0.0;
      nucleusGravity::gravity(accl,x);
      for(idim=0;idim<3;idim++) {
        scalar+=-accl[idim]*norm[idim];
        normAcclVector+=pow(accl[idim],2.0);
      }
      normAcclVector=sqrt(normAcclVector);
      angle=acos(scalar/normAcclVector)*180.0/Pi;
      gravityAngle[i]=angle;
    }
    if(gravityAngle[i]>45.0) {
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }

      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow                                                                    
        totalProductionRate+=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        totalProductionRate+=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
      }
    }
  }
  if (gravityAngleInitialized==false) gravityAngleInitialized=true;
  if (probabilityFunctionDefinedJetNASTRAN==false && PIC::ThisThread==0)  printf("totalProductionRate=%e flux(0)=%e nightsideFlux=%e \n",totalProductionRate,fluxBjorn[0],nightSideFlux);

  return totalProductionRate;
  */
#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_
  if (definedFluxDFMS[spec]==false) {
  long int File_Header = 0;
  static double **Data=NULL;
  long int Data_length;
  
  FILE *fH;
  fH = fopen("majorSpeciesFromDFMS.txt","r");

  //readDATAlength
  char str[10000];
  Data_length=0;
  rewind(fH);
  while (!feof(fH)){
    fgets(str,10000,fH);
    Data_length++;
  }
  Data_length -= (1+File_Header);

  if (PIC::ThisThread==0) printf("Data_length: %li \n",Data_length);
  
  Data = new double* [Data_length];
  Data[0] = new double [Data_length*9];
  for (int i=0;i<Data_length;i++) { 
    Data[i]=Data[0]+i*9;  
    for (int j=0;j<9;j++) Data[i][j]=0.0;
  }

  //read DATA
  long int nline,i=0;
  double f1,f2,f3,f4,f5,f6,f7,f8,f9;
  char strH[10000];

  rewind(fH);

  // Header
  for (nline=0;nline<File_Header;nline++) {
    fgets(strH,10000,fH);
  }
  // Data
  for(i=0;i<Data_length;i++,nline++) {
    fscanf(fH,"%le%le%le%le%le%le%le%le%le\n",&f1,&f2,&f3,&f4,&f5,&f6,&f7,&f8,&f9);
    Data[i][0]=f1;
    Data[i][1]=f2;
    Data[i][2]=f3;
    Data[i][3]=f4;
    Data[i][4]=f5;
    Data[i][5]=f6;
    Data[i][6]=f7;
    Data[i][7]=f8;
    Data[i][8]=f9;
  }

  int ct=0,idim;
  double scalar=0.0,norm=0.0,xSpacecraft[3],angle,r=0.0,x[3];

  //initialization to zero
  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) fluxDFMS[spec][i]=0.0;

  for (ct=0;ct<Data_length;ct++) {    
    xSpacecraft[0]=cos(Data[ct][3]*Pi/180.0)*cos(Data[ct][4]*Pi/180.0);
    xSpacecraft[1]=sin(Data[ct][3]*Pi/180.0)*cos(Data[ct][4]*Pi/180.0);
    xSpacecraft[2]=sin(Data[ct][4]*Pi/180.0);

    for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);
      scalar=0.0,norm=0.0,r=0.0;

      for (idim=0;idim<3;idim++) {
	scalar+=x[idim]*xSpacecraft[idim];
	norm+=pow(x[idim],2.0);
	r+=pow(x[idim],2.0);
      }
      norm=sqrt(norm);
      angle=acos(scalar/norm)*180.0/Pi;
      r=sqrt(r);
      
      if (spec==_H2O_SPEC_) {
	if(angle<5.0 && Data[ct][6]>0.0 && Data[ct][6]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<7.0e26) fluxDFMS[spec][i]=Data[ct][6]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
      }
      else if (spec==_CO_SPEC_) {
	if(angle<5.0 && Data[ct][7]>0.0 && Data[ct][7]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<1.0e26) fluxDFMS[spec][i]=Data[ct][7]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
      }
      else if (spec==_CO2_SPEC_) {
	if(angle<5.0 && Data[ct][8]>0.0 && Data[ct][8]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0<1.0e26) fluxDFMS[spec][i]=Data[ct][8]*1.0e6*4*Pi*pow(Data[ct][1]*1000.0,2.0)*600.0/(4*Pi*r*r);
      }
    }
  }
  
  DFMSproduction[spec]=0.0;

  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
    if (fluxDFMS[spec][i]==0.0) {
      if (spec==_H2O_SPEC_) fluxDFMS[spec][i]=1.0e14;
      if (spec==_CO_SPEC_) fluxDFMS[spec][i]=1.0e13;
      if (spec==_CO2_SPEC_) fluxDFMS[spec][i]=1.0e13;
    }
  }

  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) DFMSproduction[spec]+=fluxDFMS[spec][i]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;

  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) productionDistributionJetNASTRAN[spec][i]=fluxDFMS[spec][i]/DFMSproduction[spec]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;

  definedFluxDFMS[spec]=true;

  if (PIC::ThisThread==0) printf("DFMSproduction[%i]=%e \n",spec,DFMSproduction[spec]);
  }
  return DFMSproduction[spec];
#else
  return Comet::Jet_SourceRate[spec];
#endif

}


#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_    
bool Comet::GenerateParticlePropertiesJetNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.3*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  double totalArea=0.0;
  double totalActiveArea=0.0;

  //  if (probabilityFunctionDefinedJetNASTRAN==false) {          
  if (probabilityFunctionDefinedJetNASTRAN[spec]==false) {          
    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);
    
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    /*
    total=0.0; 
    totalArea=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);

      if(gravityAngle[i]>45.0) {

        for (c=0.0,X=0.0,idim=0;idim<3;idim++){
          c+=norm[idim]*(positionSun[idim]-x[idim]);
          X+=pow(positionSun[idim]-x[idim],2.0);
        }
        if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
          productionDistributionJetNASTRAN[i]=nightSideFlux*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
          total+=productionDistributionJetNASTRAN[i];
          totalActiveArea+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
        }else{
          double angleProd;
          int angleProdInt;
          angleProd=acos(c/sqrt(X))*180/Pi;
          angleProdInt=(int) angleProd;
          productionDistributionJetNASTRAN[i]=fluxBjorn[angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea;
          total+=productionDistributionJetNASTRAN[i];
          totalActiveArea+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
        }
      }

      totalArea+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
    }
    
    if (PIC::ThisThread==0) printf("The ratio of the active area with respect to the total area is: %e \n",total/totalArea);

    cumulativeProductionDistributionJetNASTRAN[0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionJetNASTRAN[i]+=productionDistributionJetNASTRAN[i]/total;
      }else{
	cumulativeProductionDistributionJetNASTRAN[i]=cumulativeProductionDistributionJetNASTRAN[i-1]+productionDistributionJetNASTRAN[i]/total;
      }
      
    }
    */
#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_    
    cumulativeProductionDistributionJetNASTRAN[spec][0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionJetNASTRAN[spec][i]+=productionDistributionJetNASTRAN[spec][i];
      }else{
	cumulativeProductionDistributionJetNASTRAN[spec][i]=cumulativeProductionDistributionJetNASTRAN[spec][i-1]+productionDistributionJetNASTRAN[spec][i];
      }
    }
#endif

    probabilityFunctionDefinedJetNASTRAN[spec]=true;
    }

  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_    
  while (gamma>cumulativeProductionDistributionJetNASTRAN[spec][i]){
    i++;
  }
#else
  while (gamma>cumulativeProductionDistributionJetNASTRAN[i]){
    i++;
  }
#endif

#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  Comet::SetParticleSurfaceElement(i,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh.EPS);
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,ExternalNormal,1.0);

  for (idim=0;idim<3;idim++) ExternalNormal[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};
  if(CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  double r2Tang=0.0;;
  double xFace[3];
  double vDustInit=0.01;
  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else
  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  #else
  double vInit=500;
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vInit*ExternalNormal[idim];
  #endif
#endif
  
  //init the internal degrees of freedom if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
  PIC::IDF::InitRotTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
  PIC::IDF::InitVibTemp(SurfaceTemperature,(PIC::ParticleBuffer::byte *) tempParticleData);
#endif

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}
#else

bool Comet::GenerateParticlePropertiesJetNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData){
return false;
}


#endif


double Comet::radiativeCoolingRate_Crovisier(PIC::Mesh::cDataCenterNode *CenterNode){
  double tau,r,res=0.0,dens,temp;
  int idim;
  double *position;

  const double adsorptionCrossSection=4.0E-15;

  position=CenterNode->GetX();

  for (r=0.0,idim=0;idim<DIM;idim++) r+=pow(position[idim],2.0);
  r=sqrt(r)*100.0; //convert r into cm

  dens=CenterNode->GetNumberDensity(_H2O_SPEC_)*1.0e-6; //convert density into cm^{-3}
  temp=PIC::IDF::LB::GetCellRotTemp(_H2O_SPEC_,CenterNode);

  tau=0.4*dens*r*adsorptionCrossSection;

  res=(temp<52.0) ? 4.4E-22*pow(temp,3.35) : 2.0E-20*pow(temp,2.47);
  res*=dens*exp(-tau);

  return (dens>0.0) ? res*1.0E-7/dens : 0.0; //convert energy into joule and the rate into the  rate per particle
}

void Comet::StepOverTime() {
  double LocalTimeStep,Erot;
  double radiativeCoolingRate=0.0;
  
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],FirstCellParticle,ptr;
  
  int thread,i,j,k,idim,offset,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  PIC::ParticleBuffer::byte *ParticleData;

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;
    memcpy(FirstCellParticleTable,block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (i=0;i<_BLOCK_CELLS_X_;i++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)
        for (k=0;k<_BLOCK_CELLS_Z_;k++) {
	  FirstCellParticle=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

	  if (FirstCellParticle!=-1 && cell!=NULL){
	    radiativeCoolingRate=radiativeCoolingRate_Crovisier(cell);

	  if (radiativeCoolingRate>0.0) {
	    for (ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
	      if (PIC::ParticleBuffer::GetI(ptr)==_H2O_SPEC_) {
		ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);		
	  
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[PIC::ParticleBuffer::GetI(ptr)];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=block->GetLocalTimeStep(PIC::ParticleBuffer::GetI(ptr)); //   node->Sphere->maxIntersectedNodeTimeStep[PIC::ParticleBuffer::GetI(ptr)];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

		Erot=PIC::IDF::LB::GetRotE(ParticleData);
		Erot-=radiativeCoolingRate*LocalTimeStep;
		
		if(Erot<0.0) Erot=0.0;
		PIC::IDF::LB::SetRotE(Erot,ParticleData);
	      }
	      }
	    }
	  }
	}
    }
  }
#endif
}

void Comet::PrintSurfaceTriangulationMesh(const char *fname,CutCell::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,double EPS) {
  int nface,pnode,cnt;
  bool flag;
  double *xNode,*xFace;

  list<CutCell::cNodeCoordinates> nodeCoordinates;
  list<CutCell::cNodeCoordinates>::iterator nodeitr;
  CutCell::cFaceNodeConnection *FaceNodeConnection=new CutCell::cFaceNodeConnection[nSurfaceTriangulation];

  //reconstruct the node list
  for (nface=0;nface<nSurfaceTriangulation;nface++) for (pnode=0;pnode<3;pnode++) {
    flag=false;

    switch (pnode) {
    case 0:
      xFace=(SurfaceTriangulation+nface)->x0Face;
      break;
    case 1:
      xFace=(SurfaceTriangulation+nface)->x1Face;
      break;
    case 2:
      xFace=(SurfaceTriangulation+nface)->x2Face;
      break;
    }

    for (nodeitr=nodeCoordinates.begin();nodeitr!=nodeCoordinates.end();nodeitr++) {
      xNode=nodeitr->x;

      if (pow(xFace[0]-xNode[0],2)+pow(xFace[1]-xNode[1],2)+pow(xFace[2]-xNode[2],2)<EPS*EPS) {
        flag=true;
        FaceNodeConnection[nface].node[pnode]=nodeitr;
        break;
      }
    }

    if (flag==false) {
      //the node is not found -> create new node    
      CutCell::cNodeCoordinates nd;

      nd.id=0;
      nd.x=xFace;
      nodeCoordinates.push_front(nd);

      nd.pic__shadow_attribute=(productionDistributionJetNASTRAN[nface]>0)? 1:0;

      nodeCoordinates.push_front(nd);

      FaceNodeConnection[nface].node[pnode]=nodeCoordinates.begin();
    }
  }


  //print the mesh
  FILE *fout=fopen(fname,"w");
  //  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",(int)nodeCoordinates.size(),nSurfaceTriangulation);
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\",\"Jet\"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",(int)nodeCoordinates.size(),nSurfaceTriangulation);

  for (cnt=1,nodeitr=nodeCoordinates.begin();nodeitr!=nodeCoordinates.end();nodeitr++) {
    nodeitr->id=cnt++;
    //    fprintf(fout,"%e %e %e\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2]);
    fprintf(fout,"%e %e %e %i\n",nodeitr->x[0],nodeitr->x[1],nodeitr->x[2],nodeitr->pic__shadow_attribute);
  }

  for (nface=0;nface<nSurfaceTriangulation;nface++) {
    fprintf(fout,"%i %i %i\n",FaceNodeConnection[nface].node[0]->id,FaceNodeConnection[nface].node[1]->id,FaceNodeConnection[nface].node[2]->id);
  }

  fclose(fout);
  delete [] FaceNodeConnection;
}

void Comet::GetNucleusNastranInfo(cInternalNastranSurfaceData *CG) {
  Comet::CG=CG;
}


void Comet::PrintMaxLiftableSizeSurfaceTriangulationMesh(const char *fname) {
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  const double minTemp[6]={172.0,163.0,150.0,145.0,139.0,133.0};
  long int nface,nnode,pnode;

  int rank;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);
  if (rank!=0) return;

  class cTempNodeData {
  public:
    double MaxLiftableSize;
  };

  cTempNodeData *TempNodeData=new cTempNodeData[CutCell::nBoundaryTriangleNodes];

  double x[3],accl_LOCAL[3],normGravity;
  int nd,i,j,k,idim;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL;
  double GasNumberDensity,GasBulkVelocity[3],GasMass,cr2,numerator,denominator;

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
    normGravity=0.0,numerator=0.0;
    for (idim=0;idim<DIM;idim++) accl_LOCAL[idim]=0.0;
    CutCell::BoundaryTriangleFaces[nface].GetCenterPosition(x);      

#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
    //the gravity force non spherical case
    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode,false);
    Comet::GetGravityAcceleration(accl_LOCAL,nd,startNode);
    nucleusGravity::gravity(accl_LOCAL,x);
#else
    //    the gravity force spherical case
    double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    double r=sqrt(r2);
    double mass=1.0e13;
    
    for (idim=0;idim<DIM;idim++) {
      accl_LOCAL[idim]=-GravityConstant*mass/r2*x[idim]/r;
    }
#endif
    for (idim=0;idim<DIM;idim++) normGravity+=pow(accl_LOCAL[idim],2.0);
    normGravity=sqrt(normGravity);
 
    denominator=4*ElectricallyChargedDust::MeanDustDensity*normGravity;

    double norm[3],c=0.0,X=0.0,flux,surfaceTemp,v,positionSun[3];
    double HeliocentricDistance=0.0;

    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

    for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[nface].ExternalNormal[idim];
    CutCell::BoundaryTriangleFaces[nface].GetRandomPosition(x,PIC::Mesh::mesh.EPS);
    for (c=0.0,X=0.0,idim=0;idim<3;idim++){
      c+=norm[idim]*(positionSun[idim]-x[idim]);
      X+=pow(positionSun[idim]-x[idim],2.0);
    }
    if(c<0 || CutCell::BoundaryTriangleFaces[nface].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
      flux=nightSideFlux;
      surfaceTemp=minTemp[Comet::ndist];
    }else{
      double angleProd;
      int angleProdInt;
      angleProd=acos(c/sqrt(X))*180/Pi;
      angleProdInt=(int) angleProd;
      flux=fluxBjorn[angleProdInt];
      surfaceTemp=Exosphere::GetSurfaceTemeprature(c/sqrt(X),x);
    }
    GasMass=PIC::MolecularData::GetMass(0);
    v=sqrt(2*Kbol*surfaceTemp/(Pi*GasMass));
    numerator=3*GasMass*flux*v;
     
   for (pnode=0;pnode<3;pnode++) {
      nnode=CutCell::BoundaryTriangleFaces[nface].node[pnode]-CutCell::BoundaryTriangleNodes;
      if ((nnode<0)||(nnode>=CutCell::nBoundaryTriangleNodes)) exit(__LINE__,__FILE__,"Error: out of range");

      TempNodeData[nnode].MaxLiftableSize=numerator/denominator;
    }
  }

  //print the mesh
  FILE *fout=fopen(fname,"w");
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\",\"MaxLiftableSize (m)\"");
  fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",CutCell::nBoundaryTriangleNodes,CutCell::nBoundaryTriangleFaces);

  for (nnode=0;nnode<CutCell::nBoundaryTriangleNodes;nnode++) {
    fprintf(fout,"%e %e %e %e\n",CutCell::BoundaryTriangleNodes[nnode].x[0],CutCell::BoundaryTriangleNodes[nnode].x[1],CutCell::BoundaryTriangleNodes[nnode].x[2],TempNodeData[nnode].MaxLiftableSize);
  }

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
    fprintf(fout,"%ld %ld %ld\n",1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[0]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[1]-CutCell::BoundaryTriangleNodes),1+(long int)(CutCell::BoundaryTriangleFaces[nface].node[2]-CutCell::BoundaryTriangleNodes));
  }

  fclose(fout);
  delete [] TempNodeData;
#endif
}

double PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::GetTotalCrossSection(double *v0,double *v1,int s0,int s1,PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell) {
  double T=cell->GetTranslationalTemperature(_H2O_SPEC_);

  if (s0==_H2O_SPEC_ && s1==_H2O_SPEC_) return (T>1.0) ? 1.66E-19/pow(T/300.0,0.6) : 0.0;
#if _MODEL_SOURCE_DFMS_ == _MODEL_SOURCE_DFMS_ON_
  else if ((s0==_H2O_SPEC_ && s1==_CO2_SPEC_) || (s0==_CO2_SPEC_ && s1==_H2O_SPEC_)) return 3.4E-19;
  else if ((s0==_H2O_SPEC_ && s1==_CO_SPEC_) || (s0==_CO_SPEC_ && s1==_H2O_SPEC_)) return 3.2E-19;
  else if ((s0==_CO2_SPEC_ && s1==_CO_SPEC_) || (s0==_CO_SPEC_ && s1==_CO2_SPEC_)) return 3.2E-19;
  else if (s0==_CO2_SPEC_ && s1==_CO2_SPEC_) return 3.4E-19;
  else if (s0==_CO_SPEC_ && s1==_CO_SPEC_) return 3.2E-19;
#endif
else return 0.0;
}

unsigned int Comet::GetParticleSurfaceElement(PIC::ParticleBuffer::byte *ParticleDataStart) {
  return *((int*)(ParticleDataStart+offsetSurfaceElement));
}

void Comet::SetParticleSurfaceElement(int SurfaceElement,PIC::ParticleBuffer::byte *ParticleDataStart) {
      *((int*) (ParticleDataStart+offsetSurfaceElement))=SurfaceElement;
}

