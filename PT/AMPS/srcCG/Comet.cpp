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

static bool probabilityFunctionDefinedJet=false,probabilityFunctionDefined=false,probabilityFunctionDefinedWaist=false,probabilityFunctionDefinedHartley2=false,probabilityFunctionDefinedNASTRAN=false,probabilityFunctionDefinedUniformNASTRAN=false,probabilityFunctionDefinedJetNASTRAN=false;
//static double productionDistributionJet[360][180],cumulativeProductionDistributionJet[360][180];
static double productionDistributionJet[6000],cumulativeProductionDistributionJet[6000];
static double productionDistributionWaist[6000],cumulativeProductionDistributionWaist[6000];
static double productionDistributionHartley2[6000],cumulativeProductionDistributionHartley2[6000];
static double productionDistributionNASTRAN[150000],cumulativeProductionDistributionNASTRAN[150000];
static double productionDistributionUniformNASTRAN[150000],cumulativeProductionDistributionUniformNASTRAN[150000];
static double productionDistributionJetNASTRAN[150000],cumulativeProductionDistributionJetNASTRAN[150000];
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


void Comet::Init_BeforeParser() {
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
#else
  double surfaceTemperatureConstant=180.0;
  return surfaceTemperatureConstant;
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



  /*  if (probabilityFunctionDefinedNASTRAN==false) {       
    for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      TableTotalProductionRate+=ProductionRate[i][2+Comet::ndist];
    }
    
    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth);
    positionSun[2]=0.0;
    
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    
    total=0.0;      
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x,PIC::Mesh::mesh.EPS); //I had middle element on body rotation...
      
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
	c+=norm[idim]*(positionSun[idim]-x[idim]);
	X+=pow(positionSun[idim]-x[idim],2.0);
      }
      
      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow
	productionDistributionNASTRAN[i]=NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/(2*Pi*rSphere*rSphere);
	total+=productionDistributionNASTRAN[i];
      }else{
	double angleProd;
	int angleProdInt;
	angleProd=acos(c/sqrt(X))*180/Pi;
	angleProdInt=(int) angleProd;
	productionDistributionNASTRAN[i]=ProductionRate[angleProdInt][2+Comet::ndist]/(TableTotalProductionRate/(1.0-NightSideProduction[Comet::ndist]*2))/(2*Pi*rSphere*rSphere*(cos(angleProdInt*Pi/180)-cos(angleProdInt*Pi/180+Pi/180)))*CutCell::BoundaryTriangleFaces[i].SurfaceArea+NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/(2*Pi*rSphere*rSphere);
	total+=productionDistributionNASTRAN[i];
      }
      }*/

  /*  if (probabilityFunctionDefinedNASTRAN==false) {
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
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);

      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }

      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow                                                                                                                                     
	totalSurfaceInShadow+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
      }else{
        totalSurfaceInDayLight+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
      }
    }
    if (PIC::ThisThread==0) printf("totalSurfaceinShadow=%e, totalSurfaceinDayLight=%e \n",totalSurfaceInShadow,totalSurfaceInDayLight);

    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);

      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }

      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow                                                                                                                                    
	productionDistributionNASTRAN[i]=NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/totalSurfaceInShadow;
        totalProdNightSide+=productionDistributionNASTRAN[i];
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        productionDistributionNASTRAN[i]=ProductionRate[angleProdInt][2+Comet::ndist]/(TableTotalProductionRate/(1.0-NightSideProduction[Comet::ndist]*(totalSurfaceInDayLight+totalSurfaceInShadow)/totalSurfaceInDayLight))/(totalSurfaceInDayLight*(cos(angleProdInt*Pi/180)-cos(angleProdInt*Pi/180+Pi/180)))*CutCell::BoundaryTriangleFaces[i].SurfaceArea+NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/totalSurfaceInShadow;
        totalProdDaySide+=productionDistributionNASTRAN[i];
      }
    }

    scalingFactor=NightSideProduction[Comet::ndist]/(1-NightSideProduction[Comet::ndist])/(totalProdNightSide/totalProdDaySide);

    totalProdNightSide=0.0,totalProdDaySide=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);

      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
        c+=norm[idim]*(positionSun[idim]-x[idim]);
        X+=pow(positionSun[idim]-x[idim],2.0);
      }
      
      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow                                                                                                                                       
	productionDistributionNASTRAN[i]=NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/totalSurfaceInShadow*scalingFactor;
        total+=productionDistributionNASTRAN[i];
        totalProdNightSide+=productionDistributionNASTRAN[i];
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        productionDistributionNASTRAN[i]=ProductionRate[angleProdInt][2+Comet::ndist]/(TableTotalProductionRate/(1.0-NightSideProduction[Comet::ndist]*(totalSurfaceInDayLight+totalSurfaceInShadow)/totalSurfaceInDayLight))/(totalSurfaceInDayLight*(cos(angleProdInt*Pi/180)-cos(angleProdInt*Pi/180+Pi/180)))*CutCell::BoundaryTriangleFaces[i].SurfaceArea+NightSideProduction[Comet::ndist]*CutCell::BoundaryTriangleFaces[i].SurfaceArea/totalSurfaceInShadow;
        total+=productionDistributionNASTRAN[i];
        totalProdDaySide+=productionDistributionNASTRAN[i];
      }
    }

    if (PIC::ThisThread==0) printf("ratio Production Night/Total=%e \n",totalProdNightSide/(totalProdNightSide+totalProdDaySide));
  */

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

    /*    FILE *out;
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
  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]*=-1;

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

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
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
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vInit*ExternalNormal[idim];
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
  
    //    printf("total=%e \n",total);
    
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
  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]*=-1;

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

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
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
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vInit*ExternalNormal[idim];
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
  return Comet::Jet_SourceRate[spec];
}

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
  double totalArea;

  if (probabilityFunctionDefinedJetNASTRAN==false) {          
    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);
    
    totalSurfaceElementsNumber=CutCell::nBoundaryTriangleFaces;
    
    total=0.0; 
    totalArea=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      for (idim=0;idim<3;idim++) norm[idim]=CutCell::BoundaryTriangleFaces[i].ExternalNormal[idim];
      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(x);
      /*      
      if (norm[0]<0){
      //      double azimuthAngle=25*Pi/180,zenithAngle=15.0*Pi/180,azimuthJet=135*Pi/180,zenithJet=0.0*Pi/180;
      //      double azimuthAngle=50*Pi/180,zenithAngle=20.0*Pi/180,azimuthJet=340*Pi/180,zenithJet=90.0*Pi/180;
      double azimuthAngle=70*Pi/180,zenithAngle=2.0*Pi/180,azimuthJet=340*Pi/180,zenithJet=50.0*Pi/180;
      double rxy=sqrt(x[1]*x[1]+x[0]*x[0]);
      double azimuthCenter;
      double zenithCenter=acos(x[2]/sqrt(x[2]*x[2]+rxy*rxy));

      if (asin(x[1]/rxy)>=0) azimuthCenter=acos(x[0]/rxy);
      else azimuthCenter=-acos(x[0]/rxy)+2*Pi;

      if (zenithJet-zenithAngle<0) {
	if((azimuthCenter>=(azimuthJet-azimuthAngle) && azimuthCenter<=(azimuthJet+azimuthAngle) && zenithCenter>=(zenithJet-zenithAngle) && zenithCenter<=(zenithJet+zenithAngle)) || (azimuthCenter>=(azimuthJet-Pi-azimuthAngle) && azimuthCenter<=(azimuthJet-Pi+azimuthAngle) && zenithCenter<=(zenithAngle-zenithJet) && 0.0<=(zenithAngle-zenithJet))){
	  productionDistributionJetNASTRAN[i]=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
	  total+=productionDistributionJetNASTRAN[i];
	}
      }
      else if (azimuthJet+azimuthAngle>=2*Pi){
	if((azimuthCenter>=(azimuthJet-azimuthAngle) && azimuthCenter<=(azimuthJet+azimuthAngle) && zenithCenter>=(zenithJet-zenithAngle) && zenithCenter<=(zenithJet+zenithAngle)) || (azimuthCenter<=(azimuthJet-2*Pi+azimuthAngle) && 0.0<=(azimuthJet-2*Pi+azimuthAngle) && zenithCenter>=(zenithJet-zenithAngle) && zenithCenter<=(zenithJet+zenithAngle))){
	  productionDistributionJetNASTRAN[i]=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
	  total+=productionDistributionJetNASTRAN[i];
	}
      }
      else{
	if(azimuthCenter>=(azimuthJet-azimuthAngle) && azimuthCenter<=(azimuthJet+azimuthAngle) && zenithCenter>=(zenithJet-zenithAngle) && zenithCenter<=(zenithJet+zenithAngle)){
	  productionDistributionJetNASTRAN[i]=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
	  total+=productionDistributionJetNASTRAN[i];
	}
      }
      }*/
      //      if(x[0]>-24 && x[0]<840 && x[1]>-700 && x[1]<652 && x[2]>398 && x[2]<607) {
      if(x[0]>-120 && x[0]<840 && x[1]>-1000 && x[1]<652 && x[2]>0 && x[2]<607) {
	productionDistributionJetNASTRAN[i]=CutCell::BoundaryTriangleFaces[i].SurfaceArea;
	  total+=productionDistributionJetNASTRAN[i];
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
    probabilityFunctionDefinedJetNASTRAN=true;
  }
  
  //Computation of the segment where the particle will be created
  gamma=rnd();
  i=0;
  while (gamma>cumulativeProductionDistributionJetNASTRAN[i]){
    i++;
  }
    
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
  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]*=-1;

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

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=-vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
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
  for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-vInit*ExternalNormal[idim];
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
		char ParticleData[PIC::ParticleBuffer::ParticleDataLength];
		
		memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);
	  
#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[PIC::ParticleBuffer::GetI(ptr)];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=block->GetLocalTimeStep(PIC::ParticleBuffer::GetI(ptr)); //   node->Sphere->maxIntersectedNodeTimeStep[PIC::ParticleBuffer::GetI(ptr)];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

		Erot=PIC::IDF::LB::GetRotE((PIC::ParticleBuffer::byte*)ParticleData);
		Erot-=radiativeCoolingRate*LocalTimeStep;
		
		if(Erot<0.0) Erot=0.0;
		PIC::IDF::LB::SetRotE(Erot,(PIC::ParticleBuffer::byte*)ParticleData);
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
