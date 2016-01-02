/*
 * Comet.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and vtenishe
 */

//$Id$

#include "pic.h"
#include "Dust.h"
#include "Comet.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="CHURYUMOV-GERASIMENKO";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="67P/C-G_CK";

int Comet::GravityFieldOffset=-1;

char Comet::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";

//parameters of the initial dust velocity distribution
int Comet::DustInitialVelocity::VelocityModelMode=Comet::DustInitialVelocity::Mode::ConstantVelocity;
double Comet::DustInitialVelocity::RotationPeriod=1.0E10;
double Comet::DustInitialVelocity::RotationAxis[3]={1.0,0.0,0.0};
double Comet::DustInitialVelocity::InjectionConstantVelocity=0.001;

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

static double productionDistributionJetNASTRAN[200000],cumulativeProductionDistributionJetNASTRAN[200000];
static bool probabilityFunctionDefinedJetNASTRAN=false;

static double productionDistribution[180],cumulativeProductionDistribution[180];
static double angle;
static double azimuthCenter;
static double zenithCenter;
static cInternalRotationBodyData* Nucleus;

double subSolarPointAzimuth=0.0; //53.0*Pi/180;
double subSolarPointZenith=0.0;

double DustSizeMin=1.0e-7;
double DustSizeMax=1.0e-2;
double DustTotalMassProductionRate=0.0;
int DustSampleIntervals=10;
double DustSizeDistribution=0.0;



#if _MULTISPECIES_ANALYTICAL_MODE_ == _MULTISPECIES_ANALYTICAL_MODE_ON_
static double ratioBjornSpec[13][200000];
static double productionDistributionNASTRAN[13][200000],cumulativeProductionDistributionNASTRAN[13][200000],fluxBjornANALYTICAL[13][200000];
static bool definedFluxBjorn[13],probabilityFunctionDefinedNASTRAN[13];
static double BjornProductionANALYTICAL[13];
double fluxBjorn[13][90];
double nightSideFlux[13];
#else
static double productionDistributionNASTRAN[200000],cumulativeProductionDistributionNASTRAN[200000];
static bool probabilityFunctionDefinedNASTRAN=false;
static bool definedFluxBjorn=false;
static double BjornProduction;
double fluxBjorn[90];
double nightSideFlux;
#endif

static double gravityAngle[200000];
static bool gravityAngleInitialized=false;

long int offsetSurfaceElement;

void Comet::Init_BeforeParser() {
  Exosphere::Init_SPICE();

#if _TRACKING_SURFACE_ELEMENT_MODE_ == _TRACKING_SURFACE_ELEMENT_MODE_ON_
  // Keep track of original surface element the particle was created from
  PIC::ParticleBuffer::RequestDataStorage(offsetSurfaceElement,sizeof(int));
#endif


#if _3DGRAVITY__MODE_ == _3DGRAVITY__MODE__ON_
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

  //request sampling data for calculation of the calculating of the dust brightness
  PIC::IndividualModelSampling::RequestSamplingData.push_back(Sampling::RequestSamplingData);

  //set procedure that prints the column integrals
  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SampleModelData,Sampling::PrintBrightnessMap);

  //init sampling of the dust size distribution
  Comet::Sampling::InjectedDustSizeDistribution::Init();
#endif


}

void Comet::Init_AfterParser(const char *DataFilePath) {

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
#if _3DGRAVITY__MODE_ == _3DGRAVITY__MODE__ON_
  const int nGravityVariables=3;
  int GravityVariableOffsets[nGravityVariables]={GravityFieldOffset,GravityFieldOffset+sizeof(double),GravityFieldOffset+2*sizeof(double)};

  if (PIC::Mesh::mesh.AssociatedDataFileExists("gravity",DataFilePath)==true) {
    PIC::Mesh::mesh.LoadCenterNodeAssociatedData("gravity",DataFilePath,GravityVariableOffsets,nGravityVariables);
  }
  else {
    InitGravityData();
    PIC::Mesh::mesh.SaveCenterNodeAssociatedData("gravity",GravityVariableOffsets,nGravityVariables);
  }
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
  long int ct=0;

  //evaluate the total nucleus mass
  double rTest[3]={10.0E3,0.0,0.0};
  double TotalNucleusMass=nucleusGravity::gravity(gravityAccl,rTest);
  double r;

  const double rmin=1.5*4.0E3;
  const double rmax=2.0*4.0E3;

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) {
            position=cell->GetX();
            r=sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);

            if (r<rmin) {
              nucleusGravity::gravity(gravityAccl,position);
            }
            else if (r<rmax) {
              nucleusGravity::gravity(gravityAccl,position);

              //get the interpolation weight
              double t,c;
              c=1.0-(r-rmin)/(rmax-rmin);

              //recalculate the gravity aceleration accounting for the exect gravity of the nucleus and gravity of a sphere
              t=GravityConstant*TotalNucleusMass/pow(r,3);

              for (idim=0;idim<3;idim++) gravityAccl[idim]=c*gravityAccl[idim]-(1.0-c)*t*position[idim];
            }
            else {
              //get the gravity of a sphere
              double t=GravityConstant*TotalNucleusMass/pow(r,3);

              for (idim=0;idim<3;idim++) gravityAccl[idim]=-t*position[idim];
            }


            for (idim=0;idim<3;idim++) {
              *((double*)(cell->GetAssociatedDataBufferPointer()+GravityFieldOffset+idim*sizeof(double)))=gravityAccl[idim];
              gravityAccl[idim]=0.0;
            }
          }
	     }
    }

    ct+=1;
    printf("PIC::Mesh::mesh.ThisThread=%i ct=%li \n",PIC::Mesh::mesh.ThisThread,ct);
  }

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
  int SourceProcessID,iInjectionFaceNASTRAN;

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
       flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
       GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
       GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
       SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
       if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
     }
     else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_) {
       flag=Comet::GenerateParticlePropertiesUniformNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
       GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
       GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
       SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_;
       if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
     }
     else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_) {
       flag=Comet::GenerateParticlePropertiesJetNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
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

  //apply condition of tracking the particle
  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(tempParticleData);
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x_SO_OBJECT,v_SO_OBJECT,spec,tempParticleData);
  #endif

  #if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  ElectricallyChargedDust::SetGrainCharge(0.0,(PIC::ParticleBuffer::byte*)tempParticleData);
  #endif

  ElectricallyChargedDust::SetGrainMass(GrainMass,(PIC::ParticleBuffer::byte*)tempParticleData);
  ElectricallyChargedDust::SetGrainRadius(GrainRadius,(PIC::ParticleBuffer::byte*)tempParticleData);

  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);


  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  //determine the initial charge of the dust grain
  #if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
  ElectricallyChargedDust::DustChargingProcessor_SteadyState(x_SO_OBJECT,x_SO_OBJECT,v_SO_OBJECT,spec,newParticle,newParticleData,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode);
  #endif

  nInjectedParticles++;

  //sample the injection rate
  CutCell::BoundaryTriangleFaces[iInjectionFaceNASTRAN].UserData.InjectionFlux[spec]+=
      GrainWeightCorrection*PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup]/
      startNode->block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/CutCell::BoundaryTriangleFaces[iInjectionFaceNASTRAN].SurfaceArea;

  //sample the size distribution
  Comet::Sampling::InjectedDustSizeDistribution::AddParticleData(GrainRadius,GrainWeightCorrection*PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup]);


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
    flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
    SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
    if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  }

  else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_) {
    flag=Comet::GenerateParticlePropertiesUniformNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
    SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_;
    if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__1_Uniform_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
  }  

  else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__2_Jet_) {
    flag=Comet::GenerateParticlePropertiesJetNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData,iInjectionFaceNASTRAN);
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

  //apply condition of tracking the particle
  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::InitParticleID(tempParticleData);
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x_SO_OBJECT,v_SO_OBJECT,spec,tempParticleData);
  #endif

  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);


  newParticle=PIC::ParticleBuffer::GetNewParticle();
  newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
  memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

  nInjectedParticles++;

  //sample the injection rate
  CutCell::BoundaryTriangleFaces[iInjectionFaceNASTRAN].UserData.InjectionFlux[spec]+=
      ParticleWeightCorrection*ParticleWeightCorrection/
      startNode->block->GetLocalTimeStep(spec)/CutCell::BoundaryTriangleFaces[iInjectionFaceNASTRAN].SurfaceArea;


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

double sphericalHarmonic(int i,double colatitude,double longitude){                                                                                                            
  double res;                                                                                                                                                              
                                                                                                                                                              
  switch (i) {                                                                                                                                                              
  case 0:                                                                                                                                                                     
    res=1; //Y00                                                                                                                                                              
    break;                                                                                                                                                                     
  case 1:                                                                                                                                                                      
    res=sin(colatitude)*sin(longitude);//Y1-1                                                                                                                                   
    break;                                                                                                                                                                     
  case 2:                                                                                                                                                                       
    res=cos(colatitude);                                                                                                                                                       
    break;                                                                                                                                                                     
  case 3:                                                                                                                                                                      
    res=sin(colatitude)*cos(longitude);//Y11                                                                                                                                    
    break;                                                                                                                                                                     
  case 4:                                                                                                                                                                       
    res=pow(sin(colatitude),2.0)*sin(2*longitude);//Y2-2                                                                                                                       
    break;                                                                                                                                                                      
  case 5:                                                                                                                                                                       
    res=cos(colatitude)*sin(colatitude)*sin(longitude);//Y2-1                                                                                                                   
    break;                                                                                                                                                                      
  case 6:                                                                                                                                                                      
    res=3*pow(cos(colatitude),2.0)-1.0;//Y20                                                                                                                                    
    break;                                                                                                                                                                      
  case 7:                                                                                                                                                                     
    res=cos(colatitude)*sin(colatitude)*cos(longitude);//Y21                                                                                                                    
    break;                                                                                                                                                                     
  case 8:                                                                                                                                                                       
    res=pow(sin(colatitude),2.0)*cos(2*longitude);//Y22                                                                                                                         
    break;                                                                                                                                                                     
  case 9:                                                                                                                                                                       
    res=pow(sin(colatitude),3.0)*sin(3*longitude);//Y3-3                                                                                                                       
    break;                                                                                                                                                                      
  case 10:                                                                                                                                                                      
    res=pow(sin(colatitude),2.0)*cos(colatitude)*sin(2*longitude);//Y3-2                                                                                                        
    break;                                                                                                                                                                      
  case 11:                                                                                                                                                                     
    res=5*(pow(cos(colatitude),2.0)-1)*sin(colatitude)*sin(longitude);//Y3-1                                                                                                   
    break;                                                                                                                                                                      
  case 12:                                                                                                                                                                      
    res=5*pow(cos(colatitude),3.0)-3*cos(colatitude);//Y30                                                                                                                     
    break;                                                                                                                                                                      
  case 13:                                                                                                                                                                      
    res=5*(pow(cos(colatitude),2.0)-1)*sin(colatitude)*cos(longitude);//Y31                                                                                                     
    break;                                                                                                                                                                      
  case 14:                                                                                                                                                                      
    res=pow(sin(colatitude),2.0)*cos(colatitude)*cos(2*longitude);//Y32                                                                                                         
    break;
  case 15:
    res=pow(sin(colatitude),3.0)*cos(3*longitude);//Y33                                                                                                                         
    break;
  case 16:
    res=pow(sin(colatitude),4.0)*sin(4*longitude);//Y4-4                                                                                                                        
    break;
  case 17:
    res=pow(sin(colatitude),3.0)*cos(colatitude)*sin(3*longitude);//Y4-3                                                                                                        
    break;
  case 18:
    res=pow(sin(colatitude),2.0)*(7*pow(cos(colatitude),2.0)-1)*sin(2*longitude);//Y4-2                                                                                         
    break;
  case 19:
    res=sin(colatitude)*(7*pow(cos(colatitude),3.0)-3*cos(colatitude))*sin(longitude);//Y4-1                                                                                    
    break;
  case 20:
    res=35*pow(cos(colatitude),4.0)-30*pow(cos(colatitude),2.0)+3;//Y40                                                                                                         
    break;
  case 21:
    res=sin(colatitude)*(7*pow(cos(colatitude),3.0)-3*cos(colatitude))*cos(longitude);//Y41                                                                                     
    break;
  case 22:
    res=pow(sin(colatitude),2.0)*(7*pow(cos(colatitude),2.0)-1)*cos(2*longitude);//Y42                                                                                          
    break;
  case 23:
    res=pow(sin(colatitude),3.0)*cos(colatitude)*cos(3*longitude);//Y43                                                                                                         
    break;
  case 24:
    res=pow(sin(colatitude),4.0)*cos(4*longitude);//Y44                                                                                                                         
    break;
  default:
    exit(__LINE__,__FILE__,"Error: order of spherical Harmonic too high, not implemented yet");
  }
  return res;
}


double Comet::GetTotalProductionRateBjornNASTRAN(int spec){
  double rSphere=1980,c=0.0,X=0.0,totalProductionRate=0.0;
  double positionSun[3],x[3],norm[3],xMiddle[3],lattitude,factor=1.0;
  long int totalSurfaceElementsNumber,i,j;
  double percentageActive=0.05;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  const double DistanceFromTheSun[6]={1.3,2.0,2.7,3.0,3.25,3.5};
  double HeliocentricDistance=3.3*_AU_;
  int idim;
  double ProductionRateScaleFactor,TableTotalProductionRate=0.0,BjornTotalProductionRate;
  int angle;
  double angletemp;

  int nDev=25;

  //Nucleus activity                                                                                                                                                        
  const double  Activity[3][25]={
    {1.88802089096550e+18,-3.32193190779201e+18,2.16030831636854e+18,1.16303584760745e+18,-3.48031365453629e+17,-3.97108996341047e+18,2.32187315012071e+18,2.62881801954068e+18,-1.64152743317681e+17,5.48474318492987e+16,-8.81665110610612e+16,-6.71346849527855e+17,8.17079244731431e+17,2.10263858732877e+17,-7.31447243364394e+17,1.87954830493877e+16,1.59517599584823e+16,2.22312552878431e+17,-4.12879355040244e+17,-1.37905625912140e+18,1.83112475092734e+17,1.21579175185910e+18,-2.43316081589516e+17,-4.24836863227363e+17,2.11834459021013e+16},
    {1.33147318596808e+16,-5.99325576797965e+15,-1.44574576415506e+16,-1.23844936447417e+16,-1.55154864153184e+15,-6.53313342291062e+15,1.07049479617418e+16,1.24456131751260e+16,-6.54238886353421e+15,1.12926642418814e+15,3.89604594220916e+15,-531055729734858,-398604759758765,-2.61684944191026e+15,1.41771647341777e+16,2.03706829667621e+15,-351642267595628,-1.40564295976192e+15,-2.04618374895345e+15,-6.09023703216270e+15,349833485542175,3.58729877013097e+15,-4.35629505817132e+15,-2.91104899991473e+15,1.36495458239451e+15},
    {8.24876290734347e+15,-1.15993586348543e+16,3.36505486424125e+15,-6.76013519095671e+15,-314999862632954,-1.08780416335274e+16,7.95233182311777e+15,9.16964842516085e+15,-2.81955448931900e+15,1.21059245593790e+15,-1.25443670217006e+15,-2.11455950796835e+15,1.24045282758517e+15,-1.65067535925255e+15,-5.46839069247522e+15,1.09833316361053e+15,264156854265098,1.90947201360750e+15,-892524030311892,-2.10255875207271e+15,515450866463768,3.93817676318131e+15,-2.90479115840660e+15,-5.21185256041148e+15,955141456973813}
  };

#if _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ ==  _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ANALYTICAL_

#if _MULTISPECIES_ANALYTICAL_MODE_ == _MULTISPECIES_ANALYTICAL_MODE_ON_
  if (definedFluxBjorn[spec]==false) {
    if (spec==_H2O_SPEC_) {
      double Qmin=0.02/pow(HeliocentricDistance/_AU_,4.2143229)*600,Qmax=1.0/pow(HeliocentricDistance/_AU_,4.2143229)*600;
      
      for (i=0;i<90;i++) {
        angle=(double) i;
        fluxBjorn[spec][i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
      }
      
      nightSideFlux[spec]=Qmin;
    }
    else if (spec==_CO2_SPEC_) {
      double Qmin=0.1/pow(HeliocentricDistance/_AU_,2.0)*600,Qmax=1.0/pow(HeliocentricDistance/_AU_,2.0)*600;
      for (i=0;i<90;i++) {
	angle=(double) i;
	fluxBjorn[spec][i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
      }

      nightSideFlux[spec]=Qmin;
    }
    /*    else if (spec==_CO_SPEC_) {
      double Qmin=0.1/pow(HeliocentricDistance/_AU_,2.0)*600,Qmax=1.0/pow(HeliocentricDistance/_AU_,2.0)*600;
      for (i=0;i<90;i++) {
	angle=(double) i;
	fluxBjorn[spec][i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
	}

      nightSideFlux[spec]=Qmin;
      }*/

      /*else  if (spec==_O2_SPEC_) {
	double Qmin=0.02/pow(HeliocentricDistance/_AU_,4.2143229)*60,Qmax=1.0/pow(HeliocentricDistance/_AU_,4.2143229)*60;
	
	for (i=0;i<90;i++) {
	  angle=(double) i;
	  fluxBjorn[spec][i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
	}
	
	nightSideFlux[spec]=Qmin;
	}*/
  else{
      double Qmin=0.0,Qmax=0.0;
      for (i=0;i<90;i++) {
	angle=(double) i;
	fluxBjorn[spec][i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
      }
      
      nightSideFlux[spec]=Qmin;  
    }
  }
#else //ONLY ONE SPECIES
  if (definedFluxBjorn==false) {
    double Qmin=5.0e17,Qmax=7.0e18;

     for (i=0;i<90;i++) {
       angle=(double) i;
       fluxBjorn[i]=Qmin+(Qmax-Qmin)*cos(angle*Pi/180.0);
     }
     
     nightSideFlux=Qmin;  
  }
#endif

#else
  if (definedFluxBjorn==false) {
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
  }
#endif

  

#if _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ ==  _BJORN_PRODUCTION_RATE_USERDEFINED_MODE_ON_
  return Comet::Bjorn_SourceRate[spec];

#else


#if _MULTISPECIES_ANALYTICAL_MODE_ == _MULTISPECIES_ANALYTICAL_MODE_ON_
  if (definedFluxBjorn[spec]==false) {
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

      double colatitude,longitude;
      double xCenter[3];

      double factor=0.0;


      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(xCenter);

      colatitude=acos(xCenter[2]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]+xCenter[2]*xCenter[2]));
      if (xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]==0.0) {
	longitude=0.0;
      }else if(xCenter[1]>0.0) {
        longitude=acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
      }else{
        longitude=-acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
      }


      if (spec==_H2O_SPEC_) {
	for (j=0;j<nDev;j++) factor+=Activity[0][j]*sphericalHarmonic(j,colatitude,longitude);
      }
      else if (spec==_CO2_SPEC_) {
        for (j=0;j<nDev;j++) factor+=Activity[1][j]*sphericalHarmonic(j,colatitude,longitude);
      }
      /*      else if (spec==_CO_SPEC_) {                                                                              
        for (j=0;j<nDev;j++) factor+=Activity[2][j]*sphericalHarmonic(j,colatitude,longitude);                         
        }*/
      else {
	factor=0.0;
      }

      if (factor<0.0) factor=0.0; //take care of the few cases where the constraints of positive activity is violated  

      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) { //Test Shadow
	totalProductionRate+=nightSideFlux[spec]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*factor;
      }else{
	double angleProd;
	int angleProdInt;
	angleProd=acos(c/sqrt(X))*180/Pi;
	angleProdInt=(int) angleProd;
	totalProductionRate+=fluxBjorn[spec][angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*factor;
      }
    }
    if (probabilityFunctionDefinedNASTRAN[spec]==false && PIC::ThisThread==0)  printf("spec=%i totalProductionRate=%e flux(0)=%e nightsideFlux=%e \n",spec,totalProductionRate,fluxBjorn[spec][0],nightSideFlux[spec]);
    BjornProductionANALYTICAL[spec]=totalProductionRate;
    definedFluxBjorn[spec]=true;
  }
  return BjornProductionANALYTICAL[spec];
#else
  if (definedFluxBjorn==false) {
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
    if (probabilityFunctionDefinedNASTRAN==false && PIC::ThisThread==0)  printf("spec=%i totalProductionRate=%e flux(0)=%e nightsideFlux=%e \n",spec,totalProductionRate,fluxBjorn[0],nightSideFlux);
    BjornProduction=totalProductionRate;
    definedFluxBjorn=true;
  } 
  
  return BjornProduction;

#endif

#endif

}

#if _MULTISPECIES_ANALYTICAL_MODE_ == _MULTISPECIES_ANALYTICAL_MODE_ON_

bool Comet::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData,int &iInjectionFaceNASTRAN) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[6]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0,12.7/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3],xMiddle[3],factor=1.0,lattitude;
  static double positionSun[3];
  double HeliocentricDistance=3.5*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i,j;
  double rSphere=1980.0;
  double area;
  double totalProdNightSide=0.0,totalProdDaySide=0.0,scalingFactor,scalingFactorDay,totalSurfaceInShadow=0.0,totalSurfaceInDayLight=0.0;

  int nDev=25;

  const double  Activity[3][25]={
    {1.88802089096550e+18,-3.32193190779201e+18,2.16030831636854e+18,1.16303584760745e+18,-3.48031365453629e+17,-3.97108996341047e+18,2.32187315012071e+18,2.62881801954068e+18,-1.64152743317681e+17,5.48474318492987e+16,-8.81665110610612e+16,-6.71346849527855e+17,8.17079244731431e+17,2.10263858732877e+17,-7.31447243364394e+17,1.87954830493877e+16,1.59517599584823e+16,2.22312552878431e+17,-4.12879355040244e+17,-1.37905625912140e+18,1.83112475092734e+17,1.21579175185910e+18,-2.43316081589516e+17,-4.24836863227363e+17,2.11834459021013e+16},
    {1.33147318596808e+16,-5.99325576797965e+15,-1.44574576415506e+16,-1.23844936447417e+16,-1.55154864153184e+15,-6.53313342291062e+15,1.07049479617418e+16,1.24456131751260e+16,-6.54238886353421e+15,1.12926642418814e+15,3.89604594220916e+15,-531055729734858,-398604759758765,-2.61684944191026e+15,1.41771647341777e+16,2.03706829667621e+15,-351642267595628,-1.40564295976192e+15,-2.04618374895345e+15,-6.09023703216270e+15,349833485542175,3.58729877013097e+15,-4.35629505817132e+15,-2.91104899991473e+15,1.36495458239451e+15},
    {8.24876290734347e+15,-1.15993586348543e+16,3.36505486424125e+15,-6.76013519095671e+15,-314999862632954,-1.08780416335274e+16,7.95233182311777e+15,9.16964842516085e+15,-2.81955448931900e+15,1.21059245593790e+15,-1.25443670217006e+15,-2.11455950796835e+15,1.24045282758517e+15,-1.65067535925255e+15,-5.46839069247522e+15,1.09833316361053e+15,264156854265098,1.90947201360750e+15,-892524030311892,-2.10255875207271e+15,515450866463768,3.93817676318131e+15,-2.90479115840660e+15,-5.21185256041148e+15,955141456973813}
  };



  if (probabilityFunctionDefinedNASTRAN[spec]==false) {

    if (spec>=13) exit(__LINE__,__FILE__,"Error: the max number of the species of 13 is hardwired into the comet model. Fix it!");
    if (CutCell::nBoundaryTriangleFaces>=200000) exit(__LINE__,__FILE__,"Error: the max number of the surface elements of 200000 is hardwired into the comet model. Fix it!");

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

      double colatitude,longitude;
      double xCenter[3];

      double factor=0.0;


      CutCell::BoundaryTriangleFaces[i].GetCenterPosition(xCenter);

      colatitude=acos(xCenter[2]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]+xCenter[2]*xCenter[2]));
      if (xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]==0.0) {
        longitude=0.0;
      }else if(xCenter[1]>0.0) {
        longitude=acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
      }else{
        longitude=-acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
      }


      if (spec==_H2O_SPEC_) {
        for (j=0;j<nDev;j++) factor+=Activity[0][j]*sphericalHarmonic(j,colatitude,longitude);
      }
      else if (spec==_CO2_SPEC_) {
        for (j=0;j<nDev;j++) factor+=Activity[1][j]*sphericalHarmonic(j,colatitude,longitude);
      }
      /*      else if (spec==_CO_SPEC_) {                                                                                                                              
	for (j=0;j<nDev;j++) factor+=Activity[2][j]*sphericalHarmonic(j,colatitude,longitude);                                                                
	}*/
      else {
        factor=0.0;
      }

      if (factor<0.0) factor=0.0; //take care of the few cases where the constraints of positive activity is violated                                                  


      if(c<0 || CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
        productionDistributionNASTRAN[spec][i]=nightSideFlux[spec]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*factor;
        total+=productionDistributionNASTRAN[spec][i];
      }else{
        double angleProd;
        int angleProdInt;
        angleProd=acos(c/sqrt(X))*180/Pi;
        angleProdInt=(int) angleProd;
        productionDistributionNASTRAN[spec][i]=fluxBjorn[spec][angleProdInt]*CutCell::BoundaryTriangleFaces[i].SurfaceArea*factor;
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

    //the cumulative distributino of dust has to be the same as of water
    if (_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_) {
      if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
        //water must be defined before dust
        if (probabilityFunctionDefinedNASTRAN[_H2O_SPEC_]==false) exit(__LINE__,__FILE__,"Error: water cumulative distribution has to be defined before that of the dust");

        //copy water cumulative distribution into that of the dust
        for (i=0;i<totalSurfaceElementsNumber;i++) cumulativeProductionDistributionNASTRAN[spec][i]=cumulativeProductionDistributionNASTRAN[_H2O_SPEC_][i];
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

  iInjectionFaceNASTRAN=i; //output the triangulation surface element number

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
//  double r2Tang=0.0;;
//  double xFace[3];
//  double vDustInit=0.0001;
//  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  /*for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);*/

    DustInitialVelocity::GetInitialVelocity(v_LOCAL_IAU_OBJECT,x_LOCAL_IAU_OBJECT,i);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_ 
  for (idim=0;idim<3;idim++) ExternalNormal[idim]=-ExternalNormal[idim];     
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
bool Comet::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData, int &iInjectionFaceNASTRAN) {
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

  iInjectionFaceNASTRAN=i; //output the surface element of the particle generation

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
//  double r2Tang=0.0;;
//  double xFace[3];
//  double vDustInit=0.0001;
//  double angleVelocityNormal=asin(rnd());


  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
/*
    for (idim=0;idim<3;idim++){
      v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
    }

    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);*/

    DustInitialVelocity::GetInitialVelocity(v_LOCAL_IAU_OBJECT,x_LOCAL_IAU_OBJECT,i);

    if (ExternalNormal[0]*v_LOCAL_IAU_OBJECT[0]+ ExternalNormal[1]*v_LOCAL_IAU_OBJECT[1]+ ExternalNormal[2]*v_LOCAL_IAU_OBJECT[2]<0.0) {
      exit(__LINE__,__FILE__,"Error: the initial velocity of the grains is directed inside the nucleus");
    }
  }
  else {
    for (idim=0;idim<3;idim++) ExternalNormal[idim]=-ExternalNormal[idim];

    PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  }
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]=-ExternalNormal[idim];    
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

bool Comet::GenerateParticlePropertiesUniformNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData,int &iInjectionFaceNASTRAN) {
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

  iInjectionFaceNASTRAN=i; //output the triangulation surface element number

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
//  double r2Tang=0.0;
//  double xFace[3];
//  double vDustInit=0.0001;
//  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
/*  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
    CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
    for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
    for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);*/

    DustInitialVelocity::GetInitialVelocity(v_LOCAL_IAU_OBJECT,x_LOCAL_IAU_OBJECT,i);
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else

  #if _PIC_MODEL__RADIAL_VELOCITY_MODE_ == _PIC_MODEL__RADIAL_VELOCITY_MODE__OFF_  
  for (idim=0;idim<3;idim++) ExternalNormal[idim]=-ExternalNormal[idim];    
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
  return Comet::Jet_SourceRate[spec];

}



bool Comet::GenerateParticlePropertiesJetNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData,int &iInjectionFaceNASTRAN){
return false;
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
/*
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

#if _3DGRAVITY__MODE_ == _3DGRAVITY__MODE__ON_
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
*/
}

double PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::GetTotalCrossSection(double *v0,double *v1,int s0,int s1,PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell) {
  double T=cell->GetTranslationalTemperature(_H2O_SPEC_);

  if (s0==_H2O_SPEC_ && s1==_H2O_SPEC_) return (T>1.0) ? 1.66E-19/pow(T/300.0,0.6) : 0.0;
#if  _MULTISPECIES_ANALYTICAL_MODE_ == _MULTISPECIES_ANALYTICAL_MODE_ON_
  else if ((s0==_H2O_SPEC_ && s1==_CO2_SPEC_) || (s0==_CO2_SPEC_ && s1==_H2O_SPEC_)) return 3.4E-19;
  else if (s0==_CO2_SPEC_ && s1==_CO2_SPEC_) return 3.4E-19;
  /*  else if ((s0==_H2O_SPEC_ && s1==_CO_SPEC_) || (s0==_CO_SPEC_ && s1==_H2O_SPEC_)) return 3.2E-19;
  else if ((s0==_CO2_SPEC_ && s1==_CO_SPEC_) || (s0==_CO_SPEC_ && s1==_CO2_SPEC_)) return 3.2E-19;
    else if (s0==_CO_SPEC_ && s1==_CO_SPEC_) return 3.2E-19;
  */
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  /* //else if (s0==_OH_SPEC_ && s1==_CO_SPEC_) return 3.0E-19;
  //else if (s0==_H2_SPEC_ && s1==_CO_SPEC_) return 3.0E-19;
  //else if (s0==_H_SPEC_ && s1==_CO_SPEC_) return 1.5E-19;
  //else if (s0==_O_SPEC_ && s1==_CO_SPEC_) return 1.5E-19;

  else if (s0==_H2O_SPEC_ && s1==_OH_SPEC_) return 3.2E-19;
  else if ((s0==_H2O_SPEC_ && s1==_H2_SPEC_) || (s0==_H2_SPEC_ && s1==_H2O_SPEC_)) return 3.2E-19;
  else if ((s0==_H2O_SPEC_ && s1==_H_SPEC_) || (s0==_H_SPEC_ && s1==_H2O_SPEC_)) return 1.8E-19;
  else if ((s0==_H2O_SPEC_ && s1==_O_SPEC_) || (s0==_O_SPEC_ && s1==_H2O_SPEC_)) return 1.8E-19;
  else if (s0==_OH_SPEC_ && s1==_OH_SPEC_) return 3.0E-19;
  else if ((s0==_OH_SPEC_ && s1==_H2_SPEC_) || (s0==_H2_SPEC_ && s1==_OH_SPEC_)) return 3.0E-19;
  else if ((s0==_OH_SPEC_ && s1==_H_SPEC_) || (s0==_H_SPEC_ && s1==_OH_SPEC_)) return 1.5E-19;
  else if ((s0==_OH_SPEC_ && s1==_O_SPEC_) || (s0==_O_SPEC_ && s1==_OH_SPEC_)) return 1.5E-19;
  else if (s0==_H2_SPEC_ && s1==_H2_SPEC_) return 3.0E-19;
  else if ((s0==_H2_SPEC_ && s1==_H_SPEC_) || (s0==_H_SPEC_ && s1==_H2_SPEC_)) return 1.5E-19;
  else if ((s0==_H2_SPEC_ && s1==_O_SPEC_) || (s0==_O_SPEC_ && s1==_H2_SPEC_)) return 1.5E-19;
  else if (s0==_H_SPEC_ && s1==_H_SPEC_) return 1.2E-19;
  else if ((s0==_H_SPEC_ && s1==_O_SPEC_) || (s0==_O_SPEC_ && s1==_H_SPEC_)) return 1.2E-19;
  else if (s0==_O_SPEC_ && s1==_O_SPEC_) return 1.2E-19;

  else if ((s0==_H2O_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_H2O_SPEC_)) return 3.2E-19;
  else if ((s0==_CO2_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_CO2_SPEC_)) return 3.2E-19;
  else if (s0==_O2_SPEC_ && s1==_O2_SPEC_) return 3.2E-19;
  else if ((s0==_OH_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_OH_SPEC_)) return 3.0E-19;
  else if ((s0==_H2_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_H2_SPEC_)) return 3.0E-19;
  else if ((s0==_H_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_H_SPEC_)) return 1.5E-19;
  else if ((s0==_O_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_O_SPEC_)) return 1.5E-19;
  //  else if ((s0==_CO_SPEC_ && s1==_O2_SPEC_) || (s0==_O2_SPEC_ && s1==_CO_SPEC_)) return 3.2E-19;
  */
#endif

#endif
else return 0.0;
}

unsigned int Comet::GetParticleSurfaceElement(PIC::ParticleBuffer::byte *ParticleDataStart) {
  return *((int*)(ParticleDataStart+offsetSurfaceElement));
}

void Comet::SetParticleSurfaceElement(int SurfaceElement,PIC::ParticleBuffer::byte *ParticleDataStart) {
      *((int*) (ParticleDataStart+offsetSurfaceElement))=SurfaceElement;
}

//tracing trajectories of the individual dust grains
cInternalSphericalData Comet::TrajectoryTracking::Sphere;
int Comet::TrajectoryTracking::TracedParticleNumber_GLOBAL[nZenithSurfaceElements*nAzimuthalSurfaceElements];
int Comet::TrajectoryTracking::TracedParticleNumber_LOCAL[nZenithSurfaceElements*nAzimuthalSurfaceElements];
int Comet::TrajectoryTracking::nTracedTrajectoriesPerElement;

void Comet::TrajectoryTracking::Init() {
  //set parameters of the sphere
  double x[3]={0.0,0.0,0.0};

  Sphere.SetGeneralSurfaceMeshParameters(nZenithSurfaceElements,nAzimuthalSurfaceElements);
  Sphere.SetSphereGeometricalParameters(x,TracingSurfaceRadius);

  //evaluate the total number of the traced trajectories and the number of the trajectories that originates in a particular surface element of the sphere
  nTracedTrajectoriesPerElement=nTotalTracedTrajectories/(nZenithSurfaceElements*nAzimuthalSurfaceElements);

  //init the counter of the traced trajectories
  for (int i=0;i<nZenithSurfaceElements*nAzimuthalSurfaceElements;i++) TracedParticleNumber_LOCAL[i]=0,TracedParticleNumber_GLOBAL[i]=0;
}

bool Comet::TrajectoryTracking::TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  bool res;
  long int nZenithElement,nAzimuthalElement,el;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //only those solar wind ions are traced, which trajectories are close to the surface of Mercury
  if ((spec>=_DUST_SPEC_) && (spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) {
    if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<pow(TracingSurfaceRadius,2)) return false;

    //retrive the number of the trajectories that already has been traced from the surface element
    Sphere.GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
    el=Sphere.GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

    if (TracedParticleNumber_GLOBAL[el]+TracedParticleNumber_LOCAL[el]>nTracedTrajectoriesPerElement) return false;
  }
  else return false;

  res=true; //PIC::ParticleTracker::TrajectoryTrackingCondition_default(x,v,spec,ParticleData);

  if (res==true) TracedParticleNumber_LOCAL[el]++;
  return res;
#endif
}

void Comet::TrajectoryTracking::UpdateParticleCounter() {
  int buffer[nZenithSurfaceElements*nAzimuthalSurfaceElements],i;

  MPI_Allreduce(TracedParticleNumber_LOCAL,buffer,nZenithSurfaceElements*nAzimuthalSurfaceElements,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  for (i=0;i<nZenithSurfaceElements*nAzimuthalSurfaceElements;i++) {
    TracedParticleNumber_GLOBAL[i]+=buffer[i];
    TracedParticleNumber_LOCAL[i]=0;
  }
}

double Comet::LossProcesses::PhotolyticReactionRate=0.0;
double Comet::LossProcesses::ElectronImpactRate=0.0;
double Comet::LossProcesses::ElectronTemeprature=0.0;


double Comet::LossProcesses::ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  long int nd;
  int i,j,k;
  double BackgroundPlasmaNumberDensity;
//  double PlasmaBulkVelocity[3],ElectronDensity;

  PhotolyticReactionRate=0.0;


//DEBUG -> no chemistry at all
  if ((spec!=_H2O_SPEC_) && (spec!=_H2_SPEC_) && (spec!=_H_SPEC_) && (spec!=_OH_SPEC_) && (spec!=_O_SPEC_) && (spec!=_O2_SPEC_)) {
   PhotolyticReactionAllowedFlag=false;
   return -1.0;
  }


  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
//  PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity,x,nd,node);
//  BackgroundPlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity(x,nd,node);

  PhotolyticReactionAllowedFlag=true;

  double HeliocentricDistance=0.0;

  static const double PhotolyticReactionRate_H2O=PhotolyticReactions::H2O::GetTotalReactionRate(HeliocentricDistance);
  static const double PhotolyticReactionRate_H2=PhotolyticReactions::H2::GetTotalReactionRate(HeliocentricDistance);
  static const double PhotolyticReactionRate_H=PhotolyticReactions::H::GetTotalReactionRate(HeliocentricDistance);
  static const double PhotolyticReactionRate_OH=PhotolyticReactions::OH::GetTotalReactionRate(HeliocentricDistance);
  static const double PhotolyticReactionRate_O=PhotolyticReactions::O::GetTotalReactionRate(HeliocentricDistance);
  static const double PhotolyticReactionRate_O2=PhotolyticReactions::O2::GetTotalReactionRate(HeliocentricDistance);


  switch (spec) {
    case _H2O_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_H2O;
      break;
    case _H2_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_H2;
      break;
    case _H_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_H;
      break;
    case _OH_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_OH;
      break;
    case _O_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_O;
      break;
    case _O2_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_O2;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: unknown specie");
    }


/*  //calcualte the rate due to the electron impact
  //characteristic values
  static const double ThermalElectronDensity=Europa::ElectronModel::ThermalElectronFraction;
  static const double HotElectronDensity=Europa::ElectronModel::HotElectronFraction;

  static const double HotElectronImpactRate_H2O=ElectronImpact::H2O::RateCoefficient(Europa::ElectronModel::HotElectronTemeprature)*HotElectronDensity;
  static const double ThermalElectronImpactRate_H2O=ElectronImpact::H2O::RateCoefficient(Europa::ElectronModel::ThermalElectronTemeprature)*ThermalElectronDensity;

  static const double HotElectronImpactRate_O2=ElectronImpact::O2::RateCoefficient(Europa::ElectronModel::HotElectronTemeprature)*HotElectronDensity;
  static const double ThermalElectronImpactRate_O2=ElectronImpact::O2::RateCoefficient(Europa::ElectronModel::ThermalElectronTemeprature)*ThermalElectronDensity;

  static const double HotElectronImpactRate_H2=ElectronImpact::H2::RateCoefficient(Europa::ElectronModel::HotElectronTemeprature)*HotElectronDensity;
  static const double ThermalElectronImpactRate_H2=ElectronImpact::H2::RateCoefficient(Europa::ElectronModel::ThermalElectronTemeprature)*ThermalElectronDensity;

  static const double HotElectronImpactRate_H=ElectronImpact::H::RateCoefficient(Europa::ElectronModel::HotElectronTemeprature)*HotElectronDensity;
  static const double ThermalElectronImpactRate_H=ElectronImpact::H::RateCoefficient(Europa::ElectronModel::ThermalElectronTemeprature)*ThermalElectronDensity;

  static const double HotElectronImpactRate_O=ElectronImpact::O::RateCoefficient(Europa::ElectronModel::HotElectronTemeprature)*HotElectronDensity;
  static const double ThermalElectronImpactRate_O=ElectronImpact::O::RateCoefficient(Europa::ElectronModel::ThermalElectronTemeprature)*ThermalElectronDensity;



  switch (spec){
  case _H2O_SPEC_:
    ElectronImpactRate=BackgroundPlasmaNumberDensity*(HotElectronImpactRate_H2O+ThermalElectronImpactRate_H2O);
    break;
  case _O2_SPEC_:
    ElectronImpactRate=BackgroundPlasmaNumberDensity*(HotElectronImpactRate_O2+ThermalElectronImpactRate_O2);
    break;
  case _H2_SPEC_:
    ElectronImpactRate=BackgroundPlasmaNumberDensity*(HotElectronImpactRate_H2+ThermalElectronImpactRate_H2);
    break;
  case _H_SPEC_:
    ElectronImpactRate=BackgroundPlasmaNumberDensity*(HotElectronImpactRate_H+ThermalElectronImpactRate_H);
    break;
  case _OH_SPEC_:
    ElectronImpactRate=0.0;
    break;
  case _O_SPEC_:
    ElectronImpactRate=BackgroundPlasmaNumberDensity*(HotElectronImpactRate_O+ThermalElectronImpactRate_O);
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown species");
  }
*/

  if (PhotolyticReactionRate+ElectronImpactRate<=0.0) {
    PhotolyticReactionAllowedFlag=false;
    return -1.0;
  }

  return 1.0/((PhotolyticReactionRate+ElectronImpactRate)*NumericalLossRateIncrease);  //use the "false" reaction event to increase the number of the dauter model particles. Account for this artificial correction in the ExospherePhotoionizationReactionProcessor
}


int Comet::LossProcesses::ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int *ReactionProductsList,nReactionProducts;
  double *ReactionProductVelocity;
  int ReactionChannel;
  bool PhotolyticReactionRoute;


  //init the reaction tables
  static bool initflag=false;
  static double TotalProductYeld_PhotolyticReaction[PIC::nTotalSpecies*PIC::nTotalSpecies];
  static double TotalProductYeld_ElectronImpact[PIC::nTotalSpecies*PIC::nTotalSpecies];

  /*  double HotElectronFraction=0.05;
  static const double ThermalElectronTemeprature=20.0;
  static const double HotElectronTemeprature=250.0;
  */

  if (initflag==false) {
    int iParent,iProduct;

    initflag=true;

    for (iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
      TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=0.0;
      TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=0.0;

      if (PhotolyticReactions::ModelAvailable(iParent)==true) {
        TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
      }

      /*      if (ElectronImpact::ModelAvailable(iParent)==true) {
        TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=
            Europa::ElectronModel::HotElectronFraction*ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,Europa::ElectronModel::HotElectronTemeprature) +
            Europa::ElectronModel::ThermalElectronFraction*ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,Europa::ElectronModel::ThermalElectronTemeprature);
	    }*/
    }
  }

  //determine the type of the reaction
  // PhotolyticReactionRoute=(rnd()<PhotolyticReactionRate/(PhotolyticReactionRate+ElectronImpactRate)) ? true : false;
  PhotolyticReactionRoute=true;

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
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  ParentTimeStep=node->block->GetLocalTimeStep(spec);
#else
  ParentTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif


  //account for the parent particle correction factor
  ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

  //the particle buffer used to set-up the new particle data
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  //copy the state of the initial parent particle into the new-daugher particle (just in case....)
  PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*)tempParticleData,ParticleData);

  for (int specProduct=0;specProduct<PIC::nTotalSpecies;specProduct++) {
    double ProductTimeStep,ProductParticleWeight;
    double ModelParticleInjectionRate,TimeCounter=0.0,TimeIncrement,ProductWeightCorrection=1.0/NumericalLossRateIncrease;
    int iProduct;
    long int newParticle;
    PIC::ParticleBuffer::byte *newParticleData;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
     ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[specProduct];
#else
     ProductParticleWeight=0.0;
     exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
     ProductTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[specProduct];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
     ProductTimeStep=node->block->GetLocalTimeStep(specProduct);
#else
     ProductTimeStep=0.0;
     exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

     ModelParticleInjectionRate=ParentParticleWeight/ParentTimeStep/ProductParticleWeight*((PhotolyticReactionRoute==true) ? TotalProductYeld_PhotolyticReaction[specProduct+spec*PIC::nTotalSpecies] : TotalProductYeld_ElectronImpact[specProduct+spec*PIC::nTotalSpecies]);

     //inject the product particles
     if (ModelParticleInjectionRate>0.0) {
       TimeIncrement=-log(rnd())/ModelParticleInjectionRate *rnd(); //<- *rnd() is to account for the injection of the first particle in the curent interaction

       while (TimeCounter+TimeIncrement<ProductTimeStep) {
         TimeCounter+=TimeIncrement;
         TimeIncrement=-log(rnd())/ModelParticleInjectionRate;

         //generate model particle with spec=specProduct
         bool flag=false;

         do {
           //generate a reaction channel
           if (PhotolyticReactionRoute==true) {
             PhotolyticReactions::GenerateReactionProducts(spec,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }
           else {
	     /*             if (rnd()<Europa::ElectronModel::HotElectronFraction) ElectronImpact::GenerateReactionProducts(spec,Europa::ElectronModel::HotElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
			    else ElectronImpact::GenerateReactionProducts(spec,Europa::ElectronModel::ThermalElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);*/
	   }

           //check whether the products contain species with spec=specProduct
           for (iProduct=0;iProduct<nReactionProducts;iProduct++) if (ReactionProductsList[iProduct]==specProduct) {
             flag=true;
             break;
           }
         }
         while (flag==false);


         //determine the velocity of the product specie
         double ProductParticleVelocity[3];

         for (int idim=0;idim<3;idim++) ProductParticleVelocity[idim]=vFinal[idim]+ReactionProductVelocity[idim+3*iProduct];

         //generate a particle
         PIC::ParticleBuffer::SetX(xFinal,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetV(ProductParticleVelocity,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetI(specProduct,(PIC::ParticleBuffer::byte*)tempParticleData);

         #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
         #endif

         //apply condition of tracking the particle
         #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
         PIC::ParticleTracker::InitParticleID(tempParticleData);
         PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xInit,xFinal,spec,tempParticleData);
         #endif


         //get and injection into the system the new model particle
         newParticle=PIC::ParticleBuffer::GetNewParticle();
         newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
         memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

         _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,ProductTimeStep-TimeCounter,node,true);
       }
     }

     }


  return (rnd()<1.0/NumericalLossRateIncrease) ? _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_ : _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  }

//distribute initial dust grain velocity
void Comet::DustInitialVelocity::GetInitialVelocity(double *v,double *x,int iInjectionFace) {
  double *ExternalNormal,*e0,*e1,sinTheta,cosTheta,phi,sinPhi,cosPhi;
  int idim;

  ExternalNormal=CutCell::BoundaryTriangleFaces[iInjectionFace].ExternalNormal;
  e0=CutCell::BoundaryTriangleFaces[iInjectionFace].e0Orthogonal;
  e1=CutCell::BoundaryTriangleFaces[iInjectionFace].e1Orthogonal;

  cosTheta=rnd();
  sinTheta=sqrt(1.0-cosTheta*cosTheta);
  phi=2.0*Pi*rnd();
  sinPhi=sin(phi);
  cosPhi=cos(phi);

  for (idim=0;idim<3;idim++) v[idim]=InjectionConstantVelocity*(cosTheta*ExternalNormal[idim]+sinTheta*(sinPhi*e0[idim]+cosPhi*e1[idim]));

  if (VelocityModelMode==Mode::RotationBody) {
    //add velocity due to rotation of the body
    double vRot[3],c;

    c=2.0*Pi/RotationPeriod;
    Vector3D::CrossProduct(vRot,x,RotationAxis);

    for (idim=0;idim<3;idim++) v[idim]+=c*vRot[idim];
  }
}
