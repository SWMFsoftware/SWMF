/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

#undef _PIC_MODEL__3DGRAVITY__MODE_
#define _PIC_MODEL__3DGRAVITY__MODE_ _PIC_MODEL__3DGRAVITY__MODE__OFF_

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

static bool probabilityFunctionDefinedJet=false,probabilityFunctionDefined=false,probabilityFunctionDefinedWaist=false,probabilityFunctionDefinedHartley2=false,probabilityFunctionDefinedNASTRAN=false;
//static double productionDistributionJet[360][180],cumulativeProductionDistributionJet[360][180];
static double productionDistributionJet[6000],cumulativeProductionDistributionJet[6000];
static double productionDistributionWaist[6000],cumulativeProductionDistributionWaist[6000];
static double productionDistributionHartley2[6000],cumulativeProductionDistributionHartley2[6000];
static double productionDistributionNASTRAN[30000],cumulativeProductionDistributionNASTRAN[30000];
static double productionDistribution[180],cumulativeProductionDistribution[180];
static double angle;
static double azimuthCenter;
static double zenithCenter;
static cInternalRotationBodyData* Nucleus;

double subSolarPointAzimuth=0.0; //53.0*Pi/180; //0.0;

//char Exosphere::SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_]="2008-12-20T00:00:00"; //"2011-04-13T00:00:00" ;//"2009-00-00T00:00:00";

//parameters of the soruce processes
/*double Exosphere::SourceProcesses::ThermalDesorption::uThermal=1.85*eV2J;
double Exosphere::SourceProcesses::ThermalDesorption::VibrationalFrequency=1.0E13;

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::PhotonFlux_1AU=2.0E14*1.0E4;  //Killen-2012-JGR, Yakshinskii+Madey-1999-?
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::CrossSection=3.0E-21*1.0E-4;//Satantos-2010-? ; 3.0E-20*1.0E-4;  //Killen-2012-JGR, Yakshinskii+Madey-1999-?
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;*/

/*
double Exosphere::SourceProcesses::ImpactVaporization::SourceRate=1.89E22; //////1.79e21; //Killen-2012-JGR   ;1.1e22;  2.05e22 IV for Sarantos 2010
double Exosphere::SourceProcesses::ImpactVaporization::HeliocentricDistance=1.0*_AU_;
double Exosphere::SourceProcesses::ImpactVaporization::SourceRatePowerIndex=0.0;
double Exosphere::SourceProcesses::ImpactVaporization::SourceTemeprature=6000.0; //Killen-2012-JGR ;2500.0;
*/

/*double Exosphere::SourceProcesses::SolarWindSputtering::Yield=0.1;
double Exosphere::SourceProcesses::SolarWindSputtering::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::SolarWindSputtering::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;*/


//typical parameters of solar wind
/*const double Exosphere::swVelocity_Typical[3]={-420.0E3,0.0,0.0};
const double Exosphere::swB_Typical[3]={-12.9E-9,4.71E-9,10.29E-9};
const double Exosphere::swTemperature_Typical=0.174e6,Exosphere::swNumberDensity_Typical=60.0E6;*/


//interaction with the surface
//const double Exosphere::SurfaceInteraction::AccomodationCoefficient=0.2;

//sticking probability of sodium atoms

/*
double Exosphere::SurfaceInteraction::SodiumStickingProbability(double Temp) {
  if (Temp<300.0) return 1.0;
  if (Temp<650.0) return 1.0-(Temp-300.0)/350.0;

  return 0.0;
}
*/

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
  ElectricallyChargedDust::minDustRadius=0.01*_MICROMETER_;
  ElectricallyChargedDust::maxDustRadius=100.0*_MICROMETER_;
  ElectricallyChargedDust::Sampling::SetDustSamplingIntervals(10);
  ElectricallyChargedDust::GrainVelocityGroup::minGrainVelocity=0.01;
  ElectricallyChargedDust::GrainVelocityGroup::maxGrainVelocity=1000.0;
  ElectricallyChargedDust::TotalMassDustProductionRate=1.0;
  ElectricallyChargedDust::SizeDistribution::PowerIndex=4.0;
  ElectricallyChargedDust::Init_BeforeParser();
#endif

}

void Comet::Init_AfterParser() {

  //set up the Chamberlen model
  double ExosphereEscapeRate[PIC::nTotalSpecies],ExospehreTemsprature[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) { //ExosphereEscapeRate[spec]=0.0,ExospehreTemsprature[spec]=1000.0;
    ExosphereEscapeRate[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceRate[spec];
    ExospehreTemsprature[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceTemeprature[spec];
  }

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model                                                                                                                                                                           
  ElectricallyChargedDust::Init_AfterParser();
#endif
  
  //init Gravity
#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  InitGravityData();
#endif

  /*
  Exosphere::ChamberlainExosphere::Init(ExosphereEscapeRate,ExospehreTemsprature);



  //set up the model that collected the column integrals at the subsolar point of the limn as a function of the phase angle
  Sampling::SubsolarLimbColumnIntegrals::init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Sampling::SubsolarLimbColumnIntegrals::CollectSample);

  //set up the model that prints the column integrals in the anti-solar direction
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,AntiSolarDirectionColumnMap::Print);

  //set up sampling of velocity distribution functions
  Comet::Sampling::VelocityDistribution::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Comet::Sampling::VelocityDistribution::Sampling,Comet::Sampling::VelocityDistribution::OutputSampledData);

  //call init function of the exospheric model
  Exosphere::Init_AfterParser();

  //init the model of calcualting the integrals that correspond to the Kaguya's TVIS observations
  Comet::Sampling::Kaguya::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Comet::Sampling::Kaguya::TVIS::OutputModelData);
  */
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


double SodiumStickingProbability(double& ReemissionParticleFraction,double Temp) {
  double res=-1.0;


  //define parametes of the sticking functon
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_           0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__TEMPERATURE_LIMIT_  1
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_   2

#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_     1.0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_          1.0
#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY__TEMPARATURE_LIMIT_       200.0


#define _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_



  ReemissionParticleFraction=_EXOSPHERE_SODIUM_STICKING_PROBABILITY__REEMISSION_FRACTION_;

#if _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__CONSTANT_
  res= _EXOSPHERE_SODIUM_STICKING_PROBABILITY__CONSTANT_VALUE_;

#elif _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__TEMPERATURE_LIMIT_
  res= (Temp<_EXOSPHERE_SODIUM_STICKING_PROBABILITY__TEMPARATURE_LIMIT_) ? 1.0 : 0.0;

#elif _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE_ == _EXOSPHERE_SODIUM_STICKING_PROBABILITY_MODE__YAKSHINSKY2005SS_
  static const double tmin=100.0;
  static const double tmax=500.0;
  static const double dt=5.0;

  static const int nPoints=80;

  struct cDataStruct {
    double t,s;
  };


  static const cDataStruct data[nPoints]={    //digitized from  Yakshinskiy-2005-SS
      {100.00000,0.99830}, {105.00000,0.96768}, {110.00002,0.94512}, {114.99999,0.92417}, {119.99999,0.90322}, {125.00000,0.88065}, {130.00000,0.86293}, {135.00000,0.84681}, {140.00000,0.82586}, {145.00000,0.80813},
      {150.00000,0.78718}, {154.99998,0.76623}, {160.00000,0.75011}, {165.00000,0.73077}, {170.00000,0.71305}, {175.00000,0.69371}, {179.99998,0.67759}, {184.99998,0.65986}, {189.99998,0.64375}, {195.00000,0.62763},
      {200.00000,0.61151}, {205.00000,0.59701}, {210.00002,0.57928}, {214.99998,0.56639}, {220.00000,0.55188}, {225.00000,0.53899}, {229.99998,0.52287}, {235.00000,0.51159}, {240.00000,0.49870}, {245.00000,0.48903},
      {249.99997,0.47936}, {255.00000,0.47130}, {260.00000,0.46325}, {265.00000,0.45035}, {270.00000,0.43907}, {275.00000,0.42618}, {280.00000,0.41651}, {285.00000,0.40523}, {290.00000,0.39717}, {295.00000,0.38911},
      {300.00000,0.38266}, {305.00000,0.36977}, {309.99997,0.35688}, {315.00000,0.34882}, {319.99997,0.33915}, {325.00000,0.33432}, {329.99997,0.32465}, {334.99997,0.31659}, {340.00000,0.31014}, {345.00003,0.30370},
      {350.00000,0.29564}, {354.99997,0.28919}, {360.00000,0.28274}, {365.00000,0.27791}, {370.00003,0.27146}, {375.00000,0.26502}, {379.99997,0.26018}, {384.99997,0.25535}, {390.00000,0.25051}, {395.00003,0.24568},
      {400.00000,0.24084}, {405.00000,0.23762}, {409.99997,0.23278}, {414.99997,0.22956}, {420.00000,0.22473}, {425.00000,0.22150}, {430.00000,0.21828}, {435.00003,0.21506}, {440.00000,0.21022}, {444.99997,0.20700},
      {449.99994,0.20377}, {454.99997,0.20377}, {460.00000,0.20055}, {465.00000,0.19894}, {470.00000,0.19572}, {475.00000,0.19411}, {479.99997,0.19088}, {485.00000,0.19088}, {490.00000,0.18766}, {494.99997,0.18605}
  };

  if (Temp<tmin) res=data[0].s;
  else if (Temp>tmax) res=data[nPoints-1].s;
  else {
    int n;
    double c;

    n=(int)((Temp-tmin)/dt);
    c=(Temp-tmin-dt*n)/dt;

    res=(1.0-c)*data[n].s+c*data[n+1].s;
  }


  return res;
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif

  return res;
}


double Exosphere::SurfaceInteraction::StickingProbability(int spec, double& ReemissionParticleFraction,double Temp) {
  double res=0.0;

   switch (spec) {
   case _NA_SPEC_: //case _NAPLUS_SPEC_:
     res=SodiumStickingProbability(ReemissionParticleFraction,Temp);
     break;
   default:
     exit(__LINE__,__FILE__,"the option is not implemented");
   }

   return res;
}


//surface temeprature of the planet
double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {
  /*
  //determine if the point on the night side of the Comet
  if (CosSubSolarAngle<0.0) return 100.0;

  //determine if the point is within the shadow of the Earth
  if (Comet::EarthShadowCheck(x_LOCAL_SO_OBJECT)==true) return 100.0;

  //return the day-side temeprature
  return 280*pow(CosSubSolarAngle,0.25)+100.0;*/

  const double minTemp[]={172.0,163.0,150.0,145.0,139.0};
  double res,r,zenith,azimuth;
  int angle;

  //determine if the point on the night side of the Comet
  if (CosSubSolarAngle<0.0) return minTemp[Comet::ndist];
  
  //return the day-side temeprature
  angle=(int) (acos(CosSubSolarAngle)*180.0/Pi);
  if(angle>89) angle=89;
  res=(SurfaceTemp[angle][Comet::ndist+1]>minTemp[Comet::ndist]) ?  SurfaceTemp[angle][Comet::ndist+1] : minTemp[Comet::ndist];

  return res;
  //  return 180;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  /*int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\",  \"Mean Speed Along the Line of Sight(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec),PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=2;
  }
  
  if (_NA_SPEC_>=0) Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=2*_NA_SPEC_;

  //sodium emission
  if (_NA_SPEC_>=0) {
    Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=nVariables;
    Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=nVariables+1;

    if (vlist!=NULL) sprintf(vlist,"%s,  \"Sodium Emission(5891_58A)\",  \"Sodium Emission(5897_56A)\"",vlist);
    nVariables+=2;
  }

  return nVariables;
  */}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (res[cnt]>0.0) res[cnt+1]/=res[cnt];
    cnt+=2;
  }
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  /*  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;

  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double BulkVelocity_SO[3],v_LOCAL_SO_FROZEN[3],rHeliocentric,vHeliocentric;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity_SO,spec);

    v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+BulkVelocity_SO[0]+
        Exosphere::RotationVector_SO_FROZEN[1]*x[2]-Exosphere::RotationVector_SO_FROZEN[2]*x[1];

    v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+BulkVelocity_SO[1]-
        Exosphere::RotationVector_SO_FROZEN[0]*x[2]+Exosphere::RotationVector_SO_FROZEN[2]*x[0];

    v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+BulkVelocity_SO[2]+
        Exosphere::RotationVector_SO_FROZEN[0]*x[1]-Exosphere::RotationVector_SO_FROZEN[1]*x[0];

    rHeliocentric=sqrt(pow(x[0]-Exosphere::xObjectRadial,2)+(x[1]*x[1])+(x[2]*x[2]));
    vHeliocentric=(
        (v_LOCAL_SO_FROZEN[0]*(x[0]-Exosphere::xObjectRadial))+
        (v_LOCAL_SO_FROZEN[1]*x[1])+(v_LOCAL_SO_FROZEN[2]*x[2]))/rHeliocentric;


    //brightness of the exospheric sodium
    if (spec==_NA_SPEC_) {
      //check if the point is outside of the Comet's and Earth's shadows
      if ( (Comet::EarthShadowCheck(x)==false) && ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_MOON_)*_RADIUS_(_MOON_))) ) {
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
      else res[cnt++]=0.0,res[cnt++]=0.0;
    }
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
*/}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  double res=0.0;

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceDouble State[6],ltlocal;
  double EccentricityVector[3];
  double Speed2,a,c,absEccentricity;
  const double GravitationalParameter=GravityConstant*_MASS_(_EARTH_);
  double vComet[3],xComet[3],rGeocentric=0.0;
  int idim;


  spkezr_c("Comet",EphemerisTime,"MSGR_HCI","none","Earth",State,&ltlocal);

  for (idim=0;idim<3;idim++) {
    xComet[idim]=State[idim]*1.0E3;
    vComet[idim]=State[idim+3]*1.0E3;

    rGeocentric+=pow(xComet[idim],2);
  }

  rGeocentric=sqrt(rGeocentric);
  Speed2=vComet[0]*vComet[0]+vComet[1]*vComet[1]+vComet[2]*vComet[2];
  c=xComet[0]*vComet[0]+xComet[1]*vComet[1]+xComet[2]*vComet[2];

  for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
    EccentricityVector[idim]=Speed2/GravitationalParameter*xComet[idim] - c/GravitationalParameter*vComet[idim] - xComet[idim]/rGeocentric;
    absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
    a+=EccentricityVector[idim]*xComet[idim];
  }

  absEccentricity=sqrt(absEccentricity);
  res=acos(a/(absEccentricity*rGeocentric));

  if (c<0.0) res=2.0*Pi-res;
#endif

  return res;
}


double Comet::GetTotalProductionRateBjorn(int spec,void *SphereDataPointer){
  return Comet::Bjorn_SourceRate[spec];

  //    return 1.0e26; //AT THIS STAGE, WE ARE ONLY TESTING THE FUNCTION GENERATEPARTICLEPROPERTIESBJORN BELOW
}

bool Comet::GenerateParticlePropertiesBjorn(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int i;
  double rate,TableTotalProductionRate,totalSurface,gamma,SubSolarAngle,ProjectedAngle;
  const double NightSideProduction[5]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0};
  double x[3],n[3],rSphere,*x0Sphere;

    if (probabilityFunctionDefined==false) {
      for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
      TableTotalProductionRate+=ProductionRate[i][2+Comet::ndist];
      //totalSurface+=2*Pi*pow(Comet::Rnucleus,2.0)*(cos((double) (i*Pi/180))-cos((double) ((i+1)*Pi/180)));
    }
    //totalSurface=2*totalSurface;
    
    //Computation of the probability distribution of the production
    for (i=0;i<90;i++) {
      productionDistribution[i]=ProductionRate[i][2+Comet::ndist]/(TableTotalProductionRate/(1.0-NightSideProduction[Comet::ndist]*2))+NightSideProduction[Comet::ndist]*(cos((double) (i*Pi/180))-cos(double ((i+1)*Pi/180)));
      if (i==0) {
	cumulativeProductionDistribution[i]=productionDistribution[i];
      }else{
	cumulativeProductionDistribution[i]=cumulativeProductionDistribution[i-1]+productionDistribution[i];
      }  
    }
    for (i=90;i<180;i++) {
      productionDistribution[i]=NightSideProduction[Comet::ndist]*(cos((double) (i*Pi/180))-cos(double ((i+1)*Pi/180)));
      cumulativeProductionDistribution[i]=cumulativeProductionDistribution[i-1]+productionDistribution[i];
    }
    probabilityFunctionDefined=true;
    }

  //Computation of the segment where the particle will be created
    gamma=rnd();
    i=0;
    while (gamma>cumulativeProductionDistribution[i]){
      i++;
    }


  //Computation of the angles where the particle will be created
  SubSolarAngle=acos(cos((double) (i*Pi/180))-rnd()*(cos((double) (i*Pi/180))-cos((double) ((i+1)*Pi/180))));
  ProjectedAngle=rnd()*2*Pi;

  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];

  /*  
  Sphere->GetSphereGeometricalParameters(x0Sphere,rSphere);

  ExternalNormal[0]=-cos(SubSolarAngle);
  ExternalNormal[1]=-sin(SubSolarAngle)*cos(ProjectedAngle);
  ExternalNormal[2]=-sin(SubSolarAngle)*sin(ProjectedAngle);

  
  x_LOCAL_IAU_OBJECT[0]=rSphere*cos(SubSolarAngle);
  x_LOCAL_IAU_OBJECT[1]=rSphere*sin(SubSolarAngle)*cos(ProjectedAngle);
  x_LOCAL_IAU_OBJECT[2]=rSphere*sin(SubSolarAngle)*sin(ProjectedAngle);
  */

  Sphere->GetSphereGeometricalParameters(x0Sphere,rSphere);

  n[0]=-cos(SubSolarAngle);
  n[1]=-sin(SubSolarAngle)*cos(ProjectedAngle);
  n[2]=-sin(SubSolarAngle)*sin(ProjectedAngle);
  
  x[0]=rSphere*cos(SubSolarAngle);
  x[1]=rSphere*sin(SubSolarAngle)*cos(ProjectedAngle);
  x[2]=rSphere*sin(SubSolarAngle)*sin(ProjectedAngle);
		         
  ExternalNormal[0]=cos(subSolarPointAzimuth)*n[0]+sin(subSolarPointAzimuth)*n[1];
  ExternalNormal[1]=-sin(subSolarPointAzimuth)*n[0]+cos(subSolarPointAzimuth)*n[1];
  ExternalNormal[2]=n[2];
  
  x_LOCAL_IAU_OBJECT[0]=cos(subSolarPointAzimuth)*x[0]+sin(subSolarPointAzimuth)*x[1];
  x_LOCAL_IAU_OBJECT[1]=-sin(subSolarPointAzimuth)*x[0]+cos(subSolarPointAzimuth)*x[1];
  x_LOCAL_IAU_OBJECT[2]=x[2];


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
  
  SurfaceTemperature=GetSurfaceTemeprature(cos(SubSolarAngle),x_LOCAL_SO_OBJECT);
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);

  
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

bool Comet::Radius(double &r,double x){
  double rSphere=2000.0;

  //Sphere->GetSphereGeometricalParameters(x0,rSphere);
  
  if ((x>=-rSphere-1.0e-15) && (x<=rSphere+1.0e-15)) {
    r=sqrt(rSphere*rSphere-x*x);  
    return true;
    }
  
  //Hartley 2
  /*  double boundary;
  
    if ((x>=-1165.0-1.0e-15) && (x<=1165.0+1.0e-15)) {
    boundary=sqrt(1980.0*1980.0-x*x*1980.0*2.0/2330.0*1980.0*2.0/2330.0)/1980.0*2330.0/(2.0*1.3)-(2330.0/2.0-600.0)*0.5*(1+cos((-x*1980.0*2.0/2330.0+200.0)/1800.0*Pi));
    
    r=((boundary>=0)? sqrt(boundary*boundary):0.0);
    return true;
    }
  */

  return false;
}

//bool Comet::GenerateParticlePropertiesHartley2(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalRotationBodyData* Nucleus) {
bool Comet::GenerateParticlePropertiesHartley2(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[5]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=1.064*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double l[3]={1.0,0.0,0.0};
  double x0[3]={0.0,0.0,0.0};
  double rSphere=1980.0;

  

    if (probabilityFunctionDefinedHartley2==false) {    
      Nucleus=(cInternalRotationBodyData *) Sphere;
      Nucleus= (cInternalRotationBodyData *) malloc(sizeof(cInternalRotationBodyData));
      Nucleus->SetGeneralSurfaceMeshParameters(60,100);
      //Nucleus->SetGeometricalParameters(x0,l,-rSphere,rSphere,Radius);
      Nucleus->SetGeometricalParameters(x0,l,-1165.0,1165.0,Radius);
      for (TableTotalProductionRate=0.0,i=0;i<90;i++) {
	TableTotalProductionRate+=ProductionRate[i][2+Comet::ndist];
      }
  
      positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth);
      positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth);
      positionSun[2]=0.0;

      // Nucleus->GetSphereGeometricalParameters(x0,l,xmin,xmax);
      totalSurfaceElementsNumber=Nucleus->GetTotalSurfaceElementsNumber();

      total=0.0;      
      for (i=0;i<totalSurfaceElementsNumber;i++) {
	  
	Nucleus->GetSurfaceElementNormal(norm,i);
	Nucleus->GetSurfaceElementMiddlePoint(x,i);

	  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
	    c+=norm[idim]*(positionSun[idim]-x[idim]);
	    X+=pow(positionSun[idim]-x[idim],2.0);
	  }

	  if(c<0) {
	    productionDistributionHartley2[i]=NightSideProduction[Comet::ndist]*Nucleus->GetSurfaceElementArea(i)/(2*Pi*rSphere*rSphere);
	    //productionDistributionHartley2[i]=1.0;
	    total+=productionDistributionHartley2[i];
	  }else{
	    double angleProd;
	    int angleProdInt;
	    angleProd=acos(c/sqrt(X))*180/Pi;
	    angleProdInt=(int) angleProd;
	    productionDistributionHartley2[i]=ProductionRate[angleProdInt][2+Comet::ndist]/(TableTotalProductionRate/(1.0-NightSideProduction[Comet::ndist]*2))/(2*Pi*rSphere*rSphere*(cos(angleProdInt*Pi/180)-cos(angleProdInt*Pi/180+Pi/180)))*Nucleus->GetSurfaceElementArea(i)+NightSideProduction[Comet::ndist]*Nucleus->GetSurfaceElementArea(i)/(2*Pi*rSphere*rSphere);
	    //productionDistributionHartley2[i]=0.0;
	    total+=productionDistributionHartley2[i];
	  }
      }
      printf("total=%e \n",total);
      
      cumulativeProductionDistributionHartley2[0]=0.0;
      for (i=0;i<totalSurfaceElementsNumber;i++) {
	if (i==0) {
	  cumulativeProductionDistributionHartley2[i]+=productionDistributionHartley2[i]/total;
	}else{
	  cumulativeProductionDistributionHartley2[i]=cumulativeProductionDistributionHartley2[i-1]+productionDistributionHartley2[i]/total;
	}
	
      }
    probabilityFunctionDefinedHartley2=true;
    }
    
    //Computation of the segment where the particle will be created
    gamma=rnd();
    i=0;
    while (gamma>cumulativeProductionDistributionHartley2[i]){
      i++;
    }

    Nucleus->GetSurfaceElementIndex(nAxisElement,nAzimuthalElement,i);

  
  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
  
  Nucleus->GetSurfaceElementRandomPoint(x_LOCAL_IAU_OBJECT,nAxisElement,nAzimuthalElement);
  Nucleus->GetSurfaceElementNormal(ExternalNormal,i);
  
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

  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  
  /*  // only normal velocity
  v_LOCAL_IAU_OBJECT[0]=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*SurfaceTemperature))*ExternalNormal[0];
  v_LOCAL_IAU_OBJECT[1]=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*SurfaceTemperature))*ExternalNormal[1];
  v_LOCAL_IAU_OBJECT[2]=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*SurfaceTemperature))*ExternalNormal[2];
  */

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

long int Comet::InjectionBoundaryModel_Limited() {
  int spec;
  long int res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=InjectionBoundaryModel_Limited(spec);

  return res;
}

long int Comet::InjectionBoundaryModel_Limited(int spec) {
  cInternalSphericalData *Sphere;
  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_SO_OBJECT[3],x_IAU_OBJECT[3],v_SO_OBJECT[3],v_IAU_OBJECT[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag=false;
  int SourceProcessID;

  double totalProductionRate=Comet::GetTotalProductionRateBjornNASTRAN(spec,Sphere);

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
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
  //  int _EXOSPHERE__SOURCE_MAX_ID_VALUE_=0;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_=0;

  //calcualte probabilities of each source processes                                                                 
  double TotalFlux,FluxSourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_]; //,ProbabilitySourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];                                                                                               
int iSource;

for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) FluxSourceProcess[iSource]=0.0; //,ProbabilitySourceProcess[iSource]=0.0;                                                                                         

TotalFlux=totalProductionRate;

//only Used defined source here since we only want the Bjorn model so far
//calculate the source rate due to user defined source functions                                                   
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]=Comet::GetTotalProductionRateBjornNASTRAN(spec,Sphere);

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
 if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   static double GrainInjectedMass=0.0;
   PIC::Mesh::cDataBlockAMR *block;
   double GrainRadius,GrainMass,GrainWeightCorrection;
   int GrainVelocityGroup;

   GrainInjectedMass+=ElectricallyChargedDust::TotalMassDustProductionRate*LocalTimeStep;

   while (GrainInjectedMass>0.0) {
     startNode=NULL;

     double CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];

     //generate a particle                                                                                             
     char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
     PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

     flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,tempParticleData);
     ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
     GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
     GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
     SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
     if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;

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
    SourceProcessID=(int)(rnd()*(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_));
  }
  while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());

  double CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];

  //generate a particle                                                                                             
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

  //to satisfy the compiler and fit the while structure                                                             
  if (false) {}

  //Add the user defined particle gineration                                                                        
  flag=Comet::GenerateParticlePropertiesBjornNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,tempParticleData);
  SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_;
  if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_Bjorn_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;

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

double Comet::GetTotalProductionRateBjornNASTRAN(int spec, cInternalSphericalData* Sphere){
  return 1.0e24; //AT THIS STAGE, WE ARE ONLY TESTING THE FUNCTION GENERATEPARTICLEPROPERTIESBJORN BELOW
}


bool Comet::GenerateParticlePropertiesBjornNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  const double NightSideProduction[5]={5.8/100.0,7.0/100.0,9.2/100.0,10.4/100.0,11.6/100.0};
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.3*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement;
  long int totalSurfaceElementsNumber,i;
  double rSphere=1980.0;
  
  if (probabilityFunctionDefinedNASTRAN==false) {       
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
      
      if(c<0) {
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
    }
    //    printf("total=%e \n",total);
    
    cumulativeProductionDistributionNASTRAN[0]=0.0;
    for (i=0;i<totalSurfaceElementsNumber;i++) {
      if (i==0) {
	cumulativeProductionDistributionNASTRAN[i]+=productionDistributionNASTRAN[i]/total;
      }else{
	cumulativeProductionDistributionNASTRAN[i]=cumulativeProductionDistributionNASTRAN[i-1]+productionDistributionNASTRAN[i]/total;
      }
      
    }
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
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,PIC::Mesh::mesh.EPS);
  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,1e-4);
  //  CutCell::BoundaryTriangleFaces[i].GetRandomPosition(x_LOCAL_IAU_OBJECT,1.0e-1);
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

  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=-10.0*ExternalNormal[idim];
  else PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
#else
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
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

double Comet::GetTotalProductionRateJet(int spec,void *SphereDataPointer){
  return Comet::Jet_SourceRate[spec];
  //return 1.0e27; //AT THIS STAGE, WE ARE ONLY TESTING THE FUNCTION GENERATEPARTICLEPROPERTIESBJORN BELOW
}

//NUCLEUS
bool Comet::GenerateParticlePropertiesJet(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int i,j;
  double total=0.0,TableTotalProductionRate,totalSurface,gamma,gamma2;
  double azimuth,zenith,azimuthTable[360],zenithTable[180],rSphere,*x0Sphere;
  double l[3]={1.0,0.0,0.0};
  double x0[3]={0.0,0.0,0.0};
  double c=0.0,X=0.0,cosSubSolarAngle;
  int totalSurfaceElementsNumber,nAxisElement,nAzimuthalElement,idim;
  static double positionSun[3];
  double x[3];
  double HeliocentricDistance=1.064*_AU_;

   if (probabilityFunctionDefinedJet==false) {    
  //Computation of the probability distribution of the production. ATTENTION CASE OF ZENITH GOING OVER BOUNDARIES NOT IMPLEMENTED
     Nucleus=(cInternalRotationBodyData *) Sphere;
     Nucleus= (cInternalRotationBodyData *) malloc(sizeof(cInternalRotationBodyData));
     Nucleus->SetGeneralSurfaceMeshParameters(60,100);
     Nucleus->SetGeometricalParameters(x0,l,-1165.0,1165.0,Radius);
     // Nucleus->GetSphereGeometricalParameters(x0,l,xmin,xmax);
     totalSurfaceElementsNumber=Nucleus->GetTotalSurfaceElementsNumber();
     
     positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth);
     positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth);
     positionSun[2]=0.0;

     total=0.0;      
     for (i=0;i<totalSurfaceElementsNumber;i++) {     
       Nucleus->GetSurfaceElementMiddlePoint(x,i);
       
       if (x[0]<=1165.0 && x[0]>400.0) {
	 productionDistributionJet[i]=1.0;
	 total+=1.0;
       }else{
	 productionDistributionJet[i]=0.0;
       }
     }
     cumulativeProductionDistributionJet[0]=0.0;
     for (i=0;i<totalSurfaceElementsNumber;i++) {
       if (i==0) {
	 cumulativeProductionDistributionJet[i]+=productionDistributionJet[i]/total;
       }else{
	 cumulativeProductionDistributionJet[i]=cumulativeProductionDistributionJet[i-1]+productionDistributionJet[i]/total;
       }
       
     }
     
     probabilityFunctionDefinedJet=true;
   }

   //Computation of the area where the particle will be created
    gamma=rnd();
    i=0;
    while (gamma>cumulativeProductionDistributionJet[i]){
      i++;
    }

    Nucleus->GetSurfaceElementIndex(nAxisElement,nAzimuthalElement,i);
    
    
    //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
    double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
    
  Nucleus->GetSurfaceElementRandomPoint(x_LOCAL_IAU_OBJECT,nAxisElement,nAzimuthalElement);
  Nucleus->GetSurfaceElementNormal(ExternalNormal,i);
  
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
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  
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


//NUCLEUS
bool Comet::GenerateParticlePropertiesWaist(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int i,j;
  double total=0.0,TableTotalProductionRate,totalSurface,gamma,gamma2;
  double azimuth,zenith,azimuthTable[360],zenithTable[180],rSphere,*x0Sphere;
  double l[3]={1.0,0.0,0.0};
  double x0[3]={0.0,0.0,0.0};
  double c=0.0,X=0.0,cosSubSolarAngle;
  int totalSurfaceElementsNumber,nAxisElement,nAzimuthalElement,idim;
  static double positionSun[3];
  double x[3];
  double HeliocentricDistance=1.064*_AU_;

   if (probabilityFunctionDefinedWaist==false) {    
  //Computation of the probability distribution of the production. ATTENTION CASE OF ZENITH GOING OVER BOUNDARIES NOT IMPLEMENTED
     Nucleus=(cInternalRotationBodyData *) Sphere;
     Nucleus= (cInternalRotationBodyData *) malloc(sizeof(cInternalRotationBodyData));
     Nucleus->SetGeneralSurfaceMeshParameters(60,100);
     Nucleus->SetGeometricalParameters(x0,l,-1165.0,1165.0,Radius);
     // Nucleus->GetSphereGeometricalParameters(x0,l,xmin,xmax);
     totalSurfaceElementsNumber=Nucleus->GetTotalSurfaceElementsNumber();
     
     positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth);
     positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth);
     positionSun[2]=0.0;

     total=0.0;      
     for (i=0;i<totalSurfaceElementsNumber;i++) {     
       Nucleus->GetSurfaceElementMiddlePoint(x,i);
       
       if (x[0]<=400.0 && x[0]>0.0 && x[2]<345/2.0 && x[2]>-345/2.0 && x[1]>0) {
	 productionDistributionWaist[i]=1.0;
	 total+=1.0;
       }else{
	 productionDistributionWaist[i]=0.0;
       }
     }
     cumulativeProductionDistributionWaist[0]=0.0;
     for (i=0;i<totalSurfaceElementsNumber;i++) {
       if (i==0) {
	 cumulativeProductionDistributionWaist[i]+=productionDistributionWaist[i]/total;
       }else{
	 cumulativeProductionDistributionWaist[i]=cumulativeProductionDistributionWaist[i-1]+productionDistributionWaist[i]/total;
       }
       
     }
     
     probabilityFunctionDefinedWaist=true;
   }

   //Computation of the area where the particle will be created
    gamma=rnd();
    i=0;
    while (gamma>cumulativeProductionDistributionWaist[i]){
      i++;
    }

    Nucleus->GetSurfaceElementIndex(nAxisElement,nAzimuthalElement,i);
    
    
    //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
    double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
    
  Nucleus->GetSurfaceElementRandomPoint(x_LOCAL_IAU_OBJECT,nAxisElement,nAzimuthalElement);
  Nucleus->GetSurfaceElementNormal(ExternalNormal,i);
  
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
  
  for (c=0.0,X=0.0,idim=0;idim<3;idim++){
    c+=ExternalNormal[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
    X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
  }
  cosSubSolarAngle=c/sqrt(X);
  
  SurfaceTemperature=GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  
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


/* SPHERE
bool Comet::GenerateParticlePropertiesJet(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int i,j;
  double total=0.0,TableTotalProductionRate,totalSurface,gamma,gamma2;
  double azimuth,zenith,azimuthTable[360],zenithTable[180],rSphere,*x0Sphere;

   if (probabilityFunctionDefinedJet==false) {    
  //Computation of the probability distribution of the production. ATTENTION CASE OF ZENITH GOING OVER BOUNDARIES NOT IMPLEMENTED
     
     for (i=0;i<360;i++) {
    for (j=0;j<180;j++) {
      if ((azimuthCenter+angle<360) && (azimuthCenter-angle>=0) && (zenithCenter+angle<180) && (zenithCenter-angle>=0)) { //General Case
	if ((azimuthCenter+angle>=i) && (azimuthCenter-angle<=i) && (zenithCenter+angle>=j) && (zenithCenter-angle<=j)) {
	  productionDistributionJet[i][j]=1.0;
	  total+=1.0;
	}else{
	  productionDistributionJet[i][j]=0.0;
	}
      }else if ((azimuthCenter+angle>360) && (azimuthCenter-angle>=0) && (zenithCenter+angle<180) && (zenithCenter-angle>=0)) { //Case if the azimuth angle goes above 360 degrees
	if (((azimuthCenter+angle>=i) && (azimuthCenter-angle<=i) && (zenithCenter+angle>=j) && (zenithCenter-angle<=j)) || ((azimuthCenter-360+angle>=i) && (azimuthCenter-360+angle>=0) && (zenithCenter+angle>=j) && (zenithCenter-angle<=j))) {
	  productionDistributionJet[i][j]=1.0;
	  total+=1.0;
	}else{
	  productionDistributionJet[i][j]=0.0;
	} 
      }else if ((azimuthCenter+angle<360) && (azimuthCenter-angle<0) && (zenithCenter+angle<180) && (zenithCenter-angle>=0)) { //Case if the azimuth angle goes below 0 degrees
	if (((azimuthCenter+angle>=i) && (azimuthCenter-angle<=i) && (zenithCenter+angle>=j) && (zenithCenter-angle<=j)) || ((azimuthCenter+360-angle<=i) && (azimuthCenter+360-angle<360) && (zenithCenter+angle>=j) && (zenithCenter-angle<=j))) {
	  productionDistributionJet[i][j]=1.0;
	  total+=1.0;
	}else{
	  productionDistributionJet[i][j]=0.0;
	}
      }else{
	exit(__LINE__,__FILE__,"Zenith over boundaries not implemented");
      }
    }
  }

  if (total==0.0) exit(__LINE__,__FILE__,"Error: the jet has an angle of zero");

  cumulativeProductionDistributionJet[0][0]=0.0;
  for (i=0;i<360;i++) {
    for (j=0;j<180;j++) {
      if (i==0 && j==0) {
	cumulativeProductionDistributionJet[i][j]+=productionDistributionJet[i][j]/total;
      }else if (j==0) {
	cumulativeProductionDistributionJet[i][j]=cumulativeProductionDistributionJet[i-1][179]+productionDistributionJet[i][j]/total;
      }else{
	cumulativeProductionDistributionJet[i][j]=cumulativeProductionDistributionJet[i][j-1]+productionDistributionJet[i][j]/total;
      }
    }
  }
  probabilityFunctionDefinedJet=true;
  }
  //Computation of the area where the particle will be created
    gamma=rnd();
    i=0;
    j=0;
    while (gamma>cumulativeProductionDistributioJent[i][j]){
      if (j<179) {
	j++;
      }else{
	i++;
	j=0;
      }
    }


  //Computation of the angles where the particle will be created
    azimuth=(double) i;
    zenith= (double) j;

  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];

  ExternalNormal[0]=-sin(zenith*Pi/180)*cos(azimuth*Pi/180);
  ExternalNormal[1]=-sin(zenith*Pi/180)*sin(azimuth*Pi/180);
  ExternalNormal[2]=-cos(zenith*Pi/180);

  //SPHERE
  Sphere->GetSphereGeometricalParameters(x0Sphere,rSphere);
  
  x_LOCAL_IAU_OBJECT[0]=rSphere*sin(zenith*Pi/180)*cos(azimuth*Pi/180);
  x_LOCAL_IAU_OBJECT[1]=rSphere*sin(zenith*Pi/180)*sin(azimuth*Pi/180);
  x_LOCAL_IAU_OBJECT[2]=rSphere*cos(zenith*Pi/180);
      

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
  double zenithX,azimuthX,SubSolarAngle,x[3],r;

  x[0]=cos(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[0]-sin(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[1];
  x[1]=sin(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[0]+cos(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[1];
  x[2]=x_LOCAL_IAU_OBJECT[2];

  r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  zenithX=acos(x[2]/r);
  if (x[0]!=0) azimuthX=atan(x[1]/x[0]);
  if (x[0]==0) azimuthX=Pi/2;
  
  if (x[0]>=0) {
    SubSolarAngle=acos(fabs(sin(zenithX)*cos(azimuthX)));
    if (angle!=angle) SubSolarAngle=Pi/2;
  }

  SurfaceTemperature=GetSurfaceTemeprature(cos(SubSolarAngle),x_LOCAL_SO_OBJECT);
  // SurfaceTemperature=180.0;
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  
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
*/

double Comet::GetTotalProductionRateWaist(int spec,void *SphereDataPointer){
  //  return Comet::Jet_SourceRate[spec];
  double res;
  res=(spec==0)? 9.0e26:0.0;
  return res; //AT THIS STAGE, WE ARE ONLY TESTING THE FUNCTION GENERATEPARTICLEPROPERTIESBJORN BELOW
}
/*
bool Comet::GenerateParticlePropertiesWaist(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
  double ExternalNormal[3]; 
  int i,j;
  double total=0.0,TableTotalProductionRate,totalSurface,gamma,gamma2;
  double azimuth,zenith,azimuthTable[360],zenithTable[180],rSphere,*x0Sphere;

  if (probabilityFunctionDefinedWaist==false) {    
  for (i=0;i<360;i++) {
    for (j=0;j<180;j++) {
      if (((90>=i) && (40<=i)) || ((360-40>=i) && (360-90<=i))) {
	  productionDistributionJet[i][j]=1.0;
	  total+=1.0;
	}else{
	  productionDistributionJet[i][j]=0.0;
	}
    }
  }

  if (total==0.0) exit(__LINE__,__FILE__,"Error: the jet has an angle of zero");

  cumulativeProductionDistributionJet[0][0]=0.0;
  for (i=0;i<360;i++) {
    for (j=0;j<180;j++) {
      if (i==0 && j==0) {
	cumulativeProductionDistributionJet[i][j]+=productionDistributionJet[i][j]/total;
      }else if (j==0) {
	cumulativeProductionDistributionJet[i][j]=cumulativeProductionDistributionJet[i-1][179]+productionDistributionJet[i][j]/total;
      }else{
	cumulativeProductionDistributionJet[i][j]=cumulativeProductionDistributionJet[i][j-1]+productionDistributionJet[i][j]/total;
      }
    }
  }
  probabilityFunctionDefinedWaist=true;
  }
  //Computation of the area where the particle will be created
    gamma=rnd();
    i=0;
    j=0;
    while (gamma>cumulativeProductionDistributionJet[i][j]){
      if (j<179) {
	j++;
      }else{
	i++;
	j=0;
      }
    }


  //Computation of the angles where the particle will be created
    azimuth=(double) i;
    zenith= (double) j;

  //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];

  ExternalNormal[0]=-sin(zenith*Pi/180)*cos(azimuth*Pi/180);
  ExternalNormal[1]=-sin(zenith*Pi/180)*sin(azimuth*Pi/180);
  ExternalNormal[2]=-cos(zenith*Pi/180);

  Sphere->GetSphereGeometricalParameters(x0Sphere,rSphere);
  
  x_LOCAL_IAU_OBJECT[0]=rSphere*sin(zenith*Pi/180)*cos(azimuth*Pi/180);
  x_LOCAL_IAU_OBJECT[1]=rSphere*sin(zenith*Pi/180)*sin(azimuth*Pi/180);
  x_LOCAL_IAU_OBJECT[2]=rSphere*cos(zenith*Pi/180);
  
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

    double zenithX,azimuthX,SubSolarAngle,x[3],r;

  x[0]=cos(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[0]-sin(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[1];
  x[1]=sin(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[0]+cos(subSolarPointAzimuth)*x_LOCAL_IAU_OBJECT[1];
  x[2]=x_LOCAL_IAU_OBJECT[2];

  r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  zenithX=acos(x[2]/r);
  if (x[0]!=0) azimuthX=atan(x[1]/x[0]);
  if (x[0]==0) azimuthX=Pi/2;
  
  if (x[0]>=0) {
    SubSolarAngle=acos(fabs(sin(zenithX)*cos(azimuthX)));
    if (angle!=angle) SubSolarAngle=Pi/2;
  }

   SurfaceTemperature=GetSurfaceTemeprature(cos(SubSolarAngle),x_LOCAL_SO_OBJECT);
   //SurfaceTemperature=180.0;
  PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  
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
*/

void Comet::AntiSolarDirectionColumnMap::Print(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  FILE *fout=NULL;
  char fname[300];
  SpiceDouble xform[6][6],EarthState[6],lt;
  double xEarthLSO[3];
  SpiceDouble lGSE[6]={0,0,0,0,0,0},lLSO[6]={0,0,0,0,0,0};  //only 3 first components of the vectors are used. all 6 components are needed in the definition in order SPICE routines work correctly
  int idim;

  //calculate the currect position of the Earth
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,"LSO","none","Comet",EarthState,&lt);
  for (idim=0;idim<3;idim++) xEarthLSO[idim]=1.0E3*EarthState[idim];

  //calculate the rotation matrix from 'GSE' to 'LSO'
  sxform_c("GSE","LSO",Exosphere::OrbitalMotion::et,xform);

  //determine the number of angular points
  int nZenithPoints;
  double dZ,rr,ZenithAngle,AzimuthAngle,dZenithAngle;

  dZ=dZenithAngleMin;
  rr=(maxZenithAngle+dZenithAngleMax)/(maxZenithAngle+dZenithAngleMin);
  nZenithPoints=(long int)(log(dZenithAngleMax/dZenithAngleMin)/log(rr)-2.0);
  rr=pow(dZenithAngleMax/dZenithAngleMin,1.0/(nZenithPoints+2.0));

  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {
    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }

  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];
    char vlist[_MAX_STRING_LENGTH_PIC_]="";

    //open data file
    sprintf(fname,"%s/pic.Comet.Anti-sunwardColumnIntegrals.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s\"\n",utcstr);

//    fprintf(fout,"VARIABLES=\"Angle from the anti-solar direction [degree]\", \"Angle Out of Ecpliptic Plane [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");

    ColumnIntegral::GetVariableList(vlist);
    fprintf(fout,"VARIABLES=\"l[0]\", \"l[1]\" %s \n",vlist);


    fprintf(fout,"ZONE T=\"Column Density Map\"\n");
    fprintf(fout,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nAzimuthPoints+1,nZenithPoints);
    fprintf(fout,"DATAPACKING=POINT\n");
//    fprintf(fout,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
  }


  //calculate the integrals
  //calcualte the column integrals
  int StateVectorLength=ColumnIntegral::GetVariableList(NULL);
  double StateVector[StateVectorLength];


  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {

  for (int iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints+1;iAzimuthPoint++) {
    AzimuthAngle=2.0*Pi*double(iAzimuthPoint)/double(nAzimuthPoints);

      //get the pointing vector of integration in 'GSE' frame
      lGSE[0]=-cos(ZenithAngle);
      lGSE[1]=sin(ZenithAngle)*sin(AzimuthAngle);
      lGSE[2]=sin(ZenithAngle)*cos(AzimuthAngle);

      //convert the pointing vector from 'GSE' to 'LSO'
      mxvg_c(xform,lGSE,6,6,lLSO);

      //get the integrals
      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarthLSO,lLSO,ColumnIntegral::CoulumnDensityIntegrant);
      ColumnIntegral::ProcessColumnIntegrationVector(StateVector,StateVectorLength);

//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",ZenithAngle/Pi*180.0,AzimuthAngle/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);
//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",lGSE[1]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),lGSE[2]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),StateVector[0],StateVector[1],StateVector[2]);


      if (PIC::ThisThread==0) {
        fprintf(fout,"%e   %e",lGSE[1],lGSE[2]);

        for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",StateVector[i]);
        fprintf(fout,"\n");
      }

  }

    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }


  if (PIC::ThisThread==0) fclose(fout);
#endif
}
