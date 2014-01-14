//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Moon";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MOON";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="LSO";


int Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=-1;
int Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=-1;
int Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=-1;



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

void Moon::Init_AfterParser() {

  //set up the Chamberlen model
  double ExosphereEscapeRate[PIC::nTotalSpecies],ExospehreTemsprature[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) { //ExosphereEscapeRate[spec]=0.0,ExospehreTemsprature[spec]=1000.0;
    ExosphereEscapeRate[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceRate[spec];
    ExospehreTemsprature[spec]=Exosphere::SourceProcesses::ImpactVaporization::ImpactVaporization_SourceTemeprature[spec];
  }

  Exosphere::ChamberlainExosphere::Init(ExosphereEscapeRate,ExospehreTemsprature);



  //set up the model that collected the column integrals at the subsolar point of the limn as a function of the phase angle
  Sampling::SubsolarLimbColumnIntegrals::init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Sampling::SubsolarLimbColumnIntegrals::CollectSample);

  //set up the model that prints the column integrals in the anti-solar direction
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,AntiSolarDirectionColumnMap::Print);

  //set up sampling of velocity distribution functions
  Moon::Sampling::VelocityDistribution::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Moon::Sampling::VelocityDistribution::Sampling,Moon::Sampling::VelocityDistribution::OutputSampledData);

  //call init function of the exospheric model
  Exosphere::Init_AfterParser();

  //init the model of calcualting the integrals that correspond to the Kaguya's TVIS observations
  Moon::Sampling::Kaguya::Init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Moon::Sampling::Kaguya::TVIS::OutputModelData);
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
   case _NA_SPEC_: case _NAPLUS_SPEC_:
     res=SodiumStickingProbability(ReemissionParticleFraction,Temp);
     break;
   default:
     exit(__LINE__,__FILE__,"the option is not implemented");
   }

   return res;
}


//surface temeprature of the planet
double Exosphere::GetSurfaceTemeprature(double CosSubSolarAngle,double *x_LOCAL_SO_OBJECT) {

  //determine if the point on the night side of the Moon
  if (CosSubSolarAngle<0.0) return 100.0;

  //determine if the point is within the shadow of the Earth
  if (Moon::EarthShadowCheck(x_LOCAL_SO_OBJECT)==true) return 100.0;

  //return the day-side temeprature
  return 280*pow(CosSubSolarAngle,0.25)+100.0;
}

//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\",  \"Mean Speed Along the Line of Sight(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec),PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=2;
  }

  if (_NA_SPEC_>=0) Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_=2*_NA_SPEC_;

  //sodium emission
  if (_NA_SPEC_>=0) {
    Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_=nVariables;
    Moon::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_=nVariables+1;

    if (vlist!=NULL) sprintf(vlist,"%s,  \"Sodium Emission(5891_58A)\",  \"Sodium Emission(5897_56A)\"",vlist);
    nVariables+=2;
  }

  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
  int spec,cnt=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (res[cnt]>0.0) res[cnt+1]/=res[cnt];
    cnt+=2;
  }
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
      //check if the point is outside of the Moon's and Earth's shadows
      if ( (Moon::EarthShadowCheck(x)==false) && ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_MOON_)*_RADIUS_(_MOON_))) ) {
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        res[cnt++]=1.0E-10*node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_)*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
      else res[cnt++]=0.0,res[cnt++]=0.0;
    }
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  return GetTAA("Moon","Earth",_MASS_(_EARTH_),EphemerisTime);
}


void Moon::AntiSolarDirectionColumnMap::Print(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  FILE *fout=NULL;
  char fname[300];
  SpiceDouble xform[6][6],EarthState[6],lt;
  double xEarthLSO[3];
  SpiceDouble lGSE[6]={0,0,0,0,0,0},lLSO[6]={0,0,0,0,0,0};  //only 3 first components of the vectors are used. all 6 components are needed in the definition in order SPICE routines work correctly
  int idim;

  //calculate the currect position of the Earth
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,"LSO","none","Moon",EarthState,&lt);
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
    sprintf(fname,"%s/pic.Moon.Anti-sunwardColumnIntegrals.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);

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
