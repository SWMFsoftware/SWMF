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

//parameters of the soruce processes
double Exosphere::SourceProcesses::ThermalDesorption::uThermal=1.85*eV2J;
double Exosphere::SourceProcesses::ThermalDesorption::VibrationalFrequency=1.0E13;

double Exosphere::SourceProcesses::PhotonStimulatedDesorption::PhotonFlux_1AU=2.0E14*1.0E4;  //Killen-2012-JGR, Yakshinskii+Madey-1999-?
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::CrossSection=3.0E-20*1.0E-4;  //Killen-2012-JGR, Yakshinskii+Madey-1999-?
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;

double Exosphere::SourceProcesses::ImpactVaporization::SourceRate=1.1e22;
double Exosphere::SourceProcesses::ImpactVaporization::HeliocentricDistance=1.0*_AU_;
double Exosphere::SourceProcesses::ImpactVaporization::SourceRatePowerIndex=0.0;
double Exosphere::SourceProcesses::ImpactVaporization::SourceTemeprature=2500.0;

double Exosphere::SourceProcesses::SolarWindSputtering::Yield=0.1;
double Exosphere::SourceProcesses::SolarWindSputtering::minInjectionEnergy=pow(10.0,2)*_NA__MASS_/2.0;
double Exosphere::SourceProcesses::SolarWindSputtering::maxInjectionEnergy=pow(10.0E3,2)*_NA__MASS_/2.0;


//typical parameters of solar wind
const double Exosphere::swVelocity_Typical[3]={-420.0E3,0.0,0.0};
const double Exosphere::swB_Typical[3]={-12.9E-9,4.71E-9,10.29E-9};
const double Exosphere::swTemperature_Typical=0.174e6,Exosphere::swNumberDensity_Typical=60.0E6;


//interaction with the surface
const double Exosphere::SurfaceInteraction::AccomodationCoefficient=0.2;

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

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) ExosphereEscapeRate[spec]=0.0,ExospehreTemsprature[spec]=1000.0;

  ExosphereEscapeRate[_NA_SPEC_]=Exosphere::SourceProcesses::ImpactVaporization::SourceRate;
  ExospehreTemsprature[_NA_SPEC_]=Exosphere::SourceProcesses::ImpactVaporization::SourceTemeprature;

  Exosphere::ChamberlainExosphere::Init(ExosphereEscapeRate,ExospehreTemsprature);



  //set up the model that collected the column integrals at the subsolar point of the limn as a function of the phase angle
  Sampling::SubsolarLimbColumnIntegrals::init();
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,Sampling::SubsolarLimbColumnIntegrals::CollectSample);

  //set up the model that prints the column integrals in the anti-solar direction
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SubsolarLimbColumnIntegrals::EmptyFunction,AntiSolarDirectionColumnMap::Print);

  //call init function of the exospheric model
  Exosphere::Init_AfterParser();
}

double Exosphere::SurfaceInteraction::SodiumStickingProbability(double& ReemissionParticleFraction,double Temp) {
  double res=0.0;

  ReemissionParticleFraction=0.5;

  return (Temp<200.0) ? 1.0 : 0.0;

//  return 1.0;


  static const double tmin=102.30103;
  static const double tmax=496.15054;
  static const double dt=110.66105-102.30103;

  static const int nPoints=49;

  struct cDataStruct {
    double t,s;
  };

  static const cDataStruct data[nPoints]={    //digitized from  Yakshinskiy-2005-SS
      {102.30103,0.98526},
      {110.66105,0.93862},
      {119.02106,0.90471},
      {127.38107,0.87079},
      {135.74109,0.83687},
      {143.17221,0.80861},
      {151.53223,0.77610},
      {159.89224,0.74501},
      {168.25226,0.71534},
      {176.61226,0.68707},
      {184.97227,0.65881},
      {192.40340,0.63479},
      {200.76341,0.60511},
      {209.12341,0.58250},
      {217.48343,0.55423},
      {225.84344,0.53304},
      {233.27457,0.51325},
      {241.63458,0.49629},
      {249.99460,0.47933},
      {258.35461,0.46238},
      {266.71460,0.44542},
      {275.07465,0.42846},
      {282.50577,0.41150},
      {290.86578,0.39454},
      {299.22580,0.37617},
      {307.58582,0.36063},
      {315.94583,0.34649},
      {323.37695,0.33802},
      {331.73697,0.32106},
      {340.09695,0.30693},
      {348.45700,0.29703},
      {356.81699,0.28714},
      {365.17703,0.27584},
      {372.60815,0.26736},
      {380.96814,0.25888},
      {389.32816,0.25181},
      {397.68817,0.24333},
      {406.04819,0.23344},
      {413.47931,0.23061},
      {421.83932,0.22355},
      {430.19934,0.21931},
      {438.55933,0.21224},
      {446.91937,0.20942},
      {455.27936,0.20376},
      {462.71048,0.19952},
      {471.07053,0.19670},
      {479.43051,0.19246},
      {487.79056,0.18963},
      {496.15054,0.18539}};


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
void Exosphere::SodiumCoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd;
  double NumberDensity;

  for (i=0;i<resLength;i++) res[i]=0.0;

  //get the local density number
  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(_NA_SPEC_);
  res[0]=NumberDensity;

  //get the local scattering
  if (resLength>1) {
    double BulkVelocity_SO[3],v_LOCAL_SO_FROZEN[3],rHeliocentric,vHeliocentric;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity_SO,_NA_SPEC_);

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


    if ( ((x[1]*x[1]+x[2]*x[2]<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))&&(x[0]<0.0)) || (Moon::EarthShadowCheck(x)==true) ) res[1]=0.0,res[2]=0.0;
    else {
      res[1]=1.0E-10*NumberDensity*SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      res[2]=1.0E-10*NumberDensity*SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
    }
  }
}


//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble EphemerisTime) {
  SpiceDouble State[6],ltlocal;
  double EccentricityVector[3];
  double res,Speed2,a,c,absEccentricity;
  const double GravitationalParameter=GravityConstant*_MASS_(_EARTH_);
  double vMoon[3],xMoon[3],rGeocentric=0.0;
  int idim;

  spkezr_c("Moon",EphemerisTime,"MSGR_HCI","none","Earth",State,&ltlocal);

  for (idim=0;idim<3;idim++) {
    xMoon[idim]=State[idim]*1.0E3;
    vMoon[idim]=State[idim+3]*1.0E3;

    rGeocentric+=pow(xMoon[idim],2);
  }

  rGeocentric=sqrt(rGeocentric);
  Speed2=vMoon[0]*vMoon[0]+vMoon[1]*vMoon[1]+vMoon[2]*vMoon[2];
  c=xMoon[0]*vMoon[0]+xMoon[1]*vMoon[1]+xMoon[2]*vMoon[2];

  for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
    EccentricityVector[idim]=Speed2/GravitationalParameter*xMoon[idim] - c/GravitationalParameter*vMoon[idim] - xMoon[idim]/rGeocentric;
    absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
    a+=EccentricityVector[idim]*xMoon[idim];
  }

  absEccentricity=sqrt(absEccentricity);
  res=acos(a/(absEccentricity*rGeocentric));

  if (c<0.0) res=2.0*Pi-res;

  return res;
}


void Moon::AntiSolarDirectionColumnMap::Print(int DataOutputFileNumber) {
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

    //open data file
    sprintf(fname,"pic.Moon.Anti-sunwardColumnIntegrals.out=%i.dat",DataOutputFileNumber);

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s\"\n",utcstr);

//    fprintf(fout,"VARIABLES=\"Angle from the anti-solar direction [degree]\", \"Angle Out of Ecpliptic Plane [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");

    fprintf(fout,"VARIABLES=\"l[0]\", \"l[1]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");


    fprintf(fout,"ZONE T=\"Column Density Map\"\n");
    fprintf(fout,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nAzimuthPoints+1,nZenithPoints);
    fprintf(fout,"DATAPACKING=POINT\n");
    fprintf(fout,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
  }


  //calculate the integrals
  //calcualte the column integrals
  const int StateVectorLength=3;
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
      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarthLSO,lLSO,SodiumCoulumnDensityIntegrant);

//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",ZenithAngle/Pi*180.0,AzimuthAngle/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);
//      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",lGSE[1]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),lGSE[2]/sqrt(lGSE[1]*lGSE[1]+lGSE[2]*lGSE[2]),StateVector[0],StateVector[1],StateVector[2]);
      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",lGSE[1],lGSE[2],StateVector[0],StateVector[1],StateVector[2]);

  }

    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }


  if (PIC::ThisThread==0) fclose(fout);
}
