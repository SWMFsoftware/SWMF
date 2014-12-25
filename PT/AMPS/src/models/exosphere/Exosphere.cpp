//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury.cpp
 *
 *  Created on: Feb 13, 2012
 *      Author: vtenishe
 */

//$Id$


#include "pic.h"
//#include "SpiceUsr.h"

#include "Exosphere.h"
#include "constants.h"
#include "SingleVariableDiscreteDistribution.h"

//#include "Rosetta.h"


double Exosphere::swE_Typical[3]={0.0,0.0,0.0};
double Exosphere::xObject_HCI[3]={0.0,0.0,0.0},Exosphere::vObject_HCI[3]={0.0,0.0,0.0};
double Exosphere::xEarth_HCI[3]={0.0,0.0,0.0},Exosphere::vEarth_HCI[3]={0.0,0.0,0.0};
double Exosphere::xEarth_SO[3]={0.0,0.0,0.0},Exosphere::vEarth_SO[3]={0.0,0.0,0.0};
double Exosphere::xSun_SO[3]={0.0,0.0,0.0},Exosphere::vSun_SO[3]={0.0,0.0,0.0};
double Exosphere::vObjectRadial=0.0,Exosphere::xObjectRadial=0.0;

double Exosphere::vObject_SO_FROZEN[3]={0.0,0.0,0.0};
double Exosphere::RotationVector_SO_FROZEN[3],Exosphere::RotationRate_SO_FROZEN; //the vectors of the direction of rotation and the rotation rate of the SO in SO_FROZEN

//matrix for transformation SO->HCI coordinate frame (used only when prepare data files)
SpiceDouble Exosphere::Sampling::OutputDataFile::SO_to_HCI_TransformationMartix[6][6];

//sample the source rate
double Exosphere::Sampling::CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];

//sampling offsets
int Exosphere::Sampling::SamplingDensityOffset[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
int Exosphere::Sampling::CellSamplingDataOffset=-1;
double **Exosphere::Sampling::PlanetNightSideReturnFlux=NULL;
double *Exosphere::Sampling::TotalPlanetReturnFlux=NULL,*Exosphere::Sampling::PlanetSurfaceStickingRate=NULL;

//the field in the particle's data that keeps the id of the source process due to which the particle has beed produced
long int Exosphere::Sampling::ParticleData_SourceProcessID_Offset=-1;
long int Exosphere::Sampling::ParticleData_OriginSurfaceElementNumber_Offset=-1;

//the total number of source processes
int Exosphere::nTotalSourceProcesses=0;

//the sphere that representd the planet
cInternalSphericalData *Exosphere::Planet=NULL;

//the total source rate values for specific source processes
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::SourceRate[PIC::nTotalSpecies],Exosphere::SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate[PIC::nTotalSpecies];
double Exosphere::SourceProcesses::ThermalDesorption::SourceRate[PIC::nTotalSpecies],Exosphere::SourceProcesses::ThermalDesorption::maxLocalSourceRate[PIC::nTotalSpecies];
double Exosphere::SourceProcesses::SolarWindSputtering::SourceRate[PIC::nTotalSpecies],Exosphere::SourceProcesses::SolarWindSputtering::maxLocalSourceRate[PIC::nTotalSpecies];

//evaluate nemerically the source rate
/*
double Exosphere::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Exosphere::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate=0.0;
double Exosphere::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Exosphere::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate=0.0;
*/

//user-defined additional output data
Exosphere::Sampling::fUserDefinedAdditionalData_VariableList_OutputSampledModelData Exosphere::Sampling::UserDefinedAdditionalData_VariableList_OutputSampledModelData=NULL;
Exosphere::Sampling::fUserDefinedAdditionalData_OutputSampledModelData Exosphere::Sampling::UserDefinedAdditionalData_OutputSampledModelData=NULL;

//request the sampling and particle's field data
int Exosphere::Sampling::RequestSamplingData(int offset) {
  int SamplingLength=0;

  if (CellSamplingDataOffset!=-1) exit(__LINE__,__FILE__,"Error: second request for the sampling data");

  CellSamplingDataOffset=offset;
  for (int iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) SamplingDensityOffset[iSource]=-1;

#if _EXOSPHERE_SOURCE__EXTERNAL_BOUNDARY_INJECTION_ == _EXOSPHERE_SOURCE__ON_
  SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  Exosphere::Sampling::SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
  SamplingDensityOffset[_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_]=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

  //reserve the space for sampling particles's density generated by user defined source functions
#if _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_ == _EXOSPHERE_SOURCE__ON_
$MARKER:RESERVE-CELL-SAMPLING-DATA-BUFFER$
#endif


  //reserve the particle fields
  PIC::ParticleBuffer::RequestDataStorage(ParticleData_SourceProcessID_Offset,2*sizeof(int));
  ParticleData_OriginSurfaceElementNumber_Offset=ParticleData_SourceProcessID_Offset+sizeof(int);


  return SamplingLength;
}

//init the model
void Exosphere::Init_BeforeParser() {
  int idim;

  //set the model of inhectino of the ion at the boundary of the domain due to the background plasma evnoronment
  #if _EXOSPHERE__BACKGROUND_PLASMA_ION_INJECTION_ == _PIC_MODE_ON
  PIC::ParticleWeightTimeStep::ExosphereModelExtraSourceRate=SourceProcesses::BackgroundPlasmaBoundaryIonInjection::GetTotalProductionRate;
  #endif

  //calculate the typical value of the motional electrical field
  swE_Typical[0]=-(/*Exosphere_*/swVelocity_Typical[1]*/*Exosphere_*/swB_Typical[2]-/*Exosphere_*/swVelocity_Typical[2]*/*Exosphere_*/swB_Typical[1]);
  swE_Typical[1]=+(/*Exosphere_*/swVelocity_Typical[0]*/*Exosphere_*/swB_Typical[2]-/*Exosphere_*/swVelocity_Typical[2]*/*Exosphere_*/swB_Typical[0]);
  swE_Typical[2]=-(/*Exosphere_*/swVelocity_Typical[0]*/*Exosphere_*/swB_Typical[1]-/*Exosphere_*/swVelocity_Typical[1]*/*Exosphere_*/swB_Typical[0]);


  //set the pre-processor into the ICES model
  PIC::CPLR::ICES::SWMFdataPreProcessor=SWMFdataPreProcessor;

  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SampleModelData,Sampling::OutputSampledModelData);

  //init the buffer for sampling of the source rate
  for (int s=0;s<PIC::nTotalSpecies;s++) for (int i=0;i<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;i++)  Sampling::CalculatedSourceRate[s][i]=0.0;


  //furnish the SPICE kernels
#if  _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  char str[_MAX_STRING_LENGTH_PIC_];

  for (int nKernel=0;nKernel<nFurnishedSPICEkernels;nKernel++) {
    struct stat buf;

    sprintf(str,"%s/%s",SPICE_Kernels_PATH,SPICE_Kernels[nKernel]);

    if (stat(str, &buf) != 0) {
      char f[_MAX_STRING_LENGTH_PIC_];

      sprintf(f,"SPICE kernel %s is not found",str);
      exit(__LINE__,__FILE__,f);
    }

    furnsh_c(str);
  }
#endif

  //Get the initial parameters of Mercury orbit
  SpiceDouble state[6];

  utc2et_c(SimulationStartTimeString,&OrbitalMotion::et);

  //get initial parameters of Mercury's orbit
  spkezr_c(ObjectName,OrbitalMotion::et,"J2000","none","SUN",state,&OrbitalMotion::lt);

  for (idim=0,xObjectRadial=0.0;idim<3;idim++) {
    xObject_HCI[idim]=state[idim]*1.0E3,vObject_HCI[idim]=state[idim+3]*1.0E3;
    xObjectRadial+=pow(xObject_HCI[idim],2);
  }

  xObjectRadial=sqrt(xObjectRadial);
  vObjectRadial=(xObject_HCI[0]*vObject_HCI[0]+xObject_HCI[1]*vObject_HCI[1]+xObject_HCI[2]*vObject_HCI[2])/xObjectRadial;

  //get initial position of Earth
  spkezr_c("Earth",OrbitalMotion::et,"J2000","none","SUN",state,&OrbitalMotion::lt);
  for (idim=0;idim<3;idim++) xEarth_HCI[idim]=state[idim]*1.0E3,vEarth_HCI[idim]=state[idim+3]*1.0E3;



  //calcualte the total number of source processes
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__OFF_
//do nothing
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__OFF_
//do nothing
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestSamplingData.push_back(Sampling::RequestSamplingData);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(Sampling::OutputDataFile::PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(Sampling::OutputDataFile::PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Sampling::OutputDataFile::Interpolate);

  //init the energy distribution function for specific injection processes
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) if (SourceProcesses::PhotonStimulatedDesorption::PhotonStimulatedDesorption_CrossSection[spec]>0.0) {
    SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].Init(SourceProcesses::PhotonStimulatedDesorption::PhotonStimulatedDesorption_minInjectionEnergy[spec], \
        SourceProcesses::PhotonStimulatedDesorption::PhotonStimulatedDesorption_maxInjectionEnergy[spec],SourceProcesses::PhotonStimulatedDesorption::EnergyDistributionFunction, \
        _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_,&spec);


/*    if (PIC::ThisThread==0) {
      char fname[200];

      sprintf(fname,"CumulativeEnergyDistribution-PSD.nspec=%i.%s.dat",spec,PIC::MolecularData::GetChemSymbol(spec));
      SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintCumulativeDistributionFunction(fname);

      sprintf(fname,"EnergyDistribution-PSD.nspec=%i.%s.dat",spec,PIC::MolecularData::GetChemSymbol(spec));
      SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintDistributionFunction(fname,&spec);
    }*/
  }
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_

#if _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION_ == _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION__NUMERIC_
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) if (SourceProcesses::SolarWindSputtering::SolarWindSputtering_Yield[spec]>0.0) {
    SourceProcesses::SolarWindSputtering::EnergyDistribution[spec].Init(SourceProcesses::SolarWindSputtering::SolarWindSputtering_minInjectionEnergy[spec], \
        SourceProcesses::SolarWindSputtering::SolarWindSputtering_maxInjectionEnergy[spec],SourceProcesses::SolarWindSputtering::EnergyDistributionFunction, \
        _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_,&spec);

/*
    if (PIC::ThisThread==0) {
      char fname[200];

      sprintf(fname,"CumulativeEnergyDistribution-SWS.nspec=%i.%s.dat",spec,PIC::MolecularData::GetChemSymbol(spec));
      SourceProcesses::SolarWindSputtering::EnergyDistribution[spec].fPrintCumulativeDistributionFunction(fname);

      sprintf(fname,"EnergyDistribution-SWS.nspec=%i.%s.dat",spec,PIC::MolecularData::GetChemSymbol(spec));
      SourceProcesses::SolarWindSputtering::EnergyDistribution[spec].fPrintDistributionFunction(fname,&spec);
    }
*/
  }
#elif _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION_ == _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION__USER_DEFINED_
  // do nothing
#else
  exit(__LINE__, __FILE__, "ERROR: _EXOSPHERE__ENERGY_DISTRIBUTION_INVERSION_ is not defined")
#endif

#endif


  //init the source rate buffers
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    SourceProcesses::PhotonStimulatedDesorption::SourceRate[spec]=0.0,SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate[spec]=0.0;
    SourceProcesses::ThermalDesorption::SourceRate[spec]=0.0,SourceProcesses::ThermalDesorption::maxLocalSourceRate[spec]=0.0;
    SourceProcesses::SolarWindSputtering::SourceRate[spec]=0.0,SourceProcesses::SolarWindSputtering::maxLocalSourceRate[spec]=0.0;
  }
}

void Exosphere::Init_AfterParser() {

  //init the sampling buffer of the returned flux on the night side of the planet
  int offset,s,i;

  Sampling::PlanetNightSideReturnFlux=new double *[PIC::nTotalSpecies];
  Sampling::PlanetNightSideReturnFlux[0]=new double [PIC::nTotalSpecies*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];

  Sampling::TotalPlanetReturnFlux=new double[PIC::nTotalSpecies];
  Sampling::PlanetSurfaceStickingRate=new double[PIC::nTotalSpecies];

  for (s=0,offset=0;s<PIC::nTotalSpecies;s++) {
    Sampling::PlanetNightSideReturnFlux[s]=Sampling::PlanetNightSideReturnFlux[0]+offset;
    offset+=_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;

    for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) Sampling::PlanetNightSideReturnFlux[s][i]=0.0;

    Sampling::TotalPlanetReturnFlux[s]=0.0,Sampling::PlanetSurfaceStickingRate[s]=0.0;
  }

  //remove old and create header for the new the file that contains the orbital information of the planet's motion

  if (PIC::ThisThread==0) {
    char OrbitalDataFileName[_MAX_STRING_LENGTH_PIC_];

    sprintf(OrbitalDataFileName,"rm -f %s/pic.OrbitalData.dat",PIC::OutputDataFileDirectory);
    system(OrbitalDataFileName);

    sprintf(OrbitalDataFileName,"%s/pic.OrbitalData.dat",PIC::OutputDataFileDirectory);
    FILE *fout=fopen(OrbitalDataFileName,"w");
    fprintf(fout,"Variables: \"UTC\", \
         \"xObject(J2000)\", \"yObject(J2000)\", \"zObject(J2000)\", \
         \"xEarth(J2000)\", \"yEarth(J2000)\", \"zEarth(J2000)\", \
         \"xSun(SO)\", \"ySun(SO)\", \"zSun(SO)\", \
         \"xEarth(SO)\", \"yEarth(SO)\", \"zEarth(SO)\", \
         \"UTC of the corresponding data file\", \"Output data file number\"");

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) fprintf(fout,", \"IAU2SO[%i][%i]\"",i,j);

    fprintf(fout,"\n");
    fclose(fout);
  }


  //init the model that described the available remote ground based observations
  Exosphere::Sampling::ReferenceGroundBasedObservations::init();

  //init the source models
  Exosphere::SourceProcesses::Init();
}

/*void Exosphere::Init_AfterMesh() {
  //init the model of the background plasma ion injection
  #if _EXOSPHERE__BACKGROUND_PLASMA_ION_INJECTION_ == _PIC_MODE_ON_
  SourceProcesses::BackgroundPlasmaBoundaryIonInjection::Init();
  #endif
}*/



//ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
void Exosphere::SWMFdataPreProcessor(double *x,PIC::CPLR::ICES::cDataNodeSWMF& data) {
  int i;

  if (data.status!=_PIC_ICES__STATUS_OK_) {
    for (i=0;i<3;i++) {
      data.B[i]=/*Exosphere_*/swB_Typical[i];
      data.E[i]=swE_Typical[i];
      data.swVel[i]=/*Exosphere_*/swVelocity_Typical[i];
    }

    data.swTemperature=/*Exosphere_*/swTemperature_Typical;
    data.swNumberDensity=/*Exosphere_*/swNumberDensity_Typical;

    //p=2*n*k*T assuming quasi neutrality ni=ne=n and Ti=Te=T
    data.swPressure=2.0*data.swNumberDensity*Kbol*data.swTemperature;
  }
}

//calculate the sodium column density and plot
void Exosphere::ColumnIntegral::Tail(char *fname) {
  double xEarth_new[3]={0.0,0.0,0.0};

  const int nPoints=500;
  const double IntegrationRangeBegin=-50.0E6;
  const double IntegrationRangeEnd=50.0E6;
  const double IntegrationStep=(IntegrationRangeEnd-IntegrationRangeBegin)/(nPoints-1);


  //find position of Earth
  SpiceDouble State[6],lt;

  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,State,&lt);

  xEarth_new[0]=State[0]*1.0E3;
  xEarth_new[1]=State[1]*1.0E3;
  xEarth_new[2]=State[2]*1.0E3;

  //open the output file
  FILE *fout=NULL;
  int npoint;

  if (PIC::ThisThread==0) {
    fout=fopen(fname,"w");
    char vlist[_MAX_STRING_LENGTH_PIC_]="";

    ColumnIntegral::GetVariableList(vlist);
    fprintf(fout,"VARIABLES=\"Distance from the planet\" %s\n",vlist);
  }

  for (npoint=0;npoint<nPoints;npoint++) {
     double l[3];
     int StateVectorLength=ColumnIntegral::GetVariableList(NULL);
     double StateVector[StateVectorLength];

     l[0]=-(IntegrationRangeBegin+npoint*IntegrationStep)-xEarth_new[0];
     l[1]=-xEarth_new[1];
     l[2]=-xEarth_new[2];

     PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarth_new,l,ColumnIntegral::CoulumnDensityIntegrant);
     ColumnIntegral::ProcessColumnIntegrationVector(StateVector,StateVectorLength);

     if (PIC::ThisThread==0) {
       fprintf(fout,"%e  ",IntegrationRangeBegin+npoint*IntegrationStep);

       for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",StateVector[i]);
       fprintf(fout,"\n");
     }
  }

  if (PIC::ThisThread==0) fclose(fout);
}


void Exosphere::ColumnIntegral::GetSubsolarPointDirection(double *LimbDirection,double *EarthPosition) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  SpiceInt n,nxpts;
  SpiceDouble rad[3],xpt0[3],xpt1[3],xLimb_IAU[3],xLimb_SO[3];
  SpiceDouble EarthState_IAU[6],EarthState_SO[6],lt,SunState_IAU[6];
  SpiceEllipse limb;
  SpicePlane plane;

  //find the limb
  bodvrd_c(ObjectName, "RADII", 3, &n, rad );
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,IAU_FRAME,"none",ObjectName,EarthState_IAU,&lt);
  edlimb_c(rad[0],rad[1],rad[2],EarthState_IAU,&limb);

  //find intersection of the limb with the plane hat contains position of the Earth, center of the boly and the Sun
  spkezr_c("SUN",Exosphere::OrbitalMotion::et,IAU_FRAME,"none","EARTH",SunState_IAU,&lt);
  psv2pl_c (EarthState_IAU,EarthState_IAU,SunState_IAU,&plane);
  inelpl_c (&limb,&plane,&nxpts,xpt0,xpt1);

  //choose the limb
  double cosXpt0=0.0,cosXpt1=0.0;
  int idim=0;

  if (nxpts<2) {
    exit(__LINE__,__FILE__,"Error: cannot file point of intersection of the plane that contains centers of the Earth, the Object and the Sun with the elips of the limb");
  }

  for (idim=0;idim<3;idim++) {
    cosXpt0+=(xpt0[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
    cosXpt1+=(xpt1[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
  }

  if (cosXpt0>cosXpt1) memcpy(xLimb_IAU,xpt0,3*sizeof(SpiceDouble)); //'xpt0' is closer to the Sun
  else memcpy(xLimb_IAU,xpt1,3*sizeof(SpiceDouble)); //'xpt1' is closer to the Sun


  //recalcualte position of the Earth in the SO_FRAME
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,EarthState_SO,&lt);

  //get the column integral
  double l[3],c=0.0;

  for (idim=0;idim<3;idim++) {
    //recalculate position of the limb in the 'SO_FRAME'
    xLimb_SO[idim]=
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][0]*xLimb_IAU[0])+
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][1]*xLimb_IAU[1])+
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][2]*xLimb_IAU[2]);


    //gwt the pointing vector
    l[idim]=xLimb_SO[idim]-EarthState_SO[idim];
    c+=pow(l[idim],2);

    EarthPosition[idim]=1.0E3*EarthState_SO[idim];
  }

  for (c=sqrt(c),idim=0;idim<3;idim++) LimbDirection[idim]=l[idim]/c;
#else
  for (int idim=0;idim<3;idim++) LimbDirection[idim]=0.0;
#endif
}


void Exosphere::ColumnIntegral::Limb(char *fname) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  FILE *fout=NULL,*fLimb=NULL;

  const double maxAltitude=50.0; //altititude in radii of the object
  const double dAltmin=0.001,dAltmax=0.01*maxAltitude; //the minimum and maximum resolution along the 'altitude' line

  //the 'altitude' points are distributed logarithmicaly
  //determin the factor 'rr' of the geometrical progression that defined the location od radial points of the map
  long int t;
  double dAlt,rr,R=dAltmin;

  dAlt=dAltmin;
  rr=(maxAltitude+dAltmax)/(maxAltitude+dAltmin);
  t=(long int)(log(dAltmax/dAltmin)/log(rr)-2.0);
  rr=pow(dAltmax/dAltmin,1.0/(t+2.0));


  //determine the position of the limb
  SpiceInt n,nxpts;
  SpiceDouble rad[3],xpt0[3],xpt1[3],xLimb_IAU[3],xLimb_SO[3];
  SpiceDouble EarthState_IAU[6],EarthState_SO[6],lt,SunState_IAU[6];
  SpiceEllipse limb;
  SpicePlane plane;

  //find the limb
  bodvrd_c(ObjectName, "RADII", 3, &n, rad );
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,IAU_FRAME,"none",ObjectName,EarthState_IAU,&lt);
  edlimb_c(rad[0],rad[1],rad[2],EarthState_IAU,&limb);

  //find intersection of the limb with the plane hat contains position of the Earth, center of the boly and the Sun
  spkezr_c("SUN",Exosphere::OrbitalMotion::et,IAU_FRAME,"none","EARTH",SunState_IAU,&lt);
  psv2pl_c (EarthState_IAU,EarthState_IAU,SunState_IAU,&plane);
  inelpl_c (&limb,&plane,&nxpts,xpt0,xpt1);

  //choose the limb
  double cosXpt0=0.0,cosXpt1=0.0;
  int idim=0;

  if (nxpts<2) {
    exit(__LINE__,__FILE__,"Error: cannot file point of intersection of the plane that contains centers of the Earth, the Object and the Sun with the elips of the limb");
  }

  for (idim=0;idim<3;idim++) {
    cosXpt0+=(xpt0[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
    cosXpt1+=(xpt1[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
  }

  if (cosXpt0>cosXpt1) memcpy(xLimb_IAU,xpt0,3*sizeof(SpiceDouble)); //'xpt0' is closer to the Sun
  else memcpy(xLimb_IAU,xpt1,3*sizeof(SpiceDouble)); //'xpt1' is closer to the Sun


  //recalcualte position of the Earth in the SO_FRAME
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,EarthState_SO,&lt);

  //get the column integral
  double l[3],rEarth[3],c=0.0;

  for (idim=0;idim<3;idim++) {
    //recalculate position of the limb in the 'SO_FRAME'
    xLimb_SO[idim]=
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][0]*xLimb_IAU[0])+
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][1]*xLimb_IAU[1])+
        (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][2]*xLimb_IAU[2]);

    //Convert position of the limb in meters
    xLimb_SO[idim]*=1.0E3;
    rEarth[idim]=1.0E3*EarthState_SO[idim];
  }


  //open the output file
  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];
    char vlist[_MAX_STRING_LENGTH_PIC_]="";

    ColumnIntegral::GetVariableList(vlist);

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s, Radial size of the planet=%e\"\n",utcstr,atan(_RADIUS_(_TARGET_)/sqrt(rEarth[0]*rEarth[0]+rEarth[1]*rEarth[1]+rEarth[2]*rEarth[2]))/Pi*180.0);

    fprintf(fout,"VARIABLES=\"R [Object Radii]\", \"R [m]\", \"Angle from the center of the object [degree]\" %s \n",vlist);
    fprintf(fout,"ZONE T=\"Column Density Along the Limb\"\n");


    //output emmision and column integrals at the limb at different phase angles
    char LimbIntegralsFileName[_MAX_STRING_LENGTH_PIC_];

    sprintf(LimbIntegralsFileName,"%s/pic.LimbIntegrals.dat",PIC::OutputDataFileDirectory);

    if (PIC::DataOutputFileNumber==0) {
      fLimb=fopen(LimbIntegralsFileName,"w");
      fprintf(fLimb,"VARIABLES=\"Phase Angle[degrees]\", \"Julian Date\", \"Object Radial Velocity\"  %s \n",vlist);
    }
    else fLimb=fopen(LimbIntegralsFileName,"a");
  }



  //calculate the integrals
  int StateVectorLength=ColumnIntegral::GetVariableList(NULL);
  double StateVector[StateVectorLength];
  double StateVectorMax[StateVectorLength]; //maximum values of the state variables == the values of the state variables at the limb

  for (int i=0;i<StateVectorLength;i++) StateVectorMax[i]=0.0;

  while (R<maxAltitude) {

    //get the pointing direction
    for (c=0.0,idim=0;idim<3;idim++) {
      l[idim]=(1.0+R)*xLimb_SO[idim]-rEarth[idim];
      c+=pow(l[idim],2);
    }

    for (c=sqrt(c),idim=0;idim<3;idim++) l[idim]/=c;

    PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,rEarth,l,ColumnIntegral::CoulumnDensityIntegrant);
    ColumnIntegral::ProcessColumnIntegrationVector(StateVector,StateVectorLength);

    for (int i=0;i<StateVectorLength;i++) if (StateVectorMax[i]<StateVector[i]) StateVectorMax[i]=StateVector[i];

    if (PIC::ThisThread==0) {
      fprintf(fout,"%e   %e   %e",R,R*_RADIUS_(_TARGET_),atan((1.0+R)*_RADIUS_(_TARGET_)/sqrt(rEarth[0]*rEarth[0]+rEarth[1]*rEarth[1]+rEarth[2]*rEarth[2]))/Pi*180.0);
      for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",StateVector[i]);
      fprintf(fout,"\n");
    }

    R+=dAlt;
    dAlt*=rr;
  }

  if (PIC::ThisThread==0) {
    fclose(fout);

    //output the time-dependent values of the integrals at the limb
    SpiceDouble StateMoonEarth_SO[6],StateSun_SO[6];
    static SpiceDouble lastStateMoonEarth_SO[6]={0.0,0.0,0.0,0.0,0.0,0.0};
    static int cnt=0;

    spkezr_c("Moon",Exosphere::OrbitalMotion::et,SO_FRAME,"none","Earth",StateMoonEarth_SO,&lt);
    spkezr_c("Sun",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,StateSun_SO,&lt);

    if (PIC::DataOutputFileNumber!=0) if (lastStateMoonEarth_SO[3]*StateMoonEarth_SO[3]<0.0) {
      //the radial velocity of the object has changed its direction
      fprintf(fLimb,"ZONE T=\"Pass=%i\"\n",cnt++);
    }

    memcpy(lastStateMoonEarth_SO,StateMoonEarth_SO,6*sizeof(SpiceDouble));

    //VARIABLES=\"Phase Angle[degrees]\", \"Julian Date\", \"Object Radial Velocity\", \"Column Density [m^{-2}]\", \"g-factor (5891.58A)\", \"Intensity (5891.58A) [R]\", \"g-factor (5897.56A)\", \"Intensity (5897.56A) [R]\" \n");

    fprintf(fLimb,"%e  %e  %e",Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et)/Pi*180.0,j2000_c()+Exosphere::OrbitalMotion::et/spd_c(),1.0E3*StateSun_SO[3]);
    for (int i=0;i<StateVectorLength;i++) fprintf(fLimb,"   %e",StateVectorMax[i]);
    fprintf(fLimb,"\n");
    fclose(fLimb);
  }
#endif
}


/*
void Exosphere::ColumnDensityIntegration_Limb(char *fname) {
  FILE *fout=NULL,*fLimb=NULL;
  int idim;
  double rEarth[3]={0.0,0.0,0.0},l[3]={0.0,0.0,0.0},rSun[3];

  const double maxAltitude=50.0; //altititude in radii of the object
  const double dAltmin=0.001,dAltmax=0.01*maxAltitude; //the minimum and maximum resolution along the 'altitude' line

  //the 'altitude' points are distributed logarithmicaly
  //determin the factor 'rr' of the geometrical progression that defined the location od radial points of the map
  long int t;
  double dAlt,rr,R=dAltmin;

  dAlt=dAltmin;
  rr=(maxAltitude+dAltmax)/(maxAltitude+dAltmin);
  t=(long int)(log(dAltmax/dAltmin)/log(rr)-2.0);
  rr=pow(dAltmax/dAltmin,1.0/(t+2.0));

  //determine positions of the Easth and the Sun
  SpiceDouble State[6],lt;

  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,State,&lt);
  rEarth[0]=State[0]*1.0E3;
  rEarth[1]=State[1]*1.0E3;
  rEarth[2]=State[2]*1.0E3;

  spkezr_c("Sun",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,State,&lt);
  rSun[0]=State[0]*1.0E3;
  rSun[1]=State[1]*1.0E3;
  rSun[2]=State[2]*1.0E3;



  //open the output file
  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s, Radial size of the planet=%e\"\n",utcstr,atan(_RADIUS_(_TARGET_)/sqrt(rEarth[0]*rEarth[0]+rEarth[1]*rEarth[1]+rEarth[2]*rEarth[2]))/Pi*180.0);

    fprintf(fout,"VARIABLES=\"R [Object Radii]\", \"R [m]\", \"Angle from the center of the object [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");
    fprintf(fout,"ZONE T=\"Column Density Along the Limb\"\n");


    //output emmision and column integrals at the limb at different phase angles
    if (PIC::DataOutputFileNumber==0) {
      fLimb=fopen("pic.LimbIntegrals.dat","w");
      fprintf(fLimb,"VARIABLES=\"Julian Date\", \"Phase Angle[degrees]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");
    }
    else fLimb=fopen("pic.LimbIntegrals.dat","a");
  }


  //determine the localtion of the limb
  double a,b,c,d,t4,t6,t7,t8,t10,t11,t12,t13,t14,t18;

  t4 = (rEarth[1] * rSun[1] + rEarth[0] * rSun[0] + rEarth[2] * rSun[2]);
  t6 = pow(rSun[2], 0.2e1);
  t7 = pow(rSun[0], 0.2e1);
  t8 = pow(rSun[1], 0.2e1);
  t10 = pow(rEarth[0], 0.2e1);
  t11 = _RADIUS_(_TARGET_) * _RADIUS_(_TARGET_);
  t12 = pow(rEarth[2], 0.2e1);
  t13 = pow(rEarth[1], 0.2e1);
  t14 = t10 - t11 + t12 + t13;
  t18 = -t13 - t10 - t12;

  a = (4 * t4 * t4) - 0.4e1 * (t6 + t7 + t8) * t14;
  b = 0.8e1 * t18 * (double) t4 + 0.8e1 * (double) t4 * t14;
  c = 0.4e1 * t18 * t18 + 0.4e1 * t18 * t14;

  d=b*b-4.0*a*c;
  MPI_Bcast(&d,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  if (d<0.0) {
    if (PIC::ThisThread==0) cout << "WARNINIG: The limb can not being seen for the given orbital configuration\n";
    return;
  }

  double tSunLine=0.0;

  if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

  tSunLine=(-b-sqrt(d))/(2.0*a);
  if (tSunLine<0.0) tSunLine=-2.0*c/(b+sqrt(d));


  //direction to the limb
  double lLimbObject[3];

  for (idim=0,c=0.0;idim<3;idim++) {
    l[idim]=rSun[idim]*tSunLine-rEarth[idim];
    c+=pow(l[idim],2);
  }

  for (c=sqrt(c),idim=0;idim<3;idim++) l[idim]/=c;

  //distance between the Earth and the position of the Limb
  c=sqrt(rEarth[0]*rEarth[0]+rEarth[1]*rEarth[1]+rEarth[2]*rEarth[2]-_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_));

  for (idim=0,d=0.0;idim<3;idim++) {
    lLimbObject[idim]=rEarth[idim]+c*l[idim];
    d+=pow(lLimbObject[idim],2);
  }

  for (d=sqrt(d),idim=0;idim<3;idim++) lLimbObject[idim]/=d;

  //determine the direction of the line that connects the center of the Object and the position of the limb



  //calculate the integrals
  const int StateVectorLength=3;
  double StateVector[StateVectorLength];
  double StateVectorMax[StateVectorLength]; //maximum values of the state variables == the values of the state variables at the limb

  for (int i=0;i<StateVectorLength;i++) StateVectorMax[i]=0.0;

  while (R<maxAltitude+1.0) {

    //get the pointing direction
    for (c=0.0,idim=0;idim<3;idim++) {
      l[idim]=R*lLimbObject[idim]*_RADIUS_(_TARGET_)-rEarth[idim];
      c+=pow(l[idim],2);
    }

    for (c=sqrt(c),idim=0;idim<3;idim++) l[idim]/=c;

    PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,rEarth,l,SodiumCoulumnDensityIntegrant);
    for (int i=0;i<StateVectorLength;i++) if (StateVectorMax[i]<StateVector[i]) StateVectorMax[i]=StateVector[i];

    if (PIC::ThisThread==0) {
      fprintf(fout,"%e   %e   %e   %e   %e   %e\n",R,R*_RADIUS_(_TARGET_),atan(R*_RADIUS_(_TARGET_)/sqrt(rEarth[0]*rEarth[0]+rEarth[1]*rEarth[1]+rEarth[2]*rEarth[2]))/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);
    }

    R+=dAlt;
    dAlt*=rr;
  }

  if (PIC::ThisThread==0) {
    fclose(fout);

    //output the time-dependent values of the integrals at the limb
    fprintf(fLimb,"%e %e  %e  %e  %e\n",j2000_c()+Exosphere::OrbitalMotion::et/spd_c(),Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et)/Pi*180.0,StateVectorMax[0],StateVectorMax[1],StateVectorMax[2]);
    fclose(fLimb);
  }
}
*/


/*
void Exosphere::ColumnIntegral::Map(char *fname,double dXmax,double dZmax,int nXpoints) {
  double l[3],xEarth_new[3]={0.0,0.0,0.0};

  //find position of Earth
  SpiceDouble State[6],lt;

  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,State,&lt);

  xEarth_new[0]=State[0]*1.0E3;
  xEarth_new[1]=State[1]*1.0E3;
  xEarth_new[2]=State[2]*1.0E3;


  if (dXmax<0.0) dXmax=2.0*max(fabs(PIC::Mesh::mesh.xGlobalMax[0]),fabs(PIC::Mesh::mesh.xGlobalMin[0]));
  if (dZmax<0.0) dZmax=3.0*max(fabs(PIC::Mesh::mesh.xGlobalMax[2]),fabs(PIC::Mesh::mesh.xGlobalMin[2]));

  const double maxPhiX=1.1*atan(dXmax/(1.0*_AU_));
  const double maxPhiZ=1.1*atan(dZmax/(1.0*_AU_));
  const double minPhiX=-maxPhiX,minPhiZ=-maxPhiZ;

  const double dPhi=(maxPhiX-minPhiX)/(nXpoints-1);
  const int nPointsZ=1+(int)((maxPhiZ-minPhiZ)/dPhi);

  //open the output file
  FILE *fout=NULL;
  int iZ,iX;
  double r,PhiX,PhiZ,StateVector[3];


  if (PIC::ThisThread==0) {
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Angle In Ecliptic Plane [degree]\", \"Angle Out of Ecpliptic Plane [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nPointsZ,nXpoints);
  }

  r=sqrt(xEarth_new[0]*xEarth_new[0]+xEarth_new[1]*xEarth_new[1]+xEarth_new[2]*xEarth_new[2]);

  for (iX=0;iX<nXpoints;iX++) {
    //determine the X,Y-components of the direction vector
    PhiX=maxPhiX-dPhi*iX;

    //rotate the vector
    l[0]=-(cos(PhiX)*xEarth_new[0]-sin(PhiX)*xEarth_new[1]);
    l[1]=-(sin(PhiX)*xEarth_new[0]+cos(PhiX)*xEarth_new[1]);

    for (iZ=0;iZ<nPointsZ;iZ++) {
      //determine the Z-component of the direction vector
      PhiZ=minPhiZ+dPhi*iZ;
      l[2]=r*tan(PhiZ)-xEarth_new[2];

      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,3,xEarth_new,l,ColumnIntegral::CoulumnDensityIntegrant);

      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",-PhiX/Pi*180.0,PhiZ/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);
    }

  }

  if (PIC::ThisThread==0) fclose(fout);
}
*/

/*
void Exosphere::ColumnDensityIntegration_CircularMap(char *fname,double rmax,double dRmin,double dRmax,int nAzimuthPoints,SpiceDouble EphemerisTime) {
  double l[3],xEarth_new[3]={0.0,0.0,0.0};
  double rr,dR,R=0.0;
  FILE *fout=NULL;

  //find position of Earth
  SpiceDouble State[6],lt;

  spkezr_c("Earth",EphemerisTime,SO_FRAME,"none",ObjectName,State,&lt);

  xEarth_new[0]=State[0]*1.0E3;
  xEarth_new[1]=State[1]*1.0E3;
  xEarth_new[2]=State[2]*1.0E3;


  if (rmax<0.0) rmax=2.0*max(fabs(PIC::Mesh::mesh.xGlobalMax[0]),fabs(PIC::Mesh::mesh.xGlobalMin[0]));
  if (dRmin<0.0) dRmin=0.1*_RADIUS_(_TARGET_);
  if (dRmax<0.0) dRmax=10.0*_RADIUS_(_TARGET_);

  //determin the factor 'rr' of the geometrical progression that defined the location od radial points of the map
  long int t;

  dR=dRmin;
  rr=(rmax+dRmax)/(rmax+dRmin);
  t=(long int)(log(dRmax/dRmin)/log(rr)-2.0);
  rr=pow(dRmax/dRmin,1.0/(t+2.0));

  //determine the coordinate frame with z-axis align with the Eath-Mercury direction
  //e2 -> direction from Mercury to the Earth, e0 -> projection of the z-axis ofthe SO frame, e1 -> cross product of e2 and e0
  double e0[3],e1[3],e2[3],c,c1;
  int idim;

  for (idim=0,c=0.0;idim<3;idim++) {
    e2[idim]=xEarth_new[idim];
    c+=pow(e2[idim],2);
  }

  for (c=sqrt(c),idim=0;idim<3;idim++) e2[idim]/=c;


  //determine e0
  e0[0]=-e2[2]*e2[0];
  e0[1]=-e2[2]*e2[1];
  e0[2]=1.0-e2[2]*e2[2];

  c=sqrt(e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2]);
  for (idim=0;idim<3;idim++) e0[idim]/=c;

  e1[0]=e2[1]*e0[2]-e0[1]*e2[2];
  e1[1]=e0[0]*e2[2]-e2[0]*e0[2];
  e1[2]=e2[0]*e0[1]-e0[0]*e2[1];

  //disrance between Mercury and Earth
  double rEarth2Mercury=sqrt(pow(xEarth_new[0],2)+pow(xEarth_new[1],2)+pow(xEarth_new[2],2));

  //calculate the number of the radial points
  int nRadialPoints=0;
  c1=dR,c=R;

  while (c<rmax) {
    c+=c1;
    c1*=rr;
    nRadialPoints++;
  }


  //open the output file
  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];

    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s, Radial size of the planet=%e\"\n",utcstr,atan(_RADIUS_(_TARGET_)/rEarth2Mercury)/Pi*180.0);

    fprintf(fout,"VARIABLES=\"Angle In Ecliptic Plane [degree]\", \"Angle Out of Ecpliptic Plane [degree]\", \"Column Density [m^{-2}]\", \"Intensity (5891.58A) [R]\", \"Intensity (5897.56A) [R]\" \n");
    fprintf(fout,"ZONE T=\"Column Density Map\"\n");
    fprintf(fout,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nRadialPoints,nAzimuthPoints+1);
    fprintf(fout,"DATAPACKING=POINT\n");
    fprintf(fout,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
  }


  //calcualte the column integrals
  const int StateVectorLength=3;

  double l0Earth2Mercury[3],lEarth2Mercury[3],phi;
  double PhiX,PhiZ,StateVector[StateVectorLength];

  for (int iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints+1;iAzimuthPoint++) {
    R=0.0,dR=dRmin;
    l0Earth2Mercury[0]=R;
    l0Earth2Mercury[1]=0.0;
    l0Earth2Mercury[2]=-rEarth2Mercury;

    for (int np=0;np<nRadialPoints;np++) {
      //rotate the pointing vector along the Mercury - Earth axis
      phi=2.0*Pi*double(iAzimuthPoint)/double(nAzimuthPoints);

      lEarth2Mercury[0]=cos(phi)*l0Earth2Mercury[0]-sin(phi)*l0Earth2Mercury[1];
      lEarth2Mercury[1]=sin(phi)*l0Earth2Mercury[0]+cos(phi)*l0Earth2Mercury[1];
      lEarth2Mercury[2]=l0Earth2Mercury[2];


      l[0]=lEarth2Mercury[0]*e0[0]+lEarth2Mercury[1]*e1[0]+lEarth2Mercury[2]*e2[0];
      l[1]=lEarth2Mercury[0]*e0[1]+lEarth2Mercury[1]*e1[1]+lEarth2Mercury[2]*e2[1];
      l[2]=lEarth2Mercury[0]*e0[2]+lEarth2Mercury[1]*e1[2]+lEarth2Mercury[2]*e2[2];


      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarth_new,l,SodiumCoulumnDensityIntegrant);

      PhiZ=atan(lEarth2Mercury[0]/rEarth2Mercury);
      PhiX=-atan(lEarth2Mercury[1]/rEarth2Mercury);

      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",PhiX/Pi*180.0,PhiZ/Pi*180.0,StateVector[0],StateVector[1],StateVector[2]);


    R+=dR;
    dR*=rr;

    l0Earth2Mercury[0]=R;

    }
  }


  if (PIC::ThisThread==0) fclose(fout);
}
*/
void Exosphere::ColumnIntegral::CircularMap(char *fname,double rmax,double dRmin,double dRmax,int nAzimuthPoints,SpiceDouble EphemerisTime) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  double l[3],xEarthSO[3]={0.0,0.0,0.0},lEarth;
  double rr,dR,R=0.0;
  FILE *fout=NULL;

  //find position of Earth
  SpiceDouble EarthStateSO[6],ObjectStateGSE[6],lt;

  spkezr_c("Earth",EphemerisTime,SO_FRAME,"none",ObjectName,EarthStateSO,&lt);
  spkezr_c(ObjectName,EphemerisTime,"GSE","none","Earth",ObjectStateGSE,&lt);

  xEarthSO[0]=EarthStateSO[0]*1.0E3;
  xEarthSO[1]=EarthStateSO[1]*1.0E3;
  xEarthSO[2]=EarthStateSO[2]*1.0E3;

  lEarth=sqrt(pow(xEarthSO[0],2)+pow(xEarthSO[1],2)+pow(xEarthSO[2],2));

  if (rmax<0.0) rmax=2.0*max(fabs(PIC::Mesh::mesh.xGlobalMax[0]),fabs(PIC::Mesh::mesh.xGlobalMin[0]));
  if (dRmin<0.0) dRmin=0.1*_RADIUS_(_TARGET_);
  if (dRmax<0.0) dRmax=10.0*_RADIUS_(_TARGET_);

  //determin the factor 'rr' of the geometrical progression that defined the location od radial points of the map
  long int t;

  dR=dRmin;
  rr=(rmax+dRmax)/(rmax+dRmin);
  t=(long int)(log(dRmax/dRmin)/log(rr)-2.0);
  rr=pow(dRmax/dRmin,1.0/(t+2.0));


  double e0[3],e1[3],e2[3],c,c1;
  int idim;



  //determine the coordinate frame with z-axis align with the Eath-Mercury direction
  //e2 -> direction from Mercury to the Earth, e1 -> a part of GSE z-axis ortogonal to e2, and e0 -> to get the right handed frame
/*  for (idim=0,c=0.0;idim<3;idim++) {
    e2[idim]=-EarthStateSO[idim];
    c+=pow(e2[idim],2);
  }

  for (c=sqrt(c),idim=0;idim<3;idim++) e2[idim]/=c;

  //determine e1 (projection of z-axix)
  e1[0]=-e2[2]*e2[0];
  e1[1]=-e2[2]*e2[1];
  e1[2]=1.0-e2[2]*e2[2];

  c=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
  for (idim=0;idim<3;idim++) e1[idim]/=c;

  e0[0]=e2[1]*e1[2]-e1[1]*e2[2];
  e0[1]=e1[0]*e2[2]-e2[0]*e1[2];
  e0[2]=e2[0]*e1[1]-e1[0]*e2[1];*/


  double e0GSE[3],e1GSE[3],e2GSE[3];

  for (idim=0,c=0.0;idim<3;idim++) {
    e2GSE[idim]=ObjectStateGSE[idim];
    c+=pow(e2GSE[idim],2);
  }

  for (c=sqrt(c),idim=0;idim<3;idim++) e2GSE[idim]/=c;

  //determine e1 (projection of z-axix)
  e1GSE[0]=-e2GSE[2]*e2GSE[0];
  e1GSE[1]=-e2GSE[2]*e2GSE[1];
  e1GSE[2]=1.0-e2GSE[2]*e2GSE[2];

  c=sqrt(e1GSE[0]*e1GSE[0]+e1GSE[1]*e1GSE[1]+e1GSE[2]*e1GSE[2]);
  for (idim=0;idim<3;idim++) e1GSE[idim]/=c;

  e0GSE[0]=-(e2GSE[1]*e1GSE[2]-e1GSE[1]*e2GSE[2]);
  e0GSE[1]=-(e1GSE[0]*e2GSE[2]-e2GSE[0]*e1GSE[2]);
  e0GSE[2]=-(e2GSE[0]*e1GSE[1]-e1GSE[0]*e2GSE[1]);

  SpiceDouble xform[6][6],eGSE2SO[6]={0.0,0.0,0.0,0.0,0.0,0.0},eGSE2SOres[6]={0.0,0.0,0.0,0.0,0.0,0.0};

  sxform_c("GSE",SO_FRAME,EphemerisTime,xform);

  eGSE2SO[0]=e0GSE[0],eGSE2SO[1]=e0GSE[1],eGSE2SO[2]=e0GSE[2];
  mxvg_c(xform,eGSE2SO,6,6,eGSE2SOres);
  e0[0]=eGSE2SOres[0],e0[1]=eGSE2SOres[1],e0[2]=eGSE2SOres[2];

  eGSE2SO[0]=e1GSE[0],eGSE2SO[1]=e1GSE[1],eGSE2SO[2]=e1GSE[2];
  mxvg_c(xform,eGSE2SO,6,6,eGSE2SOres);
  e1[0]=eGSE2SOres[0],e1[1]=eGSE2SOres[1],e1[2]=eGSE2SOres[2];

  eGSE2SO[0]=e2GSE[0],eGSE2SO[1]=e2GSE[1],eGSE2SO[2]=e2GSE[2];
  mxvg_c(xform,eGSE2SO,6,6,eGSE2SOres);
  e2[0]=eGSE2SOres[0],e2[1]=eGSE2SOres[1],e2[2]=eGSE2SOres[2];


  //calculate the number of the radial points
  int nRadialPoints=0;
  c1=dR,c=R;

  while (c<rmax) {
    c+=c1;
    c1*=rr;
    nRadialPoints++;
  }


  //open the output file
  if (PIC::ThisThread==0) {
    const SpiceInt lenout = 35;
    SpiceChar utcstr[lenout+2];
    char vlist[_MAX_STRING_LENGTH_PIC_]="";

    ColumnIntegral::GetVariableList(vlist);
    fout=fopen(fname,"w");

    et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
    fprintf(fout,"TITLE=\"UTC=%s, Radial size of the planet=%e\"\n",utcstr,atan(_RADIUS_(_TARGET_)/lEarth)/Pi*180.0);

    fprintf(fout,"VARIABLES=\"l0 In Ecliptic Plane []\", \"l1 Out of Ecpliptic Plane []\" %s\n",vlist);
    fprintf(fout,"ZONE T=\"Column Density Map\"\n");
    fprintf(fout,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nAzimuthPoints+1,nRadialPoints);
    fprintf(fout,"DATAPACKING=POINT\n");
  }


  //calcualte the column integrals
  int StateVectorLength=ColumnIntegral::GetVariableList(NULL);
  double StateVector[StateVectorLength];
  double AzimuthAngle,l0[3],ZenithAngle;


  R=0.0,dR=dRmin;

  for (int np=0;np<nRadialPoints;np++) {
    ZenithAngle=atan(R/lEarth);

    for (int iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints+1;iAzimuthPoint++) {
      AzimuthAngle=2.0*Pi*double(iAzimuthPoint)/double(nAzimuthPoints);

      l0[2]=cos(ZenithAngle);
      l0[1]=sin(ZenithAngle)*sin(AzimuthAngle);
      l0[0]=sin(ZenithAngle)*cos(AzimuthAngle);

      //transfer the pointing vector to 'SO_FRAME'
      l[0]=l0[0]*e0[0]+l0[1]*e1[0]+l0[2]*e2[0];
      l[1]=l0[0]*e0[1]+l0[1]*e1[1]+l0[2]*e2[1];
      l[2]=l0[0]*e0[2]+l0[1]*e1[2]+l0[2]*e2[2];

      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,xEarthSO,l,ColumnIntegral::CoulumnDensityIntegrant);
      ColumnIntegral::ProcessColumnIntegrationVector(StateVector,StateVectorLength);

      if (PIC::ThisThread==0) {
        fprintf(fout,"%e   %e",l0[0],l0[1]);
        for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",StateVector[i]);
        fprintf(fout,"\n");
      }

    }

    R+=dR;
    dR*=rr;
  }

  if (PIC::ThisThread==0) fclose(fout);
#endif
}

//=====================================================================================================
//empty sampling function
void Exosphere::Sampling::SampleModelData() {

  //sample sodium surface content
  int nZenithSurfaceElements,nAzimuthalSurfaceElements,spec;
  long int iZenith_SO,iAzimuth_SO,el_SO;
  long int iZenith_IAU,iAzimuth_IAU,el_IAU;
  double rSurfaceElement_IAU[3],rSurfaceElement_SO[3],SurfaceArea_IAU;
  SpiceDouble xform[6][6];

  nZenithSurfaceElements=Planet->nZenithSurfaceElements;
  nAzimuthalSurfaceElements=Planet->nAzimuthalSurfaceElements;

  memcpy(xform,OrbitalMotion::SO_to_IAU_TransformationMartix,36*sizeof(double));

  if (Planet!=NULL) for (iZenith_SO=0;iZenith_SO<nZenithSurfaceElements;iZenith_SO++) for (iAzimuth_SO=0;iAzimuth_SO<nAzimuthalSurfaceElements;iAzimuth_SO++)  {
    //generate position on the sphere in the coordinate frame SO (x-axis is directed to the Sun)
    Planet->GetSurfaceElementMiddlePoint(rSurfaceElement_SO,iZenith_SO,iAzimuth_SO);
    el_SO=Planet->GetLocalSurfaceElementNumber(iZenith_SO,iAzimuth_SO);

    //conver the position vector into the frane related to the planet
    rSurfaceElement_IAU[0]=xform[0][0]*rSurfaceElement_SO[0]+xform[0][1]*rSurfaceElement_SO[1]+xform[0][2]*rSurfaceElement_SO[2];
    rSurfaceElement_IAU[1]=xform[1][0]*rSurfaceElement_SO[0]+xform[1][1]*rSurfaceElement_SO[1]+xform[1][2]*rSurfaceElement_SO[2];
    rSurfaceElement_IAU[2]=xform[2][0]*rSurfaceElement_SO[0]+xform[2][1]*rSurfaceElement_SO[1]+xform[2][2]*rSurfaceElement_SO[2];

    //get the surface element of the sampled point in the frame related to the planet
    Planet->GetSurfaceElementProjectionIndex(rSurfaceElement_IAU,iZenith_IAU,iAzimuth_IAU);
    SurfaceArea_IAU=Planet->GetSurfaceElementArea(iZenith_IAU,iAzimuth_IAU);


    el_IAU=Planet->GetLocalSurfaceElementNumber(iZenith_IAU,iAzimuth_IAU);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      double SpecieSurfaceDensity;

      SpecieSurfaceDensity=Planet->SurfaceElementPopulation[spec][el_IAU]/SurfaceArea_IAU;
      Planet->SampleSpeciesSurfaceAreaDensity[spec][el_SO]+=SpecieSurfaceDensity;
    }

  }
}


//output the column integrated density
void Exosphere::Sampling::OutputSampledModelData(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  char fname[_MAX_STRING_LENGTH_PIC_];
  int ierr;

  FILE *fSource=NULL;

  //print sodium production rate
  double SourceRate=0.0,TotalSourceRate=0.0;
  int thread;
  double buffer[PIC::nTotalThreads];

  //print the simulation time
  const SpiceInt lenout = 35;
  SpiceChar utcstr[lenout+2];


  et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);

  //save the surface properties
  if (PIC::ThisThread==0) {
    char fname[300];

    sprintf(fname,"%s/pic.SurfaceProperties.%s.out=%i.dat",PIC::OutputDataFileDirectory,utcstr,DataOutputFileNumber);
    Exosphere::Planet->SaveSurfaceDensity(fname,Planet->GetTotalSurfaceElementsNumber(),PIC::nTotalSpecies);
  }

  if (PIC::ThisThread==0) {
    printf("$PREFIX:Output file number: %i\n",DataOutputFileNumber);
    printf("$PREFIX:Simulation Time: %s, et=%e,TAA=%e [deg]\n",utcstr,Exosphere::OrbitalMotion::et,Exosphere::OrbitalMotion::TAA/Pi*180.0);
    printf("$PREFIX:rHeliocenric=%e (%e AU), vRadial=%e\n",Exosphere::xObjectRadial,Exosphere::xObjectRadial/_AU_,Exosphere::vObjectRadial);
    printf("$PREFIX:Phase angle=%e [degrees]\n",Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et)*180.0/Pi);
    printf("$PREFIX:Sodium RadiationPressure Acceleration=%e\n",SodiumRadiationPressureAcceleration__Combi_1997_icarus(Exosphere::vObjectRadial,Exosphere::xObjectRadial));

    printf("\n$PREFIX:Source rates:\n");
  }


  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    SourceRate=0.0,TotalSourceRate=0.0;

    if (PIC::ThisThread==0) {
      char fname[_MAX_STRING_LENGTH_PIC_];

      sprintf(fname,"%s/pic.SourceRate.%s.spec=%i.dat",PIC::OutputDataFileDirectory,PIC::MolecularData::GetChemSymbol(spec),spec);

      if (DataOutputFileNumber==0) {
        fSource=fopen(fname,"w");
        fprintf(fSource,"VARIABLES=\"Time [JD] \", \"TAA [degrees]\", \"Phase Angle [degrees]\", \"PSD [m^{-2} s^{-1}]\", \"IV [m^{-2} s^{-1}]\", \"TD [m^{-2} s^{-1}]\", \"SWS [m^{-2} s^{-1}]\",\"Total Source Rate [m^{-2} s^{-1}]\","
          "\"Na Radiation Pressure Acceleration [m/s^2]\", \"vObject [m/s]\", \"rObject [AU]\", \"Output file number\" ");

        if (UserDefinedAdditionalData_VariableList_OutputSampledModelData!=NULL) {
          UserDefinedAdditionalData_VariableList_OutputSampledModelData(fSource);
        }

        fprintf(fSource,"\n");
      } else {
        fSource=fopen(fname,"a");
      }

      et2utc_c(Exosphere::OrbitalMotion::et,"J",5,lenout,utcstr);
      fprintf(fSource,"%s  %e  %e",utcstr+3,Exosphere::OrbitalMotion::TAA/Pi*180.0,Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et)/Pi*180.0);
    }


    //calculate the source rate
    for (int iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
      ierr=MPI_Gather(&CalculatedSourceRate[spec][iSource],1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
      if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

      for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
      if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

      if (PIC::ThisThread==0) {
        printf("$PREFIX:%s source: %s rate - %e s^{-1}\n",PIC::MolecularData::GetChemSymbol(spec),_EXOSPHERE__SOURCE_SYMBOLIC_ID_[iSource],SourceRate);
        fprintf(fSource,"  %e",SourceRate);

        TotalSourceRate+=SourceRate;
      }

      CalculatedSourceRate[spec][iSource]=0.0;
    }

    if (PIC::ThisThread==0) {
      printf("$PREFIX:Total %s production rate - %e s^{-1}\n\n",PIC::MolecularData::GetChemSymbol(spec),TotalSourceRate);

      fprintf(fSource,"  %e  %e  %e  %i\n",TotalSourceRate,Exosphere::vObjectRadial,Exosphere::xObjectRadial/_AU_,DataOutputFileNumber);

      if (UserDefinedAdditionalData_OutputSampledModelData!=NULL) {
        UserDefinedAdditionalData_OutputSampledModelData(fSource,spec);
      }

      fclose(fSource);
    }
  }



/*
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  ierr=MPI_Gather(&Exosphere::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Sodium source: Photon Stimulated Desorption rate - %e s^{-1}\n",SourceRate);
    fprintf(fSource,"  %e",SourceRate);

    TotalSourceRate+=SourceRate;
  }

  Exosphere::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate=0.0;
#else
  if (PIC::ThisThread==0) fprintf(fSource,"  0.0");
#endif


#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  ierr=MPI_Gather(&Exosphere::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Sodium source: Impact Vaporization rate - %e s^{-1}\n",SourceRate);
    fprintf(fSource,"  %e",SourceRate);

    TotalSourceRate+=SourceRate;
  }

  Exosphere::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate=0.0;
#else
  if (PIC::ThisThread==0) fprintf(fSource,"  0.0");
#endif


#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  ierr=MPI_Gather(&Exosphere::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Sodium source: Thermal Desorption rate - %e s^{-1}\n",SourceRate);
    fprintf(fSource,"  %e",SourceRate);

    TotalSourceRate+=SourceRate;
  }

  Exosphere::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate=0.0;
#else
  if (PIC::ThisThread==0) fprintf(fSource,"  0.0");
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  ierr=MPI_Gather(&Exosphere::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Sodium source: Solar Wind Sputtering rate - %e s^{-1}\n",SourceRate);
    fprintf(fSource,"  %e",SourceRate);

    TotalSourceRate+=SourceRate;
  }

  Exosphere::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate=0.0;
#else
  if (PIC::ThisThread==0) fprintf(fSource,"  0.0");
#endif
*/

/*  if (PIC::ThisThread==0) {
    printf("Total sodium production rate - %e s^{-1}\n",TotalSourceRate);

    fprintf(fSource,"  %e  %e  %e  %e  %i\n",TotalSourceRate,
        SodiumRadiationPressureAcceleration__Combi_1997_icarus(Exosphere::vObjectRadial,Exosphere::xObjectRadial),Exosphere::vObjectRadial,Exosphere::xObjectRadial/_AU_,DataOutputFileNumber);

    fclose(fSource);
    printf("\nNight side return fluxes");
  }*/

  //print the total night side flux

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    double flux[_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1],FluxAll[PIC::nTotalThreads*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];
    char ChemSymbol[_MAX_STRING_LENGTH_PIC_];
    double LocalTimeStep;
    int i;

    PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

    ierr=MPI_Gather(Sampling::PlanetNightSideReturnFlux[spec],_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1,MPI_DOUBLE,FluxAll,_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1, MPI_DOUBLE,0, MPI_GLOBAL_COMMUNICATOR);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

    for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) {
      for (flux[i]=0.0,thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) flux[i]+=FluxAll[i+thread*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];

      if (PIC::LastSampleLength!=0) flux[i]/=PIC::LastSampleLength*LocalTimeStep;
    }


    if (PIC::ThisThread==0) {
       if (spec==0) printf("$PREFIX:Night side return flux [s^{-1}]:\n");

       printf("$PREFIX:Spec=%i (%s):\n",spec,ChemSymbol);

       for (int iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
         printf("$PREFIX:Particles produced by %s: %e\n",_EXOSPHERE__SOURCE_SYMBOLIC_ID_[iSource],flux[iSource]);
       }

       printf("$PREFIX:\n");
/*

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
       printf("Particles produced by  Impact Vaporization: %e\n",flux[_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
       printf("Particles produced by Photon Stimulated Desorption: %e\n",flux[_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
       printf("Particles produced by Thermal Desorption: %e\n",flux[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
       printf("Particles produced by Solar wind sputtering: %e\n",flux[_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]);
#endif
*/

    }

    for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) Sampling::PlanetNightSideReturnFlux[spec][i]=0.0;
  }


  //print the total return flux and the total rate of sticking to the surface
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    double TotalFlux;
    char ChemSymbol[_MAX_STRING_LENGTH_PIC_];
    double LocalTimeStep;

    if ((spec==0)&&(PIC::ThisThread==0)) {
      printf("\n$PREFIX:The total return flux and sticking rate\n$PREFIX:Spec\tTotal Return Flux [1/s]\tTotal Surface Sticking Rate [1/s]\n");
    }

    PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);
    if (PIC::ThisThread==0) printf("$PREFIX:%i (%s)\t",spec,ChemSymbol);

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

    ierr=MPI_Reduce(Sampling::TotalPlanetReturnFlux+spec,&TotalFlux,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);
    if (PIC::LastSampleLength!=0) TotalFlux/=PIC::LastSampleLength*LocalTimeStep;
    if (PIC::ThisThread==0) printf("%e\t",TotalFlux);

    ierr=MPI_Reduce(Sampling::PlanetSurfaceStickingRate+spec,&TotalFlux,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);
    if (PIC::LastSampleLength!=0) TotalFlux/=PIC::LastSampleLength*LocalTimeStep;
    if (PIC::ThisThread==0) printf("%e\n",TotalFlux);

    Sampling::TotalPlanetReturnFlux[spec]=0.0,Sampling::PlanetSurfaceStickingRate[spec]=0.0;
  }

  //output the trajectory file
  if (PIC::ThisThread==0) {
    SpiceDouble etStart,dEt;
    int nTrajectoryPoint;
    FILE *fTrajectory;
    SpiceDouble State[6],lt,IAU2SO[6][6];
    char OrbitalDataFileName[_MAX_STRING_LENGTH_PIC_];

    sprintf(OrbitalDataFileName,"%s/pic.OrbitalData.dat",PIC::OutputDataFileDirectory);
    fTrajectory=fopen(OrbitalDataFileName,"a");

    dEt=PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::LastSampleLength/Exosphere::OrbitalMotion::nOrbitalPositionOutputMultiplier;
    etStart=Exosphere::OrbitalMotion::et-dEt*(Exosphere::OrbitalMotion::nOrbitalPositionOutputMultiplier-1);

    for (nTrajectoryPoint=0;nTrajectoryPoint<Exosphere::OrbitalMotion::nOrbitalPositionOutputMultiplier;nTrajectoryPoint++) {
      et2utc_c(etStart,"ISOC",0,lenout,utcstr);

      fprintf(fTrajectory,"%s",utcstr);

      spkezr_c(ObjectName,etStart,"J2000","none","SUN",State,&lt);
      fprintf(fTrajectory, "  %e  %e  %e",State[0],State[1],State[2]);

      spkezr_c("Earth",etStart,"J2000","none","SUN",State,&lt);
      fprintf(fTrajectory, "  %e  %e  %e",State[0],State[1],State[2]);

      spkezr_c("SUN",etStart,SO_FRAME,"none",ObjectName,State,&lt);
      fprintf(fTrajectory, "  %e  %e  %e",State[0],State[1],State[2]);

      spkezr_c("Earth",etStart,SO_FRAME,"none",ObjectName,State,&lt);
      fprintf(fTrajectory, "  %e  %e  %e",State[0],State[1],State[2]);

      et2utc_c(Exosphere::OrbitalMotion::et,"ISOC",0,lenout,utcstr);
      fprintf(fTrajectory, "  %s  %i",utcstr,DataOutputFileNumber);

      sxform_c(IAU_FRAME,SO_FRAME,etStart,IAU2SO);
      for (int i=0;i<3;i++) for (int j=0;j<3;j++) fprintf(fTrajectory,"  %e ",IAU2SO[i][j]);

      fprintf(fTrajectory, "\n");
      etStart+=dEt;
    }

    fclose(fTrajectory);
  }


  //the sodium column density in the tail
  sprintf(fname,"%s/pic.TailColumnDensity.%s.out=%i.dat",PIC::OutputDataFileDirectory,utcstr,DataOutputFileNumber);
  ColumnIntegral::Tail(fname);

  //column density distribution map
/*
  sprintf(fname,"pic.ColumnDensityMap.%s.out=%i.dat",utcstr,DataOutputFileNumber);
  ColumnDensityIntegration_Map(fname,-1.0,-1.0,100);

  sprintf(fname,"pic.ColumnDensityMap.close.%s.out=%i.dat",utcstr,DataOutputFileNumber);
  ColumnDensityIntegration_Map(fname,40.0*_RADIUS_(_TARGET_),40.0*_RADIUS_(_TARGET_),100);
*/

  sprintf(fname,"%s/pic.ColumnDensityMap.spherical.%s.out=%i.dat",PIC::OutputDataFileDirectory,utcstr,DataOutputFileNumber);
  double domainCharacteristicSize=0.0;
  for (int idim=0;idim<DIM;idim++) domainCharacteristicSize=max(max(fabs(PIC::Mesh::mesh.xGlobalMax[idim]),fabs(PIC::Mesh::mesh.xGlobalMin[idim])),domainCharacteristicSize);

  Exosphere::ColumnIntegral::CircularMap(fname,domainCharacteristicSize,0.05*_RADIUS_(_TARGET_),min(5*_RADIUS_(_TARGET_),0.1*domainCharacteristicSize),80,Exosphere::OrbitalMotion::et);

  //sodium column density along the limb direction
  sprintf(fname,"%s/pic.LimbColumnDensity.%s.out=%i.dat",PIC::OutputDataFileDirectory,utcstr,DataOutputFileNumber);
  ColumnIntegral::Limb(fname);


  //output data that correcponds for avalable ground based observatinos
  Exosphere::Sampling::ReferenceGroundBasedObservations::OutputSampledData(Exosphere::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::LastSampleLength,Exosphere::OrbitalMotion::et,DataOutputFileNumber);
#endif
}



/*--------------------------------- Print Output File: BEGIN  --------------------------------------*/
void Exosphere::Sampling::OutputDataFile::PrintVariableList(FILE* fout,int DataSetNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_

  //coordinate in the J2000 frame
  sxform_c(SO_FRAME,"J2000",Exosphere::OrbitalMotion::et,SO_to_HCI_TransformationMartix);
  fprintf(fout,", \"xJ2000\", \"yJ2000\", \"zJ2000\"");
#endif

  if (SamplingDensityOffset[0]!=-1) for (int iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
    fprintf(fout,", \"Number Density(spec=%s,source=%s)\"", PIC::MolecularData::GetChemSymbol(DataSetNumber),_EXOSPHERE__SOURCE_SYMBOLIC_ID_[iSource]);
  }

/*
  //macroscopic properties
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Impact Vaporization Source)\"");
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Photon Stimulated Desorption Source)\"");
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Thermal Desorption Source)\"");
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Solar Wind Sputtering Source)\"");
#endif
*/

  //add variable list from a user-defined output
#if _EXOSPHERE__USER_DEFINED_FILE_OUTPUT__MODE__ == _EXOSPHERE_SOURCE__ON_
  _EXOSPHERE__USER_DEFINED_FILE_OUTPUT__VARIABLE_LIST_(fout);
#endif

/*
  //mean particle radiation pressure acceleartion and g-factor
  fprintf(fout,", \"Radiation Pressure Acceleration [m/s^2]\", \"g-factor\", \"Radial Speed [m/s]\"");
*/

  //print variables of the Chamberlain exosphere model
  if (Exosphere::ChamberlainExosphere::ModelInitFlag==true) Exosphere::ChamberlainExosphere::PrintVariableList(fout);
//#endif
}


void Exosphere::Sampling::OutputDataFile::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double TotalMeasure=0.0,Measure=0.0;
  int i,iSource,nspec;
  char *SamplingBuffer,*CellNodeSamplingBuffer;
  double SourceRate[PIC::nTotalSpecies*(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_)];

  for (nspec=0;nspec<PIC::nTotalSpecies;nspec++) for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) SourceRate[nspec+iSource*PIC::nTotalSpecies]=0.0;

  CellNodeSamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

    Measure=InterpolationList[i]->Measure;
    if (Measure<=0.0) exit(__LINE__,__FILE__,"Error: non-positive cell volume is found");
    TotalMeasure+=Measure;

    if (SamplingDensityOffset[0]!=-1) for (nspec=0;nspec<PIC::nTotalSpecies;nspec++) {
      for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) SourceRate[nspec+iSource*PIC::nTotalSpecies]+=*(nspec+(double*)(SamplingBuffer+SamplingDensityOffset[iSource]));
    }
  }

  if (SamplingDensityOffset[0]!=-1) for (nspec=0;nspec<PIC::nTotalSpecies;nspec++) for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
    if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) SourceRate[nspec+iSource*PIC::nTotalSpecies]/=PIC::LastSampleLength*TotalMeasure;
    *(nspec+(double*)(CellNodeSamplingBuffer+SamplingDensityOffset[iSource]))=SourceRate[nspec+iSource*PIC::nTotalSpecies];
  }
}

void Exosphere::Sampling::OutputDataFile::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double t=0.0;
  char *SamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  //print position of the cell in the MSGR_HCI frame
  double x_HCI[3],*x_SO=CenterNode->GetX();

  x_HCI[0]=Exosphere::xObject_HCI[0]+
    (SO_to_HCI_TransformationMartix[0][0]*x_SO[0])+
    (SO_to_HCI_TransformationMartix[0][1]*x_SO[1])+
    (SO_to_HCI_TransformationMartix[0][2]*x_SO[2]);

  x_HCI[1]=Exosphere::xObject_HCI[1]+
    (SO_to_HCI_TransformationMartix[1][0]*x_SO[0])+
    (SO_to_HCI_TransformationMartix[1][1]*x_SO[1])+
    (SO_to_HCI_TransformationMartix[1][2]*x_SO[2]);

  x_HCI[2]=Exosphere::xObject_HCI[2]+
    (SO_to_HCI_TransformationMartix[2][0]*x_SO[0])+
    (SO_to_HCI_TransformationMartix[2][1]*x_SO[1])+
    (SO_to_HCI_TransformationMartix[2][2]*x_SO[2]);


  if (PIC::ThisThread==0) fprintf(fout,"%e  %e  %e ",x_HCI[0],x_HCI[1],x_HCI[2]);
#endif

  //output the sampled number densities
  if (SamplingDensityOffset[0]!=-1) for (int iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *(DataSetNumber+(double*)(SamplingBuffer+SamplingDensityOffset[iSource]));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }


/*
  //get the number density due to the impact evaporation source
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__ImpactVaporization_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);


#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__PhotonStimulatedDesorption_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__ThermalDesorption_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__SolarWindSputtering_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif
*/

  //user defined output in the data file of the exospehre model
#if _EXOSPHERE__USER_DEFINED_FILE_OUTPUT__MODE__ == _EXOSPHERE_SOURCE__ON_
  _EXOSPHERE__USER_DEFINED_FILE_OUTPUT__PRINT_DATA__(fout,DataSetNumber,pipe,CenterNodeThread,CenterNode);
#endif

/*  //mean particle radiation pressure acceleration and g-factor
  double gfactor=0.0,RadiationPressureAcceleration=0.0,vHeliocentric=0.0;

  if (pipe->ThisThread==CenterNodeThread) {
    if ((DataSetNumber==_NA_SPEC_)&&(x_SO[0]*x_SO[0]+x_SO[1]*x_SO[1]+x_SO[2]*x_SO[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))) {
      if ((x_SO[1]*x_SO[1]+x_SO[2]*x_SO[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x_SO[0]>0.0)) { //calculate the radiation pressure force
        double rHeliocentric;

        //calcualte velocity of the particle in the "Frozen SO Frame"; IMPORTANT: SO rotates only around z-direction!!!!
        double v_LOCAL_SO_FROZEN[3],LocalBulkVelocity_SO[3];

        CenterNode->GetBulkVelocity(LocalBulkVelocity_SO,_NA_SPEC_);

        v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+LocalBulkVelocity_SO[0]+
            Exosphere::RotationVector_SO_FROZEN[1]*x_SO[2]-Exosphere::RotationVector_SO_FROZEN[2]*x_SO[1];

        v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+LocalBulkVelocity_SO[1]-
            Exosphere::RotationVector_SO_FROZEN[0]*x_SO[2]+Exosphere::RotationVector_SO_FROZEN[2]*x_SO[0];

        v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+LocalBulkVelocity_SO[2]+
            Exosphere::RotationVector_SO_FROZEN[0]*x_SO[1]-Exosphere::RotationVector_SO_FROZEN[1]*x_SO[0];


        rHeliocentric=sqrt(pow(x_SO[0]-Exosphere::xObjectRadial,2)+(x_SO[1]*x_SO[1])+(x_SO[2]*x_SO[2]));
        vHeliocentric=(
            (v_LOCAL_SO_FROZEN[0]*(x_SO[0]-Exosphere::xObjectRadial))+
            (v_LOCAL_SO_FROZEN[1]*x_SO[1])+(v_LOCAL_SO_FROZEN[2]*x_SO[2]))/rHeliocentric;

        RadiationPressureAcceleration=SodiumRadiationPressureAcceleration__Combi_1997_icarus(vHeliocentric,rHeliocentric);
        gfactor+=SodiumGfactor__5891_58A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
        gfactor+=SodiumGfactor__5897_56A__Killen_2009_AJSS(vHeliocentric,rHeliocentric);
      }
    }
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) {
      pipe->recv(RadiationPressureAcceleration,CenterNodeThread);
      pipe->recv(gfactor,CenterNodeThread);
      pipe->recv(vHeliocentric,CenterNodeThread);
    }

    fprintf(fout,"%e  %e  %e  ",RadiationPressureAcceleration,gfactor,vHeliocentric);
  }
  else {
    pipe->send(RadiationPressureAcceleration);
    pipe->send(gfactor);
    pipe->send(vHeliocentric);
  }*/


  //Output parameters of the Chamberlain model
  if (Exosphere::ChamberlainExosphere::ModelInitFlag==true) Exosphere::ChamberlainExosphere::PrintData(fout,DataSetNumber,pipe,CenterNodeThread,CenterNode);


  //other soruce processes......
}
/*--------------------------------- Print Output File: END    --------------------------------------*/



/*--------------------------------- Source Processes: BEGIN  --------------------------------------*/
double Exosphere::SourceProcesses::totalProductionRate(int spec,int BoundaryElementType,void *BoundaryElement) {
  double res=0.0;

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
//  if (spec==_NA_SPEC_) {
  res+=Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(spec,BoundaryElementType,(cInternalSphericalData*)(BoundaryElement));
//  }
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
//  if (spec==_NA_SPEC_) {
  res+=Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
//  }
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
//  if (spec==_NA_SPEC_) {
  res+=Exosphere::SourceProcesses::ThermalDesorption::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
//  }
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
//  if (spec==_NA_SPEC_) {
  res+=Exosphere::SourceProcesses::SolarWindSputtering::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
//  }
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
//  if (spec==_NA_SPEC_) {
  res+=Exosphere::SourceProcesses::VerticalInjection::GetTotalProductionRate(spec,BoundaryElementType,(cInternalSphericalData*)(BoundaryElement));
//  }
#endif

  //add source rate due to the user defined source process
#if _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_ == _EXOSPHERE_SOURCE__ON_
  $MARKER:USER-DEFINED-TOTAL-SOURCE-RATE$
#endif

  return res;
}


long int Exosphere::SourceProcesses::InjectionBoundaryModel(int BoundaryElementType,void *BoundaryElement)  {
  int spec;
  long int res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=InjectionBoundaryModel(spec,BoundaryElementType,BoundaryElement);

  //inject the background plasma ions
  #if _EXOSPHERE__BACKGROUND_PLASMA_ION_INJECTION_ == _PIC_MODE_ON_
  res+=BackgroundPlasmaBoundaryIonInjection::ParticleInjection();
  #endif

  return res;
}

long int Exosphere::SourceProcesses::InjectionBoundaryModel(int spec,int BoundaryElementType,void *BoundaryElement) {
  cInternalSphericalData *Sphere=NULL;
//  void *SphereDataPointer=NULL;
  double *sphereX0=NULL,sphereRadius=0.0;

  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_SO_OBJECT[3],x_IAU_OBJECT[3],v_SO_OBJECT[3],v_IAU_OBJECT[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag;
  int SourceProcessID;

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;

  switch (BoundaryElementType) {
  case _INTERNAL_BOUNDARY_TYPE_SPHERE_ :
    Sphere=(cInternalSphericalData*)(BoundaryElement);
    Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

    break;
  }


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  ParticleWeight=0.0;
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
#else
  LocalTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  ModelParticlesInjectionRate=totalProductionRate(spec,BoundaryElementType,BoundaryElement)/ParticleWeight;

  if (ModelParticlesInjectionRate*ParticleWeight*LocalTimeStep<1.0E-10) return 0;

  #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (ModelParticlesInjectionRate*LocalTimeStep>nMaxInjectedParticles) {
    ParticleWeightCorrection=ModelParticlesInjectionRate*LocalTimeStep/nMaxInjectedParticles;
    ModelParticlesInjectionRate/=ParticleWeightCorrection;
  }
  #endif




  //calcualte probabilities of each source processes
  double TotalFlux,FluxSourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_]; //,ProbabilitySourceProcess[1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];
  int iSource;

  for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) FluxSourceProcess[iSource]=0.0; //,ProbabilitySourceProcess[iSource]=0.0;

/*  double TotalFlux,Flux_ImpactVaporization=0.0,Flux_PSD=0.0,Flux_TD=0.0,Flux_SW_Sputtering=0.0;
  double p,Probability_ImpactVaporization=0.0,Probability_PSD=0.0,Probability_TD=0.0,Probability_SW_Sputtering=0.0;*/

  TotalFlux=totalProductionRate(spec,BoundaryElementType,BoundaryElement);

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
//  Flux_ImpactVaporization=Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(_NA_SPEC_,SphereDataPointer);

  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]=Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate(spec,BoundaryElementType,(cInternalSphericalData*)(BoundaryElement));
  #endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
//  Flux_PSD=Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetTotalProductionRate(_NA_SPEC_,SphereDataPointer);
  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]=Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
  if (FluxSourceProcess[_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]>0.0) PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].Init(&spec);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
//  Flux_TD=Exosphere::SourceProcesses::ThermalDesorption::GetTotalProductionRate(_NA_SPEC_,SphereDataPointer);
  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]=Exosphere::SourceProcesses::ThermalDesorption::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
  if (FluxSourceProcess[_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]>0.0) ThermalDesorption::SurfaceInjectionDistribution[spec].Init(&spec);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
//  Flux_SW_Sputtering=Exosphere::SourceProcesses::SolarWindSputtering::GetTotalProductionRate(_NA_SPEC_,SphereDataPointer);
  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]=Exosphere::SourceProcesses::SolarWindSputtering::GetTotalProductionRate(spec,BoundaryElementType,BoundaryElement);
  if (FluxSourceProcess[_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]>0.0) SolarWindSputtering::SurfaceInjectionDistribution[spec].Init(&spec);
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_]=Exosphere::SourceProcesses::VerticalInjection::GetTotalProductionRate(spec,BoundaryElementType,(cInternalSphericalData*)(BoundaryElement));
#endif

#if _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_ == _EXOSPHERE_SOURCE__ON_
  //calculate the source rate due to user defined source functions
  $MARKER:CALCULATE-SOURCE-FLUX-WITH-USER-DEFINED-FUNCTIONS$
#endif

/*  for (iSource=0;iSource<1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_;iSource++) {
    ProbabilitySourceProcess[iSource]=FluxSourceProcess[iSource]/TotalFlux;
    if (iSource!=0) ProbabilitySourceProcess[iSource]+=ProbabilitySourceProcess[iSource-1];
  }*/

/*  Probability_ImpactVaporization=Flux_ImpactVaporization/TotalFlux;
  Probability_PSD=(Flux_PSD+Flux_ImpactVaporization)/TotalFlux;
  Probability_TD=(Flux_PSD+Flux_ImpactVaporization+Flux_TD)/TotalFlux;
  Probability_SW_Sputtering=(Flux_PSD+Flux_ImpactVaporization+Flux_TD+Flux_SW_Sputtering)/TotalFlux;

  //recalcualte the surface injection distributions
  if (Flux_PSD>0.0) PhotonStimulatedDesorption::SurfaceInjectionDistribution.Init();
  if (Flux_TD>0.0) ThermalDesorption::SurfaceInjectionDistribution.Init();
  if (Flux_SW_Sputtering>0.0) SolarWindSputtering::SurfaceInjectionDistribution.Init();*/


  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
    //determine the source process to generate a particle's properties
    do {
      SourceProcessID=(int)(rnd()*(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_));
    }
    while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());

    //the particle buffer used to set-up the new particle data
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

   //to satisfy the compiler and fit the while structure
   if (false) {}
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
   else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_) {
     flag=Exosphere::SourceProcesses::ImpactVaporization::GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);

//     SourceProcessID=_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_;
//     if (flag==true) Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
#endif
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
   else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_) {
     flag=Exosphere::SourceProcesses::PhotonStimulatedDesorption::GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);

//     SourceProcessID=_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_;
//     if (flag==true) Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
#endif
#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_) {
     flag=Exosphere::SourceProcesses::ThermalDesorption::GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);

//     SourceProcessID=_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_;
//     if (flag==true) Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
#endif
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_) {
     flag=Exosphere::SourceProcesses::SolarWindSputtering::GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);

//     SourceProcessID=_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_;
//     if (flag==true) Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
#endif
#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
   else if (SourceProcessID==_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_) {
     flag=Exosphere::SourceProcesses::VerticalInjection::GenerateParticleProperties(spec,(PIC::ParticleBuffer::byte*)tempParticleData,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,BoundaryElementType,BoundaryElement);

//     SourceProcessID=_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_;
//     if (flag==true) Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
#endif
#if _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_ == _EXOSPHERE_SOURCE__ON_
   //Add the user defined particle gineration
   $MARKER:GENERATE-PARTICLE-PROPERTIES-WITH-USER-DEFINED-FUNCTIONS$
#endif
   else {
     continue;
   }

   if (flag==false) continue;

   //sample the injection information
   Sampling::CalculatedSourceRate[spec][SourceProcessID]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   PIC::BC::nInjectedParticles[spec]+=1;
   PIC::BC::ParticleProductionRate[spec]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
   if (startNode->block->GetLocalTimeStep(_NA_SPEC_)/LocalTimeStep<rnd()) continue;
#endif


/*
cout << __FILE__ << "@" << __LINE__ << endl;
cout << __FILE__ << "@" << __LINE__ << "  " << x_SO_OBJECT[0] << "   " << x_SO_OBJECT[1] << "  " << x_SO_OBJECT[2] << endl;
cout << __FILE__ << "@" << __LINE__ << "  " << x_IAU_OBJECT[0] << "  " << x_IAU_OBJECT[1] << "  " << x_IAU_OBJECT[2] << endl;
*/

   //determine the surface element of the particle origin
   long int nZenithElement,nAzimuthalElement;
   int el;

#if _EXOSPHERE__SOURCE_PROCESSES__CONTROL_POSITIVE_VOLATILE_SURFACE_ABOUNDANCE_ == _PIC_MODE_ON_
   if (BoundaryElementType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
     Sphere->GetSurfaceElementProjectionIndex(x_IAU_OBJECT,nZenithElement,nAzimuthalElement);
     el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

     //check is the injection of the particle will not make the surface aboundance negative
     if (Source_DeplitSurfaceSpeciesAbundance_Flag[SourceProcessID]==true){
       if (Exosphere::Planet->SurfaceElementPopulation[spec][el]<Sphere->SurfaceElementDesorptionFluxUP[spec][el]+ParticleWeight*ParticleWeightCorrection) continue;
     }
   }
#endif

   //generate a particle
   PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
   PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
   PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

   #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
   PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
   #endif

   //apply condition of tracking the particle
   #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
   PIC::ParticleTracker::InitParticleID(tempParticleData);
   PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x_SO_OBJECT,v_SO_OBJECT,spec,tempParticleData);
   #endif

   //save the information od the particle origin: the particle origin will be sampled in SO coordinate frame
   if (BoundaryElementType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
     Sphere->GetSurfaceElementProjectionIndex(x_SO_OBJECT,nZenithElement,nAzimuthalElement);
     el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

     Exosphere::Planet->SampleSpeciesSurfaceSourceRate[spec][el][SourceProcessID]+=ParticleWeight*ParticleWeightCorrection/LocalTimeStep;

     //sample particle injection velocity
     Exosphere::Planet->SampleSpeciesSurfaceInjectionFlux[spec][el]+=ParticleWeight*ParticleWeightCorrection;
     Exosphere::Planet->SampleInjectedFluxBulkSpeed[spec][el]+=ParticleWeight*ParticleWeightCorrection*sqrt(pow(v_IAU_OBJECT[0],2)+pow(v_IAU_OBJECT[1],2)+pow(v_IAU_OBJECT[2],2));
   }
   else el=0;

   Sampling::SetParticleSourceID(SourceProcessID,(PIC::ParticleBuffer::byte*)tempParticleData);
   Sampling::SetParicleOriginSurfaceElementNumber(el,(PIC::ParticleBuffer::byte*)tempParticleData);

   newParticle=PIC::ParticleBuffer::GetNewParticle();
   newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
   memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

   nInjectedParticles++;

   //for the secondary source processes accout for the decrease of the surface density
   if (BoundaryElementType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
     if (Source_DeplitSurfaceSpeciesAbundance_Flag[SourceProcessID]==true) {
       Sphere->GetSurfaceElementProjectionIndex(x_IAU_OBJECT,nZenithElement,nAzimuthalElement);
       el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

       Sphere->SurfaceElementDesorptionFluxUP[spec][el]+=ParticleWeight*ParticleWeightCorrection;
     }
   }




   //inject the particle into the system
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
   if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
   if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is not defined");
#endif

   _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
 }

 return nInjectedParticles;
}







/*=============================== Source Processes: END  ===========================================*/



/*=============================== INTERACTION WITH THE SURFACE: BEGIN  ===========================================*/
int Exosphere::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x_SO_OBJECT,double *v_SO_OBJECT,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
  double radiusSphere,*x0Sphere,lNorm[3],rNorm,lVel[3],rVel,c;
  cInternalSphericalData *Sphere;
//  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  int idim;
  double vi,vt,vf,v_LOCAL_SO_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],x_LOCAL_IAU_OBJECT[3],cosSubsolarAngle,SurfaceTemp,beta;
  SpiceDouble xform[6][6];

  Sphere=(cInternalSphericalData*)SphereDataPointer;
//  startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

  memcpy(v_LOCAL_SO_OBJECT,v_SO_OBJECT,3*sizeof(double));
  memcpy(x_LOCAL_SO_OBJECT,x_SO_OBJECT,3*sizeof(double));

  //convert the position vector from SO_OBJECT to IAU_OBJECT coordinate frames
  memcpy(xform,OrbitalMotion::SO_to_IAU_TransformationMartix,36*sizeof(double));

  x_LOCAL_IAU_OBJECT[0]=xform[0][0]*x_LOCAL_SO_OBJECT[0]+xform[0][1]*x_LOCAL_SO_OBJECT[1]+xform[0][2]*x_LOCAL_SO_OBJECT[2];
  x_LOCAL_IAU_OBJECT[1]=xform[1][0]*x_LOCAL_SO_OBJECT[0]+xform[1][1]*x_LOCAL_SO_OBJECT[1]+xform[1][2]*x_LOCAL_SO_OBJECT[2];
  x_LOCAL_IAU_OBJECT[2]=xform[2][0]*x_LOCAL_SO_OBJECT[0]+xform[2][1]*x_LOCAL_SO_OBJECT[1]+xform[2][2]*x_LOCAL_SO_OBJECT[2];

  //get local surface temperature
  cosSubsolarAngle=Exosphere::OrbitalMotion::GetCosineSubsolarAngle(x_LOCAL_SO_OBJECT);
  SurfaceTemp=Exosphere::GetSurfaceTemeprature(cosSubsolarAngle,x_LOCAL_SO_OBJECT);


  //sample parameters of the back flux: speed is calcalate in IAU (relative to the planet) but the flux is samplined in SO (one axis is alwais directed to the Sun)
  //convert the velocity vector from SO_OBJECT to IAU_OBJECT coordinate frames
  v_LOCAL_IAU_OBJECT[0]=xform[3][0]*x_LOCAL_SO_OBJECT[0]+xform[3][1]*x_LOCAL_SO_OBJECT[1]+xform[3][2]*x_LOCAL_SO_OBJECT[2]+
      xform[3][3]*v_LOCAL_SO_OBJECT[0]+xform[3][4]*v_LOCAL_SO_OBJECT[1]+xform[3][5]*v_LOCAL_SO_OBJECT[2];

  v_LOCAL_IAU_OBJECT[1]=xform[4][0]*x_LOCAL_SO_OBJECT[0]+xform[4][1]*x_LOCAL_SO_OBJECT[1]+xform[4][2]*x_LOCAL_SO_OBJECT[2]+
      xform[4][3]*v_LOCAL_SO_OBJECT[0]+xform[4][4]*v_LOCAL_SO_OBJECT[1]+xform[4][5]*v_LOCAL_SO_OBJECT[2];

  v_LOCAL_IAU_OBJECT[2]=xform[5][0]*x_LOCAL_SO_OBJECT[0]+xform[5][1]*x_LOCAL_SO_OBJECT[1]+xform[5][2]*x_LOCAL_SO_OBJECT[2]+
      xform[5][3]*v_LOCAL_SO_OBJECT[0]+xform[5][4]*v_LOCAL_SO_OBJECT[1]+xform[5][5]*v_LOCAL_SO_OBJECT[2];


  vi=sqrt(v_LOCAL_IAU_OBJECT[0]*v_LOCAL_IAU_OBJECT[0]+v_LOCAL_IAU_OBJECT[1]*v_LOCAL_IAU_OBJECT[1]+v_LOCAL_IAU_OBJECT[2]*v_LOCAL_IAU_OBJECT[2]);

  long int nZenithElement,nAzimuthalElement;
  int el;
  double ParticleWeight;

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
#else
  ParticleWeight=0.0;
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
#endif


  Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_SO_OBJECT,nZenithElement,nAzimuthalElement);
  el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

  Sphere->SampleSpeciesSurfaceReturnFlux[spec][el]+=ParticleWeight;
  Sphere->SampleReturnFluxBulkSpeed[spec][el]+=vi*ParticleWeight;

  //total return flux
  Exosphere::Sampling::TotalPlanetReturnFlux[spec]+=ParticleWeight;

  //sample returned flux on the night side of the planet
  if (x_LOCAL_SO_OBJECT[0]<0.0) {
    //the night size
    int id;

    id=Sampling::GetParticleSourceID(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
    Exosphere::Sampling::PlanetNightSideReturnFlux[spec][id]+=ParticleWeight;
  }

  //check if the particle sticks to the surface of the planet
  double ReemissionParticleFraction=1.0;

  if (rnd()<StickingProbability(spec,ReemissionParticleFraction,SurfaceTemp)) {
    //the particle is abserbed by the surface
    Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_IAU_OBJECT,nZenithElement,nAzimuthalElement);
    el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

    Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el]+=ReemissionParticleFraction*ParticleWeight;
    Exosphere::Sampling::PlanetSurfaceStickingRate[spec]+=ParticleWeight;

//    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }









  //distribute velocity with the local surface temperature
  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*SurfaceTemp));
  for (vt=0.0,idim=0;idim<3;idim++) vt+=pow(sqrt(-log(rnd()))/beta*cos(PiTimes2*rnd()),2);

  vt=sqrt(vt);
  vf=AccomodationCoefficient[spec]*vt+(1.0-AccomodationCoefficient[spec])*vi;

  //determine the normal vector at the intersection point
  Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

  for (rNorm=0.0,rVel=0.0,c=0.0,idim=0;idim<DIM;idim++) {
    lNorm[idim]=x_LOCAL_IAU_OBJECT[idim]-x0Sphere[idim];
    rNorm+=pow(lNorm[idim],2);

    lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
    rVel+=pow(lVel[idim],2);

    c+=lNorm[idim]*lVel[idim];
  }

  rVel=sqrt(rVel);

  if (c>0.0) {
    //the distributed velocity vector is directed into the domain
    for (rVel=vf/rVel,idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=lVel[idim]*rVel;
  }
  else {
    //the distributed velocity vector is directed into the planet -> redirect it
    double c1=0.0;

    for (rNorm=sqrt(rNorm),idim=0;idim<3;idim++) {
      lVel[idim]/=rVel;
      lNorm[idim]/=rNorm;

      c1+=lVel[idim]*lNorm[idim];
    }

    for (idim=0;idim<3;idim++) {
      v_LOCAL_IAU_OBJECT[idim]=vf*(lVel[idim]-2.0*c1*lNorm[idim]);
    }
  }

  //convert velocity vector from IAU_OBJECT -> SO_OBJECT coordinate fraim
  memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

  v_LOCAL_SO_OBJECT[0]=xform[3][0]*x_LOCAL_IAU_OBJECT[0]+xform[3][1]*x_LOCAL_IAU_OBJECT[1]+xform[3][2]*x_LOCAL_IAU_OBJECT[2]+
      xform[3][3]*v_LOCAL_IAU_OBJECT[0]+xform[3][4]*v_LOCAL_IAU_OBJECT[1]+xform[3][5]*v_LOCAL_IAU_OBJECT[2];

  v_LOCAL_SO_OBJECT[1]=xform[4][0]*x_LOCAL_IAU_OBJECT[0]+xform[4][1]*x_LOCAL_IAU_OBJECT[1]+xform[4][2]*x_LOCAL_IAU_OBJECT[2]+
      xform[4][3]*v_LOCAL_IAU_OBJECT[0]+xform[4][4]*v_LOCAL_IAU_OBJECT[1]+xform[4][5]*v_LOCAL_IAU_OBJECT[2];

  v_LOCAL_SO_OBJECT[2]=xform[5][0]*x_LOCAL_IAU_OBJECT[0]+xform[5][1]*x_LOCAL_IAU_OBJECT[1]+xform[5][2]*x_LOCAL_IAU_OBJECT[2]+
      xform[5][3]*v_LOCAL_IAU_OBJECT[0]+xform[5][4]*v_LOCAL_IAU_OBJECT[1]+xform[5][5]*v_LOCAL_IAU_OBJECT[2];


  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

/*=============================== INTERACTION WITH THE SURFACE: END  ===========================================*/



/*=================================== Output surface parameters =================================================*/
void Exosphere::Sampling::OutputSurfaceDataFile::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Surface Temperature [K]\", \"Total Flux Down [s^{-1} m^{-2}]\", \"Total Flux Up [s^{-1} m^{-2}]\", \"Surface Content [m^{-2}]\", \"Returned particles' bulk speed [m/s]\"");
  fprintf(fout,", \"Injected particles' bulk speed [m/s]\", \"Sodium Sticking Coefficient\"");

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
   fprintf(fout,", \"Impact Vaporization Flux [s^{-1} m^{-2}]\"");
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
   fprintf(fout,", \"Photon Stimulated Desorption Flux [s^{-1} m^{-2}]\"");
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   fprintf(fout,", \"Thermal Desorption Flux [s^{-1} m^{-2}]\"");
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   fprintf(fout,", \"Solar Wind Sputtering Flux [s^{-1} m^{-2}]\", \"Solar Wind Flux [s^{-1} m^{-2}]\"");
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
   fprintf(fout,", \"Vertical Injection Flux [s^{-1} m^{-2}]\"");
#endif
}

void Exosphere::Sampling::OutputSurfaceDataFile::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double t,TotalFluxDown=0.0,TotalFluxUp=0.0,SurfaceContent=0.0,BulkSpeedDown=0.0,BulkSpeedUp=0.0,SampleSpeciesSurfaceInjectionFlux=0.0;


#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
   double FluxIV=0.0;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
   double FluxPDS=0.0;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   double FluxTD=0.0;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   double FluxSW=0.0,SolarWindIncidentFlux=0.0;
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
   double FluxVI=0.0;
#endif

#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
#endif


  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
    InterpolationCoefficient=Sphere->GetSurfaceElementArea(nSurfaceElement);

    BulkSpeedDown+=Sphere->SampleReturnFluxBulkSpeed[spec][nSurfaceElement];
    TotalFluxDown+=Sphere->SampleSpeciesSurfaceReturnFlux[spec][nSurfaceElement];
    SurfaceContent+=Sphere->SampleSpeciesSurfaceAreaDensity[spec][nSurfaceElement]*InterpolationCoefficient;

    BulkSpeedUp+=Sphere->SampleInjectedFluxBulkSpeed[spec][nSurfaceElement];
    SampleSpeciesSurfaceInjectionFlux+=Sphere->SampleSpeciesSurfaceInjectionFlux[spec][nSurfaceElement];


#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]/PIC::LastSampleLength : 0.0;
    FluxIV+=t;
    TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]/PIC::LastSampleLength : 0.0;
    FluxPDS+=t;
    TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]/PIC::LastSampleLength : 0.0;
   FluxTD+=t;
   TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]/PIC::LastSampleLength : 0.0;
   FluxSW+=t;
   TotalFluxUp+=t;

   SolarWindIncidentFlux+=Sphere->SolarWindSurfaceFlux[nSurfaceElement]*InterpolationCoefficient;
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_]/PIC::LastSampleLength : 0.0;
    FluxVI+=t;
    TotalFluxUp+=t;
#endif

    InterpolationNormalization+=InterpolationCoefficient;
  }



  if (ThisThread==0)  {
    //collect sampled data from all processors
    for (int thread=1;thread<nTotalThreads;thread++) {
      TotalFluxDown+=pipe->recv<double>(thread);
//      SurfaceContent+=pipe->recv<double>(thread);  All processors have the same distribution of surface content map
      TotalFluxUp+=pipe->recv<double>(thread);
      BulkSpeedDown+=pipe->recv<double>(thread);

      BulkSpeedUp+=pipe->recv<double>(thread);
      SampleSpeciesSurfaceInjectionFlux+=pipe->recv<double>(thread);

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
      FluxIV+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
      FluxPDS+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      FluxTD+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      FluxSW+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
      FluxVI+=pipe->recv<double>(thread);
#endif
    }

    if (PIC::LastSampleLength!=0) {
      if (TotalFluxDown>0.0) BulkSpeedDown/=TotalFluxDown;
      if (SampleSpeciesSurfaceInjectionFlux>0.0) BulkSpeedUp/=SampleSpeciesSurfaceInjectionFlux;

      TotalFluxDown/=PIC::LastSampleLength*PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
      SurfaceContent/=PIC::LastSampleLength;
    }

    //Print Surface Temparature
    double norm[3],CosSubSolarAngle,SurfaceTemperature,x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3];

    Exosphere::Planet->GetSurfaceNormal(norm,nZenithPoint,nAzimuthalPoint);

    x_LOCAL_IAU_OBJECT[0]=_RADIUS_(_TARGET_)*norm[0];
    x_LOCAL_IAU_OBJECT[1]=_RADIUS_(_TARGET_)*norm[1];
    x_LOCAL_IAU_OBJECT[2]=_RADIUS_(_TARGET_)*norm[2];

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

    CosSubSolarAngle=norm[0];
    SurfaceTemperature=GetSurfaceTemeprature(CosSubSolarAngle,x_LOCAL_SO_OBJECT);
    fprintf(fout," %e",SurfaceTemperature);

    //Print Sampled Surface data
    double ReemissionParticleFraction;

    fprintf(fout," %e %e %e %e %e %e ",TotalFluxDown/InterpolationNormalization,TotalFluxUp/InterpolationNormalization,SurfaceContent/InterpolationNormalization,BulkSpeedDown,BulkSpeedUp,  \
        SurfaceInteraction::StickingProbability(spec,ReemissionParticleFraction,SurfaceTemperature));

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    fprintf(fout," %e ",FluxIV/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    fprintf(fout," %e ",FluxPDS/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      fprintf(fout," %e ",FluxTD/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      fprintf(fout," %e %e ",FluxSW/InterpolationNormalization,SolarWindIncidentFlux/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
    fprintf(fout," %e ",FluxVI/InterpolationNormalization);
#endif
  }
  else {
    pipe->send(TotalFluxDown);
//    pipe->send(SurfaceContent);     All processors have the same distribution of surface content map
    pipe->send(TotalFluxUp);
    pipe->send(BulkSpeedDown);

    pipe->send(BulkSpeedUp);
    pipe->send(SampleSpeciesSurfaceInjectionFlux);

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxIV);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxPDS);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxTD);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxSW);
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxVI);
#endif
  }
}

void Exosphere::Sampling::OutputSurfaceDataFile::flushCollectingSamplingBuffer(cInternalSphericalData* Sphere) {
  int s,el,i;
  int maxSurfaceSourceID=-1;

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  if (maxSurfaceSourceID<_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_) maxSurfaceSourceID=_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
  if (maxSurfaceSourceID<_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_) maxSurfaceSourceID=_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
  if (maxSurfaceSourceID<_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_) maxSurfaceSourceID=_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  if (maxSurfaceSourceID<_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_) maxSurfaceSourceID=_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_;
#endif

#if _EXOSPHERE_SOURCE__VERTICAL_INJECTION_ == _EXOSPHERE_SOURCE__ON_
  if (maxSurfaceSourceID<_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_) maxSurfaceSourceID=_EXOSPHERE_SOURCE__ID__VERTICAL_INJECTION_;
#endif

  if (Sphere!=NULL) {
    for (s=0;s<PIC::nTotalSpecies;s++) for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      for (i=0;i<maxSurfaceSourceID+1;i++) Sphere->SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;

      Sphere->SampleSpeciesSurfaceAreaDensity[s][el]=0.0;
      Sphere->SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
      Sphere->SampleReturnFluxBulkSpeed[s][el]=0.0;
    }

    PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
  }
}


void Exosphere::Sampling::OutputSurfaceDataFile::PrintTitle(FILE* fout) {
  fprintf(fout,"TITLE=\"SurfaceData:  TAA=%e [deg]\"",Exosphere::OrbitalMotion::GetTAA(Exosphere::OrbitalMotion::et)/Pi*180.0);
}

