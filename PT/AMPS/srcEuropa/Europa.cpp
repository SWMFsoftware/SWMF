/*
 * Europa.cpp
 *
 *  Created on: Feb 13, 2012
 *      Author: vtenishe
 */

//$Id$



#include "pic.h"
#include "Europa.h"

//sputtering by the ions
#define _ION_SPUTTERING_MODE_  _PIC_MODE_ON_

//the object name and the names of the frames
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="Europa";
char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_EUROPA";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EPHIOD";


double Europa::TotalInjectionRate=0.0;

double Europa::swE_Typical[3]={0.0,0.0,0.0};
double Europa::xEuropa[3]={0.0,0.0,0.0},Europa::vEuropa[3]={0.0,0.0,0.0};
double Europa::xEarth[3]={0.0,0.0,0.0},Europa::vEarth[3]={0.0,0.0,0.0};
double Europa::vEuropaRadial=0.0,Europa::xEuropaRadial=0.0;


//sampling offsets
int Europa::Sampling::SamplingDensity__ImpactVaporization_Offset=-1;
int Europa::Sampling::SamplingDensity__PhotonStimulatedDesorption_Offset=-1;
int Europa::Sampling::SamplingDensity__ThermalDesorption_Offset=-1;
int Europa::Sampling::SamplingDensity__SolarWindSputtering_Offset=-1;
int Europa::Sampling::CellSamplingDataOffset=-1;
double **Europa::Sampling::PlanetNightSideReturnFlux=NULL;

//the field in the particle's data that keeps the id of the source process due to which the particle has beed produced
long int Europa::Sampling::ParticleData_SourceProcessID_Offset=-1;
long int Europa::Sampling::ParticleData_OriginSurfaceElementNumber_Offset=-1;

//the total number of source processes
int Europa::nTotalSourceProcesses=0;

//the sphere that represents the planet
cInternalSphericalData *Europa::Planet=NULL;

/*//the total source rate values for specific source processes
double Europa::SourceProcesses::PhotonStimulatedDesorption::SourceRate=0.0,Europa::SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate=0.0;
double Europa::SourceProcesses::ThermalDesorption::SourceRate=0.0,Europa::SourceProcesses::ThermalDesorption::maxLocalSourceRate=0.0;
double Europa::SourceProcesses::SolarWindSputtering::SourceRate=0.0,Europa::SourceProcesses::SolarWindSputtering::maxLocalSourceRate=0.0;

//evaluate numerically the source rate
double Europa::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate=0.0;*/



double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {
  exit(__LINE__,__FILE__,"not implemented");

  return 0.0;
}

//request the sampling and particle's field data
int Europa::Sampling::RequestSamplingData(int offset) {
  int SamplingLength=0;

  if (CellSamplingDataOffset!=-1) exit(__LINE__,__FILE__,"Error: second request for the sampling data");

  CellSamplingDataOffset=offset;


#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  SamplingDensity__ImpactVaporization_Offset=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  Europa::Sampling::SamplingDensity__PhotonStimulatedDesorption_Offset=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  SamplingDensity__ThermalDesorption_Offset=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  SamplingDensity__SolarWindSputtering_Offset=CellSamplingDataOffset+SamplingLength;
  SamplingLength+=sizeof(double)*PIC::nTotalSpecies;
#endif


  //reserve the particle fields
  PIC::ParticleBuffer::RequestDataStorage(ParticleData_SourceProcessID_Offset,2*sizeof(int));
  ParticleData_OriginSurfaceElementNumber_Offset=ParticleData_SourceProcessID_Offset+sizeof(int);


  return SamplingLength;
}

//init the model
void Europa::Init_BeforeParser() {
  int idim;

  //init the energy distrivution model of the injection particles
/*  SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution. Init(SourceProcesses::PhotonStimulatedDesorption::minInjectionEnergy, \
      SourceProcesses::PhotonStimulatedDesorption::maxInjectionEnergy,SourceProcesses::PhotonStimulatedDesorption::EnergyDistributionFunction, \
      _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_);

  SourceProcesses::SolarWindSputtering::EnergyDistribution.Init(SourceProcesses::SolarWindSputtering::minInjectionEnergy, \
      SourceProcesses::SolarWindSputtering::maxInjectionEnergy,SourceProcesses::SolarWindSputtering::EnergyDistributionFunction, \
      _SINGLE_VARIABLE_DISTRIBUTION__INTERPOLATION_MODE__LINEAR_);*/


  //calculate the typical value of the motional electrical field
  swE_Typical[0]=-(swVelocity_Typical[1]*swB_Typical[2]-swVelocity_Typical[2]*swB_Typical[1]);
  swE_Typical[1]=+(swVelocity_Typical[0]*swB_Typical[2]-swVelocity_Typical[2]*swB_Typical[0]);
  swE_Typical[2]=-(swVelocity_Typical[0]*swB_Typical[1]-swVelocity_Typical[1]*swB_Typical[0]);


  //set the pre-processor into the ICES model
  PIC::CPLR::ICES::SWMFdataPreProcessor=SWMFdataPreProcessor;

  //set up the model sampling procedure
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SampleModelData,Sampling::OutputSampledModelData);

  //furnish the SPICE kernels
  char str[_MAX_STRING_LENGTH_PIC_];

  for (int nKernel=0;nKernel<nFurnishedSPICEkernels;nKernel++) {
    sprintf(str,"%s/%s",SPICE_Kernels_PATH,SPICE_Kernels[nKernel]);
    furnsh_c(str);
  }

  //Get the initial parameters of Europa orbit
  SpiceDouble state[6];

  utc2et_c(SimulationStartTimeString,&OrbitalMotion::et);

  //get initial parameters of Europa's orbit
  spkezr_c("Europa",OrbitalMotion::et,"GALL_EPHIOD","none","Jupiter",state,&OrbitalMotion::lt);

  for (idim=0,xEuropaRadial=0.0;idim<3;idim++) {
    xEuropa[idim]=state[idim]*1.0E3,vEuropa[idim]=state[idim+3]*1.0E3;
    xEuropaRadial+=pow(xEuropa[idim],2);
  }

  xEuropaRadial=sqrt(xEuropaRadial);
  vEuropaRadial=(xEuropa[0]*vEuropa[0]+xEuropa[1]*vEuropa[1]+xEuropa[2]*vEuropa[2])/xEuropaRadial;

  // get initial position of GALILEO for line-of sight
  spkezr_c("GALILEO ORBITER",OrbitalMotion::et,"GALL_EPHIOD","none","Europa",state,&OrbitalMotion::lt);
  for (idim=0;idim<3;idim++) xEarth[idim]=state[idim]*1.0E3,vEarth[idim]=state[idim+3]*1.0E3;



  // calculate the total number of source processes
#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__OFF_
//do nothing
#else
  exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  ++nTotalSourceProcesses;
#elif _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__OFF_
  //do nothing
#else
    exit(__LINE__,__FILE__,"Error: the option is not found");
#endif

  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestSamplingData.push_back(Sampling::RequestSamplingData);

  //print out of the output file
  //PIC::Mesh::PrintVariableListCenterNode.push_back(Sampling::OutputDataFile::PrintVariableList);
  //PIC::Mesh::PrintDataCenterNode.push_back(Sampling::OutputDataFile::PrintData);
  //PIC::Mesh::InterpolateCenterNode.push_back(Sampling::OutputDataFile::Interpolate);

  //init the energy distribution function for specific injection processes
#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.Init();

  if (PIC::ThisThread==0) {
    SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintCumulativeDistributionFunction("CumulativeEnergyDistribution-PSD.dat");
    SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintDistributionFunction("EnergyDistribution-PSD.dat");
  }
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  SourceProcesses::SolarWindSputtering::EnergyDistribution.Init();

  if (PIC::ThisThread==0) {
    SourceProcesses::SolarWindSputtering::EnergyDistribution.fPrintCumulativeDistributionFunction("CumulativeEnergyDistribution-SWS.dat");
    SourceProcesses::SolarWindSputtering::EnergyDistribution.fPrintDistributionFunction("EnergyDistribution-SWS.dat");
  }
#endif

}

void Europa::Init_AfterParser() {

  //init the sampling buffer of the returned flux on the night side of the planet
  int offset,s,i;

  Sampling::PlanetNightSideReturnFlux=new double *[PIC::nTotalSpecies];
  Sampling::PlanetNightSideReturnFlux[0]=new double [PIC::nTotalSpecies*(_EUROPA_SOURCE_MAX_ID_VALUE_+1)];

  for (s=0,offset=0;s<PIC::nTotalSpecies;s++) {
    Sampling::PlanetNightSideReturnFlux[s]=Sampling::PlanetNightSideReturnFlux[0]+offset;
    offset+=_EUROPA_SOURCE_MAX_ID_VALUE_+1;

    for (i=0;i<_EUROPA_SOURCE_MAX_ID_VALUE_+1;i++) Sampling::PlanetNightSideReturnFlux[s][i]=0.0;
  }

}



//ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
void Europa::SWMFdataPreProcessor(double *x,PIC::CPLR::ICES::cDataNodeSWMF& data) {
  int i;

  if (data.status!=_PIC_ICES__STATUS_OK_) {
    for (i=0;i<3;i++) {
      data.B[i]=swB_Typical[i];
      data.E[i]=swE_Typical[i];
      data.swVel[i]=swVelocity_Typical[i];
    }

    data.swTemperature=swTemperature_Typical;
    data.swNumberDensity=swNumberDensity_Typical;

    //p=2*n*k*T assuming quasi neutrality ni=ne=n and Ti=Te=T
    data.swPressure=2.0*data.swNumberDensity*Kbol*data.swTemperature;
  }
}

//calculate the sodium column density and plot
void Europa::SodiumCoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd;
  double NumberDensity;

  for (i=0;i<resLength;i++) res[i]=0.0;

  //get the local density number
  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(_O2_SPEC_);
  res[0]=NumberDensity;

  //get the local scattering
  if (resLength>1) {
    double BulkVelocity[3],r;

    node->block->GetCenterNode(nd)->GetBulkVelocity(BulkVelocity,_O2_SPEC_);
    r=sqrt(pow(xEuropaRadial-x[0],2)+(x[1]*x[1])+(x[2]*x[2]));

    res[1]=1.0E-10*NumberDensity*SodiumGfactor__5891_58A__Killen_2009_AJSS(BulkVelocity[0],r);
    res[2]=1.0E-10*NumberDensity*SodiumGfactor__5897_56A__Killen_2009_AJSS(BulkVelocity[0],r);
  }

}

void Europa::ColumnDensityIntegration_Tail(char *fname) {
  double xEarth_new[3]={0.0,0.0,0.0};

  const int nPoints=500;
  const double IntegrationRangeBegin=-50.0E6;
  const double IntegrationRangeEnd=50.0E6;
  const double IntegrationStep=(IntegrationRangeEnd-IntegrationRangeBegin)/(nPoints-1);


  //find position of Galileo
  SpiceDouble State[6],lt;

  spkezr_c("GALILEO ORBITER",Europa::OrbitalMotion::et,"GALL_EPHIOD","none","Europa",State,&lt);

  xEarth_new[0]=State[0]*1.0E3;
  xEarth_new[1]=State[1]*1.0E3;
  xEarth_new[2]=State[2]*1.0E3;

  //open the output file
  FILE *fout=NULL;
  int npoint;

  if (PIC::ThisThread==0) {
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Distance from the planet\", \"Total Column Density\"\n");
  }

  for (npoint=0;npoint<nPoints;npoint++) {
     double l[3];
     double ColumnDensity;

     l[0]=-(IntegrationRangeBegin+npoint*IntegrationStep)-xEarth_new[0];
     l[1]=-xEarth_new[1];
     l[2]=-xEarth_new[2];

     PIC::ColumnIntegration::GetCoulumnIntegral(&ColumnDensity,1,xEarth_new,l,SodiumCoulumnDensityIntegrant);

     if (PIC::ThisThread==0) fprintf(fout,"%e   %e\n",IntegrationRangeBegin+npoint*IntegrationStep,ColumnDensity);
  }

  if (PIC::ThisThread==0) fclose(fout);
}

void Europa::ColumnDensityIntegration_Map(char *fname) {
  double l[3],xEarth_new[3]={0.0,0.0,0.0};

  //find position of Europa as seen from Galileo
  SpiceDouble State[6],lt;

  spkezr_c("GALILEO ORBITER",Europa::OrbitalMotion::et,"GALL_EPHIOD","none","Europa",State,&lt);

  xEarth_new[0]=State[0]*1.0E3;
  xEarth_new[1]=State[1]*1.0E3;
  xEarth_new[2]=State[2]*1.0E3;


  const int nPointsSun=1000;
  const double minPhiX=-0.003,maxPhiX=0.003;
  const double minPhiZ=-0.0005,maxPhiZ=0.0005;

  const double dPhi=(maxPhiX-minPhiX)/(nPointsSun-1);
  const int nPointsZ=1+(int)((maxPhiZ-minPhiZ)/dPhi);

  //open the output file
  FILE *fout=NULL;
  int iZ,iX;
  double r,PhiX,PhiZ,StateVector[3];


  if (PIC::ThisThread==0) {
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Angle In Ecliptic Plane\", \"Angle Out of Ecpliptic Plane\", \"Column Density\", \"Intesity (5891.58A)\", \"Intesity (5897.56A)\" \n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nPointsZ,nPointsSun);
  }

  r=sqrt(xEarth_new[0]*xEarth_new[0]+xEarth_new[1]*xEarth_new[1]+xEarth_new[2]*xEarth_new[2]);

  for (iX=0;iX<nPointsSun;iX++) {
    //determine the X,Y-components of the direction vector
    PhiX=minPhiX+dPhi*iX;

    //rotate the vector
    l[0]=-(cos(PhiX)*xEarth_new[0]-sin(PhiX)*xEarth_new[1]);
    l[1]=-(sin(PhiX)*xEarth_new[0]+cos(PhiX)*xEarth_new[1]);

    for (iZ=0;iZ<nPointsZ;iZ++) {
      //determine the Z-component of the direction vector
      PhiZ=minPhiZ+dPhi*iZ;
      l[2]=r*tan(PhiZ)-xEarth_new[2];

      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,3,xEarth_new,l,SodiumCoulumnDensityIntegrant);

      if (PIC::ThisThread==0) fprintf(fout,"%e   %e   %e   %e  %e\n",PhiX,PhiZ,StateVector[0],StateVector[1],StateVector[2]);
    }

  }

  if (PIC::ThisThread==0) fclose(fout);
}

//=====================================================================================================
//empty sampling function
void Europa::Sampling::SampleModelData() {

  if (Europa::Planet==NULL) return;


  //sample sodium surface content
  int nZenithSurfaceElements,nAzimuthalSurfaceElements,spec;
  long int iZenith_GALL_EPHIOD,iAzimuth_GALL_EPHIOD,el_GALL_EPHIOD;
  long int iZenith_IAU,iAzimuth_IAU,el_IAU;
  double rSurfaceElement_IAU[3],rSurfaceElement_GALL_EPHIOD[3],SurfaceArea_IAU;
  SpiceDouble xform[6][6];

  nZenithSurfaceElements=Planet->nZenithSurfaceElements;
  nAzimuthalSurfaceElements=Planet->nAzimuthalSurfaceElements;

  memcpy(xform,OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix,36*sizeof(double));

  for (iZenith_GALL_EPHIOD=0;iZenith_GALL_EPHIOD<nZenithSurfaceElements;iZenith_GALL_EPHIOD++) for (iAzimuth_GALL_EPHIOD=0;iAzimuth_GALL_EPHIOD<nAzimuthalSurfaceElements;iAzimuth_GALL_EPHIOD++)  {
    //generate position on the sphere in the coordinate frame GALL_EPHIOD (x-axis is directed to the Sun)
    Planet->GetSurfaceElementMiddlePoint(rSurfaceElement_GALL_EPHIOD,iZenith_GALL_EPHIOD,iAzimuth_GALL_EPHIOD);
    el_GALL_EPHIOD=Planet->GetLocalSurfaceElementNumber(iZenith_GALL_EPHIOD,iAzimuth_GALL_EPHIOD);

    //conver the position vector into the frane related to the planet
    rSurfaceElement_IAU[0]=xform[0][0]*rSurfaceElement_GALL_EPHIOD[0]+xform[0][1]*rSurfaceElement_GALL_EPHIOD[1]+xform[0][2]*rSurfaceElement_GALL_EPHIOD[2];
    rSurfaceElement_IAU[1]=xform[1][0]*rSurfaceElement_GALL_EPHIOD[0]+xform[1][1]*rSurfaceElement_GALL_EPHIOD[1]+xform[1][2]*rSurfaceElement_GALL_EPHIOD[2];
    rSurfaceElement_IAU[2]=xform[2][0]*rSurfaceElement_GALL_EPHIOD[0]+xform[2][1]*rSurfaceElement_GALL_EPHIOD[1]+xform[2][2]*rSurfaceElement_GALL_EPHIOD[2];

    //get the surface element of the sampled point in the frame related to the planet
    Planet->GetSurfaceElementProjectionIndex(rSurfaceElement_IAU,iZenith_IAU,iAzimuth_IAU);
    SurfaceArea_IAU=Planet->GetSurfaceElementArea(iZenith_IAU,iAzimuth_IAU);


    el_IAU=Planet->GetLocalSurfaceElementNumber(iZenith_IAU,iAzimuth_IAU);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      double SpecieSurfaceDensity=0.0;


//      exit(__LINE__,__FILE__,"not implemented");

      SpecieSurfaceDensity=(spec==_O2_SPEC_) ? Planet->SampleSpeciesSurfaceAreaDensity[spec][el_IAU]/SurfaceArea_IAU : 0.0;
      Planet->SurfaceElementPopulation[spec][el_GALL_EPHIOD]+=SpecieSurfaceDensity;
    }

  }
}


//output the column integrated density
void Europa::Sampling::OutputSampledModelData(int DataOutputFileNumber) {
  char fname[_MAX_STRING_LENGTH_PIC_];
  int ierr;

  //print sodium production rate
  double SourceRate=0.0,TotalSourceRate=0.0;
  int thread;
  double buffer[PIC::nTotalThreads];

  //print the simulation time
  const SpiceInt lenout = 35;
  SpiceChar utcstr[lenout+2];


  et2utc_c(Europa::OrbitalMotion::et,"C",0,lenout,utcstr);

  if (PIC::ThisThread==0) {
    printf("Simulation Time: %s, et=%e\n",utcstr,Europa::OrbitalMotion::et);
    printf("\nSource rates:\n");
  }

  if (Europa::Planet==NULL) return;

  //output total source rate
  double summSourceRate;
  MPI_Reduce(&Europa::TotalInjectionRate,&summSourceRate,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (PIC::ThisThread==0) {
    cout << "Total Injection Rate: " << summSourceRate/PIC::LastSampleLength/PIC::ParticleWeightTimeStep::GlobalTimeStep[_OPLUS_HIGH_SPEC_] << endl;
  }

  Europa::TotalInjectionRate=0.0;

  //save the surface properties
  if (PIC::ThisThread==0) {
    char fname[300];

    sprintf(fname,"pic.SurfaceProperties.%s.dat",utcstr);
//    Europa::Planet->SaveSurfaceDensity(fname,Planet->GetTotalSurfaceElementsNumber());
  }

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  ierr=MPI_Gather(&Europa::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_COMM_WORLD);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Europa sodium source: Photon Stimulated Desorption rate - %e s^{-1}\n",SourceRate);
    TotalSourceRate+=SourceRate;
  }

  Europa::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate=0.0;
#endif

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  ierr=MPI_Gather(&Europa::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_COMM_WORLD);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Europa sodium source: Impact Vaporization rate - %e s^{-1}\n",SourceRate);
    TotalSourceRate+=SourceRate;
  }

  Europa::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate=0.0;
#endif


#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  ierr=MPI_Gather(&Europa::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_COMM_WORLD);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Europa sodium source: Thermal Desorption rate - %e s^{-1}\n",SourceRate);
    TotalSourceRate+=SourceRate;
  }

  Europa::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate=0.0;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  ierr=MPI_Gather(&Europa::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate,1,MPI_DOUBLE,buffer,1, MPI_DOUBLE,0, MPI_COMM_WORLD);
  if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

  for (thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) SourceRate+=buffer[thread];
  if (PIC::LastSampleLength!=0) SourceRate/=PIC::LastSampleLength;

  if (PIC::ThisThread==0) {
    printf("Europa sodium source: Solar Wind Sputtering rate - %e s^{-1}\n",SourceRate);
    TotalSourceRate+=SourceRate;
  }

  Europa::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate=0.0;
#endif

  if (PIC::ThisThread==0) {
    printf("Europa total sodium production rate - %e s^{-1}\n",TotalSourceRate);
    printf("\nNight side return fluxes");
  }


  //print the total night side flux
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    double flux[_EUROPA_SOURCE_MAX_ID_VALUE_+1],FluxAll[PIC::nTotalThreads*(_EUROPA_SOURCE_MAX_ID_VALUE_+1)];
    char ChemSymbol[_MAX_STRING_LENGTH_PIC_];
    double LocalTimeStep;
    int i;

    PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
    exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

    ierr=MPI_Gather(Sampling::PlanetNightSideReturnFlux[spec],_EUROPA_SOURCE_MAX_ID_VALUE_+1,MPI_DOUBLE,FluxAll,_EUROPA_SOURCE_MAX_ID_VALUE_+1, MPI_DOUBLE,0, MPI_COMM_WORLD);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

    for (i=0;i<_EUROPA_SOURCE_MAX_ID_VALUE_+1;i++) {
      for (flux[i]=0.0,thread=0,SourceRate=0.0;thread<PIC::nTotalThreads;thread++) flux[i]+=FluxAll[i+thread*(_EUROPA_SOURCE_MAX_ID_VALUE_+1)];

      if (PIC::LastSampleLength!=0) flux[i]/=PIC::LastSampleLength*LocalTimeStep;
    }


    if (PIC::ThisThread==0) {
       if (spec==0) printf("Night side return flux [s^{-1}]:\n");

       printf("Spec=%i (%s):\n",spec,ChemSymbol);

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
       printf("Particles produced by  Impact Vaporization: %e\n",flux[_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_]);
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
       printf("Particles produced by Photon Stimulated Desorption: %e\n",flux[_EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
       printf("Particles produced by Thermal Desorption: %e\n",flux[_EUROPA_SOURCE__ID__THERMAL_DESORPTION_]);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
       printf("Particles produced by Solar wind sputtering: %e\n",flux[_EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_]);

#endif
    }

    for (i=0;i<_EUROPA_SOURCE_MAX_ID_VALUE_+1;i++) Sampling::PlanetNightSideReturnFlux[spec][i]=0.0;
  }




/*
  //the sodium column density in the tail
  sprintf(fname,"pic.TailColumnDensity.%s,out=%i.dat",utcstr,DataOutputFileNumber);
  ColumnDensityIntegration_Tail(fname);

  //column density distribution map
  sprintf(fname,"pic.ColumnDensityMap.%s,out=%i.dat",utcstr,DataOutputFileNumber);
  ColumnDensityIntegration_Map(fname);
*/



}



/*--------------------------------- Print Output File: BEGIN  --------------------------------------*/
void Europa::Sampling::OutputDataFile::PrintVariableList(FILE* fout,int DataSetNumber) {

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Impact Vaporization Source)\"");
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Photon Stimulated Desorption Source)\"");
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Thermal Desorption Source)\"");
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  fprintf(fout,", \"Sodium number Density(Solar Wind Sputtering Source)\"");
#endif
}


void Europa::Sampling::OutputDataFile::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double TotalMeasure=0.0,Measure=0.0,ImpactVaposizationSource=0.0;
  int i;
  char *SamplingBuffer,*CellNodeSamplingBuffer;

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  double SourcePDS=0.0;
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  double SourceTD=0.0;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  double SourceSW=0.0;
#endif

  CellNodeSamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

    Measure=InterpolationList[i]->Measure;
    if (Measure<=0.0) exit(__LINE__,__FILE__,"Error: non-positive cell volume is found");
    TotalMeasure+=Measure;

    #if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
    ImpactVaposizationSource+=*((double*)(SamplingBuffer+SamplingDensity__ImpactVaporization_Offset));
    #endif

    #if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
    SourcePDS+=*((double*)(SamplingBuffer+SamplingDensity__PhotonStimulatedDesorption_Offset));
    #endif

    #if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
    SourceTD+=*((double*)(SamplingBuffer+SamplingDensity__ThermalDesorption_Offset));
    #endif

    #if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
    SourceSW+=*((double*)(SamplingBuffer+SamplingDensity__SolarWindSputtering_Offset));
    #endif


  }

  if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) ImpactVaposizationSource/=PIC::LastSampleLength*TotalMeasure;
  *((double*)(CellNodeSamplingBuffer+SamplingDensity__ImpactVaporization_Offset))=ImpactVaposizationSource;


  #if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) SourcePDS/=PIC::LastSampleLength*TotalMeasure;
  *((double*)(CellNodeSamplingBuffer+SamplingDensity__PhotonStimulatedDesorption_Offset))=SourcePDS;
  #endif

  #if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) SourceTD/=PIC::LastSampleLength*TotalMeasure;
  *((double*)(CellNodeSamplingBuffer+SamplingDensity__ThermalDesorption_Offset))=SourceTD;
  #endif

  #if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  if ((PIC::LastSampleLength!=0)&&(nInterpolationCoeficients!=0)) SourceSW/=PIC::LastSampleLength*TotalMeasure;
  *((double*)(CellNodeSamplingBuffer+SamplingDensity__SolarWindSputtering_Offset))=SourceSW;
  #endif


}

void Europa::Sampling::OutputDataFile::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double t;
  char *SamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

  //get the number density due to the impact evaporation source
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__ImpactVaporization_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);


#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__PhotonStimulatedDesorption_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__ThermalDesorption_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(SamplingBuffer+SamplingDensity__SolarWindSputtering_Offset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);
#endif




  //other soruce processes......
}
/*--------------------------------- Print Output File: END    --------------------------------------*/



/*--------------------------------- Source Processes: BEGIN  --------------------------------------*/
/*
double Europa::SourceProcesses::totalProductionRate(int spec,void *SphereDataPointer) {
  double res=0.0;

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  if (spec==_O2_SPEC_) {
    res+=Europa::SourceProcesses::ImpactVaporization::GetTotalProductionRate(spec,SphereDataPointer);
  }
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  if (spec==_O2_SPEC_) {
    res+=Europa::SourceProcesses::PhotonStimulatedDesorption::GetTotalProductionRate(spec,SphereDataPointer);
  }
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  if (spec==_O2_SPEC_) {
    res+=Europa::SourceProcesses::ThermalDesorption::GetTotalProductionRate(spec,SphereDataPointer);
  }
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  if (spec==_O2_SPEC_) {
    res+=Europa::SourceProcesses::SolarWindSputtering::GetTotalProductionRate(spec,SphereDataPointer);
  }
#endif

  return res;
}
*/


/*long int Europa::SourceProcesses::InjectionBoundaryModel(void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_GALL_EPHIOD_EUROPA[3],x_IAU_EUROPA[3],v_GALL_EPHIOD_EUROPA[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag;
  int SourceProcessID;*/


  /*

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;

  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_O2_SPEC_];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[_O2_SPEC_];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  ModelParticlesInjectionRate=totalProductionRate(_O2_SPEC_,SphereDataPointer)/ParticleWeight;

  if (ModelParticlesInjectionRate*LocalTimeStep>nMaxInjectedParticles) {
    ParticleWeightCorrection=ModelParticlesInjectionRate*LocalTimeStep/nMaxInjectedParticles;
    ModelParticlesInjectionRate/=ParticleWeightCorrection;
  }




  //calcualte probabilities of each source processes
  double TotalFlux,Flux_ImpactVaporization=0.0,Flux_PSD=0.0,Flux_TD=0.0,Flux_SW_Sputtering=0.0;
  double p,Probability_ImpactVaporization=0.0,Probability_PSD=0.0,Probability_TD=0.0,Probability_SW_Sputtering=0.0;

  TotalFlux=totalProductionRate(_O2_SPEC_,SphereDataPointer);
  Flux_ImpactVaporization=Europa::SourceProcesses::ImpactVaporization::GetTotalProductionRate(_O2_SPEC_,SphereDataPointer);

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  Flux_PSD=Europa::SourceProcesses::PhotonStimulatedDesorption::GetTotalProductionRate(_O2_SPEC_,SphereDataPointer);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  Flux_TD=Europa::SourceProcesses::ThermalDesorption::GetTotalProductionRate(_O2_SPEC_,SphereDataPointer);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  Flux_SW_Sputtering=Europa::SourceProcesses::SolarWindSputtering::GetTotalProductionRate(_O2_SPEC_,SphereDataPointer);
#endif


  Probability_ImpactVaporization=Flux_ImpactVaporization/TotalFlux;
  Probability_PSD=(Flux_PSD+Flux_ImpactVaporization)/TotalFlux;
  Probability_TD=(Flux_PSD+Flux_ImpactVaporization+Flux_TD)/TotalFlux;
  Probability_SW_Sputtering=(Flux_PSD+Flux_ImpactVaporization+Flux_TD+Flux_SW_Sputtering)/TotalFlux;

  //recalcualte the surface injection distributions
  if (Flux_PSD>0.0) PhotonStimulatedDesorption::SurfaceInjectionDistribution.Init();
  if (Flux_TD>0.0) ThermalDesorption::SurfaceInjectionDistribution.Init();
  if (Flux_SW_Sputtering>0.0) SolarWindSputtering::SurfaceInjectionDistribution.Init();

  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {

   //Determine the source processes
   p=rnd();

   if (p<Probability_ImpactVaporization) {
     flag=Europa::SourceProcesses::ImpactVaporization::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode);
     SourceProcessID=_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_;
     if (flag==true) ImpactVaporization::CalculatedTotalSodiumSourceRate+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
   else if (p<Probability_PSD) {
     flag=Europa::SourceProcesses::PhotonStimulatedDesorption::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere);
     SourceProcessID=_EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_;
     if (flag==true) PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
   else if (p<Probability_TD) {
     flag=Europa::SourceProcesses::ThermalDesorption::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere);
     SourceProcessID=_EUROPA_SOURCE__ID__THERMAL_DESORPTION_;
     if (flag==true) ThermalDesorption::CalculatedTotalSodiumSourceRate+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }
   else {
     flag=Europa::SourceProcesses::SolarWindSputtering::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere);
     SourceProcessID=_EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_;
     if (flag==true) SolarWindSputtering::CalculatedTotalSodiumSourceRate+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
   }


   if (flag==false) continue;

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
   if (startNode->block->GetLocalTimeStep(_O2_SPEC_)/LocalTimeStep<rnd()) continue;
#endif

   //generate a particle
   char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

   PIC::ParticleBuffer::SetX(x_GALL_EPHIOD_EUROPA,(PIC::ParticleBuffer::byte*)tempParticleData);
   PIC::ParticleBuffer::SetV(v_GALL_EPHIOD_EUROPA,(PIC::ParticleBuffer::byte*)tempParticleData);
   PIC::ParticleBuffer::SetI(_O2_SPEC_,(PIC::ParticleBuffer::byte*)tempParticleData);

   PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);

   //save the information od the particle origin: the particle origin will be sampled in GALL_EPHIOD coordinate frame
   long int nZenithElement,nAzimuthalElement;
   int el;

   Sphere->GetSurfaceElementProjectionIndex(x_GALL_EPHIOD_EUROPA,nZenithElement,nAzimuthalElement);
   el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

   Europa::Planet->SampleSpeciesSurfaceSourceRate[_O2_SPEC_][el][SourceProcessID]+=ParticleWeight*ParticleWeightCorrection/LocalTimeStep;

   Sampling::SetParticleSourceID(SourceProcessID,(PIC::ParticleBuffer::byte*)tempParticleData);
   Sampling::SetParicleOriginSurfaceElementNumber(el,(PIC::ParticleBuffer::byte*)tempParticleData);

   newParticle=PIC::ParticleBuffer::GetNewParticle();
   newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
   memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

   nInjectedParticles++;

   //for the secondary source processes accout for the decrease of the surface density
   if (SourceProcessID!=_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_) {
     Sphere->GetSurfaceElementProjectionIndex(x_IAU_EUROPA,nZenithElement,nAzimuthalElement);
     el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

     Sphere->SurfaceElementDesorptionFluxUP[_O2_SPEC_][el]+=ParticleWeight*ParticleWeightCorrection;
   }




   //inject the particle into the system
   _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(_O2_SPEC_)*rnd(),startNode,true);
 }

*/


/* return nInjectedParticles;
}*/

/*--------------------------------- SOURCE: Impact Vaporization -----------------------------------*/
/*double Europa::SourceProcesses::ImpactVaporization::GetTotalProductionRate(int spec,void *SphereDataPointer) {
  return SourceRate*pow(HeliocentricDistance/Europa::xEuropaRadial,SourceRatePowerIndex);
}*/





/*=============================== Source Processes: END  ===========================================*/



/*=============================== INTERACTION WITH THE SURFACE: BEGIN  ===========================================*/
int Europa::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x_GALL_EPHIOD_EUROPA,double *v_GALL_EPHIOD_EUROPA,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
  double radiusSphere,*x0Sphere,lNorm[3],rNorm,lVel[3],rVel,c;
  cInternalSphericalData *Sphere;
//  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  int idim;
  double vi,vt,vf,v_LOCAL_GALL_EPHIOD_EUROPA[3],x_LOCAL_GALL_EPHIOD_EUROPA[3],v_LOCAL_IAU_EUROPA[3],x_LOCAL_IAU_EUROPA[3],cosSubsolarAngle,SurfaceTemp,beta;
  SpiceDouble xform[6][6];

#if  _ION_SPUTTERING_MODE_  == _PIC_MODE_ON_
  //do nothing
#elif _ION_SPUTTERING_MODE_  == _PIC_MODE_OFF_
  PIC::ParticleBuffer::DeleteParticle(ptr);
  return _PARTICLE_DELETED_ON_THE_FACE_;
#else
  exit(__LINE__,__FILE__,"the option is not defined");
#endif

  Sphere=(cInternalSphericalData*)SphereDataPointer;
//  startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

  memcpy(v_LOCAL_GALL_EPHIOD_EUROPA,v_GALL_EPHIOD_EUROPA,3*sizeof(double));
  memcpy(x_LOCAL_GALL_EPHIOD_EUROPA,x_GALL_EPHIOD_EUROPA,3*sizeof(double));

  //convert the position vector from GALL_EPHIOD_EUROPA to IAU_EUROPA coordinate frames
  memcpy(xform,OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix,36*sizeof(double));

  x_LOCAL_IAU_EUROPA[0]=xform[0][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[0][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[0][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2];
  x_LOCAL_IAU_EUROPA[1]=xform[1][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[1][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[1][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2];
  x_LOCAL_IAU_EUROPA[2]=xform[2][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[2][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[2][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2];

  //get local surface temperature
  cosSubsolarAngle=Europa::OrbitalMotion::GetCosineSubsolarAngle(x_LOCAL_GALL_EPHIOD_EUROPA);
  SurfaceTemp=Europa::GetSurfaceTemeprature(x_LOCAL_IAU_EUROPA);


  //sample parameters of the back flux: speed is calculate in IAU (relative to the planet) but the flux is sampled in GALL_EPHIOD (one axis is always directed to the Sun)
  //convert the velocity vector from GALL_EPHIOD_EUROPA to IAU_EUROPA coordinate frames
  v_LOCAL_IAU_EUROPA[0]=xform[3][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[3][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[3][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2]+
      xform[3][3]*v_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[3][4]*v_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[3][5]*v_LOCAL_GALL_EPHIOD_EUROPA[2];

  v_LOCAL_IAU_EUROPA[1]=xform[4][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[4][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[4][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2]+
      xform[4][3]*v_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[4][4]*v_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[4][5]*v_LOCAL_GALL_EPHIOD_EUROPA[2];

  v_LOCAL_IAU_EUROPA[2]=xform[5][0]*x_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[5][1]*x_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[5][2]*x_LOCAL_GALL_EPHIOD_EUROPA[2]+
      xform[5][3]*v_LOCAL_GALL_EPHIOD_EUROPA[0]+xform[5][4]*v_LOCAL_GALL_EPHIOD_EUROPA[1]+xform[5][5]*v_LOCAL_GALL_EPHIOD_EUROPA[2];


  vi=sqrt(v_LOCAL_IAU_EUROPA[0]*v_LOCAL_IAU_EUROPA[0]+v_LOCAL_IAU_EUROPA[1]*v_LOCAL_IAU_EUROPA[1]+v_LOCAL_IAU_EUROPA[2]*v_LOCAL_IAU_EUROPA[2]);

  long int nZenithElement,nAzimuthalElement;
  int el;
  double ParticleWeight;

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
#endif


  Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_GALL_EPHIOD_EUROPA,nZenithElement,nAzimuthalElement);
  el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

  Sphere->SampleSpeciesSurfaceReturnFlux[spec][el]+=ParticleWeight;
  Sphere->SampleReturnFluxBulkSpeed[spec][el]+=vi*ParticleWeight;

  //sample returned flux on the night side of the planet
  if (x_LOCAL_GALL_EPHIOD_EUROPA[0]<0.0) {
    //the night size
    int id;

    id=Sampling::GetParticleSourceID(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
    Europa::Sampling::PlanetNightSideReturnFlux[spec][id]+=ParticleWeight;
  }


  //simulate particle/surface interaction
  int ReturnCode=-1;
  double Yield,WeightCorrectionFactor,SputteringSpeed,SputteringVelocity[3],phi,theta,vSputtered_IAU[3];
  double e0[3],e1[3],e2[3],l;
  long int newParticle;

  switch (spec) {
  case _OPLUS_THERMAL_SPEC_: case _OPLUS_HIGH_SPEC_: case _O2PLUS_SPEC_:

    switch (spec) {
    case _OPLUS_THERMAL_SPEC_: case _OPLUS_HIGH_SPEC_:
      Yield=Europa::InjectEuropaMagnetosphericEPDIons::SputteringYield(vi,_MASS_(_O_),1);
      break;
    case _O2PLUS_SPEC_:
      Yield=Europa::InjectEuropaMagnetosphericEPDIons::SputteringYield(vi,_MASS_(_O2_),2);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the specie is not recognized");
    }

    Yield*=ParticleWeight/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_O2_SPEC_];

    //reconstruct the local coordinate system at the point where the ion intersects the surface of Europa
    e0[0]=x_LOCAL_IAU_EUROPA[0],e0[1]=x_LOCAL_IAU_EUROPA[1],e0[2]=x_LOCAL_IAU_EUROPA[2];
    l=sqrt((e0[0]*e0[0])+(e0[1]*e0[1])+(e0[2]*e0[2]));
    e0[0]/=l,e0[1]/=l,e0[2]/=l;

    if (fabs(e0[0])>1.0E-5) e1[0]=e0[1],e1[1]=-e0[0],e1[2]=0.0;
    else e1[0]=0.0,e1[1]=e0[2],e1[2]=-e0[1];

    l=sqrt((e1[0]*e1[0])+(e1[1]*e1[1])+(e1[2]*e1[2]));
    e1[0]/=l,e1[1]/=l,e1[2]/=l;

    e2[0]=(e0[1]*e1[2])-(e0[2]*e1[1]);
    e2[1]=(e0[2]*e1[0])-(e0[0]*e1[2]);
    e2[2]=(e0[0]*e1[1])-(e0[1]*e1[0]);

    //convert velocity vector from IAU_EUROPA -> GALL_EPHIOD_EUROPA coordinate frame
    memcpy(xform,OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix,36*sizeof(double));


    while (Yield>0.0) {
      Europa::EuropaO2Neutrals::O2SputterInjection(SputteringSpeed,WeightCorrectionFactor);
      if (WeightCorrectionFactor>Yield) WeightCorrectionFactor=Yield;
      Yield-=WeightCorrectionFactor;

      WeightCorrectionFactor*=PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_]/PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];

      //distrubute the speed
      //p=cos^2(theta)

      do {
        theta=Pi/2.0*rnd();
      }
      while (pow(cos(theta),2)<rnd());

      phi=2.0*Pi*rnd();

      //Important: e0 is directed along the normal, and e1 and e2 are tangent to the surface !!!!!
      double cosTheta,sinTheta,cosPhi,sinPhi;

      cosTheta=cos(theta);
      sinTheta=sin(theta);

      cosPhi=cos(phi);
      sinPhi=sin(phi);

      vSputtered_IAU[0]=SputteringSpeed*(cosTheta*e0[0]+sinTheta*(cosPhi*e1[0]+sinPhi*e2[0]));
      vSputtered_IAU[1]=SputteringSpeed*(cosTheta*e0[1]+sinTheta*(cosPhi*e1[1]+sinPhi*e2[1]));
      vSputtered_IAU[2]=SputteringSpeed*(cosTheta*e0[2]+sinTheta*(cosPhi*e1[2]+sinPhi*e2[2]));

      //convert velocity of the sputtered particle into the "global" coordinate frame
      v_LOCAL_GALL_EPHIOD_EUROPA[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[3][3]*vSputtered_IAU[0]+xform[3][4]*vSputtered_IAU[1]+xform[3][5]*vSputtered_IAU[2];

      v_LOCAL_GALL_EPHIOD_EUROPA[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[4][3]*vSputtered_IAU[0]+xform[4][4]*vSputtered_IAU[1]+xform[4][5]*vSputtered_IAU[2];

      v_LOCAL_GALL_EPHIOD_EUROPA[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[5][3]*vSputtered_IAU[0]+xform[5][4]*vSputtered_IAU[1]+xform[5][5]*vSputtered_IAU[2];


      //generate new particle and inject it into the system
      newParticle=PIC::ParticleBuffer::GetNewParticle();

      PIC::ParticleBuffer::SetX(x_LOCAL_GALL_EPHIOD_EUROPA,newParticle);
      PIC::ParticleBuffer::SetV(v_LOCAL_GALL_EPHIOD_EUROPA,newParticle);
      PIC::ParticleBuffer::SetI(_O2_SPEC_,newParticle);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,newParticle);

     _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_]*rnd(),(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter);

    }

    PIC::ParticleBuffer::DeleteParticle(ptr);
    ReturnCode=_PARTICLE_DELETED_ON_THE_FACE_;

    break;

  case _O2_SPEC_:
    //redistribute the particle velocity and reject it back into the domain
    //get the external normal
    e0[0]=-x_LOCAL_IAU_EUROPA[0],e0[1]=-x_LOCAL_IAU_EUROPA[1],e0[2]=-x_LOCAL_IAU_EUROPA[2];
    l=sqrt((e0[0]*e0[0])+(e0[1]*e0[1])+(e0[2]*e0[2]));
    e0[0]/=l,e0[1]/=l,e0[2]/=l;

    {
      double EmptyArray[3]={0.0,0.0,0.0};
      double SurfaceTemp=GetSurfaceTemeprature(x_LOCAL_IAU_EUROPA);

      PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_EUROPA,EmptyArray,SurfaceTemp,e0,_O2_SPEC_,-1);
    }

    //tranform the velocity into the "global" coordinate frame
    memcpy(xform,OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix,36*sizeof(double));

    v_LOCAL_GALL_EPHIOD_EUROPA[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[3][3]*v_LOCAL_IAU_EUROPA[0]+xform[3][4]*v_LOCAL_IAU_EUROPA[1]+xform[3][5]*v_LOCAL_IAU_EUROPA[2];

    v_LOCAL_GALL_EPHIOD_EUROPA[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[4][3]*v_LOCAL_IAU_EUROPA[0]+xform[4][4]*v_LOCAL_IAU_EUROPA[1]+xform[4][5]*v_LOCAL_IAU_EUROPA[2];

    v_LOCAL_GALL_EPHIOD_EUROPA[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[5][3]*v_LOCAL_IAU_EUROPA[0]+xform[5][4]*v_LOCAL_IAU_EUROPA[1]+xform[5][5]*v_LOCAL_IAU_EUROPA[2];

    memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
    ReturnCode=_PARTICLE_REJECTED_ON_THE_FACE_;
    break;


  default:
    exit(__LINE__,__FILE__,"Error: the species is not recognized");
  }

  return ReturnCode;

/*

  //check if the particle sticks to the surface of the planet
  if ((spec==_O2_SPEC_)||(spec==_O_PLUS_SPEC_)) if (rnd()<SodiumStickingProbability(SurfaceTemp)) {
    //the particle is abserbed by the surface


    Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_IAU_EUROPA,nZenithElement,nAzimuthalElement);
    el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

    Sphere->SodiumSurfaceAreaDensity_FluxDOWN[el]+=ParticleWeight;

    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }

  // Europa addition
  //check if the particle sticks to the surface of the planet
  if (spec==_O2_SPEC_) {
    //the particle is adsorbed by the surface

	Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_IAU_EUROPA,nZenithElement,nAzimuthalElement);
    el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

    Sphere->O2SurfaceAreaDensity_FluxDOWN[el]+=ParticleWeight;

    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }






  //distribute velocity with the local surface temperature
  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*SurfaceTemp));
  for (vt=0.0,idim=0;idim<3;idim++) vt+=pow(sqrt(-log(rnd()))/beta*cos(PiTimes2*rnd()),2);

  vt=sqrt(vt);
  vf=AccomodationCoefficient*vt+(1.0-AccomodationCoefficient)*vi;

  //determine the normal vector at the intersection point
  Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

  for (rNorm=0.0,rVel=0.0,c=0.0,idim=0;idim<DIM;idim++) {
    lNorm[idim]=x_LOCAL_IAU_EUROPA[idim]-x0Sphere[idim];
    rNorm+=pow(lNorm[idim],2);

    lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
    rVel+=pow(lVel[idim],2);

    c+=lNorm[idim]*lVel[idim];
  }

  rVel=sqrt(rVel);

  if (c>0.0) {
    //the distributed velocity vector is directed into the domain
    for (rVel=vf/rVel,idim=0;idim<3;idim++) v_LOCAL_IAU_EUROPA[idim]=lVel[idim]*rVel;
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
      x_LOCAL_IAU_EUROPA[idim]=vf*(lVel[idim]-2.0*c1*lNorm[idim]);
    }
  }

  //convert velocity vector from IAU_EUROPA -> GALL_EPHIOD_EUROPA coordinate frame
  memcpy(xform,OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix,36*sizeof(double));

  v_LOCAL_GALL_EPHIOD_EUROPA[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
      xform[3][3]*v_LOCAL_IAU_EUROPA[0]+xform[3][4]*v_LOCAL_IAU_EUROPA[1]+xform[3][5]*v_LOCAL_IAU_EUROPA[2];

  v_LOCAL_GALL_EPHIOD_EUROPA[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
      xform[4][3]*v_LOCAL_IAU_EUROPA[0]+xform[4][4]*v_LOCAL_IAU_EUROPA[1]+xform[4][5]*v_LOCAL_IAU_EUROPA[2];

  v_LOCAL_GALL_EPHIOD_EUROPA[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
      xform[5][3]*v_LOCAL_IAU_EUROPA[0]+xform[5][4]*v_LOCAL_IAU_EUROPA[1]+xform[5][5]*v_LOCAL_IAU_EUROPA[2];


  memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
  return _PARTICLE_REJECTED_ON_THE_FACE_;
*/
}

/*=============================== INTERACTION WITH THE SURFACE: END  ===========================================*/



/*=================================== Output surface parameters =================================================*/
void Europa::Sampling::OutputSurfaceDataFile::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Total Flux Down [s^{-1} m^{-2}]\", \"Total Flux Up [s^{-1} m^{-2}]\", \"Surface Content [s^{-1} m^{-2}]\", \"Returned particles' bulk speed [m/s]\"");

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
   fprintf(fout,", \"Impact Vaporization Flux [s^{-1} m^{-2}]\"");
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
   fprintf(fout,", \"Photon Stimulated Desorption Flux [s^{-1} m^{-2}]\"");
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
   fprintf(fout,", \"Thermal Desorption Flux [s^{-1} m^{-2}]\"");
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
   fprintf(fout,", \"Solar Wind Sputtering Flux [s^{-1} m^{-2}]\", \"Soral Wind Flux [s^{-1} m^{-2}]\"");
#endif
}

void Europa::Sampling::OutputSurfaceDataFile::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double t,TotalFluxDown=0.0,TotalFluxUp=0.0,SurfaceContent=0.0,BulkSpeed=0.0;


#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
   double FluxIV=0.0;
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
   double FluxPDS=0.0;
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
   double FluxTD=0.0;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
   double FluxSW=0.0,SolarWindIncidentFlux=0.0;
#endif




  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
    InterpolationCoefficient=Sphere->GetSurfaceElementArea(nSurfaceElement);

    BulkSpeed=Sphere->SampleReturnFluxBulkSpeed[spec][nSurfaceElement];
    TotalFluxDown=Sphere->SampleSpeciesSurfaceReturnFlux[spec][nSurfaceElement];
    SurfaceContent=Sphere->SurfaceElementPopulation[spec][nSurfaceElement]*InterpolationCoefficient;

    if (PIC::LastSampleLength!=0) {
      BulkSpeed/=PIC::LastSampleLength;
      TotalFluxDown/=PIC::LastSampleLength;
      SurfaceContent/=PIC::LastSampleLength;
    }

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_]/PIC::LastSampleLength : 0.0;
    FluxIV+=t;
    TotalFluxUp+=t;
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]/PIC::LastSampleLength : 0.0;
    FluxPDS+=t;
    TotalFluxUp+=t;
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EUROPA_SOURCE__ID__THERMAL_DESORPTION_]/PIC::LastSampleLength : 0.0;
   FluxTD+=t;
   TotalFluxUp+=t;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_]/PIC::LastSampleLength : 0.0;
   FluxSW+=t;
   TotalFluxUp+=t;

   SolarWindIncidentFlux+=Sphere->SolarWindSurfaceFlux[nSurfaceElement]*InterpolationCoefficient;
#endif

    InterpolationNormalization+=InterpolationCoefficient;
  }

  if (ThisThread==0)  {
    //collect sampled data from all processors
    for (int thread=1;thread<nTotalThreads;thread++) {
      TotalFluxDown+=pipe->recv<double>(thread);
//      SurfaceContent+=pipe->recv<double>(thread);  All processors have the same distribution of surface content map
      TotalFluxUp+=pipe->recv<double>(thread);
      BulkSpeed+=pipe->recv<double>(thread);

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
      FluxIV+=pipe->recv<double>(thread);
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
      FluxPDS+=pipe->recv<double>(thread);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
      FluxTD+=pipe->recv<double>(thread);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
      FluxSW+=pipe->recv<double>(thread);
#endif
    }

    if (TotalFluxDown>0.0) BulkSpeed/=TotalFluxDown;

    fprintf(fout," %e %e %e %e ",TotalFluxDown/InterpolationNormalization,TotalFluxUp/InterpolationNormalization,SurfaceContent/InterpolationNormalization,BulkSpeed);

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
    fprintf(fout," %e ",FluxIV/InterpolationNormalization);
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
    fprintf(fout," %e ",FluxPDS/InterpolationNormalization);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
      fprintf(fout," %e ",FluxTD/InterpolationNormalization);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
      fprintf(fout," %e %e ",FluxSW/InterpolationNormalization,SolarWindIncidentFlux/InterpolationNormalization);
#endif

  }
  else {
    pipe->send(TotalFluxDown);
//    pipe->send(SurfaceContent);     All processors have the same distribution of surface content map
    pipe->send(TotalFluxUp);
    pipe->send(BulkSpeed);

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
    pipe->send(FluxIV);
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
    pipe->send(FluxPDS);
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
    pipe->send(FluxTD);
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
    pipe->send(FluxSW);
#endif
  }
}

void Europa::Sampling::OutputSurfaceDataFile::flushCollectingSamplingBuffer(cInternalSphericalData* Sphere) {
  int s,el,i;
  int maxSurfaceSourceID=-1;

#if _EUROPA_SOURCE__IMPACT_VAPORIZATION_ == _EUROPA_SOURCE__ON_
  if (maxSurfaceSourceID<_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_) maxSurfaceSourceID=_EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_;
#endif

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
  if (maxSurfaceSourceID<_EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_) maxSurfaceSourceID=_EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_;
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
  if (maxSurfaceSourceID<_EUROPA_SOURCE__ID__THERMAL_DESORPTION_) maxSurfaceSourceID=_EUROPA_SOURCE__ID__THERMAL_DESORPTION_;
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_ == _EUROPA_SOURCE__ON_
  if (maxSurfaceSourceID<_EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_) maxSurfaceSourceID=_EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_;
#endif


  for (s=0;s<PIC::nTotalSpecies;s++) {
    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      for (i=0;i<maxSurfaceSourceID+1;i++) Sphere->SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;
    }
  }

  for (s=0;s<PIC::nTotalSpecies;s++) {
    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      Sphere->SurfaceElementPopulation[s][el]=0.0;
    }
  }

  for (s=0;s<PIC::nTotalSpecies;s++) {
    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      Sphere->SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
    }
  }

  PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
}


void Europa::Sampling::OutputSurfaceDataFile::PrintTitle(FILE* fout) {
  fprintf(fout,"TITPLE=\"SurfaceData:  TAA=?\"");
}



bool Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  //double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  //static double vSW[3]={4.0E5,000.0,000.0},nSW=5.0E6,tempSW=8.0E4;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      //startNode->GetExternalNormal(ExternalNormal,nface);
      //ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nSW,tempSW,vSW,ExternalNormal,_O_PLUS_SPEC_);

      //if (ModelParticlesInjectionRate>0.0) return true;
      return true;
    }
  }

  return false;
}


double Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double BlockSurfaceArea,ExternalNormal[3],res=0.0;
  bool ExternalFaces[6];
  int nface;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      //high energy ions (O+)
      if (spec==_OPLUS_HIGH_SPEC_) res+=GetTotalProductionRate(spec)*BlockSurfaceArea;

      //thermal ions (O+)
      if (spec==_OPLUS_THERMAL_SPEC_) res+=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Thermal_OPlus_NumberDensity,Thermal_OPlus_Temperature,Thermal_OPlus_BulkVelocity,ExternalNormal,spec)*BlockSurfaceArea;
    }
  }

  return res;
}

long int Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int s;
  long int nInjectedParticles=0;

  for (s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}

//the default sticking probability function
double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {
  ReemissionParticleFraction=1.0;

  return 0.0;
}

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

//calcualte the column integrals
//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int nVariables=0;

  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int cnt=0;

  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}

