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
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="GALL_EOJ";  //Europa-centric inertial Orbital Jupiter:
//This system has its X axis pointing in direction of Europa's motion (assumed to magnetospheric flow past Europa),
//Y is along the direction from Europa to Jupiter orthogonal to X;
//defined in galileo.tf


char Europa::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";

double Europa::TotalInjectionRate=0.0;

//double Europa::swE_Typical[3]={0.0,0.0,0.0};
double Europa::xEuropa[3]={0.0,0.0,0.0},Europa::vEuropa[3]={0.0,0.0,0.0};
double Europa::xEarth[3]={0.0,0.0,0.0},Europa::vEarth[3]={0.0,0.0,0.0};
double Europa::vEuropaRadial=0.0,Europa::xEuropaRadial=0.0;


// dust parameters
double DustSizeMin=1.0e-7;
double DustSizeMax=1.0e-2;
double DustTotalMassProductionRate=1.0e23*_H2O__MASS_;
int    DustSampleIntervals=10;
double DustSizeDistribution=0.0;

//sample the total escape and the toral return fluxes
double Europa::Sampling::TotalEscapeFlux[PIC::nTotalSpecies];
double Europa::Sampling::TotalReturnFlux[PIC::nTotalSpecies];

//sample velocity of the sputtered O2;
double Europa::Sampling::O2InjectionSpeed::SamplingBuffer[Europa::Sampling::O2InjectionSpeed::nSampleIntervals];

//the total number of source processes
//int Europa::nTotalSourceProcesses=0;

//the sphere that represents the planet
//cInternalSphericalData *Europa::Planet=NULL;

/*//the total source rate values for specific source processes
double Europa::SourceProcesses::PhotonStimulatedDesorption::SourceRate=0.0,Europa::SourceProcesses::PhotonStimulatedDesorption::maxLocalSourceRate=0.0;
double Europa::SourceProcesses::ThermalDesorption::SourceRate=0.0,Europa::SourceProcesses::ThermalDesorption::maxLocalSourceRate=0.0;
double Europa::SourceProcesses::SolarWindSputtering::SourceRate=0.0,Europa::SourceProcesses::SolarWindSputtering::maxLocalSourceRate=0.0;

//evaluate numerically the source rate
double Europa::SourceProcesses::PhotonStimulatedDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::ImpactVaporization::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::ThermalDesorption::CalculatedTotalSodiumSourceRate=0.0;
double Europa::SourceProcesses::SolarWindSputtering::CalculatedTotalSodiumSourceRate=0.0;*/


//the pair of fucntions that are used to sample model data, and output of the return amd escape fluxes
void Europa::Sampling::SamplingProcessor() {}
void Europa::Sampling::PrintOutputFile(int nfile) {
//  if (nfile!=0) return;

  if (PIC::ThisThread==0) {
    printf("$PREFIX: Escape and return fluxes\n");
    printf("$PREFIX: spec, Escape flux, Return flux\n");
  }

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    double EscapeRate=0.0,ReturnRate=0.0;
    int ierr;

    ierr=MPI_Reduce(TotalEscapeFlux+spec,&EscapeRate,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

    ierr=MPI_Reduce(TotalReturnFlux+spec,&ReturnRate,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
    if (ierr!=MPI_SUCCESS) exit(__LINE__,__FILE__);

    if (PIC::ThisThread==0) printf("$PREFIX: %s\t %e\t%e\n",PIC::MolecularData::GetChemSymbol(spec),EscapeRate/PIC::LastSampleLength,ReturnRate/PIC::LastSampleLength);

    TotalEscapeFlux[spec]=0.0;
    TotalReturnFlux[spec]=0.0;
  }
}

//process praticles when they cross boundary of the computational domain
int Europa::ParticleDomainBoundaryIntersection(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  //sample the total particle escape rate
  int spec;
  double ParticleWeight;

  spec=PIC::ParticleBuffer::GetI(ptr);
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
  Europa::Sampling::TotalEscapeFlux[spec]+=ParticleWeight;
}


//energy distribution of particles injected via ion sputtering
double Exosphere::SourceProcesses::SolarWindSputtering::EnergyDistributionFunction(double e,int *spec) {
  static const double Ee=0.015*eV2J;


  return 1.0/pow(Ee+e,2); //Brown et al., 1984;  Rubin, 2012 PATM proposal
}


//surface temperature
double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x_LOCAL_SO_OBJECT) {


  return 100.0;
}

double Europa::EnergeticIonSputteringRate(int spec) {
  double res;

  switch (spec) {
  case _O2_SPEC_ :
    res=1.0E26;
    break;
  case _H2O_SPEC_:
    res=2.0E27;
    break;
  default:
    res=0.0;
  }

  return res;
}

//init the model
void Europa::Init_BeforeParser() {
  Exosphere::Init_BeforeParser();

  //set the initial values in the flux sampling buffers
  for (int spec=0;spec<PIC::nTotalSpecies;spec++) Sampling::TotalEscapeFlux[spec]=0.0,Sampling::TotalReturnFlux[spec]=0.0;

  //set up the processoe of particles that leave the computational domain
  PIC::Mover::ProcessOutsideDomainParticles=ParticleDomainBoundaryIntersection;

  //register the sampling functions
  PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(Sampling::SamplingProcessor,Sampling::PrintOutputFile);

  //check the state of the Sputtering source
  if (_EUROPA__SPUTTERING_ION_SOURCE_ == _EUROPA__SPUTTERING_ION_SOURCE__AMPS_KINETIC_IONS_) {
    if (_EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_) {
      exit(__LINE__,__FILE__,"Error: _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ must be _EXOSPHERE_SOURCE__OFF_ when _EUROPA__SPUTTERING_ION_SOURCE_==_EUROPA__SPUTTERING_ION_SOURCE__AMPS_KINETIC_IONS_");
    }

    PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate=EnergeticIonSputteringRate;
  }

  //Get the initial parameters of Europa orbit
  SpiceDouble state[6];
  int idim;

  utc2et_c(SimulationStartTimeString,&OrbitalMotion::et);

  //get initial parameters of Europa's orbit
  spkezr_c("Europa",OrbitalMotion::et,Exosphere::SO_FRAME,"none","Jupiter",state,&OrbitalMotion::lt);

  for (idim=0,xEuropaRadial=0.0;idim<3;idim++) {
    xEuropa[idim]=state[idim]*1.0E3,vEuropa[idim]=state[idim+3]*1.0E3;
    xEuropaRadial+=pow(xEuropa[idim],2);
  }

  xEuropaRadial=sqrt(xEuropaRadial);
  vEuropaRadial=(xEuropaRadial>1.0E-5) ? (xEuropa[0]*vEuropa[0]+xEuropa[1]*vEuropa[1]+xEuropa[2]*vEuropa[2])/xEuropaRadial : 0.0;

  // get initial position of GALILEO for line-of sight
  spkezr_c("GALILEO ORBITER",OrbitalMotion::et,Exosphere::SO_FRAME,"none","Europa",state,&OrbitalMotion::lt);
  for (idim=0;idim<3;idim++) xEarth[idim]=state[idim]*1.0E3,vEarth[idim]=state[idim+3]*1.0E3;

  //init the location of the plume
  OrbitalMotion::UpdateTransformationMartix();
  Plume::SetPlumeLocation();

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  //init the dust model
  ElectricallyChargedDust::Init_BeforeParser();
#endif

}


//ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
void Europa::SWMFdataPreProcessor(double *x,PIC::CPLR::DATAFILE::ICES::cDataNodeSWMF& data) {
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

  spkezr_c("GALILEO ORBITER",Europa::OrbitalMotion::et,Exosphere::SO_FRAME,"none","Europa",State,&lt);

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

  spkezr_c("GALILEO ORBITER",Europa::OrbitalMotion::et,Exosphere::SO_FRAME,"none","Europa",State,&lt);

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
   Europa::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_,tempParticleData);

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
int Europa::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x_SO,double *v_SO,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
  double radiusSphere,*x0Sphere,lNorm[3],rNorm,lVel[3],rVel,c;
  cInternalSphericalData *Sphere;
//  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  int idim;
  double vi,vt,vf,v_LOCAL_SO[3],x_LOCAL_SO[3],v_LOCAL_IAU_EUROPA[3],x_LOCAL_IAU_EUROPA[3],SurfaceTemp,beta;
  SpiceDouble xform[6][6];



  {

     {
       long int nZenithElement,nAzimuthalElement,el;
       double ParticleWeight;
       //the particle is abserbed by the surface

       Sphere=(cInternalSphericalData*)SphereDataPointer;
       ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

      Sphere->GetSurfaceElementProjectionIndex(x_SO,nZenithElement,nAzimuthalElement);
      el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

      Sphere->SampleSpeciesSurfaceReturnFlux[spec][el]+=ParticleWeight;
      Sphere->SampleReturnFluxBulkSpeed[spec][el]+=sqrt(pow(v_SO[0],2)+pow(v_SO[1],2)+pow(v_SO[2],2))*ParticleWeight;

      //sample the total return flux
      Europa::Sampling::TotalReturnFlux[spec]+=ParticleWeight;

    }



/*  double e0[3],l,e1[3],e2[3];
  int i;

      e0[0]=x_SO[0],e0[1]=x_SO[1],e0[2]=x_SO[2];
      l=sqrt((e0[0]*e0[0])+(e0[1]*e0[1])+(e0[2]*e0[2]));
      e0[0]/=l,e0[1]/=l,e0[2]/=l;

      if (fabs(e0[0])>1.0E-5) {
        e1[0]=e0[1], e1[1]=-e0[0],e1[2]=0.0;
      }
      else {
        e1[0]=0.0,e1[1]=e0[2],e1[2]=-e0[1];
      }

      l=sqrt((e1[0]*e1[0])+(e1[1]*e1[1])+(e1[2]*e1[2]));
      e1[0]/=l,e1[1]/=l,e1[2]/=l;

      e2[0]=+(e0[1]*e1[2]-e1[1]*e0[2]);
      e2[1]=-(e0[0]*e1[2]-e1[0]*e0[2]);
      e2[2]=+(e0[0]*e1[1]-e1[0]*e0[1]);

    double c,c1,beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*400.0));
    double vnew[3];

    vnew[0]=sqrt(-log(rnd()))/beta;  //radial velocity

    c=sqrt(-log(rnd()))/beta;
    c1=rnd();

    vnew[1]=c*sin(2.0*Pi*c1);  //tangential velocity
    vnew[2]=c*cos(2.0*Pi*c1);  //tangential velocity

    for (i=0;i<3;i++) {
      v_SO[i]= vnew[0]*e0[i]+vnew[1]*e1[i]+vnew[2]*e2[i];
    }


      return _PARTICLE_REJECTED_ON_THE_FACE_;*/
  }


  Sphere=(cInternalSphericalData*)SphereDataPointer;
//  startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

  memcpy(v_LOCAL_SO,v_SO,3*sizeof(double));
  memcpy(x_LOCAL_SO,x_SO,3*sizeof(double));

  //convert the position vector from GALL_EPHIOD_EUROPA to IAU_EUROPA coordinate frames
  memcpy(xform,OrbitalMotion::SO_to_IAU_TransformationMartix,36*sizeof(double));

  x_LOCAL_IAU_EUROPA[0]=xform[0][0]*x_LOCAL_SO[0]+xform[0][1]*x_LOCAL_SO[1]+xform[0][2]*x_LOCAL_SO[2];
  x_LOCAL_IAU_EUROPA[1]=xform[1][0]*x_LOCAL_SO[0]+xform[1][1]*x_LOCAL_SO[1]+xform[1][2]*x_LOCAL_SO[2];
  x_LOCAL_IAU_EUROPA[2]=xform[2][0]*x_LOCAL_SO[0]+xform[2][1]*x_LOCAL_SO[1]+xform[2][2]*x_LOCAL_SO[2];

  //get local surface temperature
//  cosSubsolarAngle=Europa::OrbitalMotion::GetCosineSubsolarAngle(x_LOCAL_GALL_EPHIOD_EUROPA);
  SurfaceTemp=Europa::GetSurfaceTemeprature(x_LOCAL_IAU_EUROPA);


  //sample parameters of the back flux: speed is calculate in IAU (relative to the planet) but the flux is sampled in GALL_EPHIOD (one axis is always directed to the Sun)
  //convert the velocity vector from GALL_EPHIOD_EUROPA to IAU_EUROPA coordinate frames
  v_LOCAL_IAU_EUROPA[0]=xform[3][0]*x_LOCAL_SO[0]+xform[3][1]*x_LOCAL_SO[1]+xform[3][2]*x_LOCAL_SO[2]+
      xform[3][3]*v_LOCAL_SO[0]+xform[3][4]*v_LOCAL_SO[1]+xform[3][5]*v_LOCAL_SO[2];

  v_LOCAL_IAU_EUROPA[1]=xform[4][0]*x_LOCAL_SO[0]+xform[4][1]*x_LOCAL_SO[1]+xform[4][2]*x_LOCAL_SO[2]+
      xform[4][3]*v_LOCAL_SO[0]+xform[4][4]*v_LOCAL_SO[1]+xform[4][5]*v_LOCAL_SO[2];

  v_LOCAL_IAU_EUROPA[2]=xform[5][0]*x_LOCAL_SO[0]+xform[5][1]*x_LOCAL_SO[1]+xform[5][2]*x_LOCAL_SO[2]+
      xform[5][3]*v_LOCAL_SO[0]+xform[5][4]*v_LOCAL_SO[1]+xform[5][5]*v_LOCAL_SO[2];


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

  //sample the total return flux
  Exosphere::Sampling::TotalPlanetReturnFlux[spec]+=ParticleWeight;

  Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_SO,nZenithElement,nAzimuthalElement);
  el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

  Sphere->SampleSpeciesSurfaceReturnFlux[spec][el]+=ParticleWeight;
  Sphere->SampleReturnFluxBulkSpeed[spec][el]+=vi*ParticleWeight;

  //sample returned flux on the night side of the planet (only for the regular gas particles)
  if (PIC::MolecularData::GetSpecieType(spec)==_PIC_SPECIE_TYPE__GAS_) {
    if (x_LOCAL_SO[0]<0.0) {
      //the night size
      int id;

      id=Sampling::GetParticleSourceID(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
      Europa::Sampling::PlanetNightSideReturnFlux[spec][id]+=ParticleWeight;
    }
  }


  //simulate particle/surface interaction
  int ReturnCode=-1;
  double Yield,WeightCorrectionFactor,SputteringSpeed,SputteringVelocity[3],phi,theta,vSputtered_IAU[3];
  double e0[3],e1[3],e2[3],l;
  long int newParticle;

  switch (spec) {

  case _OH_PLUS_SPEC_: case _H2_PLUS_SPEC_: case _H2O_PLUS_SPEC_ :

    ReturnCode=_PARTICLE_DELETED_ON_THE_FACE_;
    break;

  case _O_PLUS_SPEC_: case _O_PLUS_THERMAL_SPEC_: case _O_PLUS_HIGH_SPEC_: case _O2_PLUS_SPEC_: case _H_PLUS_SPEC_:

#if _EUROPA__SPUTTERING_ION_SOURCE_ == _EUROPA__SPUTTERING_ION_SOURCE__AMPS_KINETIC_IONS_
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__OFF_

    int iTargetSpecies, nTargetSpecies;
    const int *TargetSpeciesTable;
    double *YieldTable;
    
    switch (spec) {
    case _O_PLUS_THERMAL_SPEC_: case _O_PLUS_HIGH_SPEC_: case _O_PLUS_SPEC_:
      //Yield=Europa::InjectEuropaMagnetosphericEPDIons::SputteringYield(vi,_MASS_(_O_),1);
      nTargetSpecies=Sputtering::GetSputteringYield(vi, _O_PLUS_THERMAL_SPEC_, TargetSpeciesTable, YieldTable);
      break;
    case _O2_PLUS_SPEC_:
      //Yield=Europa::InjectEuropaMagnetosphericEPDIons::SputteringYield(vi,_MASS_(_O2_),2);
      nTargetSpecies=Sputtering::GetSputteringYield(vi, _O2_PLUS_SPEC_, TargetSpeciesTable, YieldTable);
      break;
    case _H_PLUS_SPEC_:
      nTargetSpecies=Sputtering::GetSputteringYield(vi, _H_PLUS_SPEC_, TargetSpeciesTable, YieldTable);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the specie is not recognized");
    }


#if  _ION_SPUTTERING_MODE_  == _PIC_MODE_ON_
    Yield*=ParticleWeight/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_O2_SPEC_]*
        PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_]/PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _ION_SPUTTERING_MODE_  == _PIC_MODE_OFF_
    Yield=-1.0;
#else
    exit(__LINE__,__FILE__,"the option is not defined");
#endif

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
    memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

    for(iTargetSpecies = 0; iTargetSpecies < nTargetSpecies; iTargetSpecies++){
      while (YieldTable[iTargetSpecies]>0.0) {
	do {
	  //Europa::EuropaO2Neutrals::O2SputterInjection(SputteringSpeed,WeightCorrectionFactor);
	  Sputtering::GetSputteringSpeed(TargetSpeciesTable[iTargetSpecies], SputteringSpeed, WeightCorrectionFactor);	  
      }
      while (maxSputteredParticleVelocity<SputteringSpeed);

      if (WeightCorrectionFactor>YieldTable[iTargetSpecies]) WeightCorrectionFactor=YieldTable[iTargetSpecies];
      YieldTable[iTargetSpecies]-=WeightCorrectionFactor;

      //sample the injection speed of O2
      int vInterval=(int)(SputteringSpeed/Europa::Sampling::O2InjectionSpeed::dv);
      if (vInterval<Europa::Sampling::O2InjectionSpeed::nSampleIntervals) Europa::Sampling::O2InjectionSpeed::SamplingBuffer[vInterval]+=WeightCorrectionFactor;

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
      v_LOCAL_SO[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[3][3]*vSputtered_IAU[0]+xform[3][4]*vSputtered_IAU[1]+xform[3][5]*vSputtered_IAU[2];

      v_LOCAL_SO[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[4][3]*vSputtered_IAU[0]+xform[4][4]*vSputtered_IAU[1]+xform[4][5]*vSputtered_IAU[2];

      v_LOCAL_SO[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
          xform[5][3]*vSputtered_IAU[0]+xform[5][4]*vSputtered_IAU[1]+xform[5][5]*vSputtered_IAU[2];

      vi=sqrt(vSputtered_IAU[0]*vSputtered_IAU[0]+vSputtered_IAU[1]*vSputtered_IAU[1]+vSputtered_IAU[2]*vSputtered_IAU[2]);


      //generate new particle and inject it into the system
      newParticle=PIC::ParticleBuffer::GetNewParticle();

      PIC::ParticleBuffer::SetX(x_LOCAL_SO,newParticle);
      PIC::ParticleBuffer::SetV(v_LOCAL_SO,newParticle);
      PIC::ParticleBuffer::SetI(TargetSpeciesTable[iTargetSpecies],newParticle);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrectionFactor,newParticle);
      Europa::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_,PIC::ParticleBuffer::GetParticleDataPointer(newParticle));

      //sample the particle injection rate
#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
      ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[TargetSpeciesTable[iTargetSpecies]]*WeightCorrectionFactor;
#else
      exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

      Sphere->SampleSpeciesSurfaceInjectionFlux[TargetSpeciesTable[iTargetSpecies]][el]+=ParticleWeight;
      Sphere->SampleInjectedFluxBulkSpeed[TargetSpeciesTable[iTargetSpecies]][el]+=vi*ParticleWeight;

      //inject the particle into the system
     _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,PIC::ParticleWeightTimeStep::GlobalTimeStep[TargetSpeciesTable[iTargetSpecies]]*rnd(),(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter);

      }
    }
#endif
#endif

//    PIC::ParticleBuffer::DeleteParticle(ptr);
    ReturnCode=_PARTICLE_DELETED_ON_THE_FACE_;
    break;

  case _O_SPEC_: case _H_SPEC_: case _OH_SPEC_: case _H2O_SPEC_ : 
    //the assumptions of the boundary interaction: O,H, OH, H2O <- stick to the surface; O2 and H2 do not stick to the surface
    ReturnCode=_PARTICLE_DELETED_ON_THE_FACE_;
    break;

  case _H2_SPEC_: case _O2_SPEC_:
    //re-eject the particle into the exosphere
    //redistribute the particle velocity and reject it back into the domain
    //get the external normal
    e0[0]=-x_LOCAL_IAU_EUROPA[0],e0[1]=-x_LOCAL_IAU_EUROPA[1],e0[2]=-x_LOCAL_IAU_EUROPA[2];
    l=sqrt((e0[0]*e0[0])+(e0[1]*e0[1])+(e0[2]*e0[2]));
    e0[0]/=l,e0[1]/=l,e0[2]/=l;

    //reflect the particle velocity
    {
    double t=0.0;

    for (int i=0;i<3;i++) t+=e0[i]*v_LOCAL_IAU_EUROPA[i];
    for (int i=0;i<3;i++) v_LOCAL_IAU_EUROPA[i]-=2.0*t*e0[i];




    //REjeberate particle's velocity according to Maxwellian with the local temeprature
    double vbulk[3]={0.0,0.0,0.0};


    PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_EUROPA,vbulk,90.0,e0,spec);
    }



    //tranform the velocity into the "global" coordinate frame
    memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

    v_LOCAL_SO[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[3][3]*v_LOCAL_IAU_EUROPA[0]+xform[3][4]*v_LOCAL_IAU_EUROPA[1]+xform[3][5]*v_LOCAL_IAU_EUROPA[2];

    v_LOCAL_SO[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[4][3]*v_LOCAL_IAU_EUROPA[0]+xform[4][4]*v_LOCAL_IAU_EUROPA[1]+xform[4][5]*v_LOCAL_IAU_EUROPA[2];

    v_LOCAL_SO[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
        xform[5][3]*v_LOCAL_IAU_EUROPA[0]+xform[5][4]*v_LOCAL_IAU_EUROPA[1]+xform[5][5]*v_LOCAL_IAU_EUROPA[2];

    memcpy(v_SO,v_LOCAL_SO,3*sizeof(double));
    ReturnCode=_PARTICLE_REJECTED_ON_THE_FACE_;
    break;
  default:
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if(_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups){
      ReturnCode=_PARTICLE_DELETED_ON_THE_FACE_;
      break;
    }
#endif//_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_   
    exit(__LINE__,__FILE__,"Error: the boundary interaction model is not specified for this species");
  }


  //if a particle is deleted than sample the sticking rate
  if (ReturnCode==_PARTICLE_DELETED_ON_THE_FACE_) {
    el=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

    Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el]+=ParticleWeight;
    Exosphere::Sampling::PlanetSurfaceStickingRate[spec]+=ParticleWeight;
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
/*
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
*/



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
      if (spec==_O_PLUS_HIGH_SPEC_) res+=GetTotalProductionRate(spec)*BlockSurfaceArea;

      //thermal ions (O+)
      if (spec==_O_PLUS_THERMAL_SPEC_) res+=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Thermal_OPlus_NumberDensity,Thermal_OPlus_Temperature,Thermal_OPlus_BulkVelocity,ExternalNormal,spec)*BlockSurfaceArea;
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
  double res=1.0;
  ReemissionParticleFraction=0.0;

  if (spec==_O2_SPEC_) res=0.0;

  return res;
}

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

//calcualte the column integrals
//calculate the sodium column density and plot
int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;

  //column density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (vlist!=NULL) sprintf(vlist,"%s,  \"Column Integral(%s)\"",vlist,PIC::MolecularData::GetChemSymbol(spec));
    nVariables+=1;
  }

  //brightness of the exosphere
  if (vlist!=NULL) sprintf(vlist,"%s, \"Brightness OI(130.4)\", \"Brightness OI(135.6)\"",vlist);
  nVariables+=2;

  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
//do nothing
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity;


  nd=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node);
  for (i=0;i<resLength;i++) res[i]=0.0;

  //integrate density
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    //get the local density number
    NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(spec);
    res[cnt++]=NumberDensity;
//    res[cnt++]=NumberDensity*node->block->GetCenterNode(nd)->GetMeanParticleSpeed(spec);
  }

  //integrate brightness
  if (_O_SPEC_>=0) {
    //40.0E6 <- the characteristic value of the electron density (Hall-98-AJ); The emission rate constants are from Hall-98-AJ
    double O_NumberDensity;

    O_NumberDensity=node->block->GetCenterNode(nd)->GetNumberDensity(_O_SPEC_);
    res[cnt++]=40.0*O_NumberDensity*1.0E-6*(5.5E-9+4.9E-10)/(4.0*Pi);
    res[cnt++]=40.0*O_NumberDensity*1.0E-6*(6.0E-10+1.1E-9)/(4.0*Pi);  //*1.0E-6 is to transfer Oxygen number density to cm^{-3}
  }
  else {
    res[cnt++]=0.0;
    res[cnt++]=0.0;
  }

  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: the length of the vector is not coinsistent with the number of integrated variables");
}

double Exosphere::SourceProcesses::GetInjectionEnergy(int spec, int SourceProcessID) {
	double res,r;
	bool flag;

	//the maximum velocity of the injected particle
	static const double vmax=10.0E3;

	//parameters of the O2 sputtering distribution
  static const double U_o2=0.015*eV2J;  //the community accepted parameter of the energy distribution      //Burger 2010-SSR
  static const double Emax_o2=_O2__MASS_*vmax*vmax/2.0; //the maximum energy of the ejected particle

	//parameters of the H2O sputtring energy distribution
	static const double U_h2o=0.055*eV2J; //the parameter of the distribution (Martin's PATM proposal)     //Burger 2010-SSR
	static const double Emax_h2o=_H2O__MASS_*vmax*vmax/2.0; //the maximum energy of the ejected particle


	switch (spec) {
	case _O2_SPEC_ :
	  switch (SourceProcessID) {
	  case _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_ :
//	    r=rnd();
//	    res=U_o2*Emax_o2*(1.0-r)/(U_o2+r*Emax_o2);

      do {
        r=rnd();
        res=r*U_o2/(1.0-r);    //Burger 2010-SSR
      }
      while (res>Emax_o2);


	    break;
	  default:
	    exit(__LINE__,__FILE__,"Error: not implemented");
	  }

	  break;
	case _H2O_SPEC_:
	  switch (SourceProcessID) {
	  case _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_ :

	    do {
	      r=rnd();
	      res=U_h2o*(r+sqrt(r))/(1.0-r);    //Burger 2010-SSR
	    }
	    while (res>Emax_h2o);

	    break;
	  default:
	    exit(__LINE__,__FILE__,"Error: not implemented");
	  }

	  break;
	default:
	  exit(__LINE__,__FILE__,"Error: not implemented");
	}

/*
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
    if (isfinite(res)==false) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif
*/

	return res;
}

void Europa::Sampling::O2InjectionSpeed::flush() {
  for (int i=0;i<nSampleIntervals;i++) SamplingBuffer[i]=0.0;
}

void Europa::Sampling::O2InjectionSpeed::OutputSampledModelData(int DataOutputFileNumber) {
  double norm;
  int i;

  //collect the distribution function from all processors
  double recvBuffer[nSampleIntervals];

  MPI_Reduce(SamplingBuffer,recvBuffer,nSampleIntervals,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  flush();

  if (PIC::ThisThread==0) {
    //normalize the dsitribution function
    for (norm=0.0,i=0;i<nSampleIntervals;i++) norm+=recvBuffer[i];

    if (norm>0.0) {
      for (i=0;i<nSampleIntervals;i++) recvBuffer[i]/=norm;
    }

    //output the distribution funcion input a file
    char fname[_MAX_STRING_LENGTH_PIC_];
    FILE *fout;

    sprintf(fname,"%s/pic.O2.InjectedSpeedDistribution.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fout=fopen(fname,"w");

    fprintf(fout,"VARIABLES=\"v\", \"f(v)\"\n");
    for (i=0;i<nSampleIntervals;i++) fprintf(fout,"%e %e\n",(0.5+i)*dv,recvBuffer[i]);

    fclose(fout);
  }
}



double Europa::LossProcesses::PhotolyticReactionRate=0.0;
double Europa::LossProcesses::ElectronImpactRate=0.0;
double Europa::LossProcesses::ElectronTemeprature=0.0;

double Europa::LossProcesses::ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double BackgroundPlasmaNumberDensity;
//  double PlasmaBulkVelocity[3],ElectronDensity;

  PhotolyticReactionRate=0.0,ElectronImpactRate=0.0;


//DEBUG -> no chemistry at all
  if ((spec!=_H2O_SPEC_) && (spec!=_O2_SPEC_) && (spec!=_H2_SPEC_) && (spec!=_H_SPEC_) && (spec!=_OH_SPEC_) && (spec!=_O_SPEC_)) {
   PhotolyticReactionAllowedFlag=false;
   return -1.0;
  }


  PIC::CPLR::InitInterpolationStencil(x,node);
  BackgroundPlasmaNumberDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity();

  PhotolyticReactionAllowedFlag=true;

  static const double EuropaHeliocentricDistance=5.2*_AU_;

  static const double PhotolyticReactionRate_H2O=PhotolyticReactions::H2O::GetTotalReactionRate(EuropaHeliocentricDistance);
  static const double PhotolyticReactionRate_O2=PhotolyticReactions::O2::GetTotalReactionRate(EuropaHeliocentricDistance);
  static const double PhotolyticReactionRate_H2=PhotolyticReactions::H2::GetTotalReactionRate(EuropaHeliocentricDistance);
  static const double PhotolyticReactionRate_H=PhotolyticReactions::H::GetTotalReactionRate(EuropaHeliocentricDistance);
  static const double PhotolyticReactionRate_OH=PhotolyticReactions::OH::GetTotalReactionRate(EuropaHeliocentricDistance);
  static const double PhotolyticReactionRate_O=PhotolyticReactions::O::GetTotalReactionRate(EuropaHeliocentricDistance);


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  //determine whether the particle is in a shadow of Europa;  Jupiter will be added later
  bool shadow=false;

  //Europa:
  if (xSun_SO[0]*x[0]+xSun_SO[1]*x[1]+xSun_SO[2]*x[2]<0.0) {
    //the point is behind Euripa and could in a shadow
    //check whether the point is within the cone determined by the Sun and Europa's terminator line
    double cosThetaSquared=0.0,cosThetaSquaredMax=0.0;
    double lengthParticleSquared=0.0,lengthEuropaSquared=0.0;

    for (int i=0;i<3;i++) {
      double t=x[i]-xSun_SO[i];
      double t2=t*t;
      double xSun2=xSun_SO[i]*xSun_SO[i];

      cosThetaSquared+=t2*xSun2;
      lengthParticleSquared+=t2;
      lengthEuropaSquared+=xSun2;
    }

    cosThetaSquared/=lengthParticleSquared*lengthEuropaSquared;
    cosThetaSquaredMax=lengthEuropaSquared/(_RADIUS_(_EUROPA_)*_RADIUS_(_EUROPA_)+lengthEuropaSquared);

    if (cosThetaSquared>cosThetaSquaredMax) shadow=true;
  }

  if (shadow==false) {
    switch (spec) {
    case _H2O_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_H2O;
      break;
    case _O2_SPEC_:
      PhotolyticReactionRate=PhotolyticReactionRate_O2;
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
    default:
      exit(__LINE__,__FILE__,"Error: unknown specie");
    }
  }
#else
  switch (spec) {
  case _H2O_SPEC_:
    PhotolyticReactionRate=PhotolyticReactionRate_H2O;
    break;
  case _O2_SPEC_:
    PhotolyticReactionRate=PhotolyticReactionRate_O2;
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
  default:
    exit(__LINE__,__FILE__,"Error: unknown specie");
  }

#endif

  //calcualte the rate due to the electron impact
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

  if (PhotolyticReactionRate+ElectronImpactRate<=0.0) {
    PhotolyticReactionAllowedFlag=false;
    return -1.0;
  }

  return 1.0/((PhotolyticReactionRate+ElectronImpactRate)*NumericalLossRateIncrease);  //use the "false" reaction event to increase the number of the dauter model particles. Account for this artificial correction in the ExospherePhotoionizationReactionProcessor
}


int Europa::LossProcesses::ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int *ReactionProductsList,nReactionProducts;
  double *ReactionProductVelocity;
  int ReactionChannel;
  bool PhotolyticReactionRoute;


  //init the reaction tables
  static bool initflag=false;
  static double TotalProductYeld_PhotolyticReaction[PIC::nTotalSpecies*PIC::nTotalSpecies];
  static double TotalProductYeld_ElectronImpact[PIC::nTotalSpecies*PIC::nTotalSpecies];

  double HotElectronFraction=0.05;
  static const double ThermalElectronTemeprature=20.0;
  static const double HotElectronTemeprature=250.0;

  if (initflag==false) {
    int iParent,iProduct;

    initflag=true;

    for (iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
      TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=0.0;
      TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=0.0;

      if (PhotolyticReactions::ModelAvailable(iParent)==true) {
        TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
      }

      if (ElectronImpact::ModelAvailable(iParent)==true) {
        TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=
            Europa::ElectronModel::HotElectronFraction*ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,Europa::ElectronModel::HotElectronTemeprature) +
            Europa::ElectronModel::ThermalElectronFraction*ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,Europa::ElectronModel::ThermalElectronTemeprature);
      }
    }
  }

  //determine the type of the reaction
  PhotolyticReactionRoute=(rnd()<PhotolyticReactionRate/(PhotolyticReactionRate+ElectronImpactRate)) ? true : false;

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
             if (rnd()<Europa::ElectronModel::HotElectronFraction) ElectronImpact::GenerateReactionProducts(spec,Europa::ElectronModel::HotElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
             else ElectronImpact::GenerateReactionProducts(spec,Europa::ElectronModel::ThermalElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }

           //check whether the products contain species with spec=specProduct
           for (iProduct=0;iProduct<nReactionProducts;iProduct++) if (ReactionProductsList[iProduct]==specProduct) {
             flag=true;
             break;
           }
         }
         while (flag==false);


         //determine the velocity of the product specie
         double x[3],v[3],c=rnd();

         for (int idim=0;idim<3;idim++) {
           x[idim]=xInit[idim]+c*(xFinal[idim]-xInit[idim]);
           v[idim]=vFinal[idim]+ReactionProductVelocity[idim+3*iProduct];
         }

         //generate a particle
         PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetI(specProduct,(PIC::ParticleBuffer::byte*)tempParticleData);

         #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
         #endif

         //apply condition of tracking the particle
         #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
         PIC::ParticleTracker::InitParticleID(tempParticleData);
         PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,specProduct,tempParticleData,(void*)node);
         #endif


         //get and injection into the system the new model particle
         newParticle=PIC::ParticleBuffer::GetNewParticle();
         newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
         memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

         node=PIC::Mesh::mesh.findTreeNode(x,node);
         _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,rnd()*ProductTimeStep,node);
       }
     }

  }


  return (rnd()<1.0/NumericalLossRateIncrease) ? _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_ : _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
}






















