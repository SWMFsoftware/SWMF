//$Id$
//Evaluate density and pressure along the spacecraft trajectory using the Liouville theoreme

/*
 * RosinaMeasurements_Liouville.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: vtenishe
 */

#include <iostream>
#include <iomanip>

#include "pic.h"
#include "RosinaMeasurements.h"
#include "Comet.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

extern double *productionDistributionUniformNASTRAN;
extern double HeliocentricDistance,subSolarPointAzimuth,subSolarPointZenith;
extern bool *definedFluxBjorn;
extern double positionSun[3];
extern double **productionDistributionNASTRAN;
extern bool *probabilityFunctionDefinedNASTRAN;

//Parameters of the tests
int SphericalNucleusTest_Mode=SphericalNucleusTestMode_OFF;
static double SphericalNucleusTest_SourceRate=1.0;
static double SphericalNucleusTest_Temperature=100.0;
static double SphericalNucleusTest_Location[3]={1.0E4,0.0,0.0};
static double SphericalNucleusTest_LineOfSightRG[3]={-1.0,0.0,0.0};
static double SphericalNucleusTest_LineOfSightNG[3]={-1.0,0.0,0.0};

//step between simulated observation points
static int RosinaDataSimulationStep=1;

static bool AdjustSurfaceInjectionRate=false;
static double NudeGaugeDensitySinCorrectionFactor=1.0/3.8;

static double CutoffNudeGaugeSolidAngle=0.006;

//registering of particles that moved in the direction of the line of sight by NG
static const int SelfShadowingNG_ModeOff=0;
static const int SelfShadowingNG_ModeSimple=1;

static const int SelfShadowingNGMode=SelfShadowingNG_ModeOff;
static const double SelfShadowingNGCosAngleLimit=0.0;

//account for the species depencent sensitivity of COPS
static const bool SpeciesSensitivityMode=_PIC_MODE_OFF_;

//surface temeprature correction factor
static const double SurfaceTemperatureCorrectionFactor=1.0;

//the model of calculating of the RG pressure
static const int PressureCalculationRG_ModeCalibrationPressure=0;
static const int PressureCalculationRG_ModeFluxBalance=1;
static const int PressureCalculationRGMode=PressureCalculationRG_ModeCalibrationPressure;

//estimate the unsertanties due to the trajectory
int TrajectoryUncertanties_nTest=10;
bool TrajectoryUncertanties_Search=true;
double TrajectoryUncertanties_Radius=50.0;
double TrajectoryUncertanties_AngularLimit=5.0;

static int VelocityInjectionMode=_ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__RANDOMLY_DIRECTED_FLUX_MAXWELLAIN_;

void RosinaSample::Liouville::EvaluateLocation(int spec,double& OriginalSourceRate, double& ModifiedSourceRate,double& NudeGaugePressure,double& NudeGaugeDensity,double& NudeGaugeFlux, double& RamGaugePressure,double& RamGaugeDensity,double& RamGaugeFlux,cRosinaSamplingLocation RosinaLocation,
    double& SurfaceAreaContributedNudeGaugeMeasurements,double& SurfaceFluxContributedNudeGaugeMeasurements) {
  int iTest,iSurfaceElement;
  double c,l[3],*xLocation,rLocation,xIntersection[3],beta;
  double SampledNudeGaugeDensity=0.0,SampledRamGaugeFlux=0.0;
  double localSurfaceAreaContributedNudeGaugeMeasurements=0.0,localSurfaceFluxContributedNudeGaugeMeasurements=0.0;

  const int nTotalTests=_ROSINA_SAMPLE__LIOUVILLE__EVALUATE_LOCATION__TEST_NUMBER_;
  const double ThrehondSourceRateH2O=0.9E18;
  const double ThrehondSourceRateCO2=1.8E18;

  //species dependent COPS "bata-factors"
  const double BetaFactorH2O=0.893;
  const double BetaFactorCO2=0.704;

  //estimate the total source rate
  OriginalSourceRate=0.0,ModifiedSourceRate=0.0;


  //sample particle flux and density at the spacecraft location without accounting for shielding, and instrument orientation (used for debugging)
  bool DisregardInstrumentOrientationFlag=false;

  //change the location and orientation of the spacecraft in the case of the sperical nucleus test
  //DEBUG: BEGIN
  switch (SphericalNucleusTest_Mode) {
  case SphericalNucleusTestMode_Point:
    for (int idim=0;idim<3;idim++) {
      RosinaLocation.x[idim]=SphericalNucleusTest_Location[idim];
      RosinaLocation.RamGauge.LineOfSight[idim]=SphericalNucleusTest_LineOfSightRG[idim];
      RosinaLocation.NudeGauge.LineOfSight[idim]=SphericalNucleusTest_LineOfSightNG[idim];
    }
    break;

  case SphericalNucleusTestMode_Random:
    {
      double e1[3],rTest;
      int idim;

      //generate a random location
      rTest=Vector3D::Length(SphericalNucleusTest_Location);
      Vector3D::Distribution::Uniform(SphericalNucleusTest_Location);

      for (idim=0;idim<3;idim++) {
        SphericalNucleusTest_LineOfSightRG[idim]=-SphericalNucleusTest_Location[idim];
        SphericalNucleusTest_Location[idim]*=rTest;
      }

      Vector3D::GetRandomNormFrame(SphericalNucleusTest_LineOfSightNG,e1,SphericalNucleusTest_LineOfSightRG);

      //copy the location to all processors
      MPI_Bcast(SphericalNucleusTest_LineOfSightRG,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Bcast(SphericalNucleusTest_LineOfSightNG,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Bcast(SphericalNucleusTest_Location,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

      //output the location of the test point
      if (PIC::ThisThread==0) {
        printf("$PREFIX: Test Location=%e, %e, %e\tLine of sight RG=%e %e %e\tLine of sight NG=%e %e %e\n",
            SphericalNucleusTest_Location[0],SphericalNucleusTest_Location[1],SphericalNucleusTest_Location[2],
            SphericalNucleusTest_LineOfSightRG[0],SphericalNucleusTest_LineOfSightRG[1],SphericalNucleusTest_LineOfSightRG[2],
            SphericalNucleusTest_LineOfSightNG[0],SphericalNucleusTest_LineOfSightNG[1],SphericalNucleusTest_LineOfSightNG[2]);
      }
    }

    for (int idim=0;idim<3;idim++) {
      RosinaLocation.x[idim]=SphericalNucleusTest_Location[idim];
      RosinaLocation.RamGauge.LineOfSight[idim]=SphericalNucleusTest_LineOfSightRG[idim];
      RosinaLocation.NudeGauge.LineOfSight[idim]=SphericalNucleusTest_LineOfSightNG[idim];
    }
  }
  //DEBUG: END


  for (iSurfaceElement=0;iSurfaceElement<CutCell::nBoundaryTriangleFaces;iSurfaceElement++) {
    double ThrehondSourceRate,SourceRate=productionDistributionNASTRAN[spec][iSurfaceElement]/CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    //sample the original source rate
    OriginalSourceRate+=SourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    //change the source rate if needed
    switch (spec) {
    case _H2O_SPEC_:
      ThrehondSourceRate=ThrehondSourceRateH2O;
      break;
    case _CO2_SPEC_:
      ThrehondSourceRate=ThrehondSourceRateCO2;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not known");
    }

    if ((SourceRate<ThrehondSourceRate)&&(AdjustSurfaceInjectionRate==true)) {
      switch (spec) {
      case _H2O_SPEC_:
        productionDistributionNASTRAN[spec][iSurfaceElement]=ThrehondSourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
        SourceRate=ThrehondSourceRate;
        break;
      default:
        if (0.98*productionDistributionNASTRAN[_H2O_SPEC_][iSurfaceElement]/CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea<ThrehondSourceRate) {
          //the source rate can be increase only at the bootem of the duck: condition - no intersection with -x
          double x[3],xIntersection[3],lTestX[3]={-1.0,0.0,0.0},lTestY[3]={0.0,1.0,0.0};
          int iIntersectionFaceX,iIntersectionFaceY;

          CutCell::BoundaryTriangleFaces[iSurfaceElement].GetCenterPosition(x);

          iIntersectionFaceX=PIC::RayTracing::FindFistIntersectedFace(x,lTestX,xIntersection,CutCell::BoundaryTriangleFaces+iSurfaceElement);
          iIntersectionFaceY=PIC::RayTracing::FindFistIntersectedFace(x,lTestY,xIntersection,CutCell::BoundaryTriangleFaces+iSurfaceElement);


          if ((iIntersectionFaceX==-1)&&(iIntersectionFaceY==-1) && (x[1]<700.0)) {
            productionDistributionNASTRAN[spec][iSurfaceElement]=ThrehondSourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
            SourceRate=ThrehondSourceRate;
          }

        }
      }
    }


    //sample the modified sources rate
    ModifiedSourceRate+=SourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

    //evaluate the cosine between the external normal of the face and the pointing to the spacecraft
    double x[3],l=0.0,c=0.0;
    CutCell::BoundaryTriangleFaces[iSurfaceElement].GetCenterPosition(x);

    for (int idim=0;idim<3;idim++) {
      c+=CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal[idim]*(RosinaLocation.x[idim]-x[idim]);
      l+=pow(RosinaLocation.x[idim]-x[idim],2);
    }

    CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.ScalarProduct_FaceNormal_SpacecraftLocation=c/sqrt(l);
  }


  //set the angular limit for the ray direction generation
  xLocation=RosinaLocation.x;
  NudeGaugeDensity=0.0,NudeGaugeFlux=0.0;
  RamGaugeFlux=0.0,RamGaugeDensity=0.0;

  //loop through all surface elements
  int iStartSurfaceElement,iFinishSurfaceElement,nSurfaceElementThread;

  nSurfaceElementThread=CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads;
  iStartSurfaceElement=nSurfaceElementThread*PIC::ThisThread;
  iFinishSurfaceElement=iStartSurfaceElement+nSurfaceElementThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iFinishSurfaceElement=CutCell::nBoundaryTriangleFaces;

  double localNudeGaugeDensity=0.0,localNudeGaugeFlux=0.0,localRamGaugeFlux=0.0,localRamGaugeDensity=0.0;

  int idim;
  double x[3],r,ll[3],cosTheta,A,t;
  double tNudeGaugeDensity=0.0,tRamGaugeFlux=0.0,tRamGaugeDensity=0.0;

  if (PIC::ThisThread==0) {
    //administrator
    int iSurfaceElementStep=max(CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads/100,1);
    int thread,SignalTable[PIC::nTotalThreads];
    MPI_Request request[PIC::nTotalThreads];
    MPI_Status status;
    int SendNegativeStartElementNumberCounter=0;

    iStartSurfaceElement=0;
    iFinishSurfaceElement=iSurfaceElementStep;

    //initiate the all recieve
    for (thread=1;thread<PIC::nTotalThreads;thread++) MPI_Irecv(SignalTable+thread-1,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR,request+thread-1);

    //for the nightly tests: use predefined range of the surface elements to be processed
    if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) {
      MPI_Status status[PIC::nTotalThreads];

      iSurfaceElementStep=(CutCell::nBoundaryTriangleFaces/PIC::nTotalThreads-1);

      MPI_Waitall(PIC::nTotalThreads-1,request,status);

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        iStartSurfaceElement=iSurfaceElementStep*(thread-1);
        iFinishSurfaceElement=iStartSurfaceElement+iSurfaceElementStep;

        if (iFinishSurfaceElement>=CutCell::nBoundaryTriangleFaces) iFinishSurfaceElement=CutCell::nBoundaryTriangleFaces;

        MPI_Send(&iStartSurfaceElement,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR);
        MPI_Send(&iFinishSurfaceElement,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR);
      }

      //there is no more surface elements to be processed
      iStartSurfaceElement=-1;
      iFinishSurfaceElement=-1;

      //re-initiate recieve
      for (thread=1;thread<PIC::nTotalThreads;thread++) MPI_Irecv(SignalTable+thread-1,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR,request+thread-1);
    }

    do {
      //waite for any recieved occured
      MPI_Waitany(PIC::nTotalThreads-1,request,&thread,&status);
      thread++;

      //send the limits for the surface elements to check
      MPI_Send(&iStartSurfaceElement,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Send(&iFinishSurfaceElement,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR);

      if (iStartSurfaceElement==-1) {
        SendNegativeStartElementNumberCounter++;
      }
      else {
        //initiate recieve from thread
        MPI_Irecv(SignalTable+thread-1,1,MPI_INT,thread,0,MPI_GLOBAL_COMMUNICATOR,request+thread-1);
      }

      //update the surface element counter
      if (iStartSurfaceElement!=-1) {
        iStartSurfaceElement=iFinishSurfaceElement;
        iFinishSurfaceElement=iStartSurfaceElement+iSurfaceElementStep;

        if (iStartSurfaceElement>=CutCell::nBoundaryTriangleFaces) {
          iStartSurfaceElement=-1,iFinishSurfaceElement=-1;
        }
        else if (iFinishSurfaceElement>CutCell::nBoundaryTriangleFaces) {
          iFinishSurfaceElement=CutCell::nBoundaryTriangleFaces;
        }
      }
    }
    while (SendNegativeStartElementNumberCounter!=PIC::nTotalThreads-1);
  }
  else {
    //send the "ready" signal
    MPI_Request request;
    MPI_Status status;
    int Signal=0;

    do {
      MPI_Isend(&Signal,1,MPI_INT,0,0,MPI_GLOBAL_COMMUNICATOR,&request);
      MPI_Wait(&request,&status);

      //recieve the limits of the surface element number
      MPI_Recv(&iStartSurfaceElement,1,MPI_INT,0,0,MPI_GLOBAL_COMMUNICATOR,&status);
      MPI_Recv(&iFinishSurfaceElement,1,MPI_INT,0,0,MPI_GLOBAL_COMMUNICATOR,&status);

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp parallel for schedule(dynamic,1) default(none) shared(VelocityInjectionMode,NudeGaugeDensitySinCorrectionFactor,DisregardInstrumentOrientationFlag,SphericalNucleusTest,RosinaLocation,iStartSurfaceElement,iFinishSurfaceElement,Rosina,xLocation,CutCell::BoundaryTriangleFaces,productionDistributionNASTRAN,spec,positionSun) \
      private(xIntersection,beta,iSurfaceElement,iTest,idim,c,x,r,cosTheta,A,t,tNudeGaugeDensity,tRamGaugeFlux,tRamGaugeDensity,ll) reduction(+:localNudeGaugeDensity) reduction(+:localNudeGaugeFlux) reduction(+:localRamGaugeFlux) reduction(+:localRamGaugeDensity) \
      reduction(+:localSurfaceAreaContributedNudeGaugeMeasurements) reduction(+:localSurfaceFluxContributedNudeGaugeMeasurements)
    #endif
      for (iSurfaceElement=iStartSurfaceElement;iSurfaceElement<iFinishSurfaceElement;iSurfaceElement++) {
        bool SurfaceElementContributeNudeGaugeMeasurementsFlag=false;

        CutCell::BoundaryTriangleFaces[iSurfaceElement].GetCenterPosition(x);
        c=0.0,tNudeGaugeDensity=0.0,tRamGaugeFlux=0.0,tRamGaugeDensity=0.0;

        for (idim=0;idim<3;idim++) c+=CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal[idim]*(xLocation[idim]-x[idim]);

        if (c>=0.0) {
          //determine the surface temeprature, production rate, and the normalization coefficient
          //parameters of the primary species source at that surface element
          double SourceRate,Temperature,cosSubSolarAngle,x_LOCAL_SO_OBJECT[3]={0.0,0.0,0.0};

          SourceRate=productionDistributionNASTRAN[spec][iSurfaceElement]/CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;

          cosSubSolarAngle=Vector3D::DotProduct(CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal,positionSun)/Vector3D::Length(positionSun);
          if (CutCell::BoundaryTriangleFaces[iSurfaceElement].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow

          Temperature=SurfaceTemperatureCorrectionFactor*Comet::GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);


    //DEBUG: BEGIN
          if (SphericalNucleusTest_Mode!=SphericalNucleusTestMode_OFF) {
            Temperature=SphericalNucleusTest_Temperature;
            SourceRate=SphericalNucleusTest_SourceRate;
          }
    //DEBUG: END

          beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*Temperature));
          A=2.0*SourceRate/Pi*pow(beta,4);

          for (iTest=0;iTest<nTotalTests;iTest++) {
            CutCell::BoundaryTriangleFaces[iSurfaceElement].GetRandomPosition(x);
            for (idim=0;idim<3;idim++) ll[idim]=xLocation[idim]-x[idim];

            if (PIC::RayTracing::FindFistIntersectedFace(x,ll,xIntersection,CutCell::BoundaryTriangleFaces+iSurfaceElement)==-1) {
              //there is the direct access from the point on teh surface to the point of the observation ->  sample the number density and flux
              double ExternalNormal[3];

              switch (VelocityInjectionMode) {
              case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VERTICAL_FLUX_MAXWELLIAN_:
                memcpy(ExternalNormal,CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal,3*sizeof(double));
                break;
              case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__RANDOMLY_DIRECTED_FLUX_MAXWELLAIN_:
              case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VELOCITY_MAXWELLIAN_:
                do {
                  Vector3D::Distribution::Uniform(ExternalNormal);
                }
                while (Vector3D::DotProduct(ExternalNormal,CutCell::BoundaryTriangleFaces[iSurfaceElement].ExternalNormal)<0.0);

                Vector3D::Normalize(ExternalNormal);
                break;
              default:
                exit(__LINE__,__FILE__,"Error: the option is not found");
              }

              r=Vector3D::Length(ll);
              cosTheta=Vector3D::DotProduct(ll,ExternalNormal)/r;

              if (cosTheta<0.0) continue; //the external normal has to be in the direction of the spacecraft

              if ((Vector3D::DotProduct(ll,RosinaLocation.NudeGauge.LineOfSight)<0.0)||(DisregardInstrumentOrientationFlag==true) ||
                  ((SelfShadowingNGMode!=SelfShadowingNG_ModeOff)&&(-Vector3D::DotProduct(ll,RosinaLocation.NudeGauge.LineOfSight)/Vector3D::Length(ll)<SelfShadowingNGCosAngleLimit)) ) {
                //the particle flux can access the nude gauge
                double sinLineOfSightAngle=sqrt(1.0-pow(Vector3D::DotProduct(ll,RosinaLocation.NudeGauge.LineOfSight)/r,2));

                if (DisregardInstrumentOrientationFlag==true) {
                  sinLineOfSightAngle=0.0;
                }

                switch (VelocityInjectionMode) {
                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VERTICAL_FLUX_MAXWELLIAN_:
                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__RANDOMLY_DIRECTED_FLUX_MAXWELLAIN_:
                  tNudeGaugeDensity+=A*cosTheta/(pow(r,2)*pow(beta,3)) *sqrtPi/4.0  *   (1.0+NudeGaugeDensitySinCorrectionFactor*sinLineOfSightAngle);
                  break;

                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VELOCITY_MAXWELLIAN_:
                  A=SourceRate*pow(beta,3)/(Pi*sqrtPi);
                  tNudeGaugeDensity+=1.0/(2.0*pow(r,2)*pow(beta,2))  *   (1.0+NudeGaugeDensitySinCorrectionFactor*sinLineOfSightAngle);
                  break;

                default:
                  exit(__LINE__,__FILE__,"Error: not implemented");
                }

                //sample the nucleus surface area, and the total source rate of the surface that can contribute to the Nude
                if (SurfaceElementContributeNudeGaugeMeasurementsFlag==false) {
                  SurfaceElementContributeNudeGaugeMeasurementsFlag=true;

                  localSurfaceAreaContributedNudeGaugeMeasurements+=CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
                  localSurfaceFluxContributedNudeGaugeMeasurements+=SourceRate*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
                }

              }

              if (((c=Vector3D::DotProduct(ll,RosinaLocation.RamGauge.LineOfSight))<0.0)||(DisregardInstrumentOrientationFlag==true))  {
                 //the particle flux can bedetected by the ram gauge
                double cosLineOfSightAngle=-c/r;

                if (DisregardInstrumentOrientationFlag==true) {
                  cosLineOfSightAngle=1.0;
                }

                switch (VelocityInjectionMode) {
                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VERTICAL_FLUX_MAXWELLIAN_:
                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__RANDOMLY_DIRECTED_FLUX_MAXWELLAIN_:
                  tRamGaugeFlux+=A*cosTheta/(2.0*pow(r,2)*pow(beta,4)) * pow(cosLineOfSightAngle,1); //    *1.0/sqrt(Pi)  ;
                  tRamGaugeDensity+=A*cosTheta/(pow(r,2)*pow(beta,3)) *sqrtPi/4.0;
                  break;

                case _ROSINA_SAMPLE__LIOUVILLE__VECOLITY_DISTRIBUTION_MODE__VELOCITY_MAXWELLIAN_:
                  A=SourceRate*pow(beta,3)/(Pi*sqrtPi);
                  tRamGaugeFlux+=sqrtPi/(4.0*pow(r,2)*pow(beta,3)) * pow(cosLineOfSightAngle,2);
                  break;

                default:
                  exit(__LINE__,__FILE__,"Error: not implemented");
                }


              }
            }
          }

          //sample the contribution of the surface element to the isntrument observation
          //nude gauge density
          t=tNudeGaugeDensity/nTotalTests;
          localNudeGaugeDensity+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
          CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.NudeGaugeDensityContribution[spec]+=t;

          //nude gauge flux
          localNudeGaugeFlux=0.0;

          //ram gauge density
          t=tRamGaugeDensity/nTotalTests;
          localRamGaugeDensity+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
          CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeDensityContribution[spec]+=t;


          //ram gauge flux
          t=tRamGaugeFlux/nTotalTests;
          localRamGaugeFlux+=t*CutCell::BoundaryTriangleFaces[iSurfaceElement].SurfaceArea;
          CutCell::BoundaryTriangleFaces[iSurfaceElement].UserData.RamGaugeFluxContribution[spec]+=t;

        }
      }
    }
    while (iStartSurfaceElement!=-1);
  }

  //combine contributions from all processors
 // double t;

  MPI_Allreduce(&localNudeGaugeDensity,&NudeGaugeDensity,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localNudeGaugeFlux,&NudeGaugeFlux,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  MPI_Allreduce(&localRamGaugeFlux,&RamGaugeFlux,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localRamGaugeDensity,&RamGaugeDensity,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  MPI_Allreduce(&localSurfaceAreaContributedNudeGaugeMeasurements,&SurfaceAreaContributedNudeGaugeMeasurements,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localSurfaceFluxContributedNudeGaugeMeasurements,&SurfaceFluxContributedNudeGaugeMeasurements,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

  //convert the sampled fluxes and density into the nude and ram gauges pressures
  double BetaFactor=1.0;

  const double rgTemperature=293.0;
  const double ngTemperature=293.0;

  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*rgTemperature));

  if (SpeciesSensitivityMode==_PIC_MODE_ON_) switch (spec) {
  case _H2O_SPEC_:
    BetaFactor=BetaFactorH2O;
    break;
  case _CO2_SPEC_:
    BetaFactor=BetaFactorCO2;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the species is unknown");
  }

  switch (PressureCalculationRGMode) {
  case PressureCalculationRG_ModeCalibrationPressure:
    RamGaugePressure=RamGaugeFlux*sqrt(Kbol*rgTemperature*Pi*PIC::MolecularData::GetMass(spec)/2.0)/BetaFactor;
    break;
  case PressureCalculationRG_ModeFluxBalance:
    RamGaugePressure=4.0*Kbol*rgTemperature*RamGaugeFlux/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(spec)))/BetaFactor;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  NudeGaugePressure=NudeGaugeDensity*Kbol*ngTemperature/BetaFactor;

}

//==========================================================================================
//estimate the solud angles that the nucleus is seen by the ram and nude gauges
void RosinaSample::Liouville::GetSolidAngle(double& NudeGaugeNucleusSolidAngle,double& RamGaugeNucleusSolidAngle,int iPoint) {
  NudeGaugeNucleusSolidAngle=0.0;
  RamGaugeNucleusSolidAngle=0.0;

  //solid angle occupied by the nucleus by the nude gauge
  int t,nTotalTests,iTest,iIntersectionFace,NucleusIntersectionCounter=0;
  double l[3],xIntersection[3];

  int iStartTest,iFinishTest,nTestThread;

  nTotalTests=_ROSINA_SAMPLE__LIOUVILLE__GET_SOLID_ANGLE__TEST_NUMBER_;
  iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
  iFinishTest=(nTotalTests/PIC::nTotalThreads)*(PIC::ThisThread+1);

  nTestThread=nTotalTests/PIC::nTotalThreads;
  iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
  iFinishTest=iStartTest+nTestThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iFinishTest=nTotalTests;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic) default(none) shared(CutCell::BoundaryTriangleFaces,iStartTest,iFinishTest,Rosina,iPoint) private(iTest,l,xIntersection,iIntersectionFace) reduction(+:NucleusIntersectionCounter)
#endif
  for (iTest=iStartTest;iTest<iFinishTest;iTest++) {
    //generate a ranfom direction
    do {
      Vector3D::Distribution::Uniform(l);
    }
    while (Vector3D::DotProduct(l,Rosina[iPoint].NudeGauge.LineOfSight)<0.0);

    //determine whether an intersection with the nucleus is found
    iIntersectionFace=PIC::RayTracing::FindFistIntersectedFace(Rosina[iPoint].x,l,xIntersection,NULL);

    if (iIntersectionFace!=-1) {
      NucleusIntersectionCounter++;
      CutCell::BoundaryTriangleFaces[iIntersectionFace].UserData.FieldOfView_NudeGauge=true;
    }
  }

  MPI_Allreduce(&NucleusIntersectionCounter,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  NudeGaugeNucleusSolidAngle=2.0*Pi*double(t)/double(nTotalTests);

  //solid angle occupied by the nucleus by the ram gauge
  NucleusIntersectionCounter=0;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic) default(none) shared(CutCell::BoundaryTriangleFaces,iStartTest,iFinishTest,Rosina,iPoint) private(iTest,l,xIntersection,iIntersectionFace) reduction(+:NucleusIntersectionCounter)
#endif
  for (iTest=iStartTest;iTest<iFinishTest;iTest++) {
    //generate a ranfom direction
    do {
      Vector3D::Distribution::Uniform(l);
    }
    while (Vector3D::DotProduct(l,Rosina[iPoint].RamGauge.LineOfSight)<0.0);

    //determine whether an intersection with the nucleus is found
    iIntersectionFace=PIC::RayTracing::FindFistIntersectedFace(Rosina[iPoint].x,l,xIntersection,NULL);

    if (iIntersectionFace!=-1) {
      NucleusIntersectionCounter++;
      CutCell::BoundaryTriangleFaces[iIntersectionFace].UserData.FieldOfView_RamGauge=true;
    }
  }

  MPI_Allreduce(&NucleusIntersectionCounter,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  RamGaugeNucleusSolidAngle=2.0*Pi*double(t)/double(nTotalTests);
}


//==========================================================================================
//evaluete the fulle set of the measuremetns
void RosinaSample::Liouville::Evaluate() {
#ifndef _NO_SPICE_CALLS_  //the function will be compiled only when SPICE is used (bacause SPICE is needed to adjust the model parameters for each point of the observations
  double NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux;
  double NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle;
  double OriginalSourceRate,ModifiedSourceRate;
  double SurfaceAreaContributedNudeGaugeMeasurements,SurfaceFluxContributedNudeGaugeMeasurements;
  int iPoint,idim,spec;
  SpiceDouble lt,et,xRosetta[3],etStart;
  SpiceDouble       xform[6][6];

  FILE *fGroundTrack=NULL;
  FILE *fAllSpecies;

  const int SurfaceOutputStep=1;

  if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) RosinaDataSimulationStep=400;

  //normalize the line of sight vector of the "spherical nucleus test"
  Vector3D::Normalize(SphericalNucleusTest_LineOfSightRG);
  Vector3D::Normalize(SphericalNucleusTest_LineOfSightNG);

  if (PIC::ThisThread==0) {
    char fname[200];

    sprintf(fname,"%s/Liouville.spec=SUM.dat",PIC::OutputDataFileDirectory);
    fAllSpecies=fopen(fname,"w");

//    fprintf(fAllSpecies,"VARIABLES=\"i\", \"Nude Gauge Pressure\", \"Nude Gauge Density\", \"Nude Gauge Flux\", \
         \"Ram Gauge Pressure\", \"Ram Gauge Density\", \"Ram Gauge Flux\", \"Seconds From The First Point\", \"Nude Guage Nucleus Solid angle\", \"Ram Gauge Nucleus Solid Angle\", \"Altitude\", \"Closest Surface Element Source Rate [m^-2 s^-1]\", \"Nude Gauge COPS Measurements\", \"Ram Gauge COPS Measurements\", \"Original Total Source Rate\", \"Modified Total Source Rate\", \"Ram Gauge Pressure H2O\", \"Ram Gauge Pressure CO2\", \"Nude Gauge Pressure H2O\", \"Nude Gauge Pressure CO2\", \"Total surface area that can contribute to the Nude Gauge measurements\", \"Integrated flux from the surface that could contribute to the nude Gauge measurements\" \n");

    fprintf(fAllSpecies,"VARIABLES= ");
    fprintf(fAllSpecies,"\"Altitude\", ");
    fprintf(fAllSpecies,"\"iPoint\", ");
    fprintf(fAllSpecies,"\"TotalNudeGaugePressure\", \"TotalNudeGaugeDensity\", ");
    fprintf(fAllSpecies,"\"TotalNudeGaugeFlux\", \"TotalRamGaugePressure\", \"TotalRamGaugeDensity\", \"TotalRamGaugeFlux\", ");
    fprintf(fAllSpecies,"\"H2ORamGaugePressure\", \"H2ONudeGaugePressure\", ");
    fprintf(fAllSpecies,"\"CO2RamGaugePressure\", \"CO2NudeGaugePressure\", ");
    fprintf(fAllSpecies,"\"SecondsFromBegining\", ");
    fprintf(fAllSpecies,"\"NudeGaugeNucleusSolidAngle\", \"RamGaugeNucleusSolidAngle\", ");
    fprintf(fAllSpecies,"\"productionDistributionNASTRAN iNucleusClosestFace\", ");
    fprintf(fAllSpecies,"\"100.0*NudeGaugeReferenceData[iPoint]\", \"100.0*RamGaugeReferenceData[iPoint]\", ");
    fprintf(fAllSpecies,"\"125.0*NudeGaugeReferenceData[iPoint]\", \"125.0*RamGaugeReferenceData[iPoint]\", ");
    fprintf(fAllSpecies,"\"75.0*NudeGaugeReferenceData[iPoint]\", \"75.0*RamGaugeReferenceData[iPoint]\", ");
    fprintf(fAllSpecies,"\"OriginalSourceRate\", \"ModifiedSourceRate\", ");
    fprintf(fAllSpecies,"\"SurfaceAreaContributedNudeGaugeMeasurements\", \"SurfaceFluxContributedNudeGaugeMeasurements\"");
    fprintf(fAllSpecies,"\"minTotalRamGaugePressure\" ,\"maxTotalRamGaugePressure\", \"minTotalNudeGaugePressure\", \"maxTotalNudeGaugePressure\"");
    fprintf(fAllSpecies,"\n");


    fGroundTrack=fopen("GroundTracks.dat","w");
    printf("VARIABLES=\"i\", \"Nude Gauge Pressure\", \"Nude Gauge Density\", \"Nude Gauge Flux\",  \"Ram Gauge Pressure\", \"Ram Gauge Density\", \"Ram Gauge Flux\", \"Seconds From The First Point\", \"Nude Guage Nucleus Solid angle\", \"Ram Gauge Nucleus Solid Angle\", \"Altitude\", \"Closest Surface Element Source Rate [m^-2 s^-1]\", \"Nude Gauge COPS Measurements\", \"Ram Gauge COPS Measurements\", \"Original Total Source Rate\", \"Modified Total Source Rate\" \n");
  }

  for (iPoint=0*950;iPoint<RosinaSample::nPoints;iPoint+=RosinaDataSimulationStep) {
    //process the location only if the Nude gauge solid angle is above the cutoff 'CutoffNudeGaugeSolidAngle'
    int ProcessPoinFlag=true;

    //simulate the solid angles and the measuremetns
    GetSolidAngle(NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,iPoint);

    if ((NudeGaugeNucleusSolidAngle<CutoffNudeGaugeSolidAngle)||(Rosina[iPoint].Altitude<0.0)) {
      ProcessPoinFlag=false;
    }

    MPI_Bcast(&ProcessPoinFlag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    if (ProcessPoinFlag==false) continue;

    //cange Step with the altitude
    if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_OFF_) {
      if (Rosina[iPoint].Altitude<500) RosinaDataSimulationStep=12;
      if (Rosina[iPoint].Altitude<100) RosinaDataSimulationStep=6;
      if (Rosina[iPoint].Altitude<50) RosinaDataSimulationStep=3;
      if (Rosina[iPoint].Altitude<20) RosinaDataSimulationStep=1;
    }

    //set the vectors, location of the Sub for the observation
    //init line-of-sight vectors
    //set the location of the spacecraft
    utc2et_c(RosinaSample::ObservationTime[iPoint],&et);
    spkpos_c("ROSETTA",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xRosetta,&lt);
    for (idim=0;idim<3;idim++) xRosetta[idim]*=1.0E3;

    if (iPoint==0) etStart=et;

    Rosina[iPoint].SecondsFromBegining=et-etStart;
    for (idim=0;idim<3;idim++) Rosina[iPoint].x[idim]=xRosetta[idim];

    sxform_c("ROS_SPACECRAFT","67P/C-G_CK",et,xform);

    //set pointing of COPS nude gauge: negative y
    SpiceDouble NudeGauge[6]={0.0,-1.0,0.0 ,0.0,0.0,0.0};
    SpiceDouble l[6];

    mxvg_c(xform,NudeGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[iPoint].NudeGauge.LineOfSight[idim]=l[idim];


    //set pointing of COPS ram gauge: positive z
    SpiceDouble RamGauge[6]={0.0,0.0,1.0 ,0.0,0.0,0.0};

    mxvg_c(xform,RamGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[iPoint].RamGauge.LineOfSight[idim]=l[idim];

    //update the location of the Sun and the nucleus shadowing
    SpiceDouble xSun[3];

    spkpos_c("SUN",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xSun,&lt);
    reclat_c(xSun,&HeliocentricDistance,&subSolarPointAzimuth,&subSolarPointZenith);

    HeliocentricDistance*=1.0E3;

    positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth)*sin(subSolarPointZenith);
    positionSun[2]=HeliocentricDistance*cos(subSolarPointZenith);

    PIC::RayTracing::SetCutCellShadowAttribute(positionSun,true);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      definedFluxBjorn[spec]=false,probabilityFunctionDefinedNASTRAN[spec]=false;

      Comet::GetTotalProductionRateBjornNASTRAN(spec);
    }


    Comet::BjornNASTRAN::Init();


    //reset the fiew of view flags
    for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
      CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_NudeGauge=false;
      CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_RamGauge=false;
    }


    //save surface source rate corresponding to the point of the obervation
    if ((PIC::ThisThread==0)&&((iPoint/RosinaDataSimulationStep)%SurfaceOutputStep==0)) {
       FILE *fSource;
       char fname[1000];
       int iface,spec;

       for (spec=0;spec<PIC::nTotalSpecies;spec++) {
         sprintf(fname,"%s/SoruceRate.spec=%i.iPoint=%i.bin",PIC::OutputDataFileDirectory,spec,iPoint);
         fSource=fopen(fname,"w");

         for (iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
           double t=productionDistributionNASTRAN[spec][iface];

           fwrite(&t,sizeof(double),1,fSource);
         }

         fclose(fSource);
       }
    }


    //correct the source rate
    bool FluxCorrectionFlag=false;

 
    if (_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_OFF_) FluxCorrectionFlag=false; 

    if (FluxCorrectionFlag==true) { //(false) {
      static bool InitFlag=false;
      static double *CorrectionMaskTable=NULL;

      if (InitFlag==false) {
        //init the CorrectionMaskTable
        InitFlag=true;
        CorrectionMaskTable=new double [CutCell::nBoundaryTriangleFaces];


        class cCorectionEntry {
        public:
          double refMax,refMean,refMin;
          double *refData;

          void Set(const char *fname) {
            FILE *fRef=fopen(fname,"r");

            refData=new double [CutCell::nBoundaryTriangleFaces];
            fread(refData,CutCell::nBoundaryTriangleFaces,sizeof(double),fRef);

            refMax=refData[0];
            refMin=refData[0];
            refMean=0.0;

            for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
              if (refMax<refData[iface]) refMax=refData[iface];
              if (refMin>refData[iface]) refMin=refData[iface];
              refMean+=refData[iface];
            }

            refMean/=CutCell::nBoundaryTriangleFaces;

            std::cout << std::endl << "Ref Data File=" << fname << ", refMin=" << refMin << ", refMean=" << refMean << ", refMax=" << refMax << std::endl << std::flush;
            fclose(fRef);
          }

          cCorectionEntry() {
            refMax=-1.0,refMin=-1.0,refMean=-1.0;
            refData=NULL;
          }

          cCorectionEntry(const char* fname) {
            refMax=-1.0,refMin=-1.0,refMean=-1.0;
            refData=NULL;

            Set(fname);
          }
        };


        char fname[1000];
        int iface,spec;
        double RefMaxValue=0.0;
        double RefThrehold=0.75;
        double RefMultiplier=0.5;

        //open the reference files
       cCorectionEntry p1046spec0("NudeGaugeDensityContribution.spec=0.iPoint=1046.bin");
       cCorectionEntry p1058spec0("NudeGaugeDensityContribution.spec=0.iPoint=1058.bin");
       cCorectionEntry p948spec0("NudeGaugeDensityContribution.spec=0.iPoint=948.bin");
       cCorectionEntry p1028spec0("NudeGaugeDensityContribution.spec=0.iPoint=1028.bin");
       cCorectionEntry p1028spec1("NudeGaugeDensityContribution.spec=1.iPoint=1028.bin");
       cCorectionEntry p1043spec0("NudeGaugeDensityContribution.spec=0.iPoint=1043.bin");

       bool ApplyRefCorrection_p1046spec0=false;
       bool ApplyRefCorrection_p1058spec0=false;
       bool ApplyRefCorrection_p948spec0=false;
       bool ApplyRefCorrection_p1028spec0=false;
       bool ApplyRefCorrection_p1028spec1=false;
       bool ApplyRefCorrection_p1043spec0=false;

        //create the correction mask table
        for (iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
          CorrectionMaskTable[iface]=1.0; //0.0001*ref3Mean*ref7Mean;
          //if (Reference7[iface]*Reference5[iface]>ref5Mean*ref7Mean) CorrectionMaskTable[iface]*=Reference7[iface]*Reference3[iface];

          if (ApplyRefCorrection_p1028spec0==true) {
            if ( (p1028spec0.refData[iface]/p1028spec0.refMax>5.0E-2) && (p1043spec0.refData[iface]/p1043spec0.refMax<1.0E-1) )  CorrectionMaskTable[iface]/=2.0E0;
          }

          //if (Reference3[iface]>0.0) CorrectionMaskTable[iface]*=Reference3[iface]/ref3Max;
          //if (Reference7[iface]/ref7Max>8.0E-1) CorrectionMaskTable[iface]*=2.0;
        }
      }

      //apply the correction mask
      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
         productionDistributionNASTRAN[spec][iface]*=CorrectionMaskTable[iface];
        }
      }
    }


    //evaluate the location
    RosinaSample::cRosinaSamplingLocation RosinaLocation=Rosina[iPoint];

    //test the actual location of the spacecraft as well as its vicinity
    double TotalOriginalSourceRate=0.0,TotalModifiedSourceRate=0.0,TotalNudeGaugePressure=0.0;
    double OriginalSourceRateH2O=0.0,ModifiedSourceRateH2O=0.0;
    double OriginalSourceRateCO2=0.0,ModifiedSourceRateCO2=0.0;
    double TotalNudeGaugeDensity=0.0,TotalNudeGaugeFlux=0.0,TotalRamGaugePressure=0.0,TotalRamGaugeDensity=0.0,TotalRamGaugeFlux=0.0;
    double H2ORamGaugePressure=0.0,H2ONudeGaugePressure=0.0;
    double CO2RamGaugePressure=0.0,CO2NudeGaugePressure=0.0;

    double trajectoryTotalNudeGaugePressure=0.0,trajectoryTotalRamGaugePressure=0.0;
    double minTotalRamGaugePressure=-1.0,maxTotalRamGaugePressure=-1.0;
    double minTotalNudeGaugePressure=-1.0,maxTotalNudeGaugePressure=-1.0;

    int iLocationTest,nLocationTestPoints=(TrajectoryUncertanties_Search==false) ? 1 :1+TrajectoryUncertanties_nTest;


    //evaluate the location
    TotalOriginalSourceRate=0.0,TotalModifiedSourceRate=0.0,TotalNudeGaugePressure=0.0;
    TotalNudeGaugeDensity=0.0,TotalNudeGaugeFlux=0.0,TotalRamGaugePressure=0.0,TotalRamGaugeDensity=0.0,TotalRamGaugeFlux=0.0;
    H2ORamGaugePressure=0.0,H2ONudeGaugePressure=0.0;
    CO2RamGaugePressure=0.0,CO2NudeGaugePressure=0.0;

    trajectoryTotalNudeGaugePressure=0.0,trajectoryTotalRamGaugePressure=0.0;


    //loop through all species
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      //save the original source rate of the species
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++)  {
        CutCell::BoundaryTriangleFaces[iface].UserData.OriginalSourceRate[spec]=productionDistributionNASTRAN[spec][iface]/CutCell::BoundaryTriangleFaces[iface].SurfaceArea;
      }
    }

    //start loop over possible spacecraft offsetfs from the trijectory, and variations of the pointing
    for (iLocationTest=0;iLocationTest<nLocationTestPoints;iLocationTest++) {
      //generate the search location
      double e0[3],e1[3],e2[3];

      RosinaLocation=Rosina[iPoint];

      switch (iLocationTest) {
      case 0:
        //do nothing
        break;
      default:
        //generate a location of the shifted from that provided by SPICE by the distance TrajectoryUnsertanties::Radius
        memcpy(e2,Rosina[iPoint].v,3*sizeof(double));

        Vector3D::Normalize(e2);
        Vector3D::GetNormFrame(e0,e1,e2);
        Vector3D::Distribution::Circle::Uniform(RosinaLocation.x,e0,e1,Rosina[iPoint].x,min(0.1*Rosina[iPoint].Altitude,TrajectoryUncertanties_Radius));
        MPI_Bcast(RosinaLocation.x,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

        if (TrajectoryUncertanties_Search==true) {
          //generate a new coordinate frame of the spacecraft, and then recalculate pointing directions of the gauges
          int idim;
          double phi,CosPhi,SinPhi,e0new[3],e1new[3],e2new[3],e0[3]={1.0,0.0,0.0},e1[3]={0.0,1.0,0.0},e2[3]={0.0,0.0,1.0};

          if (PIC::ThisThread==0) {
            //rotation in e0-e1 plane
            phi=TrajectoryUncertanties_AngularLimit/180.0*Pi*(2.0*rnd()-1.0);
            CosPhi=cos(phi);
            SinPhi=sin(phi);

            for (idim=0;idim<3;idim++) {
              e0new[idim]=CosPhi*e0[idim]+SinPhi*e1[idim];
              e1new[idim]=-SinPhi*e0[idim]+CosPhi*e1[idim];
            }

            Vector3D::Copy(e0,e0new);
            Vector3D::Copy(e1,e1new);

            //rotation ins the e0-e2 plane
            phi=TrajectoryUncertanties_AngularLimit/180.0*Pi*(2.0*rnd()-1.0);
            CosPhi=cos(phi);
            SinPhi=sin(phi);

            for (idim=0;idim<3;idim++) {
              e0new[idim]=CosPhi*e0[idim]+SinPhi*e2[idim];
              e2new[idim]=-SinPhi*e0[idim]+CosPhi*e2[idim];
            }

            Vector3D::Copy(e0,e0new);
            Vector3D::Copy(e2,e2new);

            //rotation ins the e1-e2 plane
            phi=TrajectoryUncertanties_AngularLimit/180.0*Pi*(2.0*rnd()-1.0);
            CosPhi=cos(phi);
            SinPhi=sin(phi);

            for (idim=0;idim<3;idim++) {
              e1new[idim]=CosPhi*e1[idim]+SinPhi*e2[idim];
              e2new[idim]=-SinPhi*e1[idim]+CosPhi*e2[idim];
            }

            Vector3D::Copy(e1,e1new);
            Vector3D::Copy(e2,e2new);

            //recalculate pointing directions of the nude and ram gauges
            double NudeGaugeNewPointing[3],RamGaugeNewPointing[3];

            for (int i=0;i<3;i++) {
              NudeGaugeNewPointing[i]=RosinaLocation.NudeGauge.LineOfSight[0]*e0[i]+RosinaLocation.NudeGauge.LineOfSight[1]*e1[i]+RosinaLocation.NudeGauge.LineOfSight[2]*e2[i];
              RamGaugeNewPointing[i]=RosinaLocation.RamGauge.LineOfSight[0]*e0[i]+RosinaLocation.RamGauge.LineOfSight[1]*e1[i]+RosinaLocation.RamGauge.LineOfSight[2]*e2[i];
            }

            memcpy(RosinaLocation.NudeGauge.LineOfSight,NudeGaugeNewPointing,3*sizeof(double));
            memcpy(RosinaLocation.RamGauge.LineOfSight,RamGaugeNewPointing,3*sizeof(double));
          }

          MPI_Bcast(RosinaLocation.NudeGauge.LineOfSight,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
          MPI_Bcast(RosinaLocation.RamGauge.LineOfSight,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
        }
      }

      //reset the summators
      TotalOriginalSourceRate=0.0,TotalModifiedSourceRate=0.0,TotalNudeGaugePressure=0.0;
      TotalNudeGaugeDensity=0.0,TotalNudeGaugeFlux=0.0,TotalRamGaugePressure=0.0,TotalRamGaugeDensity=0.0,TotalRamGaugeFlux=0.0;
      H2ORamGaugePressure=0.0,H2ONudeGaugePressure=0.0;
      CO2RamGaugePressure=0.0,CO2NudeGaugePressure=0.0;


      //perform evaluation of the location
      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        EvaluateLocation(spec,OriginalSourceRate,ModifiedSourceRate,NudeGaugePressure,NudeGaugeDensity,NudeGaugeFlux,RamGaugePressure,RamGaugeDensity,RamGaugeFlux,RosinaLocation,
          SurfaceAreaContributedNudeGaugeMeasurements,SurfaceFluxContributedNudeGaugeMeasurements);

        TotalOriginalSourceRate+=OriginalSourceRate;
        TotalModifiedSourceRate+=ModifiedSourceRate;
        TotalNudeGaugePressure+=NudeGaugePressure;
        TotalNudeGaugeDensity+=NudeGaugeDensity;
        TotalNudeGaugeFlux+=NudeGaugeFlux;
        TotalRamGaugePressure+=RamGaugePressure;
        TotalRamGaugeDensity+=RamGaugeDensity;
        TotalRamGaugeFlux+=RamGaugeFlux;

        switch (spec) {
        case _H2O_SPEC_:
          H2ORamGaugePressure=RamGaugePressure,H2ONudeGaugePressure=NudeGaugePressure;
          if (iLocationTest==0) OriginalSourceRateH2O=OriginalSourceRate,ModifiedSourceRateH2O=ModifiedSourceRate;
          break;
        case _CO2_SPEC_:
          CO2RamGaugePressure=RamGaugePressure,CO2NudeGaugePressure=NudeGaugePressure;
          if (iLocationTest==0) OriginalSourceRateCO2=OriginalSourceRate,ModifiedSourceRateCO2=ModifiedSourceRate;
          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is not found");
        }

        //save the modified source rate of the species
        for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++)  {
          CutCell::BoundaryTriangleFaces[iface].UserData.ModifiedSourceRate[spec]=productionDistributionNASTRAN[spec][iface]/CutCell::BoundaryTriangleFaces[iface].SurfaceArea;
        }
      } //end of the loop over all species for a given location and pointing direction

      switch (iLocationTest) {
      case 0:
        //the location on the spacecraft trajectory
        //output solution for this location
        if (PIC::ThisThread==0) {
          maxTotalRamGaugePressure=TotalRamGaugePressure,minTotalRamGaugePressure=TotalRamGaugePressure;
          maxTotalNudeGaugePressure=TotalNudeGaugePressure,minTotalNudeGaugePressure=TotalNudeGaugePressure;

          trajectoryTotalNudeGaugePressure=TotalNudeGaugePressure,trajectoryTotalRamGaugePressure=TotalRamGaugePressure;

          fprintf(fGroundTrack,"%e %e %e\n",Rosina[iPoint].xNucleusClosestPoint[0],Rosina[iPoint].xNucleusClosestPoint[1],Rosina[iPoint].xNucleusClosestPoint[2]);


          fprintf(fAllSpecies,"%e   %i    %e %e   %e %e %e %e   %e %e   %e %e   %e  %e %e   %e  %e %e   %e %e   %e %e  %e %e  %e %e  ",
              Rosina[iPoint].Altitude,
              iPoint,
              TotalNudeGaugePressure,TotalNudeGaugeDensity,
              TotalNudeGaugeFlux,TotalRamGaugePressure,TotalRamGaugeDensity,TotalRamGaugeFlux,
              H2ORamGaugePressure,H2ONudeGaugePressure,
              CO2RamGaugePressure,CO2NudeGaugePressure,
              Rosina[iPoint].SecondsFromBegining,
              NudeGaugeNucleusSolidAngle,RamGaugeNucleusSolidAngle,
              ((Rosina[iPoint].Altitude>0.0) ? (productionDistributionNASTRAN[_H2O_SPEC_][Rosina[iPoint].iNucleusClosestFace]+productionDistributionNASTRAN[_CO2_SPEC_][Rosina[iPoint].iNucleusClosestFace])/CutCell::BoundaryTriangleFaces[Rosina[iPoint].iNucleusClosestFace].SurfaceArea : 0.0),
              100.0*RosinaSample::NudeGaugeReferenceData[iPoint],100.0*RosinaSample::RamGaugeReferenceData[iPoint],
              125.0*RosinaSample::NudeGaugeReferenceData[iPoint],125.0*RosinaSample::RamGaugeReferenceData[iPoint],
              75.0*RosinaSample::NudeGaugeReferenceData[iPoint],75.0*RosinaSample::RamGaugeReferenceData[iPoint],
              OriginalSourceRate,ModifiedSourceRate,
              SurfaceAreaContributedNudeGaugeMeasurements,SurfaceFluxContributedNudeGaugeMeasurements);

          fflush(fGroundTrack);
          fflush(fAllSpecies);
        }

        break;
      default:
        //sample max and min values
        if ((maxTotalRamGaugePressure==0.0)||((TotalRamGaugePressure>0.0)&&(maxTotalRamGaugePressure<TotalRamGaugePressure))) maxTotalRamGaugePressure=TotalRamGaugePressure;
        if ((minTotalRamGaugePressure==0.0)||((TotalRamGaugePressure>0.0)&&(minTotalRamGaugePressure>TotalRamGaugePressure))) minTotalRamGaugePressure=TotalRamGaugePressure;

        if ((maxTotalNudeGaugePressure==0.0)||((TotalNudeGaugePressure>0.0)&&(maxTotalNudeGaugePressure<TotalNudeGaugePressure))) maxTotalNudeGaugePressure=TotalNudeGaugePressure;
        if ((minTotalNudeGaugePressure==0.0)||((TotalNudeGaugePressure>0.0)&&(minTotalNudeGaugePressure>TotalNudeGaugePressure))) minTotalNudeGaugePressure=TotalNudeGaugePressure;
      }


      if (PIC::ThisThread==0) {
        double CosNudeGaugeAngle,CosRamGaugeAngle,NudeGaugeAngle,RamGaugeAngle;

        CosNudeGaugeAngle=Vector3D::DotProduct(Rosina[iPoint].NudeGauge.LineOfSight,RosinaLocation.NudeGauge.LineOfSight);
        NudeGaugeAngle=acos((CosNudeGaugeAngle<1.0)? CosNudeGaugeAngle : 1.0-1.0E-15)/Pi*180.0;

        CosRamGaugeAngle=Vector3D::DotProduct(Rosina[iPoint].RamGauge.LineOfSight,RosinaLocation.RamGauge.LineOfSight);
        RamGaugeAngle=acos((CosRamGaugeAngle<1.0)? CosRamGaugeAngle : 1.0-1.0E-15)/Pi*180.0;



        std::cout << "s/c location/orientation perturbation: iLocationTest=" <<  iLocationTest <<
            ",\t TotalRamGaugePressure=" << minTotalRamGaugePressure << " < " << trajectoryTotalRamGaugePressure << " <" << maxTotalRamGaugePressure <<
            ", TotalNudeGaugePressure=" << minTotalNudeGaugePressure << " < " << trajectoryTotalNudeGaugePressure << " <" << maxTotalNudeGaugePressure <<
            ", d=" << std::setprecision(5) << sqrt(pow(Rosina[iPoint].x[0]-RosinaLocation.x[0],2)+pow(Rosina[iPoint].x[1]-RosinaLocation.x[1],2)+pow(Rosina[iPoint].x[2]-RosinaLocation.x[2],2)) <<
            ",\t angles=" << NudeGaugeAngle << ",  " << RamGaugeAngle <<
            ", cos(angles)=" << CosNudeGaugeAngle << ",  " << CosRamGaugeAngle <<
            " (" << __LINE__ << "@" << __FILE__ << ")" << std::endl << std::flush;
      }


    }


    //output pressure max and min values
    if (PIC::ThisThread==0) {
      fprintf(fAllSpecies,"%e %e    %e %e\n",minTotalRamGaugePressure,maxTotalRamGaugePressure,minTotalNudeGaugePressure,maxTotalNudeGaugePressure);
      fflush(fAllSpecies);

      printf("%e |1  %i |2  %e %e |3  %e %e |4  %e %e |5  %e %e |6  %e %e |7 %e %e\n",
/*1*/     Rosina[iPoint].Altitude,
/*2*/     iPoint,
/*3*/     trajectoryTotalNudeGaugePressure,trajectoryTotalRamGaugePressure,
/*4*/     minTotalRamGaugePressure,maxTotalRamGaugePressure,
/*5*/     minTotalNudeGaugePressure,maxTotalNudeGaugePressure,
/*6*/     OriginalSourceRateH2O,ModifiedSourceRateH2O,
/*7*/     OriginalSourceRateCO2,ModifiedSourceRateCO2,
/*8*/     SurfaceAreaContributedNudeGaugeMeasurements,SurfaceFluxContributedNudeGaugeMeasurements);
   }


    //save the surface properties
    if ((iPoint/RosinaDataSimulationStep)%SurfaceOutputStep==0) {
      //create the surface output file, and clean the buffers
      char fname[200];

//      sprintf(fname,"%s/SurfaceContributionParameters.iPoint=%i.dat",PIC::OutputDataFileDirectory,iPoint);
//      CutCell::PrintSurfaceData(fname);

      //save binary file of the surface contribution
      FILE *fSource;
      int iface,spec;


      double *t=new double [CutCell::nBoundaryTriangleFaces];
      double *tall=new double [CutCell::nBoundaryTriangleFaces];

      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        for (iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++)  t[iface]=CutCell::BoundaryTriangleFaces[iface].UserData.NudeGaugeDensityContribution[spec];

        MPI_Reduce(t,tall,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

        if (PIC::ThisThread==0) {
          double minNudeGaugeDensityContribution=-1.0,maxNudeGaugeDensityContribution=-1.0;

          for (iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
            if ((minNudeGaugeDensityContribution<0.0)||(minNudeGaugeDensityContribution>t[iface])) minNudeGaugeDensityContribution=t[iface];
            if ((maxNudeGaugeDensityContribution<0.0)||(maxNudeGaugeDensityContribution<tall[iface])) maxNudeGaugeDensityContribution=tall[iface];
          }

          std::cout << "minNudeGaugeDensityContribution(spec=" << spec <<")=" << minNudeGaugeDensityContribution << std::endl << std::flush;
          std::cout << "maxNudeGaugeDensityContribution(spec=" << spec <<")=" << maxNudeGaugeDensityContribution << std::endl << std::flush;

          sprintf(fname,"%s/NudeGaugeDensityContribution.spec=%i.iPoint=%i.bin",PIC::OutputDataFileDirectory,spec,iPoint);
          fSource=fopen(fname,"w");
          fwrite(tall,sizeof(double),CutCell::nBoundaryTriangleFaces,fSource);
          fclose(fSource);
        }
      }

      delete [] t;
      delete [] tall;

      //determine the "field of view map" for the given configuration and output the surface data file
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
        bool FieldOfViewRG=false,FieldOfViewNG=false;
        double x[3],ll[3],xIntersection[3];

        CutCell::BoundaryTriangleFaces[iface].GetRandomPosition(x);
        for (idim=0;idim<3;idim++) ll[idim]=Rosina[iPoint].x[idim]-x[idim];

        if (PIC::RayTracing::FindFistIntersectedFace(x,ll,xIntersection,CutCell::BoundaryTriangleFaces+iface)==-1) {
          //the location can be seen from the spacecraft.
          //veryfy whether the location can be seen by the instruments
          if (Vector3D::DotProduct(Rosina[iPoint].RamGauge.LineOfSight,ll)<0.0) FieldOfViewRG=true;
          if (Vector3D::DotProduct(Rosina[iPoint].NudeGauge.LineOfSight,ll)<0.0) FieldOfViewNG=true;
        }

        for (int i=0;i<3;i++) {
          if (FieldOfViewRG==true) CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_RamGauge=true;
          if (FieldOfViewNG==true) CutCell::BoundaryTriangleFaces[iface].UserData.FieldOfView_NudeGauge=true;
        }
      }

      sprintf(fname,"%s/SurfaceContributionParameters.iPoint=%i.dat",PIC::OutputDataFileDirectory,iPoint);
      CutCell::PrintSurfaceData(fname);

      //clear the samplign buffers
      for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        CutCell::BoundaryTriangleFaces[iface].UserData.NudeGaugeDensityContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.NudeGaugeFluxContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.RamGaugeDensityContribution[spec]=0.0;
        CutCell::BoundaryTriangleFaces[iface].UserData.RamGaugeFluxContribution[spec]=0.0;
      }
    }
  }

  if (PIC::ThisThread==0) {
    fclose(fGroundTrack);
    fclose(fAllSpecies);
  }
#endif //_NO_SPICE_CALLS_
}



