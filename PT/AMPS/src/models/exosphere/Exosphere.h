/*
 * Mercury.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#ifndef _EXOSPHERE_
#define _EXOSPHERE_

#include "Na.h"

#include "SpiceUsr.h"
#include "SingleVariableDistribution.h"
#include "SingleVariableDiscreteDistribution.h"
#include "constants.h"

/*//path to the SPICE Kernels directory
const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/SPICE/Kernels/MESSENGER/kernels";

//SPICE Kernels to be loaded
const int nFurnishedSPICEkernels=5+12;
const char SPICE_Kernels[nFurnishedSPICEkernels][_MAX_STRING_LENGTH_PIC_]={"spk/msgr_de405_de423s.bsp","fk/msgr_dyn_v600.tf","../../NAIF/naif0010.tls","pck/pck00009_MSGR_v10.tpc","fk/msgr_v210.tf",
    "ik/msgr_epps_v100.ti","ck/msgr20110413.bc","ck/msgr20110414.bc","ck/msgr20110415.bc","ck/msgr20110416.bc","ck/msgr20110417.bc","ck/msgr20110418.bc","ck/msgr20110419.bc","ck/msgr20110420.bc","ck/msgr20110421.bc",
    "sclk/messenger_1486.tsc","spk/msgr_20040803_20140823_od266sc_0.bsp"};



//time stams for sampling remote column density observations
const int nReferenceGroundBasedObservations=5; //25;
const char ReferenceGroundBasedObservationTime[nReferenceGroundBasedObservations][_MAX_STRING_LENGTH_PIC_]={
    "2008-05-18T00:00:00","2008-07-06T00:00:00","2008-11-07T00:00:00","2007-11-12T00:00:00","2007-06-03T00:00:00"};


/ *    "2001-05-25T00:00:00","2006-01-12T00:00:00","2005-12-18T00:00:00","2006-06-17T00:00:00","2003-10-04T00:00:00",
    "2006-10-21T00:00:00","2003-05-07T00:00:00","1999-04-27T00:00:00","1998-05-28T00:00:00","1990-12-10T00:00:00",
    "1990-12-04T00:00:00","1988-01-11T00:00:00","2000-06-05T00:00:00","2003-02-06T00:00:00","2002-12-16T00:00:00",
    "2002-08-21T00:00:00","2003-12-13T00:00:00","2006-06-12T00:00:00","2008-07-13T00:00:00","2006-11-09T00:00:00"};* /*/


//default setting of the exospheric model
//descriptors of the source processes
#define _EXOSPHERE_SOURCE__ON_    0
#define _EXOSPHERE_SOURCE__OFF_   1

#define _EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_              0
#define _EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_     1
#define _EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_               2
#define _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_            3

//define wich of the source processes are active
#define _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_              _EXOSPHERE_SOURCE__ON_
#define _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_     _EXOSPHERE_SOURCE__ON_
#define _EXOSPHERE_SOURCE__THERMAL_DESORPTION_               _EXOSPHERE_SOURCE__ON_
#define _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_            _EXOSPHERE_SOURCE__ON_


//integration mode: steady state/time dependent
#define _EXOSPHERE_INTEGRATION_MODE__STEADY_STATE_    0
#define _EXOSPHERE_INTEGRATION_MODE__TIME_DEPENDENT_  1


#define _EXOSPHERE_INTEGRATION_MODE_ _EXOSPHERE_INTEGRATION_MODE__TIME_DEPENDENT_

//redistribute the surface density of the exospheric component
#define _EXOSPHERE__SURFACE_CONTENT__BALANCE_FLUXES_         0
#define _EXOSPHERE__SURFACE_CONTENT__UNIFORM_                1
#define _EXOSPHERE__SURFACE_CONTENT__RADIAL_DISTRIBUTION_    2
#define _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_           3

//default macro that is called when '_EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_' is used
#define _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(el)   (0.0)

//#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__UNIFORM_
#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__BALANCE_FLUXES_

//user defined settings of the exospheric model
#include "UserDefinition.Exosphere.h"



//set up "dependent" model parameters
//maximum value of the source ID number
#define _EXOSPHERE_SOURCE_MAX_ID_VALUE_ 0

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
#undef _EXOSPHERE_SOURCE_MAX_ID_VALUE_
#define _EXOSPHERE_SOURCE_MAX_ID_VALUE_ _EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
#undef _EXOSPHERE_SOURCE_MAX_ID_VALUE_
#define _EXOSPHERE_SOURCE_MAX_ID_VALUE_ _EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_  == _EXOSPHERE_SOURCE__ON_
#undef _EXOSPHERE_SOURCE_MAX_ID_VALUE_
#define _EXOSPHERE_SOURCE_MAX_ID_VALUE_ _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_
#endif

namespace Exosphere {

  //the name os the simulated object and related coordinate frames
  extern char ObjectName[_MAX_STRING_LENGTH_PIC_],IAU_FRAME[_MAX_STRING_LENGTH_PIC_],SO_FRAME[_MAX_STRING_LENGTH_PIC_];


  //simulation date and position of Object at the time of simulations
  //const char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_]="2011-04-13T00:00:00"; //"2001-03-01T00:00:00";  ////"2011-01-01T00:00:00";
  extern char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_];
  extern double xObject_HCI[3],vObject_HCI[3],xEarth_HCI[3],vEarth_HCI[3],xEarth_SO[3],vEarth_SO[3],xSun_SO[3],vSun_SO[3];
  extern double vObjectRadial,xObjectRadial;

  extern double vObject_SO_FROZEN[3]; //Velocity of Object in inertial frame, which axis coinsides with that of SO at a given time
  extern double RotationVector_SO_FROZEN[3],RotationRate_SO_FROZEN; //the vectors of the direction of rotation and the rotation rate of the SO in SO_FROZEN


  //typical solar wind conditions far from the planet
  extern const double swVelocity_Typical[3];
  extern const double swB_Typical[3];
  extern const double swTemperature_Typical,swNumberDensity_Typical;
  extern double swE_Typical[3];

  //the total number of source processes
  extern int nTotalSourceProcesses;

  //the sphere that representd the planet
  extern cInternalSphericalData *Planet;

  //init the model
  void Init_BeforeParser();
  void Init_AfterParser();

  //ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
  void SWMFdataPreProcessor(double *x,PIC::ICES::cDataNodeSWMF& data);


  //make coulumn integration
  void SodiumCoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
  void GetColumnIntegralLimbSubsolarPoint(double* ColumnIntegralVector,int ColumnIntegralVectorLength);

  void ColumnDensityIntegration_Tail(char *name);
  void ColumnDensityIntegration_Limb(char *name);
  void ColumnDensityIntegration_Map(char *fname,double dXmax,double dZmax,int nXpoints);
  void ColumnDensityIntegration_CircularMap(char *fname,double rmax,double dRmin,double dRmax,int nAzimuthPoints,SpiceDouble EphemerisTime);


  //Chamberlain Exosphere Model (Shen-1963-JAS,Valeille-2009)
  namespace ChamberlainExosphere {
    extern double *SpecieExobaseEscapeRate,*SpecieExobaseTemperature;
    extern bool ModelInitFlag;

    void Init(double *ExobaseEscapeRate,double *ExobaseTemperature);

    void PrintVariableList(FILE* fout);
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  }



  //Sampling
  namespace Sampling {

    //reference ground based observations
    namespace ReferenceGroundBasedObservations {

       struct cObservationTag {
         char TimeStamp[_MAX_STRING_LENGTH_PIC_];
         double TAA,PhaseAngle;
         SpiceDouble et;
         int nOutputFile;
       };

       extern cObservationTag RemoteObservationList[nReferenceGroundBasedObservations];
       void init();
       void OutputSampledData(SpiceDouble etStartInterval,SpiceDouble etFinishInterval,int nObjectOutputFile);
    }

    //offsets for sampling densities that are due to different sampling processes
    extern int SamplingDensity__ImpactVaporization_Offset;
    extern int SamplingDensity__PhotonStimulatedDesorption_Offset;
    extern int SamplingDensity__ThermalDesorption_Offset;
    extern int SamplingDensity__SolarWindSputtering_Offset;
    extern int CellSamplingDataOffset;

    //the field in the particle's data that keeps the id of the source process due to which the particle has beed produced
    extern long int ParticleData_SourceProcessID_Offset;
    extern long int ParticleData_OriginSurfaceElementNumber_Offset;

    //sample the planet's night side return flux
    extern double **PlanetNightSideReturnFlux;

    //total return flux; the rate of sticking to the planet's surface
    extern double *TotalPlanetReturnFlux,*PlanetSurfaceStickingRate;

    namespace OutputDataFile {
      //matrix for transformation SO->HCI coordinate frame (used only when prepare data files)
      extern SpiceDouble SO_to_HCI_TransformationMartix[6][6];

      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
    }

    namespace OutputSurfaceDataFile {
      void PrintVariableList(FILE* fout);
      void PrintTitle(FILE* fout);
      void PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);
    }

    //clean the model sampling buffers after output of the data file
    inline void FlushSamplingDataBuffers() {OutputSurfaceDataFile::flushCollectingSamplingBuffer(Planet);}

    //request the sampling buffer
    int RequestSamplingData(int offset);

    //set and get the source id
    inline int GetParticleSourceID(PIC::ParticleBuffer::byte *pData) {return *(int*)(pData+ParticleData_SourceProcessID_Offset);}
    inline void SetParticleSourceID(int id,PIC::ParticleBuffer::byte *pData) {*(int*)(pData+ParticleData_SourceProcessID_Offset)=id;}

    //set and get the surface elemetn number of the particle's origine
    inline int GetParicleOriginSurfaceElementNumber(PIC::ParticleBuffer::byte *pData) {return *(int*)(pData+ParticleData_OriginSurfaceElementNumber_Offset);}
    inline void SetParicleOriginSurfaceElementNumber(int el,PIC::ParticleBuffer::byte *pData) {*(int*)(pData+ParticleData_OriginSurfaceElementNumber_Offset)=el;}

    //sample particle's data
    void inline SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec) {
      int id;

      //determine if the particle is a sodium atom
      if (spec!=_NA_SPEC_) return;

      id=GetParticleSourceID((PIC::ParticleBuffer::byte*)ParticleData);

      switch (id) {
      case _EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_:
        *((double*)(SamplingBuffer+SamplingDensity__ImpactVaporization_Offset))+=LocalParticleWeight;
        break;
      case _EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_:
        *((double*)(SamplingBuffer+SamplingDensity__PhotonStimulatedDesorption_Offset))+=LocalParticleWeight;
        break;
      case _EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_:
        *((double*)(SamplingBuffer+SamplingDensity__ThermalDesorption_Offset))+=LocalParticleWeight;
        break;
      case _EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_:
        *((double*)(SamplingBuffer+SamplingDensity__SolarWindSputtering_Offset))+=LocalParticleWeight;
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is not found");
      }
    }


    void SampleModelData();
    void OutputSampledModelData(int DataOutputFileNumber);
  }


  //orbital motion of the Planet
  namespace OrbitalMotion {
    extern double AccumulatedPlanetRotation,TotalSimulationTime,TAA;

    //SPICE ephemeris time
    extern SpiceDouble et,lt;

    //direction to the Sun and the angle of the rotaio berween planetary axises and the direction to the Sun on the Z-plane
    extern double SunDirection_IAU_OBJECT[3];

    //matrixes for tranformation SO->IAU and IAU->SO coordinate frames
    extern SpiceDouble SO_to_IAU_TransformationMartix[6][6],IAU_to_SO_TransformationMartix[6][6];

    //parameters of orbital motion of Object
    extern double CoordinateFrameRotationRate;

    //the number of sub-intervals printed for a single interval between outputs of model data files in the pic.OrbitalData.dat
    extern int nOrbitalPositionOutputMultiplier;

    inline double GetCosineSubsolarAngle(double *x) { return x[0]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}

    //calculate TAA
    double GetTAA(SpiceDouble);
    double GetTAA(const char*);

    //calculate the phase angle
    double GetPhaseAngle(SpiceDouble);
    double GetPhaseAngle(const char*);

    namespace FrameRotation {
       void GetRotationAxis(double *RotationAxis,double &RotationAngle,const char *FrameName,double etStartRotation,double etFinishRotation);
       double GetRotationVector(double *RotationVector,const char *FrameName,double etStartRotation,double etFinishRotation);
       void GetRotationMatrix(double RotationMatrix[3][3],double *RotationAxis,double RotationAngle);
    }
  }

  //surface temeprature of the planet
  double GetSurfaceTemeprature(double CosSubSolarAngle,double *x_SO);



  //sources
  namespace SourceProcesses {

    //generate particle position and velocity
  //generate particle properties
  inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,cSingleVariableDiscreteDistribution *SurfaceInjectionDistribution,cSingleVariableDistribution *EnergyDistribution) {
    double ExternalNormal[3];

    //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
    double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
    int nZenithElement,nAzimuthalElement;
    unsigned int idim,el;

    el=SurfaceInjectionDistribution->DistributeVariable();
    Exosphere::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
    Exosphere::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

    x_LOCAL_IAU_OBJECT[0]=sphereRadius*ExternalNormal[0];
    x_LOCAL_IAU_OBJECT[1]=sphereRadius*ExternalNormal[1];
    x_LOCAL_IAU_OBJECT[2]=sphereRadius*ExternalNormal[2];

    /*
    ExternalNormal[0]*=-1.0;
    ExternalNormal[1]*=-1.0;
    ExternalNormal[2]*=-1.0;
    */

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
    double c=0.0,rVel=0.0,lVel[3];
    double Speed=sqrt(EnergyDistribution->DistributeVariable()*2.0/_NA__MASS_);


    for (idim=0;idim<3;idim++) {
      lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      rVel+=pow(lVel[idim],2);

      c+=ExternalNormal[idim]*lVel[idim];
    }

    rVel=Speed/sqrt(rVel);

    if (c>0.0) {
      //the distributed velocity vector is directed into the domain
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=lVel[idim]*rVel;
    }
    else {
      //the distributed velocity vector is directed into the planet -> redirect it
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=(lVel[idim]-2.0*c*ExternalNormal[idim])*rVel;
    }

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

    namespace ThermalDesorption {
      extern double uThermal;
      extern double VibrationalFrequency;

      extern double SourceRate,maxLocalSourceRate;
      extern double CalculatedTotalSodiumSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],temp,SodiumSurfaceElementPopulation;

        SodiumSurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SodiumSurfaceElementPopulation[SurfaceElement];
        if (SodiumSurfaceElementPopulation<0.0) return 0.0;

        memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];

        //calculate local temeprature
        //get the local position of the surface element in SO_FRAME
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3];

        x_LOCAL_IAU_OBJECT[0]=_RADIUS_(_TARGET_)*norm[0];
        x_LOCAL_IAU_OBJECT[1]=_RADIUS_(_TARGET_)*norm[1];
        x_LOCAL_IAU_OBJECT[2]=_RADIUS_(_TARGET_)*norm[2];

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


        temp=GetSurfaceTemeprature(res,x_LOCAL_SO_OBJECT);
        res=VibrationalFrequency*SodiumSurfaceElementPopulation*exp(-uThermal/(Kbol*temp));

        return res;
      }


      inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double* v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        double ExternalNormal[3],CosSubSolarAngle;

        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
        int nZenithElement,nAzimuthalElement;
        unsigned int el;

        el=SurfaceInjectionDistribution.DistributeVariable();
        Exosphere::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
        Exosphere::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

        CosSubSolarAngle=(Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[0]*ExternalNormal[0])+
            (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[1]*ExternalNormal[1])+
            (Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[2]*ExternalNormal[2]);

        x_LOCAL_IAU_OBJECT[0]=sphereRadius*ExternalNormal[0];
        x_LOCAL_IAU_OBJECT[1]=sphereRadius*ExternalNormal[1];
        x_LOCAL_IAU_OBJECT[2]=sphereRadius*ExternalNormal[2];

        ExternalNormal[0]*=-1.0;
        ExternalNormal[1]*=-1.0;
        ExternalNormal[2]*=-1.0;

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

        SurfaceTemperature=GetSurfaceTemeprature(CosSubSolarAngle,x_LOCAL_SO_OBJECT);
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,_NA_SPEC_);


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



    }

    namespace PhotonStimulatedDesorption {
      extern double PhotonFlux_1AU;
      extern double CrossSection;

      extern double minInjectionEnergy;
      extern double maxInjectionEnergy;

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],SodiumSurfaceElementPopulation;

        SodiumSurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SodiumSurfaceElementPopulation[SurfaceElement];
        if (SodiumSurfaceElementPopulation<0.0) return 0.0;

        memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];


        res*=(res>0.0) ? SodiumSurfaceElementPopulation*CrossSection*PhotonFlux_1AU*pow(_AU_/Exosphere::xObjectRadial,2.0) : 0.0;
        return res;
      }

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        return Exosphere::SourceProcesses::GenerateParticleProperties(x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }
    }

    namespace ImpactVaporization {
      extern double SourceRate;
      extern double HeliocentricDistance;
      extern double SourceRatePowerIndex;
      extern double SourceTemeprature;

      double GetTotalProductionRate(int spec,void *SphereDataPointer);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode) {
        unsigned int idim;
        double r=0.0,vbulk[3]={0.0,0.0,0.0},ExternalNormal[3];


        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
        SpiceDouble xform[6][6];

        memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

        //Geenrate new particle position
        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=pow(ExternalNormal[idim],2);
        }

        r=sqrt(r);

        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]/=r;
          x_LOCAL_IAU_OBJECT[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
        }

        //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
        x_LOCAL_SO_OBJECT[0]=xform[0][0]*x_LOCAL_IAU_OBJECT[0]+xform[0][1]*x_LOCAL_IAU_OBJECT[1]+xform[0][2]*x_LOCAL_IAU_OBJECT[2];
        x_LOCAL_SO_OBJECT[1]=xform[1][0]*x_LOCAL_IAU_OBJECT[0]+xform[1][1]*x_LOCAL_IAU_OBJECT[1]+xform[1][2]*x_LOCAL_IAU_OBJECT[2];
        x_LOCAL_SO_OBJECT[2]=xform[2][0]*x_LOCAL_IAU_OBJECT[0]+xform[2][1]*x_LOCAL_IAU_OBJECT[1]+xform[2][2]*x_LOCAL_IAU_OBJECT[2];


        //determine if the particle belongs to this processor
        startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
        if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

        //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SourceTemeprature,ExternalNormal,_NA_SPEC_);


/*
//DEBUG -> injected velocity is normal to the surface

for (int i=0;i<3;i++)  v_LOCAL_IAU_OBJECT[i]=-ExternalNormal[i]*4.0E3;
//END DEBUG
*/


        //transform the velocity vector to the coordinate frame 'MSGR_SO'
        v_LOCAL_SO_OBJECT[0]=xform[3][0]*x_LOCAL_IAU_OBJECT[0]+xform[3][1]*x_LOCAL_IAU_OBJECT[1]+xform[3][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[3][3]*v_LOCAL_IAU_OBJECT[0]+xform[3][4]*v_LOCAL_IAU_OBJECT[1]+xform[3][5]*v_LOCAL_IAU_OBJECT[2];

        v_LOCAL_SO_OBJECT[1]=xform[4][0]*x_LOCAL_IAU_OBJECT[0]+xform[4][1]*x_LOCAL_IAU_OBJECT[1]+xform[4][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[4][3]*v_LOCAL_IAU_OBJECT[0]+xform[4][4]*v_LOCAL_IAU_OBJECT[1]+xform[4][5]*v_LOCAL_IAU_OBJECT[2];

        v_LOCAL_SO_OBJECT[2]=xform[5][0]*x_LOCAL_IAU_OBJECT[0]+xform[5][1]*x_LOCAL_IAU_OBJECT[1]+xform[5][2]*x_LOCAL_IAU_OBJECT[2]+
            xform[5][3]*v_LOCAL_IAU_OBJECT[0]+xform[5][4]*v_LOCAL_IAU_OBJECT[1]+xform[5][5]*v_LOCAL_IAU_OBJECT[2];

        memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
        memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
        memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
        memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));

        return true;
      }

    }

    namespace SolarWindSputtering {
      extern double Yield;
      extern double minInjectionEnergy;
      extern double maxInjectionEnergy;

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double norm_IAU_OBJECT[3],norm_SO_OBJECT[3];

        if (((cInternalSphericalData*)SphereDataPointer)->SodiumSurfaceElementPopulation[SurfaceElement]<=0.0) return 0.0;

        //get the normal to the surface element in 'IAU' frame
        memcpy(norm_IAU_OBJECT,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        //convert the normal vector to the 'SO' frame
        norm_SO_OBJECT[0]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*norm_IAU_OBJECT[2]);

        norm_SO_OBJECT[1]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*norm_IAU_OBJECT[2]);

        norm_SO_OBJECT[2]=
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*norm_IAU_OBJECT[0])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*norm_IAU_OBJECT[1])+
            (OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*norm_IAU_OBJECT[2]);

        //get the surface element that is pointer by the vectorm norm_SO_OBJECT
        long int nZenithElement,nAzimuthalElement,nd;

        ((cInternalSphericalData*)SphereDataPointer)->GetSurfaceElementProjectionIndex(norm_SO_OBJECT,nZenithElement,nAzimuthalElement);
        nd=((cInternalSphericalData*)SphereDataPointer)->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

        //return the local source rate
        return ((cInternalSphericalData*)SphereDataPointer)->SolarWindSurfaceFlux[nd]*Yield*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement];
      }

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        return Exosphere::SourceProcesses::GenerateParticleProperties(x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }
    }


    double totalProductionRate(int spec,void *SphereDataPointer);
    long int InjectionBoundaryModel(void *SphereDataPointer);
  }



  //combine together the surface area densities and recalcualte the
  void ExchangeSurfaceAreaDensity();


  //Interaction of the particles with the surface
  namespace SurfaceInteraction {


    //sampling the surface area density of the sticking species
#define _OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_             3

#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_DOWN_           0
#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_UP_             1
#define _OBJECT_SURFACE_STICKING_SAMPLING_OFFSET__AREA_NUMBER_DENSITY_ 2


    //sodium/surface interaction model
    extern const double AccomodationCoefficient;

    //sticking probability of sodium atoms
    double SodiumStickingProbability(double& ReemissionParticleFraction,double Temp);

    //model of the interaction between particles and the planetary surface
    int ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);
  }










}



#endif /* OBJECT_H_ */
