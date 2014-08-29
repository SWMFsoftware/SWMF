//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * Mercury.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#ifndef _EXOSPHERE_
#define _EXOSPHERE_



#include "SingleVariableDistribution.h"
#include "SingleVariableDiscreteDistribution.h"
#include "constants.h"

#include "Exosphere.dfn"






//user defined settings of the exospheric model
//#include "UserDefinition.Exosphere.h"
#include "Na.h"

//define the symbolic id of source processes
static const char _EXOSPHERE__SOURCE_SYMBOLIC_ID_[][100]={"ImpactVaposization","PhotonStimulatedDesorption","ThermalDesorption","SolarWindSputtering"};

//the default value for for the list of the SPICE kernels that will be loaded
static const int nFurnishedSPICEkernels=0;
static const char SPICE_Kernels[][_MAX_STRING_LENGTH_PIC_]={""};

//the default values of the list of the referenced ground based observations
static const int nReferenceGroundBasedObservations=0;
static const char ReferenceGroundBasedObservationTime[][_MAX_STRING_LENGTH_PIC_]={""};

//the default location of the SPICE kernels
static const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="";

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif


namespace Exosphere {

  //the name os the simulated object and related coordinate frames
  extern char ObjectName[_MAX_STRING_LENGTH_PIC_],IAU_FRAME[_MAX_STRING_LENGTH_PIC_],SO_FRAME[_MAX_STRING_LENGTH_PIC_];


  //simulation date and position of Object at the time of simulations
  //const char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_]="2011-04-13T00:00:00"; //"2001-03-01T00:00:00";  ////"2011-01-01T00:00:00";
  //extern char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_];

  static const char SimulationStartTimeString[]="2000-01-01T00:00:00";

  extern double xObject_HCI[3],vObject_HCI[3],xEarth_HCI[3],vEarth_HCI[3],xEarth_SO[3],vEarth_SO[3],xSun_SO[3],vSun_SO[3];
  extern double vObjectRadial,xObjectRadial;

  extern double vObject_SO_FROZEN[3]; //Velocity of Object in inertial frame, which axis coinsides with that of SO at a given time
  extern double RotationVector_SO_FROZEN[3],RotationRate_SO_FROZEN; //the vectors of the direction of rotation and the rotation rate of the SO in SO_FROZEN


  //typical solar wind conditions far from the planet
  static const double /*Exosphere_*/swVelocity_Typical[]={0.0,0.0,0.0};
  static const double /*Exosphere_*/swB_Typical[]={0.0,0.0,0.0};
  static const double /*Exosphere_*/swTemperature_Typical=0.0;
  static const double /*Exosphere_*/swNumberDensity_Typical=0.0;
  extern double swE_Typical[3];

  //the total number of source processes
  extern int nTotalSourceProcesses;

  //the array of flags that defines wether the source process change the surface aboundance of the volatile
  static const bool Source_DeplitSurfaceSpeciesAbundance_Flag[]={true};


  //the sphere that representd the planet
  extern cInternalSphericalData *Planet;

  //init the model
  void Init_BeforeParser();
  void Init_AfterParser();

  //ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
  void SWMFdataPreProcessor(double *x,PIC::CPLR::ICES::cDataNodeSWMF& data);


  //make coulumn integration
  namespace ColumnIntegral {
    void CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
    int GetVariableList(char *vlist=NULL);
    void ProcessColumnIntegrationVector(double *res,int resLength);



    void GetSubsolarPointDirection(double *LimbDirection,double *EarthPosition);
    void Tail(char *name);
    void Limb(char *name);
//    void Map(char *fname,double dXmax,double dZmax,int nXpoints);
    void CircularMap(char *fname,double rmax,double dRmin,double dRmax,int nAzimuthPoints,SpiceDouble EphemerisTime);
  }

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

    //sample the source rates
    extern double CalculatedSourceRate[PIC::nTotalSpecies][1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_];

    //offsets for sampling densities that are due to different sampling processes
    extern int SamplingDensityOffset[PIC::nTotalSpecies*(1+_EXOSPHERE__SOURCE_MAX_ID_VALUE_)];
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
//      if (spec!=_NA_SPEC_) return;

      id=GetParticleSourceID((PIC::ParticleBuffer::byte*)ParticleData);
      *((double*)(SamplingBuffer+spec+SamplingDensityOffset[id]))+=LocalParticleWeight;
/*
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
      */
    }


    void SampleModelData();
    void OutputSampledModelData(int DataOutputFileNumber);

    typedef void (*fUserDefinedAdditionalData_VariableList_OutputSampledModelData)(FILE*);
    typedef void (*fUserDefinedAdditionalData_OutputSampledModelData)(FILE*,int);
    extern fUserDefinedAdditionalData_VariableList_OutputSampledModelData UserDefinedAdditionalData_VariableList_OutputSampledModelData;
    extern fUserDefinedAdditionalData_OutputSampledModelData UserDefinedAdditionalData_OutputSampledModelData;

    void inline SetUserDefinedAdditionalOutputSampledModelDataFunctions(fUserDefinedAdditionalData_VariableList_OutputSampledModelData t0,fUserDefinedAdditionalData_OutputSampledModelData t1) {
      UserDefinedAdditionalData_VariableList_OutputSampledModelData=t0;
      UserDefinedAdditionalData_OutputSampledModelData=t1;
    }
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
    double GetTAA(const char* TargetName, const char* CenterBodyName, double CenterBodyMass, SpiceDouble EphemerisTime);

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
  inline bool GenerateParticleProperties(double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,cSingleVariableDiscreteDistribution<int> *SurfaceInjectionDistribution,cSingleVariableDistribution<int> *EnergyDistribution) {
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
      static const double ThermalDesorption_uThermal[]={0.0};
      static const double ThermalDesorption_VibrationalFrequency[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];
//      extern double CalculatedTotalSodiumSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return SourceRate[spec];}

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],temp,SodiumSurfaceElementPopulation;

        SodiumSurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement];
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
        res=ThermalDesorption_VibrationalFrequency[spec]*SodiumSurfaceElementPopulation*exp(-ThermalDesorption_uThermal[spec]/(Kbol*temp));

        return res;
      }


      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double* v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
        double ExternalNormal[3],CosSubSolarAngle;

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif


        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
        double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
        int nZenithElement,nAzimuthalElement;
        unsigned int el;

        el=SurfaceInjectionDistribution[spec].DistributeVariable();
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



    }

    namespace PhotonStimulatedDesorption {
      static const double PhotonStimulatedDesorption_PhotonFlux_1AU=0.0;
      static const double PhotonStimulatedDesorption_CrossSection[]={0.0};

      static const double PhotonStimulatedDesorption_minInjectionEnergy[]={0.0};
      static const double PhotonStimulatedDesorption_maxInjectionEnergy[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],SurfaceElementPopulation;

        SurfaceElementPopulation=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement];
        if (SurfaceElementPopulation<0.0) return 0.0;

        memcpy(sun,Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];


        res*=(res>0.0) ? SurfaceElementPopulation*PhotonStimulatedDesorption_CrossSection[spec]*PhotonStimulatedDesorption_PhotonFlux_1AU*pow(_AU_/Exosphere::xObjectRadial,2.0) : 0.0;
        return res;
      }

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return SourceRate[spec];}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution[PIC::nTotalSpecies];
      double EnergyDistributionFunction(double e,int *spec);

      //generate particle properties
      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif

        if (BoundaryElementType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: particle ejection from a non-spehtical body is not implemented");

        cInternalSphericalData* Sphere=(cInternalSphericalData*)BoundaryElement;

        return Exosphere::SourceProcesses::GenerateParticleProperties(x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,SurfaceInjectionDistribution+spec,EnergyDistribution+spec);
      }
    }

    namespace ImpactVaporization {
      static const double ImpactVaporization_SourceRate[]={0.0};
      static const double ImpactVaporization_HeliocentricDistance=1.0*_AU_;
      static const double ImpactVaporization_SourceRatePowerIndex=0.0;
      static const double ImpactVaporization_SourceTemeprature[]={0.0};

      double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
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
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,ImpactVaporization_SourceTemeprature[spec],ExternalNormal,spec);


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

        //set up the intermal energy if needed
#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_

  #if _PIC_INTERNAL_DEGREES_OF_FREEDOM__TR_RELAXATION_MODE_  == _PIC_MODE_ON_
        PIC::IDF::InitRotTemp(ImpactVaporization_SourceTemeprature[spec],tempParticleData);
  #endif

  #if _PIC_INTERNAL_DEGREES_OF_FREEDOM__VT_RELAXATION_MODE_  == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
  #endif

#endif

        return true;
      }

    }

    namespace SolarWindSputtering {
      static const double SolarWindSputtering_Yield[]={0.0};
      static const double SolarWindSputtering_minInjectionEnergy[]={0.0};
      static const double SolarWindSputtering_maxInjectionEnergy[]={0.0};

      extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
      double GetSurfaceElementProductionRate(int nElement,int *spec);

      //evaluate nemerically the source rate
//      extern double CalculatedTotalSodiumSourceRate;

      inline double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double norm_IAU_OBJECT[3],norm_SO_OBJECT[3];

        if (((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[spec][SurfaceElement]<=0.0) return 0.0;

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
        return ((cInternalSphericalData*)SphereDataPointer)->SolarWindSurfaceFlux[nd]*SolarWindSputtering_Yield[spec]*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement];
      }

      inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return SourceRate[spec];}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution[PIC::nTotalSpecies];
      double EnergyDistributionFunction(double e,int *spec);

      //generate particle properties
      inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {

#if _PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_
        exit(__LINE__,__FILE__,"Error: not implemented");
#endif

        if (BoundaryElementType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: particle ejection from a non-spehtical body is not implemented");

        cInternalSphericalData* Sphere=(cInternalSphericalData*)BoundaryElement;

        return Exosphere::SourceProcesses::GenerateParticleProperties(x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,Sphere,SurfaceInjectionDistribution+spec,EnergyDistribution+spec);
      }
    }


    double totalProductionRate(int spec,int BoundaryElementType,void *BoundaryElement);
    long int InjectionBoundaryModel(int BoundaryElementType,void *BoundaryElement);
    long int InjectionBoundaryModel(int spec,int BoundaryElementType,void *BoundaryElement);

    long int InjectionBoundaryModelLimited(void *SphereDataPointer);
    long int InjectionBoundaryModelLimited(int spec,void *SphereDataPointer);
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
    static const double AccomodationCoefficient[]={0.0};

    //sticking probability of sodium atoms
    double StickingProbability(int spec,double& ReemissionParticleFraction,double Temp);

    //model of the interaction between particles and the planetary surface
    int ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);
  }










}



#endif /* OBJECT_H_ */
