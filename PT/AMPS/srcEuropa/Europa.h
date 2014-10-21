/*
 * Europa.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#ifndef EUROPA_H_
#define EUROPA_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


#include "Exosphere.h"


#include "Na.h"

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif

#include "SingleVariableDistribution.h"
#include "SingleVariableDiscreteDistribution.h"

#include "constants.h"

/*
//path to the SPICE Kernels directory
const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/SPICE/Kernels"; //"/Users/rubinmar/Codes/CSPICE/Kernels/GALILEO";

//SPICE Kernels to be loaded
const int nFurnishedSPICEkernels=6;
const char SPICE_Kernels[nFurnishedSPICEkernels][_MAX_STRING_LENGTH_PIC_]={"GALILEO/MK00062B.TSC","NAIF/naif0010.tls","PCK/PCK00006.TPC","GALILEO/s980326B.bsp","GALILEO/pk96030a.tpc","GALILEO/galileo.tf"};
*/



//descriptors of the source processes
#define _EUROPA_SOURCE__ON_    0
#define _EUROPA_SOURCE__OFF_   1

#define _EUROPA_SOURCE__ID__IMPACT_VAPORIZATION_              0
#define _EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_     1
#define _EUROPA_SOURCE__ID__THERMAL_DESORPTION_               2
#define _EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_            3
#define _EUROPA_SOURCE__ID__MAGNETOSPHERIC_PLS                4
#define _EUROPA_SOURCE__ID__MAGNETOSPHERIC_EPD                5

//define wich of the source processes are active
#define _EUROPA_SOURCE__IMPACT_VAPORIZATION_              _EUROPA_SOURCE__OFF_
#define _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_     _EUROPA_SOURCE__OFF_
#define _EUROPA_SOURCE__THERMAL_DESORPTION_               _EUROPA_SOURCE__OFF_
#define _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_            _EUROPA_SOURCE__OFF_


//integration mode: steady state/time dependent
#define _EUROPA_INTEGRATION_MODE__STEADY_STATE_    0
#define _EUROPA_INTEGRATION_MODE__TIME_DEPENDENT_  1


#define _EUROPA_INTEGRATION_MODE_ _EUROPA_INTEGRATION_MODE__TIME_DEPENDENT_

//maximum value of the source ID number
#define _EUROPA_SOURCE_MAX_ID_VALUE_ 0

#if _EUROPA_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EUROPA_SOURCE__ON_
#undef _EUROPA_SOURCE_MAX_ID_VALUE_
#define _EUROPA_SOURCE_MAX_ID_VALUE_ _EUROPA_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_
#endif

#if _EUROPA_SOURCE__THERMAL_DESORPTION_ == _EUROPA_SOURCE__ON_
#undef _EUROPA_SOURCE_MAX_ID_VALUE_
#define _EUROPA_SOURCE_MAX_ID_VALUE_ _EUROPA_SOURCE__ID__THERMAL_DESORPTION_
#endif

#if _EUROPA_SOURCE__SOLAR_WIND_SPUTTERING_  == _EUROPA_SOURCE__ON_
#undef _EUROPA_SOURCE_MAX_ID_VALUE_
#define _EUROPA_SOURCE_MAX_ID_VALUE_ _EUROPA_SOURCE__ID__SOLAR_WIND_SPUTTERING_
#endif


#define _EUROPA__SPUTTERING_ION_SOURCE__SWMF_PLASMA_FLOW_   0
#define _EUROPA__SPUTTERING_ION_SOURCE__AMPS_KINETIC_IONS_  1

#define _EUROPA__SPUTTERING_ION_SOURCE_ _EUROPA__SPUTTERING_ION_SOURCE__SWMF_PLASMA_FLOW_

namespace Europa {
  using namespace Exosphere;


  //the total injection rate
  extern double TotalInjectionRate;


  //simulation date and position of Europa at the time of simulations
//  const char SimulationStartTimeString[_MAX_STRING_LENGTH_PIC_]="1996-12-19T06:00:00";  //closest approach "1996-12-19T06:21:00";
  extern double xEuropa[3],vEuropa[3],xEarth[3],vEarth[3];
  extern double vEuropaRadial,xEuropaRadial;

/*  //typical solar wind conditions far from the planet
  static const double swVelocity_Typical[3]={-420.0E3,0.0,0.0};
  static const double swB_Typical[3]={-12.9E-9,4.71E-9,10.29E-9};
  static const double swTemperature_Typical=0.174e6,swNumberDensity_Typical=60.0E6;
  extern double swE_Typical[3];*/

  //the total number of source processes
  //extern int nTotalSourceProcesses;

  //the sphere that representd the planet
  //extern cInternalSphericalData *Planet;

  //init the model
  void Init_BeforeParser();
//  void Init_AfterParser();

  //ICES data preprocessor -> set up typical values of the solar wind in the regions where the SWMF values have not been found
  void SWMFdataPreProcessor(double *x,PIC::CPLR::ICES::cDataNodeSWMF& data);


  //make coulumn integration
  void SodiumCoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
  void ColumnDensityIntegration_Tail(char *name);
  void ColumnDensityIntegration_Map(char *name);


  //injection parameters of the thermal ions
  const double Thermal_OPlus_NumberDensity=35.0E6;
  const double Thermal_OPlus_Temperature=1.044E6;
  const double Thermal_OPlus_BulkVelocity[3]={90300.0,0.0,0.0};

/*  namespace ThermalIon {
    namespace O {
      const double NumberDensity=35.0E6;
      const double Temeprature=1.044E6;
      extern double BulkVelocity[3];
    }
  }*/

  namespace Sampling {
    using namespace Exosphere::Sampling;

    namespace O2InjectionSpeed {
      const int nSampleIntervals=1000;
      const double vMax=10.0E3;
      const double dv=vMax/nSampleIntervals;

      extern double SamplingBuffer[nSampleIntervals];

      void flush();
      void OutputSampledModelData(int DataOutputFileNumber);
    }
  }


  //orbital motion of the Planet
  namespace OrbitalMotion {
    extern double AccumulatedPlanetRotation,TotalSimulationTime,TAA;

    //SPICE ephemeris time
    extern SpiceDouble et,lt;

    //direction to the Sun and the angle of the rotation between planetary axes and the direction to the Sun on the Z-plane
    extern double SunDirection_IAU_EUROPA[3],PlanetAxisToSunRotationAngle;

    //matrixes for transformation GALL_EPHIOD->IAU and IAU->GALL_EPHIOD coordinate frames
    extern SpiceDouble GALL_EPHIOD_to_IAU_TransformationMartix[6][6],IAU_to_GALL_EPHIOD_TransformationMartix[6][6];
    extern double IAU_to_LS_TransformationMatrix[3][3];


    //transformation matrix: IAU -> GALL_ESOM (used for calcualtion of the surface temeprature)
    extern SpiceDouble IAU_to_GALL_ESOM_TransformationMatrix[6][6];

    //parameters of orbital motion of Europa
    extern double CoordinateFrameRotationRate;

    inline double GetCosineSubsolarAngle(double *x) { return x[0]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}
  }

  //surface temeprature of the planet
  /*
  inline double GetSurfaceTemeprature(double CosSubSolarAngle) {
    static const double Tn=110.0;
    static const double Td0_Aphelion=590.0,Td0_Perihelion=725.0;
    static const double TAA_Aphelion=Pi,TAA_Perihelion=0.0;
    static const double alpha=(Td0_Aphelion-Td0_Perihelion)/(TAA_Aphelion-TAA_Perihelion);

    double Td,Angle;

    Angle=(OrbitalMotion::TAA<Pi) ? OrbitalMotion::TAA : 2.0*Pi-OrbitalMotion::TAA;
    Td=Td0_Perihelion+alpha*(Angle-TAA_Perihelion);

    return (CosSubSolarAngle>0.0) ? Tn+(Td-Tn)*pow(CosSubSolarAngle,0.25) : Tn;
  }
  */

  //Surface temeprature
  inline double GetSurfaceTemeprature(double *x_IAU) {
/*    double x_GALL_ESOM[3];
    int i,j;

    //transform coordinates into x_GALL_ESOM
    for (i=0;i<3;i++) {
      x_GALL_ESOM[i]=0.0;

      for (j=0;j<3;j++) {
        x_GALL_ESOM[i]+=Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[i][j]*x_IAU[j];
      }
    }*/



    //calculate the surface itself
    double cosSubSolarAngle,r,cosLat,Temp;

    r=sqrt(x_IAU[0]*x_IAU[0]+x_IAU[1]*x_IAU[1]+x_IAU[2]*x_IAU[2]);

    cosSubSolarAngle=(Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0]*x_IAU[0]+
        Europa::OrbitalMotion::SunDirection_IAU_EUROPA[1]*x_IAU[1]+
        Europa::OrbitalMotion::SunDirection_IAU_EUROPA[2]*x_IAU[2])/r;

    cosLat=sqrt(x_IAU[0]*x_IAU[0]+x_IAU[1]*x_IAU[1])/r;

    Temp=40.0+50.0*pow(cosLat,0.25)+40.0*cosSubSolarAngle;

    return Temp;
  }



  //sources
  /*namespace SourceProcesses {

    //generate particle position and velocity
  //generate particle properties
  inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,cSingleVariableDiscreteDistribution<int> *SurfaceInjectionDistribution,cSingleVariableDistribution<int> *EnergyDistribution) {
    double ExternalNormal[3];

    //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_EUROPA'
    double x_LOCAL_IAU_EUROPA[3],x_LOCAL_GALL_EPHIOD_EUROPA[3],v_LOCAL_IAU_EUROPA[3],v_LOCAL_GALL_EPHIOD_EUROPA[3];
    int nZenithElement,nAzimuthalElement;
    unsigned int idim,el;

    el=SurfaceInjectionDistribution->DistributeVariable();
    Europa::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
    Europa::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

    x_LOCAL_IAU_EUROPA[0]=sphereX0[0]+sphereRadius*ExternalNormal[0];
    x_LOCAL_IAU_EUROPA[1]=sphereX0[1]+sphereRadius*ExternalNormal[1];
    x_LOCAL_IAU_EUROPA[2]=sphereX0[2]+sphereRadius*ExternalNormal[2];

    ExternalNormal[0]*=-1.0;
    ExternalNormal[1]*=-1.0;
    ExternalNormal[2]*=-1.0;

    //transfer the position into the coordinate frame related to the rotating coordinate frame 'GALL_EPHIOD'
    x_LOCAL_GALL_EPHIOD_EUROPA[0]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][2]*x_LOCAL_IAU_EUROPA[2]);

    x_LOCAL_GALL_EPHIOD_EUROPA[1]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][2]*x_LOCAL_IAU_EUROPA[2]);

    x_LOCAL_GALL_EPHIOD_EUROPA[2]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][2]*x_LOCAL_IAU_EUROPA[2]);


    //determine if the particle belongs to this processor
    startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_GALL_EPHIOD_EUROPA,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

    //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_EUROPA'
    double c=0.0,rVel=0.0,lVel[3];
    double Speed=sqrt(EnergyDistribution->DistributeVariable()*2.0/_O2__MASS_);


    for (idim=0;idim<3;idim++) {
      lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      rVel+=pow(lVel[idim],2);

      c+=ExternalNormal[idim]*lVel[idim];
    }

    rVel=Speed/sqrt(rVel);

    if (c>0.0) {
      //the distributed velocity vector is directed into the domain
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_EUROPA[idim]=lVel[idim]*rVel;
    }
    else {
      //the distributed velocity vector is directed into the planet -> redirect it
      for (idim=0;idim<3;idim++) v_LOCAL_IAU_EUROPA[idim]=(lVel[idim]-2.0*c*ExternalNormal[idim])*rVel;
    }

    //transform the velocity vector to the coordinate frame 'GALL_EPHIOD'
    v_LOCAL_GALL_EPHIOD_EUROPA[0]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][2]*x_LOCAL_IAU_EUROPA[2])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][3]*v_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][4]*v_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][5]*v_LOCAL_IAU_EUROPA[2]);

    v_LOCAL_GALL_EPHIOD_EUROPA[1]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][2]*x_LOCAL_IAU_EUROPA[2])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][3]*v_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][4]*v_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][5]*v_LOCAL_IAU_EUROPA[2]);

    v_LOCAL_GALL_EPHIOD_EUROPA[2]=
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][0]*x_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][1]*x_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][2]*x_LOCAL_IAU_EUROPA[2])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][3]*v_LOCAL_IAU_EUROPA[0])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][4]*v_LOCAL_IAU_EUROPA[1])+
        (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][5]*v_LOCAL_IAU_EUROPA[2]);

    memcpy(x_GALL_EPHIOD_EUROPA,x_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
    memcpy(x_IAU_EUROPA,x_LOCAL_IAU_EUROPA,3*sizeof(double));
    memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));

    return true;
  }

    namespace ThermalDesorption {
      const double uThermal=1.85*eV2J;
      const double VibrationalFrequency=1.0E13;

      extern double SourceRate,maxLocalSourceRate;
      extern double CalculatedTotalSodiumSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      inline double GetLocalProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],temp,SurfaceDensity;

        SurfaceDensity=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement];
        if (SurfaceDensity<0.0) return 0.0;

        memcpy(sun,Europa::OrbitalMotion::SunDirection_IAU_EUROPA,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];


        //get the coordinate of the middle point of the surface element
        double xMiddle[3];
        ((cInternalSphericalData*)SphereDataPointer)->GetSurfaceElementMiddlePoint(xMiddle,SurfaceElement);


        temp=GetSurfaceTemeprature(xMiddle);
        res=VibrationalFrequency*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement]*exp(-uThermal/(Kbol*temp));

        return res;
      }


      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        double ExternalNormal[3];

        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_EUROPA'
        double x_LOCAL_IAU_EUROPA[3],x_LOCAL_GALL_EPHIOD_EUROPA[3],v_LOCAL_IAU_EUROPA[3],v_LOCAL_GALL_EPHIOD_EUROPA[3];
        int nZenithElement,nAzimuthalElement;
        unsigned int el;

        el=SurfaceInjectionDistribution.DistributeVariable();
        Europa::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
        Europa::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

        x_LOCAL_IAU_EUROPA[0]=sphereX0[0]+sphereRadius*ExternalNormal[0];
        x_LOCAL_IAU_EUROPA[1]=sphereX0[1]+sphereRadius*ExternalNormal[1];
        x_LOCAL_IAU_EUROPA[2]=sphereX0[2]+sphereRadius*ExternalNormal[2];

        ExternalNormal[0]*=-1.0;
        ExternalNormal[1]*=-1.0;
        ExternalNormal[2]*=-1.0;

        //transfer the position into the coordinate frame related to the rotating coordinate frame 'GALL_EPHIOD'
        x_LOCAL_GALL_EPHIOD_EUROPA[0]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][2]*x_LOCAL_IAU_EUROPA[2]);

        x_LOCAL_GALL_EPHIOD_EUROPA[1]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][2]*x_LOCAL_IAU_EUROPA[2]);

        x_LOCAL_GALL_EPHIOD_EUROPA[2]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][2]*x_LOCAL_IAU_EUROPA[2]);


        //determine if the particle belongs to this processor
        startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_GALL_EPHIOD_EUROPA,startNode);
        if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

        //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_EUROPA'
        double CosSubSolarAngle,SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};

        CosSubSolarAngle=sqrt(1.0-ExternalNormal[2]*ExternalNormal[2]);
        SurfaceTemperature=GetSurfaceTemeprature(x_LOCAL_IAU_EUROPA);
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_EUROPA,vbulk,SurfaceTemperature,ExternalNormal,_O2_SPEC_);


        //transform the velocity vector to the coordinate frame 'GALL_EPHIOD'
        v_LOCAL_GALL_EPHIOD_EUROPA[0]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][2]*x_LOCAL_IAU_EUROPA[2])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][3]*v_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][4]*v_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][5]*v_LOCAL_IAU_EUROPA[2]);

        v_LOCAL_GALL_EPHIOD_EUROPA[1]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][2]*x_LOCAL_IAU_EUROPA[2])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][3]*v_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][4]*v_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][5]*v_LOCAL_IAU_EUROPA[2]);

        v_LOCAL_GALL_EPHIOD_EUROPA[2]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][0]*x_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][1]*x_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][2]*x_LOCAL_IAU_EUROPA[2])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][3]*v_LOCAL_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][4]*v_LOCAL_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][5]*v_LOCAL_IAU_EUROPA[2]);

        memcpy(x_GALL_EPHIOD_EUROPA,x_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
        memcpy(x_IAU_EUROPA,x_LOCAL_IAU_EUROPA,3*sizeof(double));
        memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));

        return true;
      }



    }

    namespace PhotonStimulatedDesorption {
      const double PhotonFlux_1AU=3.0E15*1.0E+4;
      const double CrossSection=2.0E-20*1.0E-4;

      const double minInjectionEnergy=pow(10.0,2)*_O2__MASS_/2.0;
      const double maxInjectionEnergy=pow(10.0E3,2)*_O2__MASS_/2.0;

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetLocalProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double res=0.0,norm[3],sun[3],SurfaceDensity;

        SurfaceDensity=((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement];
        if (SurfaceDensity<0.0) return 0.0;

        memcpy(sun,Europa::OrbitalMotion::SunDirection_IAU_EUROPA,3*sizeof(double));
        memcpy(norm,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));

        for (int idim=0;idim<DIM;idim++) res+=sun[idim]*norm[idim];


        res*=(res>0.0) ? ((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement]*CrossSection*PhotonFlux_1AU*pow(_AU_/Europa::xEuropaRadial,2.0) : 0.0;
        return res;
      }

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        return Europa::SourceProcesses::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }
    }

    namespace ImpactVaporization {
      const double SourceRate=2.6E23;
      const double HeliocentricDistance=0.387098*_AU_;
      const double SourceRatePowerIndex=1.9;

      const double SourceTemeprature=2500.0;

      double GetTotalProductionRate(int spec,void *SphereDataPointer);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode) {
        unsigned int idim;
        double r=0.0,vbulk[3]={0.0,0.0,0.0},ExternalNormal[3];


        //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_EUROPA'
        double x_LOCAL_IAU_EUROPA[3],x_LOCAL_GALL_EPHIOD_EUROPA[3],v_LOCAL_IAU_EUROPA[3],v_LOCAL_GALL_EPHIOD_EUROPA[3];
        SpiceDouble xform[6][6];

        memcpy(xform,OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix,36*sizeof(double));

        //Geenrate new particle position
        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=pow(ExternalNormal[idim],2);
        }

        r=sqrt(r);

        for (idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]/=r;
          x_LOCAL_IAU_EUROPA[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
        }

        //transfer the position into the coordinate frame related to the rotating coordinate frame 'GALL_EPHIOD'
        x_LOCAL_GALL_EPHIOD_EUROPA[0]=xform[0][0]*x_LOCAL_IAU_EUROPA[0]+xform[0][1]*x_LOCAL_IAU_EUROPA[1]+xform[0][2]*x_LOCAL_IAU_EUROPA[2];
        x_LOCAL_GALL_EPHIOD_EUROPA[1]=xform[1][0]*x_LOCAL_IAU_EUROPA[0]+xform[1][1]*x_LOCAL_IAU_EUROPA[1]+xform[1][2]*x_LOCAL_IAU_EUROPA[2];
        x_LOCAL_GALL_EPHIOD_EUROPA[2]=xform[2][0]*x_LOCAL_IAU_EUROPA[0]+xform[2][1]*x_LOCAL_IAU_EUROPA[1]+xform[2][2]*x_LOCAL_IAU_EUROPA[2];


        //determine if the particle belongs to this processor
        startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_GALL_EPHIOD_EUROPA,startNode);
        if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

        //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_EUROPA'
        PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_EUROPA,vbulk,SourceTemeprature,ExternalNormal,_O2_SPEC_);

        //transform the velocity vector to the coordinate frame 'GALL_EPHIOD'
        v_LOCAL_GALL_EPHIOD_EUROPA[0]=xform[3][0]*x_LOCAL_IAU_EUROPA[0]+xform[3][1]*x_LOCAL_IAU_EUROPA[1]+xform[3][2]*x_LOCAL_IAU_EUROPA[2]+
            xform[3][3]*v_LOCAL_IAU_EUROPA[0]+xform[3][4]*v_LOCAL_IAU_EUROPA[1]+xform[3][5]*v_LOCAL_IAU_EUROPA[2];

        v_LOCAL_GALL_EPHIOD_EUROPA[1]=xform[4][0]*x_LOCAL_IAU_EUROPA[0]+xform[4][1]*x_LOCAL_IAU_EUROPA[1]+xform[4][2]*x_LOCAL_IAU_EUROPA[2]+
            xform[4][3]*v_LOCAL_IAU_EUROPA[0]+xform[4][4]*v_LOCAL_IAU_EUROPA[1]+xform[4][5]*v_LOCAL_IAU_EUROPA[2];

        v_LOCAL_GALL_EPHIOD_EUROPA[2]=xform[5][0]*x_LOCAL_IAU_EUROPA[0]+xform[5][1]*x_LOCAL_IAU_EUROPA[1]+xform[5][2]*x_LOCAL_IAU_EUROPA[2]+
            xform[5][3]*v_LOCAL_IAU_EUROPA[0]+xform[5][4]*v_LOCAL_IAU_EUROPA[1]+xform[5][5]*v_LOCAL_IAU_EUROPA[2];

        memcpy(x_GALL_EPHIOD_EUROPA,x_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
        memcpy(x_IAU_EUROPA,x_LOCAL_IAU_EUROPA,3*sizeof(double));
        memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));

        return true;
      }

    }

    namespace SolarWindSputtering {
      const double Yield=0.1;

      const double minInjectionEnergy=pow(10.0,2)*_O2__MASS_/2.0;
      const double maxInjectionEnergy=pow(10.0E3,2)*_O2__MASS_/2.0;

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetLocalProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double norm_IAU_EUROPA[3],norm_GALL_EPHIOD_EUROPA[3];
        int idim;

        if (((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement]<=0.0) return 0.0;

        //get the normal to the surface element in 'IAU' frame
        memcpy(norm_IAU_EUROPA,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));
        for (idim=0;idim<DIM;idim++) norm_IAU_EUROPA[idim]*=-1.0;

        //convert the normal vector to the 'GALL_EPHIOD' frame
        norm_GALL_EPHIOD_EUROPA[0]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][2]*norm_IAU_EUROPA[2]);

        norm_GALL_EPHIOD_EUROPA[1]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][2]*norm_IAU_EUROPA[2]);

        norm_GALL_EPHIOD_EUROPA[2]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][2]*norm_IAU_EUROPA[2]);

        //get the surface element that is pointer by the vectorm norm_GALL_EPHIOD_EUROPA
        long int nZenithElement,nAzimuthalElement,nd;

        ((cInternalSphericalData*)SphereDataPointer)->GetSurfaceElementProjectionIndex(norm_GALL_EPHIOD_EUROPA,nZenithElement,nAzimuthalElement);
        nd=((cInternalSphericalData*)SphereDataPointer)->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

        //return the local source rate
        return ((cInternalSphericalData*)SphereDataPointer)->SolarWindSurfaceFlux[nd]*Yield*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement];
      }

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {
        return Europa::SourceProcesses::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }
    }


    double totalProductionRate(int spec,void *SphereDataPointer);
    long int InjectionBoundaryModel(void *SphereDataPointer);
  }
*/

  namespace EuropaO2Neutrals {
      const double Yield=0.1;

      const double minInjectionEnergy=pow(1E1,2)*_MASS_(_O2_)/2.0;
      const double maxInjectionEnergy=pow(1E5,2)*_MASS_(_O2_)/2.0;

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection positino on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate nemerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetLocalProductionRate(int spec,int SurfaceElement,void *SphereDataPointer) {
        double norm_IAU_EUROPA[3],norm_GALL_EPHIOD_EUROPA[3];
        int idim;

        if (((cInternalSphericalData*)SphereDataPointer)->SurfaceElementPopulation[_NA_SPEC_][SurfaceElement]<=0.0) return 0.0;

        //get the normal to the surface element in 'IAU' frame
        memcpy(norm_IAU_EUROPA,(((cInternalSphericalData*)SphereDataPointer)->SurfaceElementExternalNormal+SurfaceElement)->norm,3*sizeof(double));
        for (idim=0;idim<DIM;idim++) norm_IAU_EUROPA[idim]*=-1.0;

        //convert the normal vector to the 'GALL_EPHIOD' frame
        norm_GALL_EPHIOD_EUROPA[0]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][2]*norm_IAU_EUROPA[2]);

        norm_GALL_EPHIOD_EUROPA[1]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][2]*norm_IAU_EUROPA[2]);

        norm_GALL_EPHIOD_EUROPA[2]=
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][0]*norm_IAU_EUROPA[0])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][1]*norm_IAU_EUROPA[1])+
            (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][2]*norm_IAU_EUROPA[2]);

        //get the surface element that is pointer by the vectorm norm_GALL_EPHIOD_EUROPA
        long int nZenithElement,nAzimuthalElement,nd;

        ((cInternalSphericalData*)SphereDataPointer)->GetSurfaceElementProjectionIndex(norm_GALL_EPHIOD_EUROPA,nZenithElement,nAzimuthalElement);
        nd=((cInternalSphericalData*)SphereDataPointer)->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);

        //return the local source rate
        return ((cInternalSphericalData*)SphereDataPointer)->SolarWindSurfaceFlux[nd]*Yield*((cInternalSphericalData*)SphereDataPointer)->SurfaceElementArea[SurfaceElement];
      }

      inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate;}

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {

        exit(__LINE__,__FILE__,"not implemented");
        return false;

        //        return Europa::SourceProcesses::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }


    inline int O2SputterInjection(double& Vel, double& RelWeight){
      static double vmin = sqrt(2.0*minInjectionEnergy/_MASS_(_O2_));
      static double vmax = sqrt(2.0*maxInjectionEnergy/_MASS_(_O2_));
      static double F, U, P, v;
      static int initialize = 0;

      if (initialize == 0) {
        U = 0.015*ElectronCharge;
        F = log(vmax)-log(vmin); // Integral of used injection dF ~ 1/v
        // Integral of the injection flux to be modeled
        v = vmax;
        P=U*(-2*v*sqrt(U*_MASS_(_O2_))+sqrt(2)*atan(_MASS_(_O2_)*v*sqrt(2)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1)/2)*_MASS_(_O2_)*v*v+
        		2*U*sqrt(2)*atan(_MASS_(_O2_)*v*sqrt(2)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1)/2))/(_MASS_(_O2_)*v*v+2*U)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1);
        v= vmin;
        P-=U*(-2*v*sqrt(U*_MASS_(_O2_))+sqrt(2)*atan(_MASS_(_O2_)*v*sqrt(2)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1)/2)*_MASS_(_O2_)*v*v+2*U*sqrt(2)*atan(_MASS_(_O2_)*v*sqrt(2)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1)/2))/(_MASS_(_O2_)*v*v+2*U)*pow(U*_MASS_(_O2_),-0.1e1/0.2e1);

        //printf("P = %e in modeled interval vmin=%.2e to vmax=%.2e\n", P, vmin, vmax);
        //printf("%1.4f times flux in interval vmin=0 to vmax=infinity\n", P/(3.14159265359*sqrt(U/(2.0*m))));

        initialize = 1;
      }

      Vel = rnd() * (log(vmax)-log(vmin)) + log(vmin);
      Vel = exp(Vel);
      RelWeight = U * _MASS_(_O2_) * Vel * Vel / pow(_MASS_(_O2_) * Vel * Vel / 2.0 + U, 2.0); // Flux density distribution to be modeled
      RelWeight *= F / P * Vel; // RelWeight = F/dF*dP/P

      return 0;
    }


    inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere,cSingleVariableDiscreteDistribution<int> *SurfaceInjectionDistribution,cSingleVariableDistribution<int> *EnergyDistribution) {
       double ExternalNormal[3];
       double RelWeight;

       //'x' is the position of a particle in the coordinate frame related to the planet 'IAU_EUROPA'
       double x_LOCAL_IAU_EUROPA[3],x_LOCAL_GALL_EPHIOD_EUROPA[3],v_LOCAL_IAU_EUROPA[3],v_LOCAL_GALL_EPHIOD_EUROPA[3];
       int nZenithElement,nAzimuthalElement;
       unsigned int idim,el;

       el=SurfaceInjectionDistribution->DistributeVariable();
       Europa::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,el);
       Europa::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);

       x_LOCAL_IAU_EUROPA[0]=sphereX0[0]+sphereRadius*ExternalNormal[0];
       x_LOCAL_IAU_EUROPA[1]=sphereX0[1]+sphereRadius*ExternalNormal[1];
       x_LOCAL_IAU_EUROPA[2]=sphereX0[2]+sphereRadius*ExternalNormal[2];

       ExternalNormal[0]*=-1.0;
       ExternalNormal[1]*=-1.0;
       ExternalNormal[2]*=-1.0;

       //transfer the position into the coordinate frame related to the rotating coordinate frame 'GALL_EPHIOD'
       x_LOCAL_GALL_EPHIOD_EUROPA[0]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[0][2]*x_LOCAL_IAU_EUROPA[2]);

       x_LOCAL_GALL_EPHIOD_EUROPA[1]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[1][2]*x_LOCAL_IAU_EUROPA[2]);

       x_LOCAL_GALL_EPHIOD_EUROPA[2]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[2][2]*x_LOCAL_IAU_EUROPA[2]);


       //determine if the particle belongs to this processor
       startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_GALL_EPHIOD_EUROPA,startNode);
       if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

       //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_EUROPA'
       double c=0.0,rVel=0.0,lVel[3];
       double Speed;

       O2SputterInjection(Speed, RelWeight);

// check???
       for (idim=0;idim<3;idim++) {
         lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
         rVel+=pow(lVel[idim],2);

         c+=ExternalNormal[idim]*lVel[idim];
       }

       rVel=Speed/sqrt(rVel);

       if (c>0.0) {
         //the distributed velocity vector is directed into the domain
         for (idim=0;idim<3;idim++) v_LOCAL_IAU_EUROPA[idim]=lVel[idim]*rVel;
       }
       else {
         //the distributed velocity vector is directed into the planet -> redirect it
         for (idim=0;idim<3;idim++) v_LOCAL_IAU_EUROPA[idim]=(lVel[idim]-2.0*c*ExternalNormal[idim])*rVel;
       }

       //transform the velocity vector to the coordinate frame 'GALL_EPHIOD'
       v_LOCAL_GALL_EPHIOD_EUROPA[0]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][2]*x_LOCAL_IAU_EUROPA[2])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][3]*v_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][4]*v_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[3][5]*v_LOCAL_IAU_EUROPA[2]);

       v_LOCAL_GALL_EPHIOD_EUROPA[1]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][2]*x_LOCAL_IAU_EUROPA[2])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][3]*v_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][4]*v_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[4][5]*v_LOCAL_IAU_EUROPA[2]);

       v_LOCAL_GALL_EPHIOD_EUROPA[2]=
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][0]*x_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][1]*x_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][2]*x_LOCAL_IAU_EUROPA[2])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][3]*v_LOCAL_IAU_EUROPA[0])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][4]*v_LOCAL_IAU_EUROPA[1])+
           (OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[5][5]*v_LOCAL_IAU_EUROPA[2]);

       memcpy(x_GALL_EPHIOD_EUROPA,x_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));
       memcpy(x_IAU_EUROPA,x_LOCAL_IAU_EUROPA,3*sizeof(double));
       memcpy(v_GALL_EPHIOD_EUROPA,v_LOCAL_GALL_EPHIOD_EUROPA,3*sizeof(double));

       return true;
     }


    inline int GetIAU_to_LS_TransformationMatrix(){
    	// evaluates rotation matrix from IAU system to local solar coordinate system
    	// x-axis pointing towards sun, y perpendicular to x and spin axis, z completes right-handed system
    	double xLocalSun[6], yLocalSun[3],zLocalSun[3], length;
    	SpiceDouble lt;
    	spkezr_c("SUN",Europa::OrbitalMotion::et,"IAU_Europa","none","Europa",xLocalSun,&lt);
    	length = 0.0;
    	for(int i=0;i<3;i++) length+=pow(xLocalSun[i],2);
    	length = sqrt(length);
    	for(int i=0;i<3;i++) xLocalSun[i]/=length;
    	// y is perpendicular to plane of xLocalSun*RotationAxis ([0,0,1])
    	yLocalSun[0] = xLocalSun[1];
    	yLocalSun[1] = -xLocalSun[0];
    	yLocalSun[2] = 0.0;
    	length = 0.0;
    	for(int i=0;i<3;i++) length+=pow(yLocalSun[i],2);
    	length = sqrt(length);
    	for(int i=0;i<3;i++) yLocalSun[i]/=length;
    	// z is perpendicular to x-y plane
    	zLocalSun[0] = -xLocalSun[2]*yLocalSun[1];
    	zLocalSun[1] = xLocalSun[2]*yLocalSun[0];
    	zLocalSun[2] = xLocalSun[0]*yLocalSun[1]-xLocalSun[1]*yLocalSun[0];
    	for(int i=0;i<3;i++) Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[0][i]=xLocalSun[i];
    	for(int i=0;i<3;i++) Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[1][i]=yLocalSun[i];
    	for(int i=0;i<3;i++) Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[2][i]=zLocalSun[i];
    	return 0;
    }

    inline int GetLatLon_vs_Sun(double IAUCoord[3], double &LatVsSun, double &LonVsSun) {
    	double xLocal[3], length;
    	for(int i=0;i<3;i++) xLocal[i] = 0.0;
    	for(int j=0;j<3;j++) for(int i=0;i<3;i++) xLocal[j]+=Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[j][i]*IAUCoord[i];
    	for(int i=0;i<3;i++) length+=pow(xLocal[i],2);
    	length = sqrt(length);
    	for(int i=0;i<3;i++) xLocal[i]/=length;
    	LatVsSun = asin(xLocal[2]);
    	if (xLocal[1] != 0.0)
    		LonVsSun = xLocal[1]/fabs(xLocal[1])*acos(xLocal[0]/cos(LatVsSun));
    	else if (xLocal[0] > 0.0)
    		LonVsSun = 0.0;
    	else
    		LonVsSun = Pi;

    	if (LonVsSun < 0.0) LonVsSun += 2.0*Pi;

    	return 0;
    }

    inline double GetSurfaceTemperature(double IAUCoord[3]) {
    	 double subSolarAngle, Temp, LatVsSun, LonVsSun;

    	 GetLatLon_vs_Sun(IAUCoord, LatVsSun, LonVsSun);
    	 subSolarAngle=acos(cos(LonVsSun)*cos(LatVsSun));

    	 Temp = 40.0+50.0*pow(fabs(cos(LatVsSun)),0.25);
    	 if (subSolarAngle<Pi/2.0) Temp += 40.0*cos(subSolarAngle);

    	 return Temp;
    }

    inline double AdsorptionTimeO2(double IAUCoord[3]) {
    	// Surface sticking time for O2 in seconds
    	return 3.80738E16*pow(GetSurfaceTemperature(IAUCoord),-8.80346);
    }

  double totalProductionRate(int spec,void *SphereDataPointer);
  long int InjectionBoundaryModel(void *SphereDataPointer);
}

  namespace InjectEuropaMagnetosphericEPDIons {
      const double Yield=0.1;

      const double minInjectionEnergy=1E1*1000*ElectronCharge; // Energy range minimum 10 keV
      const double maxInjectionEnergy=1E5*1000*ElectronCharge; // Energy range minimum 100'000 keV

      extern double SourceRate,maxLocalSourceRate;

      //the object for distribution of injection position on the planet's surface
      extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution;
      double GetSurfaceElementSodiumProductionRate(int nElement);

      //evaluate numerically the source rate
      extern double CalculatedTotalSodiumSourceRate;

      inline double GetTotalProductionRate(int spec) {
    	  //double TotalSurfaceArea;
    	  const static double EPD_Flux = 1.243919681E6*2*Pi*1e4; // Total Flux [(m^2*s)^-1]

    	  return EPD_Flux;
      }

      //energy distribution function of injected particles
      extern cSingleVariableDistribution<int> EnergyDistribution;
      double EnergyDistributionFunction(double e);

      //generate particle properties
      inline bool GenerateParticleProperties(double *x_GALL_EPHIOD_EUROPA,double *x_IAU_EUROPA,double *v_GALL_EPHIOD_EUROPA,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere) {

        exit(__LINE__,__FILE__,"not implemented");
        return false;

        //        return Europa::SourceProcesses::GenerateParticleProperties(x_GALL_EPHIOD_EUROPA,x_IAU_EUROPA,v_GALL_EPHIOD_EUROPA,sphereX0,sphereRadius,startNode,Sphere,&SurfaceInjectionDistribution,&EnergyDistribution);
      }


    inline int MagnetosphericInjection(double& Vel, double& RelWeight){
      //const static double EPD_Flux = 1.243919681E6*1E4; // Total Flux [(m^2*sr*s)^-1]
      const static double vmin = sqrt(2*minInjectionEnergy/_MASS_(_S_)); // Energy range minimum 10 keV
      const static double vmax = sqrt(2*maxInjectionEnergy/_MASS_(_S_)); // Energy range minimum 100'000 keV

      Vel = rnd() * (log(vmax)-log(vmin)) + log(vmin); // Create random speed with ~1/v probability
      Vel = exp(Vel);

      // Relative weight from ratio of injected distribution (~1/v) and measured distribution from E12 flyby (Paranicas 2002), calculated by Maple
      RelWeight=0.7452471071e-13*pow(Vel,4)*pow(0.1658282218e-9*Vel*Vel+0.5442256265e3,-0.374873e1)/(1+0.6619630708e-44*pow(Vel*Vel,0.31176e1));

      return 0;
    }


    bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    double BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);


    //injection of model particles through the faces of the bounding box
    inline long int  BoundingBoxInjection(int Spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
      bool ExternalFaces[6];
      double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
      int nface,idim;
    //  long int nInjectedParticles;
      long int newParticle;
      PIC::ParticleBuffer::byte *newParticleData;
      long int nInjectedParticles=0;

      if (Spec!=_OPLUS_HIGH_SPEC_) return 0; //inject only O ions







      //static double swVel[3]={4.0E5,000.0,000.0}, nSW=5.0E6;
      //static const double tempSW=8.0E4;
      double v[3],lVel[3],rVel,Speed, RelWeight,c;


      double dTimeCounter;

      double ModelParticlesInjectionRate;

      if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
        for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
          startNode->GetExternalNormal(ExternalNormal,nface);



//if ((nface!=4)&&(nface!=5)) continue;


          //inject thermal O ions
          ParticleWeight=startNode->block->GetLocalParticleWeight(_OPLUS_THERMAL_SPEC_);
          LocalTimeStep=startNode->block->GetLocalTimeStep(_OPLUS_THERMAL_SPEC_);
          TimeCounter=0.0;
          ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(Thermal_OPlus_NumberDensity,Thermal_OPlus_Temperature,Thermal_OPlus_BulkVelocity,ExternalNormal,_OPLUS_THERMAL_SPEC_);

          if (ModelParticlesInjectionRate>0.0) {
            ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;
            PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

            while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
              //generate the new particle position on the face
              for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim]-ExternalNormal[idim]*PIC::Mesh::mesh.EPS;

              //generate particles' velocity
              PIC::Distribution::InjectMaxwellianDistribution(v,Thermal_OPlus_BulkVelocity,Thermal_OPlus_Temperature,ExternalNormal,_OPLUS_THERMAL_SPEC_,-1);

              //generate new particle and inject it into the system
              newParticle=PIC::ParticleBuffer::GetNewParticle();
              newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
              nInjectedParticles++;

              PIC::ParticleBuffer::SetX(x,newParticleData);
              PIC::ParticleBuffer::SetV(v,newParticleData);
              PIC::ParticleBuffer::SetI(_OPLUS_THERMAL_SPEC_,newParticleData);
              PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);
              Europa::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_,newParticleData);

              _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
            }

          }

          //inject the high energy O ions
          ParticleWeight=startNode->block->GetLocalParticleWeight(_OPLUS_HIGH_SPEC_);
          LocalTimeStep=startNode->block->GetLocalTimeStep(_OPLUS_HIGH_SPEC_);
          TimeCounter=0.0;
          ModelParticlesInjectionRate=GetTotalProductionRate(_OPLUS_HIGH_SPEC_);

          if (ModelParticlesInjectionRate>0.0) {

            ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;


            PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);
            while (true) {    /////((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
              //generate the new particle position on the face
              for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim]-ExternalNormal[idim]*PIC::Mesh::mesh.EPS;


              MagnetosphericInjection(Speed, RelWeight);


//==========  DEBUG  ============

//Speed=1.0E6;
//RelWeight=1.0;
 //==========  END DEBUG  ============



              dTimeCounter=-log(rnd())/ModelParticlesInjectionRate*RelWeight;

              if (TimeCounter+dTimeCounter>LocalTimeStep) break;
              TimeCounter+=dTimeCounter;


              //sample the injection rate
              Europa::TotalInjectionRate+=ParticleWeight*RelWeight;

              //generate a particle
              newParticle=PIC::ParticleBuffer::GetNewParticle();
              newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
              nInjectedParticles++;


              //         PIC::Distribution::MaxwellianVelocityDistribution(v,swVel,tempSW,SW);

          //    MagnetosphericInjection(Speed, RelWeight);

              rVel=0.0,c=0.0;




              for (idim=0;idim<3;idim++) {
                 lVel[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
                 rVel+=pow(lVel[idim],2);

                 c+=ExternalNormal[idim]*lVel[idim];
               }
               rVel=Speed/sqrt(rVel);

               if (c<0.0) {
                 //the distributed velocity vector is directed into the domain
                 for (idim=0;idim<3;idim++) v[idim]=lVel[idim]*rVel;
               }
               else {
                 //the distributed velocity vector is directed into the planet -> redirect it
                 for (idim=0;idim<3;idim++) v[idim]=(lVel[idim]-2.0*c*ExternalNormal[idim])*rVel;
               }



               //cosine distribution of the injection angle
               double theta,phi,sinTheta,cosTheta,sinPhi,cosPhi,le0,le1;

               phi=2.0*Pi*rnd();
               theta=asin(rnd());

//               theta=0.0;

//               RelWeight=cos(theta);

/*               double p;

               do {
                 theta=Pi/2.0*rnd();

                 p=sin(theta)*cos(theta);
               }
               while (p<rnd());*/


               theta=asin(sqrt(rnd()));


               sinPhi=sin(phi);
               cosPhi=cos(phi);

               sinTheta=sin(theta);
               cosTheta=cos(theta);

               le0=sqrt(e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2]);
               le1=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);

               for (idim=0;idim<3;idim++) v[idim]=Speed*(-cosTheta*ExternalNormal[idim]+sinTheta*(cosPhi*e0[idim]/le0+sinPhi*e1[idim]/le1));

//============= DEBUG ==============
/*
if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>5.0E7*5.0E7) {
  cout << "ERROR: Velocity exceeds the limit, " << __FILE__ << "%" << __LINE__ << endl;
}

v[0]=Speed,v[1]=0.0,v[2]=0.0;
*/
//============= END DEBUG ==============


              PIC::ParticleBuffer::SetX(x,newParticleData);
              PIC::ParticleBuffer::SetV(v,newParticleData);
              PIC::ParticleBuffer::SetI(_OPLUS_HIGH_SPEC_,newParticleData);

              PIC::ParticleBuffer::SetIndividualStatWeightCorrection(RelWeight,newParticleData);
              Europa::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_,newParticleData);


              //inject the particle into the system
              _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
            }
          }


        }
      }

      return nInjectedParticles;
    }

    long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);


   //Europa addition
    inline double SputteringYield(double v, double m, int spec){
      // Matrix contains sputtering yields for
      // Column 1 H^+, use spec = 0
      // Column 2 O^(n+), use spec = 1
      // Column 3 O2^(n+) , use spec = 2
      // Entry is calculated from E(eV) by int((log10(E)-1.47)/0.1)
      // Sputtering yield is interpolated by log10(Y)=b*log10(E)+a
      const static double YieldData[77][3]={
        {0.0000E+00, 0.0000E+00, 0.0000E+00},
        {1.2189E+00, 0.0000E+00, 0.0000E+00},
        {1.6346E+00, 0.0000E+00, 0.0000E+00},
        {2.0413E+00, 1.3084E+00, 0.0000E+00},
        {2.4445E+00, 1.7076E+00, 0.0000E+00},
        {2.9274E+00, 2.2285E+00, 0.0000E+00},
        {3.5057E+00, 2.9084E+00, 1.3049E+00},
        {4.1424E+00, 3.7957E+00, 1.7029E+00},
        {4.4843E+00, 4.5900E+00, 2.2224E+00},
        {4.8545E+00, 4.9956E+00, 2.9004E+00},
        {5.2552E+00, 5.4370E+00, 3.7853E+00},
        {5.5834E+00, 5.9174E+00, 4.6910E+00},
        {5.6050E+00, 6.1035E+00, 5.3956E+00},
        {5.6267E+00, 6.2297E+00, 6.2061E+00},
        {5.6486E+00, 6.3585E+00, 7.1383E+00},
        {5.6705E+00, 6.4899E+00, 7.7580E+00},
        {5.6924E+00, 6.8240E+00, 8.3513E+00},
        {5.7145E+00, 7.2673E+00, 8.9992E+00},
        {5.7367E+00, 7.7394E+00, 9.6974E+00},
        {5.7589E+00, 8.2421E+00, 1.0388E+01},
        {6.2521E+00, 8.8610E+00, 1.1014E+01},
        {6.8697E+00, 9.7362E+00, 1.1678E+01},
        {7.5482E+00, 1.0698E+01, 1.2255E+01},
        {8.4838E+00, 1.1755E+01, 1.2742E+01},
        {9.6829E+00, 1.2916E+01, 1.3331E+01},
        {1.1051E+01, 1.4191E+01, 1.4073E+01},
        {1.2626E+01, 1.5721E+01, 1.4857E+01},
        {1.4436E+01, 1.8377E+01, 1.8273E+01},
        {1.6507E+01, 2.1482E+01, 2.3093E+01},
        {1.8874E+01, 2.5111E+01, 2.9116E+01},
        {2.1581E+01, 3.1285E+01, 3.6559E+01},
        {2.4676E+01, 4.0769E+01, 4.5904E+01},
        {2.7509E+01, 4.7429E+01, 5.8892E+01},
        {3.0155E+01, 5.1168E+01, 7.5930E+01},
        {3.0883E+01, 5.5574E+01, 9.7896E+01},
        {2.9704E+01, 6.9086E+01, 1.2629E+02},
        {2.6026E+01, 8.5882E+01, 1.6294E+02},
        {2.1287E+01, 1.0676E+02, 2.1024E+02},
        {1.5623E+01, 1.3320E+02, 2.7022E+02},
        {1.1466E+01, 1.6724E+02, 3.4406E+02},
        {8.4154E+00, 2.0999E+02, 4.3807E+02},
        {6.3132E+00, 2.6936E+02, 5.6099E+02},
        {4.9850E+00, 3.4728E+02, 7.2768E+02},
        {3.9363E+00, 4.4775E+02, 9.4389E+02},
        {3.1082E+00, 5.4918E+02, 1.2243E+03},
        {2.4543E+00, 6.6287E+02, 1.5532E+03},
        {0.0000E+00, 8.0008E+02, 1.9659E+03},
        {0.0000E+00, 9.6570E+02, 2.4881E+03},
        {0.0000E+00, 1.0561E+03, 3.1346E+03},
        {0.0000E+00, 1.1225E+03, 3.8439E+03},
        {0.0000E+00, 1.0649E+03, 4.7136E+03},
        {0.0000E+00, 1.0102E+03, 5.7683E+03},
        {0.0000E+00, 9.0521E+02, 6.5421E+03},
        {0.0000E+00, 7.4138E+02, 7.4197E+03},
        {0.0000E+00, 6.0720E+02, 7.5824E+03},
        {0.0000E+00, 4.8969E+02, 7.5824E+03},
        {0.0000E+00, 3.8318E+02, 7.0401E+03},
        {0.0000E+00, 2.9983E+02, 6.3214E+03},
        {0.0000E+00, 2.3245E+02, 5.1415E+03},
        {0.0000E+00, 1.7994E+02, 4.1819E+03},
        {0.0000E+00, 1.3895E+02, 3.3329E+03},
        {0.0000E+00, 1.0712E+02, 2.6428E+03},
        {0.0000E+00, 8.2586E+01, 2.0899E+03},
        {0.0000E+00, 6.3903E+01, 1.6317E+03},
        {0.0000E+00, 4.9725E+01, 1.2740E+03},
        {0.0000E+00, 3.8693E+01, 9.9500E+02},
        {0.0000E+00, 3.0108E+01, 7.7766E+02},
        {0.0000E+00, 2.3363E+01, 6.0779E+02},
        {0.0000E+00, 1.8062E+01, 4.6551E+02},
        {0.0000E+00, 1.3963E+01, 3.5118E+02},
        {0.0000E+00, 1.0799E+01, 2.6510E+02},
        {0.0000E+00, 8.3838E+00, 2.0681E+02},
        {0.0000E+00, 6.5090E+00, 1.6133E+02},
        {0.0000E+00, 5.0624E+00, 0.0000E+00},
        {0.0000E+00, 3.9453E+00, 0.0000E+00},
        {0.0000E+00, 3.0748E+00, 0.0000E+00},
        {0.0000E+00, 0.0000E+00, 0.0000E+00},
      };
      double E = 0.5*m*pow(v,2.0)/ElectronCharge;  // Energy in eV
      int counter = int((log10(E)-1.47)/0.1);      // Derive array entry from energy
      if ((E > 1E9) || (YieldData[counter][spec] < 1.0) || (counter < 0)) return 0.0;
      //double Emin = pow(10,1.47+counter*0.1);      // Upper limit in energy interval
      //double Emax = pow(10,1.47+(counter+1)*0.1);  // Lower limit in energy interval
      double b = (log10(YieldData[counter+1][spec])-log10(YieldData[counter][spec]))/0.1;
      double a = log10(YieldData[counter][spec])-b*(1.47+counter*0.1);
      //double b = (log10(YieldData[counter+1][spec])-log10(YieldData[counter][spec]))/(log10(Emax)-log10(Emin));
      //double a = log10(YieldData[counter][spec])-b*log10(Emin);
      //printf("counter = %d\n", counter);
      //printf("Emin = %e eV, E = %e eV, Emax = %e eV\n",Emin,E,Emax);
      //printf("Yield(Emin) = %f, Yield(Emax) = %f\n",YieldData[counter][spec],YieldData[counter+1][2]);
      //printf("a = %.2e,b = %.2e\n",a,b);
      return pow(10,b*log10(E)+a);
    }



  double totalProductionRate(int spec,void *SphereDataPointer);
  long int InjectionBoundaryModel(void *SphereDataPointer);
}


  //combine together the surface area densities and recalculate the
  void ExchangeSurfaceAreaDensity();


  //Interaction of the particles with the surface
  namespace SurfaceInteraction {


    //sampling the surface area density of the sticking species
#define _EUROPA_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_             3

#define _EUROPA_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_DOWN_           0
#define _EUROPA_SURFACE_STICKING_SAMPLING_OFFSET__FLUX_UP_             1
#define _EUROPA_SURFACE_STICKING_SAMPLING_OFFSET__AREA_NUMBER_DENSITY_ 2


    //sodium/surface interaction model
    const double AccomodationCoefficient=0.1;

    //sticking probability of sodium atoms
    inline double SodiumStickingProbability(double Temp) {
      if (Temp<300.0) return 1.0;
      if (Temp<650.0) return 1.0-(Temp-300.0)/350.0;

      return 0.0;
    }


    //model of the interaction between particles and the planetary surface
    int ParticleSphereInteraction_SurfaceAccomodation(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);
  }








  //the total acceleration acting on a particle
//double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);


//defined the forces that acts upon a particle on
#define _FORCE_GRAVITY_MODE_ _PIC_MODE_OFF_
#define _FORCE_LORENTZ_MODE_ _PIC_MODE_OFF_
#define _FORCE_FRAMEROTATION_MODE_ _PIC_MODE_OFF_

void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};



  //Mass and heliocentric velocity of Europa

//  static const double EuropaHeliocentricVelocity=-10.0E3;
//  static const double EuropaHeliocentricDistance=57909100.0E3;   // 0.387 098 AU


  memcpy(x_LOCAL,x,3*sizeof(double));
  memcpy(v_LOCAL,v,3*sizeof(double));

  accl[0]=0.0,accl[1]=0.0,accl[2]=0.0;
//  return;




  //get the radiation pressure acceleration
/*
  if (spec==_O2_SPEC_) {
    double r,r2,vHeliocentric,radiationPressureAcceleration;

    r2=x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];

    if ((r2>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x_LOCAL[0]>0.0)) { //calculate the radiation pressure force
      r2+=pow(xEuropaRadial+x_LOCAL[0],2);
      r=sqrt(r2);

      vHeliocentric=((vEuropaRadial+v_LOCAL[0])*(xEuropaRadial+x_LOCAL[0])+v_LOCAL[1]*x_LOCAL[1]+v_LOCAL[2]*x_LOCAL[2])/r;
      radiationPressureAcceleration=SodiumRadiationPressureAcceleration__Combi_1997_icarus(vHeliocentric,r);

      accl_LOCAL[0]-=radiationPressureAcceleration*(xEuropaRadial+x_LOCAL[0])/r;
      accl_LOCAL[1]-=radiationPressureAcceleration*x_LOCAL[1]/r,accl_LOCAL[2]-=radiationPressureAcceleration*x_LOCAL[2]/r;
    }
  }
  else if (spec==_O_PLUS_SPEC_) { //the Lorent force
*/


#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
    long int nd;
    char *offset;
    int i,j,k;
    PIC::Mesh::cDataCenterNode *CenterNode;
    double E[3],B[3];

    if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
      startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL,startNode);

      if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
#endif

    CenterNode=startNode->block->GetCenterNode(nd);
    offset=CenterNode->GetAssociatedDataBufferPointer();

/*    if (*((int*)(offset+PIC::CPLR::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
      memcpy(E,offset+PIC::CPLR::ICES::ElectricFieldOffset,3*sizeof(double));
      memcpy(B,offset+PIC::CPLR::ICES::MagneticFieldOffset,3*sizeof(double));
    }
    else {
      memcpy(E,swE_Typical,3*sizeof(double));
      memcpy(B,swB_Typical,3*sizeof(double));
    }*/

    PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
    PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);

    double ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    double mass=PIC::MolecularData::GetMass(spec);

    accl_LOCAL[0]+=ElectricCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
    accl_LOCAL[1]+=ElectricCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
    accl_LOCAL[2]+=ElectricCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;

//    exit(__LINE__,__FILE__,"error: not implemented");
#endif

//  }


#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
  //the gravity force
  double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
  }
#endif


/*  //correct the gravity acceleration: accout for solar gravity of the particle location
  double rSolarVector[3],r2Solar,rSolar;

  rSolarVector[0]=xEuropaRadial-x_LOCAL[0];
  rSolarVector[1]=x_LOCAL[1],rSolarVector[2]=x_LOCAL[2];

  r2Solar=(rSolarVector[0]*rSolarVector[0])+(rSolarVector[1]*rSolarVector[1])+(rSolarVector[2]*rSolarVector[2]);
  rSolar=sqrt(r2Solar);

  accl_LOCAL[0]+=GravityConstant*_MASS_(_SUN_)*(1.0/r2Solar*rSolarVector[0]/rSolar-1.0/pow(xEuropaRadial,2)); //x-axis in solar coordinate frame and GALL_EPHIOD have opposite direction
  accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)/r2Solar*rSolarVector[1]/rSolar;
  accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)/r2Solar*rSolarVector[2]/rSolar;*/


  //account for the planetary rotation around the Sun
#if _FORCE_FRAMEROTATION_MODE_ == _PIC_MODE_ON_
  double aCen[3],aCorr[3],t3,t7,t12;

  t3 = RotationVector_SO_FROZEN[0] * x_LOCAL[1] - RotationVector_SO_FROZEN[1] * x_LOCAL[0];
  t7 = RotationVector_SO_FROZEN[2] * x_LOCAL[0] - RotationVector_SO_FROZEN[0] * x_LOCAL[2];
  t12 = RotationVector_SO_FROZEN[1] * x_LOCAL[2] - RotationVector_SO_FROZEN[2] * x_LOCAL[1];

  aCen[0] = -RotationVector_SO_FROZEN[1] * t3 + RotationVector_SO_FROZEN[2] * t7;
  aCen[1] = -RotationVector_SO_FROZEN[2] * t12 + RotationVector_SO_FROZEN[0] * t3;
  aCen[2] = -RotationVector_SO_FROZEN[0] * t7 + RotationVector_SO_FROZEN[1] * t12;


  aCorr[0] = -2.0*(RotationVector_SO_FROZEN[1] * v_LOCAL[2] - RotationVector_SO_FROZEN[2] * v_LOCAL[1]);
  aCorr[1] = -2.0*(RotationVector_SO_FROZEN[2] * v_LOCAL[0] - RotationVector_SO_FROZEN[0] * v_LOCAL[2]);
  aCorr[2] = -2.0*(RotationVector_SO_FROZEN[0] * v_LOCAL[1] - RotationVector_SO_FROZEN[1] * v_LOCAL[0]);

  accl_LOCAL[0]+=aCen[0]+aCorr[0];
  accl_LOCAL[1]+=aCen[1]+aCorr[1];
  accl_LOCAL[2]+=aCen[2]+aCorr[2];

#endif



  //copy the local value of the acceleration to the global one
  memcpy(accl,accl_LOCAL,3*sizeof(double));
}

inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {

  long int nd;
  int i,j,k;
  double PlasmaBulkVelocity[3],ElectronDensity;

  nd=PIC::Mesh::mesh.findCenterNodeIndex(x,i,j,k,node);
  PIC::CPLR::GetBackgroundPlasmaVelocity(PlasmaBulkVelocity,x,nd,node);
  ElectronDensity=PIC::CPLR::GetBackgroundPlasmaNumberDensity(x,nd,node);

  //only sodium can be ionized
  if (spec!=_O2_SPEC_) {
    PhotolyticReactionAllowedFlag=false;
    return -1.0;
  }

  PhotolyticReactionAllowedFlag=true;
  return 3.8e5;


  static const double LifeTime=3600.0*5.8/pow(0.4,2);

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  double res,r2=x[1]*x[1]+x[2]*x[2];

  //check if the particle is outside of the Earth and lunar shadows
  if ( ((r2>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x[0]<0.0)) /*&& (Mercury::EarthShadowCheck(x)==false)*/ ) {
    res=LifeTime,PhotolyticReactionAllowedFlag=true;
  }
  else {
    res=-1.0,PhotolyticReactionAllowedFlag=false;
  }

/*    //check if the particle intersect the surface of the Earth
  if (pow(x[0]-Mercury::xEarth_SO[0],2)+pow(x[1]-Mercury::xEarth_SO[1],2)+pow(x[2]-Mercury::xEarth_SO[2],2)<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) {
    res=1.0E-10*LifeTime,PhotolyticReactionAllowedFlag=true;
  }*/

#else
  double res=LifeTime;
  PhotolyticReactionAllowedFlag=true;
#endif

  return res;
}

inline int ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {


  PIC::ParticleBuffer::SetI(_O2PLUS_SPEC_,ParticleData);
  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;

//  return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
}


//xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode
inline int GenericUnimolecularReactionProcessor(double *xInit,double *xFinal,double *vFinal, int &spec, long int ptr,PIC::ParticleBuffer::byte *ParticleData, double TimeInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double ParentSpeciesLifeTime;
  bool PhotolyticReactionAllowedFlag;
  int ResCode=_GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_;

  if (spec!=_O2_SPEC_) return _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_;


  ParentSpeciesLifeTime=ExospherePhotoionizationLifeTime(xFinal,spec,ptr,PhotolyticReactionAllowedFlag,node);

  //determine if the parent particle sould be removed
  double c,p;

  c=exp(-TimeInterval/ParentSpeciesLifeTime);
  p=1.0-c; //the probability for reaction to occur

  if (rnd()<p) {
    //the reaction has occured -> the original partiucle should be removed
    ResCode=_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_;
  }


  //determine if a daugher particle should be generated
  if (_O2PLUS_SPEC_!=-1) {
    double TimeIntervalProduct;
    PIC::Mesh::cDataBlockAMR *block=node->block;
    long int newParticle,nDaugherParticles;

    TimeIntervalProduct=TimeInterval*block->GetLocalTimeStep(_O2PLUS_SPEC_)/block->GetLocalTimeStep(_O2_SPEC_);
    c=(1.0-exp(-TimeIntervalProduct/ParentSpeciesLifeTime))*block->GetLocalParticleWeight(_O2_SPEC_)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData)/block->GetLocalParticleWeight(_O2PLUS_SPEC_); //the number of model O2+ particles generated during time interval 'TimeIntervalProduct'

    nDaugherParticles=(int)c;
    c-=nDaugherParticles;

    if (rnd()<c) nDaugherParticles++;

    //generate new particles and inject them into the system
    while (--nDaugherParticles>0) {
      newParticle=PIC::ParticleBuffer::GetNewParticle();

      double x[3];
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newParticleNode;

      for (int idim=0;idim<3;idim++) x[idim]=xInit[idim]+rnd()*(xFinal[idim]-xInit[idim]);
      newParticleNode=PIC::Mesh::mesh.findTreeNode(x,node);


      PIC::ParticleBuffer::CloneParticle(newParticle,ptr);
      PIC::ParticleBuffer::SetI(_O2PLUS_SPEC_,newParticle);
      PIC::ParticleBuffer::SetV(vFinal,newParticle);
      PIC::ParticleBuffer::SetX(x,newParticle);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticle);

      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,0.0,newParticleNode);
    }
  }

  return ResCode;

}


}



#endif /* EUROPA_H_ */
