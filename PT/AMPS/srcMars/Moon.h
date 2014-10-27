/*
 * Mercury.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

#ifndef MERCURY_H_
#define MERCURY_H_


#include "Exosphere.h"

namespace Moon {
  using namespace Exosphere;



  //init the model
  void Init_AfterParser();

  //output the column integrals in the anti-solar direction
  namespace AntiSolarDirectionColumnMap {
    const int nAzimuthPoints=150;
    const double maxZenithAngle=Pi/4.0;
    const double dZenithAngleMin=0.001*maxZenithAngle,dZenithAngleMax=0.02*maxZenithAngle;


    void Print(int DataOutputFileNumber);
  }

  //sampling procedure uniqueck for the model of the lunar exosphere
  namespace Sampling {
    using namespace Exosphere::Sampling;

    namespace SubsolarLimbColumnIntegrals {
      //sample values of the column integrals at the subsolar limb as a function of the phase angle
      const int nPhaseAngleIntervals=50;
      const double dPhaseAngle=Pi/nPhaseAngleIntervals;

      //sample the altitude variatin of the column integral in the solar direction from the limb
      const double maxAltitude=50.0; //altititude in radii of the object
      const double dAltmin=0.001,dAltmax=0.001*maxAltitude; //the minimum and maximum resolution along the 'altitude' line
      extern double rr;
      extern long int nSampleAltitudeDistributionPoints;

      struct cSampleAltitudeDistrubutionBufferElement {
        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
      };

      extern cSampleAltitudeDistrubutionBufferElement *SampleAltitudeDistrubutionBuffer;


      extern SpiceDouble etSampleBegin;
      extern int SamplingPhase;
      extern int firstPhaseRadialVelocityDirection;

      extern int nOutputFile;

      struct cSampleBufferElement {
        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
        cSampleAltitudeDistrubutionBufferElement *AltitudeDistrubutionBuffer;

        int nSamples;
        double JulianDate;
      };

      //sample buffer for the current sampling cycle
      extern cSampleBufferElement SampleBuffer_AntiSunwardMotion[nPhaseAngleIntervals];
      extern cSampleBufferElement SampleBuffer_SunwardMotion[nPhaseAngleIntervals];

      //sample buffer for the lifetime of the simulation
      extern cSampleBufferElement SampleBuffer_AntiSunwardMotion__TotalModelRun[nPhaseAngleIntervals];
      extern cSampleBufferElement SampleBuffer_SunwardMotion__TotalModelRun[nPhaseAngleIntervals];

      void EmptyFunction();

      void init();
      void CollectSample(int DataOutputFileNumber); //the function will ba called by that part of the core that prints output files
      void PrintDataFile();
    }

  }


  //check if the point is in the Earth shadow: the function returns 'true' if it is in the shadow, and 'false' if the point is outside of the Earth shadow
  bool inline EarthShadowCheck(double *x_LOCAL_SO_OBJECT) {
    double lPerp[3],lSun2=0.0,c=0.0,xSun_LOCAL[3],xEarth_LOCAL[3];
    int idim;

    memcpy(xSun_LOCAL,xSun_SO,3*sizeof(double));
    memcpy(xEarth_LOCAL,xEarth_SO,3*sizeof(double));

    for (idim=0;idim<3;idim++) {
      xSun_LOCAL[idim]-=x_LOCAL_SO_OBJECT[idim];
      xEarth_LOCAL[idim]-=x_LOCAL_SO_OBJECT[idim];

      lSun2+=pow(xSun_LOCAL[idim],2);
      c+=xSun_LOCAL[idim]*xEarth_LOCAL[idim];
    }

    for (idim=0;idim<3;idim++) {
      lPerp[idim]=xEarth_LOCAL[idim]-c*xSun_LOCAL[idim]/lSun2;
    }

    return (lPerp[0]*lPerp[0]+lPerp[1]*lPerp[1]+lPerp[2]*lPerp[2]<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) ? true : false;
  }


    //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

    //get the radiation pressure acceleration
    if (spec==_NA_SPEC_) {
      if ( ((x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x_LOCAL[0]>0.0)) && (EarthShadowCheck(x_LOCAL)==false) ) { //calculate the radiation pressure force
        double rHeliocentric,vHeliocentric,radiationPressureAcceleration;

        //calcualte velocity of the particle in the "Frozen SO Frame";
        double v_LOCAL_SO_FROZEN[3];

        v_LOCAL_SO_FROZEN[0]=Exosphere::vObject_SO_FROZEN[0]+v_LOCAL[0]+
            Exosphere::RotationVector_SO_FROZEN[1]*x_LOCAL[2]-Exosphere::RotationVector_SO_FROZEN[2]*x_LOCAL[1];

        v_LOCAL_SO_FROZEN[1]=Exosphere::vObject_SO_FROZEN[1]+v_LOCAL[1]-
            Exosphere::RotationVector_SO_FROZEN[0]*x_LOCAL[2]+Exosphere::RotationVector_SO_FROZEN[2]*x_LOCAL[0];

        v_LOCAL_SO_FROZEN[2]=Exosphere::vObject_SO_FROZEN[2]+v_LOCAL[2]+
            Exosphere::RotationVector_SO_FROZEN[0]*x_LOCAL[1]-Exosphere::RotationVector_SO_FROZEN[1]*x_LOCAL[0];


        rHeliocentric=sqrt(pow(x_LOCAL[0]-Exosphere::xObjectRadial,2)+(x_LOCAL[1]*x_LOCAL[1])+(x_LOCAL[2]*x_LOCAL[2]));
        vHeliocentric=(
            (v_LOCAL_SO_FROZEN[0]*(x_LOCAL[0]-Exosphere::xObjectRadial))+
            (v_LOCAL_SO_FROZEN[1]*x_LOCAL[1])+(v_LOCAL_SO_FROZEN[2]*x_LOCAL[2]))/rHeliocentric;

        radiationPressureAcceleration=SodiumRadiationPressureAcceleration__Combi_1997_icarus(vHeliocentric,rHeliocentric);

        accl_LOCAL[0]+=radiationPressureAcceleration*(x_LOCAL[0]-Exosphere::xObjectRadial)/rHeliocentric;
        accl_LOCAL[1]+=radiationPressureAcceleration*x_LOCAL[1]/rHeliocentric;
        accl_LOCAL[2]+=radiationPressureAcceleration*x_LOCAL[2]/rHeliocentric;
      }
    }
    else if (spec==_NA_PLUS_SPEC_) { //the Lorentz force
      long int nd;
      char *offset;
      int i,j,k;
      PIC::Mesh::cDataCenterNode *CenterNode;
      double E[3],B[3];

      exit(__LINE__,__FILE__,"check the numbers!");

      if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  #endif

      CenterNode=startNode->block->GetCenterNode(nd);
      offset=CenterNode->GetAssociatedDataBufferPointer();

      if (*((int*)(offset+PIC::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
        memcpy(E,offset+PIC::ICES::ElectricFieldOffset,3*sizeof(double));
        memcpy(B,offset+PIC::ICES::MagneticFieldOffset,3*sizeof(double));
      }
      else {
        memcpy(E,swE_Typical,3*sizeof(double));
        memcpy(B,swB_Typical,3*sizeof(double));
      }


      accl_LOCAL[0]+=ElectronCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/_MASS_(_NA_);
      accl_LOCAL[1]+=ElectronCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/_MASS_(_NA_);
      accl_LOCAL[2]+=ElectronCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/_MASS_(_NA_);

    }


    //the gravity force
    double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
    double r=sqrt(r2);
    int idim;

    for (idim=0;idim<DIM;idim++) {
      accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
    }


    //correct the gravity acceleration: accout for solar gravity of the particle location
    //correction of the solar gravity

    double rSun2Moon,rSun2Particle;

    rSun2Moon=xObjectRadial;
    rSun2Particle=sqrt(pow(x_LOCAL[0]-xObjectRadial,2)+pow(x_LOCAL[1],2)+pow(x_LOCAL[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_SUN_)*((x_LOCAL[0]-xObjectRadial)/pow(rSun2Particle,3)+xObjectRadial/pow(rSun2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[1]/pow(rSun2Particle,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[2]/pow(rSun2Particle,3));

    //correction of the Earth gravity
    double rEarth2Moon,rEarth2Particle;

    rEarth2Moon=sqrt(pow(xEarth_SO[0],2)+pow(xEarth_SO[1],2)+pow(xEarth_SO[2],2));
    rEarth2Particle=sqrt(pow(x_LOCAL[0]-xEarth_SO[0],2)+pow(x_LOCAL[1]-xEarth_SO[1],2)+pow(x_LOCAL[2]-xEarth_SO[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[0]-xEarth_SO[0])/pow(rEarth2Particle,3)+xEarth_SO[0]/pow(rEarth2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[1]-xEarth_SO[1])/pow(rEarth2Particle,3)+xEarth_SO[1]/pow(rEarth2Moon,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_EARTH_)*((x_LOCAL[2]-xEarth_SO[2])/pow(rEarth2Particle,3)+xEarth_SO[2]/pow(rEarth2Moon,3));

    //account for the planetary rotation around the Sun
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

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }


}


#endif /* MERCURY_H_ */
