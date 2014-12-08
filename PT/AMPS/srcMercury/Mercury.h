/*
 * Mercury.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */

#ifndef MERCURY_H_
#define MERCURY_H_

#include <math.h>

#include "pic.h"

#include "SpiceUsr.h"
#include "Exosphere.h"
#include "constants.h"

//path to the SPICE Kernels directory
//const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/SPICE/Kernels"; 

//SPICE Kernels to be loaded
/*
const int nFurnishedSPICEkernels=5+12+1;
const char SPICE_Kernels[nFurnishedSPICEkernels][_MAX_STRING_LENGTH_PIC_]={"MESSENGER/kernels/spk/msgr_de405_de423s.bsp","MESSENGER/kernels/fk/msgr_dyn_v600.tf","NAIF/naif0010.tls","MESSENGER/kernels/pck/pck00009_MSGR_v10.tpc","MESSENGER/kernels/fk/msgr_v210.tf",
   "MESSENGER/kernels/ik/msgr_epps_v100.ti","MESSENGER/kernels/ck/msgr20110413.bc","MESSENGER/kernels/ck/msgr20110414.bc","MESSENGER/kernels/ck/msgr20110415.bc","MESSENGER/kernels/ck/msgr20110416.bc","MESSENGER/kernels/ck/msgr20110417.bc","MESSENGER/kernels/ck/msgr20110418.bc","MESSENGER/kernels/ck/msgr20110419.bc","MESSENGER/kernels/ck/msgr20110420.bc","MESSENGER/kernels/ck/msgr20110421.bc",
     "MESSENGER/kernels/sclk/messenger_1486.tsc","MESSENGER/kernels/spk/msgr_20040803_20140823_od266sc_0.bsp","OTHER/GSE.tf"};
*/

namespace Mercury {
  using namespace Exosphere;

  namespace Sampling {
    using namespace Exosphere::Sampling;

    //calculate the column integrals along the limb direction
    namespace SubsolarLimbColumnIntegrals {
      //sample values of the column integrals at the subsolar limb as a function of the phase angle
      const int nPhaseAngleIntervals=50;
      const double dPhaseAngle=Pi/nPhaseAngleIntervals;

      //sample the altitude variatin of the column integral in the solar direction from the limb
      const double maxAltitude=50.0; //altititude in radii of the object
      const double dAltmin=0.001,dAltmax=0.001*maxAltitude; //the minimum and maximum resolution along the 'altitude' line
      extern double rr;
      extern long int nSampleAltitudeDistributionPoints;

      //offsets in the sample for individual parameters of separate species
      extern int _NA_EMISSION_5891_58A_SAMPLE_OFFSET_,_NA_EMISSION_5897_56A_SAMPLE_OFFSET_,_NA_COLUMN_DENSITY_OFFSET_;

/*
      struct cSampleAltitudeDistrubutionBufferElement {
        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
      };

      extern cSampleAltitudeDistrubutionBufferElement *SampleAltitudeDistrubutionBuffer;
*/


      extern SpiceDouble etSampleBegin;
      extern int SamplingPhase;
      extern int firstPhaseRadialVelocityDirection;

      extern int nOutputFile;

      class cSampleBufferElement {
      public:
//        double NA_ColumnDensity,NA_EmissionIntensity__5891_58A,NA_EmissionIntensity__5897_56A;
//        cSampleAltitudeDistrubutionBufferElement *AltitudeDistrubutionBuffer;

        double *SampleColumnIntegrals;
        double **AltitudeDistributionColumnIntegrals;

        int nSamples;
        double JulianDate;

        cSampleBufferElement() {
          SampleColumnIntegrals=NULL,AltitudeDistributionColumnIntegrals=NULL,nSamples=0,JulianDate=0;
        }
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


    //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

    //get the radiation pressure acceleration
    if (spec==_NA_SPEC_) {
      if ((x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2]>_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_))||(x_LOCAL[0]>0.0)) { //calculate the radiation pressure force
        double rHeliocentric,vHeliocentric,radiationPressureAcceleration;

        //calcualte velocity of the particle in the "Frozen SO Frame"; IMPORTANT: SO rotates only around z-direction!!!!
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
//      char *offset;
      int i,j,k;
//      PIC::Mesh::cDataCenterNode *CenterNode;
      double E[3],B[3];

//      exit(__LINE__,__FILE__,"check the numbers!");

      if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  #endif

/*
      CenterNode=startNode->block->GetCenterNode(nd);
      offset=CenterNode->GetAssociatedDataBufferPointer();

      if (*((int*)(offset+PIC::CPLR::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
        memcpy(E,offset+PIC::CPLR::ICES::ElectricFieldOffset,3*sizeof(double));
        memcpy(B,offset+PIC::CPLR::ICES::MagneticFieldOffset,3*sizeof(double));
      }
      else {
        memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
        memcpy(B,Exosphere_swB_Typical,3*sizeof(double));
      }
*/

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere_swB_Typical,3*sizeof(double));
#else 
      PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);
      PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
#endif


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
    double rSun2Moon,rSun2Particle;

    rSun2Moon=xObjectRadial;
    rSun2Particle=sqrt(pow(x_LOCAL[0]-xObjectRadial,2)+pow(x_LOCAL[1],2)+pow(x_LOCAL[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_SUN_)*((x_LOCAL[0]-xObjectRadial)/pow(rSun2Particle,3)+xObjectRadial/pow(rSun2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[1]/pow(rSun2Particle,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[2]/pow(rSun2Particle,3));


    if (isnan(accl_LOCAL[0])||isnan(accl_LOCAL[1])||isnan(accl_LOCAL[2])) exit(__LINE__,__FILE__,"Error in calculation of the acceleration");


/*    double rSolarVector[3],r2Solar,rSolar;


    exit(__LINE__,__FILE__,"Check the gravity correction");

    rSolarVector[0]=xObjectRadial-x_LOCAL[0];
    rSolarVector[1]=x_LOCAL[1],rSolarVector[2]=x_LOCAL[2];

    r2Solar=(rSolarVector[0]*rSolarVector[0])+(rSolarVector[1]*rSolarVector[1])+(rSolarVector[2]*rSolarVector[2]);
    rSolar=sqrt(r2Solar);

    accl_LOCAL[0]+=GravityConstant*_MASS_(_SUN_)*(1.0/r2Solar*rSolarVector[0]/rSolar-1.0/pow(xObjectRadial,2)); //x-axis in solar coordinate frame and SO have opposite direction
    accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)/r2Solar*rSolarVector[1]/rSolar;
    accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)/r2Solar*rSolarVector[2]/rSolar;*/

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


  inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    static const double LifeTime=3600.0*5.8/pow(0.4,2);


    //only sodium can be ionized
    if (spec!=_NA_SPEC_) {
      PhotolyticReactionAllowedFlag=false;
      return -1.0;
    }

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
    spec=_NA_PLUS_SPEC_;

    PIC::ParticleBuffer::SetI(spec,ParticleData);
    if (_NA_PLUS_SPEC_>=0) return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;

    return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
  }

}


#endif /* MERCURY_H_ */
