
#ifndef CCMC_H
#define CCMC_H

//$Id$

#include <math.h>

#include "pic.h"
#include "Exosphere.h"
#include "constants.h"

#include "ccmc.dfn"

namespace CCMC {
  using namespace Exosphere;

  //the size of the computational domain
  namespace Domain {
    static const double xmin[]={0.0,0.0,0.0};
    static const double xmax[]={0.0,0.0,0.0};
  }

  //the condition of the particle trajectory tracking
  namespace ParticleTracker {
    inline bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {

/*      //only those solar wind ions are traced, which trajectories are close to the surface of Mercury
      if (spec==_H_PLUS_SPEC_) {
        if (x[1]*x[1]+x[2]*x[2]>pow(2.0*_RADIUS_(_TARGET_),2)) return false;
      }*/

      return PIC::ParticleTracker::TrajectoryTrackingCondition_default(x,v,spec,ParticleData);
    }
  }

  //the total acceleration acting on a particle
  //double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

/*
    memcpy(accl,accl_LOCAL,3*sizeof(double));
    return;
*/


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

    //the Lorentz force
    double elCharge;

    if ((elCharge=PIC::MolecularData::GetElectricCharge(spec))>0.0) {
      long int nd;
      int i,j,k;
      double E[3],B[3];

      if ((nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

      #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
      #endif


#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere_swB_Typical,3*sizeof(double));
#else 
      PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);
      PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
#endif

      elCharge/=PIC::MolecularData::GetMass(spec);

      accl_LOCAL[0]+=elCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1]);
      accl_LOCAL[1]+=elCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0]);
      accl_LOCAL[2]+=elCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0]);

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

}


#endif
