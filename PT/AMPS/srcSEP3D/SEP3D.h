//$Id$

#ifndef _SEP3D_H_
#define _SEP3D_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

// SPICE is NOT used for this application: include empty substitutes
#include "SpiceEmptyDefinitions.h"


namespace SEP3D {
  // it is based on a generic application named Exosphere
  using namespace Exosphere;

  //--------------------------------------------------------------------------
  //particle tracking condition
  namespace ParticleTracker {
    inline bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
      //only those particles are traced, which are close to the magnetic island
      if (x[0]<9.0E8||x[2]>0.8E7||x[2]<-0.8E7) return false;

      return PIC::ParticleTracker::TrajectoryTrackingCondition_default(x,v,spec,ParticleData);
    }
  }

  //--------------------------------------------------------------------------
  // specific initialization procedures
  void Init_BeforeParser();

  //--------------------------------------------------------------------------
  // methods for the reaction processing specific for the application
  namespace Physics {
    //------------------------------------------------------------------------
    // particle mover
    int Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode, bool FirstBoundaryFlag);
    int Mover_Axisymmetric_SecondOrder(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    //------------------------------------------------------------------------
    // self-explanatory
    void inline TotalParticleAcceleration(double *accl, int spec, long int ptr,
					  double *x, double *v, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
      //......................................................................
      // local values of position, velocity and acceleration
      double x_LOCAL[3], v_LOCAL[3], accl_LOCAL_CLASSIC[3]={0.0};
      double accl_LOCAL_RELATIVISTIC[3]={0.0};
      memcpy(x_LOCAL,x,3*sizeof(double));
      memcpy(v_LOCAL,v,3*sizeof(double));

      //......................................................................
      // calculate Lorentz force if needed
//#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
      // find electro-magnetic field
      double E[3],B[3];
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      // if no coupler is used -> use characteristic values
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere::swB_Typical,3*sizeof(double));
#else 
      // if coupler is used -> get values from it
      PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
      PIC::CPLR::GetBackgroundElectricField(E);
      PIC::CPLR::GetBackgroundMagneticField(B);
#endif//_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 

      //......................................................................
      // calculate acceleraton due to Lorentz force
      double ElCharge=PIC::MolecularData::GetElectricCharge(spec);
      double mass=PIC::MolecularData::GetMass(spec);
      accl_LOCAL_CLASSIC[0]+=ElCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
      accl_LOCAL_CLASSIC[1]+=ElCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
      accl_LOCAL_CLASSIC[2]+=ElCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;
//#endif//_FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

      //......................................................................
      // calculate gravitational force if needed
//#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
      // find the distance to the Sun's center, here it's {0, 0, 0} 
      double r, r2;
      r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
      r=sqrt(r2);
      // calculate acceleration due to gravity
      for (int idim=0;idim<DIM;idim++) {
	accl_LOCAL_CLASSIC[idim]-=GravityConstant*_SUN__MASS_/r2*x_LOCAL[idim]/r;
      }
//#endif//_FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_

      //......................................................................
      // convert classic acceleration to relativistic
      double speed2=0.0, VelDotAccl=0.0;
      double c2=SpeedOfLight*SpeedOfLight;
      for(int idim=0; idim<DIM;idim++){
	speed2    += v_LOCAL[idim]*v_LOCAL[idim];
	VelDotAccl+= v_LOCAL[idim]*accl_LOCAL_CLASSIC[idim]; 
      }
      if(speed2 > c2) exit(__LINE__,__FILE__,"Error: superluminous particle");
      double RelGamma = pow(1.0 - speed2/c2, -0.5);
      for(int idim=0; idim<DIM;idim++)
	accl_LOCAL_RELATIVISTIC[idim]= 
	  (accl_LOCAL_CLASSIC[idim]-VelDotAccl*v_LOCAL[idim]/c2)/RelGamma;

      //copy the local values of the acceleration to the global ones
      memcpy(accl,accl_LOCAL_RELATIVISTIC,3*sizeof(double));
    }
  }
  

  inline long int inject_particle_onto_field_line(int spec){
    //namespace aliases
    namespace PB = PIC::ParticleBuffer;
    namespace FL = PIC::FieldLine;
    // pointer to the particle to be injected
    long int ptr;
    PB::byte* ptrData;
    ptrData = new PB::byte[PB::GetParticleDataLength()];
    
    // set species
    PB::SetI(spec, ptrData);

    // particle mass
    double m0 = PIC::MolecularData::GetMass(spec);
    
    // pick a random field line and set field line id
    int iFieldLine = (int)(FL::nFieldLine * rnd());
    PB::SetFieldLineId(iFieldLine, ptrData);
    
    // get a segment to inject particle
    int iSegment;
    double WeightCorrection;
    FL::FieldLinesAll[iFieldLine].GetSegmentRandom(iSegment,
						   WeightCorrection, spec);
    FL::cFieldLineSegment* Segment = 
      FL::FieldLinesAll[iFieldLine].GetSegment(iSegment);

    // particle is injected onto random location along the segment
    double S = iSegment + rnd();
    PB::SetFieldLineCoord(S, ptrData);

    // set the cartesian coordinate for the particle
    double x[3], v[3];
    FL::FieldLinesAll[iFieldLine].GetCartesian(x, S);
    PB::SetX(x, ptrData);
    
    
    // the numerical approach used in this simulation (guiding motion center)
    // applies to energetic particles only;
    // inject particles with Maxwellain velocity but cut center out
    
    // envelope distribution is uniform on a spherical shell
    const double SpeedMin  = 5e+5;
    const double SpeedMin3 = pow(SpeedMin,3);
    const double SpeedMax  = 1.2e+6;
    const double SpeedMax3 = pow(SpeedMax,3);
    
    // generate random speed
    double Speed = pow(SpeedMin3 + rnd() * (SpeedMax3-SpeedMin3), 1.0/3);
    
    // random direction
    double cosTheta =  1 - 2*rnd();
    double sinTheta =  pow(1-cosTheta*cosTheta, 0.5);
    double Phi = 2 * Pi * rnd();
    
    // the full initial speed
    v[0] = Speed * sinTheta * cos(Phi);
    v[1] = Speed * sinTheta * sin(Phi);
    v[2] = Speed * cosTheta;
    
    // characteristic temperature
    double Tc = 6e+3;

    // local plasma parameters (temperature and velocity)
    double T;
    Segment->GetPlasmaTemperature(S-(int)S, T);
    double U[3];
    Segment->GetPlasmaVelocity(S-(int)S, U);

    // compute the weight correction
    double misc = 0;
    for(int i=0; i<3; i++) misc += (v[i]-U[i])*(v[i]-U[i]);
    double omega = 
      pow(Tc/T, 1.5) *
      exp( m0 / (2*Kbol*Tc) * SpeedMin*SpeedMin - m0 / (2*Kbol*T) * misc);

    WeightCorrection *= omega;
    PB::SetIndividualStatWeightCorrection(WeightCorrection,ptrData);
    
    // find component of the velocity parallel to the B field
    double dir[3];
    Segment->GetDir(dir);
    double vpar = v[0] * dir[0] + v[1] * dir[1] + v[2] * dir[2];
    for(int i=0; i<3; i++) v[i] = vpar * dir[i];
    
    // set only parallel component
    PB::SetV(v, ptrData);
    
    // perpendicular motiion is characterized by particles mag moment
    double KinEnergyPerp = 0.5 * m0 * (Speed*Speed - vpar*vpar);      

    //magnetic field
    double B[3], AbsB;
    Segment->GetMagneticField(S-(int)S, B);
    AbsB = pow(B[0]*B[0]+B[1]*B[1]+B[2]*B[2], 0.5)+1E-15;
    
    //magnetic moment
    PB::SetMagneticMoment(KinEnergyPerp/AbsB, ptrData);
    
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    node=PIC::Mesh::mesh.findTreeNode(x);
    
    long int res=PB::InitiateParticle(x,v,NULL,NULL,
				      ptrData,
				      _PIC_INIT_PARTICLE_MODE__ADD2LIST_,
				      (void*)node);

    delete [] ptrData;
    // double misc =    PB::GetFieldLineCoord(res);
    return res;
    
  }
}

#endif//_SEP3D_H_
