/*
 * Titan.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */
//$Id$

#ifndef _TITAN_H_
#define _TITAN_H_

#include <math.h>

#include "pic.h"

#include <stdexcept>
#include <cstdlib>
#include "Exosphere.h"
#include "constants.h"


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif



namespace Titan {
  using namespace Exosphere;

  namespace tgitm_exobase {
	
	const int ndes=5184,nclmns=7,nintp=nclmns-3,Tnum = 800.0;
	const int Snum = 1;
	extern double tgitm_grid[ndes][nclmns],interp_val[nintp],
	maxflx[PIC::nTotalSpecies],totalflux[PIC::nTotalSpecies];
	
	void read_tgitm();
	
	void tgitm_interpolate(double polar, double azimuth);
	
	inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
			return totalflux[spec];
		}

	inline bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
	double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,
	double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
		
		unsigned int idim;
		long int nZenithElement,nAzimuthalElement;
		int el;
		double r,vbulk[3]={0.0,0.0,0.0},ExternalNormal[3];
		//'x' is the position of a particle in the coordinate frame related to the planet 'IAU_OBJECT'
		double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];
		SpiceDouble xform[6][6];
		double pol,azi,rrr,r4,r3;
		double  flx=0.0, ParticleWeight;
		double ParticleWeightCorrection=1.0;
		const double _N2__MASS_ = 4.64950898e-26;

		memcpy(xform,OrbitalMotion::IAU_to_SO_TransformationMartix,36*sizeof(double));

//Geenrate new particle position
//First the generated particle position must be consistent with local flux
		while( rnd() > flx/maxflx[spec] ){ 
			r=0.0;
			for (idim=0;idim<DIM;idim++) {
				ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
				r+=pow(ExternalNormal[idim],2);
			}
			r=sqrt(r);
			azi = atan2(ExternalNormal[1],ExternalNormal[0]);
			if(azi<0.0)azi=azi+2.0*Pi; 
			pol = acos(ExternalNormal[2]/r);
			tgitm_interpolate(pol, azi);
			flx=interp_val[spec];
		}

//		
//###### Particle has been accepted
			
		for (idim=0;idim<DIM;idim++) {
			ExternalNormal[idim]/=r;
			x_LOCAL_IAU_OBJECT[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
		}
		//transfer the position into the coordinate frame related to the rotting coordinate frame 'MSGR_SO'
		x_LOCAL_SO_OBJECT[0]=xform[0][0]*x_LOCAL_IAU_OBJECT[0]+xform[0][1]*x_LOCAL_IAU_OBJECT[1]+xform[0][2]*x_LOCAL_IAU_OBJECT[2];
		x_LOCAL_SO_OBJECT[1]=xform[1][0]*x_LOCAL_IAU_OBJECT[0]+xform[1][1]*x_LOCAL_IAU_OBJECT[1]+xform[1][2]*x_LOCAL_IAU_OBJECT[2];
		x_LOCAL_SO_OBJECT[2]=xform[2][0]*x_LOCAL_IAU_OBJECT[0]+xform[2][1]*x_LOCAL_IAU_OBJECT[1]+xform[2][2]*x_LOCAL_IAU_OBJECT[2];

		//determine if the particle belongs to this processor
		startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
		if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
		
		//generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
		
		//Switch statement used to implement numerical velocity distribution for better statistics
		double speed=0.0;
		PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,interp_val[3],ExternalNormal,spec);
		/*switch (spec) {
		case _N2_SPEC_:
		
				PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,Tnum,ExternalNormal,spec);
				for (idim=0;idim<DIM;idim++) {
					speed+=pow(v_LOCAL_IAU_OBJECT[idim],2);
				}
				speed=sqrt(speed);
				
				//correction factor based on ratios of maxwelliams Tnum temperature of hot distribution and 
				//interp_val true temperature of surface element 
				ParticleWeightCorrection=exp(-_N2__MASS_*speed*speed/2.0/Kbol/interp_val[3])/exp(-_N2__MASS_*speed*speed/2.0/Kbol/Tnum);
				  
				  ParticleWeightCorrection=pow((1.0/interp_val[3])/(1.0/Tnum),1.5)
				  *exp(-_N2__MASS_*speed*speed/2.0/Kbol/interp_val[3])/exp(-_N2__MASS_*speed*speed/2.0/Kbol/Tnum);
				//cout<<speed<<'\t'<<ParticleWeightCorrection<<endl;
				PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
		  break;
		
		case _CH4_SPEC_:
			PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,interp_val[3],ExternalNormal,spec);
			break;
		case _H2_SPEC_:
			PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,interp_val[3],ExternalNormal,spec);
			break;
		default:
		  exit(__LINE__,__FILE__,"Error: The speed for the species is not defined");
		}*/
		


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
	
  namespace Sampling {
    using namespace Exosphere::Sampling;

  }

  namespace UserdefinedSoruce {
    inline double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {
      return 1.0E20;
    }

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
//   PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,ImpactVaporization_SourceTemeprature[spec],ExternalNormal,spec);



//DEBUG -> injected velocity is normal to the surface

for (int i=0;i<3;i++)  v_LOCAL_IAU_OBJECT[i]=-ExternalNormal[i]*1.0E3;
//END DEBUG



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

  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

/*    memcpy(accl,accl_LOCAL,3*sizeof(double));
    return;*/

    //the gravity force
    double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
    double r=sqrt(r2);
    int idim;

    for (idim=0;idim<DIM;idim++) {
      accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
    }

    //the Lorentz force
    double ElectricCharge, mass,B[3],E[3];

    ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    mass          =PIC::MolecularData::GetMass(spec);

    if (ElectricCharge!=0.0) {
      PIC::CPLR::InitInterpolationStencil(x_LOCAL,startNode);
      PIC::CPLR::GetBackgroundMagneticField(B);
      PIC::CPLR::GetBackgroundElectricField(E);

      accl_LOCAL[0]+=ElectricCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
      accl_LOCAL[1]+=ElectricCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
      accl_LOCAL[2]+=ElectricCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;
    }


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    //correct the gravity acceleration: accout for solar gravity of the particle location
    double rSun2Moon,rSun2Particle;

    rSun2Moon=xObjectRadial;
    rSun2Particle=sqrt(pow(x_LOCAL[0]-xObjectRadial,2)+pow(x_LOCAL[1],2)+pow(x_LOCAL[2],2));

    accl_LOCAL[0]-=GravityConstant*_MASS_(_SUN_)*((x_LOCAL[0]-xObjectRadial)/pow(rSun2Particle,3)+xObjectRadial/pow(rSun2Moon,3));
    accl_LOCAL[1]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[1]/pow(rSun2Particle,3));
    accl_LOCAL[2]-=GravityConstant*_MASS_(_SUN_)*(x_LOCAL[2]/pow(rSun2Particle,3));


    if (isnan(accl_LOCAL[0])||isnan(accl_LOCAL[1])||isnan(accl_LOCAL[2])) exit(__LINE__,__FILE__,"Error in calculation of the acceleration");


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
#endif //_EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }


  inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {
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
