/*
 * ProtostellarNebula.h
 *
 *  Created on: Jun 21, 2012
 *      Author: vtenishe
 */
//$Id$

#ifndef _PROTOSTELLARNEBULA_H_
#define _PROTOSTELLARNEBULA_H_

#include <math.h>

#include "pic.h"


#include "Exosphere.h"
#include "constants.h"


#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif



namespace ProtostellarNebula {
  using namespace Exosphere;

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
//      PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,ImpactVaporization_SourceTemeprature[spec],ExternalNormal,spec);



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



//defined the forces that acts upon a particle on
#define _FORCE_GRAVITY_MODE_ _PIC_MODE_OFF_
#define _FORCE_LORENTZ_MODE_ _PIC_MODE_OFF_
#define _FORCE_FRAMEROTATION_MODE_ _PIC_MODE_OFF_

  void inline TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
    double x_LOCAL[3],v_LOCAL[3],accl_LOCAL[3]={0.0,0.0,0.0};

    memcpy(x_LOCAL,x,3*sizeof(double));
    memcpy(v_LOCAL,v,3*sizeof(double));

    accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

    exit(__LINE__,__FILE__,"ERROR: Lorentz force not implemented");

    //********************************************************
    // the following code is copied from srcEuropa/Europa.h
    //********************************************************
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
    
    PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
    PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);

    double ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    double mass=PIC::MolecularData::GetMass(spec);

    accl_LOCAL[0]+=ElectricCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
    accl_LOCAL[1]+=ElectricCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
    accl_LOCAL[2]+=ElectricCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;

#endif


#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
  //the gravity force
  double r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x_LOCAL[idim]/r;
  }
#endif

#if _FORCE_FRAMEROTATION_MODE_ == _PIC_MODE_ON_
  // by default rotation period is ~25 days (freq = 4.63E-7 sec^-1)
  // frame angular velocity
  static const double Omega        = 2.0*Pi*4.63E-7;//rad/sec 
  // frame angular velocity x 2
  static const double TwoOmega     = 2.0*Omega;
  // frame angular velocity squared
  static const double SquaredOmega = Omega*Omega;
  double aCen[3],aCorr[3];

  aCen[0] = SquaredOmega * x_LOCAL[0];
  aCen[1] = SquaredOmega * x_LOCAL[1];
  aCen[2] = 0.0;


  aCorr[0] =   TwoOmega * v_LOCAL[1];
  aCorr[1] = - TwoOmega * v_LOCAL[0];
  aCorr[2] =   0.0;

  accl_LOCAL[0]+=aCen[0]+aCorr[0];
  accl_LOCAL[1]+=aCen[1]+aCorr[1];
  accl_LOCAL[2]+=aCen[2]+aCorr[2];

#endif

    //copy the local value of the acceleration to the global one
    memcpy(accl,accl_LOCAL,3*sizeof(double));
  }


  inline double ExospherePhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {
    static const double LifeTime=3600.0*5.8/pow(0.4,2);

    // no photoionization for now
    PhotolyticReactionAllowedFlag=false;
    return -1.0;

  }

  inline int ExospherePhotoionizationReactionProcessor(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {

    // no photoionization for now
    return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  }

}


#endif /* PROTOSTELLARNEBULA_H_ */
