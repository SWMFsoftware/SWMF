//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>


#include <sys/time.h>
#include <sys/resource.h>


#include "pic.h"
#include "constants.h"
#include "Europa.h"
#include "ElectronImpact.h"

//the modeling case
#define _EUROPA_MODEL_MODE__NEAR_PLANET_    0
#define _EUROPA_MODEL_MODE__TAIL_MODEL_     1
#define _EUROPA_MODEL_MODE_ _EUROPA_MODEL_MODE__TAIL_MODEL_


//defiend the modes for included models
#define _EUROPA_MODE_ON_    0
#define _EUROPA_MODE_OFF_   1

#define _EUROPA_IMPACT_VAPORIZATION_MODE_ _EUROPA_MODE_ON_
#define _EUROPA_PSD_MODE_ _EUROPA_MODE_ON_
#define _EUROPA_THERMAL_DESORPTION_MODE_ _EUROPA_MODE_ON_

//sample particles that enter EPD
#define _GALILEO_EPD_SAMPLING_  _EUROPA_MODE_ON_


//the parameters of the domain and the sphere
const double DebugRunMultiplier=4.0;


const double rSphere=_RADIUS_(_TARGET_);

//const double xMaxDomain=400.0; //modeling of the tail
double xMaxDomain=5; //modeling the vicinity of the planet
double yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


//const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
double dxMinSphere=DebugRunMultiplier*10.0/100,dxMaxSphere=DebugRunMultiplier*1.0/10.0;
double dxMinGlobal=1,dxMaxGlobal=1;


//the codes for the soruces processes
#define _ALL_SOURCE_PROCESSES_                -1
#define _IMPACT_VAPORIZATION_SOURCE_PROCESS_   0
#define _PSD_SOURCE_PROCESS_                   1

//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
	double res=0.0;

#if _EUROPA_IMPACT_VAPORIZATION_MODE_ == _EUROPA_MODE_ON_
	if ((SourceProcessCode==-1)||(SourceProcessCode==_IMPACT_VAPORIZATION_SOURCE_PROCESS_)) { //the total impact vaporization flux
		res+=2.6E23;
	}
#endif

#if _EUROPA_PSD_MODE_ == _EUROPA_MODE_ON_
	if ((SourceProcessCode==-1)||(SourceProcessCode==_PSD_SOURCE_PROCESS_)) { //the photon stimulated desorption flux
		res+=3.6E24;
	}
#endif

	return res;
}





//impact vaporization
//const double ProductionRateIVsodium=2.6E23;
//const double TemperatureInjectionIVsodium=2000.0;



//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
        double res,r,l[3] = {1.0,0.0,0.0};
	int idim;
	double SubsolarAngle;

	for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
	if (r > 0.8 * rSphere) for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

	SubsolarAngle=acos(l[0]);

	SubsolarAngle=0.0;
	res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;

  res/=2.1;

  if (r>0.95) {
    if (strcmp(Europa::Mesh::sign,"0x3030203cdedcf30")==0) { //reduced mesh
      //do nothing
    }
    else if (strcmp(Europa::Mesh::sign,"0x203009b6e27a9")==0) { //full mesh
      res/=4.1*2.1;
    }
    else if (strcmp(Europa::Mesh::sign,"new")==0) { //reduced mesh
      res=20.0E3/rSphere;
    }
    else exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

	return rSphere*res;
}


double localResolution(double *x) {
  double res;

  if (strcmp(Europa::Mesh::sign,"0x3030203cdedcf30")==0) { //reduced mesh
    int idim;
    double lnR,r=0.0;

    for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

    r=sqrt(r);

    if (r<2.0*rSphere) return localSphericalSurfaceResolution(x);

    if (r>dxMinGlobal*rSphere) {
      lnR=log(r);
      res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
    }
    else res=dxMinGlobal;
  }
  else if (strcmp(Europa::Mesh::sign,"0x203009b6e27a9")==0) { //full mesh
    int idim;
    double lnR,r=0.0, d1,d2,d3,d=0.0;

    for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);
    r=sqrt(r);
    // accomodation for the O2Plus tail's shape
    d1 = x[0];
    d2 = 0.5*(           x[1] + sqrt(3.0)*x[2]);
    d3 = 0.5*(-sqrt(3.0)*x[1] +           x[2]);
    d  = - rSphere*d1 + 0.2*d2*d2 + 1.6*d3*d3;
    d  = ( (d>0.0) ? 1. : -1.) * sqrt(fabs(d));

    if (r<rSphere) return rSphere*dxMaxGlobal;
    if (r<1.03*rSphere) return localSphericalSurfaceResolution(x);

    if (r>dxMinGlobal*rSphere && d > 1.2*rSphere) {
      lnR=log(r);
      res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
    }
    else res=dxMinGlobal;
  }
  else if (strcmp(Europa::Mesh::sign,"new")==0) { //reduced mesh
    int idim;
    double lnR,r=0.0;

    for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

    r=sqrt(r);

    if (r<0.98*rSphere) return rSphere;
    else if (r<1.05*rSphere) return localSphericalSurfaceResolution(x);

    if (r>dxMinGlobal*rSphere) {
      lnR=log(r);
      res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
    }
    else res=dxMinGlobal;
  }
  else exit(__LINE__,__FILE__,"Error: unknown option");

  return rSphere*res;
}


//set up the local time step
#if  _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_OFF_ 
//use the time step distribution for real model runs
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	double CellSize;

	const double maxInjectionEnergy=1E5*1000*ElectronCharge;
	double CharacteristicSpeed = sqrt(2*maxInjectionEnergy/_MASS_(_O_));

	CellSize=startNode->GetCharacteristicCellSize();

   switch (spec) {
   case _O_PLUS_HIGH_SPEC_:
     CharacteristicSpeed=10.0*1.6E6;
     break;
   case _O_PLUS_THERMAL_SPEC_: case _O_PLUS_SPEC_:
     CharacteristicSpeed=9.6E4;
     break;
 
   case _O2_SPEC_:case _H2O_SPEC_:case _H2_SPEC_:case _H_SPEC_:case _OH_SPEC_:case _O_SPEC_:
     CharacteristicSpeed=1.0e5;
     break;
   case _O2_PLUS_SPEC_:case _H_PLUS_SPEC_:case _H2_PLUS_SPEC_:case _H2O_PLUS_SPEC_:case _OH_PLUS_SPEC_:
     CharacteristicSpeed=1.0e6;
     break;
   default:
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
     if(_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups){
       CharacteristicSpeed=2.0e4;
       break;
     }
#endif//_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
     char error_message[300];
     sprintf(error_message,"unknown species %i", spec);
     exit(__LINE__,__FILE__,error_message);
    }

  return 0.2* 0.3*CellSize/CharacteristicSpeed;
}
#elif _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_ 
//use the time step distribution for the nightly tests
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
         double CellSize;
 
 //      double CharacteristicSpeed_NA=5.0E3;
 
         const double maxInjectionEnergy=1E5*1000*ElectronCharge;
         double CharacteristicSpeed = sqrt(2*maxInjectionEnergy/_MASS_(_O_));
 
 
         //  CharacteristicSpeed*=sqrt(PIC::MolecularData::GetMass(NA)/PIC::MolecularData::GetMass(spec));
 
         CellSize=startNode->GetCharacteristicCellSize();
 
 
 
 
 //      CharacteristicSpeed=1.0E6;
 
 /*if (spec==_O_PLUS_HIGH_SPEC_) CharacteristicSpeed=10.0*1.6E6;
 if (spec==_O_PLUS_THERMAL_SPEC_) CharacteristicSpeed=10.0*9.6E4;*/
 
   switch (spec) {
   case _O_PLUS_HIGH_SPEC_:
     CharacteristicSpeed=10.0*1.6E6;
     break;
   case _O_PLUS_THERMAL_SPEC_: case _O_PLUS_SPEC_:
     CharacteristicSpeed=10.0*9.6E4;
     break;
 
   case _O2_SPEC_:case _H2O_SPEC_:case _H2_SPEC_:case _H_SPEC_:case _OH_SPEC_:case _O_SPEC_:
     CharacteristicSpeed=1.0e4;
     break;
   case _O2_PLUS_SPEC_:case _H_PLUS_SPEC_:case _H2_PLUS_SPEC_:case _H2O_PLUS_SPEC_:case _OH_PLUS_SPEC_:
     CharacteristicSpeed=10*1.0e4;
     break;
   default:
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
     if(_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups){
       CharacteristicSpeed=20.0e3;
       break;
     }
#endif//_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
     char error_message[300];
     sprintf(error_message,"unknown species %i", spec);
     exit(__LINE__,__FILE__,error_message);
    }
 
    return 0.3*CellSize/CharacteristicSpeed;
 }

#endif




double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

	/*
  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
	 */

	return 0.0;

	/*
  if (spec!=NA) return 0.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
	 */
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	bool ExternalFaces[6];
	double ExternalNormal[3],ModelParticlesInjectionRate;
	int nface;

	static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;

	if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
		for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
			startNode->GetExternalNormal(ExternalNormal,nface);
			ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O_PLUS_HIGH_SPEC_);

			if (ModelParticlesInjectionRate>0.0) return true;
		}
	}

	return false;
}



//injection of model particles through the faces of the bounding box
long int  BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  //  long int nInjectedParticles;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;
  
  if ((spec!=_O_PLUS_HIGH_SPEC_)&&(spec!=_O_PLUS_THERMAL_SPEC_)) return 0; //inject only spec=0
  
  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
  double v[3];
  
  
  double ModelParticlesInjectionRate;
  
  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);
    
    
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
	startNode->GetExternalNormal(ExternalNormal,nface);
	TimeCounter=0.0;
	
	ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O_PLUS_THERMAL_SPEC_);
	
	
	if (ModelParticlesInjectionRate>0.0) {
	  ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;
	  
	  PIC::Mesh::mesh.GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);
	  
	  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
	    //generate the new particle position on the face
	    for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];
	    
	    //generate a particle
	    newParticle=PIC::ParticleBuffer::GetNewParticle();
	    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
	    nInjectedParticles++;
	    
	    PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O_PLUS_THERMAL_SPEC_);
	    
	    PIC::ParticleBuffer::SetX(x,newParticleData);
	    PIC::ParticleBuffer::SetV(v,newParticleData);
	    PIC::ParticleBuffer::SetI(spec,newParticleData);
	    
	    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);
	    
	    
	    //inject the particle into the system
	    //          PIC::Mover::MoveParticleTimeStep[spec](newParticle,LocalTimeStep-TimeCounter,startNode);
	    
	  }
	}
	
	
      }
  }
  
  return nInjectedParticles;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
	long int nInjectedParticles=0;

	for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

	return nInjectedParticles;
}

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_

double GetTotalProductionRateUniformNASTRAN(int spec){
  return 1.0e25;
}

bool GenerateParticlePropertiesUniformNASTRAN(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,char* tempParticleData) {
  double ExternalNormal[3]; 
  int idim;
  double rate,TableTotalProductionRate,totalSurface,gamma,cosSubSolarAngle,ProjectedAngle,elementSubSolarAngle[180],r;
  double x[3],n[3],c=0.0,X,total,xmin,xmax,*x0Sphere,norm[3];
  static double positionSun[3];
  double HeliocentricDistance=3.3*_AU_;
  int nAzimuthalSurfaceElements,nAxisSurfaceElements,nAxisElement,nAzimuthalElement, nZenithElement;
  long int totalSurfaceElementsNumber,i;
//  double rSphere=1980.0;
  double area;
  static double productionDistributionUniformNASTRAN[200000];
  static double cumulativeProductionDistributionUniformNASTRAN[200000];
  static bool probabilityFunctionDefinedUniformNASTRAN;
  unsigned int el;
  double x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3],v_LOCAL_IAU_OBJECT[3],v_LOCAL_SO_OBJECT[3];

  { // get random position on Europa's surface
    /**************************************************
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * CHANGE DISTRIBUTION OF INJECION OVER SURFACE !!!
     * ACCOUNT FOR DIFFERENT AREAS OF ELEMENTS !!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     **************************************************/
    el=(int)(rnd()*Europa::Planet->nZenithSurfaceElements*Europa::Planet->nAzimuthalSurfaceElements);
    Europa::Planet->GetSurfaceElementIndex(nZenithElement,nAzimuthalElement, el);
    Europa::Planet->GetSurfaceElementRandomDirection(ExternalNormal,nZenithElement,nAzimuthalElement);
    
    x_LOCAL_IAU_OBJECT[0]=sphereX0[0]+sphereRadius*ExternalNormal[0];
    x_LOCAL_IAU_OBJECT[1]=sphereX0[1]+sphereRadius*ExternalNormal[1];
    x_LOCAL_IAU_OBJECT[2]=sphereX0[2]+sphereRadius*ExternalNormal[2];
    
    //    ExternalNormal[0]*=-1.0;
    //    ExternalNormal[1]*=-1.0;
    //    ExternalNormal[2]*=-1.0;
  }

  
  //transfer the position into the coordinate frame related to the rotating coordinate frame 'MSGR_SO'
  x_LOCAL_SO_OBJECT[0]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[0][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[0][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[1]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[1][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[1][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[1][2]*x_LOCAL_IAU_OBJECT[2]);
  
  x_LOCAL_SO_OBJECT[2]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[2][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[2][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[2][2]*x_LOCAL_IAU_OBJECT[2]);
    

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;
  
  //generate particle's velocity vector in the coordinate frame related to the planet 'IAU_OBJECT'
  double SurfaceTemperature,vbulk[3]={0.0,0.0,0.0};
  //  if(CutCell::BoundaryTriangleFaces[i].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) cosSubSolarAngle=-1; //Get Temperature from night side if in the shadow
  SurfaceTemperature=100;//GetSurfaceTemeprature(cosSubSolarAngle,x_LOCAL_SO_OBJECT);

  double r2Tang=0.0;
  double xFace[3];
//  double vDustInit=2500;
  double angleVelocityNormal=asin(rnd());

  if (spec>=_DUST_SPEC_ && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) { //for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]=vDustInit*ExternalNormal[idim];
  for (idim=0;idim<3;idim++){
        v_LOCAL_IAU_OBJECT[idim]=ElectricallyChargedDust::InitialGrainSpeed*ExternalNormal[idim]*cos(angleVelocityNormal);
     
    }
  /*CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
   *while(xFace[0]==x_LOCAL_SO_OBJECT[0] && xFace[1]==x_LOCAL_SO_OBJECT[1] && xFace[2]==x_LOCAL_SO_OBJECT[2]) CutCell::BoundaryTriangleFaces[i].GetRandomPosition(xFace,1e-4);
   *for (idim=0;idim<3;idim++) r2Tang+=pow(x_LOCAL_SO_OBJECT[idim]-xFace[idim],2.0);
   *for (idim=0;idim<3;idim++) v_LOCAL_IAU_OBJECT[idim]+=vDustInit*sin(angleVelocityNormal)*(x_LOCAL_SO_OBJECT[idim]-xFace[idim])/sqrt(r2Tang);
   */
  }
  else for (idim=0;idim<3;idim++) PIC::Distribution::InjectMaxwellianDistribution(v_LOCAL_IAU_OBJECT,vbulk,SurfaceTemperature,ExternalNormal,spec);
  

  //transform the velocity vector to the coordinate frame 'MSGR_SO'
  v_LOCAL_SO_OBJECT[0]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][2]*x_LOCAL_IAU_OBJECT[2])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][3]*v_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][4]*v_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[3][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[1]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][2]*x_LOCAL_IAU_OBJECT[2])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][3]*v_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][4]*v_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[4][5]*v_LOCAL_IAU_OBJECT[2]);
  
  v_LOCAL_SO_OBJECT[2]=
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][0]*x_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][1]*x_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][2]*x_LOCAL_IAU_OBJECT[2])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][3]*v_LOCAL_IAU_OBJECT[0])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][4]*v_LOCAL_IAU_OBJECT[1])+
    (Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[5][5]*v_LOCAL_IAU_OBJECT[2]);
  
  memcpy(x_SO_OBJECT,x_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(x_IAU_OBJECT,x_LOCAL_IAU_OBJECT,3*sizeof(double));
  memcpy(v_SO_OBJECT,v_LOCAL_SO_OBJECT,3*sizeof(double));
  memcpy(v_IAU_OBJECT,v_LOCAL_IAU_OBJECT,3*sizeof(double));
  
  return true;
}



long int DustInjection(int spec) {
  double ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,TimeCounter=0.0,x_SO_OBJECT[3],x_IAU_OBJECT[3],v_SO_OBJECT[3],v_IAU_OBJECT[3],sphereX0[3]={0.0},sphereRadius=_EUROPA__RADIUS_;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  double ParticleWeightCorrection=1.0;
  bool flag=false;
  int SourceProcessID;

  // ignore non-dust species
  if (!(_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) return 0;

  double totalProductionRate=GetTotalProductionRateUniformNASTRAN(spec);

  const int nMaxInjectedParticles=10*PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  exit(__LINE__,__FILE__,"Error: not implemented!");
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  ModelParticlesInjectionRate=totalProductionRate/ParticleWeight;

  if (ModelParticlesInjectionRate*LocalTimeStep>nMaxInjectedParticles) {
    ParticleWeightCorrection=ModelParticlesInjectionRate*LocalTimeStep/nMaxInjectedParticles;
    ModelParticlesInjectionRate/=ParticleWeightCorrection;
  }

  //definition of indexes TEMPORARY!!!!!
  //int _EXOSPHERE__SOURCE_MAX_ID_VALUE_=2;
  int _EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_=0;
  int _exosphere__SOURCE_MAX_ID_VALUE_=0;


  //calcualte probabilities of each source processes
  double TotalFlux,FluxSourceProcess[1+_exosphere__SOURCE_MAX_ID_VALUE_]; 
  for (int iSource=0;iSource<1+_exosphere__SOURCE_MAX_ID_VALUE_;iSource++) FluxSourceProcess[iSource]=0.0; 
  
  TotalFlux=totalProductionRate;
  
  //only Used defined source here since we only want the Bjorn model so far
  //calculate the source rate due to user defined source functions                                                   
  //Distribution of dust injection correlated with water
  
  FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_]=GetTotalProductionRateUniformNASTRAN(_H2O_SPEC_);

  TotalFlux=GetTotalProductionRateUniformNASTRAN(_H2O_SPEC_);
  
  double CalculatedSourceRate[PIC::nTotalSpecies][1+_exosphere__SOURCE_MAX_ID_VALUE_];
  CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_]=0.0;
  
  
  static double GrainInjectedMass=0.0;
  PIC::Mesh::cDataBlockAMR *block;
  double GrainRadius,GrainMass,GrainWeightCorrection;
  int GrainVelocityGroup;
  
  GrainInjectedMass+=ElectricallyChargedDust::TotalMassDustProductionRate*LocalTimeStep;
  
  while (GrainInjectedMass>0.0) {
    startNode=NULL;
    
    do {
      SourceProcessID=(int)(rnd()*(1+_exosphere__SOURCE_MAX_ID_VALUE_));
    }
    while (FluxSourceProcess[SourceProcessID]/TotalFlux<rnd());
    
    //generate a particle                                                                                             
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
    PIC::ParticleBuffer::SetI(spec,(PIC::ParticleBuffer::byte*)tempParticleData);

    if (SourceProcessID==_EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_) {
      flag=GenerateParticlePropertiesUniformNASTRAN(spec,x_SO_OBJECT,x_IAU_OBJECT,v_SO_OBJECT,v_IAU_OBJECT,sphereX0,sphereRadius,startNode,tempParticleData);
      ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(GrainRadius,GrainWeightCorrection);
      GrainMass=4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity*pow(GrainRadius,3);
      GrainInjectedMass-=GrainMass*ParticleWeight*GrainWeightCorrection;
      SourceProcessID=_EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_;
      if (flag==true) CalculatedSourceRate[spec][_EXOSPHERE_SOURCE__ID__USER_DEFINED__Uniform_]+=ParticleWeightCorrection*ParticleWeight/LocalTimeStep;
    }
    else {
      continue;
    }
    if (flag==false) continue;
    if ((block=startNode->block)->GetLocalTimeStep(_DUST_SPEC_)/LocalTimeStep<rnd()) continue;
    
    //determine the velocity group of the injected grain;
    //calculate additional particle weight correction because the particle will be placed in a different weight group
    GrainVelocityGroup=ElectricallyChargedDust::GrainVelocityGroup::GetGroupNumber(v_SO_OBJECT);
    GrainWeightCorrection*=block->GetLocalTimeStep(_DUST_SPEC_+GrainVelocityGroup)/block->GetLocalTimeStep(_DUST_SPEC_);
    
#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    GrainWeightCorrection*=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_]/PIC::ParticleWeightTimeStep::GlobalParticleWeight[_DUST_SPEC_+GrainVelocityGroup];
#else
    exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif
    
    //determine the surface element of the particle origin                                                            
    PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

    PIC::ParticleBuffer::SetX(x_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::SetV(v_SO_OBJECT,(PIC::ParticleBuffer::byte*)tempParticleData);
    PIC::ParticleBuffer::SetI(_DUST_SPEC_+GrainVelocityGroup,(PIC::ParticleBuffer::byte*)tempParticleData);
    
    Europa::Sampling::SetParticleSourceID(_EXOSPHERE_SOURCE__ID__EXTERNAL_BOUNDARY_INJECTION_,(PIC::ParticleBuffer::byte*)tempParticleData);
    
    ElectricallyChargedDust::SetGrainCharge(0.0,(PIC::ParticleBuffer::byte*)tempParticleData);
    ElectricallyChargedDust::SetGrainMass(GrainMass,(PIC::ParticleBuffer::byte*)tempParticleData);
    ElectricallyChargedDust::SetGrainRadius(GrainRadius,(PIC::ParticleBuffer::byte*)tempParticleData);
    
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(GrainWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
    
    
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);
    
    //determine the initial charge of the dust grain
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
    _PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(x_SO_OBJECT,x_SO_OBJECT,v_SO_OBJECT,spec,newParticle,newParticleData,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode);
#endif
    
    nInjectedParticles++;
    
    //apply condition of tracking the particle
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    PIC::ParticleTracker::InitParticleID(newParticleData);
    PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x_SO_OBJECT,v_SO_OBJECT,spec,newParticleData);
#endif

    //inject the particle into the system                                                                             
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
    if (startNode==NULL) exit(__LINE__,__FILE__,"Error: the node is not defined");
    if ((startNode->Thread!=PIC::ThisThread)||(startNode->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is n\
ot defined");
#endif
    
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode,true);
  }
  
  return nInjectedParticles;
}

long int DustInjection(){
  int spec;
  long int res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=DustInjection(spec);

  return res;
}


#endif


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
	double res=1.0;

	// for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

	return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
	double radiusSphere,*x0Sphere,l[3],r,vNorm,c;
	cInternalSphericalData *Sphere;
	cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
	int idim;

	//   long int newParticle;
	//   PIC::ParticleBuffer::byte *newParticleData;
	//   double ParticleStatWeight,WeightCorrection;


	Sphere=(cInternalSphericalData*)SphereDataPointer;
	startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

	Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

	for (r=0.0,idim=0;idim<DIM;idim++) {
		l[idim]=x[idim]-x0Sphere[idim];
		r+=pow(l[idim],2);
	}

	for (r=sqrt(r),vNorm=0.0,idim=0;idim<DIM;idim++) vNorm+=v[idim]*l[idim]/r;
	if (vNorm<0.0) for (c=2.0*vNorm/r,idim=0;idim<DIM;idim++) v[idim]-=c*l[idim];

	//sample the particle data
	double *SampleData;
	long int nSurfaceElement,nZenithElement,nAzimuthalElement;

	Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
	nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
	SampleData=Sphere->SamplingBuffer+PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


	SampleData[PIC::BC::InternalBoundary::Sphere::sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);


	//   if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<0.9*rSphere*rSphere) exit(__LINE__,__FILE__,"Particle inside the sphere");


	r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);


	//particle-surface interaction
	if (false) { /////(spec!=NA) { //surface reactiona
		exit(__LINE__,__FILE__,"no BC for the space is implemented");

		/*
     ParticleStatWeight=startNode->block->GetLocalParticleWeight(0);
     ParticleStatWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

     //model the rejected H+ and neutralized H
     //1. model rejected SW protons (SPEC=1)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

exit(__LINE__,__FILE__,"ERROR: SetI(1,2) -> looks very strangle");

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(1,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.01*ParticleStatWeight/startNode->block->GetLocalParticleWeight(1);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[1](newParticle,dtTotal,startNode);

     //1. model rejected SW NEUTRALIZED protons (SPEC=2)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(2,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.2*ParticleStatWeight/startNode->block->GetLocalParticleWeight(2);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[2](newParticle,dtTotal,startNode);

     //remove the oroginal particle (s=0)
     PIC::ParticleBuffer::DeleteParticle(ptr);
     return _PARTICLE_DELETED_ON_THE_FACE_;
		 */
	}


	//delete all particles that was not reflected on the surface
//	PIC::ParticleBuffer::DeleteParticle(ptr);
	return _PARTICLE_DELETED_ON_THE_FACE_;
}




double sphereInjectionRate(int spec,void *SphereDataPointer) {
	double res=0.0;

	/*
	if (spec==NA) res=sodiumTotalProductionRate();
	else if (spec==NAPLUS) res=0.0;
	else exit(__LINE__,__FILE__,"Error: the source rate for the species is not determined");
*/

	return res;
}



void amps_init_mesh() {
//	MPI_Init(&argc,&argv);

  //set up the resolution levels of the mesh
  if (strcmp(Europa::Mesh::sign,"0x3030203cdedcf30")==0) { //reduced mesh
   dxMinGlobal=1,dxMaxGlobal=1;

    xMaxDomain=5; //modeling the vicinity of the planet
    yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


   //const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
    dxMinSphere=DebugRunMultiplier*10.0/100,dxMaxSphere=DebugRunMultiplier*1.0/10.0;
  }
  else if (strcmp(Europa::Mesh::sign,"0x203009b6e27a9")==0) { //full mesh
    dxMinGlobal=0.4/2.1,dxMaxGlobal=1;

    xMaxDomain=5; //modeling the vicinity of the planet
    yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


   //const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
    dxMinSphere=DebugRunMultiplier*10.0/100,dxMaxSphere=DebugRunMultiplier*1.0/10.0;
  }
  else if (strcmp(Europa::Mesh::sign,"new")==0) { //full mesh
    dxMinGlobal=0.4/2.1,dxMaxGlobal=1;

    xMaxDomain=1+200E3/_EUROPA__RADIUS_; //modeling the vicinity of the planet
    yMaxDomain=1+200E3/_EUROPA__RADIUS_; //the minimum size of the domain in the direction perpendicular to the direction to the sun


   //const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
    dxMinSphere=20E3,dxMaxSphere=100E3;
  }
  else exit(__LINE__,__FILE__,"Error: unknown option");


PIC::InitMPI();


	//SetUp the alarm
	//  PIC::Alarm::SetAlarm(2000);



	//  VT_OFF();
	//  VT_traceoff();


#ifdef _MAIN_PROFILE_
	initSaturn ("");
#endif


	rnd_seed();



	char inputFile[]="Europa.input";



	MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


	//init the Europa model
	//Exosphere::Init_BeforeParser();
	Europa::Init_BeforeParser();

	//init the particle solver
	PIC::Init_BeforeParser();


	//output the parameters of the implemented physical models
	if (PIC::ThisThread==0) {
	  ElectronImpact::H2O::Print("H2O-ElectronImpact.dat",PIC::OutputDataFileDirectory);
	  ElectronImpact::O2::Print("O2-ElectronImpact.dat",PIC::OutputDataFileDirectory);
	  ElectronImpact::H2::Print("H2-ElectronImpact.dat",PIC::OutputDataFileDirectory);
	  ElectronImpact::O::Print("O-ElectronImpact.dat",PIC::OutputDataFileDirectory);
	  ElectronImpact::H::Print("H-ElectronImpact.dat",PIC::OutputDataFileDirectory);
	}


//	PIC::Parser::Run(inputFile);


	const int InitialSampleLength=600;

	//PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=1600;
	//PIC::RequiredSampleLength=InitialSampleLength; //0;

	Europa::Init_AfterParser();




	//output the PDS energy distribution function
	if (PIC::ThisThread==0) {
		//Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintCumulativeDistributionFunction("CumulativeEnergyDistribution-PSD.dat");
		//Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintDistributionFunction("EnergyDistribution-PSD.dat");

		//cout << Europa::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
	}



	//register the sphere
	static const bool SphereInsideDomain=true;

	if (SphereInsideDomain==true) {
		double sx0[3]={0.0,0.0,0.0};
		cInternalBoundaryConditionsDescriptor SphereDescriptor;
		cInternalSphericalData *Sphere;


		//reserve memory for sampling of the surface balance of sticking species
		long int ReserveSamplingSpace[PIC::nTotalSpecies];

		for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_EUROPA_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


		cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);



		PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
		SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
		Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
		Sphere->SetSphereGeometricalParameters(sx0,rSphere);


		//init the object for distribution of the injection surface elements
    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
      Europa::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


    #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      Europa::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
    #endif

    #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      Europa::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
          Europa::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
    #endif
    }






    Sphere->Radius=_RADIUS_(_TARGET_);
		Sphere->PrintSurfaceMesh("Sphere.dat");
		Sphere->PrintSurfaceData("SpheraData.dat",0);
		Sphere->localResolution=localSphericalSurfaceResolution;
		Sphere->InjectionRate=Europa::SourceProcesses::totalProductionRate;
		Sphere->faceat=0;
		Sphere->ParticleSphereInteraction=Europa::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
		Sphere->InjectionBoundaryCondition=Europa::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;

		Sphere->PrintTitle=Europa::Sampling::OutputSurfaceDataFile::PrintTitle;
		Sphere->PrintVariableList=Europa::Sampling::OutputSurfaceDataFile::PrintVariableList;
		Sphere->PrintDataStateVector=Europa::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

		//set up the planet pointer in Europa model
		Europa::Planet=Sphere;

		Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);
	}

	//init the solver
	PIC::Mesh::initCellSamplingDataBuffer();

	//init the mesh
	cout << "Init the mesh" << endl;

	int maxBlockCellsnumber,minBlockCellsnumber,idim;

	maxBlockCellsnumber=_BLOCK_CELLS_X_;
	if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
	if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);

	minBlockCellsnumber=_BLOCK_CELLS_X_;
	if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
	if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);

	double DomainLength[3],DomainCenterOffset[3],xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

	if (maxBlockCellsnumber==minBlockCellsnumber) {
		for (idim=0;idim<DIM;idim++) {
			DomainLength[idim]=2.0*xMaxDomain*rSphere;
			DomainCenterOffset[idim]=-xMaxDomain*rSphere;
		}
	}
	else {
		if (maxBlockCellsnumber!=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__);
		if (minBlockCellsnumber!=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__);
		if (minBlockCellsnumber!=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__);

		DomainLength[0]=xMaxDomain*rSphere*(1.0+double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_);
		DomainLength[1]=DomainLength[0]*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
		DomainLength[2]=DomainLength[0]*double(_BLOCK_CELLS_Z_)/_BLOCK_CELLS_X_;

		if (DomainLength[1]<2.01*yMaxDomain*_RADIUS_(_TARGET_)) {
			double r;

			printf("Size of the domain is smaller that the radius of the body: the size of the domain is readjusted\n");
			r=2.01*yMaxDomain*_RADIUS_(_TARGET_)/DomainLength[1];

			for (idim=0;idim<DIM;idim++) DomainLength[idim]*=r;
		}

		DomainCenterOffset[0]=-yMaxDomain*rSphere;////-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
		DomainCenterOffset[1]=-DomainLength[1]/2.0;
		DomainCenterOffset[2]=-DomainLength[2]/2.0;
	}

	for (idim=0;idim<DIM;idim++) {
		xmax[idim]=-DomainCenterOffset[idim];
		xmin[idim]=-(DomainLength[idim]+DomainCenterOffset[idim]);
	}


	//generate only the tree
	PIC::Mesh::mesh.AllowBlockAllocation=false;
	PIC::Mesh::mesh.init(xmin,xmax,localResolution);
	PIC::Mesh::mesh.memoryAllocationReport();


	char mesh[200]="";
  bool NewMeshGeneratedFlag=false;
  FILE *fmesh=NULL;

  sprintf(mesh,"amr.sig=%s.mesh.bin",Europa::Mesh::sign);
  fmesh=fopen(mesh,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(mesh);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh.ThisThread==0) {
       PIC::Mesh::mesh.buildMesh();
       PIC::Mesh::mesh.saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }
  }


	cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

 if (NewMeshGeneratedFlag==true) PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");

	PIC::Mesh::mesh.memoryAllocationReport();
	PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
	PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

	PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
	PIC::Mesh::mesh.CreateNewParallelDistributionLists();

	//initialize the blocks
	PIC::Mesh::mesh.AllowBlockAllocation=true;
	PIC::Mesh::mesh.AllocateTreeBlocks();

	PIC::Mesh::mesh.memoryAllocationReport();
	PIC::Mesh::mesh.GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
	PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);
#endif

	//init the volume of the cells'
	PIC::Mesh::mesh.InitCellMeasure();

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


	if (PIC::ThisThread==0) cout << "AMPS' Initialization is complete" << endl;

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


}

void amps_init() {
   int idim;

   //init the PIC solver
   PIC::Init_AfterParser ();
   PIC::Mover::Init();
   

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
   //init the dust model
   ElectricallyChargedDust::Init_AfterParser();
#endif

   //set up the time step
   PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
   PIC::ParticleWeightTimeStep::initTimeStep();
   
   //set up the particle weight
   PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjectionRate;
   if (_O_PLUS_HIGH_SPEC_>=0) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O_PLUS_HIGH_SPEC_);
   if (_O_PLUS_THERMAL_SPEC_>=0) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O_PLUS_THERMAL_SPEC_);
   
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_  ==  _EXOSPHERE_SOURCE__ON_
   //evaluate the weight from the parameters of the sputtering (in the Exosphere sputtering model)
   if (_O_PLUS_THERMAL_SPEC_>=0) {
     PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O2_SPEC_,_O_PLUS_THERMAL_SPEC_,1.0e4);
   }
   else {
     PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O2_SPEC_);
   }
#else
   //copy the weight distribution from that of the thermal ions
   PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O2_SPEC_,_O_PLUS_THERMAL_SPEC_,1.0e4);
#endif
   
   if (_H2O_SPEC_>=0) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H2O_SPEC_);
   if (_H2O_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_H2O_PLUS_SPEC_, _H2O_SPEC_, 1.0);

   if (_H2_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_H2_SPEC_, _H2O_SPEC_, 1.0);
   if (_H2_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_H2_PLUS_SPEC_, _H2_SPEC_, 1.0);

   if (_O_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O_SPEC_, _H2O_SPEC_, 1.0);
   if (_O_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O_PLUS_SPEC_, _O_SPEC_, 1.0);

   if (_OH_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O_SPEC_, _OH_SPEC_, 1.0);
   if (_OH_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_OH_PLUS_SPEC_, _OH_SPEC_, 1.0);

   if (_H_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_H_SPEC_, _H2O_SPEC_, 1.0);
   if (_H_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_H_PLUS_SPEC_, _H_SPEC_, 1.0);


   if (_O2_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_O2_PLUS_SPEC_,_O2_SPEC_,1.0E10*1.0E-7);
//   if (_O2_PLUS_SPEC_>=0) PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(_O2_PLUS_SPEC_,_O_PLUS_THERMAL_SPEC_,1.0);

   #if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
   for (int s=0;s<PIC::nTotalSpecies;s++)
     if (_DUST_SPEC_<=s && s<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)
       PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(s,_H2O_SPEC_,1e-7);
   #endif // _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_



     //PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);
   
   //init weight of the daugter products of the photolytic and electron impact reactions
   for (int spec=0;spec<PIC::nTotalSpecies;spec++) if (PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]<0.0) {
       double yield=0.0;
       
       yield+=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_]*
	 (PhotolyticReactions::H2O::GetSpeciesReactionYield(spec)+ElectronImpact::H2O::GetSpeciesReactionYield(spec,20.0));
       
       yield+=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_O2_SPEC_]*
	 (PhotolyticReactions::O2::GetSpeciesReactionYield(spec)+ElectronImpact::O2::GetSpeciesReactionYield(spec,20.0));
       
       yield/=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_];
       PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(spec,_H2O_SPEC_,yield);
     }
   
   //set photolytic reactions
   //PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(sodiumPhotoionizationReactionProcessor,_O2_SPEC_);
   //PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(sodiumPhotoionizationLifeTime,_O2_SPEC_);
   
   
   //	PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
   //	PIC::Mesh::mesh.outputMeshDataTECPLOT("mesh.data.dat",0);
   
   MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
   if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;
   
   
   
   
   //output final data
   //  PIC::Mesh::mesh.outputMeshDataTECPLOT("final.data.dat",0);
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
   //set the User Definted injection function for dust
   PIC::BC::UserDefinedParticleInjectionFunction=DustInjection;
#endif
   //create the list of mesh nodes where the injection boundary conditions are applied
   PIC::BC::BlockInjectionBCindicatior=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxParticleInjectionIndicator;
   PIC::BC::userDefinedBoundingBlockInjectionFunction=Europa::InjectEuropaMagnetosphericEPDIons::BoundingBoxInjection;
   PIC::BC::InitBoundingBoxInjectionBlockList();
   
   
   
   
   //init the particle buffer
   PIC::ParticleBuffer::Init(10000000);
   //  double TimeCounter=time(NULL);
   int LastDataOutputFileNumber=-1;
   
   
   //init the sampling of the particls' distribution functions
   //const int nSamplePoints=3;
   //double SampleLocations[nSamplePoints][DIM]={{2.0E6,0.0,0.0}, {0.0,2.0E6,0.0}, {-2.0E6,0.0,0.0}};
   
   /* THE DEFINITION OF THE SAMPLE LOCATIONS IS IN THE INPUT FILE
      PIC::DistributionFunctionSample::vMin=-40.0E3;
      PIC::DistributionFunctionSample::vMax=40.0E3;
      PIC::DistributionFunctionSample::nSampledFunctionPoints=500;
      
      PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);
   */
   
   //also init the sampling of the particles' pitch angle distribution functions
   //PIC::PitchAngleDistributionSample::nSampledFunctionPoints=101;
   
   //PIC::PitchAngleDistributionSample::Init(SampleLocations,nSamplePoints);
   
   
   
#if _GALILEO_EPD_SAMPLING_ == _EUROPA_MODE_ON_
   //init sampling points along the s/c trajectory
   const int nFlybySamplePasses=1;
   const double FlybySamplingInterval=120.0*60.0,FlybySamplingIntervalStep=120.0; //in seconds
   const int nSampleSteps=(int)(FlybySamplingInterval/FlybySamplingIntervalStep);
   
   const char *FlybySamplePassesUTC[nFlybySamplePasses]={"1996-12-19T06:00:00"};
   
   SpiceDouble et,lt;
   SpiceDouble state[6];
   int nFlybyPass,n;
   
   /* FIPS POINTING
      INS-236720_FOV_FRAME       = 'MSGR_EPPS_FIPS'
      INS-236720_FOV_SHAPE       = 'CIRCLE'
      INS-236720_BORESIGHT       = ( 0.0, 0.0, 1.0 )
   */
   
   SpiceDouble pointing[3],bsight[3],bsight_INIT[3]={0.0,0.0,1.0};
   SpiceDouble rotate[3][3];
   
   const SpiceInt lenout = 35;
   SpiceChar utcstr[lenout+2];
   
   
   int nFluxSamplePoint=0;
   
   int nTotalFluxSamplePoints=nFlybySamplePasses*nSampleSteps;
   double FluxSampleLocations[nTotalFluxSamplePoints][3];
   double FluxSampleDirections[nTotalFluxSamplePoints][3];
   double Dist;
   
   
   for (nFlybyPass=0;nFlybyPass<nFlybySamplePasses;nFlybyPass++) {
     utc2et_c(FlybySamplePassesUTC[nFlybyPass],&et);
     
     if (PIC::ThisThread==0) {
       cout << "S/C Flyby Sampling: Pass=" << nFlybyPass << ":" << endl;
       //???      cout << "Flux Sample Point\tUTS\t\t\t x[km]\t\ty[km]\t\tz[km]\t\t\t lx\t\tly\t\tlz\t" << endl;
       cout << "Flux Sample Point\tUTS\t\t\tx[km]\t\ty[km]\t\tz[km]\tr[rTarget]\t" << endl;
     }
     
     for (n=0;n<nSampleSteps;n++) {
       //position of the s/c
       spkezr_c("GALILEO ORBITER",et,Exosphere::SO_FRAME,"NONE","EUROPA",state,&lt);
       
       
       //get the pointing vector in the 'MSO' frame
       //memcpy(bsight,bsight_INIT,3*sizeof(double));
       
       //???      pxform_c ("MSGR_EPPS_FIPS",Exosphere::SO_FRAME,et,rotate);
       //???      mxv_c(rotate,bsight,pointing);
       
       //print the pointing information
       if (PIC::ThisThread==0) {
	 et2utc_c(et,"C",0,lenout,utcstr);
	 printf("%i\t\t\t%s\t",nFluxSamplePoint,utcstr);
	 for (idim=0;idim<3;idim++) printf("%e\t",state[idim]);
	 Dist = 0.0;
	 for (idim=0;idim<3;idim++) Dist+=pow(state[idim],2);
	 printf("%.2f\t",sqrt(Dist)/_RADIUS_(_TARGET_)*1E3);
	 cout << "\t";
	 
	 //???        for (idim=0;idim<3;idim++) printf("%e\t",pointing[idim]);
	 cout << endl;
       }
       
       //save the sample pointing information
       for (idim=0;idim<3;idim++) {
	 FluxSampleLocations[nFluxSamplePoint][idim]=state[idim]*1.0E3;
	 //???  FluxSampleDirections[nFluxSamplePoint][idim]=pointing[idim];
       }
       
       //increment the flyby time
       et+=FlybySamplingIntervalStep;
       ++nFluxSamplePoint;
     }
     
     cout << endl;
   }
   
   // PIC::ParticleFluxDistributionSample::Init(FluxSampleLocations,FluxSampleDirections,30.0/180.0*Pi,nTotalFluxSamplePoints);
   
#elif _GALILEO_EPD_SAMPLING_ == _EUROPA_MODE_OFF_
   printf("No Galileo EPD sampling\n");
#else
   exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif
   
   
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
   
   
#if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__ICES_
   //init ICES
   
#ifdef _ICES_CREATE_COORDINATE_LIST_
   /*
    *const char IcesLocationPath[]="";//"/Users/dborovik/MyICES/ICES";
    *const char IcesModelCase[]="";//"Europa09";
    *
    *
    *PIC::CPLR::ICES::createCellCenterCoordinateList();
    *PIC::CPLR::ICES::SetLocationICES(IcesLocationPath);
    *PIC::CPLR::ICES::retriveSWMFdata(IcesModelCase);  ////("EUROPA_RESTART_n070001");
    */
   
   PIC::CPLR::ICES::createCellCenterCoordinateList();
   //PIC::CPLR::ICES::SetLocationICES("/Users/vtenishe/ices/ICES/Models"); //("/Users/vtenishe/CODES/ICES/Models");
   PIC::CPLR::ICES::retriveSWMFdata(); //"Europa09"); //("RESTART_t001.52m"); //("MERCURY_RESTART_n070100");  ////("MERCURY_RESTART_n070001");
#endif
   

  //#ifdef _ICES_LOAD_DATA_
  
  PIC::CPLR::ICES::readSWMFdata(1.0);
  //  PIC::Mesh::mesh.outputMeshDataTECPLOT("ices.data.dat",0);
  

  //output the solar wind ion flux at the palnet's surface
  PIC::CPLR::ICES::PrintSphereSurfaceIonFlux("SurfaceIonFlux.dat",1.05*_RADIUS_(_TARGET_));
  PIC::CPLR::ICES::EvaluateSurfaceIonFlux(1.05);
  
  PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test.dat",0);
  
  Exosphere::SourceProcesses::Init();
  PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test-1.dat",0);
  

  #elif _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_
  //TECPLOT
  //read the background data
    if (PIC::CPLR::DATAFILE::BinaryFileExists("EUROPA-BATSRUS")==true)  {
      PIC::CPLR::DATAFILE::LoadBinaryFile("EUROPA-BATSRUS");
    }
    else {
      double xminTECPLOT[3]={-5.1,-5.1,-5.1},xmaxTECPLOT[3]={5.1,5.1,5.1};

      double RotationMatrix_BATSRUS2AMPS[3][3]={ { 1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      
      //  1  0  0
      //  0  1  0
      //  0  0  1
      
      PIC::CPLR::DATAFILE::TECPLOT::SetRotationMatrix_DATAFILE2LocalFrame(RotationMatrix_BATSRUS2AMPS);
      
      PIC::CPLR::DATAFILE::TECPLOT::UnitLength=_EUROPA__RADIUS_;
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsXYZ(xminTECPLOT,xmaxTECPLOT);
      PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsSPHERICAL(1.001,10.0);

      PIC::CPLR::DATAFILE::TECPLOT::DataMode=PIC::CPLR::DATAFILE::TECPLOT::DataMode_SPHERICAL;
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedVelocityVariableData(5,1.0E3);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedIonPressureVariableData(11,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedMagneticFieldVariableData(8,1.0E-9);
      PIC::CPLR::DATAFILE::TECPLOT::SetLoadedDensityVariableData(4,1.0E6);
      PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=11;
      PIC::CPLR::DATAFILE::TECPLOT::ImportData("./3d__mhd_3_n00045039-extracted.plt"); 
      
      PIC::CPLR::DATAFILE::SaveBinaryFile("EUROPA-BATSRUS");
    }

  #else
    exit(__LINE__,__FILE__,"ERROR: unrecognized datafile reader mode");

  #endif //_PIC_COUPLER_DATAFILE_READER_MODE_

#endif //_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
	

    //init particle weight of neutral species that primary source is sputtering


#if  _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE_ == _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_MODE__USER_SOURCE_RATE_ 
	if(_O2_SPEC_>=0)
	  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O2_SPEC_);
#endif
#endif


	// initialize sputtering physical model
#if _SPUTTERING__MODE_ == _SPUTTERING__ON_
	Sputtering::Init();
#endif



  if (_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE_==_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE__TECPLOT_ASCII_) {
    PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);
  }



}

	//time step

void amps_time_step () {

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int idim;

  //determine the parameters of the orbital motion of Europa
  SpiceDouble StateBegin[6],StateEnd[6],lt;
  double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;
  //    int idim;
  
  spkezr_c("Europa",Europa::OrbitalMotion::et,Exosphere::SO_FRAME,"none","Jupiter",StateBegin,&lt);
  Europa::OrbitalMotion::et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_];
  spkezr_c("Europa",Europa::OrbitalMotion::et,Exosphere::SO_FRAME,"none","Jupiter",StateEnd,&lt);
  
  for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
    StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
    StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;
    
    rBegin+=pow(StateBegin[idim],2);
    rEnd+=pow(StateEnd[idim],2);
    
    Europa::xEuropa[idim]=StateBegin[idim];
    Europa::vEuropa[idim]=StateBegin[3+idim];
  }
  
  rBegin=sqrt(rBegin);
  rEnd=sqrt(rEnd);
  
  for (idim=0;idim<3;idim++) {
    lBegin[idim]=StateBegin[idim]/rBegin;
    lEnd[idim]=StateEnd[idim]/rEnd;
    
    c0+=StateBegin[3+idim]*lBegin[idim];
    c1+=StateEnd[3+idim]*lEnd[idim];
  }
  
  Europa::xEuropaRadial=rBegin;
  Europa::vEuropaRadial=0.5*(c0+c1);
  
  //calculate TAA
  {
    double EccentricityVector[3];
    double Speed2,a,c,absEccentricity,cosTAA;
    
    const double GravitationalParameter=GravityConstant*_MASS_(_SUN_);
    
    Speed2=Europa::vEuropa[0]*Europa::vEuropa[0]+Europa::vEuropa[1]*Europa::vEuropa[1]+Europa::vEuropa[2]*Europa::vEuropa[2];
    c=Europa::xEuropa[0]*Europa::vEuropa[0]+Europa::xEuropa[1]*Europa::vEuropa[1]+Europa::xEuropa[2]*Europa::vEuropa[2];
    
    for (idim=0,absEccentricity=0.0,a=0.0;idim<3;idim++) {
      EccentricityVector[idim]=Speed2/GravitationalParameter*Europa::xEuropa[idim] - c/GravitationalParameter*Europa::vEuropa[idim] - Europa::xEuropa[idim]/rBegin;
      absEccentricity+=EccentricityVector[idim]*EccentricityVector[idim];
      a+=EccentricityVector[idim]*Europa::xEuropa[idim];
    }
    
    absEccentricity=sqrt(absEccentricity);
    cosTAA=a/(absEccentricity*rBegin);
    
    Europa::OrbitalMotion::TAA=acos(cosTAA);
    if (c<0.0) Europa::OrbitalMotion::TAA=2.0*Pi-Europa::OrbitalMotion::TAA;
  }
  
  
  //the index varying from 0 to 1!!!!!! because the velocity must be perpendicular to z-axis, which is the axis of rotation
  for (idim=0;idim<2;idim++) {
    vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
    vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
  }
  
  vTangentialBegin=sqrt(vTangentialBegin);
  vTangentialEnd=sqrt(vTangentialEnd);
  
  Europa::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);
  
  
  //determine direction to the Sun and rotation angle in the coordinate frame related to Europa
  SpiceDouble state[6],l=0.0;
  
  spkezr_c("SUN",Europa::OrbitalMotion::et-0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_O2_SPEC_],Exosphere::IAU_FRAME,"none","EUROPA",state,&lt);

  for (idim=0;idim<3;idim++) l+=pow(state[idim],2);
  
  for (l=sqrt(l),idim=0;idim<3;idim++) {
    Europa::OrbitalMotion::SunDirection_IAU_EUROPA[idim]=state[idim]/l;
  }
  
  Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=acos(Europa::OrbitalMotion::SunDirection_IAU_EUROPA[1]);
  if (Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0]<0.0) Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=2.0*Pi-Europa::OrbitalMotion::PlanetAxisToSunRotationAngle;
  
  
  //Update the transformation matrixes and rotation vector
  Europa::OrbitalMotion::UpdateTransformationMartix();
  Europa::OrbitalMotion::UpdateRotationVector_SO_J2000();
  Europa::OrbitalMotion::UpdateSunJupiterLocation();
  
#else
		double xEuropaTest[3] = {-7.27596e-09, - 6.76112e+08, - 2.76134e+06}; 
		std::copy(&xEuropaTest[0], &xEuropaTest[0]+3, &Europa::xEuropa[0]);
		
		double vEuropaTest[3] = { 0,            77.5556,       92.079};
		std::copy(&vEuropaTest[0], &vEuropaTest[0]+3, &Europa::vEuropa[0]);

		double SunDirTest[3]  = {-0.364833, 0.931007, -0.0111228};
		std::copy(&SunDirTest[0], &SunDirTest[0]+3, &Europa::OrbitalMotion::SunDirection_IAU_EUROPA[0]);
		
		Europa::OrbitalMotion::et  = -9.57527e+07;
		Europa::OrbitalMotion::TAA =  3.14159;
		Europa::OrbitalMotion::CoordinateFrameRotationRate = 5.5426e-10;
		Europa::OrbitalMotion::PlanetAxisToSunRotationAngle= 5.90955;

/*
		SpiceDouble TransMatrix1[6][6] = {
		  {-0.0307892,   0.999518,    0.00388983,  0,         0,         0}, 
		  {-0.999503,   -0.0307619,  -0.00689729,  0,         0,         0},
		  {-0.0067743,  -0.00410026,  0.999969,    0,         0,         0},
		  {-3.05046e-07,-8.84681e-09,-1.41289e-07,-0.0307892, 0.999518,  0.00388983},
		  { 9.95762e-09,-3.05672e-07,-7.9686e-08, -0.999503, -0.0307619,-0.00689729},
		  {-8.27439e-08, 1.367e-07,  -2.56459e-14,-0.0067743,-0.00410026,0.999969}};
		std::copy(&TransMatrix1[0][0],&TransMatrix1[0][0]+36, &Europa::OrbitalMotion::SO_to_IAU_TransformationMartix[0][0]);

		SpiceDouble TransMatrix2[6][6] = {
		  {-0.0307892,  -0.999503,   -0.0067743,   0,          0,         0},
		  { 0.999518,   -0.0307619,  -0.00410026,  0,          0,         0},
		  { 0.00388983, -0.00689729,  0.999969,    0,          0,         0},
		  {-3.05046e-07, 9.95762e-09,-8.27439e-08,-0.0307892, -0.999503, -0.0067743},
		  {-8.84681e-09,-3.05672e-07, 1.367e-07,   0.999518,  -0.0307619,-0.00410026},
		  {-1.41289e-07,-7.9686e-08, -2.56459e-14, 0.00388983,-0.00689729,0.999969}};
		std::copy(&TransMatrix2[0][0],&TransMatrix2[0][0]+36,&Europa::OrbitalMotion::IAU_to_SO_TransformationMartix[0][0]);

		SpiceDouble TransMatrix3[6][6] = {
		  {-0.364833,   0.931007,  -0.0111228,   0,         0,        0},
		  {-0.931064,  -0.364856,  -2.77556e-17, 0,         0,        0},
		  {-0.00405822, 0.0103561,  0.999938,    0,         0,        0},
		  { 1.90558e-05,7.46742e-06,1.08988e-09,-0.364833,  0.931007,-0.0111228},
		  {-7.46787e-06,1.9057e-05,-1.07033e-24,-0.931064, -0.364856,-2.77556e-17},
		  { 2.12366e-07,8.2049e-08, 1.21233e-11,-0.00405822,0.0103561,0.999938}};
		std::copy(&TransMatrix3[0][0],&TransMatrix3[0][0]+36,&Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[0][0]);
*/

#endif
  
  //make the time advance
  static int LastDataOutputFileNumber=0;
  
  
  PIC::TimeStep();
  
  //     cell=(node->block!=NULL) ? node->block->GetCenterNode(nd) : NULL;
  
  
  if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
    PIC::RequiredSampleLength*=2;
    if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;
    
    
    LastDataOutputFileNumber=PIC::DataOutputFileNumber;
    if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
  }
  
  

}
