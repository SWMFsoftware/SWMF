//the particle class
#include "pic.h"
#include "constants.h"

//#include "CG-BC.h"

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

//$Id$


//#include <saturn.h>

//#define _MAIN_PROFILE_



//#include "vt_user.h"
//#include <VT.h>




//the modeling case
#define _MERCURY_MODEL_MODE__NEAR_PLANET_    0
#define _MERCURY_MODEL_MODE__TAIL_MODEL_     1

#define _MERCURY_MODEL_MODE_ _MERCURY_MODEL_MODE__TAIL_MODEL_

//Check the mesh consistency
//#define _CHECK_MESH_CONSISTENCY_
//#define _ICES_CREATE_COORDINATE_LIST_
//#define _ICES_LOAD_DATA_


//defiend the modes for included models
#define _MERCURY_MODE_ON_    0
#define _MERCURY_MODE_OFF_   1

#define _MERCURY_IMPACT_VAPORIZATION_MODE_ _MERCURY_MODE_ON_
#define _MERCURY_PSD_MODE_ _MERCURY_MODE_ON_
#define _MERCURY_THERMAL_DESORPTION_MODE_ _MERCURY_MODE_ON_

//sample particles that enter FIPS
#define _MERCURY_FIPS_SAMPLING_  _MERCURY_MODE_OFF_

//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;


const double rSphere=1980.0; //_RADIUS_(_TARGET_);


const double xMaxDomain=3.0; //modeling of the tail
//const double xMaxDomain=20; //modeling the vicinity of the planet
const double yMaxDomain=3.0; //the minimum size of the domain in the direction perpendicular to the direction to the sun



/*const double xMaxDomain=4.0; //modeling of the tail
//const double xMaxDomain=20; //modeling the vicinity of the planet
const double yMaxDomain=4.0; //the minimum size of the domain in the direction perpendicular to the direction to the sun
*/


const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100,dxMaxSphere=DebugRunMultiplier*2.0/10.0;

/*---------------------------------------------------------------------------------*/
//the user defined additionas output data in the exosphery model outputdata file
void ExospherUserDefinedOutputVariableList(FILE *fout) {
  fprintf(fout, ",\"Radiation Pressure Acceleration [m/s^2]\"");
}

void ExosphereUserDefinedOutputData(FILE *fout,int spec) {
}

/*---------------------------------------------------------------------------------*/


//photoionization lifetime of sodium
/*double sodiumPhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {

  static const double LifeTime=3600.0*5.8/pow(0.4,2);

  double res,r2=x[1]*x[1]+x[2]*x[2];

  //check if the particle is outside of the Earth and lunar shadows
  if ( ((r2>rSphere*rSphere)||(x[0]<0.0)) && (Comet::EarthShadowCheck(x)==false) ) {
    res=LifeTime,PhotolyticReactionAllowedFlag=true;
  }
  else {
    res=-1.0,PhotolyticReactionAllowedFlag=false;
  }

  //check if the particle intersect the surface of the Earth
  if (pow(x[0]-Comet::xEarth_SO[0],2)+pow(x[1]-Comet::xEarth_SO[1],2)+pow(x[2]-Comet::xEarth_SO[2],2)<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) {
    res=1.0E-10*LifeTime,PhotolyticReactionAllowedFlag=true;
  }

  return res;
}

int sodiumPhotoionizationReactionProcessor(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {
  spec=_NAPLUS_SPEC_22_;

  PIC::ParticleBuffer::SetI(spec,ParticleData);
//  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;

  return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
}*/


//the codes for the soruces processes
#define _ALL_SOURCE_PROCESSES_                -1
#define _IMPACT_VAPORIZATION_SOURCE_PROCESS_   0
#define _PSD_SOURCE_PROCESS_                   1

//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
  double res=0.0;

#if _MERCURY_IMPACT_VAPORIZATION_MODE_ == _MERCURY_MODE_ON_
  if ((SourceProcessCode==-1)||(SourceProcessCode==_IMPACT_VAPORIZATION_SOURCE_PROCESS_)) { //the total impact vaporization flux
    res+=2.6E23;
  }
#endif

#if _MERCURY_PSD_MODE_ == _MERCURY_MODE_ON_
  if ((SourceProcessCode==-1)||(SourceProcessCode==_PSD_SOURCE_PROCESS_)) { //the photon stimulated desorption flux
    res+=3.6E24;
  }
#endif

  return res;
}


/*
bool sodiumDistributeParticleParameters(double *x,double *v, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,double &ParticleWeightCorrection) {

  int const nIntervalsLoockupTable=10000;

  //determine throught which process the particle will be generated
  static bool initflag=false;


  static double injectionProbabilityIV_LimitMin=-1.0,injectionProbabilityIV_LimitMax=-1.0;
  static double injectionProbabilityPSD_LimitMin=-1.0,injectionProbabilityPSD_LimitMax=-1.0;

  static double EnergyDistributionPSDTable[nIntervalsLoockupTable];

  ParticleWeightCorrection=1.0;


  if (initflag==false) {
    initflag=true;

    //init the probability limits
    double upperProbabilityUpperLimit=0.0,totalPtoductionRate=sodiumTotalProductionRate();

#if _MERCURY_IMPACT_VAPORIZATION_MODE_ == _MERCURY_MODE_ON_
    injectionProbabilityIV_LimitMin=upperProbabilityUpperLimit;
    upperProbabilityUpperLimit+=sodiumTotalProductionRate(_IMPACT_VAPORIZATION_SOURCE_PROCESS_)/totalPtoductionRate;
    injectionProbabilityIV_LimitMax=upperProbabilityUpperLimit;
#endif

#if _MERCURY_PSD_MODE_ == _MERCURY_MODE_ON_
    injectionProbabilityPSD_LimitMin=upperProbabilityUpperLimit;
    upperProbabilityUpperLimit+=sodiumTotalProductionRate(_PSD_SOURCE_PROCESS_)/totalPtoductionRate;
    injectionProbabilityPSD_LimitMax=upperProbabilityUpperLimit;


    //generte the loock-up table for the generation of the energy of injected particle
    double emin=0.0,emax=0.5*PIC::MolecularData::GetMass(_NA_SPEC_)*pow(2000.0,2);

    const int iIntegratedIntervals=10.0*nIntervalsLoockupTable;
    double e,de=(emax-emin)/iIntegratedIntervals,A=0.0,dA;
    int i,iTable;

    const double beta=0.7,U=0.052*ElectronCharge;

    for (i=0;i<iIntegratedIntervals;i++) {
      e=emin+((double)i+0.5)*de;
      A+=e/pow(e+U,2.0+beta);
    }

    dA=A/nIntervalsLoockupTable;
    A=0.0;
    EnergyDistributionPSDTable[0]=emin;

    for (i=0,iTable=0;i<iIntegratedIntervals;i++) {
      e=emin+((double)i+0.5)*de;
      A+=e/pow(e+U,2.0+beta);

      if (A>dA) {
        A-=dA;
        ++iTable;

        if (iTable==nIntervalsLoockupTable) break;

        EnergyDistributionPSDTable[iTable]=e;
      }
    }

#endif
  }


  //when the process is determines: defined the position and the velocty of the newly generated particle
  double r=rnd();

  if  ((injectionProbabilityIV_LimitMin<=r)&&(r<=injectionProbabilityIV_LimitMax)) { //injection though the impact vaporization process
    static const double InjectionTemperature=2000.0;  //the temeprature of the injected sodium atoms

    int idim;
    double ExternalNormal[3];
    double vbulk[3]={0.0,0.0,0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    for (r=0.0,idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      r+=ExternalNormal[idim]*ExternalNormal[idim];
    }



    r=-sqrt(r);

    for (idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]/=r;
      x[idim]=-rSphere*ExternalNormal[idim];
    }

    node=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (node->Thread!=PIC::Mesh::mesh.ThisThread) return false;
    startNode=node;

    //generate the particle velocity
    static const double MinParticleWeightCorrectionFactor=1.0E-5;
    bool flag;

    do {
      ParticleWeightCorrection=PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,InjectionTemperature,ExternalNormal,_NA_SPEC_,_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_);

      if (ParticleWeightCorrection>MinParticleWeightCorrectionFactor) flag=true;
      else {
        double p; // <= the probability to accept the particle

        p=ParticleWeightCorrection/MinParticleWeightCorrectionFactor;

        if (p>rnd()) flag=true,ParticleWeightCorrection=MinParticleWeightCorrectionFactor;
        else flag=false;
      }
    }
    while (flag==false);
  }
  else if ((injectionProbabilityPSD_LimitMin<=r)&&(r<=injectionProbabilityPSD_LimitMax)) { //injection through the PSD process

    int idim;
    double ExternalNormal[3]={0.0,0.0,0.0},p,l;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    //1. determine the position on the sphere where the new particle is generated
    do {  //calculate the probability of choosing the positions

      do { //generate the random position on the day side
        for (r=0.0,idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=ExternalNormal[idim]*ExternalNormal[idim];
        }
      }
      while (ExternalNormal[0]>0.0);

      r=sqrt(r);

      for (idim=0;idim<DIM;idim++) {
        ExternalNormal[idim]/=r;
        x[idim]=-(rSphere+2.0*PIC::Mesh::mesh.EPS)*ExternalNormal[idim];
      }

      //r -> the random position in teh dayside of the sphere
      //p -> the probability to generate particle at the location 'r'
      //calculate the Longitude and the Latitude of the position (Lon, Lat)

      p=x[0]/rSphere;
    }
    while (p<rnd());

    node=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (node->Thread!=PIC::Mesh::mesh.ThisThread) return false;
    startNode=node;

    //2. determine the velocity  vector of the newly generated particle
    double e,speed,emin,emax;
    int eLevel;

    eLevel=(int)(rnd()*nIntervalsLoockupTable);
    emin=EnergyDistributionPSDTable[eLevel];
    emax=(eLevel<nIntervalsLoockupTable-1) ? EnergyDistributionPSDTable[eLevel+1] : EnergyDistributionPSDTable[eLevel];

    e=emin+rnd()*(emax-emin);
    speed=sqrt(2.0*e/PIC::MolecularData::GetMass(_NA_SPEC_));

    //destribute the direction of the particle velocity vector
    do {
      l=0.0,p=0.0;

      for (idim=0;idim<3;idim++) {
        v[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());

        l+=v[idim]*v[idim];
        p+=v[idim]*ExternalNormal[idim];
      }
    }
    while (p>0.0);

    for (l=speed/sqrt(l),idim=0;idim<3;idim++) v[idim]*=l;




  }
  else exit(__LINE__,__FILE__,"Error: the injection process is not determined");


  return true;
}
*/



double SurfaceResolution(CutCell::cTriangleFace* t) {
  return 1.0 ;//max(1.0,t->CharacteristicSize()*18.0)/4.0;
}

//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double res,r,l[3];
  int idim;
  double SubsolarAngle;

  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  SubsolarAngle=acos(l[0]);


  SubsolarAngle=0.0;

  res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;



  return rSphere*res;
}


double localResolution(double *x) {
  int idim;
  double lnR,res,r=0.0;

  for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

  r=sqrt(r);

  //  if (r>dxMinGlobal*rSphere) {
  if (r>dxMinGlobal*110.0) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;

//  if ((x[0]>0.0)&&(sqrt(x[1]*x[1]+x[2]*x[2])<3.0*rSphere)) res=(res<0.4) ? res : 0.4;  ///min(res,0.4);

//  return rSphere*res;
  return rSphere*res/200; 
}


/*
double localResolution(double *x) {
  int idim;
  double lnR,res,r=0.0;

  for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

  r=sqrt(r);

  if (r>dxMinGlobal*rSphere) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;

//  if ((x[0]>0.0)&&(sqrt(x[1]*x[1]+x[2]*x[2])<3.0*rSphere)) res=(res<0.4) ? res : 0.4;  ///min(res,0.4);

 return rSphere*res;
}
*/

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;

  double CharacteristicSpeed_NA=5.0E2;

//  CharacteristicSpeed*=sqrt(PIC::MolecularData::GetMass(NA)/PIC::MolecularData::GetMass(spec));

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed_NA;


}

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

  /*
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_NA_SPEC_11_);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }*/

  return false;
}




//injection of model particles through the faces of the bounding box
/*
long int  BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
//  long int nInjectedParticles;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;

  if (spec!=_NA_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_NA_SPEC_);


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

          PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_NA_SPEC_);

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
*/

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;



//  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}






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
   PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}


/*
void prePopulateSWprotons(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double x[3],v[3],aNpart,NodeMeasure=1.0;
  int nPart,idim,i,j,k;

  static const double nSW=5.0E6; //solar wind number denstiy
  static const double tempSW=8.0E4;
  static double swVel[3]={4.0E5,0.0,0.0};

  long int newParticle,nd;
  PIC::ParticleBuffer::byte *newParticleData;

  static long int nTotalGeneratedParticles=0,nTotalProcessorBlocks=0;
  static double GlobalParticleWeight=0.0,aNpartTotal=0.0,TotalDomainVolume=0.0;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->Thread==PIC::Mesh::mesh.ThisThread) {
      nTotalProcessorBlocks++;


      //place particles in the blocks
      for (idim=0;idim<DIM;idim++) NodeMeasure*=(startNode->xmax[idim]-startNode->xmin[idim]);
      aNpart=NodeMeasure*nSW/startNode->block->GetLocalParticleWeight(SW);

      aNpartTotal+=aNpart;
      TotalDomainVolume+=NodeMeasure;

      GlobalParticleWeight=startNode->block->GetLocalParticleWeight(SW);

      nPart=(int)aNpart;
      aNpart-=nPart;
      if (aNpart>rnd()) nPart++;


      //generate particles
      for (;nPart>0;nPart--) {

        //generate a random particle's position
        for (idim=0;idim<DIM;idim++) x[idim]=startNode->xmin[idim]+(startNode->xmax[idim]-startNode->xmin[idim])*rnd();
        if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<rSphere*rSphere) continue;

        nTotalGeneratedParticles++;

        PIC::Mesh::mesh.fingCellIndex(x,i,j,k,startNode);
        nd=startNode->block->getCenterNodeLocalNumber(i,j,k);

        newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->GetCenterNode(nd)->FirstCellParticle);
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

        PIC::Distribution::MaxwellianVelocityDistribution(v,swVel,tempSW,SW);

        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetX(x,newParticleData);
        PIC::ParticleBuffer::SetI(SW,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

      }


    }

  }
  else for (int nDownNode=0;nDownNode<8;nDownNode++) prePopulateSWprotons(startNode->downNode[nDownNode]);


  if (startNode==PIC::Mesh::mesh.rootTree) {
    long int *GeneratedParticle=new long int [PIC::Mesh::mesh.nTotalThreads];
    long int *GeneratedNodes=new long int [PIC::Mesh::mesh.nTotalThreads];
    double *anpart=new double [PIC::Mesh::mesh.nTotalThreads];
    double *volume=new double [PIC::Mesh::mesh.nTotalThreads];

    MPI_Gather(&nTotalGeneratedParticles,1,MPI_LONG,GeneratedParticle,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&nTotalProcessorBlocks,1,MPI_LONG,GeneratedNodes,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&aNpartTotal,1,MPI_DOUBLE,anpart,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&TotalDomainVolume,1,MPI_DOUBLE,volume,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);


    if (PIC::Mesh::mesh.ThisThread==0) {
      cout << "Pre Population of the domain:\n Thread\tGenerated Particles\tDomain Block Number\t aNpartTotal\tSubDomain Volume" << endl;

      for (int t=0;t<PIC::Mesh::mesh.nTotalThreads;t++) cout << t << "\t" << GeneratedParticle[t] << "\t" << GeneratedNodes[t] << "\t" << anpart[t] << "\t" << volume[t] << endl;

      cout << "Global Particle's weight=" << GlobalParticleWeight << endl;
    }

    delete [] GeneratedNodes;
    delete [] GeneratedParticle;
    delete [] anpart;
    delete [] volume;
  }

}
*/

/*double sphereInjectionRate(int spec,void *SphereDataPointer) {
  double res=0.0;

  if (spec==_NA_SPEC_11_) res=sodiumTotalProductionRate();
  else if (spec==_NAPLUS_SPEC_11_) res=0.0;
  else exit(__LINE__,__FILE__,"Error: the source rate for the species is not determined");


  return res;
}*/

/*
long int sphereParticleInjection(void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ParticleWeight,LocalTimeStep,x[3],v[3],*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
//  int idim;

  double ParticleWeightCorrection;

//  static const double Temp=200.0;
//  double vbulk[3]={0.0,0.0,0.0};


//  return 0;

//====================  DEBUG ===========================
//  static bool FirstPArticleGenerated=false;
//====================  END DEBUG ===================================


  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_NA_SPEC_];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[_NA_SPEC_];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif



  static double InjectParticles=0.0;

  bool flag;

  InjectParticles+=sphereInjectionRate(_NA_SPEC_,SphereDataPointer)*LocalTimeStep;

  while (InjectParticles>0.0) {



    startNode=NULL;


    //if (sodiumDistributeParticleParameters(x,v,startNode,ParticleWeightCorrection)==false) continue;


    flag=sodiumDistributeParticleParameters(x,v,startNode,ParticleWeightCorrection);
    InjectParticles-=ParticleWeight*ParticleWeightCorrection;
    if (flag==false) continue;

//====================  DEBUG ===========================
    {
static double InjectionRadialVelocity=0.0,InjectionTangentionalSpeed=0.0;
static long int nTotalInjectedParticles=0;

double l[3],r=0.0,v0=0.0,v1=0.0;
int idim;

for (idim=0;idim<3;idim++) r+=pow(x[idim],2);
r=sqrt(r);
for (idim=0;idim<3;idim++) {
  l[idim]=x[idim]/r;

  v0+=v[idim]*l[idim];
}

for (idim=0;idim<3;idim++) v1+=pow(v[idim]-v0*l[idim],2);

nTotalInjectedParticles++;
InjectionRadialVelocity+=v0;
InjectionTangentionalSpeed+=sqrt(v1);
    }
//====================  END DEBUG ===================================




    if (startNode->block->GetLocalTimeStep(_NA_SPEC_)/LocalTimeStep<rnd()) continue;

    //generate a particle
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    nInjectedParticles++;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(v,newParticleData);
    PIC::ParticleBuffer::SetI(_NA_SPEC_,newParticleData);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);


    //inject the particle into the system
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(_NA_SPEC_)*rnd() LocalTimeStep-TimeCounter,startNode,true);
//    PIC::Mover::MoveParticleBoundaryInjection[NA](newParticle,startNode->block->GetLocalTimeStep(NA)*rnd() / *LocalTimeStep-TimeCounter* /,startNode,true);
  }

  return nInjectedParticles;
}
*/


bool Radius(double &r,double x) {


 // r=1.0;
 // return true;

  /*
  if ((x>=-Pi-1.0E-15)&&(x<=Pi+1.0E-15)) {
    r=fabs(sin(x));
    return true;
    }*/

  double boundary;
  
    if ((x>=-1165.0-1.0e-15) && (x<=1165.0+1.0e-15)) {
    boundary=sqrt(1980.0*1980.0-x*x*1980.0*2.0/2330.0*1980.0*2.0/2330.0)/1980.0*2330.0/(2.0*1.3)-(2330.0/2.0-600.0)*0.5*(1+cos((-x*1980.0*2.0/2330.0+200.0)/1800.0*Pi));
    
    r=((boundary>=0)? sqrt(boundary*boundary):0.0);
    return true;
    }
    /*
  if ((x>=-rSphere-1.0e-15) && (x<=rSphere+1.0e-15)) {
    r=sqrt(rSphere*rSphere-x*x);  
    return true;
    }*/

  return false;
}

void testPrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalRotationBodyData *Nucleus,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
 int nInterpolationElement,nSurfaceElement;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double t,TotalFluxDown=0.0,TotalFluxUp=0.0,SurfaceContent=0.0,BulkSpeedDown=0.0,BulkSpeedUp=0.0,SampleSpeciesSurfaceInjectionFlux=0.0;


#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
   double FluxIV=0.0;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
   double FluxPDS=0.0;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   double FluxTD=0.0;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   double FluxSW=0.0,SolarWindIncidentFlux=0.0;
#endif



#if _SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  exit(__LINE__,__FILE__"Error: the model is implemeted only for _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_");
#endif

  /*
  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
    InterpolationCoefficient=Nucleus->GetSurfaceElementArea(nSurfaceElement);

    BulkSpeedDown+=Sphere->SampleReturnFluxBulkSpeed[spec][nSurfaceElement];
    TotalFluxDown+=Sphere->SampleSpeciesSurfaceReturnFlux[spec][nSurfaceElement];
    SurfaceContent+=Sphere->SampleSpeciesSurfaceAreaDensity[spec][nSurfaceElement]*InterpolationCoefficient;

    BulkSpeedUp+=Sphere->SampleInjectedFluxBulkSpeed[spec][nSurfaceElement];
    SampleSpeciesSurfaceInjectionFlux+=Sphere->SampleSpeciesSurfaceInjectionFlux[spec][nSurfaceElement];


#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_]/PIC::LastSampleLength : 0.0;
    FluxIV+=t;
    TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__PHOTON_STIMULATED_DESPRPTION_]/PIC::LastSampleLength : 0.0;
    FluxPDS+=t;
    TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__THERMAL_DESORPTION_]/PIC::LastSampleLength : 0.0;
   FluxTD+=t;
   TotalFluxUp+=t;
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
   t=(PIC::LastSampleLength!=0) ? Sphere->SampleSpeciesSurfaceSourceRate[spec][nSurfaceElement][_EXOSPHERE_SOURCE__ID__SOLAR_WIND_SPUTTERING_]/PIC::LastSampleLength : 0.0;
   FluxSW+=t;
   TotalFluxUp+=t;

   SolarWindIncidentFlux+=Sphere->SolarWindSurfaceFlux[nSurfaceElement]*InterpolationCoefficient;
#endif

    InterpolationNormalization+=InterpolationCoefficient;
  }
  */


  if (ThisThread==0)  {
    /* //collect sampled data from all processors
    for (int thread=1;thread<nTotalThreads;thread++) {
      TotalFluxDown+=pipe->recv<double>(thread);
//      SurfaceContent+=pipe->recv<double>(thread);  All processors have the same distribution of surface content map
      TotalFluxUp+=pipe->recv<double>(thread);
      BulkSpeedDown+=pipe->recv<double>(thread);

      BulkSpeedUp+=pipe->recv<double>(thread);
      SampleSpeciesSurfaceInjectionFlux+=pipe->recv<double>(thread);

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
      FluxIV+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
      FluxPDS+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      FluxTD+=pipe->recv<double>(thread);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      FluxSW+=pipe->recv<double>(thread);
#endif
    }

    if (PIC::LastSampleLength!=0) {
      if (TotalFluxDown>0.0) BulkSpeedDown/=TotalFluxDown;
      if (SampleSpeciesSurfaceInjectionFlux>0.0) BulkSpeedUp/=SampleSpeciesSurfaceInjectionFlux;

      TotalFluxDown/=PIC::LastSampleLength*PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
      SurfaceContent/=PIC::LastSampleLength;
    }
  */

    //Print Surface Temparature
    double norm[3],CosSubSolarAngle,SurfaceTemperature=0.0,x_LOCAL_IAU_OBJECT[3],x_LOCAL_SO_OBJECT[3];
    int i;
    double totalInterpolationCoefficient=0.0;
    for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
      i=SurfaceElementsInterpolationList[nInterpolationElement];
      InterpolationCoefficient=Nucleus->GetSurfaceElementArea(i);
      Nucleus->GetSurfaceElementMiddlePoint(x_LOCAL_IAU_OBJECT,i);
      Nucleus->GetSurfaceElementNormal(norm,i);	
      
      double c,X;
      int idim;
      
      double positionSun[3],HeliocentricDistance=1.064*_AU_,subSolarPointAzimuth=53.0*Pi/180;
      
      for (idim=0;idim<3;idim++) x_LOCAL_SO_OBJECT[idim]=x_LOCAL_IAU_OBJECT[idim];
      
      positionSun[0]=HeliocentricDistance*cos(subSolarPointAzimuth);
      positionSun[1]=HeliocentricDistance*sin(subSolarPointAzimuth);
      positionSun[2]=0.0;
      
      for (c=0.0,X=0.0,idim=0;idim<3;idim++){
	c+=norm[idim]*(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim]);
	X+=pow(positionSun[idim]-x_LOCAL_IAU_OBJECT[idim],2.0);
      }
      CosSubSolarAngle=c/sqrt(X);
            
      SurfaceTemperature+=InterpolationCoefficient*Comet::GetSurfaceTemeprature(CosSubSolarAngle,x_LOCAL_SO_OBJECT);
      totalInterpolationCoefficient+=InterpolationCoefficient;
    }
    fprintf(fout," %e",SurfaceTemperature/totalInterpolationCoefficient);

    //Print Sampled Surface data
    double ReemissionParticleFraction;
    
    //    fprintf(fout," %e %e %e %e %e %e ",TotalFluxDown/InterpolationNormalization,TotalFluxUp/InterpolationNormalization,SurfaceContent/InterpolationNormalization,BulkSpeedDown,BulkSpeedUp, \
    //    	    Comet::SurfaceInteraction::StickingProbability(spec,ReemissionParticleFraction,SurfaceTemperature));
    fprintf(fout," %e %e %e %e %e %e ",0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    fprintf(fout," %e ",FluxIV/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    fprintf(fout," %e ",FluxPDS/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
      fprintf(fout," %e ",FluxTD/InterpolationNormalization);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
      fprintf(fout," %e %e ",FluxSW/InterpolationNormalization,SolarWindIncidentFlux/InterpolationNormalization);
#endif

}
  else {
    /*    pipe->send(TotalFluxDown);
//    pipe->send(SurfaceContent);     All processors have the same distribution of surface content map
    pipe->send(TotalFluxUp);
    pipe->send(BulkSpeedDown);

    pipe->send(BulkSpeedUp);
    pipe->send(SampleSpeciesSurfaceInjectionFlux);

#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxIV);
#endif

#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxPDS);
#endif

#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxTD);
#endif

#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
    pipe->send(FluxSW);
#endif
    */ }
}


void amps_init() {
  //  MPI_Init(&argc,&argv);
    PIC::InitMPI();


    //SetUp the alarm
  //  PIC::Alarm::SetAlarm(2000);



  //  VT_OFF();
  //  VT_traceoff();


  #ifdef _MAIN_PROFILE_
    initSaturn ("");
  #endif


    rnd_seed();



    //char inputFile[]="mercury.input";



    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


    //init the Mercury model
    Comet::Init_BeforeParser();

    //init the particle solver
    PIC::Init_BeforeParser();
  //  PIC::Parser::Run(inputFile);



    Comet::Sampling::SetUserDefinedAdditionalOutputSampledModelDataFunctions(ExospherUserDefinedOutputVariableList,ExosphereUserDefinedOutputData);


    PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=500; //0; //00; //*10;
    PIC::RequiredSampleLength=10; //00; //0; //0;

    Comet::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;

    Comet::Init_AfterParser();


    //output the PDS enerfy distribution function
    if (PIC::ThisThread==0) {
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        char fname[200];

        sprintf(fname,"%s/CumulativeEnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
        Comet::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintCumulativeDistributionFunction(fname);

        sprintf(fname,"%s/EnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
        Comet::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintDistributionFunction(fname,&spec);
      }

  //    cout << Comet::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
    }



    //register the sphere
    if (false) {
      double sx0[3]={0.0,0.0,0.0};
      cInternalBoundaryConditionsDescriptor SphereDescriptor;
      cInternalSphericalData *Sphere;


      //reserve memory for sampling of the surface balance of sticking species
      long int ReserveSamplingSpace[PIC::nTotalSpecies];

      for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


      cInternalSphericalData::SetGeneralSurfaceMeshParameters(50,50);



      PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
      SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
      Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
      Sphere->SetSphereGeometricalParameters(sx0,rSphere);


      //init the object for distribution of the injection surface elements
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        Comet::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


      #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
        Comet::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
      #endif

      #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
        Comet::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
      #endif
      }


      char fname[_MAX_STRING_LENGTH_PIC_];

      sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
      Sphere->PrintSurfaceMesh(fname);

      sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
      Sphere->PrintSurfaceData(fname,0);

      Sphere->localResolution=localSphericalSurfaceResolution;
      Sphere->InjectionRate=Comet::SourceProcesses::totalProductionRate;
      Sphere->faceat=0;
      Sphere->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
      Sphere->InjectionBoundaryCondition=Comet::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;

      Sphere->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
      Sphere->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;
      Sphere->PrintDataStateVector=Comet::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

      //set up the planet pointer in Mercury model
      Comet::Planet=Sphere;

      //allocate the buffers for collecting the sodium surface density
      Sphere->SurfaceElementDesorptionFluxUP=new double*[PIC::nTotalSpecies];
      Sphere->SurfaceElementAdsorptionFluxDOWN=new double*[PIC::nTotalSpecies];
      Sphere->SurfaceElementPopulation=new double*[PIC::nTotalSpecies];

      Sphere->SurfaceElementDesorptionFluxUP[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
      Sphere->SurfaceElementAdsorptionFluxDOWN[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
      Sphere->SurfaceElementPopulation[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
        Sphere->SurfaceElementDesorptionFluxUP[spec]=Sphere->SurfaceElementDesorptionFluxUP[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
        Sphere->SurfaceElementAdsorptionFluxDOWN[spec]=Sphere->SurfaceElementAdsorptionFluxDOWN[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
        Sphere->SurfaceElementPopulation[spec]=Sphere->SurfaceElementPopulation[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
      }

      Sphere->SurfaceElementExternalNormal=new cInternalSphericalData_UserDefined::cSurfaceElementExternalNormal[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
      Sphere->SurfaceElementArea=new double[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];


      Sphere->ElementSourceRate=new cInternalSphericalData_UserDefined::cElementSourceRate*[PIC::nTotalSpecies];
      Sphere->ElementSourceRate[0]=new cInternalSphericalData_UserDefined::cElementSourceRate[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
        Sphere->ElementSourceRate[spec]=Sphere->ElementSourceRate[spec-1]+PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
      }


      Sphere->SolarWindSurfaceFlux=new double[PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      for (int el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
          Sphere->SurfaceElementDesorptionFluxUP[spec][el]=0.0;
          Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el]=0.0;
          Sphere->SurfaceElementPopulation[spec][el]=0.0;
        }

        Sphere->SolarWindSurfaceFlux[el]=-1.0;

        Sphere->SurfaceElementArea[el]=Sphere->GetSurfaceElementArea(el);
        Sphere->GetSurfaceElementNormal((Sphere->SurfaceElementExternalNormal+el)->norm,el);
      }

      //allocate buffers for sampling surface sodium source rates and sodikum surface content
      int offsetSpecie,offsetElement,s,el,i;

      Sphere->SampleSpeciesSurfaceSourceRate=new double** [PIC::nTotalSpecies];
      Sphere->SampleSpeciesSurfaceSourceRate[0]=new double *[PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];
      Sphere->SampleSpeciesSurfaceSourceRate[0][0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];

      for (offsetSpecie=0,s=0,offsetElement=0;s<PIC::nTotalSpecies;s++) {
        Sphere->SampleSpeciesSurfaceSourceRate[s]=Sphere->SampleSpeciesSurfaceSourceRate[0]+offsetSpecie;
        offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
          Sphere->SampleSpeciesSurfaceSourceRate[s][el]=Sphere->SampleSpeciesSurfaceSourceRate[0][0]+offsetElement;
          offsetElement+=_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;

          for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) Sphere->SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;
        }
      }

      Sphere->SampleSpeciesSurfaceAreaDensity=new double* [PIC::nTotalSpecies];
      Sphere->SampleSpeciesSurfaceAreaDensity[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
        Sphere->SampleSpeciesSurfaceAreaDensity[s]=Sphere->SampleSpeciesSurfaceAreaDensity[0]+offsetSpecie;
        offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
          Sphere->SampleSpeciesSurfaceAreaDensity[s][el]=0.0;
        }
      }

      Sphere->SampleSpeciesSurfaceReturnFlux=new double* [PIC::nTotalSpecies];
      Sphere->SampleSpeciesSurfaceReturnFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      Sphere->SampleSpeciesSurfaceInjectionFlux=new double* [PIC::nTotalSpecies];
      Sphere->SampleSpeciesSurfaceInjectionFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      Sphere->SampleReturnFluxBulkSpeed=new double* [PIC::nTotalSpecies];
      Sphere->SampleReturnFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      Sphere->SampleInjectedFluxBulkSpeed=new double* [PIC::nTotalSpecies];
      Sphere->SampleInjectedFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber];

      for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
        Sphere->SampleSpeciesSurfaceReturnFlux[s]=Sphere->SampleSpeciesSurfaceReturnFlux[0]+offsetSpecie;
        Sphere->SampleSpeciesSurfaceInjectionFlux[s]=Sphere->SampleSpeciesSurfaceInjectionFlux[0]+offsetSpecie;

        Sphere->SampleReturnFluxBulkSpeed[s]=Sphere->SampleReturnFluxBulkSpeed[0]+offsetSpecie;
        Sphere->SampleInjectedFluxBulkSpeed[s]=Sphere->SampleInjectedFluxBulkSpeed[0]+offsetSpecie;

        offsetSpecie+=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
          Sphere->SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
          Sphere->SampleSpeciesSurfaceInjectionFlux[s][el]=0.0;

          Sphere->SampleReturnFluxBulkSpeed[s][el]=0.0;
          Sphere->SampleInjectedFluxBulkSpeed[s][el]=0.0;
        }
      }

    }


    //+++++++++++++++++++++++++++++  NUCLEUS +++++++++++++++++++++++++++++++++
    //register the nucleus (the body of rotation)
    if (true) {
      double sx0[3]={0.0,0.0,0.0};
      cInternalBoundaryConditionsDescriptor RotationBodyDescriptor;
      cInternalRotationBodyData *Nucleus;


      //reserve memory for sampling of the surface balance of sticking species
      long int ReserveSamplingSpace[PIC::nTotalSpecies];

      for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


//      cInternalRotationBodyData::SetGeneralSurfaceMeshParameters(60,100);


      double x0[3]={0.0,0.0,0.0},l[3]={1.0,0.0,0.0};
      //cInternalRotationBodyData::SetGeneralSurfaceMeshParameters(60,100);


      PIC::BC::InternalBoundary::RotationBody::Init(ReserveSamplingSpace,NULL);
      RotationBodyDescriptor=PIC::BC::InternalBoundary::RotationBody::RegisterInternalRotationBody();
      Nucleus=(cInternalRotationBodyData*) RotationBodyDescriptor.BoundaryElement;

      Nucleus->SetGeneralSurfaceMeshParameters(60,100);
      //Nucleus->SetGeometricalParameters(x0,l,-rSphere,rSphere,Radius);
      Nucleus->SetGeometricalParameters(x0,l,-1165.0,1165.0,Radius);


      //init the object for distribution of the injection surface elements
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        Comet::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


      #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
        Comet::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
      #endif

      #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
        Comet::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber-1,
            Comet::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
      #endif
      }


      char fname[_MAX_STRING_LENGTH_PIC_];

      Nucleus->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
      Nucleus->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;
      Nucleus->PrintDataStateVector=testPrintDataStateVector; ///Comet::Sampling::OutputSurfaceDataFile::PrintDataStateVector;


      sprintf(fname,"%s/RotationBodyNucleus.dat",PIC::OutputDataFileDirectory);
      Nucleus->PrintSurfaceMesh(fname);

      //sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
      //Nucleus->PrintSurfaceData(fname,0);

      Nucleus->localResolution=localSphericalSurfaceResolution;
      Nucleus->InjectionRate=Comet::SourceProcesses::totalProductionRate;
      Nucleus->faceat=0;
      Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
      //Nucleus->InjectionBoundaryCondition=Comet::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
      //      Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
      PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;

      //set up the planet pointer in Mercury model


      //allocate the buffers for collecting the sodium surface density
      Nucleus->SurfaceElementDesorptionFluxUP=new double*[PIC::nTotalSpecies];
      Nucleus->SurfaceElementAdsorptionFluxDOWN=new double*[PIC::nTotalSpecies];
      Nucleus->SurfaceElementPopulation=new double*[PIC::nTotalSpecies];

      Nucleus->SurfaceElementDesorptionFluxUP[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];
      Nucleus->SurfaceElementAdsorptionFluxDOWN[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];
      Nucleus->SurfaceElementPopulation[0]=new double[PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
        Nucleus->SurfaceElementDesorptionFluxUP[spec]=Nucleus->SurfaceElementDesorptionFluxUP[spec-1]+PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;
        Nucleus->SurfaceElementAdsorptionFluxDOWN[spec]=Nucleus->SurfaceElementAdsorptionFluxDOWN[spec-1]+PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;
        Nucleus->SurfaceElementPopulation[spec]=Nucleus->SurfaceElementPopulation[spec-1]+PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;
      }

      Nucleus->SurfaceElementExternalNormal=new cInternalSphericalData_UserDefined::cSurfaceElementExternalNormal[PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];
      Nucleus->SurfaceElementArea=new double[PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];


      Nucleus->ElementSourceRate=new cInternalSphericalData_UserDefined::cElementSourceRate*[PIC::nTotalSpecies];
      Nucleus->ElementSourceRate[0]=new cInternalSphericalData_UserDefined::cElementSourceRate[PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
        Nucleus->ElementSourceRate[spec]=Nucleus->ElementSourceRate[spec-1]+PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;
      }


      Nucleus->SolarWindSurfaceFlux=new double[PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      for (int el=0;el<PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;el++) {
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
          Nucleus->SurfaceElementDesorptionFluxUP[spec][el]=0.0;
          Nucleus->SurfaceElementAdsorptionFluxDOWN[spec][el]=0.0;
          Nucleus->SurfaceElementPopulation[spec][el]=0.0;
        }

        Nucleus->SolarWindSurfaceFlux[el]=-1.0;

        Nucleus->SurfaceElementArea[el]=Nucleus->GetSurfaceElementArea(el);
        Nucleus->GetSurfaceElementNormal((Nucleus->SurfaceElementExternalNormal+el)->norm,el);
      }

      //allocate buffers for sampling surface sodium source rates and sodikum surface content
      int offsetSpecie,offsetElement,s,el,i;

      Nucleus->SampleSpeciesSurfaceSourceRate=new double** [PIC::nTotalSpecies];
      Nucleus->SampleSpeciesSurfaceSourceRate[0]=new double *[PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];
      Nucleus->SampleSpeciesSurfaceSourceRate[0][0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber*(_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1)];

      for (offsetSpecie=0,s=0,offsetElement=0;s<PIC::nTotalSpecies;s++) {
        Nucleus->SampleSpeciesSurfaceSourceRate[s]=Nucleus->SampleSpeciesSurfaceSourceRate[0]+offsetSpecie;
        offsetSpecie+=PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;el++) {
          Nucleus->SampleSpeciesSurfaceSourceRate[s][el]=Nucleus->SampleSpeciesSurfaceSourceRate[0][0]+offsetElement;
          offsetElement+=_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;

          for (i=0;i<_EXOSPHERE__SOURCE_MAX_ID_VALUE_+1;i++) Nucleus->SampleSpeciesSurfaceSourceRate[s][el][i]=0.0;
        }
      }

      Nucleus->SampleSpeciesSurfaceAreaDensity=new double* [PIC::nTotalSpecies];
      Nucleus->SampleSpeciesSurfaceAreaDensity[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
        Nucleus->SampleSpeciesSurfaceAreaDensity[s]=Nucleus->SampleSpeciesSurfaceAreaDensity[0]+offsetSpecie;
        offsetSpecie+=PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;el++) {
          Nucleus->SampleSpeciesSurfaceAreaDensity[s][el]=0.0;
        }
      }

      Nucleus->SampleSpeciesSurfaceReturnFlux=new double* [PIC::nTotalSpecies];
      Nucleus->SampleSpeciesSurfaceReturnFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      Nucleus->SampleSpeciesSurfaceInjectionFlux=new double* [PIC::nTotalSpecies];
      Nucleus->SampleSpeciesSurfaceInjectionFlux[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      Nucleus->SampleReturnFluxBulkSpeed=new double* [PIC::nTotalSpecies];
      Nucleus->SampleReturnFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      Nucleus->SampleInjectedFluxBulkSpeed=new double* [PIC::nTotalSpecies];
      Nucleus->SampleInjectedFluxBulkSpeed[0]=new double [PIC::nTotalSpecies*PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber];

      for (offsetSpecie=0,s=0;s<PIC::nTotalSpecies;s++) {
        Nucleus->SampleSpeciesSurfaceReturnFlux[s]=Nucleus->SampleSpeciesSurfaceReturnFlux[0]+offsetSpecie;
        Nucleus->SampleSpeciesSurfaceInjectionFlux[s]=Nucleus->SampleSpeciesSurfaceInjectionFlux[0]+offsetSpecie;

        Nucleus->SampleReturnFluxBulkSpeed[s]=Nucleus->SampleReturnFluxBulkSpeed[0]+offsetSpecie;
        Nucleus->SampleInjectedFluxBulkSpeed[s]=Nucleus->SampleInjectedFluxBulkSpeed[0]+offsetSpecie;

        offsetSpecie+=PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;

        for (el=0;el<PIC::BC::InternalBoundary::RotationBody::TotalSurfaceElementNumber;el++) {
          Nucleus->SampleSpeciesSurfaceReturnFlux[s][el]=0.0;
          Nucleus->SampleSpeciesSurfaceInjectionFlux[s][el]=0.0;

          Nucleus->SampleReturnFluxBulkSpeed[s][el]=0.0;
          Nucleus->SampleInjectedFluxBulkSpeed[s][el]=0.0;
        }
      }

    }

  //=============   DEBUG ==============

    /*
    double xProbeMin[3]={-2172234,-1.236913E-10,2000742};
    double dx=-2172234-(-2229398),v=1.0,t;
    double xProbeMax[3];
    int IntersectionStatus;

    for (int i=0;i<3;i++) {
      xProbeMax[i]=xProbeMin[i]+dx;
      v*=dx;
    }

    t=Sphere->GetRemainedBlockVolume(xProbeMin,xProbeMax,0.0,&IntersectionStatus);
  */

  //============= END DEBUG =============

    //init the solver
    PIC::Mesh::initCellSamplingDataBuffer();

    //init the mesh
  //  cout << "Init the mesh" << endl;

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

        fprintf(PIC::DiagnospticMessageStream,"Size of the domain is smaller that the radius of the body: the size of the domain is readjusted\n");
        r=2.01*yMaxDomain*_RADIUS_(_TARGET_)/DomainLength[1];

        for (idim=0;idim<DIM;idim++) DomainLength[idim]*=r;
      }

      DomainCenterOffset[0]=-yMaxDomain*rSphere;////-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
      DomainCenterOffset[1]=-DomainLength[1]/2.0;
      DomainCenterOffset[2]=-DomainLength[2]/2.0;
    }

    /*
    for (idim=0;idim<DIM;idim++) {
      xmin[idim]=DomainCenterOffset[idim];
      xmax[idim]=DomainLength[idim]+DomainCenterOffset[idim];
    }
    */

    for (idim=0;idim<DIM;idim++) {
      xmax[idim]=-DomainCenterOffset[idim];
      xmin[idim]=-(DomainLength[idim]+DomainCenterOffset[idim]);
    }

  //  double xmin[3]={-xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};
  //  double xmax[3]={xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};


    //generate only the tree
    //    PIC::Mesh::mesh.CutCellSurfaceLocalResolution=SurfaceResolution;
    PIC::Mesh::mesh.AllowBlockAllocation=false;
    PIC::Mesh::mesh.init(xmin,xmax,localSphericalSurfaceResolution);

  //  VT_ON();
  //  VT_USER_START("name");
  //  VT_ON();

  //  {
  //    VT_TRACER("name");

   char fname[_MAX_STRING_LENGTH_PIC_];

   /**/
   sprintf(fname,"%s/mesh.msh",PIC::OutputDataFileDirectory);
    if (PIC::Mesh::mesh.ThisThread==0) {
      PIC::Mesh::mesh.buildMesh();
      PIC::Mesh::mesh.saveMeshFile(fname);
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Mesh::mesh.readMeshFile(fname);
    }
   
    sprintf(fname,"%s/mesh.dat",PIC::OutputDataFileDirectory);
    PIC::Mesh::mesh.outputMeshTECPLOT(fname);
   

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

  //init the PIC solver
    PIC::Init_AfterParser ();
    PIC::Mover::Init();
  //  PIC::Mover::TotalParticleAcceleration=TotalParticleAcceleration;

  //  for (int s=0;s<PIC::nTotalSpecies;s++) {
  //    PIC::Mover::MoveParticleTimeStep[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder; ///UniformWeight_UniformTimeStep_SecondOrder;
  //    PIC::Mover::MoveParticleBoundaryInjection[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder;
  //  }



    //set up the time step
    PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
    PIC::ParticleWeightTimeStep::initTimeStep();


    PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
    PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H2O_SPEC_);
    PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_CO2_SPEC_);
    


    //set photolytic reactions
  /*
    PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(sodiumPhotoionizationReactionProcessor,_NA_SPEC_11_);
    PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(sodiumPhotoionizationLifeTime,_NA_SPEC_11_);
  */


  //  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
  //  PIC::Mesh::mesh.outputMeshDataTECPLOT("mesh.data.dat",0);

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  //  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;




    //output final data
  //  PIC::Mesh::mesh.outputMeshDataTECPLOT("final.data.dat",0);

    //create the list of mesh nodes where the injection boundary conditinos are applied
    PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
    PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
    PIC::BC::InitBoundingBoxInjectionBlockList();


    //init the particle buffer
    PIC::ParticleBuffer::Init(1000000);
  //  double TimeCounter=time(NULL);
  //  int LastDataOutputFileNumber=-1;


    //init the sampling of the particls' distribution functions
    /*    const int nSamplePoints=3;
    double SampleLocations[nSamplePoints][DIM]={{3.0E3,3.0E3,0.0}, {2.8E3,5.6E3,0.0}, {-2.3E3,3.0E3,0.0}};

    PIC::DistributionFunctionSample::vMin=-40.0E3;
    PIC::DistributionFunctionSample::vMax=40.0E3;
    PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

    PIC::DistributionFunctionSample::Init();
    */

    //time step
    double SimulationTimeStep=-1.0;
    
    //    for (int spec=0;spec<PIC::nTotalSpecies;spec++) if ((SimulationTimeStep<0.0)||(SimulationTimeStep>PIC::ParticleWeightTimeStep::GlobalTimeStep[spec])) {
    //SimulationTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
      // }

    //set up the time step                                                                                    
    PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
    PIC::ParticleWeightTimeStep::initTimeStep();

    //    SimulationTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_];
    //SimulationTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[_CO2_SPEC_];
  }

void amps_time_step() {



    //make the time advance
     PIC::TimeStep();

     // cout << __LINE__ << __FILE__ << endl;
}


