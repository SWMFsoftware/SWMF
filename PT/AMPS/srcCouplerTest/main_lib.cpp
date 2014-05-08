//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

//the particle class
#include "pic.h"
#include "constants.h"

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

#define _TARGET_ _EARTH_

const double rSphere=_RADIUS_(_TARGET_);


const double xMaxDomain=5.0; //modeling of the tail
//const double xMaxDomain=20; //modeling the vicinity of the planet
const double yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun



/*const double xMaxDomain=4.0; //modeling of the tail
//const double xMaxDomain=20; //modeling the vicinity of the planet
const double yMaxDomain=4.0; //the minimum size of the domain in the direction perpendicular to the direction to the sun
*/


const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100,dxMaxSphere=DebugRunMultiplier*2.0/10.0;



//the species
/*int NA=0;
int NAPLUS=1;*/


//the total acceleration of the particles
#include "Na.h"



/*---------------------------------------------------------------------------------*/
//the user defined additionas output data in the exosphery model outputdata file
void ExospherUserDefinedOutputVariableList(FILE *fout) {
  fprintf(fout, ",\"Radiation Pressure Acceleration [m/s^2]\"");
}

void ExosphereUserDefinedOutputData(FILE *fout,int spec) {

  if (PIC::ThisThread==0) {
    double t;

    t=(spec==_NA_SPEC_) ? SodiumRadiationPressureAcceleration__Combi_1997_icarus(Moon::vObjectRadial,Moon::xObjectRadial) : 0.0;
    fprintf(fout,"%e  ",t);
  }
}

/*---------------------------------------------------------------------------------*/


//photoionization lifetime of sodium
/*double sodiumPhotoionizationLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag) {

  static const double LifeTime=3600.0*5.8/pow(0.4,2);

  double res,r2=x[1]*x[1]+x[2]*x[2];

  //check if the particle is outside of the Earth and lunar shadows
  if ( ((r2>rSphere*rSphere)||(x[0]<0.0)) && (Moon::EarthShadowCheck(x)==false) ) {
    res=LifeTime,PhotolyticReactionAllowedFlag=true;
  }
  else {
    res=-1.0,PhotolyticReactionAllowedFlag=false;
  }

  //check if the particle intersect the surface of the Earth
  if (pow(x[0]-Moon::xEarth_SO[0],2)+pow(x[1]-Moon::xEarth_SO[1],2)+pow(x[2]-Moon::xEarth_SO[2],2)<_RADIUS_(_EARTH_)*_RADIUS_(_EARTH_)) {
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



  //return rSphere*res;
  return 1e9;
}

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


//  return rSphere*res;

  //return max(pow(r/rSphere,0.25)*5e8,1e7);

  return max(r*400,1e7);

//  return 5e7;
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;

  double CharacteristicSpeed_NA=2.0E3;

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
//    Moon::Init_BeforeParser();

    //init the particle solver
    PIC::Init_BeforeParser();
  //  PIC::Parser::Run(inputFile);



    Moon::Sampling::SetUserDefinedAdditionalOutputSampledModelDataFunctions(ExospherUserDefinedOutputVariableList,ExosphereUserDefinedOutputData);


    const int InitialSampleLength=100;
    PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=20; //0; //00; //*10;
    PIC::RequiredSampleLength=InitialSampleLength; //00; //0; //0;

    Moon::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;

//    Moon::Init_AfterParser();


    //output the PDS enerfy distribution function
    if (PIC::ThisThread==0) {
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        char fname[200];

        sprintf(fname,"%s/CumulativeEnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
        Moon::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintCumulativeDistributionFunction(fname);

        sprintf(fname,"%s/EnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
        Moon::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintDistributionFunction(fname,&spec);
      }

  //    cout << Moon::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
    }



    //register the sphere
    if (false) {
      double sx0[3]={0.0,0.0,0.0};
      cInternalBoundaryConditionsDescriptor SphereDescriptor;
      cInternalSphericalData *Sphere;


      //reserve memory for sampling of the surface balance of sticking species
      long int ReserveSamplingSpace[PIC::nTotalSpecies];

      for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


      cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);



      PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
      SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
      Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
      Sphere->SetSphereGeometricalParameters(sx0,rSphere);


      //init the object for distribution of the injection surface elements
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
        Moon::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Moon::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


      #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
        Moon::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Moon::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
      #endif

      #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
        Moon::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
            Moon::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
      #endif
      }


      char fname[_MAX_STRING_LENGTH_PIC_];

      sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
      Sphere->PrintSurfaceMesh(fname);

      sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
      Sphere->PrintSurfaceData(fname,0);

      Sphere->localResolution=localSphericalSurfaceResolution;
      Sphere->InjectionRate=Moon::SourceProcesses::totalProductionRate;
      Sphere->faceat=0;
      Sphere->ParticleSphereInteraction=Moon::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
      Sphere->InjectionBoundaryCondition=Moon::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;

      Sphere->PrintTitle=Moon::Sampling::OutputSurfaceDataFile::PrintTitle;
      Sphere->PrintVariableList=Moon::Sampling::OutputSurfaceDataFile::PrintVariableList;
      Sphere->PrintDataStateVector=Moon::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

      //set up the planet pointer in Mercury model
      Moon::Planet=Sphere;

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


//the domain is the cube (-100,100)
    

    std::cout<<"Setting grid bounds"<<std::endl;
    for (idim=0;idim<DIM;idim++) xmin[idim]=-16.0*rSphere,xmax[idim]=16.0*rSphere;

    std::cout<<"xmin: "<<xmin[0]<<","<<xmin[1]<<","<<xmin[2]<<std::endl;
    std::cout<<"xmax: "<<xmax[0]<<","<<xmax[1]<<","<<xmax[2]<<std::endl;

    //generate only the tree
    PIC::Mesh::mesh.AllowBlockAllocation=false;
    std::cout<<"Initializing mesh"<<std::endl;
    PIC::Mesh::mesh.init(xmin,xmax,localResolution);
    PIC::Mesh::mesh.memoryAllocationReport();


  //  VT_ON();
  //  VT_USER_START("name");
  //  VT_ON();

  //  {
  //    VT_TRACER("name");

   char fname[_MAX_STRING_LENGTH_PIC_];


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

  //  }
  //  VT_USER_END("name");
  //  MPI_Finalize();
  //  return EXIT_SUCCESS;

  //  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

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

    //set up the particle weight
    if (_NA_SPEC_<0) exit(__LINE__,__FILE__,"Error: the species that is usd for defining the weight is used used in the simulation");



    PIC::ParticleWeightTimeStep::GlobalParticleWeight=new double [PIC::nTotalSpecies];
    for (int i=0;i<PIC::nTotalSpecies;i++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[i]=1.0;


/*    PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
    PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_NA_SPEC_);

    //copy the weight and time step from Na neutra to Na ions
    if (_NAPLUS_SPEC_!=-1) {
      PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_NAPLUS_SPEC_,_NA_SPEC_,5.0E3/800.0E3);
      PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(_NAPLUS_SPEC_,_NA_SPEC_,5.0E3/800.0E3);
    }
*/

    //set photolytic reactions
  /*
    PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(sodiumPhotoionizationReactionProcessor,_NA_SPEC_11_);
    PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(sodiumPhotoionizationLifeTime,_NA_SPEC_11_);
  */


  //  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
  //  PIC::Mesh::mesh.outputMeshDataTECPLOT("mesh.data.dat",0);

//    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
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
    double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

    PIC::DistributionFunctionSample::vMin=-40.0E3;
    PIC::DistributionFunctionSample::vMax=40.0E3;
    PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

    PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);*/


  #if _MERCURY_FIPS_SAMPLING_ == _MERCURY_MODE_ON_
    //init sampling points along the s/c trajectory
    const int nFlybySamplePasses=6;
    const double FlybySamplingInterval=20.0*60.0,FlybySamplingIntervalStep=60.0; //in seconds
    const int nSampleSteps=(int)(FlybySamplingInterval/FlybySamplingIntervalStep);

    const char *FlybySamplePassesUTC[nFlybySamplePasses]={"2011-04-13T16:15:00","2011-04-14T16:30:00","2011-04-16T04:35:00",
        "2011-04-13T04:40:00","2011-04-15T04:55:00","2011-04-21T17:50:00"};

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



    for (nFlybyPass=0;nFlybyPass<nFlybySamplePasses;nFlybyPass++) {
      utc2et_c(FlybySamplePassesUTC[nFlybyPass],&et);

      if (PIC::ThisThread==0) {
        cout << "S/C Flyby Sampling: Pass=" << nFlybyPass << ":" << endl;
        cout << "Flux Sample Point\tUTS\t\t\t x[km]\t\ty[km]\t\tz[km]\t\t\t lx\t\tly\t\tlz\t" << endl;
      }

      for (n=0;n<nSampleSteps;n++) {
        //position of the s/c
        spkezr_c("MESSENGER",et,"LSO","NONE","MERCURY",state,&lt);


        //get the pointing vector in the 'SO' frame
        memcpy(bsight,bsight_INIT,3*sizeof(double));

        pxform_c ("MSGR_EPPS_FIPS","LSO",et,rotate);
        mxv_c(rotate,bsight,pointing);

        //print the pointing information
        if (PIC::ThisThread==0) {
          et2utc_c(et,"C",0,lenout,utcstr);
          fprintf(PIC::DiagnospticMessageStream,"%i\t\t\t%s\t",nFluxSamplePoint,utcstr);
          for (idim=0;idim<3;idim++) fprintf(PIC::DiagnospticMessageStream,"%e\t",state[idim]);

          cout << "\t";

          for (idim=0;idim<3;idim++) fprintf(PIC::DiagnospticMessageStream,"%e\t",pointing[idim]);
          cout << endl;
        }

        //save the samlpe pointing information
        for (idim=0;idim<3;idim++) {
          FluxSampleLocations[nFluxSamplePoint][idim]=state[idim]*1.0E3;
          FluxSampleDirections[nFluxSamplePoint][idim]=pointing[idim];
        }

        //increment the flyby time
        et+=FlybySamplingIntervalStep;
        ++nFluxSamplePoint;
      }

      if (PIC::ThisThread==0) cout << endl;
    }

    PIC::ParticleFluxDistributionSample::Init(FluxSampleLocations,FluxSampleDirections,30.0/180.0*Pi,nTotalFluxSamplePoints);

  #elif _MERCURY_FIPS_SAMPLING_ == _MERCURY_MODE_OFF_
    fprintf(PIC::DiagnospticMessageStream,"No FIPS sampling\n");
  #else
    exit(__LINE__,__FILE__,"Error: the option is not recognized");
  #endif

    //init ICES

  #ifdef _ICES_CREATE_COORDINATE_LIST_
    PIC::CPLR::ICES::createCellCenterCoordinateList();
    PIC::CPLR::ICES::SetLocationICES("/left/ices/ICES");
    PIC::CPLR::ICES::retriveSWMFdata("MERCURY_RESTART_n070100");  ////("MERCURY_RESTART_n070001");
  #endif


  #ifdef _ICES_LOAD_DATA_
    PIC::CPLR::ICES::readSWMFdata(1.0);
    PIC::Mesh::mesh.outputMeshDataTECPLOT("ices.data.dat",0);


    //create the map of the solar wind flux
    int el;

    for (el=0;el<PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;el++) {
      int i,j,k;
      long int nd;
      double FaceCenterPoint[3],PlasmaVelocity[3],PlasmaNumberDensity,FaceElementNormal[3],c;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
      PIC::Mesh::cDataCenterNode *CenterNode;
      char *offset;

      Moon::Planet->GetSurfaceElementMiddlePoint(FaceCenterPoint,el);
      Moon::Planet->GetSurfaceElementNormal(FaceElementNormal,el);

      if (FaceElementNormal[0]<0.0) continue;

      node=PIC::Mesh::mesh.findTreeNode(FaceCenterPoint);

      if ((nd=PIC::Mesh::mesh.fingCellIndex(FaceCenterPoint,i,j,k,node,false))==-1) {
        exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

  #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
      if ((node->Thread==PIC::ThisThread)&&(node->block==NULL)) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  #endif

      if (node->Thread==PIC::ThisThread) {
        CenterNode=node->block->GetCenterNode(nd);
        offset=CenterNode->GetAssociatedDataBufferPointer();

        if (*((int*)(offset+PIC::CPLR::ICES::DataStatusOffsetSWMF))==_PIC_ICES__STATUS_OK_) {
          memcpy(PlasmaVelocity,offset+PIC::CPLR::ICES::PlasmaBulkVelocityOffset,3*sizeof(double));
          memcpy(&PlasmaNumberDensity,offset+PIC::CPLR::ICES::PlasmaNumberDensityOffset,sizeof(double));
        }
        else {
          double EmptyArray[3]={0.0,0.0,0.0};

          memcpy(PlasmaVelocity,EmptyArray,3*sizeof(double));
          memcpy(&PlasmaNumberDensity,EmptyArray,sizeof(double));
        }

        c=-(PlasmaVelocity[0]*FaceElementNormal[0]+PlasmaVelocity[1]*FaceElementNormal[1]+PlasmaVelocity[2]*FaceElementNormal[2]);
        if (c<0.0) c=0.0;

        Moon::Planet->SolarWindSurfaceFlux[el]=c*PlasmaNumberDensity;
      }

      MPI_Bcast(Moon::Planet->SolarWindSurfaceFlux+el,1,MPI_DOUBLE,node->Thread,MPI_GLOBAL_COMMUNICATOR);
    }


  #endif

    //prepopulate the solar wind protons
  //  prePopulateSWprotons(PIC::Mesh::mesh.rootTree);


    /*
    //check the volume of local cells
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
    while (node!=NULL) {
      int i,j,k,di,dj,dk,idim;
      long int LocalCellNumber;
      double r,rmin=0.0,rmax=0.0,rprobe[3]={0.0,0.0,0.0};
      PIC::Mesh::cDataCenterNode *cell;

      if (node->Temp_ID==1862) {
        cout << __FILE__ << "@" << __LINE__ << endl;
      }

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
              cell=node->block->GetCenterNode(LocalCellNumber);
              rmin=-1,rmax=-1;

              if (cell->Measure<=0.0) {
                for (dk=0;dk<=((DIM==3) ? 1 : 0);dk++) for (dj=0;dj<=((DIM>1) ? 1 : 0);dj++) for (di=0;di<=1;di++) {
                  node->GetCornerNodePosition(rprobe,i+di,j+dj,k+dk);

                  for (idim=0,r=0.0;idim<DIM;idim++) r+=pow(rprobe[idim],2);
                  r=sqrt(r);

                  if ((rmin<0.0)||(rmin>r)) rmin=r;
                  if ((rmax<0.0)||(rmax<r)) rmax=r;
                }

                if ((rmin<rSphere)&&(rmax>rSphere)) {
                  cout << "Node ("<< i+di << "," << j+dj << "," << k+dk << "): r=" << rprobe[0] << "," << rprobe[1] << "," << rprobe[2] << ", |r|=" << r << endl;
                }

                PIC::Mesh::mesh.InitCellMeasure(node);

              }


            }
         }
      }

      node=node->nextNodeThisThread;
    }
  */

  //  VT_TRACER("main");









    //time step
    double SimulationTimeStep=-1.0;

    for (int spec=0;spec<PIC::nTotalSpecies;spec++) if ((SimulationTimeStep<0.0)||(SimulationTimeStep>PIC::ParticleWeightTimeStep::GlobalTimeStep[spec])) {
      SimulationTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
    }

  }

void amps_time_step() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
    //determine the parameters of the orbital motion of Mercury
    SpiceDouble StateBegin[6],StateEnd[6],lt,StateSun[6],StateMiddle[6];
    double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;
    int idim;

    SpiceDouble HCI_to_MSO_TransformationMartix[6][6];

    spkezr_c("Moon",Moon::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateBegin,&lt);
    spkezr_c("SUN",Moon::OrbitalMotion::et,"LSO","none","Moon",StateSun,&lt);

    double SimulationTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];

    //calculate position of the Earth ar the middle of the interation
    SpiceDouble StateEarth_SO[6],StateEarth_HCI[6];
    spkezr_c("Earth",Moon::OrbitalMotion::et+0.5*SimulationTimeStep,Moon::SO_FRAME,"none","MOON",StateEarth_SO,&lt);
    spkezr_c("Earth",Moon::OrbitalMotion::et+0.5*SimulationTimeStep,"MSGR_HCI","none","MOON",StateEarth_HCI,&lt);

    //calculate lunar velocity in an itertial frame, which have dirtectional vectors that coinsides with that of LSO
    sxform_c("MSGR_HCI","LSO",Moon::OrbitalMotion::et+0.5*SimulationTimeStep,HCI_to_MSO_TransformationMartix);
    spkezr_c("Moon",Moon::OrbitalMotion::et+0.5*SimulationTimeStep,"MSGR_HCI","none","SUN",StateMiddle,&lt);


    Moon::OrbitalMotion::et+=SimulationTimeStep;
    spkezr_c("Moon",Moon::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateEnd,&lt);


    for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
      StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
      StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;
      StateMiddle[idim]*=1.0E3,StateMiddle[3+idim]*=1.0E3;

      rBegin+=pow(StateBegin[idim],2);
      rEnd+=pow(StateEnd[idim],2);

      Moon::xObject_HCI[idim]=StateBegin[idim];
      Moon::vObject_HCI[idim]=StateBegin[3+idim];

      Moon::xSun_SO[idim]=1.0E3*StateSun[idim];
      Moon::vSun_SO[idim]=1.0E3*StateSun[3+idim];

      Moon::xEarth_SO[idim]=1.0E3*StateEarth_SO[idim];
      Moon::vEarth_SO[idim]=1.0E3*StateEarth_SO[3+idim];

      Moon::xEarth_HCI[idim]=1.0E3*StateEarth_HCI[idim];
      Moon::vEarth_HCI[idim]=1.0E3*StateEarth_HCI[3+idim];
    }

    //calculate parameters of SO_FROZEN
    //velocity of the coordinate frame
    for (idim=0;idim<3;idim++) {
      Moon::vObject_SO_FROZEN[idim]=
          HCI_to_MSO_TransformationMartix[idim][0]*StateMiddle[3+0]+
          HCI_to_MSO_TransformationMartix[idim][1]*StateMiddle[3+1]+
          HCI_to_MSO_TransformationMartix[idim][2]*StateMiddle[3+2];
    }

    //the axis of rotation of the MSO fraim in MSO_FROZEN during the next time step
    //get pointing direction to the Sun at the end of the current iteration in MSO_FROZEN
    SpiceDouble fmatrix[6][6];
    double SunPointingDirectionEnd[3],SunPointingDirectionEnd_MSO_FROZEN[3];

    //calculate Sun pointing at the end of the iteration in HCI frame (et is already incremented!!!!!!)
    sxform_c("LSO","MSGR_HCI",Moon::OrbitalMotion::et,fmatrix);

    SunPointingDirectionEnd[0]=fmatrix[0][0];
    SunPointingDirectionEnd[1]=fmatrix[1][0];
    SunPointingDirectionEnd[2]=fmatrix[2][0];

    //convert the pointing direction vector into MSO_FROZEN frame
    sxform_c("MSGR_HCI","LSO",Moon::OrbitalMotion::et-SimulationTimeStep,fmatrix);

    for (idim=0;idim<3;idim++) {
      SunPointingDirectionEnd_MSO_FROZEN[idim]=
          fmatrix[idim][0]*SunPointingDirectionEnd[0]+
          fmatrix[idim][1]*SunPointingDirectionEnd[1]+
          fmatrix[idim][2]*SunPointingDirectionEnd[2];
    }

    //calculate the rate of rotation in MSO_FROZEN
    Moon::RotationRate_SO_FROZEN=acos(SunPointingDirectionEnd_MSO_FROZEN[0])/SimulationTimeStep;


    //calculate the direction of rotation
    double c=sqrt(pow(SunPointingDirectionEnd_MSO_FROZEN[1],2)+pow(SunPointingDirectionEnd_MSO_FROZEN[2],2));

    if (c>0.0) {
      Moon::RotationVector_SO_FROZEN[0]=0.0;
      Moon::RotationVector_SO_FROZEN[1]=-SunPointingDirectionEnd_MSO_FROZEN[2]/c*Moon::RotationRate_SO_FROZEN;
      Moon::RotationVector_SO_FROZEN[2]=SunPointingDirectionEnd_MSO_FROZEN[1]/c*Moon::RotationRate_SO_FROZEN;
    }
    else {
      Moon::RotationVector_SO_FROZEN[0]=0.0;
      Moon::RotationVector_SO_FROZEN[1]=0.0;
      Moon::RotationVector_SO_FROZEN[2]=0.0;
    }


    //RECALCUALTE THE ROTATION VECTOR USING THE TRANSOFRMATON MARTICX FROM MSO_FROSEN at the time step (n) to the MSO_FROZEN at the time step (n+1)
    //the rotation vector is the eigrnvector of the transformation matrix
    //Zhuravlev, Osnovy teoreticheskoi mehaniki, Chapter 2, paragraph 6.2 (sposoby zadaniya orientacii tverdogo tela)

    //get the transformation matrix T(LSO[n]->LSO[n+1])=T1(LSO[n]->MSGR_HCI)*T2(MSGR_HCI->LSO[n+1])
    SpiceDouble T1[6][6],T2[6][6];
    double T[3][3];
    int i,j,k;


    sxform_c("LSO","MSGR_HCI",Moon::OrbitalMotion::et-SimulationTimeStep,T1);
    sxform_c("MSGR_HCI","LSO",Moon::OrbitalMotion::et,T2);


    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T[i][j]=0.0;

      for (k=0;k<3;k++) T[i][j]+=T1[i][k]*T2[k][j];
    }

    //determine the rate and the vectrot of the rotation
    double RotationAngle,t,RotationVector[3],RotationRate;

    RotationAngle=acos((T[0][0]+T[1][1]+T[2][2]-1.0)/2.0);

    t=2.0*sin(RotationAngle);
    RotationVector[0]=(T[2][1]-T[1][2])/t;
    RotationVector[1]=(T[0][2]-T[2][0])/t;
    RotationVector[2]=(T[1][0]-T[0][1])/t;

    RotationRate=RotationAngle/SimulationTimeStep;

    t=RotationRate/sqrt(RotationVector[0]*RotationVector[0]+RotationVector[1]*RotationVector[1]+RotationVector[2]*RotationVector[2]);
    RotationVector[0]*=t,RotationVector[1]*=t,RotationVector[2]*=t;


    rBegin=sqrt(rBegin);
    rEnd=sqrt(rEnd);

    for (idim=0;idim<3;idim++) {
      lBegin[idim]=StateBegin[idim]/rBegin;
      lEnd[idim]=StateEnd[idim]/rEnd;

      c0+=StateBegin[3+idim]*lBegin[idim];
      c1+=StateEnd[3+idim]*lEnd[idim];
    }

    Moon::xObjectRadial=0.5*(rBegin+rEnd);
    Moon::vObjectRadial=0.5*(c0+c1);

    //calculate TAA
    Moon::OrbitalMotion::TAA=Moon::OrbitalMotion::GetTAA(Moon::OrbitalMotion::et);

    for (idim=0;idim<3;idim++) {
      vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
      vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
    }

    vTangentialBegin=sqrt(vTangentialBegin);
    vTangentialEnd=sqrt(vTangentialEnd);

    Moon::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);


    //determine direction to the Sun and rotation angle in the coordiname frame related to the Moon
    SpiceDouble state[6],l=0.0;

    spkezr_c("SUN",Moon::OrbitalMotion::et,"IAU_MOON","none","MOON",state,&lt);

    for (idim=0;idim<3;idim++) l+=pow(state[idim],2);

    for (l=sqrt(l),idim=0;idim<3;idim++) {
      Moon::OrbitalMotion::SunDirection_IAU_OBJECT[idim]=state[idim]/l;
    }

    //matrixes for tranformation LSO->IAU and IAU->LSO coordinate frames
    sxform_c("LSO","IAU_MOON",Moon::OrbitalMotion::et,Moon::OrbitalMotion::SO_to_IAU_TransformationMartix);
    sxform_c("IAU_MOON","LSO",Moon::OrbitalMotion::et,Moon::OrbitalMotion::IAU_to_SO_TransformationMartix);
#endif



    //make the time advance
     PIC::TimeStep();
}


