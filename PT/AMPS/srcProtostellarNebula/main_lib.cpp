

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


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "ProtostellarNebula.h"




//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;
const double rSphere=_RADIUS_(_TARGET_);


const double xMaxDomain=30; //modeling the vicinity of the planet
const double yMaxDomain=30; //the minimum size of the domain in the direction perpendicular to the direction to the sun

const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100/2.5,dxMaxSphere=DebugRunMultiplier*2.0/10.0;





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

  if (r>dxMinGlobal*rSphere) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;

  return rSphere*res;
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;

  double CharacteristicSpeed_NA=1.0E6;

//  CharacteristicSpeed*=sqrt(PIC::MolecularData::GetMass(NA)/PIC::MolecularData::GetMass(spec));

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed_NA;


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

    MPI_Gather(&nTotalGeneratedParticles,1,MPI_LONG,GeneratedParticle,1,MPI_LONG,0,MPI_COMM_WORLD);
    MPI_Gather(&nTotalProcessorBlocks,1,MPI_LONG,GeneratedNodes,1,MPI_LONG,0,MPI_COMM_WORLD);
    MPI_Gather(&aNpartTotal,1,MPI_DOUBLE,anpart,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(&TotalDomainVolume,1,MPI_DOUBLE,volume,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


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

double sphereInjectionRate(int spec,int BoundaryElementType,void *BoundaryElement) {

  double res=1.0E20;

  return res;
}



long int sphereParticleInjection(int spec,int BoundaryElementType,void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ParticleWeight,LocalTimeStep,/*ExternalNormal[3],*/x[3],v[3],/*r,*/*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
//  int idim;

  double ParticleWeightCorrection=1.0;

//  static const double Temp=200.0;
//  double vbulk[3]={0.0,0.0,0.0};


//  return 0;

//====================  DEBUG ===========================
//  static bool FirstPArticleGenerated=false;
//====================  END DEBUG ===================================


  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  /*
  double TimeCounter=0.0;
  double ModelParticlesInjectionRate=sphereInjectionRate(NA,SphereDataPointer)/ParticleWeight;

  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {

*/

  static double InjectParticles=0.0;

  bool flag;
  int idim;
  double r,ExternalNormal[3];

  InjectParticles+=sphereInjectionRate(spec,BoundaryElementType,SphereDataPointer)*LocalTimeStep;

  while (InjectParticles>0.0) {

    //generate the particle position
    for (r=0.0,idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      r+=ExternalNormal[idim]*ExternalNormal[idim];
    }

    r=-sqrt(r);

    for (idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]/=r;
      x[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
    }

    InjectParticles-=ParticleWeight;

    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) continue;

    //generate the particle velocity
//    PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,Temp,ExternalNormal,NA);


    for (idim=0;idim<3;idim++) v[idim]=-ExternalNormal[idim]*3.0e3;


//    InjectParticles-=ParticleWeight*ParticleWeightCorrection;
//    if (flag==false) continue;

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




    if (startNode->block->GetLocalTimeStep(spec)/LocalTimeStep<rnd()) continue;

    //generate a particle
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    nInjectedParticles++;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(v,newParticleData);
    PIC::ParticleBuffer::SetI(spec,newParticleData);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);


    //inject the particle into the system
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd() /*LocalTimeStep-TimeCounter*/,startNode,true);
//    PIC::Mover::MoveParticleBoundaryInjection[NA](newParticle,startNode->block->GetLocalTimeStep(NA)*rnd() /*LocalTimeStep-TimeCounter*/,startNode,true);
  }

  return nInjectedParticles;
}

long int sphereParticleInjection(int BoundaryElementType,void *BoundaryElement) {
  long int spec,res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=sphereParticleInjection(spec,BoundaryElementType,BoundaryElement);

  return res;
}

long int sphereParticleInjection(void *SphereDataPointer)  {
  long int res=0.0;
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=sphereParticleInjection(spec,SphereDataPointer);
  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,0.0,0.0},nNA=5.0E6,tempNA=1.0E5;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O_SPEC_);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }

  return false;
}

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;

  if (spec!=_O_SPEC_ && spec!=_H_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={0.0,0.0,0.0},nNA=5.0E6,tempNA=1.0E5;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,spec);


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

          //generate particles' velocity
          PIC::Distribution::InjectMaxwellianDistribution(v,vNA,tempNA,ExternalNormal,spec,-1);

          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

          //inject the particle into the system
          _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
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

double BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  bool ExternalFaces[6];
  double ExternalNormal[3],BlockSurfaceArea;
  int nface;

  if (spec!=_O_SPEC_ && spec!=_H_SPEC_) return 0; //inject only spec=0

  double ModelParticlesInjectionRate=0.0;
  static double vNA[3]={0.0,0.0,0.0},nNA=5.0E6,tempNA=1.0E5;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
	startNode->GetExternalNormal(ExternalNormal,nface);
	BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);
	ModelParticlesInjectionRate+=BlockSurfaceArea*PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,spec);
      }
  }

  return ModelParticlesInjectionRate;
}


void amps_init() {
  PIC::InitMPI();


  //SetUp the alarm
//  PIC::Alarm::SetAlarm(2000);


  rnd_seed();
  MPI_Barrier(MPI_COMM_WORLD);


  //init the Mercury model
  ProtostellarNebula::Init_BeforeParser();
  PIC::Init_BeforeParser();

  ProtostellarNebula::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;
  ProtostellarNebula::Init_AfterParser();



  //register the sphere
  {
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


    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceMesh(fname);

    sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceData(fname,0);


    Sphere->localResolution=localSphericalSurfaceResolution;
    Sphere->InjectionRate=ProtostellarNebula::SourceProcesses::totalProductionRate; 
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=ProtostellarNebula::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
    Sphere->InjectionBoundaryCondition=ProtostellarNebula::SourceProcesses::InjectionBoundaryModel; 

    Sphere->PrintTitle=ProtostellarNebula::Sampling::OutputSurfaceDataFile::PrintTitle;
    Sphere->PrintVariableList=ProtostellarNebula::Sampling::OutputSurfaceDataFile::PrintVariableList;
    Sphere->PrintDataStateVector=ProtostellarNebula::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

    //set up the planet pointer in Mercury model
    ProtostellarNebula::Planet=Sphere;
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


  if (PIC::Mesh::mesh.ThisThread==0) {
    PIC::Mesh::mesh.buildMesh();
    PIC::Mesh::mesh.saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh.readMeshFile("mesh.msh");
  }

  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");

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


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  //set up the particle weight
  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=BoundingBoxInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O_SPEC_);
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H_SPEC_);




  MPI_Barrier(MPI_COMM_WORLD);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
//  double TimeCounter=time(NULL);
//  int LastDataOutputFileNumber=-1;


/*  //init the sampling of the particls' distribution functions: THE DECLARATION IS MOVED INTO THE INPUT FILE
  const int nSamplePoints=3;
  double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

  PIC::DistributionFunctionSample::vMin=-40.0E3;
  PIC::DistributionFunctionSample::vMax=40.0E3;
  PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

  PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);*/
}







  //time step
// for (long int niter=0;niter<100000001;niter++) {
void amps_time_step(){

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  
  int idim;

    //determine the parameters of the orbital motion of Mercury
    SpiceDouble StateBegin[6],StateEnd[6],lt,StateSun[6],StateMiddle[6];
    double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;

    SpiceDouble HCI_to_MSO_TransformationMartix[6][6];

    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateBegin,&lt);
    spkezr_c("SUN",Mercury::OrbitalMotion::et,"MSGR_MSO","none","Mercury",StateSun,&lt);

    //calculate Mercury's velocity in an itertial frame, which have dirtectional vectors that coinsides with that of MSO
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],HCI_to_MSO_TransformationMartix);
    spkezr_c("Mercury",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],"MSGR_HCI","none","SUN",StateMiddle,&lt);


    Mercury::OrbitalMotion::et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateEnd,&lt);


    for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
      StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
      StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;
      StateMiddle[idim]*=1.0E3,StateMiddle[3+idim]*=1.0E3;

      rBegin+=pow(StateBegin[idim],2);
      rEnd+=pow(StateEnd[idim],2);

      Mercury::xObject_HCI[idim]=StateBegin[idim];
      Mercury::vObject_HCI[idim]=StateBegin[3+idim];

      Mercury::xSun_SO[idim]=1.0E3*StateSun[idim];
      Mercury::vSun_SO[idim]=1.0E3*StateSun[3+idim];
    }

    //calculate parameters of SO_FROZEN
    //velocity of the coordinate frame
    for (idim=0;idim<3;idim++) {
      Mercury::vObject_SO_FROZEN[idim]=
          HCI_to_MSO_TransformationMartix[idim][0]*StateMiddle[3+0]+
          HCI_to_MSO_TransformationMartix[idim][1]*StateMiddle[3+1]+
          HCI_to_MSO_TransformationMartix[idim][2]*StateMiddle[3+2];
    }

    //the axis of rotation of the MSO fraim in MSO_FROZEN during the next time step
    //get pointing direction to the Sun at the end of the current iteration in MSO_FROZEN
    SpiceDouble fmatrix[6][6];
    double SunPointingDirectionEnd[3],SunPointingDirectionEnd_MSO_FROZEN[3];

    //calculate Sun pointing at the end of the iteration in HCI frame (et is already incremented!!!!!!)
    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,fmatrix);

    SunPointingDirectionEnd[0]=fmatrix[0][0];
    SunPointingDirectionEnd[1]=fmatrix[1][0];
    SunPointingDirectionEnd[2]=fmatrix[2][0];

    //convert the pointing direction vector into MSO_FROZEN frame
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],fmatrix);

    for (idim=0;idim<3;idim++) {
      SunPointingDirectionEnd_MSO_FROZEN[idim]=
          fmatrix[idim][0]*SunPointingDirectionEnd[0]+
          fmatrix[idim][1]*SunPointingDirectionEnd[1]+
          fmatrix[idim][2]*SunPointingDirectionEnd[2];
    }

    //calculate the rate of rotation in MSO_FROZEN
    Mercury::RotationRate_SO_FROZEN=acos(SunPointingDirectionEnd_MSO_FROZEN[0])/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];


    //calculate the direction of rotation
    double c=sqrt(pow(SunPointingDirectionEnd_MSO_FROZEN[1],2)+pow(SunPointingDirectionEnd_MSO_FROZEN[2],2));

    if (c>0.0) {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=-SunPointingDirectionEnd_MSO_FROZEN[2]/c*Mercury::RotationRate_SO_FROZEN;
      Mercury::RotationVector_SO_FROZEN[2]=SunPointingDirectionEnd_MSO_FROZEN[1]/c*Mercury::RotationRate_SO_FROZEN;
    }
    else {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=0.0;
      Mercury::RotationVector_SO_FROZEN[2]=0.0;
    }


    //RECALCUALTE THE ROTATION VECTOR USING THE TRANSOFRMATON MARTICX FROM MSO_FROSEN at the time step (n) to the MSO_FROZEN at the time step (n+1)
    //the rotation vector is the eigrnvector of the transformation matrix
    //Zhuravlev, Osnovy teoreticheskoi mehaniki, Chapter 2, paragraph 6.2 (sposoby zadaniya orientacii tverdogo tela)

    //get the transformation matrix T(LSO[n]->LSO[n+1])=T1(LSO[n]->MSGR_HCI)*T2(MSGR_HCI->LSO[n+1])
    SpiceDouble T1[6][6],T2[6][6];
    double T[3][3];
    int i,j,k;



    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,T1);
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);




    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T[i][j]=0.0;

      for (k=0;k<3;k++) T[i][j]+=T2[i][k]*T1[k][j];
    }

    //determine the rate and the vectrot of the rotation
    double RotationAngle,t,RotationVector[3],RotationRate;

    RotationAngle=acos((T[0][0]+T[1][1]+T[2][2]-1.0)/2.0);

    t=2.0*sin(RotationAngle);
    RotationVector[0]=(T[2][1]-T[1][2])/t;
    RotationVector[1]=(T[0][2]-T[2][0])/t;
    RotationVector[2]=(T[1][0]-T[0][1])/t;

    RotationRate=RotationAngle/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];

//    t=RotationRate/sqrt(RotationVector[0]*RotationVector[0]+RotationVector[1]*RotationVector[1]+RotationVector[2]*RotationVector[2]);
//    RotationVector[0]*=t,RotationVector[1]*=t,RotationVector[2]*=t;

t=1.0;

    //TEST THE ROTATION RATE AND THE ROTATION VECTOR
    double testRoptationMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double cosRotationAngle,sinRotationAngle;

    cosRotationAngle=cos(RotationAngle);
    sinRotationAngle=sin(RotationAngle);

    for (i=0;i<3;i++) testRoptationMatrix[i][i]+=cosRotationAngle;
    for (i=0;i<3;i++) for (j=0;j<3;j++) testRoptationMatrix[i][j]+=(1.0-cosRotationAngle)*RotationVector[i]*RotationVector[j]/pow(t,2);

    testRoptationMatrix[0][1]-=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[0][2]+=sinRotationAngle*RotationVector[1]/t;
    testRoptationMatrix[1][0]+=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[1][2]-=sinRotationAngle*RotationVector[0]/t;
    testRoptationMatrix[2][0]-=sinRotationAngle*RotationVector[1]/t,testRoptationMatrix[2][1]+=sinRotationAngle*RotationVector[0]/t;


    //CALCULATE THE EFECT OF THE TRANSFORMATION AND TRANSTER THE TRANSFORMED VEWCTROS TO HCI
    double T3[3][3];

    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T3[i][j]=0.0;

      for (k=0;k<3;k++) T3[i][j]+=T2[i][k]*testRoptationMatrix[k][j];
    }


    //GET THE ROTATION MATRIX FROM SPICE
//    SpiceDouble rVect1[3],rVect2[3],rVect[3],rot[3][3],xform[6][6];

//    sxform_c (  "MSGR_HCI","MSGR_MSO", et, tsipm ) ;

    //RECALCULATE THE MATRIX AGIN
    double newRotationVector[3],newRate;

    newRate=Exosphere::OrbitalMotion::FrameRotation::GetRotationVector(newRotationVector,"MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],Mercury::OrbitalMotion::et);


    //END OF TRANSFORMATIUON TESTS ------------------

    rBegin=sqrt(rBegin);
    rEnd=sqrt(rEnd);

    for (idim=0;idim<3;idim++) {
      lBegin[idim]=StateBegin[idim]/rBegin;
      lEnd[idim]=StateEnd[idim]/rEnd;

      c0+=StateBegin[3+idim]*lBegin[idim];
      c1+=StateEnd[3+idim]*lEnd[idim];
    }

    Mercury::xObjectRadial=0.5*(rBegin+rEnd);
    Mercury::vObjectRadial=0.5*(c0+c1);

    //calculate TAA
    Mercury::OrbitalMotion::TAA=Mercury::OrbitalMotion::GetTAA(Mercury::OrbitalMotion::et);

    for (idim=0;idim<3;idim++) {
      vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
      vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
    }

    vTangentialBegin=sqrt(vTangentialBegin);
    vTangentialEnd=sqrt(vTangentialEnd);

    Mercury::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);


    //determine direction to the Sun and rotation angle in the coordiname frame related to Mercury
    SpiceDouble state[6],l=0.0;

    spkezr_c("SUN",Mercury::OrbitalMotion::et,"IAU_MERCURY","none","MERCURY",state,&lt);

    for (idim=0;idim<3;idim++) l+=pow(state[idim],2);

    for (l=sqrt(l),idim=0;idim<3;idim++) {
      Mercury::OrbitalMotion::SunDirection_IAU_OBJECT[idim]=state[idim]/l;
    }
    
    //matrixes for tranformation MSO->IAU and IAU->MSO coordinate frames
    sxform_c("MSGR_MSO","IAU_MERCURY",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::SO_to_IAU_TransformationMartix);
    sxform_c("IAU_MERCURY","MSGR_MSO",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::IAU_to_SO_TransformationMartix);
 
#endif

    //make the time advance
    static int LastDataOutputFileNumber=0;

    //make the time advance
     PIC::TimeStep();

     // write output file
     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;
       
       
       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
     
}



