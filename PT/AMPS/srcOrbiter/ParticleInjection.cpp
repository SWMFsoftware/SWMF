//$Id$
//particle injection procedures

/*
 *  ParticleInjection.cpp
 *
 *  Created on: Feb 11, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "SingleVariableDiscreteDistribution.h"

int Orbiter::InjectionModel::PointSource::InjectionDataTableLength=0;
Orbiter::InjectionModel::PointSource::cInjectionData Orbiter::InjectionModel::PointSource::InjectionDataTable[]={{0.0,0.0,{0.0,0.0,0.0},0}};


int Orbiter::InjectionModel::FaceEjection::InjectionDataTableLength=0;
Orbiter::InjectionModel::FaceEjection::cInjectionData Orbiter::InjectionModel::FaceEjection::InjectionDataTable[]={{0,0.0,0.0,0,NULL}};
bool Orbiter::InjectionModel::FaceEjection::SourceInitFlag=false;

int  Orbiter::InjectionModel::FaceEjection::InjectionFaceTableLength=0;
Orbiter::InjectionModel::FaceEjection::cFaceTableEntry *Orbiter::InjectionModel::FaceEjection::InjectionFaceTable=NULL;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inject particle from faces that has the same face attribute
void Orbiter::InjectionModel::FaceEjection::Init() {
  int i,j,iAttribute;

  if (SourceInitFlag==true) {
    //the model is already initialized -> exit
    return;
  }
  else SourceInitFlag=true;

  //determine the length of the face table that is needed (how many face attributed are set for the particle injection
  bool flag;
  vector<int> FaceAttributeTable;
  int nTotalInjectionFaceAttributes=0;

  for (i=0;i<InjectionDataTableLength;i++) {
    flag=true;

    for (j=0;j<i;j++) if (InjectionDataTable[i].faceat==InjectionDataTable[j].faceat) {
      //the face attribute is already registered
      flag=false;
      break;
    }

    if (flag==true) {
      //this is a new face attribute
      FaceAttributeTable.push_back(InjectionDataTable[i].faceat);
      nTotalInjectionFaceAttributes++;
    }
  }

  //allocate and populate 'InjectionFaceTable'
  InjectionFaceTable=new cFaceTableEntry[nTotalInjectionFaceAttributes];

  for (iAttribute=0;iAttribute<nTotalInjectionFaceAttributes;iAttribute++) {
    InjectionFaceTable[iAttribute].FaceTableLength=0;
    InjectionFaceTable[iAttribute].faceat=FaceAttributeTable[iAttribute];
  }

  //determine how many faces exists for each attribute
  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) for (iAttribute=0;iAttribute<nTotalInjectionFaceAttributes;iAttribute++) {
    if (InjectionFaceTable[iAttribute].faceat==CutCell::BoundaryTriangleFaces[i].attribute) {
      InjectionFaceTable[iAttribute].FaceTableLength++;
      break;
    }
  }

  //populate the FaceTable
  for (iAttribute=0;iAttribute<nTotalInjectionFaceAttributes;iAttribute++) {
    InjectionFaceTable[iAttribute].FaceTable=new int [InjectionFaceTable[iAttribute].FaceTableLength];
    InjectionFaceTable[iAttribute].FaceTableLength=0;
    InjectionFaceTable[iAttribute].TotalArea=0.0;
  }

  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) for (iAttribute=0;iAttribute<nTotalInjectionFaceAttributes;iAttribute++) {
    if (InjectionFaceTable[iAttribute].faceat==CutCell::BoundaryTriangleFaces[i].attribute) {

      InjectionFaceTable[iAttribute].FaceTable[InjectionFaceTable[iAttribute].FaceTableLength]=i;
      InjectionFaceTable[iAttribute].TotalArea+=CutCell::BoundaryTriangleFaces[i].SurfaceArea;

      InjectionFaceTable[iAttribute].FaceTableLength++;
      break;
    }
  }

  //initialize the object used for distributing the face where injection of the next particle will occur
  for (iAttribute=0;iAttribute<nTotalInjectionFaceAttributes;iAttribute++) {
    double FaceAreaTable[InjectionFaceTable[iAttribute].FaceTableLength];

    for (i=0;i<InjectionFaceTable[iAttribute].FaceTableLength;i++) FaceAreaTable[i]=CutCell::BoundaryTriangleFaces[InjectionFaceTable[iAttribute].FaceTable[i]].SurfaceArea;

    //init the surface distribution module
    InjectionFaceTable[iAttribute].InjectionFaceGenerator.InitArray(FaceAreaTable,InjectionFaceTable[iAttribute].FaceTableLength,10*InjectionFaceTable[iAttribute].FaceTableLength);

    //add pointer to InjectionFaceTable[].FaceTable to all appropriate sources
    for (i=0;i<InjectionDataTableLength;i++) if (InjectionDataTable[i].faceat==InjectionFaceTable[iAttribute].faceat) InjectionDataTable[i].FaceTable=InjectionFaceTable+iAttribute;
  }
}

//inject model particles from faces
long int Orbiter::InjectionModel::FaceEjection::InjectParticles() {
  int spec,iSource;
  double ParticleWeight,TimeStep,ModelParticlesInjectionRate,TimeCounter;
  double x[3],v[3];
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int nTotalInjectedParticles=0;

  //check whether the model is initialized
  if (SourceInitFlag==false) Init();

  //the model is developed for the case of the global time step and global particle weight
  if (_SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_) exit(__LINE__,__FILE__,"Error: the model should be used with the species dependent global time step");
  if (_SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_) exit(__LINE__,__FILE__,"Error: the model should be used with the specied dependent global particle weight");


  for (iSource=0;iSource<InjectionDataTableLength;iSource++) {
    spec=InjectionDataTable[iSource].Species;

    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    TimeStep=PIC::ParticleWeightTimeStep::GetGlobalTimeStep(spec);

    //determine the model particle injection rate
    ModelParticlesInjectionRate=InjectionDataTable[iSource].SourceRate*InjectionDataTable[iSource].FaceTable->TotalArea/ParticleWeight;
    TimeCounter-=-log(rnd())/ModelParticlesInjectionRate*rnd();

    while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<TimeStep) {
      //distribute the particle velocity
      double BulkFlowVelocity[3]={0.0,0.0,0.0};
      int i,iFace;

      //generate location of the new particle
      i=InjectionDataTable[iSource].FaceTable->InjectionFaceGenerator.DistributeVariable();
      iFace=InjectionDataTable[iSource].FaceTable->FaceTable[i];
      CutCell::BoundaryTriangleFaces[iFace].GetRandomPosition(x);

      startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
      if (startNode->Thread!=PIC::ThisThread) continue;

      //generate velocity of the new particle and inject it into the system
      do {
        PIC::Distribution::MaxwellianVelocityDistribution(v,BulkFlowVelocity,InjectionDataTable[iSource].Temperature,spec);
      }
      while (Vector3D::DotProduct(v,CutCell::BoundaryTriangleFaces[iFace].ExternalNormal)<=0.0);

      //generate the new particle, and inject it into the system
      newParticle=PIC::ParticleBuffer::GetNewParticle();
      newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
      nTotalInjectedParticles++;

      //sample the source rate
      CutCell::BoundaryTriangleFaces[iFace].UserData.SourceRate[spec]+=ParticleWeight/TimeStep/CutCell::BoundaryTriangleFaces[iFace].SurfaceArea;

      PIC::ParticleBuffer::SetX(x,newParticleData);
      PIC::ParticleBuffer::SetV(v,newParticleData);
      PIC::ParticleBuffer::SetI(spec,newParticleData);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

      //apply condition of tracking the particle
      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      PIC::ParticleTracker::InitParticleID(newParticleData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
      #endif

       //inject the particle into the system
      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_FACE_(newParticle,TimeStep*rnd(),startNode,true,CutCell::BoundaryTriangleFaces+iFace);
    }
  }

  return nTotalInjectedParticles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//inject particle from a point source
long int Orbiter::InjectionModel::PointSource::InjectParticles() {
  double TimeCounter,ModelParticlesInjectionRate;
  int spec,idim;
  double *x,v[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nTotalInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
  int iSourceLocation;

/*
  for (iSourceLocation=0;iSourceLocation<InjectionDataTableLength;iSourceLocation) {
    x=InjectionDataTable[iSourceLocation].Location;
    spec=InjectionDataTable[iSourceLocation].Species;

    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::ThisThread) continue;

    double ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    double TimeStep=startNode->block->GetLocalTimeStep(spec);

    //determine the model particle injection rate
    ModelParticlesInjectionRate=InjectionDataTable[iSourceLocation].SourceRate/ParticleWeight;
    TimeCounter-=-log(rnd())/ModelParticlesInjectionRate*rnd();

    while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<TimeStep) {
      //distribute the particle velocity
      double BulkFlowVelocity[3]={0.0,0.0,0.0};

      PIC::Distribution::MaxwellianVelocityDistribution(v,BulkFlowVelocity,InjectionDataTable[iSourceLocation].Temperature,spec);

      //generate the new particle, and inject it into the system
      newParticle=PIC::ParticleBuffer::GetNewParticle();
      newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
      nTotalInjectedParticles++;

      PIC::ParticleBuffer::SetX(x,newParticleData);
      PIC::ParticleBuffer::SetV(v,newParticleData);
      PIC::ParticleBuffer::SetI(spec,newParticleData);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

      //apply condition of tracking the particle
      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      PIC::ParticleTracker::InitParticleID(newParticleData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
      #endif

       //inject the particle into the system
      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,TimeStep*rnd(),startNode,true);
    }
  }*/

  exit(__LINE__,__FILE__,"Error: implementations is not compelete");

  return nTotalInjectedParticles;
}

//inject particle from a "ring"
long int Orbiter::InjectionModel::Ring::InjectParticles() {
  double TimeCounter,ModelParticlesInjectionRate;
  int spec,idim;
  double x[3],v[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nTotalInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;

  //the produre works only for the case of the "global" time step and "global" particle weight
  if (_SIMULATION_PARTICLE_WEIGHT_MODE_ != _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_) exit(__LINE__,__FILE__,"Error: the procedure is developed only for the case of global particle weight");
  if (_SIMULATION_TIME_STEP_MODE_ != _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_) exit(__LINE__,__FILE__,"Error: the procedure is developed only for the case of a global time step");

/*
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    double ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
    double TimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];

    //determine the model particle injection rate
    ModelParticlesInjectionRate=Orbiter::InjectionModel::Ring::SourceRate[spec]/ParticleWeight;
    TimeCounter-=-log(rnd())/ModelParticlesInjectionRate*rnd();

    while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<TimeStep) {
      //distribute the location of a point
      double phi=PiTimes2*rnd();
      double RsinPhi,RcosPhi;

      RsinPhi=Radius*sin(phi);
      RcosPhi=Radius*cos(phi);

      for (idim=0;idim<3;idim++) x[idim]=x0[idim]+RsinPhi*e0[idim]+RcosPhi*e1[idim];

      startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
      if (startNode->Thread!=PIC::ThisThread) continue;

      //distribute the particle velocity
      double BulkFlowVelocity[3]={0.0,0.0,0.0};

      PIC::Distribution::MaxwellianVelocityDistribution(v,BulkFlowVelocity,SourceTemperature[spec],spec);

      //generate the new particle, and inject it into the system
      newParticle=PIC::ParticleBuffer::GetNewParticle();
      newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
      nTotalInjectedParticles++;

      PIC::ParticleBuffer::SetX(x,newParticleData);
      PIC::ParticleBuffer::SetV(v,newParticleData);
      PIC::ParticleBuffer::SetI(spec,newParticleData);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

      //apply condition of tracking the particle
      #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
      PIC::ParticleTracker::InitParticleID(newParticleData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
      #endif

       //inject the particle into the system
      _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,TimeStep*rnd(),startNode,true);
    }
  }
*/

  exit(__LINE__,__FILE__,"Error: implementation is not compelte");

  return nTotalInjectedParticles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//the manager of the particle injection module
long int  Orbiter::InjectionModel::InjectParticles() {
  long int res=0;

  if (FaceEjection::InjectionDataTableLength!=0) res+=FaceEjection::InjectParticles();
  if (PointSource::InjectionDataTableLength!=0) res+=PointSource::InjectParticles();

  return res;
}

double Orbiter::InjectionModel::GetTotalInjectionRate(int spec) {
  double res=0.0;

  if (FaceEjection::InjectionDataTableLength!=0) {
    for (int iSource=0;iSource<FaceEjection::InjectionDataTableLength;iSource++) if (spec==FaceEjection::InjectionDataTable[iSource].Species) {
      if (FaceEjection::SourceInitFlag==false) FaceEjection::Init();

      res+=FaceEjection::InjectionDataTable[iSource].SourceRate*FaceEjection::InjectionDataTable[iSource].FaceTable->TotalArea;
    }
  }

  return res;
}




