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


//  #define _ICES_CREATE_COORDINATE_LIST_

//the species
//unsigned int _DUST_SPEC_=0;

//the particle class
//#include "pic.h"
#include "constants.h"

#include "Enceladus.h"


//the parameters of the domain and the sphere
const double DebugRunMultiplier=3.0;

const double rSphere=_RADIUS_(_TARGET_);
const double xMaxDomain=5.0;
const double dxMinGlobal=1.0,dxMaxGlobal=5.0;
const double dxMinSphere=DebugRunMultiplier*2.0/100,dxMaxSphere=DebugRunMultiplier*4.0/100.0;


//Preprocessor of the SWMF data the obtained through ICES -> rotate the vectors
//Need to rotate: E[3],B[3],swVel[3]
void SWMFdataPreProcessor(double *x,PIC::CPLR::ICES::cDataNodeSWMF& data) {
  VectorRotation::Along_Z_direction(data.B,-270.0/180.0*Pi);
  VectorRotation::Along_Z_direction(data.E,-270.0/180.0*Pi);
  VectorRotation::Along_Z_direction(data.swVel,-270.0/180.0*Pi);
}




//double SodiumRadiationPressureAcceleration_Combi_1997_icarus(double HeliocentricVelocity,double HeliocentricDistance);
void TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {


  accl[0]=0.0,accl[1]=0.0,accl[2]=0.0;

//  return;

  //the gravity force
  double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
  }



}






//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double r,l[3];
  int idim;
  double SubsolarAngle;

  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  SubsolarAngle=acos(l[0]);


  SubsolarAngle=0.0;

//  res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;



//  return rSphere*res;
  double Multiplier=0.25;

  if (x[2]>0.0) Multiplier=1.0;

  return Multiplier*DebugRunMultiplier*0.5*(dxMaxSphere+dxMinSphere)*rSphere*4; //*4
}

double localResolution(double *x) {
  int idim;
  double res=0.0,r=0.0;


  if ((x[2]<0.0)&&(asin(sqrt((x[0]*x[0]+x[1]*x[1])/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])))/Pi*180.0<45.0)) return DebugRunMultiplier*0.25*(dxMaxSphere+dxMinSphere)*rSphere*4; //*4

  return DebugRunMultiplier*rSphere*dxMaxGlobal*4; //*4


  for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

  r=sqrt(r);

  if ((3.5E6<r)&&(r<3.8E6)) return DebugRunMultiplier*100.0E3;
  else return DebugRunMultiplier*rSphere*dxMaxGlobal;


  /*
  if (r>dxMinGlobal*rSphere) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;

*/

  return DebugRunMultiplier*rSphere*res; 
}

//set up the local time step
//determine the limits on the "geometrical (flight time) times step

void GetTimeStepLimits(double &minTimeStep,double &maxTimeStep,int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double dt;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {

    if (ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode)==false) exit(__LINE__,__FILE__,"Error: not only the dust species are present in the system");

    if (dt>maxTimeStep) maxTimeStep=dt;
    if ((minTimeStep<0.0)||(dt<minTimeStep)) minTimeStep=dt;
  }
  else for (int i=0;i<(1<<_MESH_DIMENSION_);i++) if (startNode->downNode[i]!=NULL) GetTimeStepLimits(minTimeStep,maxTimeStep,spec,startNode->downNode[i]);
}



double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double dt;


  /*
  const double TimeStepMinValue=minParticleMovingInterval/CharacteristicSpeed_NA;



  static double *minTimeStep=NULL,*maxTimeStep=NULL;

  if (minTimeStep==NULL) {
    //init the buffers

    if (PIC::ThisThread==0) cout << "Determine the time step limits:" << endl;

    minTimeStep=new double[PIC::nTotalSpecies];
    maxTimeStep=new double[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      minTimeStep[s]=-1.0,maxTimeStep[s]=-1.0;
      GetTimeStepLimits(minTimeStep[s],maxTimeStep[s],s,PIC::Mesh::mesh.rootTree);

      if (PIC::ThisThread==0){
        cout << "spec=" << s << ", minTimeStep=" << minTimeStep[s] << ", maxTimeStep=" << maxTimeStep[s] << ": The time step range that will be actually used is (" << TimeStepMinValue << "," << maxTimeStep[s] << ")" << endl;
      }
    }
  }

  //get the local "geometrical" time step
  if (ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode)==false) exit(__LINE__,__FILE__,"Error: not only the dust species are present in the system");


  //modify the "geometrical" time step
  dt=TimeStepMinValue+(maxTimeStep[spec]-TimeStepMinValue)/(maxTimeStep[spec]-minTimeStep[spec])*(dt-minTimeStep[spec]);
*/

  //get the local "geometrical" time step
  if (ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode)==false) exit(__LINE__,__FILE__,"Error: not only the dust species are present in the system");

  return dt;
}


bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {


  return false;
}








double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
  double r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  double c=v[0]*x[0]/r+v[1]*x[1]/r+v[2]*x[2]/r;

  v[0]-=2.0*c*x[0]/r;
  v[1]-=2.0*c*x[1]/r;
  v[2]-=2.0*c*x[2]/r;

  return _PARTICLE_REJECTED_ON_THE_FACE_;

  /*
   //delete all particles that was not reflected on the surface
   PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;*/
}




int main(int argc,char **argv) {
  PIC::InitMPI();

  //  MPI_Init(&argc,&argv);

  rnd_seed();






  //==================  set up the PIC solver

  //  char inputFile[]="enceladus.input";


  //  MPI_Barrier(MPI_COMM_WORLD);

  //set up the alarm
  PIC::Alarm::SetAlarm(8.0*3600.0-15*60);

  //init the particle solver
  PIC::Init_BeforeParser();
  Enceladus::Init_BeforeParser();

  //  PIC::Parser::Run(inputFile);



  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=100;
  PIC::RequiredSampleLength=10; //200; //0;


  //register the sphere
  double sx0[3]={0.0,0.0,0.0};
  cInternalBoundaryConditionsDescriptor SphereDescriptor;
  cInternalSphericalData *Sphere;

  SetGeneralSphereSurfaceMeshParameters(200,200);

  PIC::BC::InternalBoundary::Sphere::Init();
  SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
  Sphere->SetSphereGeometricalParameters(sx0,rSphere);

  Sphere->localResolution=localSphericalSurfaceResolution;
  Sphere->InjectionRate=Enceladus::sphereInjectionRate;
  Sphere->faceat=0;
  Sphere->ParticleSphereInteraction=Enceladus::ParticleSphereInteraction;
  Sphere->InjectionBoundaryCondition=Enceladus::DustInjection__Sphere;

  Sphere->PrintSurfaceMesh("Sphere.dat");
  Sphere->PrintSurfaceData("SpheraData.dat",0);

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

    DomainCenterOffset[0]=-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
    DomainCenterOffset[1]=-DomainLength[1]/2.0;
    DomainCenterOffset[2]=-DomainLength[2]/2.0;
  }

  for (idim=0;idim<DIM;idim++) {
    xmin[idim]=DomainCenterOffset[idim];
    xmax[idim]=DomainLength[idim]+DomainCenterOffset[idim];
  }



  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,localResolution);
  PIC::Mesh::mesh.memoryAllocationReport();



  if (PIC::Mesh::mesh.ThisThread==0) {
    PIC::Mesh::mesh.buildMesh();
    PIC::Mesh::mesh.saveMeshFile("mesh.msh");
    //    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    // MPI_Barrier(MPI_COMM_WORLD);
    PIC::Mesh::mesh.readMeshFile("mesh.msh");
  }


  cout << __LINE__ << " rnd=" << rnd() << " " << PIC::Mesh::mesh.ThisThread << endl;

  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();

//  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();

  //initialize the blocks
  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();

//  PIC::Mesh::mesh.checkMeshConsistency(PIC::Mesh::mesh.rootTree);

  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();


  //============================== DEBUG =========================
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

    while (node!=NULL) {
      if (node->Temp_ID==2478) {
        cout << __FILE__<< "@" << __LINE__ << endl;

        PIC::Mesh::mesh.InitCellMeasure(node);
      }

      node=node->nextNodeThisThread;
    }

  //============================== END DEBUG ====================================



    /*
  //set up the output of the mars model production rate
  PIC::Mesh::PrintVariableListCenterNode.push_back(newMars::PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(newMars::PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(newMars::Interpolate);
*/


//init the PIC solver
  PIC::Init_AfterParser();
  PIC::Mover::Init();
  Enceladus::Init_AfterParser();
//  PIC::Mover::TotalParticleAcceleration=ElectricallyChargedDust::TotalGrainAcceleration; ///TotalParticleAcceleration;


/*
  for (int s=0;s<PIC::nTotalSpecies;s++) {
    PIC::Mover::MoveParticleTimeStep[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder; ///UniformWeight_UniformTimeStep_SecondOrder;
    PIC::Mover::MoveParticleBoundaryInjection[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder;
  }
*/

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //set up the particle weight
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=20*PIC::nTotalThreads;
  if (PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber>5000) PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=5000;

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=NULL;
//  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_DUST_SPEC_);
  for (int s=0;s<PIC::nTotalSpecies;s++)  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);


  //setup the dust charging functions
  for (int s=_DUST_SPEC_;s<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups;s++) {
    PIC::ChemicalReactions::GenericParticleTranformation::SetSpeciesModel(ElectricallyChargedDust::DustChargingProcessorIndicator,ElectricallyChargedDust::DustChargingProcessor_Implicit_SecondOrder,s);
  }


//  PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
  PIC::Mesh::mesh.outputMeshDataTECPLOT("mesh.data.dat",_DUST_SPEC_);

  //  MPI_Barrier(MPI_COMM_WORLD);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;





  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=NULL;
  PIC::BC::InitBoundingBoxInjectionBlockList();


  //init the particle buffer
  PIC::ParticleBuffer::Init(2000000);

//  double TimeCounter=time(NULL);
  int LastDataOutputFileNumber=-1;



  //dust grains size distribution function:
  //the distribution function is sampled along E7 trajectory at (-100,0,100) sec from CA
  const int nSamplingPoints=3;
//  double ProbeLocations[nSamplingPoints][DIM]={{0.0,0.0,-5.0E5}};

  double ProbeLocations[nSamplingPoints][DIM]={
  {700.835E3,268.886E3,-358.365E3},
  {-20.8542E3,-10.7278E3,-346.632E3},
  {-745.461E3,-282.658E3,-334.544E3},
  };

  ElectricallyChargedDust::Sampling::SampleSizeDistributionFucntion::Init(ProbeLocations,nSamplingPoints,200);





//============================== DEBUG =========================
  /* cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * */  node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  while (node!=NULL) {
    if (node->Temp_ID==2478) {
      cout << __FILE__<< "@" << __LINE__ << endl;


    }

    node=node->nextNodeThisThread;
  }

//============================== END DEBUG ====================================


/*-------------------------  ICES --------------------------------*/
  //Load the plasma parameters from ICES
  //init ICES

  PIC::CPLR::ICES::SWMFdataPreProcessor=SWMFdataPreProcessor;
#ifdef _ICES_CREATE_COORDINATE_LIST_
  PIC::CPLR::ICES::createCellCenterCoordinateList();
  PIC::CPLR::ICES::SetLocationICES("/left/ices/ICES");
  PIC::CPLR::ICES::retriveSWMFdata("Enceladus");
#endif

  PIC::CPLR::ICES::readSWMFdata(1.0);
//  PIC::CPLR::ICES::readDSMCdata();
      cout << __FILE__<< "@" << __LINE__ << endl;


  PIC::Mesh::mesh.outputMeshDataTECPLOT("ices.data.dat",0);
/*-------------------------  END ICES --------------------------------*/

  //time step
  for (long int niter=0;niter<100000001;niter++) {



     PIC::TimeStep();






     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh.ThisThread==0) cout << "The new lample length is " << PIC::RequiredSampleLength << endl;
     }


  }


  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}





