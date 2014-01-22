//$Id$


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


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"

 static double totalTime=0;
 static bool updateTotalTime=false;
 static bool computeInjectedBulkVelocity=true;
 static double averageBulkVelocity;
 static double averageBulkVelocity2;

double BulletLocalResolution(double *x) {
  int idim;



//  return  ((fabs(x[0])<100.0)||(x[1]*x[1]+x[2]*x[2]<40.0*40.0)) ? 5.0 : 100.0;

  return 0.01;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace) {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  return max(6.0,t->CharacteristicSize());
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed=4.0E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}

int ExternalBoundaryConditions(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  if (nIntersectionFace==0 || nIntersectionFace==1)   return _PARTICLE_DELETED_ON_THE_FACE_;

  static const int ExternalNormal[6][3]={ 
       {-1.0,0.0,0.0}, {1.0,0.0,0.0}, 
       {0.0,-1.0,0.0}, {0.0,1.0,0.0}, 
       {0.0,0.0,-1.0}, {0.0,0.0,1.0}};

  double c=vInit[0]*ExternalNormal[nIntersectionFace][0]+vInit[1]*ExternalNormal[nIntersectionFace][1]+vInit[2]*ExternalNormal[nIntersectionFace][2];

  vInit[0]-=2.0*c*ExternalNormal[nIntersectionFace][0];
  vInit[1]-=2.0*c*ExternalNormal[nIntersectionFace][1];
  vInit[2]-=2.0*c*ExternalNormal[nIntersectionFace][2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  double gamma=5.0/3.0;
  static double M1=1.4,temp=2.93e2,n=1.0E20,M2,n2,temp2;

  double cSon=sqrt(gamma*Kbol*temp/PIC::MolecularData::GetMass(spec));
  printf("spec=%i cSon=%e gamma=%e mass=%e\n",spec,cSon,gamma,PIC::MolecularData::GetMass(spec));

  static double v[3]={-M1*cSon,000.0,000.0};

  M2=sqrt((pow(M1,2.0)*(gamma-1.0)+2.0)/(2.0*gamma*pow(M1,2.0)-(gamma-1.0)));
  n2=(gamma+1)*pow(M1,2.0)/((gamma-1.0)*pow(M1,2.0)+2.0)*n;
  temp2=(1+(gamma-1.0)/2.0*pow(M1,2.0))*(2.0*gamma/(gamma-1.0)*pow(M1,2.0)-1.0)/(pow(M1,2.0)*(2.0*gamma/(gamma-1.0)+(gamma-1.0)/2.0))*temp;

  static double v2[3]={-M2*cSon,000.0,000.0};

  printf("M1=%e v[0]=%e n1=%e temp=%e \n",M1,v[0],n,temp);
  printf("M2=%e v2[0]=%e n2=%e temp2=%e \n",M2,v2[0],n2,temp2);

  /*
  //Mach 8
  static double v[3]={-1930.73,000.0,000.0},n=1.0E20,temp=2.0e2;
  static double v2[3]={-692.98,000.0,000.0},n2=2.786e20,temp2=temp*20.87;
  */

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      //      if (spec!=_O2_SPEC_) return 0.0;

      if (ExternalNormal[0]>0) {
	ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,0);
      }else if (ExternalNormal[0]<0) {
        ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n2,temp2,v2,ExternalNormal,0);
      }

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  double gamma=5.0/3.0;
  static double M1=1.4,temp=2.93e2,n=1.0E20,M2,n2,temp2;
  
  double cSon=sqrt(gamma*Kbol*temp/PIC::MolecularData::GetMass(0));
 
  static double v[3]={-M1*cSon,000.0,000.0};

  M2=sqrt((pow(M1,2.0)*(gamma-1.0)+2.0)/(2.0*gamma*pow(M1,2.0)-(gamma-1.0)));
  n2=(gamma+1)*pow(M1,2.0)/((gamma-1.0)*pow(M1,2.0)+2.0)*n;
  temp2=(1+(gamma-1.0)/2.0*pow(M1,2.0))*(2.0*gamma/(gamma-1.0)*pow(M1,2.0)-1.0)/(pow(M1,2.0)*(2.0*gamma/(gamma-1.0)+(gamma-1.0)/2.0))*temp;

  static double v2[3]={-M2*cSon,000.0,000.0};

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      if (ExternalNormal[0]>0) {
	ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,0);
      }else if (ExternalNormal[0]<0) {
        ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n2,temp2,v2,ExternalNormal,0);
      }

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

  //  if (spec!=_O2_SPEC_) return 0; //inject only spec=0


  double gamma=5.0/3.0;
  static double M1=1.4,temp=2.93e2,n=1.0E20,M2,n2,temp2;

  double cSon=sqrt(gamma*Kbol*temp/PIC::MolecularData::GetMass(spec));

  static double v[3]={-M1*cSon,000.0,000.0};

  M2=sqrt((pow(M1,2.0)*(gamma-1.0)+2.0)/(2.0*gamma*pow(M1,2.0)-(gamma-1.0)));
  n2=(gamma+1)*pow(M1,2.0)/((gamma-1.0)*pow(M1,2.0)+2.0)*n;
  temp2=(1+(gamma-1.0)/2.0*pow(M1,2.0))*(2.0*gamma/(gamma-1.0)*pow(M1,2.0)-1.0)/(pow(M1,2.0)*(2.0*gamma/(gamma-1.0)+(gamma-1.0)/2.0))*temp;

  static double v2[3]={-M2*cSon,000.0,000.0};

  double ModelParticlesInjectionRate;

  double vel[3];

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);

    if (updateTotalTime==true) {
      totalTime+=LocalTimeStep;
      updateTotalTime=false;
    }

    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      if (ExternalNormal[0]>0) {
	ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,0);
      }else if (ExternalNormal[0]<0) {
        ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n2,temp2,v2,ExternalNormal,0);
      }


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
	  if (ExternalNormal[0]>0) {
	    //	    PIC::Distribution::InjectMaxwellianDistribution(vel,v,temp,ExternalNormal,_O2_SPEC_,-1);
	    PIC::Distribution::InjectMaxwellianDistribution(vel,v,temp,ExternalNormal,0);
          }else if (ExternalNormal[0]<0) {
	    PIC::Distribution::InjectMaxwellianDistribution(vel,v2,temp2,ExternalNormal,0);
          }

	  PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(vel,newParticleData);
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

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int main(int argc,char **argv) {

  //init the particle solver
  PIC::InitMPI();
  PIC::Init_BeforeParser();



  double xmin[3]={-0.23,-0.0025,-0.0025};
  double xmax[3]={0.23,0.0025,0.0025};



  /*    //load the NASTRAN mesh
  CutCell::ReadNastranSurfaceMeshLongFormat("bullet.surface.nas",CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,xmin,xmax,1.0E-8);
  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);
  }
  */

  for (int i=0;i<3;i++) xmin[i]*=2.0,xmax[i]*=2.0;


  PIC::Mesh::mesh.CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.buildMesh();

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();



  //initialize the blocks

  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();


  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();



  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=100000; //0; //00; //*10;
  PIC::RequiredSampleLength=1000; //00; //0; //0;


  PIC::Init_AfterParser ();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_O2_SPEC_);

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  //init the particle buffer
  PIC::ParticleBuffer::Init(1000000);

  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;
  PIC::Mover::ProcessOutsideDomainParticles=ExternalBoundaryConditions;


  char fname[_MAX_STRING_LENGTH_PIC_];
  sprintf(fname,"%s/VolumeMesh.dat",PIC::OutputDataFileDirectory);
  PIC::Mesh::mesh.outputMeshTECPLOT(fname);



  for (long int niter=0;niter<100001;niter++) {
    if (niter>0) updateTotalTime=true;
    PIC::TimeStep();
  }

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
