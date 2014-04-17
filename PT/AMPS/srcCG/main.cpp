
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
#include "Comet.h"

double BulletLocalResolution(double *x) {
  int idim,i;
  double res,r,l[3];
  double SubsolarAngle;



  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

    //3.3 AU
  if (r<1600) return 300;

  SubsolarAngle=acos(l[0])*180.0/Pi;
  if (SubsolarAngle>=89.5) return 2000.0*pow(r/1600.0,2.0);

 //  return 300*pow(r/1600.0,2.0)*(1+5*SubsolarAngle/180.0);

 return 300*pow(r/1600.0,2.0);

 /*  //2.7 AU
  if (r<1600) return 300;

  SubsolarAngle=acos(l[0])*180.0/Pi;
  if (SubsolarAngle>=89.5) return 700.0*pow(r/1600.0,2.0);

 //  return 300*pow(r/1600.0,2.0)*(1+5*SubsolarAngle/180.0);

 return 70*pow(r/1600.0,2.0);
 */

//  return 1200.0; //Hartley 2
//  return 200.0;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace) {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  return max(1.0,t->CharacteristicSize());
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed=4.0E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  /*
  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double v[3]={000.0,2.0e3,000.0},n=5.0E6,temp=20.0;


  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (spec!=_O2_SPEC_) return 0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,_O2_SPEC_);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;*/
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  /*  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={000.0,2.0e3,0.0},nNA=5.0E6,tempNA=20.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O2_SPEC_);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }

  return false;*/
}

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  /*  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;

  if (spec!=_O2_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={000.0,2.0e3,000.0},nNA=5.0E6,tempNA=20.0;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O2_SPEC_);


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
          PIC::Distribution::InjectMaxwellianDistribution(v,vNA,tempNA,ExternalNormal,_O2_SPEC_,-1);

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

  return nInjectedParticles;*/
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  /*  long int nInjectedParticles=0;

  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;*/
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



  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  //  CutCell::ReadNastranSurfaceMeshLongFormat("surface_Thomas_elements.nas",CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,xmin,xmax,1.0E-8);
  CutCell::ReadNastranSurfaceMeshLongFormat("cg.Lamy-surface.nas",CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,xmin,xmax,1.0E-8);
  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);
  }

  //Computation of the gravity field for an irregular nucleus shape
  nucleusGravity::readMesh("cg.Lamy.nas");
  nucleusGravity::setDensity(300);
  

  //  for (int i=0;i<3;i++) xmin[i]*=6.0,xmax[i]*=6.0;
  for (int i=0;i<3;i++) xmin[i]=-1.0e4,xmax[i]=1.0e4;


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


  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=300; //0; //00; //*10;
  PIC::RequiredSampleLength=10; //00; //0; //0;


  PIC::Init_AfterParser ();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H2O_SPEC_);

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  //init the particle buffer
  PIC::ParticleBuffer::Init(1000000);

  //set the model of the boundary conditinos
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceBoundaryCondition;

  //set the User Definted injection function
  PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;


  //output the volume mesh
  char fname[_MAX_STRING_LENGTH_PIC_];
  sprintf(fname,"%s/VolumeMesh.dat",PIC::OutputDataFileDirectory);
  PIC::Mesh::mesh.outputMeshTECPLOT(fname);


  for (long int niter=0;niter<100000001;niter++) {
    PIC::TimeStep();
  }

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return 1;
}
