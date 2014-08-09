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

double BulletLocalResolution(double *x) {
  int idim;



//  return  ((fabs(x[0])<100.0)||(x[1]*x[1]+x[2]*x[2]<40.0*40.0)) ? 5.0 : 100.0;

  return 10.0;
}

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace) {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  double res,size;

  size=t->CharacteristicSize();

  if (size<0.02) res=0.02;
  else if (size>1.0) res=1.0;
  else res=size;


  return res;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize;
  double CharacteristicSpeed=1.0E3;

  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double v[3]={2.0e3,000.0,000.0},n=5.0E6,temp=20.0;


  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      if (spec!=_O2_SPEC_) return 0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(n,temp,v,ExternalNormal,_O2_SPEC_);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={2.0e3,0.0,0.0},nNA=5.0E6,tempNA=20.0;

  if (PIC::Mesh::mesh.ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,_O2_SPEC_);

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

  if (spec!=_O2_SPEC_) return 0; //inject only spec=0

  static double vNA[3]={2.0e3,000.0,000.0},nNA=5.0E6,tempNA=20.0;
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




  double Rosetta::GetTotalProduction(int spec,void *BoundaryElement) {
    return 1.0E20;
  }

  double Rosetta::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
    return GetTotalProduction(spec,BoundaryElement);
  }


  bool Rosetta::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement) {
//  bool Rosetta::GenerateParticleProperties(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, char *tempParticleData,int BoundaryElementType,void *BoundaryElement) {

    static long int ProbabilityTableLengh=0;
    static double ProbabilityTableIncrement=0.0;

    class cFaceDescriptor {
    public:
      double weight;
      int nFace;
      cFaceDescriptor* next;

      cFaceDescriptor() {
        weight=0.0,nFace=-1,next=NULL;
      }
    };

    class cProbabilityTable {
    public:
      int nTotalFaces;
      cFaceDescriptor *firstFaceDescriptor;

      cProbabilityTable() {
        nTotalFaces=0;
        firstFaceDescriptor=NULL;
      }
    };


    static cFaceDescriptor *FaceDescriptorTable=NULL;
    static cProbabilityTable *ProbabilityTable=NULL;

    static bool initflag=false;

    if (initflag==false) {
      initflag=true;

      double t,minSurfaceArea=-1.0,totalSurfaceArea=0.0;
      int nt,cnt,next;

      const int defaultProbabilityTableLengh=10000;

      //calculate the length of the probability table
      for (nt=0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        if ((minSurfaceArea<0.0)||(minSurfaceArea>CutCell::BoundaryTriangleFaces[nt].SurfaceArea)) minSurfaceArea=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        totalSurfaceArea+=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
      }

      ProbabilityTableIncrement=totalSurfaceArea/defaultProbabilityTableLengh;
      ProbabilityTableLengh=defaultProbabilityTableLengh;

      ProbabilityTable=new cProbabilityTable [ProbabilityTableLengh];

      //calculate the number of the face descriptors that is needed for the mesh
      int nFaceDescriptor=0,start=0,finish=0.0;

      for (nt=0,t=0.0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        t=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        finish=(int)(t/ProbabilityTableIncrement);

        nFaceDescriptor+=finish+1;
      }

      FaceDescriptorTable=new cFaceDescriptor [nFaceDescriptor];

      //init the face descriptor table
      double tstart=0.0,tfinish=0.0;

      for (cnt=0,nt=0,start=0;nt<CutCell::nBoundaryTriangleFaces;nt++) {
        tfinish=tstart+CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
        finish=(int)(tfinish/ProbabilityTableIncrement);

        if (finish>=CutCell::nBoundaryTriangleFaces) finish=CutCell::nBoundaryTriangleFaces-1;

        for (int ii=start;ii<=finish;ii++) {
          if (start==finish) FaceDescriptorTable[cnt].weight=CutCell::BoundaryTriangleFaces[nt].SurfaceArea;
          else {
        	if (ii==start) FaceDescriptorTable[cnt].weight=(start+1)*ProbabilityTableIncrement-tstart;
        	else if (ii==finish) FaceDescriptorTable[cnt].weight=tfinish-finish*ProbabilityTableIncrement;
        	else FaceDescriptorTable[cnt].weight=ProbabilityTableIncrement;
          }

          FaceDescriptorTable[cnt].nFace=nt;
          FaceDescriptorTable[cnt].next=ProbabilityTable[ii].firstFaceDescriptor;

          ProbabilityTable[ii].firstFaceDescriptor=FaceDescriptorTable+cnt;
          ProbabilityTable[ii].nTotalFaces++;

          cnt++;
        }

        start=finish;
        tstart=tfinish;
      }

      //normalize weights
      for (int np=0;np<ProbabilityTableLengh;np++) if (ProbabilityTable[np].firstFaceDescriptor!=NULL) {
        double summ=0.0;
        cFaceDescriptor *face,*prev;

        for (face=ProbabilityTable[np].firstFaceDescriptor;face!=NULL;face=face->next) summ+=face->weight;
        for (face=ProbabilityTable[np].firstFaceDescriptor;face!=NULL;face=face->next) face->weight/=summ;

        //Convert the weight into a cumulative distribution
        prev=ProbabilityTable[np].firstFaceDescriptor;
        face=prev->next;

        for (;face!=NULL;face=face->next) {
          face->weight+=prev->weight;
          prev=face;
        }
      }
    }


    //Determine the face number
    int nt,nface=-1;
    double x[3],v[3];
    double weight=rnd();
    cFaceDescriptor *face;

    nt=(int)(rnd()*ProbabilityTableLengh);

    for (face=ProbabilityTable[nt].firstFaceDescriptor;face!=NULL;face=face->next) if (face->weight>=weight){
	  nface=face->nFace;
	  break;
    }

    bool PositionGenerated;

    do {
    	PositionGenerated=true;

		CutCell::BoundaryTriangleFaces[nface].GetRandomPosition(x);

		//place the point inside the domain
		for (int idim=0;idim<3;idim++) x[idim]+=1.2*PIC::Mesh::mesh.EPS*CutCell::BoundaryTriangleFaces[nface].ExternalNormal[idim];


		//check if the point is inside the domain
		if (CutCell::CheckPointInsideDomain(x,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,PIC::Mesh::mesh.EPS)==false) {
			//exit(__LINE__,__FILE__,"The point is outside of the domain");
			PositionGenerated=false;
		}
    }
    while (PositionGenerated==false);


    //determine if the particle belongs to this processor
    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

    //get the velocity vector
    for (int idim=0;idim<3;idim++) v[idim]=500.0*CutCell::BoundaryTriangleFaces[nface].ExternalNormal[idim];

    memcpy(x_SO_OBJECT,x,3*sizeof(double));
    memcpy(x_IAU_OBJECT,x,3*sizeof(double));
    memcpy(v_SO_OBJECT,v,3*sizeof(double));
    memcpy(v_IAU_OBJECT,v,3*sizeof(double));

    return true;

  }



int main(int argc,char **argv) {

  //init the particle solver
  PIC::InitMPI();
  PIC::Init_BeforeParser();

  Rosetta::Init_BeforeParser();


  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};



  //load the NASTRAN mesh
  CutCell::ReadNastranSurfaceMeshLongFormat("rosetta.surface.reduced.nas",xmin,xmax,1.0E-8);
  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);
    CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);
  }


  //set up the s/c object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor RosettaSurfaceDescriptor;
  cInternalNastranSurfaceData *Rosetta;

//  PIC::BC::InternalBoundary::RotationBody::Init(ReserveSamplingSpace,NULL);
  RosettaSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  Rosetta=(cInternalNastranSurfaceData*) RosettaSurfaceDescriptor.BoundaryElement;

//  Rosetta->PrintTitle=Comet::Sampling::OutputSurfaceDataFile::PrintTitle;
//  Rosetta->PrintVariableList=Comet::Sampling::OutputSurfaceDataFile::PrintVariableList;

//  Nucleus->localResolution=Comet::localSphericalSurfaceResolution;
  Rosetta->InjectionRate=Rosetta::GetTotalProduction;
  Rosetta->faceat=0;
//  Nucleus->ParticleSphereInteraction=Comet::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
  Rosetta->InjectionBoundaryCondition=Rosetta::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;
//  Nucleus->InjectionBoundaryCondition=Comet::InjectionBoundaryModel_Limited; ///sphereParticleInjection;
  //PIC::BC::UserDefinedParticleInjectionFunction=Comet::InjectionBoundaryModel_Limited;




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


  //test the shadow procedure
  double xLightSource[3]={200.0,0.0,100.0}; //{6000.0e3,1.5e6,0.0};
  PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,false);

  if (PIC::ThisThread==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/SurfaceTriangulation-shadow.dat",PIC::OutputDataFileDirectory);
    CutCell::PrintSurfaceTriangulationMesh(fname,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,1.0E-8);
  }

  //init the volume of the cells'
  PIC::Mesh::mesh.InitCellMeasure();



  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=2000; //0; //00; //*10;
  PIC::RequiredSampleLength=100; //00; //0; //0;


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
