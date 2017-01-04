
//$Id$
//functionality for calculating of the cutogg rigidity

#include "pic.h"
#include "Earth.h"
#include "specfunc.h"

/*
 * CutoffRigidity.cpp
 *
 *  Created on: Dec 25, 2016
 *      Author: vtenishe
 */

bool Earth::CutoffRigidity::SampleRigidityMode=false;
long int Earth::CutoffRigidity::InitialRigidityOffset=-1;
long int Earth::CutoffRigidity::InitialLocationOffset=-1;
double*** Earth::CutoffRigidity::CutoffRigidityTable=NULL;


void Earth::CutoffRigidity::Init_BeforeParser() {
  if (SampleRigidityMode==true) {
    //request data for sampling of the cutoff rigidity in the particle state vector
    PIC::ParticleBuffer::RequestDataStorage(InitialRigidityOffset,sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(InitialLocationOffset,3*sizeof(double));
  }
}

void Earth::CutoffRigidity::AllocateCutoffRigidityTable() {
  int i,j,k,offset;

  if (SampleRigidityMode==true) {
    //allocate the cutoff rigidity table
    //access pattern CutoffRigidityTable[spec][iZenith][iAzimuthal]
    CutoffRigidityTable=new double** [PIC::nTotalSpecies];
    CutoffRigidityTable[0]=new double *[PIC::nTotalSpecies*Earth::Planet->nZenithSurfaceElements];
    CutoffRigidityTable[0][0]=new double [PIC::nTotalSpecies*Earth::Planet->nZenithSurfaceElements*Earth::Planet->nAzimuthalSurfaceElements];

    for (i=0,offset=0;i<PIC::nTotalSpecies;i++) {
      CutoffRigidityTable[i]=CutoffRigidityTable[0]+offset;
      offset+=Earth::Planet->nZenithSurfaceElements;
    }

    for (i=0,offset=0;i<PIC::nTotalSpecies;i++) for (j=0;j<Earth::Planet->nZenithSurfaceElements;j++) {
      CutoffRigidityTable[i][j]=CutoffRigidityTable[0][0]+offset;
      offset+=Earth::Planet->nAzimuthalSurfaceElements;
    }

    //init the default value of the cutoff rigidity table
    for (i=0;i<PIC::nTotalSpecies;i++) for (j=0;j<Earth::Planet->nZenithSurfaceElements;j++) for (k=0;k<Earth::Planet->nAzimuthalSurfaceElements;k++) {
      CutoffRigidityTable[i][j][k]=-1.0;
    }
  }
}

int Earth::CutoffRigidity::ReversedTimeRelativisticBoris(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  return PIC::Mover::Relativistic::Boris(ptr,-dtTotal,startNode);
}

//process model particles that leaves the computational domain
int Earth::CutoffRigidity::ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double *x,Rigidity;
  long int iAzimuth,iZenith;
  int spec;
  PIC::ParticleBuffer::byte *ParticleData;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  spec=PIC::ParticleBuffer::GetI(ParticleData);
  x=(double*)(ParticleData+InitialLocationOffset);
  Rigidity=*((double*)(ParticleData+InitialRigidityOffset));

  //get dooedinates of the point of origin of the particle
  Earth::Planet->GetSurfaceElementProjectionIndex(x,iZenith,iAzimuth);

  //update the rigidity data
  if ((CutoffRigidityTable[spec][iZenith][iAzimuth]<0.0)||(CutoffRigidityTable[spec][iZenith][iAzimuth]>Rigidity)) CutoffRigidityTable[spec][iZenith][iAzimuth]=Rigidity;

  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//output sampled cutoff rigidity map
void Earth::CutoffRigidity::OutputDataFile::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Cutoff Rigidity\"");
}

void Earth::CutoffRigidity::OutputDataFile::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement,iZenith,iAzimuth;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double t,CutoffRigidity=0.0,SurfaceElementCutoffRigidity;

  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
    InterpolationCoefficient=Sphere->GetSurfaceElementArea(nSurfaceElement);

    Sphere->GetSurfaceElementIndex(iZenith,iAzimuth,nSurfaceElement);

    if (PIC::ThisThread!=0) {
      pipe->send(CutoffRigidityTable[spec][iZenith][iAzimuth]);
    }
    else {
      double t;
      int thread;

      SurfaceElementCutoffRigidity=CutoffRigidityTable[spec][iZenith][iAzimuth];

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        t=pipe->recv<double>(thread);
        if ((SurfaceElementCutoffRigidity<0.0)||(SurfaceElementCutoffRigidity>t)) SurfaceElementCutoffRigidity=t;
      }

      if (SurfaceElementCutoffRigidity>0.0) {
        CutoffRigidity+=SurfaceElementCutoffRigidity*InterpolationCoefficient;
        InterpolationNormalization+=InterpolationCoefficient;
      }
    }
  }

  if (PIC::ThisThread==0) fprintf(fout," %e ",((InterpolationNormalization>0.0) ? CutoffRigidity/InterpolationNormalization : -1));
}

//injection rate of the test particles when calculate the cut-off rigidity
double Earth::CutoffRigidity::ParticleInjector::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return 1.0;}

//generate a new particle
bool Earth::CutoffRigidity::ParticleInjector::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
                                       double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,
                                       double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,
                                       int BoundaryElementType,void *BoundaryElement) {

  int idim;
  double ExternalNormal[3];


  //Generate new particle position
  Vector3D::Distribution::Uniform(ExternalNormal);

  for (idim=0;idim<DIM;idim++) {
    x_IAU_OBJECT[idim]=-RigidityTestRadiusVector*ExternalNormal[idim];
    x_SO_OBJECT[idim]=x_IAU_OBJECT[idim];
  }

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh.findTreeNode(x_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) return false;

  //generate velocity of the injected particle
  //IMPORTANT! The velocity vector is directed inside the sphere; and the trajectory is traced backward

  //direction of the new particle velocity
  do {
    Vector3D::Distribution::Uniform(v_IAU_OBJECT);
  }
  while (Vector3D::DotProduct(v_IAU_OBJECT,ExternalNormal)<=0.0);

  //energy and rigidity of the new particles
  double mass,speed,energy,rigidity,momentum,charge;

  mass=PIC::MolecularData::GetMass(spec);
  charge=PIC::MolecularData::GetElectricCharge(spec);

  energy=Earth::CutoffRigidity::RigidityTestMinEnergy+rnd()*(RigidityTestMaxEnergy-RigidityTestMinEnergy);
  speed=Relativistic::E2Speed(energy,mass);
  momentum=Relativistic::Speed2Momentum(speed,mass);

  rigidity=(charge>0.0) ? momentum/charge : 0.0;

  for (idim=0;idim<3;idim++) {
    v_IAU_OBJECT[idim]*=speed;
    v_SO_OBJECT[idim]=v_IAU_OBJECT[idim];
  }

  //save the initial location, and rigidity of the particle
  *((double*)(tempParticleData+InitialRigidityOffset))=rigidity;
  memcpy((tempParticleData+InitialLocationOffset),v_IAU_OBJECT,3*sizeof(double));

  return true;
}

