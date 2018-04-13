
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
long int Earth::CutoffRigidity::IntegratedPathLengthOffset=-1;
double*** Earth::CutoffRigidity::CutoffRigidityTable=NULL;

//enable/disable the particle injection procedure
bool Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=true;

//the number of the locations where the rigidity and particle flux is simulated
int Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength=0;
double** Earth::CutoffRigidity::IndividualLocations::xTestLocationTable=NULL;
double Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit=0.0;
double Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit=0.0;
double **Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable=NULL;

double CutoffRigidityTestLocationTable[][3]={{0.0,0.0,0.0}};
int CutoffRigidityTestLocationTableLength=0;

long int Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex=-1;
long int Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed=-1;


int Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations=0;  //the total number of model particles ejected from a test location
int Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations=1;

void Earth::CutoffRigidity::Init_BeforeParser() {
  if ((SampleRigidityMode==true)&&(InitialRigidityOffset==-1)) {
    //request data for sampling of the cutoff rigidity in the particle state vector
    PIC::ParticleBuffer::RequestDataStorage(InitialRigidityOffset,sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(InitialLocationOffset,3*sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(IntegratedPathLengthOffset,sizeof(double));
  }

  //allocate 'IndividualLocations::xTestLocationTable'
  if (CutoffRigidityTestLocationTableLength!=0) {
    //de allocate previously allocated buffers
    if (IndividualLocations::xTestLocationTable!=NULL) {
      delete [] IndividualLocations::xTestLocationTable[0];
      delete [] IndividualLocations::xTestLocationTable;

      delete [] IndividualLocations::CutoffRigidityTable[0];
      delete [] IndividualLocations::CutoffRigidityTable;
    }


    //allocate new sampling buffers
    IndividualLocations::xTestLocationTableLength=CutoffRigidityTestLocationTableLength;

    IndividualLocations::xTestLocationTable=new double* [CutoffRigidityTestLocationTableLength];
    IndividualLocations::xTestLocationTable[0]=new double [3*CutoffRigidityTestLocationTableLength];

    for (int i=1;i<CutoffRigidityTestLocationTableLength;i++) IndividualLocations::xTestLocationTable[i]=IndividualLocations::xTestLocationTable[i-1]+3;

    IndividualLocations::CutoffRigidityTable=new double* [PIC::nTotalSpecies];
    IndividualLocations::CutoffRigidityTable[0]=new double [PIC::nTotalSpecies*CutoffRigidityTestLocationTableLength];

    for (int spec=1;spec<PIC::nTotalSpecies;spec++) IndividualLocations::CutoffRigidityTable[spec]=IndividualLocations::CutoffRigidityTable[spec-1]+CutoffRigidityTestLocationTableLength;


    for (int i=0;i<CutoffRigidityTestLocationTableLength;i++) {
      for (int j=0;j<3;j++) IndividualLocations::xTestLocationTable[i][j]=CutoffRigidityTestLocationTable[i][j];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) IndividualLocations::CutoffRigidityTable[spec][i]=-1.0;
    }

    //request particle storage if needed
    if (Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex==-1) {
      PIC::ParticleBuffer::RequestDataStorage(Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex,sizeof(int));
      PIC::ParticleBuffer::RequestDataStorage(Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed,sizeof(double));
    }
  }
}

void Earth::CutoffRigidity::AllocateCutoffRigidityTable() {
  int i,j,k,offset;

  if (CutoffRigidityTestLocationTableLength!=0) {
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


//process model particles that leaves the computational domain
int Earth::CutoffRigidity::ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double *x,Rigidity;
  long int iAzimuth,iZenith;
  int spec;
  PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  spec=PIC::ParticleBuffer::GetI(ParticleData);

  //update the rigidity data
  if (SampleRigidityMode==true) {
    x=(double*)(ParticleData+InitialLocationOffset);

    Rigidity=*((double*)(ParticleData+InitialRigidityOffset));

    //get coordinates of the point of origin of the particle
    if (Earth::Planet!=NULL) {
      Earth::Planet->GetSurfaceElementProjectionIndex(x,iZenith,iAzimuth);
      if ((CutoffRigidityTable[spec][iZenith][iAzimuth]<0.0)||(CutoffRigidityTable[spec][iZenith][iAzimuth]>Rigidity)) CutoffRigidityTable[spec][iZenith][iAzimuth]=Rigidity;
    }
  }

  //register the particle velocity vector
  int iOriginIndex=0;

  if (Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0) {
    iOriginIndex=*((int*)(ParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex));

    DomainBoundaryParticleProperty::RegisterParticleProperties(PIC::ParticleBuffer::GetI(ptr),xInit,vInit,iOriginIndex,nIntersectionFace);

    if (SampleRigidityMode==true) {
      //sample the cutoff rigidity
      if ((IndividualLocations::CutoffRigidityTable[spec][iOriginIndex]<0.0)||(Rigidity<IndividualLocations::CutoffRigidityTable[spec][iOriginIndex])) {
        IndividualLocations::CutoffRigidityTable[spec][iOriginIndex]=Rigidity;
      }
    }
  }


  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//output sampled cutoff rigidity map
void Earth::CutoffRigidity::OutputDataFile::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Cutoff Rigidity\"");
}

void Earth::CutoffRigidity::OutputDataFile::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement,iZenith,iAzimuth;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double CutoffRigidity=0.0,SurfaceElementCutoffRigidity;

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
double Earth::CutoffRigidity::ParticleInjector::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return 1.0E14;}

//generate a new particle
bool Earth::CutoffRigidity::ParticleInjector::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
  double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,
  double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,int BoundaryElementType,void *BoundaryElement) {

  int idim;
  double ExternalNormal[3];

  //if the injection model is disabled => exit from the procedure
  if (ParticleInjectionMode==false) return false;

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

  static double logRigidityTestMinEnergy=log(RigidityTestMinEnergy);
  static double logRigidityTestMaxEnergy=log(RigidityTestMaxEnergy);

  energy=exp(logRigidityTestMinEnergy+rnd()*(logRigidityTestMaxEnergy-logRigidityTestMinEnergy));

  speed=Relativistic::E2Speed(energy,mass);
  momentum=Relativistic::Speed2Momentum(speed,mass);

  rigidity=(charge>0.0) ? momentum/charge : 0.0;

  for (idim=0;idim<3;idim++) {
    v_IAU_OBJECT[idim]*=speed;
    v_SO_OBJECT[idim]=v_IAU_OBJECT[idim];
  }

  //save the initial location, and rigidity of the particle
  if (SampleRigidityMode==true) {
    *((double*)(tempParticleData+InitialRigidityOffset))=rigidity;
    memcpy((tempParticleData+InitialLocationOffset),v_IAU_OBJECT,3*sizeof(double));
  }

  return true;
}

