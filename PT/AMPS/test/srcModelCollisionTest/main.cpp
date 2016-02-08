//$Id$
//modeling of the collsiion procedure, and the models of the internal degrees of freedom

#include "pic.h"

//parametes of the test model
namespace H2O {
  const double Density=1.0E18;
  const double Temperature=200;
}

namespace O {
  const double Density=1.0E10;
  const double Temperature=400.0;
}

namespace H2 {
  const double Density=1.0E10;
  const double Temperature=800.0;
}

//the desired number of the model particles per cell
const int nParticlePerCell=25;

//size of the domain and cell resolution
const double DomainLength=10.0;
const double dxDomain=1.0;

//functions that returns the local resolution and time step
double localResolution(double *x) {
  return dxDomain;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CollFreqPerParticle=1.0E10*400.0*1.0E-20;

  return 1.0/CollFreqPerParticle;
}

//distribute the blocks between processors
double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

//initialize AMPS
void amps_init() {
  PIC::InitMPI();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //init the particle solver
  PIC::Init_BeforeParser();

  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  double xmax[3]={DomainLength,DomainLength,DomainLength};
  double xmin[3]={-DomainLength,-DomainLength,-DomainLength};
  int idim;


  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,localResolution);
  PIC::Mesh::mesh.memoryAllocationReport();

  if (PIC::Mesh::mesh.ThisThread==0) {
    PIC::Mesh::mesh.buildMesh();
    PIC::Mesh::mesh.saveMeshFile("mesh.msh");
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh.readMeshFile("mesh.msh");
  }

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
  PIC::Init_AfterParser();
  PIC::Mover::Init();


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //init particle weight
  for (int s=0;s<PIC::nTotalSpecies;s++) {
    double weight;

    weight=H2O::Density*pow(2.0*DomainLength/dxDomain,3)/nParticlePerCell;
    PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,weight,PIC::Mesh::mesh.rootTree);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
}


int main(int argc,char **argv) {
  amps_init();

  //generate the new population of the model particles
  double v[3]={0.0,0.0,0.0};

  PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,H2O::Temperature);
  PIC::InitialCondition::PrepopulateDomain(_O_SPEC_,O::Density,v,O::Temperature);
  PIC::InitialCondition::PrepopulateDomain(_H2_SPEC_,H2::Density,v,H2::Temperature);

  PIC::MolecularCollisions::ParticleCollisionModel::ntc();






  //finish execution of the test
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
