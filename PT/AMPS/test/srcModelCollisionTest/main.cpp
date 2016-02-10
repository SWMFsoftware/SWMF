//$Id$
//modeling of the collsiion procedure, and the models of the internal degrees of freedom

#include "pic.h"

//parametes of the test model
namespace H2O {
  const double Density=1.0E18;
  const double Temperature=200;
}

namespace O {
  const double Density=1.0E15;
  const double Temperature=400.0;
}

namespace H2 {
  const double Density=1.0E12;
  const double Temperature=800.0;
}

//the desired number of the model particles per cell
const int nParticlePerCell=25;

//size of the domain and cell resolution
const double DomainLength=10.0;
const double dxDomain=1.0;

double TimeStepMultiplierTable[PIC::nTotalSpecies]={1.0,1.0,1.0};

//functions that returns the local resolution and time step
double localResolution(double *x) {
  return dxDomain;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CollFreqPerParticle=1.0E10*400.0*1.0E-20;

  static bool initflag=false;
  static double DensityTable[PIC::nTotalSpecies],SummDensity=0.0;

  if (initflag==false) {
    initflag=true;

    DensityTable[_H2O_SPEC_]=H2O::Density;
    DensityTable[_O_SPEC_]=O::Density;
    DensityTable[_H2_SPEC_]=H2::Density;

    for (int s=0;s<PIC::nTotalSpecies;s++) SummDensity+=DensityTable[s];
  }

  CollFreqPerParticle=SummDensity*400.0*1.0E-18;


  return TimeStepMultiplierTable[spec]/CollFreqPerParticle;
}

//distribute the blocks between processors
double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

//calculate the tital number of the simulation cells
int GetTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int res=0.0;
  int i;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    res+=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      res+=GetTotalCellNumber(downNode);
    }
  }

  return res;
}

//sample relative velocity
void SampleRelativeSpeed(double RelativeSpeed[PIC::nTotalSpecies][PIC::nTotalSpecies],int RelativeSpeedCouter[PIC::nTotalSpecies][PIC::nTotalSpecies],cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int i,j,k,FirstCellParticle,p0,p1,s0,s1;
  double v0[3],v1[3],cr;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    //evaluate the mean relative speed
    if (startNode->block!=NULL) {

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              FirstCellParticle=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              for (p0=FirstCellParticle;p0!=-1;p0=PIC::ParticleBuffer::GetNext(p0)) {
                s0=PIC::ParticleBuffer::GetI(p0);
                PIC::ParticleBuffer::GetV(v0,p0);

                for (p1=PIC::ParticleBuffer::GetNext(p0);p1!=-1;p1=PIC::ParticleBuffer::GetNext(p1)) {
                  s1=PIC::ParticleBuffer::GetI(p1);
                  PIC::ParticleBuffer::GetV(v1,p1);

                  cr=sqrt(pow(v1[0]-v0[0],2)+pow(v1[1]-v0[1],2)+pow(v1[2]-v0[2],2));

                  RelativeSpeed[s0][s1]+=cr;
                  RelativeSpeedCouter[s0][s1]++;

                  RelativeSpeed[s1][s0]+=cr;
                  RelativeSpeedCouter[s1][s0]++;
                }
              }

            }
         }
      }
    }

  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,downNode);
    }
  }
}

//delete all particles
void DeleteAllParticles(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int i,j,k,FirstCellParticle,p,pnext;
  double v0[3],v1[3],cr;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    //evaluate the mean relative speed
    if (startNode->block!=NULL) {

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              FirstCellParticle=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              if (FirstCellParticle!=-1) {
                p=FirstCellParticle;

                do {
                  pnext=PIC::ParticleBuffer::GetNext(p);
                  PIC::ParticleBuffer::DeleteParticle(p);
                  p=pnext;
                }
                while (p!=-1);

                startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
              }
            }
         }
      }
    }

  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      DeleteAllParticles(downNode);
    }
  }

}

//initialize AMPS
void amps_init() {
  PIC::InitMPI();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  rnd_seed(-1);

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
    double weight,density;

    switch(s) {
    case _H2O_SPEC_ :
      density=H2O::Density;
      break;
    case _O_SPEC_:
      density=O::Density;
      break;
    case _H2_SPEC_:
      density=H2::Density;
      break;
    }

    weight=density*pow(2.0*DomainLength,3)/(GetTotalCellNumber(PIC::Mesh::mesh.rootTree)*nParticlePerCell);
    PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,weight,PIC::Mesh::mesh.rootTree);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
}


void GetTotalCollisionFreq(double CollisionFreq[PIC::nTotalSpecies][PIC::nTotalSpecies],cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int i,j,k,s0,s1;
  char* SamplingData;

  if (PIC::Mesh::mesh.rootTree==startNode) {
    for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) CollisionFreq[s0][s1]=0.0;
  }

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      SamplingData=startNode->block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+
        PIC::Mesh::collectingCellSampleDataPointerOffset;

      for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++)  {
        int CollFreqOffset=PIC::MolecularCollisions::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset+
            sizeof(double)*PIC::MolecularCollisions::ParticleCollisionModel::CollsionFrequentcySampling::Offset(s0,s1);

        CollisionFreq[s0][s1]+=(*((double*)(SamplingData+CollFreqOffset))); //+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData)*LocalParticleWeight_s0/LocalTimeStep_s0/cellMeasure*CollisionLimitingFactor; // calculate collision frequency taking into account the collision limiting factor
      }
    }
  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      GetTotalCollisionFreq(CollisionFreq,downNode);
    }
  }
}

int main(int argc,char **argv) {
  int s0,s1;
  amps_init();

  double RelativeSpeed[PIC::nTotalSpecies][PIC::nTotalSpecies];
  int RelativeSpeedCouter[PIC::nTotalSpecies][PIC::nTotalSpecies];
  for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) RelativeSpeed[s0][s1]=0.0,RelativeSpeedCouter[s0][s1]=0;

  //======================= TEST 1 BEGINS: COLLISION FREQUENCY AND GENERATED MEAN RELATIVE VELOCITY ===============================
  //generate the new population of the model particles
  double v[3]={0.0,0.0,0.0};
  char fname[400];
  std::fstream fout;
  int n,nTotalTestIterations=50;

  for (n=0;n<nTotalTestIterations;n++) {
    //populate the domain with partiucles
    PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,H2O::Temperature);
    PIC::InitialCondition::PrepopulateDomain(_O_SPEC_,O::Density,v,O::Temperature);
    PIC::InitialCondition::PrepopulateDomain(_H2_SPEC_,H2::Density,v,H2::Temperature);

    //sample relative speed
    if (n==0) {
      //the relative velocity sampling length is sufficient even for a single iteration
      SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

      sprintf(fname,"%s/test_ModelCollisionTest--test1.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
    }

    //callc particle collision model
    PIC::MolecularCollisions::ParticleCollisionModel::ntc();

    //remove all particles
    DeleteAllParticles(PIC::Mesh::mesh.rootTree);
  }

  //collect the collision frequentcy from all processors and output into a file
  double CollisionFreq[PIC::nTotalSpecies][PIC::nTotalSpecies];

  GetTotalCollisionFreq(CollisionFreq,PIC::Mesh::mesh.rootTree);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/test_ModelCollisionTest--test1-1.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    cout << "Test 1: collision frequency: \ns\t";
    fout << "Test 1: collision frequency: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;

      tLocal=CollisionFreq[s0][s1];
      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nTotalTestIterations/GetTotalCellNumber(PIC::Mesh::mesh.rootTree) << "  ";
        fout << tGlobal/nTotalTestIterations/GetTotalCellNumber(PIC::Mesh::mesh.rootTree) << "  ";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) {
    cout << "\nTest 1: the relative speed generated by the domain prepopulation procedure: \ns\t";
    fout << "\nTest 1: the relative speed generated by the domain prepopulation procedure: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }


  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) fout.close();

  //======================= TEST 1 ENDS: COLLISION FREQUENCY AND GENERATED MEAN RELATIVE VELOCITY ===============================


  //======================= TEST 2 BEGINS: RELATIVE VELOCITY AFTER MULTIPLE COLLISIONS WHEN INITAL TEMEPRATURE OF ALL SPECIES IS THE SAME =======
  double Temp=H2O::Temperature;

  PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,Temp);
  PIC::InitialCondition::PrepopulateDomain(_O_SPEC_,O::Density,v,Temp);
  PIC::InitialCondition::PrepopulateDomain(_H2_SPEC_,H2::Density,v,Temp);

  for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) RelativeSpeed[s0][s1]=0.0,RelativeSpeedCouter[s0][s1]=0;
  SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/test_ModelCollisionTest--test2-1.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    cout << "\nTest 2: the relative before collisions: \ns\t";
    fout << "\nTest 2: the relative before collisions: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  //collision loop
  for (n=0;n<nTotalTestIterations;n++) {
    PIC::MolecularCollisions::ParticleCollisionModel::ntc();
  }

  for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) RelativeSpeed[s0][s1]=0.0,RelativeSpeedCouter[s0][s1]=0;
  SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

  sprintf(fname,"%s/test_ModelCollisionTest--test2.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  if (PIC::ThisThread==0) {
    cout << "\nTest 2: the relative after collisions: \ns\t";
    fout << "\nTest 2: the relative after collisions: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) fout.close();

  //remove all particles
  DeleteAllParticles(PIC::Mesh::mesh.rootTree);

  //======================= TEST 2 ENDS: RELATIVE VELOCITY AFTER MULTIPLE COLLISIONS WHEN INITAL TEMEPRATURE OF ALL SPECIES IS THE SAME =======


  //======================= TEST 3 BEGINS: COLLISION FREQUENCY AND GENERATED MEAN RELATIVE VELOCITY (SPECIES DEPENDENT TIME STEP)  =============
  //generate the new population of the model particles
  TimeStepMultiplierTable[_O_SPEC_]=2.0;
  TimeStepMultiplierTable[_H2_SPEC_]=4.0;

  for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  PIC::ParticleWeightTimeStep::initTimeStep();

  for (n=0;n<nTotalTestIterations;n++) {
    //populate the domain with partiucles
    PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,H2O::Temperature);
    PIC::InitialCondition::PrepopulateDomain(_O_SPEC_,O::Density,v,O::Temperature);
    PIC::InitialCondition::PrepopulateDomain(_H2_SPEC_,H2::Density,v,H2::Temperature);

    //sample relative speed
    if (n==0) {
      //the relative velocity sampling length is sufficient even for a single iteration
      SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

      sprintf(fname,"%s/test_ModelCollisionTest--test3.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
    }

    //callc particle collision model
    PIC::MolecularCollisions::ParticleCollisionModel::ntc();

    //remove all particles
    DeleteAllParticles(PIC::Mesh::mesh.rootTree);
  }

  //collect the collision frequentcy from all processors and output into a file
  GetTotalCollisionFreq(CollisionFreq,PIC::Mesh::mesh.rootTree);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/test_ModelCollisionTest--test3-1.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    cout << "Test 3: collision frequency: \ns\t";
    fout << "Test 3: collision frequency: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;

      tLocal=CollisionFreq[s0][s1];
      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nTotalTestIterations/GetTotalCellNumber(PIC::Mesh::mesh.rootTree) << "  ";
        fout << tGlobal/nTotalTestIterations/GetTotalCellNumber(PIC::Mesh::mesh.rootTree) << "  ";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) {
    cout << "\nTest 3: the relative speed generated by the domain prepopulation procedure: \ns\t";
    fout << "\nTest 3: the relative speed generated by the domain prepopulation procedure: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }


  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) fout.close();

  //======================= TEST 3 ENDS: COLLISION FREQUENCY AND GENERATED MEAN RELATIVE VELOCITY  (SPECIES DEPENDENT TIME STEP) ===============================


  //======================= TEST 4 BEGINS: RELATIVE VELOCITY AFTER MULTIPLE COLLISIONS WHEN INITAL TEMEPRATURE OF ALL SPECIES IS THE SAME(SPECIES DEPENDENT TIME STEP)  =======
  Temp=H2O::Temperature;

  PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,Temp);
  PIC::InitialCondition::PrepopulateDomain(_O_SPEC_,O::Density,v,Temp);
  PIC::InitialCondition::PrepopulateDomain(_H2_SPEC_,H2::Density,v,Temp);

  for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) RelativeSpeed[s0][s1]=0.0,RelativeSpeedCouter[s0][s1]=0;
  SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/test_ModelCollisionTest--test4-1.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    cout << "\nTest 4: the relative before collisions: \ns\t";
    fout << "\nTest 4: the relative before collisions: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  //collision loop
  for (n=0;n<nTotalTestIterations;n++) {
    PIC::MolecularCollisions::ParticleCollisionModel::ntc();
  }

  for (s0=0;s0<PIC::nTotalSpecies;s0++) for (s1=0;s1<PIC::nTotalSpecies;s1++) RelativeSpeed[s0][s1]=0.0,RelativeSpeedCouter[s0][s1]=0;
  SampleRelativeSpeed(RelativeSpeed,RelativeSpeedCouter,PIC::Mesh::mesh.rootTree);

  sprintf(fname,"%s/test_ModelCollisionTest--test4.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  if (PIC::ThisThread==0) {
    cout << "\nTest 4: the relative after collisions: \ns\t";
    fout << "\nTest 4: the relative after collisions: \ns\t";

    for (s0=0;s0<PIC::nTotalSpecies;s0++) {
      cout << "s=" << s0 << "\t";
      fout << "s=" << s0 << "\t";
    }

    cout << endl;
    fout << endl;
  }

  //output the relative speed
  for (s0=0;s0<PIC::nTotalSpecies;s0++) {
    if (PIC::ThisThread==0) cout << s0 << "  ";

    for (s1=0;s1<PIC::nTotalSpecies;s1++) {
      double tLocal,tGlobal;
      double nLocal,nGlobal;

      tLocal=RelativeSpeed[s0][s1];
      nLocal=RelativeSpeedCouter[s0][s1];

      MPI_Reduce(&tLocal,&tGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
      MPI_Reduce(&nLocal,&nGlobal,1,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread==0) {
        cout << tGlobal/nGlobal << "\t";
        fout << tGlobal/nGlobal << "\t";
      }
    }

    if (PIC::ThisThread==0) {
      cout << endl;
      fout << endl;
    }
  }

  if (PIC::ThisThread==0) fout.close();

  //remove all particles
  DeleteAllParticles(PIC::Mesh::mesh.rootTree);

  //======================= TEST 4 ENDS: RELATIVE VELOCITY AFTER MULTIPLE COLLISIONS WHEN INITAL TEMEPRATURE OF ALL SPECIES IS THE SAME (SPECIES DEPENDENT TIME STEP) =======

  //finish execution of the test
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
