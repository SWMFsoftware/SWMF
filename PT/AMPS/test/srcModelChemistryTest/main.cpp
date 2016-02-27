//$Id$
//modeling of the collsiion procedure, and the models of the internal degrees of freedom

#include "pic.h"
#include "PhotolyticReactions.h"
#include "ElectronImpact.h"

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
const double DomainLength=10000.0;
const double dxDomain=1000.0;

double TimeStepMultiplierTable[PIC::nTotalSpecies]={1.0,1.0,1.0};

//functions that returns the local resolution and time step
double localResolution(double *x) {
  return dxDomain;
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  return TimeStepMultiplierTable[spec]*startNode->GetCharacteristicCellSize()/1500.0;
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

//////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////
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

    density=(s==_H2O_SPEC_) ? H2O::Density : 1.0E-2*H2O::Density;

//    density=H2O::Density;

    weight=density*pow(2.0*DomainLength,3)/(GetTotalCellNumber(PIC::Mesh::mesh.rootTree)*nParticlePerCell);
    PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,weight,PIC::Mesh::mesh.rootTree);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;

  //init the particle buffer
  PIC::ParticleBuffer::Init(10000000);
}



//================================================================================
//delete all particles
void CountParticles(double *ParticleCouter,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int i,j,k,FirstCellParticle;
  double v0[3],v1[3],cr;
  double static v[3]={0.0,0.0,0.0};
  long int p;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    //evaluate the mean relative speed
    if (startNode->block!=NULL) {

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              FirstCellParticle=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              for (p=FirstCellParticle;p!=-1;p=PIC::ParticleBuffer::GetNext(p)) {
                if (PIC::ParticleBuffer::IsParticleAllocated(p)==false) exit(__LINE__,__FILE__,"Error: the particle is re-deleted");

                ParticleCouter[PIC::ParticleBuffer::GetI(p)]+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
              }
           }
         }
      }
    }

  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      CountParticles(ParticleCouter,downNode);
    }
  }

}


//================================================================================================
//reflection of a particle from the boundaries of the compukrational domain
int ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int idim;
  double c=0.0,e[3]={0.0,0.0,0.0};

  //determine the normal to the face of intersection
  switch (nIntersectionFace) {
  case 0:case 1:
    e[0]=1.0;
    break;
  case 2:case 3:
    e[1]=1.0;
    break;
  case 4:case 5:
    e[2]=1.0;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the case number is wrong");
  }

  for (idim=0;idim<3;idim++) c+=vInit[idim]*e[idim];
  for (idim=0;idim<3;idim++) vInit[idim]-=2.0*c*e[idim];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//================================================================================================
//The photolitic reaction rate and reaction processor
double h2oTheoreticalLifeTime=0.0;
double ProductionYieldTable[PIC::nTotalSpecies][PIC::nTotalSpecies];

double TheoreticalLifeTime(double *x,int spec,long int ptr,bool &PhotolyticReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double res=0.0;

  PhotolyticReactionAllowedFlag=(spec==_H2O_SPEC_) ? true : false;

  switch (spec) {
  case _H2O_SPEC_:
    res=-PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_]/log(0.1);
    h2oTheoreticalLifeTime=res;
    break;
  }

  return res;
}

int ReactionProcessor(double *xInit,double *xFinal,double *vFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData, double ReactionTimeInterval,double ReactionTimeIntervalLeft,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int *ReactionProductsList,nReactionProducts;
  double *ReactionProductVelocity;
  int ReactionChannel;
  bool PhotolyticReactionRoute;


  //init the reaction tables
  static bool initflag=false;
  static double TotalProductYeld_PhotolyticReaction[PIC::nTotalSpecies*PIC::nTotalSpecies];
  static double TotalProductYeld_ElectronImpact[PIC::nTotalSpecies*PIC::nTotalSpecies];

  double HotElectronFraction=0.05;
  static const double ThermalElectronTemeprature=20.0;
  static const double HotElectronTemeprature=250.0;

  if (initflag==false) {
    int iParent,iProduct;

    initflag=true;

    for (iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
      TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=0.0;
      TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=0.0;

      if (PhotolyticReactions::ModelAvailable(iParent)==true) {
        TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
      }

      if (ElectronImpact::ModelAvailable(iParent)==true) {
        TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=
            ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,HotElectronTemeprature);
      }

      ProductionYieldTable[iParent][iProduct]=TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]+
          TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies];
    }
  }

  //determine the type of the reaction
  double PhotolyticReactionRate=1.0,ElectronImpactRate=0.0;


  PhotolyticReactionRoute=(rnd()<PhotolyticReactionRate/(PhotolyticReactionRate+ElectronImpactRate)) ? true : false;

  //inject the products of the reaction
  double ParentTimeStep,ParentParticleWeight;

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  ParentParticleWeight=0.0;
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
  ParentTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif


  //account for the parent particle correction factor
  ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

  //the particle buffer used to set-up the new particle data
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  //copy the state of the initial parent particle into the new-daugher particle (just in case....)
  PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*)tempParticleData,ParticleData);

  for (int specProduct=0;specProduct<PIC::nTotalSpecies;specProduct++) {
    double ProductTimeStep,ProductParticleWeight;
    double ModelParticleInjectionRate,TimeCounter=0.0,TimeIncrement,ProductWeightCorrection=1.0;
    int iProduct;
    long int newParticle;
    PIC::ParticleBuffer::byte *newParticleData;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
     ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[specProduct];
#else
     ProductParticleWeight=0.0;
     exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
     ProductTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[specProduct];
#else
     ProductTimeStep=0.0;
     exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

//     ModelParticleInjectionRate=ParentParticleWeight/ParentTimeStep/ProductParticleWeight*((PhotolyticReactionRoute==true) ? TotalProductYeld_PhotolyticReaction[specProduct+spec*PIC::nTotalSpecies] : TotalProductYeld_ElectronImpact[specProduct+spec*PIC::nTotalSpecies]);

     ModelParticleInjectionRate=ParentParticleWeight/h2oTheoreticalLifeTime/ProductParticleWeight*((PhotolyticReactionRoute==true) ? TotalProductYeld_PhotolyticReaction[specProduct+spec*PIC::nTotalSpecies] : TotalProductYeld_ElectronImpact[specProduct+spec*PIC::nTotalSpecies]);

     //inject the product particles
     if (ModelParticleInjectionRate>0.0) {
       TimeIncrement=-log(rnd())/ModelParticleInjectionRate *rnd(); //<- *rnd() is to account for the injection of the first particle in the curent interaction

       //the time interval during which the nre partilces will be injected
       double ProductInjectionTimeInterval,ProductTimeIntervalLeft;

       ProductInjectionTimeInterval=ReactionTimeInterval*ProductTimeStep/ParentTimeStep;
       ProductTimeIntervalLeft=ReactionTimeIntervalLeft*ProductTimeStep/ParentTimeStep;

       TimeCounter=ProductTimeStep-ProductTimeIntervalLeft;

       while (TimeCounter+TimeIncrement<ProductTimeStep) {
         TimeCounter+=TimeIncrement;
         TimeIncrement=-log(rnd())/ModelParticleInjectionRate;

         //generate model particle with spec=specProduct
         bool flag=false;

         do {
           //generate a reaction channel
           if (PhotolyticReactionRoute==true) {
             PhotolyticReactions::GenerateReactionProducts(spec,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }
           else {
             ElectronImpact::GenerateReactionProducts(spec,HotElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }

           //check whether the products contain species with spec=specProduct
           for (iProduct=0;iProduct<nReactionProducts;iProduct++) if (ReactionProductsList[iProduct]==specProduct) {
             flag=true;
             break;
           }
         }
         while (flag==false);


         //determine the velocity of the product specie
         double x[3],v[3],c=rnd();

         for (int idim=0;idim<3;idim++) {
           x[idim]=xInit[idim]+c*(xFinal[idim]-xInit[idim]);
           v[idim]=vFinal[idim]+ReactionProductVelocity[idim+3*iProduct];
         }

         //generate a particle
         PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetI(specProduct,(PIC::ParticleBuffer::byte*)tempParticleData);

         #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
         #endif

         //apply condition of tracking the particle
         #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
         PIC::ParticleTracker::InitParticleID(tempParticleData);
         PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,specProduct,tempParticleData);
         #endif


         //get and injection into the system the new model particle
         newParticle=PIC::ParticleBuffer::GetNewParticle();
         newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
         memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);

         node=PIC::Mesh::mesh.findTreeNode(x,node);
         _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,ProductTimeStep-TimeCounter,node);
       }
     }

  }


  return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
}


void PhotochemicalModel(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int *ReactionProductsList,nReactionProducts;
  double *ReactionProductVelocity;
  int ReactionChannel,spec;
  bool PhotolyticReactionRoute;
  PIC::ParticleBuffer::byte *ParticleData;
  double vParent[3],xParent[3];

  //get the particle data
  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  PIC::ParticleBuffer::GetV(vParent,ParticleData);
  PIC::ParticleBuffer::GetX(xParent,ParticleData);

  bool PhotolyticReactionAllowedFlag;

  TheoreticalLifeTime(NULL,spec,ptr,PhotolyticReactionAllowedFlag,node);

  if (PhotolyticReactionAllowedFlag==false) {
    //no reaction occurs -> add the particle to the list
    //the particle is remain in the system
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
    return;
  }

  //init the reaction tables
  static bool initflag=false;
  static double TotalProductYeld_PhotolyticReaction[PIC::nTotalSpecies*PIC::nTotalSpecies];
  static double TotalProductYeld_ElectronImpact[PIC::nTotalSpecies*PIC::nTotalSpecies];

  double HotElectronFraction=0.05;
  static const double ThermalElectronTemeprature=20.0;
  static const double HotElectronTemeprature=250.0;

  static double MeanLog=0.0;

  if (initflag==false) {
    int iParent,iProduct;

    initflag=true;

    for (iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
      TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=0.0;
      TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=0.0;

      if (PhotolyticReactions::ModelAvailable(iParent)==true) {
        TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
      }

      if (ElectronImpact::ModelAvailable(iParent)==true) {
        TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies]=
            ElectronImpact::GetSpeciesReactionYield(iProduct,iParent,HotElectronTemeprature);
      }

      ProductionYieldTable[iParent][iProduct]=TotalProductYeld_PhotolyticReaction[iProduct+iParent*PIC::nTotalSpecies]+
          TotalProductYeld_ElectronImpact[iProduct+iParent*PIC::nTotalSpecies];
    }


    int ntest=1000000;
    for (int i=0;i<ntest;i++) MeanLog+=log(rnd());
    MeanLog/=ntest;
  }

  //determine the type of the reaction
  double PhotolyticReactionRate=1.0,ElectronImpactRate=0.0;


  PhotolyticReactionRoute=(rnd()<PhotolyticReactionRate/(PhotolyticReactionRate+ElectronImpactRate)) ? true : false;

  //inject the products of the reaction
  double ParentTimeStep,ParentParticleWeight;

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  ParentParticleWeight=0.0;
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  ParentTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#else
  ParentTimeStep=0.0;
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif


  //account for the parent particle correction factor
  ParentParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

  //the particle buffer used to set-up the new particle data
  char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::SetParticleAllocated((PIC::ParticleBuffer::byte*)tempParticleData);

  //copy the state of the initial parent particle into the new-daugher particle (just in case....)
  PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*)tempParticleData,ParticleData);

  for (int specProduct=0;specProduct<PIC::nTotalSpecies;specProduct++) if (specProduct!=spec) {
    double ProductTimeStep,ProductParticleWeight;
    double ModelParticleInjectionRate,TimeCounter=0.0,TimeIncrement,ProductWeightCorrection=1.0;
    int iProduct;
    long int newParticle;
    PIC::ParticleBuffer::byte *newParticleData;


#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
     ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[specProduct];
#else
     ProductParticleWeight=0.0;
     exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif

#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
     ProductTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[specProduct];
#else
     ProductTimeStep=0.0;
     exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

//     ModelParticleInjectionRate=ParentParticleWeight/ParentTimeStep/ProductParticleWeight*((PhotolyticReactionRoute==true) ? TotalProductYeld_PhotolyticReaction[specProduct+spec*PIC::nTotalSpecies] : TotalProductYeld_ElectronImpact[specProduct+spec*PIC::nTotalSpecies]);

     ModelParticleInjectionRate=ParentParticleWeight/h2oTheoreticalLifeTime/ProductParticleWeight*((PhotolyticReactionRoute==true) ? TotalProductYeld_PhotolyticReaction[specProduct+spec*PIC::nTotalSpecies] : TotalProductYeld_ElectronImpact[specProduct+spec*PIC::nTotalSpecies]);


  //   ModelParticleInjectionRate=1.0/h2oTheoreticalLifeTime/MeanLog;



     //inject the product particles
     if (ModelParticleInjectionRate>0.0) {
       TimeIncrement=-log(rnd())/ModelParticleInjectionRate *rnd(); //<- *rnd() is to account for the injection of the first particle in the curent interaction

//       while (TimeCounter+TimeIncrement<ProductTimeStep) {


       double anpart=ProductionYieldTable[spec][specProduct]*(1.0-exp(-ProductTimeStep/h2oTheoreticalLifeTime))*ParentParticleWeight/ProductParticleWeight;
       int npart=(int)anpart;
       if (anpart-npart>rnd()) npart+=1;

       for (int n=0;n<npart;n++) {

         TimeCounter+=TimeIncrement;
         TimeIncrement=-log(rnd())/ModelParticleInjectionRate;

         //generate model particle with spec=specProduct
         bool flag=false;

         do {
           //generate a reaction channel
           if (PhotolyticReactionRoute==true) {
             PhotolyticReactions::GenerateReactionProducts(spec,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }
           else {
             ElectronImpact::GenerateReactionProducts(spec,HotElectronTemeprature,ReactionChannel,nReactionProducts,ReactionProductsList,ReactionProductVelocity);
           }

           //check whether the products contain species with spec=specProduct
           for (iProduct=0;iProduct<nReactionProducts;iProduct++) if (ReactionProductsList[iProduct]==specProduct) {
             flag=true;
             break;
           }
         }
         while (flag==false);


         //determine the velocity of the product specie
         double x[3],v[3],c=rnd();

         for (int idim=0;idim<3;idim++) {
           x[idim]=xParent[idim];
           v[idim]=vParent[idim]+ReactionProductVelocity[idim+3*iProduct];
         }

         //generate a particle
         PIC::ParticleBuffer::SetX(x,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetV(v,(PIC::ParticleBuffer::byte*)tempParticleData);
         PIC::ParticleBuffer::SetI(specProduct,(PIC::ParticleBuffer::byte*)tempParticleData);

         #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
         PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ProductWeightCorrection,(PIC::ParticleBuffer::byte*)tempParticleData);
         #endif

         //apply condition of tracking the particle
         #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
         PIC::ParticleTracker::InitParticleID(tempParticleData);
         PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,specProduct,tempParticleData);
         #endif


         //get and injection into the system the new model particle
         newParticle=PIC::ParticleBuffer::GetNewParticle(FirstParticleCell);
         newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

         PIC::ParticleBuffer::CloneParticle(newParticleData,(PIC::ParticleBuffer::byte *)tempParticleData);

//         memcpy((void*)newParticleData,(void*)tempParticleData,PIC::ParticleBuffer::ParticleDataLength);
       }
     }

  }

  //determine whether the parent particle is removed
  if (rnd()<exp(-ParentTimeStep/h2oTheoreticalLifeTime)) {
    //the particle is remain in the system
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
  else {
    //the particle is removed from the system
    PIC::ParticleBuffer::DeleteParticle(ptr);
  }

}

void PhotochemicalModelWrapper(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int i,j,k;
  double v0[3],v1[3],cr;
  double static v[3]={0.0,0.0,0.0};
  long int p,pnext,oldFirstCellParticle,newFirstCellParticle;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    //evaluate the mean relative speed
    if (startNode->block!=NULL) {

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              oldFirstCellParticle=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
              newFirstCellParticle=-1;

              if (oldFirstCellParticle!=-1) {
                p=oldFirstCellParticle;

                while (p!=-1) {

                static int cnt=0;

                cnt++;

                  if (PIC::ParticleBuffer::IsParticleAllocated(p)==false) exit(__LINE__,__FILE__,"Error: the particle is re-deleted");


                  pnext=PIC::ParticleBuffer::GetNext(p);
                  PhotochemicalModel(p,newFirstCellParticle,startNode);
                  p=pnext;


                  for (long int pp=newFirstCellParticle;pp!=-1;pp=PIC::ParticleBuffer::GetNext(pp)) {
                    if (PIC::ParticleBuffer::IsParticleAllocated(pp)==false) exit(__LINE__,__FILE__,"Error: the particle is re-deleted");


                  }
                }

                startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=newFirstCellParticle;
              }
           }
         }
      }
    }

  }
  else {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      PhotochemicalModelWrapper(downNode);
    }
  }

}


//================================================================================================
//test 1:
//run the photochemistry model directly calling the model processor from PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
void Test1(double ReactionYieldTable[PIC::nTotalSpecies][PIC::nTotalSpecies]) {
  double v[3]={0.0,0.0,0.0};
  char fname[400];
  std::fstream fout;
  int s,n,nTotalTestIterations=100;

  double ProductParticleCounter[PIC::nTotalSpecies],InitialParticleCounter[PIC::nTotalSpecies];
  for (s=0;s<PIC::nTotalSpecies;s++) ProductParticleCounter[s]=0.0,InitialParticleCounter[s]=0.0;

  for (n=0;n<nTotalTestIterations;n++) {
    //populate the domain with partiucles
    PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,H2O::Temperature);
    CountParticles(InitialParticleCounter,PIC::Mesh::mesh.rootTree);

    //apply the chemical model
//    PhotochemicalModelWrapper(PIC::Mesh::mesh.rootTree);
    PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();

    //count and remove all particles
    CountParticles(ProductParticleCounter,PIC::Mesh::mesh.rootTree);
    DeleteAllParticles(PIC::Mesh::mesh.rootTree);
  }

  //cpllecte the data from all processors and determine the lifetile
  double GlobalProductParticleCounter[PIC::nTotalSpecies],GlobalInitialParticleCounter[PIC::nTotalSpecies];

  MPI_Reduce(ProductParticleCounter,GlobalProductParticleCounter,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(InitialParticleCounter,GlobalInitialParticleCounter,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    //H2O loss rate
    double LifeTimeH2O;
    bool h2oRectionAllowedFlag;

    sprintf(fname,"%s/test_Chemistry-Test1.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    LifeTimeH2O=-PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_]/
        log(GlobalProductParticleCounter[_H2O_SPEC_]/GlobalInitialParticleCounter[_H2O_SPEC_]);

    h2oTheoreticalLifeTime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(NULL,_H2O_SPEC_,-1,h2oRectionAllowedFlag,NULL);

    cout << "Test 1: H2O lifetime: Numerical\tTheoretical\tRelative Error\n" << LifeTimeH2O << "\t" << h2oTheoreticalLifeTime << "\t" << fabs(LifeTimeH2O-h2oTheoreticalLifeTime)/max(LifeTimeH2O,h2oTheoreticalLifeTime) << endl << endl;
    fout << "Test 1: H2O lifetime: Numerical\tTheoretical\tRelative Error\n" << LifeTimeH2O << "\t" << h2oTheoreticalLifeTime << "\t" << fabs(LifeTimeH2O-h2oTheoreticalLifeTime)/max(LifeTimeH2O,h2oTheoreticalLifeTime) << endl << endl;

    //get the source rate of the products
    double LossRateH2O,SourceRate;

    LossRateH2O=-(GlobalProductParticleCounter[_H2O_SPEC_]-GlobalInitialParticleCounter[_H2O_SPEC_])/PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_];

    cout << "Test 1: H2O lifetime derived from the reaction daughter product: Numerical\nspec\tH2O LifeTime\n";
    fout << "Test 1: H2O lifetime derived from the reaction daughter product: Numerical\nspec\tH2O LifeTime\n";


    for (s=0;s<PIC::nTotalSpecies;s++) if (s!=_H2O_SPEC_) {
      double ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_];
      double ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[s];

      LifeTimeH2O=-PIC::ParticleWeightTimeStep::GlobalTimeStep[s]/
          log(1.0-GlobalProductParticleCounter[s]/GlobalInitialParticleCounter[_H2O_SPEC_]/ReactionYieldTable[_H2O_SPEC_][s]*ProductParticleWeight/ParentParticleWeight);

      cout << s << "\t" << LifeTimeH2O << endl;
      fout << s << "\t" << LifeTimeH2O << endl;

    }

    cout << "Test 1: Particle Number:\n spec\tInitial Particle Number\tFinal Particle Number\tReaction Yield\n";
    fout << "Test 1: Particle Number:\n spec\tInitial Particle Number\tFinal Particle Number\tReaction Yield\n";

    for (s=0;s<PIC::nTotalSpecies;s++) {
      cout << s << "\t" << GlobalInitialParticleCounter[s] << "\t" << GlobalProductParticleCounter[s] << "\t" << ReactionYieldTable[_H2O_SPEC_][s] << endl;
      fout << s << "\t" << GlobalInitialParticleCounter[s] << "\t" << GlobalProductParticleCounter[s] << "\t" << ReactionYieldTable[_H2O_SPEC_][s] << endl;
    }

    fout.close();
  }
}

//================================================================================================
//test 2:
//run the photochemistry model. Similar to Test1. But the model processor is called from inside of the AMPS' time step procedure
void Test2(double ReactionYieldTable[PIC::nTotalSpecies][PIC::nTotalSpecies]) {
  double v[3]={0.0,0.0,0.0};
  char fname[400];
  std::fstream fout;
  int s,n,nTotalTestIterations=100;

  double ProductParticleCounter[PIC::nTotalSpecies],InitialParticleCounter[PIC::nTotalSpecies];
  for (s=0;s<PIC::nTotalSpecies;s++) ProductParticleCounter[s]=0.0,InitialParticleCounter[s]=0.0;

  for (n=0;n<nTotalTestIterations;n++) {
    //populate the domain with partiucles
    PIC::InitialCondition::PrepopulateDomain(_H2O_SPEC_,H2O::Density,v,H2O::Temperature);
    CountParticles(InitialParticleCounter,PIC::Mesh::mesh.rootTree);

    //apply the chemical model
//    PhotochemicalModelWrapper(PIC::Mesh::mesh.rootTree);
//    PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel();
    PIC::TimeStep();

    //count and remove all particles
    CountParticles(ProductParticleCounter,PIC::Mesh::mesh.rootTree);
    DeleteAllParticles(PIC::Mesh::mesh.rootTree);
  }

  //cpllecte the data from all processors and determine the lifetile
  double GlobalProductParticleCounter[PIC::nTotalSpecies],GlobalInitialParticleCounter[PIC::nTotalSpecies];

  MPI_Reduce(ProductParticleCounter,GlobalProductParticleCounter,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(InitialParticleCounter,GlobalInitialParticleCounter,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    //H2O loss rate
    double LifeTimeH2O;
    bool h2oRectionAllowedFlag;

    sprintf(fname,"%s/test_Chemistry-Test2.dat",PIC::OutputDataFileDirectory);
    fout.open(fname,std::fstream::out);

    LifeTimeH2O=-PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_]/
        log(GlobalProductParticleCounter[_H2O_SPEC_]/GlobalInitialParticleCounter[_H2O_SPEC_]);

    h2oTheoreticalLifeTime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(NULL,_H2O_SPEC_,-1,h2oRectionAllowedFlag,NULL);

    cout << "Test 2: H2O lifetime: Numerical\tTheoretical\tRelative Error\n" << LifeTimeH2O << "\t" << h2oTheoreticalLifeTime << "\t" << fabs(LifeTimeH2O-h2oTheoreticalLifeTime)/max(LifeTimeH2O,h2oTheoreticalLifeTime) << endl << endl;
    fout << "Test 2: H2O lifetime: Numerical\tTheoretical\tRelative Error\n" << LifeTimeH2O << "\t" << h2oTheoreticalLifeTime << "\t" << fabs(LifeTimeH2O-h2oTheoreticalLifeTime)/max(LifeTimeH2O,h2oTheoreticalLifeTime) << endl << endl;

    //get the source rate of the products
    double LossRateH2O,SourceRate;

    LossRateH2O=-(GlobalProductParticleCounter[_H2O_SPEC_]-GlobalInitialParticleCounter[_H2O_SPEC_])/PIC::ParticleWeightTimeStep::GlobalTimeStep[_H2O_SPEC_];

    cout << "Test 2: H2O lifetime derived from the reaction daughter product: Numerical\nspec\tH2O LifeTime\n";
    fout << "Test 2: H2O lifetime derived from the reaction daughter product: Numerical\nspec\tH2O LifeTime\n";


    for (s=0;s<PIC::nTotalSpecies;s++) if (s!=_H2O_SPEC_) {
      double ParentParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[_H2O_SPEC_];
      double ProductParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[s];

      LifeTimeH2O=-PIC::ParticleWeightTimeStep::GlobalTimeStep[s]/
          log(1.0-GlobalProductParticleCounter[s]/GlobalInitialParticleCounter[_H2O_SPEC_]/ReactionYieldTable[_H2O_SPEC_][s]*ProductParticleWeight/ParentParticleWeight);

      cout << s << "\t" << LifeTimeH2O << endl;
      fout << s << "\t" << LifeTimeH2O << endl;

    }

    cout << "Test 2: Particle Number:\n spec\tInitial Particle Number\tFinal Particle Number\tReaction Yield\n";
    fout << "Test 2: Particle Number:\n spec\tInitial Particle Number\tFinal Particle Number\tReaction Yield\n";

    for (s=0;s<PIC::nTotalSpecies;s++) {
      cout << s << "\t" << GlobalInitialParticleCounter[s] << "\t" << GlobalProductParticleCounter[s] << "\t" << ReactionYieldTable[_H2O_SPEC_][s] << endl;
      fout << s << "\t" << GlobalInitialParticleCounter[s] << "\t" << GlobalProductParticleCounter[s] << "\t" << ReactionYieldTable[_H2O_SPEC_][s] << endl;
    }


    fout.close();
  }
}

//================================================================================================
//test the chemistry model
int main(int argc,char **argv) {
  int s;
  amps_init();

  //set time step different for each species
  TimeStepMultiplierTable[_O_SPEC_]=1.0/2.0;
  TimeStepMultiplierTable[_H2_SPEC_]=1.0/4.0;

  for (s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //set the reflective boundary conditions on the boundary of the computational domain
  PIC::Mover::ProcessOutsideDomainParticles=ProcessOutsideDomainParticles;

  //init the reaction ProductionYieldTable table
  double TotalProductYeld_PhotolyticReaction[PIC::nTotalSpecies][PIC::nTotalSpecies];
  double ProductionYieldTable[PIC::nTotalSpecies][PIC::nTotalSpecies];

  for (int iParent=0;iParent<PIC::nTotalSpecies;iParent++) for (int iProduct=0;iProduct<PIC::nTotalSpecies;iProduct++) {
    TotalProductYeld_PhotolyticReaction[iParent][iProduct]=0.0;

    if (PhotolyticReactions::ModelAvailable(iParent)==true) {
      TotalProductYeld_PhotolyticReaction[iParent][iProduct]=PhotolyticReactions::GetSpeciesReactionYield(iProduct,iParent);
    }

    ProductionYieldTable[iParent][iProduct]=TotalProductYeld_PhotolyticReaction[iParent][iProduct];
  }


  //excute test sequence
  Test1(ProductionYieldTable);
  Test2(ProductionYieldTable);

  //finish execution
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}
