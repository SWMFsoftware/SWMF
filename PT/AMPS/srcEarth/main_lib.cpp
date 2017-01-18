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


#include "pic.h"
#include "constants.h"
#include "Earth.h"

#include "GeopackInterface.h"
#include "T96Interface.h"


const double rSphere=_RADIUS_(_TARGET_);

double xMaxDomain=5; //modeling the vicinity of the planet
double yMaxDomain=5; //the minimum size of the domain in the direction perpendicular to the direction to the sun


double dxMinSphere=0.5,dxMaxSphere=0.5;
double dxMinGlobal=1,dxMaxGlobal=1;


//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
  double res=0.0;
  return res;
}


//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double res,r,l[3] = {1.0,0.0,0.0};


  if ( (strcmp(Earth::Mesh::sign,"0x301020156361a50")!=0)) {
    //test mesh
    return 0.1*_RADIUS_(_TARGET_);
  }
  else {
    return 0.5*_RADIUS_(_TARGET_);
  }





  res=dxMinSphere;
  res/=2.1;

  res*=2.1;
  
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  return 5.5* 2.5*rSphere*res;
  #endif

  return 0.3* 5.5* 2.5*rSphere*res;
}


double localResolution(double *x) {
  double res;

  if ( (strcmp(Earth::Mesh::sign,"0x301020156361a50")!=0)) {
    //test mesh
    res=0.5*_RADIUS_(_TARGET_);
  }
  else {
    res=0.5*_RADIUS_(_TARGET_);
  }

  //  if (strcmp(Earth::Mesh::sign,"new")==0) 
{ // new mesh
    int idim;
    double r=0.0;

    for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

    r=sqrt(r);

    #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
    if (r<0.98*rSphere) res=rSphere;
    else if (r<1.05*rSphere) res=localSphericalSurfaceResolution(x);
    else if (r<2.0*rSphere) res=2.5* rSphere * dxMinGlobal;
    else res=6.0*rSphere*dxMinGlobal*max(1.0+(5.0-1.0)/((6.0-2.0)*_RADIUS_(_TARGET_))*(r-2.0*_RADIUS_(_TARGET_)),1.0);
    #else
    if (r<0.98*rSphere) res=rSphere;
    else if (r<2*1.05*rSphere) res=localSphericalSurfaceResolution(x);
    else if (r<3.0*rSphere) res=2.5* rSphere * dxMinGlobal;
    else res=6.0*rSphere*dxMinGlobal*max(1.0+(5.0-1.0)/((6.0-2.0)*_RADIUS_(_TARGET_))*(r-2.0*_RADIUS_(_TARGET_)),1.0);
    #endif
  }

  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  return 2.5*res;
  #endif

  return 0.3*  2.5*res;
}

//set up the local time step
double localTimeStep(int spec, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double mass,maxSpeed,CellSize=startNode->GetCharacteristicCellSize();
  int nCompositionGroup;

  nCompositionGroup=Earth::CompositionGroupTableIndex[spec];
  maxSpeed=Earth::CompositionGroupTable[nCompositionGroup].GetMaxVelocity(spec);

  
/*  //evaluate the maximum particle speed with the energy limit used in the Earth magnetosphere model
  mass=PIC::MolecularData::GetMass(spec);
  maxSpeed=SpeedOfLight*sqrt(Earth::BoundingBoxInjection::maxEnergy/(Earth::BoundingBoxInjection::maxEnergy+mass*SpeedOfLight*SpeedOfLight));*/

  return 0.2*CellSize/maxSpeed;
}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;
  return res;
}









void amps_init_mesh() {

  //if (strcmp(Earth::Mesh::sign,"new")==0) 
{ //full mesh

    dxMinGlobal=2.2*0.4/2.1,dxMaxGlobal=1;

    //modeling the vicinity of the planet
    xMaxDomain=16 *_RADIUS_(_TARGET_); 
    //the minimum size of the domain in the direction perpendicular
    //to the direction to the sun
    yMaxDomain=16 * _RADIUS_(_TARGET_); 

    dxMinSphere=20E3,dxMaxSphere=100E3;
  }
//  else exit(__LINE__,__FILE__,"Error: unknown option");


 PIC::InitMPI();
 rnd_seed();

 MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


 Earth::Init_BeforeParser();

 //init the particle solver
 PIC::Init_BeforeParser();



 const int InitialSampleLength=600;
 
 //register the sphere
 static const bool SphereInsideDomain=true;
 
 if (SphereInsideDomain) {
   double sx0[3]={0.0,0.0,0.0};
   cInternalBoundaryConditionsDescriptor SphereDescriptor;
   cInternalSphericalData *Sphere;
   
   //reserve memory for sampling of the surface balance of sticking species
   long int ReserveSamplingSpace[PIC::nTotalSpecies];
   for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=0;
   
   cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);
   
   
   
   PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
   SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
   Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
   Sphere->SetSphereGeometricalParameters(sx0,rSphere);
   Sphere->ParticleSphereInteraction=Earth::BC::ParticleSphereInteraction;
   
   Earth::Planet=Sphere;
   
   Sphere->Radius=_RADIUS_(_TARGET_);
   Sphere->PrintSurfaceMesh("Sphere.dat");
   Sphere->PrintSurfaceData("SpheraData.dat",0);
   Sphere->localResolution=localSphericalSurfaceResolution;
   Sphere->faceat=0;
   
   Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);

    if (Earth::CutoffRigidity::SampleRigidityMode==true) {
      Earth::CutoffRigidity::AllocateCutoffRigidityTable();

      Sphere->PrintDataStateVector=Earth::CutoffRigidity::OutputDataFile::PrintDataStateVector;
      Sphere->PrintVariableList=Earth::CutoffRigidity::OutputDataFile::PrintVariableList;
    }
 }
 
 //Init the spherical shells used for sampling of the energetic particle flux
 Earth::Sampling::Init();

 //init the solver
 PIC::Mesh::initCellSamplingDataBuffer();
 
 //init the mesh
 if (PIC::ThisThread==0)
   cout << "Init the mesh" << endl;
 
 int maxBlockCellsnumber,minBlockCellsnumber,idim;
 
 maxBlockCellsnumber=_BLOCK_CELLS_X_;
 if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
 if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);
 
 minBlockCellsnumber=_BLOCK_CELLS_X_;
 if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
 if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);
 
 double xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

 for (idim=0;idim<DIM;idim++) {
   xmax[idim]= 29 * _RADIUS_(_TARGET_);
   xmin[idim]=-29 * _RADIUS_(_TARGET_);
 }
 
 
 //generate only the tree
 PIC::Mesh::mesh.AllowBlockAllocation=false;
 PIC::Mesh::mesh.init(xmin,xmax,localResolution);
 PIC::Mesh::mesh.memoryAllocationReport();
 
 
 char mesh[200]="";
 bool NewMeshGeneratedFlag=false;
 FILE *fmesh=NULL;
 
 sprintf(mesh,"amr.sig=%s.mesh.bin",Earth::Mesh::sign);
 fmesh=fopen(mesh,"r");
 
 if (fmesh!=NULL) {
   fclose(fmesh);
   PIC::Mesh::mesh.readMeshFile(mesh);
 }
 else {
   NewMeshGeneratedFlag=true;
   
   if (PIC::Mesh::mesh.ThisThread==0) {
     PIC::Mesh::mesh.buildMesh();
     PIC::Mesh::mesh.saveMeshFile("mesh.msh");
     MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
   }
   else {
     MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
     PIC::Mesh::mesh.readMeshFile("mesh.msh");
   }
 }
 
 if (NewMeshGeneratedFlag==true) PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
 
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

/*
 //print out the mesh file
 PIC::Mesh::mesh.outputMeshTECPLOT("mesh.dat");
*/
 
 //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
 if (NewMeshGeneratedFlag==true) {
   unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();
   
   if (PIC::Mesh::mesh.ThisThread==0) {
     char command[300];
     
     sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
     system(command);
   }
 }
 
 MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
 
 
 if (PIC::ThisThread==0) cout << "AMPS' Initialization is complete" << endl;
 
 MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
 
 
}
 
 void amps_init() {
   int idim;
   
   //init the PIC solver
   PIC::Init_AfterParser ();
   PIC::Mover::Init();
   
   //set up the sampling routine
   //combine all species into the same distribution function
   vector<int> SpeciesTable;

   for (int s=0;s<PIC::nTotalSpecies;s++) SpeciesTable.push_back(s);

/*   const int nSamplePoints=7;
   double SampleLocations[nSamplePoints][DIM]={{1.8E8,0.0,0.0}, {1.3E8,0.0,0.0}, {9.4E7,0.0,0.0}, {5.7E7,0.0,0.0}, {3.6E7,0.0,0.0}, {2.3E7,0.0,0.0}, {7.3E6,0.0,0.0}};

   PIC::EnergyDistributionSampleRelativistic::eMin=1.0E6;
   PIC::EnergyDistributionSampleRelativistic::eMax=1.0E10;
   PIC::EnergyDistributionSampleRelativistic::AddCombinedCombinedParticleDistributionList(SpeciesTable);

   PIC::EnergyDistributionSampleRelativistic::nSamleLocations=nSamplePoints;
   PIC::EnergyDistributionSampleRelativistic::SamplingLocations=SampleLocations;*/

   PIC::EnergyDistributionSampleRelativistic::AddCombinedCombinedParticleDistributionList(SpeciesTable);

   //set up the time step
   PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
   PIC::ParticleWeightTimeStep::initTimeStep();
   

   //set up the particle weight
   if (Earth::ImpulseSource::Mode==true) {
     Earth::ImpulseSource::InitParticleWeight();
   }
   else {
     PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=Earth::BoundingBoxInjection::InjectionRate;
     for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(s);
   }
   
   MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
   if (PIC::Mesh::mesh.ThisThread==0) cout << "The mesh is generated" << endl;
   
   //output final data
   //create the list of mesh nodes where the injection boundary conditions are applied
   
   if (Earth::ImpulseSource::Mode==true) {
     PIC::BC::UserDefinedParticleInjectionFunction=Earth::ImpulseSource::InjectParticles;
   }
   else {
     PIC::BC::BlockInjectionBCindicatior=Earth::BoundingBoxInjection::InjectionIndicator;
     PIC::BC::userDefinedBoundingBlockInjectionFunction=Earth::BoundingBoxInjection::InjectionProcessor;
     PIC::BC::InitBoundingBoxInjectionBlockList();
   }


   //init the particle buffer
//   PIC::ParticleBuffer::Init(10000000);

   int LastDataOutputFileNumber=-1;
   

   //init the nackground magnetic field
   switch (_PIC_COUPLER_MODE_) {
   case _PIC_COUPLER_MODE__DATAFILE_ :

     if (PIC::CPLR::DATAFILE::BinaryFileExists("EARTH-BATSRUS")==true)  {
       PIC::CPLR::DATAFILE::LoadBinaryFile("EARTH-BATSRUS");
     }
     else if (_PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_) {
       // BATL reader

       if (PIC::CPLR::DATAFILE::BinaryFileExists("EARTH-BATSRUS")==true)  {
         PIC::CPLR::DATAFILE::LoadBinaryFile("EARTH-BATSRUS");
       }
       else {
         // initialize the reader
         #if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_
         PIC::CPLR::DATAFILE::BATSRUS::Init("3d__ful_2_t00000000_n00020000.idl");
         PIC::CPLR::DATAFILE::BATSRUS::LoadDataFile();
         #endif

         //initialize derived data
         if (PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate==true) {
           #if _PIC_COUPLER__INTERPOLATION_MODE_==_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_
           exit(__LINE__,__FILE__,"ERROR: magnetic field gradient can't be computed with 0th order interpolation method");
           #endif

           for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
             PIC::CPLR::DATAFILE::GenerateMagneticFieldGradient(node);
           }

           //Exchange derived data betwenn the boundary nodes
           PIC::Mesh::mesh.ParallelBlockDataExchange();
         }

         PIC::CPLR::DATAFILE::SaveBinaryFile("EARTH-BATSRUS");
       }
     }
     else {
       exit(__LINE__,__FILE__,"ERROR: the background importing procedure is not defined");
     }

     break;
   case _PIC_COUPLER_MODE__T96_ :
     if (PIC::CPLR::DATAFILE::BinaryFileExists("EARTH-T96")==true)  {
       PIC::CPLR::DATAFILE::LoadBinaryFile("EARTH-T96");
     }
     else {
       //calculate the geomegnetic filed
       T96::Init(Exosphere::SimulationStartTimeString,NULL);

       //set the magnetic field;
       //set default == 0 electric field


       class cSetBackgroundMagneticField {
       public:
         void Set(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

           const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
           const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
           const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

           if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
             int ii,S=(kMax-kMin+1)*(jMax-jMin+1)*(iMax-iMin+1);

             if (startNode->block!=NULL) {
               #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
               #pragma omp parallel for schedule(dynamic,1) default (none) shared (PIC::Mesh::mesh,iMin,jMin,kMin,S,PIC::CPLR::DATAFILE::Offset::MagneticField,PIC::CPLR::DATAFILE::Offset::ElectricField,startNode)
               #endif

               for (ii=0;ii<S;ii++) {
                 int i,j,k;
                 double *xNodeMin=startNode->xmin;
                 double *xNodeMax=startNode->xmax;
                 double x[3],B[3],xCell[3];
                 PIC::Mesh::cDataCenterNode *CenterNode;

                 //set the value of the geomagnetic field calculated at the centers of the cells
                 int nd,idim;
                 char *offset;

                 //determine the coordinates of the cell
                 int S1=ii;

                 i=iMin+S1/((kMax-kMin+1)*(jMax-jMin+1));
                 S1=S1%((kMax-kMin+1)*(jMax-jMin+1));

                 j=jMin+S1/(kMax-kMin+1);
                 k=kMin+S1%(kMax-kMin+1);

                 //locate the cell
                 nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
                 if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
                 offset=CenterNode->GetAssociatedDataBufferPointer();

                 //the interpolation location
                 xCell[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i));
                 xCell[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j));
                 xCell[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k));

                 //calculate the geomagnetic field
                 T96::GetMagneticField(B,xCell);

                 //save E and B
                 for (idim=0;idim<3;idim++) {
                   if (PIC::CPLR::DATAFILE::Offset::MagneticField.active==true) {
                     *((double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.offset+idim*sizeof(double)))=B[idim];
                   }

                   if (PIC::CPLR::DATAFILE::Offset::ElectricField.active==true) {
                     *((double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.offset+idim*sizeof(double)))=0.0;
                   }
                 }
               }

             }
           }
           else {
             int i;
             cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

             for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) Set(downNode);
           }


         }
       } SetBackgroundMagneticField;

       SetBackgroundMagneticField.Set(PIC::Mesh::mesh.rootTree);
       PIC::CPLR::DATAFILE::SaveBinaryFile("EARTH-T96");
     }


     break;
   default:
     exit(__LINE__,__FILE__,"Error: the option is unknown");
   }



    //init particle weight of neutral species that primary source is sputtering



//  if (_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE_==_PIC_OUTPUT_MACROSCOPIC_FLOW_DATA_MODE__TECPLOT_ASCII_) {
    PIC::Mesh::mesh.outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);
//  }



}

 //time step

void amps_time_step () {
  
  PIC::TimeStep();
  
}
