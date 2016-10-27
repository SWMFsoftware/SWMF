//$Id$
//modeling of the collsiion procedure, and the models of the internal degrees of freedom

#include "pic.h"
#include <fstream>

//the desired number of the model particles per cell
const int nCase=256;
int iCase;

//size of the domain and cell resolution
const double dxSubMesh=1.0;
const double dxCell=0.5*dxSubMesh / _BLOCK_CELLS_X_;
double XMin[3] = {            0,             0,             0};
double XMax[3] = {    dxSubMesh,     dxSubMesh,     dxSubMesh};
double XMid[3] = {0.5*dxSubMesh, 0.5*dxSubMesh, 0.5*dxSubMesh};

// refinement level at given location
int get_level(double *x){
  int iLocation = 0;
  for(int i=0; i<3; i++)
    iLocation |= ((x[i] > 0.5*dxSubMesh) ? 1 : 0) << i;
  return (iCase >> iLocation) & 1;
}
//-----------------------------------------------------------------------------
//functions that returns the local resolution and time step
double localResolution(double *x) {
  if(x[0] >= 0.25*dxSubMesh && x[0] < 0.75*dxSubMesh ||
     x[1] >= 0.25*dxSubMesh && x[1] < 0.75*dxSubMesh ||
     x[2] >= 0.25*dxSubMesh && x[2] < 0.75*dxSubMesh)
    return 4*dxCell;
  int iLevel = get_level(x);
  return 4*dxCell / (1 + 3*iLevel);
}
//-----------------------------------------------------------------------------
//distribute the blocks between processors
double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  return 1.0;
}
//-----------------------------------------------------------------------------
//initialize AMPS
void amps_init() {
  static bool FirstCall = true;
  if(FirstCall){
    PIC::InitMPI();
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    rnd_seed(-1);
    
    //init the particle solver
    PIC::Init_BeforeParser();
    FirstCall = false;
  }

  //generate only the tree
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(XMin,XMax,localResolution);
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
  
  char filename[100];
  sprintf(filename,"%s%i%s","mesh.case=",iCase,".dat");
  
  PIC::Mesh::mesh.outputMeshTECPLOT(filename);
  
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
}
//-----------------------------------------------------------------------------
//if x is a ghost cell center, may need to fix its coordinates 
//so they correspond to the proper refinement level
void fix_coords(double* x){
    //find refinement level at the location
  int iLevel = get_level(x);
  if(iLevel == 1)
    return;
  double xShift[3] = {2*(x[0]-XMid[0])/dxCell,
		      2*(x[1]-XMid[1])/dxCell,
		      2*(x[2]-XMid[2])/dxCell};
  if( fabs( fabs(xShift[0]) - (int)fabs(xShift[0]) ) < 0.1)
    return;
  for(int i=0; i<3; i++)
    x[i] = 0.5+(2*floor(0.5*xShift[i]) + 1) * 0.5*dxCell;
}
//-----------------------------------------------------------------------------
int main(int argc,char **argv) {
  // alias
  namespace Interpolation = PIC::InterpolationRoutines::CellCentered;

  //init the interpolation routines
  PIC::InterpolationRoutines::Init();

  //init mesh
  //  amps_init();

  // output file name
  char fname[400];
  sprintf(fname,"%s/test_InterpolateAMR-test.dat",
	  PIC::OutputDataFileDirectory);
  fstream fout;
  bool open_file = true;
  bool fail = false;

  Interpolation::cStencil* Stencil;
  if (PIC::ThisThread==0) for(iCase=0; iCase<nCase && !fail; iCase++){
      amps_init();
    if(open_file){
      fout.open(fname, fstream::out);
      open_file=false;
    }

    int nPoint = 100000;
    double x0[3], x1[3], *xcell, diff;
    for(int iPoint=0; iPoint < nPoint && !fail; iPoint++){
      //      cout <<"\riCase="<<iCase<<" iPoint="<<iPoint<<"       ";
      // generate test point
      for(int i=0; i<3; i++){
	x0[i] = XMin[i] + 0.5*dxCell + rnd() * (XMax[i] - XMin[i] - dxCell);
	x1[i] = 0;
      }
      // get interpolation stencil
      Stencil = Interpolation::Linear::InitStencil(x0,NULL);
      //check order of interpolation
      for(int icell=0; icell<Stencil->Length; icell++){
	xcell = Stencil -> cell[icell] -> GetX();
	fix_coords(xcell);
	for(int i=0; i<3; i++)
	  x1[i] += xcell[i] * Stencil->Weight[icell];
      }
      diff=0;
      for(int i=0; i<3; i++)
	diff+= pow(x1[i]-x0[i],2);
      diff = pow(diff, 0.5);
      if(diff > 0.001){
	fout << "iCase="   << iCase
	     << " iPoint=" << iPoint
	     << " x=" << x0[0] << " " << x0[1] << " " << x0[2]
	     << " diff=" << diff << "\nStencil:\n";
	for(int icell=0; icell<Stencil->Length; icell++){
	  xcell = Stencil->cell[icell] -> GetX();
	  fix_coords(xcell);
	  fout << "Cell " << icell << ": "
	       << xcell[0] << " " << xcell[1] << " " << xcell[2]
	       << " Weight=" << Stencil->Weight[icell] << endl;
	}
	fail = true;
      }
    }
  }

  //finish execution of the test
  MPI_Finalize();
  cout << "\nEnd of the run" << endl;

  fout.close();
  return EXIT_SUCCESS;
}
