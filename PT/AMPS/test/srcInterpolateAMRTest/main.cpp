//$Id$
//modeling of the collsiion procedure, and the models of the internal degrees of freedom

#include "pic.h"
#include <fstream>

//the desired number of the model particles per cell
const int nCase=256;
int iCase;

//size of the domain and cell resolution
const double dxSubMesh=1.0e16;
const double dxCell=0.5*dxSubMesh / _BLOCK_CELLS_X_;
double XMin[3] = {            0,             0,             0};
double XMax[3] = {    dxSubMesh,     dxSubMesh,     dxSubMesh};
double XMid[3] = {0.5*dxSubMesh, 0.5*dxSubMesh, 0.5*dxSubMesh};

// refinement level at given location
int get_level(double *x){
  int iLocation = 0;

  for(int i=0; i<3; i++) iLocation |= ((x[i] > 0.5*dxSubMesh) ? 1 : 0) << i;
  return (iCase >> iLocation) & 1;
}
//-----------------------------------------------------------------------------
//functions that returns the local resolution and time step
double localResolution(double *x) {
  if ((x[0] >= 0.25*dxSubMesh && x[0] < 0.75*dxSubMesh) || (x[1] >= 0.25*dxSubMesh && x[1] < 0.75*dxSubMesh) || (x[2] >= 0.25*dxSubMesh && x[2] < 0.75*dxSubMesh)) return 4*dxCell;

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
  if (FirstCall==true) {
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

  if (iLevel == 1) return;

  double xShift[3] = {2*(x[0]-XMid[0])/dxCell, 2*(x[1]-XMid[1])/dxCell,2*(x[2]-XMid[2])/dxCell};

  if( fabs( fabs(xShift[0]) - (int)fabs(xShift[0]) ) < 0.1) return;

  for(int i=0; i<3; i++) x[i] = 0.5*dxSubMesh+(2*floor(0.5*xShift[i]) + 1) * 0.5*dxCell;
}
//-----------------------------------------------------------------------------
int main(int argc,char **argv) {
  //execute the test only on ONE processor
  PIC::InitMPI();

  MPI_Comm New_Comm;
  int color=PIC::ThisThread,newID=0;

  MPI_Comm_split(MPI_COMM_WORLD,color,newID,&New_Comm);

  if (PIC::ThisThread==0) {
    MPI_GLOBAL_COMMUNICATOR=New_Comm;

    //init the interpolation routines
    PIC::InterpolationRoutines::Init();

    //init mesh
    //  amps_init();

    // output file name
    char fname[400];
    fstream foutCenterBased,foutCornerBased;
    bool open_file = true;

    PIC::InterpolationRoutines::CellCentered::cStencil *CenterBasedStencil;
    PIC::InterpolationRoutines::CornerBased::cStencil *CornerBasedStencil;

    if (PIC::ThisThread==0) for (iCase=0;iCase<nCase;iCase++) {
      amps_init();

      if (open_file==true) {
        sprintf(fname,"%s/test_InterpolateAMR-CenterBased-test.dat",PIC::OutputDataFileDirectory);
        foutCenterBased.open(fname, fstream::out);

        sprintf(fname,"%s/test_InterpolateAMR-CornerBased-test.dat",PIC::OutputDataFileDirectory);
        foutCornerBased.open(fname, fstream::out);

        open_file=false;
      }

      int i,nPoint=100000;
      double x0[3],x1CornerBased[3],x1CenterBased[3],*xcell,diff;

      for (int iPoint=0;iPoint<nPoint;iPoint++) {
        //      cout <<"\riCase="<<iCase<<" iPoint="<<iPoint<<"       ";
        // generate test point
        for (i=0;i<3;i++) {
          x0[i]=XMin[i]+0.5*dxCell+rnd()*(XMax[i]-XMin[i]-dxCell);
          x1CornerBased[i]=0.0,x1CenterBased[i]=0.0;
        }

        //get the interpolation stencil
        CenterBasedStencil=PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x0,NULL);
        CornerBasedStencil=PIC::InterpolationRoutines::CornerBased::InitStencil(x0,NULL);

        //check order of interpolation
        for (int icell=0; icell<CenterBasedStencil->Length; icell++) {
          xcell=CenterBasedStencil->cell[icell]->GetX();
          fix_coords(xcell);

          for (i=0; i<3; i++) x1CenterBased[i]+=xcell[i]*CenterBasedStencil->Weight[icell];
        }

        for (diff=0.0,i=0;i<3;i++) diff+=pow((x1CenterBased[i]-x0[i])/(0.5*(x1CenterBased[i]+x0[i])),2);

        diff=pow(diff,0.5);

        if (diff>0.001) {
          foutCenterBased << "iCase="   << iCase << " iPoint=" << iPoint  << " x=" << x0[0] << " " << x0[1] << " " << x0[2] << " diff=" << diff << "\nStencil:\n";

          for (int icell=0; icell<CenterBasedStencil->Length; icell++) {
            xcell = CenterBasedStencil->cell[icell]->GetX();
            fix_coords(xcell);

            foutCenterBased << "Cell " << icell << ": " << xcell[0] << " " << xcell[1] << " " << xcell[2] << " Weight=" << CenterBasedStencil->Weight[icell] << endl;
          }
        }

        //check the corner-based interpolation procedure
        for (int iCorner=0; iCorner<CornerBasedStencil->Length; iCorner++) {
          xcell=CornerBasedStencil->cell[iCorner]->GetX();
          fix_coords(xcell);

          for (i=0; i<3; i++) x1CornerBased[i]+=xcell[i]*CornerBasedStencil->Weight[iCorner];
        }

        for (diff=0.0,i=0;i<3;i++) diff+=pow((x1CornerBased[i]-x0[i])/(0.5*(x1CornerBased[i]+x0[i])),2);

        diff=pow(diff,0.5);

        if (diff>0.001) {
          foutCornerBased << "iCase="   << iCase << " iPoint=" << iPoint  << " x=" << x0[0] << " " << x0[1] << " " << x0[2] << " diff=" << diff << "\nStencil:\n";

          for (int iCorner=0; iCorner<CornerBasedStencil->Length; iCorner++) {
            xcell = CornerBasedStencil->cell[iCorner]->GetX();
            fix_coords(xcell);

            foutCornerBased << "Cell " << iCorner << ": " << xcell[0] << " " << xcell[1] << " " << xcell[2] << " Weight=" << CornerBasedStencil->Weight[iCorner] << endl;
          }
        }


      }
    }

    foutCenterBased.close();
    foutCornerBased.close();
  }

  //finish execution of the test
  MPI_Finalize();
  cout << "\nEnd of the run" << endl;


  return EXIT_SUCCESS;
}
