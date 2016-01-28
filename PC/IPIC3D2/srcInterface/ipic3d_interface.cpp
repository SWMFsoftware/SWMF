#include <iomanip>
#include "ipic3d_interface.h"
#include "multi_ipic3d_domain.h"
#include "MPIdata.h"
#include <iostream>
#include <cstdlib>
#include "timing_SWMF.h"

using namespace iPic3D;
using namespace std;

int IPIC3D_init_mpi(MPI_Comm iComm, signed int* iProc,signed int* nProc){

  // At this time we do not have a proper timstep 
  // so we only do steping internaly iSimCycle
  
  //std::cout<<"ipic3dComm pre :: "<<ipic3dComm<<endl;
  if(MPI_SUCCESS != MPI_Comm_dup(iComm, &ipic3dComm)){
    cout<<"IPIC3D_init_mpi :: Can not copy  MPI comunicator, arborting!"<<endl;
    cout.flush();
    abort();
  } 
  
  //std::cout<<"ipic3dComm init :: "<<ipic3dComm<<endl;
  myrank     = *iProc;
  numProc    = *nProc; 
  param      = NULL;
  nIPIC      = 0;
  iIPIC      = 0;
  SimRun  = new iPic3D::c_Solver *[nmaxIPIC];
  iSimCycle = new int[nmaxIPIC];
  
  for(int i = 0; i < nmaxIPIC; i++){
    firstcall[i] = true;
    starttime[i] = 0.0;
    SimRun[i]    = NULL;
    //std::cout<<" SimRun["<<i<<"]  = "<<SimRun[i]<<std::endl;
  }
  
  return(0);
}

int IPIC3D_read_paramin(std::stringstream *paramin){    
  param = paramin;
  iIPIC++;

  return(0);
}


int IPIC3D_init(double inittime){

  //  Called by  PC_init_session
  //  Not used as this time as 
  //  most of the grid are set up when
  //  we read PARAM.in, and the rest 
  //  when IPIC3D recive data for the 
  //  first time. 

  //std::cout<<"IPIC3D_init"<<std::endl;

  timenow = inittime;

  return(0);

}

int IPIC3D_from_gm_init(int *paramint, double *paramreal, std::stringstream *ss){
 
  char **dummy = NULL; // dummy argument
  int firstIPIC = nIPIC;
  string nameFunc="PC: IPIC3D_from_gm_init";
  timing_start(nameFunc);

  // number of dimensions in GM
  nDim = paramint[0];

  // number of PIC regions (added?)
  nIPIC +=paramint[1];
  
  // the number of PIC fluids ns = nFluid + nSpecies - 1
  int ns = paramint[3] + paramint[4] - 1;

  MPIdata::init(ipic3dComm, &myrank, &numProc);
  for(int i = firstIPIC; i < nIPIC; i++){ 
    param->clear();
    param->seekg(0,param->beg);
    starttime[i] = timenow;
    SimRun[i]  = new iPic3D::c_Solver;
    SimRun[i]->Init(0, dummy, timenow, param, i, paramint,
		    &paramreal[(i - firstIPIC)*9], 
		    &paramreal[(nIPIC - firstIPIC)*9], ss,true);
  }
  timing_stop(nameFunc);
  return(0);
}

int IPIC3D_finalize_init(){

  // This function is called by the coupler the first time 
  // it want to couple from GM -> PC. At this point we should have all
  // information to finnish the initialization of SimRun[i]

  char **dummy = NULL; // dummy argument

  string nameFunc="PC: IPIC3D_finalize_init";
  timing_start(nameFunc);

  // now we should have all the the infomation
  for(int i = 0; i < nIPIC; i++){     
    SimRun[i]->Init(0, dummy, timenow);
    SimRun[i]->CalculateMoments();
    SimRun[i]->WriteOutput(0);
    iSimCycle[i] = SimRun[i]->FirstCycle();
  }
  timing_stop(nameFunc);
  return(0);
}


int IPIC3D_run(double time){
  bool b_err = false;
 
  timenow = time;
  string nameFunc="PC: IPIC3D_run";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++){ 
 
   if (SimRun[i]->get_myrank() == 0)
       cout << " ======= Cycle " << iSimCycle[i] << ", dt=" << SimRun[i]->getDt() << " ======= " << endl;

   timing_start("PC: SyncWithFluid");  
   SimRun[i]->SyncWithFluid(iSimCycle[i]);
   timing_stop("PC: SyncWithFluid");

   timing_start("PC: CalculateField");  
   SimRun[i]->CalculateField();
   timing_stop("PC: CalculateField");

   timing_start("PC: ParticlesMover");  
   b_err = SimRun[i]->ParticlesMover();
   timing_stop("PC: ParticlesMover");

   timing_start("PC: CalculateBField");  
   if (!b_err) SimRun[i]->CalculateB();
   timing_stop("PC: CalculateBField");
   
   timing_start("PC: GatherMoments");  
   if (!b_err) SimRun[i]->CalculateMoments();
   timing_stop("PC: GatherMoments");  
   
    if (b_err) {
      iSimCycle[i] = SimRun[i]->LastCycle() + 1;
    }

    timing_start("PC: WriteOutput");  
    SimRun[i]->WriteOutput(iSimCycle[i]+1);
    timing_stop("PC: WriteOutput");
    
    SimRun[i]->WriteConserved(iSimCycle[i]);

    iSimCycle[i]++;    
  }
  timing_stop(nameFunc);
  return(0);
}

int IPIC3D_save_restart(){
  string nameFunc="PC: IPIC3D_save_restart";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++) 
    SimRun[i]->WriteRestart(iSimCycle[i]);
  
  timing_stop(nameFunc);  
  return(0);
}

int IPIC3D_end(){
  try {
    for(int i = 0; i < nIPIC; i++){
      SimRun[i]->Finalize();
    }
    delete[] param;
    delete[] SimRun;
    delete[] iSimCycle;
    delete param;    
  }
  catch( int e) {
    cout<<" IPIC3D have not a clean dealocation, this is a error chatch : "<< e <<endl;
    cout.flush();
  }
  return 0;
}

double IPIC3D_getSItime(){
  
  // All simulations are in the same time zone
  return(SimRun[0]->getSItime());
}

int IPIC3D_get_ngridpoints(int *nPoint){
  
  *nPoint = 0;
  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->GetNgridPnt(&nGridPntSim[i]);
    nShiftGridPntSim[i] = *nPoint;
    *nPoint += nGridPntSim[i];
  }
  // Fortran operates on an 2D array [ndim,nPoint]
  *nPoint = *nPoint/nDim;
  return(0);
}

int IPIC3D_get_grid(double *Pos_DI, int n){

  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->GetGridPnt(&Pos_DI[nShiftGridPntSim[i]]);
  
  return(0);
}

int IPIC3D_set_state_var(double *Data_VI, int *iPoint_I){
  string nameFunc="PC: IPIC3D_set_state_var";
  timing_start(nameFunc);

  // WARNING  iPoint_I is a reindexing of state var array.
  
  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->setStateVar(Data_VI,
                           &iPoint_I[(nShiftGridPntSim[i]/nDim)]);
  }
  timing_stop(nameFunc);
  return(0);
}

int IPIC3D_get_state_var(int nDim, int nPoint, double *Xyz_I, double *data_I, int nVar, string *NameVar){
  string nameFunc="PC: IPIC3D_get_state_var";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->getStateVar(nDim, nPoint, Xyz_I, data_I, nVar, NameVar);    
  }
  timing_stop(nameFunc);  
  return(0);
}

int IPIC3D_find_point(int nPoint, double *Xyz_I, int *iProc_I){
  string nameFunc="PC: IPIC3D_find_point";
  timing_start(nameFunc);
  
  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->findProcForPoint(nPoint, Xyz_I, iProc_I);
  }

  timing_stop(nameFunc);
  return(0);
}



int IPIC3D_set_dt( double DtSi){
  
  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->setSIDt(DtSi);
  return(0);
}



