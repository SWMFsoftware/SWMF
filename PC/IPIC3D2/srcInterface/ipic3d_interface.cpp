#include <iomanip>
#include <iostream>
#include <cstdlib>
#include "ipic3d_interface.h"
#include "multi_ipic3d_domain.h"
#include "MPIdata.h"
#include "timing_SWMF.h"
#include "ReadParam.h"

using namespace iPic3D;
using namespace std;

int ipic3d_init_mpi_(MPI_Fint *iComm,signed int* iProc,signed int* nProc){
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);

  // At this time we do not have a proper timstep 
  // so we only do steping internaly iSimCycle
  
  if(MPI_SUCCESS != MPI_Comm_dup(c_iComm, &ipic3dComm)){
    cout<<"IPIC3D_init_mpi :: Can not copy  MPI comunicator, arborting!"<<endl;
    cout.flush();
    abort();
  } 
  
  myrank     = *iProc;
  numProc    = *nProc; 
  param      = NULL;
  nIPIC      = 0;
  iIPIC      = 0;
  SimRun  = new iPic3D::c_Solver *[nmaxIPIC];
  iSimCycle = new int[nmaxIPIC];

  isFirstSession = true;
  isInitilized   = false;
  for(int i = 0; i < nmaxIPIC; i++){
    firstcall[i] = true;
    starttime[i] = 0.0;
    SimRun[i]    = NULL;
  }
  
  myProc = *iProc;
  
  MPIdata::init(ipic3dComm, &myrank, &numProc);
    
  return(0);
}

int ipic3d_init_(double *inittime){
  //  Called by  PC_init_session
  //  Not used as this time as 
  //  most of the grid are set up when
  //  we read PARAM.in, and the rest 
  //  when IPIC3D recive data for the 
  //  first time. 
  timenow = *inittime;
  return(0);
}

int ipic3d_read_param_(char *paramIn, int *nlines, int *ncharline, int *iProc){
  // So far, iPIC3D does not support multi sessions.
  if(!isFirstSession) return(0);

  // convert character array to string stream object
  myProc = *iProc;
  std::stringstream  *ss;
  ss = char_to_stringstream(paramIn, (*nlines)*(*ncharline), *ncharline, *iProc); 
  param = ss;
  iIPIC++;
  isFirstSession = false;
  return(0);
}

int ipic3d_from_gm_init_(int *paramint, double *paramreal, char *NameVar){
  if(isInitilized) return(0);
  std::stringstream  *ss;
  ss = new std::stringstream;
  (*ss)<<NameVar;

  char **dummy = NULL; // dummy argument
  int firstIPIC;
  string nameFunc="PC: ipic3d_from_gm_init";
  timing_start(nameFunc);

  // number of dimensions in GM
  nDim = paramint[0];

  firstIPIC = 0; 
  nIPIC =paramint[1];
  
  // the number of PIC fluids ns = nFluid + nSpecies - 1
  int ns = paramint[3] + paramint[4] - 1;

  for(int i = firstIPIC; i < nIPIC; i++){ 
    param->clear();
    param->seekg(0,param->beg);
    starttime[i] = timenow;
    SimRun[i]  = new iPic3D::c_Solver;
    SimRun[i]->Init(0, dummy, timenow, param, i, paramint,
		    &paramreal[(i - firstIPIC)*9], 
		    &paramreal[(nIPIC - firstIPIC)*9], ss,true);
    SimRun[i]->SetCycle(0);
  }
  timing_stop(nameFunc);
  isInitilized = true;
  return(0);
}

int ipic3d_finalize_init_(){

  // This function is called by the coupler the first time 
  // it want to couple from GM -> PC. At this point we should have all
  // information to finnish the initialization of SimRun[i]

  char **dummy = NULL; // dummy argument

  string nameFunc="PC: ipic3d_finalize_init_";
  timing_start(nameFunc);

  // now we should have all the the infomation
  for(int i = 0; i < nIPIC; i++){     
    SimRun[i]->Init(0, dummy, timenow);
    SimRun[i]->CalculateMoments();   
    iSimCycle[i] = SimRun[i]->FirstCycle();
    SimRun[i]->WriteOutput(iSimCycle[i]);
  }
  timing_stop(nameFunc);
  return(0);
}

int ipic3d_run_(double *time){
  bool b_err = false;
 
  timenow = *time;
  string nameFunc="PC: IPIC3D_run";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++){ 
 
   if (SimRun[i]->get_myrank() == 0)
       cout << " ======= Cycle " << iSimCycle[i] << ", dt=" << SimRun[i]->getDt() << " ======= " << endl;

   SimRun[i]->SetCycle(iSimCycle[i]);
   
   timing_start("PC: SyncWithFluid");  
   SimRun[i]->SyncWithFluid(iSimCycle[i]);
   timing_stop("PC: SyncWithFluid");

   timing_start("PC: CalculateField");  
   SimRun[i]->CalculateField(iSimCycle[i]);
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

    SimRun[i]->updateSItime();
    
    timing_start("PC: WriteOutput");  
    SimRun[i]->WriteOutput(iSimCycle[i]+1);
    timing_stop("PC: WriteOutput");
    
    SimRun[i]->WriteConserved(iSimCycle[i]+1);
 
    iSimCycle[i]++;    
  }
  
  // All simulations are in the same time zone
  *time = SimRun[0]->getSItime();
  timing_stop(nameFunc);
  return(0);
}

int ipic3d_save_restart_(){
  string nameFunc="PC: ipic3d_save_restart";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++) 
    SimRun[i]->WriteRestart(iSimCycle[i]-1);
  
  timing_stop(nameFunc);  
  return(0);
}

int ipic3d_get_ngridpoints_(int *nPoint){
  
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

int ipic3d_get_grid_(double *Pos_DI, int *n){
  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->GetGridPnt(&Pos_DI[nShiftGridPntSim[i]]);
  
  return(0);
}

int ipic3d_set_state_var_(double *Data_VI, int *iPoint_I){
  string nameFunc="PC: ipic3d_set_state_var_";
  timing_start(nameFunc);

  // WARNING  iPoint_I is a reindexing of state var array.
  
  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->setStateVar(Data_VI,
                           &iPoint_I[(nShiftGridPntSim[i]/nDim)]);
  }
  timing_stop(nameFunc);
  return(0);
}

int ipic3d_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I, int *nVar){
  string nameFunc="PC: IPIC3D_get_state_var";
  timing_start(nameFunc);

  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->getStateVar(*nDim, *nPoint, Xyz_I, data_I, *nVar);    
  }
  timing_stop(nameFunc);  
  return(0);
}

int ipic3d_find_points_(int *nPoint, double *Xyz_I, int *iProc_I){
  string nameFunc="PC: ipic3d_find_points_";
  timing_start(nameFunc);
  
  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->findProcForPoint(*nPoint, Xyz_I, iProc_I);
  }
  timing_stop(nameFunc);
  return(0);
}

int ipic3d_set_swmf_dt_( double *DtSi){  
  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->setSIDt(*DtSi,true);
  return(0);
}

int ipic3d_set_dt_( double *DtSi){  
  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->setSIDt(*DtSi,false);
  return(0);
}

int ipic3d_cal_dt_(double *dtOut){
  // Each PIC domain should be treated seperately in the future. Now, all the
  // domains use the same time step.

  double dt = 1e10, dt0;
  
  for(int i = 0; i < nIPIC; i++){
    dt0 = SimRun[i]->calSIDt();
    if(dt0 < dt) dt = dt0;
  }
  *dtOut = dt;
  return(0);
}

int ipic3d_end_(){
  try {
    for(int i = 0; i < nIPIC; i++){
      //SimRun[i]->Finalize();
      delete SimRun[i];
    }
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
