#include <iomanip>
#include "ipic3d_interface.h"
#include "multi_ipic3d_domain.h"
#include "MPIdata.h"
#include <iostream>
#include <cstdlib>

using namespace iPic3D;
using namespace std;

int IPIC3D_init_mpi(MPI_Comm iComm, signed int* iProc,signed int* nProc){

  // At this time we do not have a proper timstep 
  // so we only do steping internaly iSimCycle
  
  //std::cout<<"ipic3dComm pre :: "<<ipic3dComm<<endl;
  if(MPI_SUCCESS != MPI_Comm_dup(iComm, &ipic3dComm)){
    cout<<"IPIC3D_init_mpi :: Can not copy  MPI comunicator, arborting!"<<endl;
    cout.flush();
    //abort();
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
    // ([ndim, nRegions, nVar, nFluid, nSpecies...], [min X, max X, dx, ....], 
    // [Charge*ns, Mass*ns, PePerPtotal, NormX, NormU, NormM]
  }
  return(0);
}

int IPIC3D_finalize_init(){

  // This function is called by the coupler the first time 
  // it want to couple from GM -> PC. At this point we should have all
  // information to finnish the initialization of SimRun[i]

  char **dummy = NULL; // dummy argument
  
  // now we should have all the the infomation
  for(int i = 0; i < nIPIC; i++){     
    SimRun[i]->Init(0, dummy, timenow);
    SimRun[i]->CalculateMoments();
    SimRun[i]->WriteOutput(0);
    iSimCycle[i] = SimRun[i]->FirstCycle();
  }
  return(0);
}


int IPIC3D_run(double time){
  bool b_err = false;
  //static bool firstcall = true;
  //char **dummy = NULL; // dummy argument
  
  //std::cout<<"IPIC3D_run :: "<<time<<std::endl;
 
  timenow = time;
  /*
  // initial setup may come from the coupler  
  if(firstcall){
    // now we should have all the the infomation
    for(int i = 0; i < nIPIC; i++){ 
      if(SimRun[i] == NULL){
        SimRun[i]  = new iPic3D::c_Solver;
	SimRun[i]->InitMpi(ipic3dComm, &myrank, &numProc);
      }
      SimRun[i]->Init(0, dummy, starttime[i], param); 
      SimRun[i]->GatherMoments();
      SimRun[i]->WriteOutput(0);
      iSimCycle[i] = SimRun[i]->FirstCycle();
    }
    firstcall = false;
  }
  */
  for(int i = 0; i < nIPIC; i++){ 
 
   if (SimRun[i]->get_myrank() == 0)
       cout << " ======= Cycle " << iSimCycle[i] << ", dt=" << SimRun[i]->getDt() << " ======= " << endl;
    
   SimRun[i]->SyncWithFluid(iSimCycle[i]);
 
   SimRun[i]->CalculateField();
   
   b_err = SimRun[i]->ParticlesMover();

   if (!b_err) SimRun[i]->CalculateB();
   if (!b_err) SimRun[i]->CalculateMoments();

    if (b_err) {
      iSimCycle[i] = SimRun[i]->LastCycle() + 1;
    }

    SimRun[i]->WriteOutput(iSimCycle[i]+1);
    SimRun[i]->WriteConserved(iSimCycle[i]);

    iSimCycle[i]++;    
  }

  return(0);
}

int IPIC3D_save_restart(){

  for(int i = 0; i < nIPIC; i++) 
    SimRun[i]->WriteRestart(iSimCycle[i]);
  
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
    //cout<<"nGridPntSim["<<i<<"] = "<<nGridPntSim[i]<<endl;
    //cout<<"nShiftGridPntSim["<<i<<"] = "<<nShiftGridPntSim[i]<<endl;
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

int IPIC3D_set_state_var(std::stringstream *ss, int nVar, double *Data_VI, int *iPoint_I){

  // WARNING  iPoint_I is a reindexing of state var array.
  
  //cout<<"IPIC3D_set_state_var "<<endl;

  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->setStateVar(ss, nVar, Data_VI,
                           &iPoint_I[(nShiftGridPntSim[i]/nDim)]);
  }
  return(0);
}

int IPIC3D_get_state_var(int nDim, int nPoint, double *Xyz_I, double *data_I, int nVar, string *NameVar){
  //cout<<"IPIC3D_get_state_var"<<endl;

  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->getStateVar(nDim, nPoint, Xyz_I, data_I, nVar, NameVar);    
  }

  return(0);
}

int IPIC3D_find_point(int nPoint, double *Xyz_I, int *iProc_I){
  
  //cout<<"IPIC3D_find_point "<<endl;

  for(int i = 0; i < nIPIC; i++){
    SimRun[i]->findProcForPoint(nPoint, Xyz_I, iProc_I);
  }

  return(0);
}



int IPIC3D_set_dt( double DtSi){
  
  for(int i = 0; i < nIPIC; i++)
    SimRun[i]->setSIDt(DtSi);

  cout<<"setdt = "<<DtSi<<endl;
  return(0);
}



