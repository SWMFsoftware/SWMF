#include "PC_wrapper_c.h"
#include "ipic3d_interface.h"
#include "ReadParam.h"

#include <iostream>

using std::cout;
using std::endl;
/* 
  Concerting data and structured used by C/Frotran and convert them
  into C++ constructs. With ecseption of MPI where we use the 
  C structure.
*/


int pc_wrapper_c_init_mpi_(MPI_Fint *iComm,signed int* iProc,signed int* nProc){
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);
  IPIC3D_init_mpi(c_iComm,iProc,nProc);
  myProc = *iProc;
  return(0);
}

int pc_wrapper_c_init_(double *inittime){
  IPIC3D_init(*inittime);
  return(0);
}

int pc_wrapper_c_finalize_init_(){
  IPIC3D_finalize_init();
  return(0);
}

int pc_wrapper_c_run_(double *time){
  IPIC3D_run(*time);
  *time = IPIC3D_getSItime();
  return(0);
}

int pc_wrapper_c_save_restart_(){
  IPIC3D_save_restart();
  return(0);
}

int pc_wrapper_c_end_(){
  IPIC3D_end();
  return(0);
}

int pc_wrapper_c_read_param_(char *param, int *nlines, int *ncharline, int *iProc){
  // convert character array to string stream object
  myProc = *iProc;
  std::stringstream  *ss;
  ss = char_to_stringstream(param, (*nlines)*(*ncharline), *ncharline, *iProc); 
  IPIC3D_read_paramin(ss);
  return(0);
}

int pc_wrapper_c_from_gm_init_(int *paramint, double *paramreal, char *NameVar){
  std::stringstream  *ss;
  ss = new std::stringstream;
  //ss->write(NameVar, 1);
  (*ss)<<NameVar;
  IPIC3D_from_gm_init(paramint, paramreal, ss);
  return(0);
}

int pc_wrapper_c_get_ngridpoints_(int *nPoint){
  IPIC3D_get_ngridpoints(nPoint);
  return(0);
}

int pc_wrapper_c_get_grid_(double *Pos_DI, int *n){
  IPIC3D_get_grid(Pos_DI, *n);
  return(0);  
}

int pc_wrapper_c_set_state_var_(char *NameVar, int *nNameVar, int *nVar, double *Data_VI, int *iPoint_I){
  std::stringstream  *ss;
  ss = new std::stringstream;
  ss->write(NameVar, *nNameVar); 
  IPIC3D_set_state_var(ss, *nVar, Data_VI, iPoint_I);
  return(0);
}

int pc_wrapper_c_set_dt_(double *DtSi){
  IPIC3D_set_dt(*DtSi);
  return(0);
}

int pc_wrapper_c_find_points_(int *nDim, int *nPoint, double *Xyz_DI, int *iProc_I){
  IPIC3D_find_point(*nPoint, Xyz_DI, iProc_I);
  return(0);
}

int pc_wrapper_c_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,int *nVar, char *NameVarChar, int *nChar){
  std::string NameVar;

  IPIC3D_get_state_var(*nDim, *nPoint, Xyz_I, data_I, *nVar, &NameVar);
  return(0);
}



