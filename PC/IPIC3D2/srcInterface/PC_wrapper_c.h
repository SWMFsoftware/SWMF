#ifndef PC_WRAPPER_C_H
#define PC_WRAPPER_C_H

#include <mpi.h>

// practical for debuging
int myProc;

extern "C" {
  
  int pc_wrapper_c_init_mpi_(MPI_Fint *iComm, signed int* iProc, signed int* nProc);
  int pc_wrapper_c_init_(double *inittime);
  int pc_wrapper_c_finalize_init_();
  int pc_wrapper_c_run_(double *timetime);
  int pc_wrapper_c_save_restart_();
  int pc_wrapper_c_end_();
  int pc_wrapper_c_read_param_(char *param, int *nlines, int *ncharline, int *iProc);
  int pc_wrapper_c_from_gm_init_(int *paramint, double *paramreal, char *NameVar);
  int pc_wrapper_c_get_ngridpoints_(int *nPoint);
  int pc_wrapper_c_get_grid_(double *Pos_DI, int *n);
  int pc_wrapper_c_set_state_var_(char *NameVar, int *nNameVar, int *nVar, double *Data_VI, int *iPoint_I);
  int pc_wrapper_c_set_dt_(double *DtSi);
  int pc_wrapper_c_find_points_(int *nDim, int *nPoint, double *Xyz_DI, int *iProc_I);
  int pc_wrapper_c_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I, int *nVar, char *NameVarChar, int *nChar);
}

#endif
