#ifndef IPIC3D_INTERFACE_HPP
#define IPIC3D_INTERFACE_HPP

#include <sstream>
#include <mpi.h>

int myProc;

extern "C" {  
  int ipic3d_init_mpi_(MPI_Fint *iComm, signed int* iProc, signed int* nProc);
  int ipic3d_init_(double *inittime); 
  int ipic3d_finalize_init_();
  int ipic3d_run_(double *time); 
  int ipic3d_save_restart_();
  int ipic3d_end_(); 
  int ipic3d_read_param_(char *param, int *nlines, int *ncharline, int *iProc);
  int ipic3d_from_gm_init_(int *paramint, double *paramreal, char *NameVar);
  int ipic3d_get_ngridpoints_(int *nPoint);
  int ipic3d_get_grid_(double *Pos_DI, int *n);
  int ipic3d_set_state_var_(double *Data_VI, int *iPoint_I);
  int ipic3d_set_dt_(double *DtSi);
  int ipic3d_set_swmf_dt_(double *DtSi);
  int ipic3d_cal_dt_(double *dt);
  int ipic3d_find_points_(int *nPoint, double *Xyz_DI, int *iProc_I);
  int ipic3d_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I, int *nVar);
}
#endif
