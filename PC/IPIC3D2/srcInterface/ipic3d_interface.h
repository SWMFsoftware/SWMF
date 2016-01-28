#ifndef IPIC3D_INTERFACE_HPP
#define IPIC3D_INTERFACE_HPP


#include <sstream>
#include <mpi.h>

int IPIC3D_init_mpi(MPI_Comm iComm, signed int* iProc, signed int* nProc); 
int IPIC3D_from_gm_init(int *paramint, double *paramreal, std::stringstream *ss);
int IPIC3D_init(double inittime); 
int IPIC3D_finalize_init();
int IPIC3D_run(double time); 
int IPIC3D_save_restart();
int IPIC3D_end(); 
double IPIC3D_getSItime();
int IPIC3D_read_paramin(std::stringstream *param);
int IPIC3D_get_ngridpoints(int *nPoint);
int IPIC3D_get_grid(double *Pos_DI, int n);
int IPIC3D_set_state_var(double *Data_VI, int *iPoint_I);
int IPIC3D_get_state_var(int nDim, int nPoint, double *Xyz_I, double *data_I, int nVar, std::string *NameVar);
int IPIC3D_find_point(int nPoint, double *Xyz_I, int *iProc_I);
int IPIC3D_set_dt(double DtSi);

#endif
