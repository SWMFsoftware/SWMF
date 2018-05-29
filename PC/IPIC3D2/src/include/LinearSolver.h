#ifndef LINEARSOLVER
#define LINEARSOLVER

#include "linear_solver_wrapper_c.h"
#include "../../srcInterface/multi_ipic3d_domain.h"

typedef void (*MATVEC) (double *, double *, int);
void iPIC3D_MaxwellImage(double *vecIn, double *vecOut, int n);
void iPIC3D_PoissonImage(double *vecIn, double *vecOut, int n);
void iPIC3D_matvec_particle_correction(double *vecIn, double *vecOut, int n);


#endif
