#include "linear_solver_wrapper_c.h"

//matvec function interface for user
void (*linear_solver_matvec_c)(double* VecIn, double * VecOut, int n);

void linear_wrapper_matvec_c(double* VecIn, double * VecOut, int n) {
  linear_solver_matvec_c(VecIn,VecOut,n);
}

