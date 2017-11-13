
extern "C"{
  void linear_solver_matvec_c(double * VecIn, double * VecOut, int n );
  void  linear_solver_matvec(double * VecIn, double * VecOut, int *n);
  void  linear_solver_gmres( double * Rhs_I,double * x_I, int * lInit, int * n, int * nKrylov,
			     double * Tol, int * nIter, int * iError, int * lTest);
}
