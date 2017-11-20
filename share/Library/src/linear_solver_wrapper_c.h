// function interface for the fortran code ModLinearSolver
extern "C"{
  void linear_wrapper_matvec_c(double * VecIn, double * VecOut, int n );
  void linear_solver_gmres(double * Rhs_I, double * x_I, 
			   int * lInit, int * n, int * nKrylov,
			   double * Tol, int * nIter, int * iError, int * lTest, int * iComm);  
}

//matvec function interface for user
extern void (*linear_solver_matvec_c)(double* VecIn, double * VecOut, int n);

