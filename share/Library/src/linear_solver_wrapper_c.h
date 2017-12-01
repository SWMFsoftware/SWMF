// function interface for the fortran code ModLinearSolver
extern "C"{
  void linear_wrapper_matvec_c(double * VecIn, double * VecOut, int n );
  void linear_solver_wrapper(const char * solverType,double * tolerance, int * nIteration, 
			     int* nVar,int* nDim,int* nI,int* nJ,int* nK,int* nBlock,
			     int* iComm, double* Rhs_I, double* x_I, 
			     double* PrecondParam,double *precond_matrix, int* lTest);
}

//matvec function interface for user
extern void (*linear_solver_matvec_c)(double* VecIn, double * VecOut, int n);

