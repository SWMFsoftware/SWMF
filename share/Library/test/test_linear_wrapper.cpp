# include <cstdlib>
# include <iostream>
# include <cmath>
# include <mpi.h>
# include "linear_solver_matvec_c.h"
# include <stdio.h>

int main(){

  // a test of the backward implicit solver for 1D diffusion equation
  // first type of boundary condition is used, i.e., T(at two ends) =0
  int n=101; // number of true cells
 
  //alpha = 0.1; // the coefficient (dt)/(dx)^2, defined in MatVec_C.cpp
  //dt =1e-5;
  double  dx=1.0/((double)n-1);
  double nTotalStep=1000;
  int lInit=1;
  int nKrylov=100;
  double Tol=1e-5;
  int nIter=200;
  int iError;
  int lTest=0;

  //set the init condition as a sin function
  // the Sol_I contains the first initial guess the same as Rhs_I
  double Rhs_I[n];
  double Sol_I[n];
  double Pi=3.1415926;
  for (int i=1; i<n-1; i++) Sol_I[i] = 100.*sin(Pi*i/(n-1.));
  //set boundary condition
  Sol_I[0]=0.0;
  Sol_I[n-1]=0.0;


  for (int iStep=1; iStep<=nTotalStep; iStep++){   
 
    for (int i=0; i<n; i++) Rhs_I[i] = Sol_I[i];
    
    linear_solver_gmres(Rhs_I, Sol_I, &lInit, &n , &nKrylov, &Tol,
			&nIter, &iError, &lTest);
      
  }
  //the analytic solution is 100*exp(-Pi^2*t)sin(Pi*x)
  // dt=1e-5; the peak at t=1e-2 is 90.6

  /*
  FILE *fout= fopen("test_linear_wrapper.tmp","w");
  for (int i=0; i<n; i++) fprintf(fout, "%f\n", Sol_I[i]);
  fclose(fout);
  */
  for (int i=0; i<n; i++) printf("%f\n", Sol_I[i]);    

  return 0;

}
