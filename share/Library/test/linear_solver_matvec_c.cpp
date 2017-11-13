# include <cstdlib>
# include <iostream>
# include <cmath>
# include <mpi.h>
# include "linear_solver_matvec_c.h"


void linear_solver_matvec_c(double* VecIn, double * VecOut, int n){
  
  double alpha=0.1;
  for (int i=1; i<n-1; i++){    
    VecOut[i] = -alpha*(VecIn[i-1]+VecIn[i+1])+(1+2*alpha)*VecIn[i];
  }
  VecOut[0]=0.0;
  VecOut[n-1]=0.0;

}


