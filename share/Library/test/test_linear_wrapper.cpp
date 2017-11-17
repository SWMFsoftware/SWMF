# include <cstdlib>
# include <iostream>
# include <cmath>
# include <mpi.h>
# include "test_linear_wrapper.h"
# include <stdio.h>

int rank, size;

void pass_cell_info(double * Vec, int n, double  *head, double *end){
  double  buffer_forward, buffer_backward;
  MPI_Status  status;
  
  if (rank != size-1){
    MPI_Recv(&buffer_forward,1,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
    (*end)=buffer_forward;
  }else{
    (*end) = 0.0;
  }

  if (rank != 0) {
    buffer_forward = Vec[0];
    MPI_Send(&buffer_forward,1,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
  }
  
  if (rank != 0){
    MPI_Recv(&buffer_backward,1,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
    (*head)=buffer_backward;
    // printf("buffer_backward:%f,head:%f\n",buffer_backward,*head);
  }else{
    (*head) = 0.0;
  }
  
  if (rank != size-1){
    buffer_backward = Vec[n-1];
    MPI_Send(&buffer_backward,1,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
  }
}

void linear_solver_matvec_c(double* VecIn, double * VecOut, int n){
  
  double  head;
  double  end;

  pass_cell_info(VecIn,n,&head,&end);//pass ghost cell info
  
  double alpha=0.1;

  // printf("rank:%d,head:%f,end:%f\n",rank,head,end);
  for (int i=1; i<n-1; i++){ 
    VecOut[i] = -alpha*(VecIn[i-1]+VecIn[i+1])+(1+2.*alpha)*VecIn[i];
  }

  if(rank != 0){
    VecOut[0] = -alpha*((head)+VecIn[1])+(1+2.*alpha)*VecIn[0];
  }else{
    VecOut[0] = 0.0; //set bc
  }

  if(rank != size-1){
    VecOut[n-1] = -alpha*(VecIn[n-2]+(end))+(1+2.*alpha)*VecIn[n-1];
  }else{
    VecOut[n-1] = 0.0; //set bc
  }

}

//========================================================================

int main(){

  // a test of the backward implicit solver for 1D diffusion equation
  // first type of boundary condition is used, i.e., T(at two ends) =0
  MPI_Init(NULL,NULL);
 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n=101; // number of true cells
  int nLen =n/size; //length of first size-1 procs
  int nCellProc=nLen;//number of cell for this proc
  if (rank==size-1){
    nCellProc =(n%nLen==0)?nLen:n-(size-1)*nLen;//nCellProc not includes ghost cells  
  }
 
  //alpha = 0.1; // the coefficient (dt)/(dx)^2, defined in MatVec_C.cpp
  //dt =1e-5;
  double  dx=1.0/((double)n-1);
  double nTotalStep=1000;

  int iMatvec=1;
  int lInit=1;
  int nKrylov=200;
  double Tol=1e-7;
  int nIter=200;
  int iError;
  int lTest=0;

  // set the init condition as a sine function
  // the Sol_I contains the first initial guess the same as Rhs_I
  double Rhs_I[nCellProc];
  double Sol_I[nCellProc];
  double Pi = 3.14159265;
  
  for (int i=0; i<nCellProc; i++) {
    Sol_I[i] = 100.*sin(Pi*(i+rank*nLen)/(n-1));
  }

  // Convert C MPI communicator to Fortran integer
  MPI_Fint iComm = MPI_Comm_c2f(MPI_COMM_WORLD);
  MPI_Status  status;
  int nTotalLenProc;// cells including ghost cells

  for (int iStep=1; iStep<=nTotalStep; iStep++){   
    if (rank==0)      Sol_I[0]           = 0.0;
    if (rank==size-1) Sol_I[nCellProc-1] = 0.0;
    
    for (int i=0; i<nCellProc; i++) Rhs_I[i] = Sol_I[i];
    /* test input
    for (int iProc=0;iProc<size;iProc++){
      if (rank==iProc) {
	printf("rank:%d,nCellProc:%d\n", rank,nCellProc);
	for (int i=0; i<nCellProc; i++) {
	  printf("rank:%d input i:%d, %f\n",rank,i,Rhs_I[i]);
	}
      }
    }
    */ 
    linear_solver_gmres(&iMatvec, Rhs_I, Sol_I, &lInit, &nCellProc, &nKrylov, &Tol, &nIter, &iError, &lTest, &iComm);
  }
  
  // the analytic solution is 100*exp(-Pi^2*t)sin(Pi*x)
  // dt=1e-5; the peak at t=1e-2 is 90.6
  
  int tag=0, nLength;
  double buffer[nCellProc], output[n];

  // send solution to proc 0 
  if (rank != 0){
    for (int i=0; i<nCellProc; i++) buffer[i]=Sol_I[i];
    MPI_Send(buffer, nCellProc, MPI_DOUBLE, 0,  tag,
             MPI_COMM_WORLD);
  }

  // recieve solution and write it out
  if (rank==0){

    // copy local solution into output array
    for (int i=0; i<nCellProc; i++) {
      output[i]=Sol_I[i];
    }

    // recieve from others
    for (int i=1; i<size; i++) {
      // number of elements to recieve
      nLength = nLen;

      // Last processor can have different number of elements
      if (i == size-1) nLength = (n%nLen==0)?nLen:n-nLen*(size-1);

      MPI_Recv(buffer, nLength, MPI_DOUBLE, i,  tag, MPI_COMM_WORLD,  &status);

      // print result out
      for (int j=0; j<nLength; j++) output[j+i*nLen] = buffer[j];
    }
    for (int i=0; i<n; i++) printf("%f\n",output[i]);
  }

  MPI_Finalize();
  return 0;
}
