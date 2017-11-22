# include <cstdlib>
# include <iostream>
# include <cmath>
# include <mpi.h>
# include "linear_solver_wrapper_c.h"
# include <stdio.h>

int rank, nProc;

void pass_cell_info(double * Vec, int n, double  *head, double *end){
  double  buffer_forward, buffer_backward;
  MPI_Status  status;
  
  if (rank != nProc-1){
    MPI_Recv(&buffer_forward,1,MPI_DOUBLE,rank+1,1,MPI_COMM_WORLD,&status);
    (*end)=buffer_forward;
  }else{
    (*end) = n-1;
  }

  if (rank != 0) {
    buffer_forward = Vec[0];
    MPI_Send(&buffer_forward,1,MPI_DOUBLE,rank-1,1,MPI_COMM_WORLD);
  }
  
  if (rank != 0){
    MPI_Recv(&buffer_backward,1,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,&status);
    (*head)=buffer_backward;
  }else{
    (*head) = 0.0;
  }
  
  if (rank != nProc-1){
    buffer_backward = Vec[n-1];
    MPI_Send(&buffer_backward,1,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
  }
}


void linear_solver_matvec_c_default(double* VecIn, double * VecOut, int n){
  
  double  head;
  double  end;

  pass_cell_info(VecIn,n,&head,&end);//pass ghost cell info

  for (int i=1; i<n-1; i++){ 
    VecOut[i] = VecIn[i-1]+VecIn[i+1]-2.*VecIn[i];
  }

  if(rank != 0){
    VecOut[0] = head+VecIn[1]-2*VecIn[0];
  }else{
    VecOut[0] = VecIn[1]-2*VecIn[0]; //set bc
  }

  if(rank != nProc-1){
    VecOut[n-1] = VecIn[n-2]+end-2.*VecIn[n-1];
  }else{
    VecOut[n-1] = VecIn[n-2]-2.*VecIn[n-1]; //set bc
  }

}

double ** GetPrecondMatrix(int nVectorLen, int nDim){
  double ** matrix = new double * [2*nDim+1];
  matrix[0] = new double [(2*nDim+1)*nVectorLen];
  for (int i=1;i<2*nDim+1;i++) matrix[i]=matrix[i-1]+nVectorLen;

  for (int i=0;i<nVectorLen;i++){
    // set the diagonal vector   
    matrix[0][i] = -2;
    // set the sub diagonal vector 
    matrix[1][i] = (i!=0)?1:0;   
    // set the super diagonal vector 
    matrix[2][i] = (i!=nVectorLen-1)?1:0;
    
  }

  return matrix;

}

//========================================================================

int main(){
  // set C matvec function pointer to the actual function 
  // the function pointer is declared and used in 
  // share/Library/src/linear_solver_wrapper_c.cpp
  linear_solver_matvec_c = linear_solver_matvec_c_default;

  // a test of the backward implicit solver for 1D diffusion equation
  // first type of boundary condition is used, i.e., T(at two ends) =0
  MPI_Init(NULL,NULL);
 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);

  int n=200; // number of true cells excluding boundary cells
  int nLen =n/nProc; //length of first nProc-1 procs
  int nCellProc=nLen;//number of cell for this proc
  if (rank==nProc-1){
    nCellProc =(n%nLen==0)?nLen:n-(nProc-1)*nLen;//nCellProc not includes ghost cells  
  }
 
  //alpha = 0.1; // the coefficient (dt)/(dx)^2, defined in MatVec_C.cpp
  //dt =1e-5;
  double  dx=1.0/((double)n-1);
  double nTotalStep=2;
  double Tol=1e-5;
  int nIter=200;
  int iError;
  int lTest = (rank == 0);

  // set the init condition as a sine function
  // the Sol_I contains the first initial guess the same as Rhs_I
  double Rhs_I[nCellProc];
  double Sol_I[nCellProc];
  
  for (int i=0; i<nCellProc; i++) {
    Rhs_I[i] = 0;
  }

  //set boundary condition 
  if (rank==nProc-1) Rhs_I[nCellProc-1] = -n-1; 

  // Convert C MPI communicator to Fortran integer
  MPI_Fint iComm = MPI_Comm_c2f(MPI_COMM_WORLD);
  MPI_Status  status;
  int nTotalLenProc;// cells including ghost cells
    
  double ** precond_matrix_II;
  // parameter to choose preconditioner types
  //0:No precondition; 1: BILU; 2:DILU;
    //[-1,0): MBILU; 
  double PrecondParam=2;
  int nVar=1, nDim=1, nI=nCellProc, nJ=1, nK=1, nBlock=1;

  // get precond matrix for 1D case 
  precond_matrix_II = GetPrecondMatrix(nVar*nVar*nI*nJ*nK,nDim); 
  linear_solver_wrapper("GMRES", &Tol,&nIter, &nVar, &nDim,
			&nI, &nJ, &nK, &nBlock, &iComm, Rhs_I,
			Sol_I, &PrecondParam, precond_matrix_II[0], 
			&lTest);
      
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
    for (int i=1; i<nProc; i++) {
      // number of elements to recieve
      nLength = nLen;

      // Last processor can have different number of elements
      if (i == nProc-1) nLength = (n%nLen==0)?nLen:n-nLen*(nProc-1);

      MPI_Recv(buffer, nLength, MPI_DOUBLE, i,  tag, MPI_COMM_WORLD,  &status);

      // print result out
      for (int j=0; j<nLength; j++) output[j+i*nLen] = buffer[j];
    }
    for (int i=0; i<n; i++) {
      printf("%f\n",output[i]);
      if (fabs(output[i]-(i+1))>1e-3) printf("ERROR!\n");
    }
  }

  MPI_Finalize();
  return 0;
}
