
#include <mpi.h>
#include "GMRES.h"
#include "Basic.h"
#include "parallel.h"
#include "errors.h"
#include "Alloc.h"
#include "TimeTasks.h"
#include "ipicdefs.h"

void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen,
  const double *b, int m, int max_iter, double tol, Field * field)
{
  if (m > xkrylovlen) {
    // m need not be the same for all processes,
    // we cannot restrict this test to the main process,
    // (although we could probably restrict it to the
    // process with the highest cartesian rank).
    eprintf("In GMRES the dimension of Krylov space(m) "
      "can't be > (length of krylov vector)/(# processors)\n");
  }
  bool GMRESVERBOSE = false;
  double initial_error, rho_tol;
  double *r = new double[xkrylovlen];
  double *im = new double[xkrylovlen];

  double *s = new double[m + 1];
  double *cs = new double[m + 1];
  double *sn = new double[m + 1];
  double *y = new double[m + 3];
  eqValue(0.0, s, m + 1);
  eqValue(0.0, cs, m + 1);
  eqValue(0.0, sn, m + 1);
  eqValue(0.0, y, m + 3);


  // allocate H for storing the results from decomposition
  double **H = newArr2(double, m + 1, m);
  for (int ii = 0; ii < m + 1; ii++)
    for (int jj = 0; jj < m; jj++)
      H[ii][jj] = 0;
  // allocate V
  double **V = newArr2(double, m+1, xkrylovlen);
  for (int ii = 0; ii < m+1; ii++)
    for (int jj = 0; jj < xkrylovlen; jj++)
      V[ii][jj] = 0;



  if (GMRESVERBOSE && is_output_thread()) {
    printf( "------------------------------------\n"
            "-             GMRES                -\n"
            "------------------------------------\n\n");
  }

  double normb = normP(b, xkrylovlen);
  if (normb == 0.0)
    normb = 1.0;

  int itr=0;
  for (itr = 0; itr < max_iter; itr++)
  {
    // r = b - A*x
    (field->*FunctionImage) (im, xkrylov);
    sub(r, b, im, xkrylovlen);
    initial_error = normP(r, xkrylovlen);

    if (itr == 0) {
      if (is_output_thread())
        printf("Initial residual: %g norm b vector (source) = %g\n",
          initial_error, normb);
        //cout << "Initial residual: " << initial_error << " norm b vector (source) = " << normb << endl;
      rho_tol = initial_error * tol;

      if ((initial_error / normb) <= tol) {
        if (is_output_thread())
          printf("GMRES converged without iterations: initial error < tolerance\n");
          //cout << "GMRES converged without iterations: initial error < tolerance" << endl;
        break;
      }
    }

    scale(V[0], r, (1.0 / initial_error), xkrylovlen);
    eqValue(0.0, s, m + 1);
    s[0] = initial_error;
    int k = 0;
    while (rho_tol < initial_error && k < m) {

      // w= A*V(:,k)
      double *w = V[k+1];
      (field->*FunctionImage) (w, V[k]);
      // old code (many MPI_Allreduce calls)
      //
      //const double av = normP(w, xkrylovlen);
      //for (register int j = 0; j <= k; j++) {
      //  H[j][k] = dotP(w, V[j], xkrylovlen);
      //  addscale(-H[j][k], w, V[j], xkrylovlen);
      //}

      // new code to make a single MPI_Allreduce call
      for (int j = 0; j <= k; j++)
      {
        y[j] = dot(w, V[j], xkrylovlen);
      }
      y[k+1] = norm2(w,xkrylovlen);
      {
        
        MPI_Allreduce(MPI_IN_PLACE, y, (k+2),
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      for (int j = 0; j <= k; j++) {
        H[j][k] = y[j];
        addscale(-H[j][k], V[k+1], V[j], xkrylovlen);
      }
      // Is there a numerically stable way to
      // eliminate this second all-reduce all?
      H[k+1][k] = normP(V[k+1], xkrylovlen);
      //
      // check that vectors are orthogonal
      //
      //for (register int j = 0; j <= k; j++) {
      //  dprint(dotP(w, V[j], xkrylovlen));
      //}

      double av = sqrt(y[k+1]);
      // why are we testing floating point numbers
      // for equality?  Is this supposed to say
      //if (av < delta * fabs(H[k + 1][k]))
      const double delta=0.001;
      if (av + delta * H[k + 1][k] == av)
      {
        for (int j = 0; j <= k; j++) {
          const double htmp = dotP(w, V[j], xkrylovlen);
          H[j][k] = H[j][k] + htmp;
          addscale(-htmp, w, V[j], xkrylovlen);
        }
        H[k + 1][k] = normP(w, xkrylovlen);
      }
      // normalize the new vector
      scale(w, (1.0 / H[k + 1][k]), xkrylovlen);

      if (0 < k) {

        for (int j = 0; j < k; j++)
          ApplyPlaneRotation(H[j + 1][k], H[j][k], cs[j], sn[j]);

        getColumn(y, H, k, m + 1);
      }

      const double mu = sqrt(H[k][k] * H[k][k] + H[k + 1][k] * H[k + 1][k]);
      cs[k] = H[k][k] / mu;
      sn[k] = -H[k + 1][k] / mu;
      H[k][k] = cs[k] * H[k][k] - sn[k] * H[k + 1][k];
      H[k + 1][k] = 0.0;

      ApplyPlaneRotation(s[k + 1], s[k], cs[k], sn[k]);
      initial_error = fabs(s[k]);
      k++;
    }

    k--;
    y[k] = s[k] / H[k][k];

    for (int i = k - 1; i >= 0; i--) {
      double tmp = 0.0;
      for (int l = i + 1; l <= k; l++)
        tmp += H[i][l] * y[l];
      y[i] = (s[i] - tmp) / H[i][i];

    }


    for (int j = 0; j < k; j++)
    {
      const double yj = y[j];
      double* Vj = V[j];
      for (int i = 0; i < xkrylovlen; i++)
        xkrylov[i] += yj * Vj[i];
    }

    if (initial_error <= rho_tol) {
      if (is_output_thread())
      {
        printf("GMRES converged at restart # %d; iteration #%d with error: %g\n",
          itr, k,  initial_error / rho_tol * tol);
        //cout << "GMRES converged at restart # " << itr << "; iteration #" << k << " with error: " << initial_error / rho_tol * tol << endl;
      }
      break;
    }
    if (is_output_thread() && GMRESVERBOSE)
    {
      printf("Restart: %d error: %g\n", itr,  initial_error / rho_tol * tol);
      //cout << "Restart: " << itr << " error: " << initial_error / rho_tol * tol << endl;
    }

  }
  if(itr==max_iter && is_output_thread())
  {
    printf("GMRES not converged !! Final error: %g\n",
      initial_error / rho_tol * tol);
    //cout << "GMRES not converged !! Final error: " << initial_error / rho_tol * tol << endl;
  }

  delete[]r;
  delete[]im;
  delete[]s;
  delete[]cs;
  delete[]sn;
  delete[]y;
  delArr2(H, m + 1);
  delArr2(V, m + 1);
  return;
}


void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn) {
  double temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}
