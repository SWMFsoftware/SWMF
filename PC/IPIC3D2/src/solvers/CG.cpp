
#include <mpi.h>
#include "CG.h"
#include "Basic.h"
#include "parallel.h"

/**
 * 
 * Electrostatic Field with 3 components(x,y,z) defined on a 2-D grid
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
 * @version 2.0
 *
 */

bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, FIELD_IMAGE FunctionImage, Field * field) {
  // allocate residual, image, p, b, calculated on central points
  double *r = new double[xkrylovlen];
  double *v = new double[xkrylovlen];
  double *z = new double[xkrylovlen];
  double *im = new double[xkrylovlen];
  double c, t, d, initial_error;
  eqValue(0.0, r, xkrylovlen);
  eqValue(0.0, v, xkrylovlen);
  eqValue(0.0, z, xkrylovlen);
  eqValue(0.0, im, xkrylovlen);

  int i = 0;
  bool CONVERGED = false;
  bool CGVERBOSE = false;
  // initial guess for x: all the components are equal to 0
  eqValue(0, xkrylov, xkrylovlen);
  // Compute r = b -Ax
  (field->*FunctionImage) (im, xkrylov);
  sub(r, b, im, xkrylovlen);
  // v = r
  eq(v, r, xkrylovlen);
  c = dotP(r, r, xkrylovlen);
  initial_error = sqrt(c);
  if (is_output_thread())
    printf("CG Initial error: %g\n", initial_error);
    //cout << "CG Initial error: " << initial_error << endl;
  if (initial_error < 1E-16)
    return (true);
  while (i < maxit) {
    (field->*FunctionImage) (z, v);
    t = c / dotP(v, z, xkrylovlen);
    // x(i+1) = x + t*v
    addscale(t, xkrylov, v, xkrylovlen);
    // r(i+1) = r - t*z
    addscale(-t, r, z, xkrylovlen);
    d = dotP(r, r, xkrylovlen);
    if (CGVERBOSE && is_output_thread())
    {
      printf("Iteration # %d - norm of residual relative to initial error %g\n",
        i, sqrt(d) / initial_error);
      //cout << "Iteration # " << i << " - norm of residual relative to initial error " << sqrt(d) / initial_error << endl;
    }
    if (sqrt(d) < tol * initial_error) {
      // if (d < tol){

      if (is_output_thread())
      {
        printf("CG converged at iteration # %d with error %g\n", i, sqrt(d));
        //cout << "CG converged at iteration # " << i << " with error " << sqrt(d) << endl;
      }
      CONVERGED = true;
      break;
    }
    else if (sqrt(d) > 10E8 * initial_error) {
      if (is_output_thread()) {
        printf("CG not converged after %d iterations\n", maxit);
        //cerr << "CG not converging" << endl;
        //cerr << "CG stopped" << endl;

      }
      CONVERGED = false;
      return (CONVERGED);
      break;

    }

    // calculate the new v
    addscale(1, d / c, v, r, xkrylovlen);
    c = d;
    i++;

  }
  if (i == maxit)
  {
    printf("CG not converged after %d iterations\n", maxit);
    //cout << "CG not converged after " << maxit << " iterations" << endl;
  }
  // deallocate
  delete[]r;
  delete[]im;
  delete[]v;
  delete[]z;
  return (CONVERGED);
}

