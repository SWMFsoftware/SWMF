/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/***************************************************************************
  Bessel.h  -  
  -------------------
begin                : Tue Jan 30 2007

 ***************************************************************************/

#ifndef Bessel_H
#define Bessel_H

#include <math.h>
#include <errors.h>
// #include "gsl/gsl_sf_bessel.h"


/*! \brief Calculate Bessel functions of the first kind J_n \param lambda_s Argument of function \param nmax Sequence J_n(x) for n=0 ... nmax+1 is calculated \param bessel_Jn_array Array of values for J_n corresponding to n=0 ... nmax+1 \param bessel_Jn_prime_array Array of values for J_n_prime ((ie derivate of J_n) for n=0 ... nmax \note J_n array goes to n=nmax + 1, but the J_n_prime array goes to nmax only. \note This function uses the corresponding call in the GSL. */
void calc_bessel_Jn_seq(double lambda, int nmax, double bessel_Jn_array[], double bessel_Jn_prime_array[]) {
  int gsl_fn_result = 0.0;
  // CHANGE IT HERE: PUT HERE THE BESSEL FUNCTION
  // gsl_fn_result = gsl_sf_bessel_Jn_array( 0, nmax+1, lambda, bessel_Jn_array );

  if (gsl_fn_result)
    eprintf("Error in calc_bessel_Jn_seq !!!!");

  bessel_Jn_prime_array[0] = -bessel_Jn_array[1];

  for (int i = 1; i <= nmax; ++i)
    bessel_Jn_prime_array[i] = 0.5 * (bessel_Jn_array[i - 1] - bessel_Jn_array[i + 1]);

}


#endif
