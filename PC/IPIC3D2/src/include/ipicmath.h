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

#ifndef _ipicmath_h_
#define _ipicmath_h_
#include "assert.h"
#include "math.h"
#include "stdlib.h" // for rand

// valid if roundup power is representable.
inline int
pow2roundup (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

// does not work if highest non-sign bit is set
inline int
pow2rounddown (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;

    // set all bits below highest bit
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    // set the bit higher than the highest bit
    x++;
    // shift it down and return it
    return (x >> 1);
}

// return integer ceiling of n over m
inline int ceiling_of_ratio(int n, int m)
{
  // probably this way is cheapest
  return ((n-1)/m+1);
  //return (n+m-1)/m;
  //const int out = ceil(n/double(m));
  //return out;
}

// round n up to next multiple of m
inline int roundup_to_multiple(int n, int m)
{
  //return ((n-1)/m+1)*m;
  return (n+m-1)/m*m;
}

// sample from clopen unit interval (0,1]
inline double sample_clopen_u_double()
{
  // old way (retained for sake of bit-wise code agreement)
  const double harvest = rand() / (double) RAND_MAX;
  return 1.0 - .999999 * harvest;
  // better way
  const double max_inv = 1./(double(RAND_MAX)+1);
  return (double(rand())+1)*max_inv;
}

// sample from open unit interval (0,1)
static inline double sample_open_u_double()
{
  const double max_inv = 1./(double(RAND_MAX)+2);
  return (double(rand())+1)*max_inv;
}

// sample from unit interval [0,1]
inline double sample_u_double()
{
  // old way
  return rand()/double(RAND_MAX);
  // faster way
  const double max_inv = 1./(double(RAND_MAX));
  return double(rand())*max_inv;
}

inline void sample_standard_maxwellian(double& u)
{
  // we sample a single component by pretending that it
  // is part of a two-dimensional joint distribution.
  const double prob = sqrt(-2.0 * log(sample_clopen_u_double()));
  const double theta = 2.0 * M_PI * sample_u_double();
  u = prob * cos(theta);
}

inline void sample_standard_maxwellian(double& u, double& v)
{
  // the distribution of the magnitude of (u,v)
  // can be integrated analytically
  const double prob = sqrt(-2.0 * log(sample_clopen_u_double()));
  const double theta = 2.0 * M_PI * sample_u_double();
  u = prob * cos(theta);
  v = prob * sin(theta);
}

inline void sample_standard_maxwellian(double& u, double& v, double& w)
{
  sample_standard_maxwellian(u,v);
  sample_standard_maxwellian(w);
}

inline void sample_maxwellian(double& u, double ut)
{
  sample_standard_maxwellian(u);
  u *= ut;
}

inline void sample_maxwellian(double& u, double& v, double& w,
  double ut, double vt, double wt)
{
  sample_standard_maxwellian(u,v,w);
  u *= ut;
  v *= vt;
  w *= wt;
}

inline void sample_maxwellian(
  double& u, double& v, double& w,
  double ut, double vt, double wt,
  double u0, double v0, double w0)
{
  sample_standard_maxwellian(u,v,w);
  u = u0 + ut*u;
  v = v0 + vt*v;
  w = w0 + wt*w;
}

// add or subtract multiples of L until x lies
// between 0 and L.
//
/** RIFAI QUESTA PARTE questo e' la tomba delle performance*/
//inline void MODULO(double *x, double L)
//{
//  *x = *x - floor(*x / L) * L;
//}
// version of previous method that assumes Linv = 1/L
// (faster if 1/L is precomputed)
inline double modulo(double x, double L, double Linv)
{
  return x - floor(x * Linv) * L;
}

#endif
