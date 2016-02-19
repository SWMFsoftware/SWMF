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

//#ifdef _OPENMP
#include <omp.h>
//#endif
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h> // for uint64_t

// icpc testv.cpp -mmic -fopenmp -vec-report6 -o test.mic
// icpc testv.cpp -fopenmp -vec-report6 -o test

#define ALIGNMENT (64)
#define ALIGNED(X) __assume_aligned(X, ALIGNMENT)
#define AlignedAlloc(T, NUM) \
(T *const __restrict__)(_mm_malloc(sizeof(T)*NUM, ALIGNMENT))
#define AlignedFree(S) (_mm_free(S))

class SpeciesParticle
{
  double x[3];
  double q;
  double u[3];
  long long ID;
 public:
  // accessors
  long long get_ID()const{ return ID; }
  double get_x(int i)const{ return x[i]; }
  double get_u(int i)const{ return u[i]; }
  double get_q()const{ return q; }
  void set_ID(long long in){ ID=in; }
  void set_x(int i, double in) { x[i] = in; }
  void set_u(int i, double in) { u[i] = in; }
  void set_q(double in) { q = in; }
  // alternative accessors
  double get_x()const{ return x[0]; }
  double get_y()const{ return x[1]; }
  double get_z()const{ return x[2]; }
  double get_u()const{ return u[0]; }
  double get_v()const{ return u[1]; }
  double get_w()const{ return u[2]; }
  void set_x(double in){ x[0]=in; }
  void set_y(double in){ x[1]=in; }
  void set_z(double in){ x[2]=in; }
  void set_u(double in){ u[0]=in; }
  void set_v(double in){ u[1]=in; }
  void set_w(double in){ u[2]=in; }
  void set(long long _ID,
    double _x, double _y, double _z,
    double _u, double _v, double _w,
    double _q)
  {
    ID = _ID;
    x[0] = _x; x[1] = _y; x[2] = _z;
    u[0] = _u; u[1] = _v; u[2] = _w;
    q = _q;
  }
};

class BEfield
{
  double BE[6];
 public:
  void set(
    double Bx,
    double By,
    double Bz,
    double Ex,
    double Ey,
    double Ez)
  {
    BE[0] = Bx;
    BE[1] = By;
    BE[2] = Bz;
    BE[3] = Ex;
    BE[4] = Ey;
    BE[5] = Ez;
  }
  double getB(int j)
  {
    return BE[j];
  }
  double getE(int j)
  {
    return BE[j+3];
  }
};

#define NUMPCLS 1000000

// fast random number generator
// http://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
//inline int fastrand() { 
//  long long g_seed = (214013*g_seed+2531011); 
//  return (g_seed>>16)&0x7FFF; 
//} 
// see http://www.digicortex.net/node/22
// for fast generator using SSE4
// http://en.wikipedia.org/wiki/SSE4

double usample()
{
  const double RAND_MAX_inv = 1./RAND_MAX;
  return rand()*RAND_MAX_inv;
}

void test_mover()
{
  // initialize particles and field data
  //
  SpeciesParticle* pcls = new SpeciesParticle[NUMPCLS];
  BEfield* flds = new BEfield[NUMPCLS];
  ALIGNED(pcls);
  //
  // initialize with random data
  //
  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].set( i,
      usample(), usample(), usample(),
      usample(), usample(), usample(),
      usample());
    flds[i].set(
      usample(), usample(), usample(),
      usample(), usample(), usample());
  }

  // global constants
  const double dt = 1+usample();
  const double dto2 = .5 * dt;
  double qdto2mc = 1+usample();

  const double start = omp_get_wtime();

  const int NUM_PCLS_MOVED_AT_A_TIME = 2;
  const int D=4;
  assert(NUMPCLS%NUM_PCLS_MOVED_AT_A_TIME==0);
  for (int pidx = 0; pidx < NUMPCLS; pidx+=NUM_PCLS_MOVED_AT_A_TIME)
  {
    // copy the particles
    SpeciesParticle* pcl[NUM_PCLS_MOVED_AT_A_TIME];
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    {
      pcl[i] = &pcls[pidx+i];
    }
    // actually, all the particles are aligned,
    // but the compiler should be able to see that.
    ALIGNED(pcl[0]);
    double xorig[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
    double uorig[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
    double  xavg[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
    double  uavg[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
    // gather data into vectors
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    for(int j=0;j<3;j++)
    {
      xavg[i][j] = xorig[i][j] = pcl[i]->get_x(j);
      uorig[i][j] = pcl[i]->get_u(j);
    }
    const int NiterMover = 1;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {

      // // compute weights for field components
      // //
      // double weights[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(ALIGNMENT)));
      // int cx[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      // for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      // {
      //   grid->get_safe_cell_and_weights(xavg[i],cx[i],weights[i]);
      // }

      // arr1_double_get field_components[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(ALIGNMENT)));
      // for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      // {
      //   get_field_components_for_cell(field_components[i],fieldForPcls,
      //     cx[i][0],cx[i][1],cx[i][2]);
      // }

      double E[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      double B[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      //// could do this with memset
      //for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      //for(int j=0;j<3;j++)
      //{
      //  E[i][j]=0;
      //  B[i][j]=0;
      //}
      //for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      //for(int j=0;j<3;j++)
      //for(int c=0; c<8; c++)
      //{
      //  B[i][j] += weights[i][c] * field_components[i][c][j];
      //  E[i][j] += weights[i][c] * field_components[i][c][j+3];
      //}
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        B[i][j] = flds[pidx+i].getB(j);
        E[i][j] = flds[pidx+i].getE(j);
      }
      double Om[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        Om[i][j] = qdto2mc*B[i][j];
      }

      // can these dot products vectorize if
      // NUM_PCLS_MOVED_AT_A_TIME is large enough?
      double omsq[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(ALIGNMENT)));
      double denom[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        omsq[i] = Om[i][0] * Om[i][0]
                + Om[i][1] * Om[i][1]
                + Om[i][2] * Om[i][2];
        denom[i] = 1.0 / (1.0 + omsq[i]);
      }
      // solve the position equation
      double ut[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        ut[i][j] = uorig[i][j] + qdto2mc * E[i][j];
      }
      double udotOm[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        udotOm[i] = ut[i][0] * Om[i][0]
                  + ut[i][1] * Om[i][1]
                  + ut[i][2] * Om[i][2];
      }
      // solve the velocity equation 
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        // With D=4, the cross product can be handled with 64-bit cross-product swizzles.
        // With D=4, the scalar-vector product can be handled with 64-bit broadcast
        // with 4-element granularity.  See pages 27, 29, and 30 in the
        // "Intel Xeon Phi Coprocessor Instruction Set Reference Manual".
        uavg[i][0] = (ut[i][0] + (ut[i][1] * Om[i][2] - ut[i][2] * Om[i][1] + udotOm[i] * Om[i][0])) * denom[i];
        uavg[i][1] = (ut[i][1] + (ut[i][2] * Om[i][0] - ut[i][0] * Om[i][2] + udotOm[i] * Om[i][1])) * denom[i];
        uavg[i][2] = (ut[i][2] + (ut[i][0] * Om[i][1] - ut[i][1] * Om[i][0] + udotOm[i] * Om[i][2])) * denom[i];
      }
      // update average position
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        xavg[i][j] = xorig[i][j] + uavg[i][j] * dto2;
      }
    } // end of iteration
    // update the final position and velocity (scatter)
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    for(int j=0;j<3;j++)
    {
      pcl[i]->set_x(j, xorig[i][j] + uavg[i][j] * dt);
      pcl[i]->set_u(j, 2.*uavg[i][j] - uorig[i][j]);
    }
  }
  const double end = omp_get_wtime();
  printf(" moving in test_mover_AoS() took %g s.\n", end - start);
}

void test_mover_AoS_vec()
{
  // initialize particles and field data
  //
  SpeciesParticle* pcls = new SpeciesParticle[NUMPCLS];
  BEfield* flds = new BEfield[NUMPCLS];
  ALIGNED(pcls);
  //
  // initialize with random data
  //
  for(int i=0;i<NUMPCLS;i++)
  {
    pcls[i].set( i,
      usample(), usample(), usample(),
      usample(), usample(), usample(),
      usample());
    flds[i].set(
      usample(), usample(), usample(),
      usample(), usample(), usample());
  }

  // global constants
  const double dt = 1+usample();
  const double dto2 = .5 * dt;
  double qdto2mc = 1+usample();

  const double start = omp_get_wtime();

  const int NUM_PCLS_MOVED_AT_A_TIME = 2;
  const int D=3; // dimensions of space
  assert(NUMPCLS%NUM_PCLS_MOVED_AT_A_TIME==0);
  for (int pidx = 0; pidx < NUMPCLS; pidx+=NUM_PCLS_MOVED_AT_A_TIME)
  {
    // copy the particles
    SpeciesParticle* pcl[NUM_PCLS_MOVED_AT_A_TIME];
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    {
      pcl[i] = &pcls[pidx+i];
    }
    // actually, all the particles are aligned,
    // but the compiler should be able to see that.
    ALIGNED(pcl[0]);
    double xorig[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
    double uorig[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
    double  xavg[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
    double  uavg[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
    // gather data into vectors
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    {
      double* xorigi = xorig+D*i;
      double* uorigi = uorig+D*i;
      SpeciesParticle* pcli = pcl[i];
      ALIGNED(pcli);
      // #pragma simd
      for(int j=0;j<D;j++)
      {
        xorigi[j] = pcli->get_x(j);
        uorigi[j] = pcli->get_u(j);
      }
    }
    for(int k=0;k<NUM_PCLS_MOVED_AT_A_TIME*D;k++)
    {
      xavg[k] = xorig[k];
    }
    const int NiterMover = 1;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {

      // // compute weights for field components
      // //
      // double weights[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(ALIGNMENT)));
      // int cx[NUM_PCLS_MOVED_AT_A_TIME][D] __attribute__((aligned(ALIGNMENT)));
      // for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      // {
      //   grid->get_safe_cell_and_weights(xavg[i],cx[i],weights[i]);
      // }

      // arr1_double_get field_components[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(ALIGNMENT)));
      // for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      // {
      //   get_field_components_for_cell(field_components[i],fieldForPcls,
      //     cx[i][0],cx[i][1],cx[i][2]);
      // }

      double E[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      double B[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      //// could do this with memset
      //for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      //for(int j=0;j<D;j++)
      //{
      //  E[i][j]=0;
      //  B[i][j]=0;
      //}
      //for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      //for(int c=0; c<8; c++)
      //for(int j=0;j<D;j++)
      //{
      //  B[i][j] += weights[i][c] * field_components[i][c][j];
      //  E[i][j] += weights[i][c] * field_components[i][c][j+D];
      //}
      //#pragma simd collapse(2)
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<D;j++)
      {
        B[i*D+j] = flds[pidx+i].getB(j);
        E[i*D+j] = flds[pidx+i].getE(j);
      }
      double Om[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      for(int k=0; k<NUM_PCLS_MOVED_AT_A_TIME*D; k++)
      {
        Om[k] = qdto2mc*B[k];
      }

      double denom[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        double omsqp1 = 1.0 + Om[i*D+0] * Om[i*D+0]
                + Om[i*D+1] * Om[i*D+1]
                + Om[i*D+2] * Om[i*D+2];
        denom[i] = 1.0 / omsqp1;
      }
      // solve the position equation
      double ut[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      for(int k=0; k<NUM_PCLS_MOVED_AT_A_TIME*D;k++)
      {
        ut[k] = uorig[k] + qdto2mc * E[k];
      }
      double udotOm[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        udotOm[i] = ut[i*D+0] * Om[i*D+0]
                  + ut[i*D+1] * Om[i*D+1]
                  + ut[i*D+2] * Om[i*D+2];
      }
      // solve the velocity equation 
      //
      // cross-products are difficult to vectorize,
      // so we handle the bad stuff separately.
      double uxOm_A[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      double uxOm_B[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        uxOm_A[i*D+0] = ut[i*D+1] * Om[i*D+2];
        uxOm_A[i*D+1] = ut[i*D+2] * Om[i*D+0];
        uxOm_A[i*D+2] = ut[i*D+0] * Om[i*D+1];
      }
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        uxOm_B[i*D+0] =  ut[i*D+2] * Om[i*D+1];
        uxOm_B[i*D+1] =  ut[i*D+0] * Om[i*D+2];
        uxOm_B[i*D+2] =  ut[i*D+1] * Om[i*D+0];
      }
      // with D=4 this can be done with a broadcast with 4-element granularity
      double denomv[NUM_PCLS_MOVED_AT_A_TIME*D] __attribute__((aligned(ALIGNMENT)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        denomv[i*D+0] = denom[i];
        denomv[i*D+1] = denom[i];
        denomv[i*D+2] = denom[i];
      }
      for(int k=0; k<NUM_PCLS_MOVED_AT_A_TIME*D;k++)
      {
        int i=k/D; //k=i*D+j;
        uavg[k] = ut[k] + (uxOm_A[k] - uxOm_B[k] + udotOm[i]*Om[k]) * denom[i];
      }
      // update average position
      for(int k=0; k<NUM_PCLS_MOVED_AT_A_TIME*D;k++)
      {
        xavg[k] = xorig[k] + uavg[k] * dto2;
      }
    } // end of iteration
    // update the final position and velocity (scatter)
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    {
      double* xorigi = xorig+D*i;
      double* uorigi = uorig+D*i;
      double* uavgi = uavg+D*i;
      SpeciesParticle* pcli = pcl[i];
      ALIGNED(pcli);
      // compiler thinks inefficient to vectorize this
      for(int j=0;j<D;j++)
      {
        pcli->set_x(j, xorigi[j] + uavgi[j] * dt);
        pcli->set_u(j, 2.*uavgi[j] - uorigi[j]);
      }
    }
  }
  const double end = omp_get_wtime();
  printf(" moving in test_mover_AoS_vec() took %g s.\n", end - start);
}

//void time_test_mover()
//{
//  const double start = omp_get_wtime();
//  test_mover_AoS();
//  const double end = omp_get_wtime();
//  printf("test_mover_AoS() took %g s.\n", end - start);
//}

void vmul ( int n, double *a, double *b, double *c )
{
  // #pragma omp parallel for simd
  //#pragma omp simd
  for ( int i = 0; i < n ; i ++ )
     a[i] += b[i] * c[i];
}

#define N 16*100000

double a[N] __attribute__((aligned(ALIGNMENT)));
double b[N] __attribute__((aligned(ALIGNMENT)));
double c[N] __attribute__((aligned(ALIGNMENT)));

void test_vectorization()
{
  for ( int i = 0; i < N ; i++ ) {
    a[i] = 0;
    b[i] = i;
    c[i] = i*2;
  }


  double start,end;

  start = omp_get_wtime();
  vmul( N, a, b, c );
  end = omp_get_wtime();

  printf(" vmul took %g s.\n", end - start);

}

int main()
{
  #ifdef _OPENMP
    printf("_OPENMP is set.\n");
  #else
    printf("_OPENMP is not set.\n");
  #endif
  //printf("omp_get_thread_num==%d\n",omp_get_thread_num());
  printf("omp_get_max_threads==%d\n",omp_get_max_threads());
  //test_vectorization();
  test_mover_AoS();
  test_mover_AoS_vec();
}

