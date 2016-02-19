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

#ifndef __IPIC_DEFS_H__
#define __IPIC_DEFS_H__

typedef unsigned long long longid;
//typedef uint64_t longid; // requires #include <stdint.h>

// comment this out if OpenMP is not installed on your system.
#define USING_OMP

// uncomment the following line to use parallel hdf5
#define USING_PARALLEL_HDF5

// use precprocessor to remove former MPI_Barrier() calls.
//#define MPI_Barrier(args...)
#define former_MPI_Barrier(args...)

#define flds_MPI_Allreduce(first, args...) \
{ \
  assert(timeTasks.is_active(TimeTasks::REDUCE_FIELDS)); \
  timeTasks_set_task(TimeTasks::FLDS_COMM); \
  { \
    timeTasks_set_task(TimeTasks::FLDS_MPI_ALLREDUCE); \
    MPI_Allreduce(first, ## args); \
  } \
}

#define flds_MPI_Sendrecv_replace(first, args...) \
  { \
    assert(timeTasks.is_active(TimeTasks::FLDS_COMM)); \
    timeTasks_set_task(TimeTasks::FLDS_MPI_SENDRECV); \
    MPI_Sendrecv_replace(first, ## args); \
  }

// updates time for appropriate XXXX_MPI_SENDRECV task
#define ipic_MPI_Sendrecv_replace(first, args...) \
  { \
    double start_time = MPI_Wtime(); \
    MPI_Sendrecv_replace(first, ## args); \
    timeTasks.end_sendrecv(start_time); \
  }

#define pcls_MPI_Isend(first, args...) \
  { \
    assert(timeTasks.is_active(TimeTasks::PCLS_COMM)); \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Isend); \
    MPI_Isend(first, ## args); \
  }

#define pcls_MPI_Irecv(first, args...) \
  { \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Irecv); \
    MPI_Irecv(first, ## args); \
  }

#define pcls_MPI_Wait(first, args...) \
  { \
    assert(timeTasks.is_active(TimeTasks::PCLS_COMM)); \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Wait); \
    MPI_Wait(first, ## args); \
  }

#define pcls_MPI_Cancel(first, args...) \
  { \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Cancel); \
    MPI_Cancel(first, ## args); \
  }

#define pcls_MPI_Request_free(first, args...) \
  { \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Request_free); \
    MPI_Request_free(first, ## args); \
  }

#define pcls_MPI_Test(first, args...) \
  { \
    assert(timeTasks.is_active(TimeTasks::PCLS_COMM)); \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Test); \
    MPI_Test(first, ## args); \
  }

#define pcls_MPI_Waitany(first, args...) \
  { \
    assert(timeTasks.is_active(TimeTasks::PCLS_COMM)); \
    timeTasks_set_task(TimeTasks::PCLS_MPI_Waitany); \
    MPI_Waitany(first, ## args); \
  }

// determine the width of the vector unit
//
#if defined(__MIC__)
  const int VECBITS = 512;
#elif defined(__AVX__)
  const int VECBITS = 256;
#elif defined(__SSE_)
  const int VECBITS = 128;
#else
  const int VECBITS = 64;
#endif
const int VECBYTES = VECBITS/8;

// the number of doubles that fill a vector
const int DVECWIDTH = VECBYTES/sizeof(double);
const int SVECWIDTH = VECBYTES/sizeof(float);
//#define SINGLE_PRECISION_PCLS
//
// single precision does not seem to help on the MIC
typedef double pfloat;
//#ifdef SINGLE_PRECISION_PCLS
//  typedef float pfloat;
//  #ifdef __MIC__
//    #define VECTOR_WIDTH 16
//  #else
//    #define VECTOR_WIDTH 8
//  #endif
//#else
//  #ifdef __MIC__
//    #define VECTOR_WIDTH 8
//  #else
//    #define VECTOR_WIDTH 4
//  #endif
//  typedef double pfloat;
//#endif

#endif
