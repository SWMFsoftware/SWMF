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

#include "IDgenerator.h"
#include "MPIdata.h"
#include "math.h"
#include "debug.h"

// Increasing the number of unique IDs to 2^64 would be possible
// using uint64_t or via large double precision IDs created by
// grouping processors e.g. into 2048 different groups each
// having a unique exponent.
//
// Idea of how to assign unique double precision IDs
// using most of the double precision bits:
//
//   MPI_Comm_rank returns int, so assume that
//   there are at most 2^31 processors.
//   Putting 2^30 processors in 2^10 groups
//   means that there are 2^20 processors in a group.
//   let group_idx range from 1 to 2^10.
//   let ID_increment = group_idx.
//   set counter as if there are only 2^20 processors.
//   then multiply counter by group_idx.
//
// But I need to be careful about the sign bits
// in the mantissa and exponent.
//
// We can work out the details when we really need more
// than 2^53 unique particle IDs.
//
// Points relevant to use of float versus integer keys in code:
// * In particle sorting and communication, performance
//   imperatives dictate against special handling of the ID
//   field.  Casting allows use of a double field to store a
//   64-bit integer.  Correct communication would require the
//   expense of calls to MPI_Pack, though I think that this would
//   only be an issue in the unlikely event of communicating
//   between big-endian and little-endian architectures.
// * double supports consecutive integers up to 2^53.
//   (uint64_t supports consecutive integers up to 2^64.)
// * double can still represent up to 2^64 distinct IDs
//   if one gives up on the IDs being consecutive integers,
//   and 2^62 of these values are positive integers.
// * There is wider support for double precision floats
//   than for 64-bit integers (e.g. among legacy codes,
//   architectures, and databases, which apparently allow use of
//   floating point numbers as keys in key-value pairs).  For
//   example, on my laptop (which has a 32-bit architecture),
//   int64_t and long long are simply mapped to int, but double
//   precision is supported -- though it's unlikely that we would
//   want more than 4 billion particles on a 32-bit architecture.
//
// Considerations relevant to use of float versus integer keys
// in saved data (a separable issue):
// * what kind of support is provided by scientific databases
//   for integer versus double-precision keys?
// * use of double means that there is no need for a separate
//   HDF5 method to save the ID field.

void doubleIDgenerator::reserve_particles_in_range(double lowest, double highest)
{
  int num_threads_in_this_proc = omp_get_max_threads();
  delete [] counter;
  counter = new double[num_threads_in_this_proc];

  double num_procs = MPIdata::get_nprocs();
  // technically should add 1 to the difference.
  // in integer case could overflow.
  // after division or floor it doesn't matter anyway.
  double max_pcls = highest - lowest;
  // allocate range of IDs unique to this processor
  double max_pcls_per_proc = floor(max_pcls / num_procs);
  double counter_offset_for_this_proc
    = lowest + max_pcls_per_proc * MPIdata::get_rank();
  // allocate range of IDs unique to each thread
  double max_pcls_per_thread = floor(max_pcls_per_proc / num_threads_in_this_proc);
  for(int i=0;i<num_threads_in_this_proc;i++)
  {
    counter[i] = counter_offset_for_this_proc + i*max_pcls_per_thread;
  }
}

// assumes that all processors will generate at most nop particles
// and that only the master thread will be used to generate particles
// (used on startup)
//
void doubleIDgenerator::reserve_num_particles(int nop)
{
  num_threads_in_this_proc = 1;
  counter = new double[num_threads_in_this_proc];
  counter[0] = nop * MPIdata::get_rank();
  // dprintf("initialized first particle ID to %g", counter[0]);
}

