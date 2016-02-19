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

#ifndef ompdefs_H
#define ompdefs_H

#include <stdio.h>
#include "asserts.h"
// the compiler sets _OPENMP if the -openmp flag is used
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0;}
inline int omp_get_max_threads(){ return 1;}
#define omp_set_num_threads(num_threads)
#endif

class Caller_to_SetMaxThreadsForScope{
 int max_threads;
 public:
  Caller_to_SetMaxThreadsForScope(int i)
  {
    max_threads = omp_get_max_threads();
    // omp_set_num_threads should have been
    // called omp_set_max_threads
    omp_set_num_threads(i);
  }
  ~Caller_to_SetMaxThreadsForScope()
  {
    // restore the original maximum number of threads
    omp_set_num_threads(max_threads);
  }
};

#define set_max_threads_for_scope(num_threads) \
  Caller_to_SetMaxThreadsForScope \
  instanceOfCaller_to_SetMaxThreadsForScope(num_threads);

#endif
