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

#ifndef IDgenerator_h
#define IDgenerator_h

#include "ompdefs.h"

// Class to generate unique double-precision IDs
//
// I use double rather than 64-bit integers for reasons
// discussed in .cpp file.
//
static const double DINTMAX = 0x20000000000000p0; // 2^53
class doubleIDgenerator
{
  double * counter;
  int num_threads_in_this_proc;
  // largest consecutive positive integer representable by double
 public:
  doubleIDgenerator():counter(0),num_threads_in_this_proc(0){};
  void reserve_num_particles(int nop);
  void reserve_particles_in_range(double lowest, double highest=DINTMAX);
  ~doubleIDgenerator() { delete [] counter; }
 public: 
  double generateID()
  {
    return counter[omp_get_thread_num()]++;
  }
};

#endif // IDgenerator_h
