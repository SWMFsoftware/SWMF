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

#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

// class to accumulate node-centered species moments
// 
class Moments10
{
  private:
    arr4_double arr;
    int nx;
    int ny;
    int nz;
  public:
    void set_to_zero();

    // fetch accessors (write access)
    arr4_double fetch_arr() { return arr; }

    Moments10(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      arr (nxn, nyn, nzn,10)
    {
    };
    ~Moments10(){};
};

#endif
