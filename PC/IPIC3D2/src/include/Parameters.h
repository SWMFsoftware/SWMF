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

#ifndef _Parameters_h_
#define _Parameters_h_

// namespace provides a more flexible, succinct singleton via "using Parameters"
//
namespace Parameters
{
  enum Enum
  {
    SoA=0, // struct of arrays
    AoS, // array of structs
    // for moments type
    AoSvec,
    SoAvec,
    AoSintr,
    // for mover type
    SoA_vec_onesort,
    AoS_vec_onesort,
    SoA_vec_resort,
    AoS_vec_resort,
    AoS_Relativistic
  };

  void init_parameters();

  bool get_USING_AOS();
  bool get_SORTING_SOA();
  bool get_SORTING_PARTICLES();
  // for resorting particles with each iteration of mover
  bool get_RESORTING_PARTICLES();
  inline bool get_USING_XAVG() { return get_RESORTING_PARTICLES(); }
  bool get_VECTORIZE_MOMENTS();
  //bool get_VECTORIZE_MOVER();
  Enum get_MOVER_TYPE();
  Enum get_MOMENTS_TYPE();

  int get_multiple_of_vector_width_in_doubles();
  // blocksize and numblocks for use in BlockCommunicator
  int get_blockSize();
  int get_numBlocks();
  bool get_doWriteOutput();
}
#endif
