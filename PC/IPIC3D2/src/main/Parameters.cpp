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

#include "Parameters.h"

using namespace Parameters;

//********** edit these parameters *********
// parameters relevant to parallelization
//
// This should be a multiple of the width of the vector unit
// as measured in double precision floating point numbers
int Parameters::get_multiple_of_vector_width_in_doubles(){ return 8; }
//
// parameters relevant to parallel IO
//bool Parameters::call_H5Block3dSetChunk() { return true; }
//
// parameters relevant to particle communication
//
// int Parameters::get_blockSize() { return 64; }
int Parameters::get_blockSize() { return 8192; }
int Parameters::get_numBlocks() { return 4; }
//
// parameters relevant to parallelization
//
bool Parameters::get_VECTORIZE_MOMENTS() { return false; }
// supported options: SoA AoS
Parameters::Enum Parameters::get_MOMENTS_TYPE() { return AoS; }
// supported options: SoA AoS AoSvec AoSintr AoS_vec_onesort SoA_vec_resort
Parameters::Enum Parameters::get_MOVER_TYPE() { return AoS; }
//********** derived parameters *********

static bool SORTING_PARTICLES;
static bool RESORTING_PARTICLES;
static bool USING_AOS;
static bool SORTING_SOA;

void Parameters::init_parameters()
{
  RESORTING_PARTICLES = 
       get_MOVER_TYPE()==SoA_vec_resort
    || get_MOVER_TYPE()==AoS_vec_resort;
  SORTING_PARTICLES = get_VECTORIZE_MOMENTS()
    || get_MOVER_TYPE()==SoA_vec_onesort
    || get_MOVER_TYPE()==AoS_vec_onesort
    || get_MOVER_TYPE()==SoA_vec_resort
    || get_MOVER_TYPE()==AoS_vec_resort;
  SORTING_SOA = get_VECTORIZE_MOMENTS()
    || get_MOVER_TYPE()==SoA_vec_onesort
    || get_MOVER_TYPE()==SoA_vec_resort;
  USING_AOS =
       get_MOMENTS_TYPE()==AoS
    || get_MOVER_TYPE()==AoS
    || get_MOVER_TYPE()==AoSintr
    || get_MOVER_TYPE()==AoS_vec_onesort
    || get_MOVER_TYPE()==AoS_vec_resort;
}

bool Parameters::get_RESORTING_PARTICLES() { return RESORTING_PARTICLES; }
bool Parameters::get_SORTING_PARTICLES() { return SORTING_PARTICLES; }
bool Parameters::get_SORTING_SOA() { return SORTING_SOA; }
bool Parameters::get_USING_AOS() { return USING_AOS; }

//bool Parameters::get_RESORTING_PARTICLES() { return true; }
//bool Parameters::get_SORTING_PARTICLES() { return true; }
//bool Parameters::get_SORTING_SOA() { return true; }
//bool Parameters::get_USING_AOS() { return true; }
//
bool Parameters::get_doWriteOutput() { return true; }
