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


// this uses std::vector, which
// must include about 8000 lines
//
//#include "aligned_allocator.h"
//#include <vector> // needed for aligned_vector
//#define aligned_vector(type) std::vector<type, aligned_allocator<type, 64> >
//
// this approximate implementation of std::vector 
// includes only about 2800 lines
//
#include "Larray.h"
// using declaration macro confuses ctags
#define aligned_vector(type) Larray<type>

// canonical workaround for lack of support in C++ for templated typedef
//template <typename T>
//struct aligned_vector
//{
//    typedef Larray<T> type;
//    //typedef std::vector<type, aligned_allocator<type, 64> > type;
//};
//
//// and yet another layer of indirection to avoid template brackets...
//class SpeciesParticle;
//typedef aligned_vector<SpeciesParticle>::type SpeciesParticleVector;
class SpeciesParticle;
typedef aligned_vector(SpeciesParticle) vector_SpeciesParticle;
typedef aligned_vector(double) vector_double;

