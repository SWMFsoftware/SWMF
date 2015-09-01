
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

