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

// downloaded from https://gist.github.com/donny-dont/1471329
#ifndef aligned_allocator_h
#define aligned_allocator_h

//#include <cstdint>
//#include <stdint.h>
#include <cstddef> // for std::size_t
#include <stdexcept> // for std::length_error
//#include <stdlib.h> // for posix_memalign or aligned_alloc
//#include <malloc.h> // for memalign(alignment, size)
//#include <stdio.h> // for printf
#include "errors.h" // for eprintf

/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
template <typename T, std::size_t Alignment>
class aligned_allocator
{
  public:
 
    // The following will be the same for virtually all allocators.
    typedef T * pointer;
    typedef const T * const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    typedef std::size_t size_type;
    typedef ptrdiff_t difference_type;
 
    T * address(T& r) const
    {
      return &r;
    }
 
    const T * address(const T& s) const
    {
      return &s;
    }
 
    std::size_t max_size() const
    {
      // The following has been carefully written to be independent of
      // the definition of size_t and to avoid signed/unsigned warnings.
      return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
    }
 
 
    // The following must be the same for all allocators.
    template <typename U>
    struct rebind
    {
      typedef aligned_allocator<U, Alignment> other;
    } ;
 
    bool operator!=(const aligned_allocator& other) const
    {
      return !(*this == other);
    }
 
    void construct(T * const p, const T& t) const
    {
      void * const pv = static_cast<void *>(p);
 
      new (pv) T(t);
    }
 
    void destroy(T * const p) const
    {
      p->~T();
    }
 
    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const aligned_allocator& other) const
    {
      return true;
    }
 
 
    // Default constructor, copy constructor, rebinding constructor, and destructor.
    // Empty for stateless allocators.
    aligned_allocator() { }
 
    aligned_allocator(const aligned_allocator&) { }
 
    template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
 
    ~aligned_allocator() { }
 
 
    // The following will be different for each allocator.
    T * allocate(const std::size_t n) const
    {
      // The return value of allocate(0) is unspecified.
      // Mallocator returns NULL in order to avoid depending
      // on malloc(0)'s implementation-defined behavior
      // (the implementation can define malloc(0) to return NULL,
      // in which case the bad_alloc check below would fire).
      // All allocators can return NULL in this case.
      if (n == 0) {
        return NULL;
      }
 
      // All allocators should contain an integer overflow check.
      // The Standardization Committee recommends that std::length_error
      // be thrown in the case of integer overflow.
      if (n > max_size())
      {
        //throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
        eprintf("aligned_allocator<T>::allocate() - Integer overflow.");
        //abort();
      }
 
      // Mallocator wraps malloc().
#ifdef __INTEL_COMPILER
      void * const pv = _mm_malloc(n * sizeof(T), Alignment);
#else // is there a way to test whether posix_memalign is available?
      //void * const pv = new T[n];
      // this is supposed to be aligned with SSE types anyway...
      void * const pv = malloc(n*sizeof(T));
      // should I use posix_memalign()?
      //void * const pv = memalign(Alignment, n*sizeof(T));
      //void * pv;
      //posix_memalign(&pv, Alignment, n*sizeof(T));
#endif
 
      // Allocators should throw std::bad_alloc in the case of memory allocation failure.
      if (pv == NULL)
      {
        //throw std::bad_alloc();
        eprintf("std::bad_alloc()");
        //abort();
      }
 
      return static_cast<T *>(pv);
    }
 
    void deallocate(T * const p, const std::size_t n) const
    {
#ifdef __INTEL_COMPILER
      _mm_free(p);
#else
     free(p); // for memalign
#endif
    }
 
 
    // The following will be the same for all allocators that ignore hints.
    template <typename U>
    T * allocate(const std::size_t n, const U * /* const hint */) const
    {
      return allocate(n);
    }
 
 
    // Allocators are not required to be assignable, so
    // all allocators should have a private unimplemented
    // assignment operator. Note that this will trigger the
    // off-by-default (enabled under /Wall) warning C4626
    // "assignment operator could not be generated because a
    // base class assignment operator is inaccessible" within
    // the STL headers, but that warning is useless.
  private:
    aligned_allocator& operator=(const aligned_allocator&);
};

// unfortunately, C++ seems not to support typedefs for templates,
// so I use declaration macros.
// #include <vector> // needed if using aligned_vector
//#define aligned_vector(type) std::vector<type, aligned_allocator<type, 64> >

#endif // aligned_allocator_h
