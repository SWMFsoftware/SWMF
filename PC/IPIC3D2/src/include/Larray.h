#ifndef Larray_h
#define Larray_h

#include "Alloc.h" // for ALIGNED
#include "ipicmath.h" // for pow2roundup
#include "asserts.h"
#include "string.h" // for memcpy

// linear array (e.g. of particles)
// 
// This basically extends aligned_vector(type),
// and the interface should be maintained in line with std::vector.
// Unfortunately, the mechanisms that C++ provides for extending
// the interface of a class are deficient, particularly in
// the case of templates.  Also, by avoiding implementation
// via inheritance from std::vector we also can specify the
// code at the lowest level, which is important for this
// performance-critical class.
//
template<class type>
class Larray
{
  static const int num_elem_in_block = 8;
 private: // members
  type* list;
  int _size; // number of particles in list
  int _capacity; // maximum number of particles
 private:
  void check_index(int i)const
  {
    #ifdef CHECK_BOUNDS
      assert_le(0, i);
      assert_lt(i, _size);
    #endif
  }
 public: // access on list
  int size()const { return _size; }
  int capacity()const { return _capacity; }

 public: // operations
  type& back(){ return list[_size-1]; }
  void pop_back()
  {
    assert(_size>0);
    _size--;
  }
  void clear()
  {
    _size = 0;
  }
  // this modifies std::vector::resize() by not initializing added elements
  //void resize(int newsize, value_type val = value_type())
  void resize(int newsize)
  {
    if(newsize <= _size)
    {
      _size = newsize;
    }
    else
    {
      reserve(newsize);
      _size = newsize;
    }
  }
  // unsafe version that assumes sufficient capacity
  // (extends std::vector)
  void fast_push_back(const type& element)
  {
    // remove this assertion to make this fast but unsafe
    assert_ge(_size,_capacity);
    list[_size] = element;
    _size++;
  }
  void push_back(const type& element)
  {
    if(__builtin_expect(_size>=_capacity,false))
    {
      int newcapacity = pow2roundup(_size+1);
      reserve(newcapacity);
    }
    list[_size] = element;
    _size++;
  }
  inline const type& operator[](int i)const
  {
    check_index(i);
    ALIGNED(list);
    return list[i];
  }
  inline type& operator[](int i)
  {
    check_index(i);
    ALIGNED(list);
    return list[i];
  }
  // this extends std::vector
  void delete_element(int i)
  {
    // works even for last particle,
    // though pop would be faster in that case.
    // 
    // why doesn't this compile?
    //this->operator[i] = list[--_size];
    check_index(i);
    ALIGNED(list);
    list[i] = list[--_size];
  }
 public: // memory
  ~Larray()
  {
    AlignedFree(list);
  }
  Larray():
    list(0),
    _size(0),
    _capacity(0)
  {}
  Larray(int requested_size):
    list(0),
    _size(0),
    _capacity(0)
  { if(requested_size > 0) reserve(requested_size); }
  // exchange content of this class with content of x
  void swap(Larray<type>& x)
  {
    // could do this with std::swap if willing to include
    // <utility> (since C++11) or <algorithm> (until C++11):
    //
    //std::swap(list,x.list);
    //std::swap(_size,x._size);
    //std::swap(_capacity,x._capacity);
    //
    type* tmp_list = list;
    int tmp_size = _size;
    int tmp_capacity = _capacity;
    list = x.list;
    _size = x._size;
    _capacity = x._capacity;
    x.list = tmp_list;
    x._size = tmp_size;
    x._capacity = tmp_capacity;
  }
  // request capacity to be at least newcapacity without deleting elements
  //
  // to bring this into conformity with std::vector, this should
  // be change so as never to shrink the capacity.
  void reserve(int newcapacity)
  {
    // ignore request if requested size is too small
    if(_size > newcapacity) return;

    // round up size to a multiple of num_elem_in_block
    //newcapacity = roundup_to_multiple(newcapacity,num_elem_in_block);
    newcapacity = ((newcapacity-1)/num_elem_in_block+1)*num_elem_in_block;
    if(newcapacity != _capacity)
    {
      _capacity = newcapacity;
      type* oldList = list;
      // the next two lines assume that type has no indirection
      list = AlignedAlloc(type,_capacity);
      memcpy(list,oldList,sizeof(type)*_size);
      // for(int i=0;i<_size;i++) list[i] = oldList[i];
      AlignedFree(oldList);
    }
  }
  // should rename this function as shrink_to_fit to conform to std::vector.
  void realloc_if_smaller_than(int required_max_size)
  {
    if(_size < required_max_size)
      reserve(required_max_size);
  }
  void shrink()
  {
    // shrink _capacity by a factor of two if elements will fit.
    int proposed_size = pow2rounddown(_capacity/2);
    if( _size <= proposed_size && proposed_size < _capacity)
    {
      reserve(proposed_size);
    }
  }
};

#endif
