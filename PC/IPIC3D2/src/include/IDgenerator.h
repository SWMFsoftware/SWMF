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
