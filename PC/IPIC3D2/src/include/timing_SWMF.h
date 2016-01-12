#include <string>
#ifndef TIMING_SWMF_H
#define TIMING_SWMF_H

extern "C" {
  // In SWMF/util/TIMING/src/timing.f90
  void timing_start_c(size_t* nameLen, char* name);
  void timing_stop_c(size_t* nameLen, char* name);
}

inline void timing_start(string name){
#ifdef BATSRUS
  size_t nameLen;
  nameLen = name.length();
  char *nameChar=new char[nameLen+1];
  name.copy(nameChar,nameLen);
  timing_start_c(&nameLen, nameChar);
#endif
}


inline void timing_stop(string name){
#ifdef BATSRUS
  size_t nameLen;
  nameLen = name.length();
  char *nameChar=new char[nameLen+1];
  name.copy(nameChar,nameLen);
  timing_stop_c(&nameLen, nameChar);
#endif
}

#endif
