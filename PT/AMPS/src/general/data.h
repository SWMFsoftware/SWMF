#ifndef DATA
#define DATA

#include "global.dfn"

#define DSMCTARGET 0
#define PICTARGET  1
#define EULERTARGET 2

extern int DIM;
extern double tau,tmax,GeneralTime;
extern unsigned char NS;
extern int ThisThread;
extern int TotalThreadsNumber;

extern bool dsmc_flag;
extern bool chem_flag;
extern bool idf_flag;

extern int SymmetryMode;

extern bool ExternalSpeciesUsingFlag;


#ifndef CompilationTarget
#define CompilationTarget DSMCTARGET 
#endif

#if CompilationTarget==DSMCTARGET 
  #include "dsmc.h"
  #include "mol.h"

  extern Cdsmc dsmc;
  extern Cmol mol;
#endif

#if CompilationTarget==PICTARGET 
  #include "dsmc.h"
  #include "mol.h"
  #include "pic.h"

  extern Cpic pic;
  extern Cmol mol;
#endif

#if CompilationTarget==EULERTARGET
  #include "euler.h"

  extern CEuler euler;
#endif 



#endif
