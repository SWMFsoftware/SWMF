//===================================================
//$Id$
//===================================================


#ifdef MPI_ON
#include "mpi.h"
#endif


#include "global.dfn"
#include "data.h"

int DIM=3;
double tau,tmax,GeneralTime=0.0;
unsigned char NS=1;
int ThisThread=0;
int TotalThreadsNumber=1;

bool dsmc_flag=true;
bool chem_flag=false;
bool idf_flag=false;

int SymmetryMode=no_symmetry;

//External species
bool ExternalSpeciesUsingFlag=false;

//===================================================
//global objects

#if CompilationTarget==DSMCTARGET 
  #include "dsmc.h"
  #include "mol.h"
  Cdsmc dsmc;
  Cmol mol;
#endif 

#if CompilationTarget==PICTARGET  
  #include "dsmc.h"
  #include "mol.h"
  #include "pic.h"
  Cpic pic;
  Cmol mol;
#endif

#if CompilationTarget==EULERTARGET
  #include "euler.h"
  CEuler euler;
#endif

#if CompilationTarget==HYBRIDTARGET
  #include "dsmc.h"
  #include "mol.h"
  #include "euler.h"
  #include "hybrid.h"

  Cmol mol;
  Chybrid hybrid;
#endif




