#ifndef DATA
#define DATA

#include "global.dfn"

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

//=================



#include "dsmc.h"
#include "mol.h"
extern Cdsmc dsmc;
extern Cmol mol;
#endif
