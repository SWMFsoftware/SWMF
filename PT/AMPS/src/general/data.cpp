
#include "global.dfn"

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
#include "dsmc.h"
#include "mol.h"

Cdsmc dsmc;
Cmol mol;
