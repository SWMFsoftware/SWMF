//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//define the global variables for the whole execution code

#include "mpi.h"

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "global.dfn"

#ifndef _GLOBAL_VARIABLES_
#define _GLOBAL_VARIABLES_

extern MPI_Comm MPI_GLOBAL_COMMUNICATOR;

#endif

//compilation mode
#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_

//inlcude settings of the general block
#include "../../.general.conf"


#endif

