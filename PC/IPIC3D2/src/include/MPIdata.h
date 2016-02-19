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

/***************************************************************************
  MPIdata.h  -  MPI data and methods wrapper
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef MPIDATA_H
#define MPIDATA_H

#ifndef NO_MPI
#include <mpi.h>
#include "my_mpi.h"
#endif

/**
 * MPI Data Structure. This class contains:
 *
 * - rank of process
 * - number of processors
 * - size of communication buffer
 *
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 *
 * I made this class a singleton.  It should only be created once,
 * since MPI_Init should be called only once. -Alec
 */
class MPIdata {
public:
  static MPIdata& instance();
private:
  // disable constructor and destructor of this singleton
  // by making them private.
  ~MPIdata(){}
  MPIdata(){}
public:
  /** initialize MPI environment */
  static void init(int *, char ***);
  /** Initialize for MHD-IPIC coupling */
  static void init(MPI_Comm iComm, signed int* iProc, signed int* nProc);
  /** close MPI environment */
  static void finalize_mpi();
  /** finalize and exit with error code */
  static void exit(int code);
  /** print MPI data structure */
  void Print(void);
  /** MPI status during the communication */
  //MPI_Status status;
public:
  static int get_rank(){return instance().rank;}
  static int get_nprocs(){return instance().nprocs;}
private:
  /** rank of the process */
  static int rank;
  /** number of processes */
  static int nprocs;

  // evidently unused...
  //char *buffer;
  //int buffer_size;
};

#define printf0(first, args...) \
  if(!MPIdata::get_rank()) printf(first, ## args);
#endif
