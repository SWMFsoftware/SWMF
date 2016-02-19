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
  VirtualTopology.h - Abstract Base class for virtual process topologies
  -------------------
begin                : May 2008
copyright            : KUL Leuven
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef VirtualTopology3D_H
#define VirtualTopology3D_H

//#include "mpi.h"
/**
 *  Abstract base class for virtual process topologies
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
//class VirtualTopology3D {
//public:
//  /** Find the neighbors in the new communicator  */
//  virtual void setup_vctopology(MPI_Comm comm_old) = 0;
//  /** Print topology info */
//  virtual void Print() = 0;
//  /** Print the mapping of topology */
//  virtual void PrintMapping() = 0;
//  /** get XLEN */
//  virtual int getXLEN() = 0;
//  /** get YLEN */
//  virtual int getYLEN() = 0;
//  /** get ZLEN */
//  virtual int getZLEN() = 0;
//  /** get nprocs */
//  virtual int getNprocs() = 0;
//  /** get periodicity on boundaries - DIRECTION X*/
//  virtual bool getPERIODICX() = 0;
//  /** get periodicity on boundaries - DIRECTION Y*/
//  virtual bool getPERIODICY() = 0;
//  /** get periodicity on boundaries - DIRECTION Z*/
//  virtual bool getPERIODICZ() = 0;
//  /** get the cartesian rank of the process */
//  virtual int getCartesian_rank() = 0;
//  /** get the cartesian rank of XLEFT neighbor */
//  virtual int getXleft_neighbor() = 0;
//  /** get the cartesian rank of XRIGHT neighbor */
//  virtual int getXright_neighbor() = 0;
//  /** get the cartesian rank of YLEFT neighbor */
//  virtual int getYleft_neighbor() = 0;
//  /** get the cartesian rank of YRIGHT neighbor */
//  virtual int getYright_neighbor() = 0;
//  /** get the cartesian rank of ZLEFT neighbor */
//  virtual int getZleft_neighbor() = 0;
//  /** get the cartesian rank of ZRIGHT neighbor */
//  virtual int getZright_neighbor() = 0;
//  /** get the PARTICLES cartesian rank of XLEFT neighbor */
//  virtual int getXleft_neighbor_P() = 0;
//  /** get the PARTICLES cartesian rank of XRIGHT neighbor */
//  virtual int getXright_neighbor_P() = 0;
//  /** get the PARTICLES cartesian rank of YLEFT neighbor */
//  virtual int getYleft_neighbor_P() = 0;
//  /** get the PARTICLES cartesian rank of YRIGHT neighbor */
//  virtual int getYright_neighbor_P() = 0;
//  /** get the PARTICLES cartesian rank of ZLEFT neighbor */
//  virtual int getZleft_neighbor_P() = 0;
//  /** get the PARTICLES cartesian rank of ZRIGHT neighbor */
//  virtual int getZright_neighbor_P() = 0;
//  /** get the coordinates in dir direction of process*/
//  virtual int getCoordinates(int dir) = 0;
//  /** get the coordinates of process*/
//  virtual int *getCoordinates() = 0;
//  /** get Periodicity condition in dir direction */
//  virtual int getPeriods(int dir) = 0;
//  /** if cVERBOSE == true, print to the screen all the comunication */
//  virtual bool getcVERBOSE() = 0;
//  /** get the MPI communicator */
//  virtual inline MPI_Comm getComm() = 0;
//};
#include "VCtopology3D.h"
#endif
