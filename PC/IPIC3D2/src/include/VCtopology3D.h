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
  VCtopology3D.h  -  a 3D Virtual cartesian topology
  A virtual topology is a mechanism for naming the processes
  in a communicator in a way that fits the communication
  pattern better. Since our processes will communicate mainly
  with the nearest neighbours after the fashion of a two-dimensional
  grid, we create a virtual topology to reflect this fact
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL Luveun
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef VCtopology3D_H
#define VCtopology3D_H

#include "mpi.h"

/**
 *  
 * Virtual cartesian topology
 * A virtual topology is a mechanism for naming the processes
 * in a communicator in a way that fits the communication
 * pattern better. Since our processes will communicate mainly
 * with the nearest neighbours after the fashion of a two-dimensional
 * grid, we create a virtual topology to reflect this fact
 * @version 2.0
 */

class Collective;

class VCtopology3D //:public VirtualTopology3D
{
public:
  /** constructor: Define topology parameters: dimension, domain decomposition,... */
  VCtopology3D(const Collective& col);
  /** destructor */
  ~VCtopology3D();
  /** Find the neighbors in the new communicator  */
  void setup_vctopology(MPI_Comm comm_old);
  /** Print topology info */
  void Print();
  /** Print the mapping of topology */
  void PrintMapping();

  int getXLEN()const{ return (XLEN); }
  int getYLEN()const{ return (YLEN); }
  int getZLEN()const{ return (ZLEN); }
  int getNprocs()const{ return (nprocs); }
  bool getPERIODICX()const{ return (PERIODICX); }
  bool getPERIODICY()const{ return (PERIODICY); }
  bool getPERIODICZ()const{ return (PERIODICZ); }

  bool getPERIODICX_P()const{ return (PERIODICX_P); }
  bool getPERIODICY_P()const{ return (PERIODICY_P); }
  bool getPERIODICZ_P()const{ return (PERIODICZ_P); }

  //the below is for field communicator
  int getCartesian_rank()const{ return (cartesian_rank); }
  int getXleft_neighbor()const{ return (xleft_neighbor); }
  int getXright_neighbor()const{ return (xright_neighbor); }
  int getYleft_neighbor()const{ return (yleft_neighbor); }
  int getYright_neighbor()const{ return (yright_neighbor); }
  int getZleft_neighbor()const{ return (zleft_neighbor); }
  int getZright_neighbor()const{ return (zright_neighbor); }

//  bool hasXleftNeighbor()const{ return !_noXleftNeighbor; }
//  bool hasXrghtNeighbor()const{ return !_noXrghtNeighbor; }
//  bool hasYleftNeighbor()const{ return !_noYleftNeighbor; }
//  bool hasYrghtNeighbor()const{ return !_noYrghtNeighbor; }
//  bool hasZleftNeighbor()const{ return !_noZleftNeighbor; }
//  bool hasZrghtNeighbor()const{ return !_noZrghtNeighbor; }

  bool isXupper()const{ return (coordinates[0]==dims[0]-1); }
  bool isYupper()const{ return (coordinates[1]==dims[1]-1); }
  bool isZupper()const{ return (coordinates[2]==dims[2]-1); }

  //the below only called by particle
  int getXleft_neighbor_P()const{ return (xleft_neighbor_P); }
  int getXright_neighbor_P()const{ return (xright_neighbor_P); }
  int getYleft_neighbor_P()const{ return (yleft_neighbor_P); }
  int getYright_neighbor_P()const{ return (yright_neighbor_P); }
  int getZleft_neighbor_P()const{ return (zleft_neighbor_P); }
  int getZright_neighbor_P()const{ return (zright_neighbor_P); }

  bool isPeriodicXlower_P()const{ return _isPeriodicXlower_P; }
  bool isPeriodicXupper_P()const{ return _isPeriodicXupper_P; }
  bool isPeriodicYlower_P()const{ return _isPeriodicYlower_P; }
  bool isPeriodicYupper_P()const{ return _isPeriodicYupper_P; }
  bool isPeriodicZlower_P()const{ return _isPeriodicZlower_P; }
  bool isPeriodicZupper_P()const{ return _isPeriodicZupper_P; }

  bool noXleftNeighbor_P()const{ return _noXleftNeighbor_P; }
  bool noXrghtNeighbor_P()const{ return _noXrghtNeighbor_P; }
  bool noYleftNeighbor_P()const{ return _noYleftNeighbor_P; }
  bool noYrghtNeighbor_P()const{ return _noYrghtNeighbor_P; }
  bool noZleftNeighbor_P()const{ return _noZleftNeighbor_P; }
  bool noZrghtNeighbor_P()const{ return _noZrghtNeighbor_P; }

  bool hasXleftNeighbor_P()const{ return !_noXleftNeighbor_P; }
  bool hasXrghtNeighbor_P()const{ return !_noXrghtNeighbor_P; }
  bool hasYleftNeighbor_P()const{ return !_noYleftNeighbor_P; }
  bool hasYrghtNeighbor_P()const{ return !_noYrghtNeighbor_P; }
  bool hasZleftNeighbor_P()const{ return !_noZleftNeighbor_P; }
  bool hasZrghtNeighbor_P()const{ return !_noZrghtNeighbor_P; }

  bool isBoundaryProcess_P()const{ return _isBoundaryProcess_P; }


  bool getcVERBOSE()const{ return (cVERBOSE); }
  int getCoordinates(int dir)const{ return (coordinates[dir]); }
  const int *getCoordinates()const{ return (coordinates); }
  const int *getDims()const{ return dims; }
  int getPeriods(int dir)const{ return (periods[dir]); }
  MPI_Comm getFieldComm()const{ return (CART_COMM); }
  MPI_Comm getParticleComm()const{ return (CART_COMM_P); }


private:
  /** New communicator with virtual cartesian topology */
  MPI_Comm CART_COMM;
  /** New communicator with virtual cartesian topology for Particles*/
  MPI_Comm CART_COMM_P;
  /** MPI status during sending and receiving communication */
  MPI_Status status;
  /** Direction X for shift MPI_Cart_Shift*/
  int XDIR;
  /** Direction Y for shift MPI_Cart_Shift*/
  int YDIR;
  /** Direction Z for shift MPI_Cart_Shift*/
  int ZDIR;
  /** RIGHT = +    upwards   shift */
  int RIGHT;
  /** LEFT  = -    downwards shift */
  int LEFT;
  /** dimension of virtual topology */
  int PROCDIM;
  /** number of subdomains - Direction X */
  int XLEN;
  /** number of subdomains - Direction Y */
  int YLEN;
  /** number of subdomains - Direction Z */
  int ZLEN;
  /** nprocs = number of processors */
  int nprocs;
  /** periodicity on boundaries - DIRECTION X*/
  bool PERIODICX;
  bool PERIODICX_P;
  /** periodicity on boundaries - DIRECTION Y*/
  bool PERIODICY;
  bool PERIODICY_P;
  /** periodicity on boundaries - DIRECTION Z*/
  bool PERIODICZ;
  bool PERIODICZ_P;
  /** periodicity on boundaries - DIRECTION X*/
  //bool PERIODICX_P;
  /** periodicity on boundaries - DIRECTION Y*/
  //bool PERIODICY_P;
  /** periodicity on boundaries - DIRECTION Z*/
  //bool PERIODICZ_P;
  /** rank may be reordered     */
  int reorder;
  /** arrays for Create_Cart_create  */
  int dims[3]; // i.e. divisions
  /** periodicity */
  int periods[3];
  int periods_P[3];
  /** coordinates on processors grid */
  int coordinates[3];
  /** cartesian rank */
  int cartesian_rank;
  /** cartesian rank of XLEFT neighbor */
  int xleft_neighbor;
  int xleft_neighbor_P;
  /** cartesian rank of XRIGHT neighbor */
  int xright_neighbor;
  int xright_neighbor_P;
  /** cartesian rank of YLEFT neighbor */
  int yleft_neighbor;
  int yleft_neighbor_P;
  /** cartesian rank of YRIGHT neighbor */
  int yright_neighbor;
  int yright_neighbor_P;
  /** cartesian rank of ZRIGHT neighbor */
  int zleft_neighbor;
  int zleft_neighbor_P;
  /** cartesian rank of ZLEFT neighbor */
  int zright_neighbor;
  int zright_neighbor_P;

  
  /**  for Field Communicator **/
  bool _noXrghtNeighbor;
  bool _noXleftNeighbor;
  bool _noYrghtNeighbor;
  bool _noYleftNeighbor;
  bool _noZrghtNeighbor;
  bool _noZleftNeighbor;

  /**  for Particle Communicator **/
  bool _isPeriodicXlower_P;
  bool _isPeriodicXupper_P;
  bool _isPeriodicYlower_P;
  bool _isPeriodicYupper_P;
  bool _isPeriodicZlower_P;
  bool _isPeriodicZupper_P;

  bool _noXrghtNeighbor_P;
  bool _noXleftNeighbor_P;
  bool _noYrghtNeighbor_P;
  bool _noYleftNeighbor_P;
  bool _noZrghtNeighbor_P;
  bool _noZleftNeighbor_P;

  int _isBoundaryProcess_P;

  /** if cVERBOSE == true, print to the screen all the comunication */
  bool cVERBOSE;
};

typedef VCtopology3D VirtualTopology3D;

#endif
