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

#include "mpi.h"
#include "Alloc.h"
#include "Collective.h"
#include "VCtopology3D.h"
#include <iostream>
#include "MPIdata.h"
#include "debug.h"

using std::cout;
using std::endl;

/** DEFINE THE Topology HERE, setting XLEN,YLEN,ZLEN */
VCtopology3D::VCtopology3D(const Collective& col) {
  // *******************************************
  // *******************************************
  // change these values to change the topology
  XLEN = col.getXLEN();
  YLEN = col.getYLEN();
  ZLEN = col.getZLEN();
  nprocs = XLEN * YLEN * ZLEN;
  // here you have to set the topology
  PERIODICX = col.getPERIODICX();
  PERIODICY = col.getPERIODICY();
  PERIODICZ = col.getPERIODICZ();

  PERIODICX_P = col.getPERIODICX_P();
  PERIODICY_P = col.getPERIODICY_P();
  PERIODICZ_P = col.getPERIODICZ_P();

  // *******************************************
  // *******************************************
  XDIR = 0;
  YDIR = 1;
  ZDIR = 2;
  RIGHT = 1;
  LEFT = -1;

  reorder = 1;

  dims[0] = XLEN;
  dims[1] = YLEN;
  dims[2] = ZLEN;

  periods[0] = PERIODICX;
  periods[1] = PERIODICY;
  periods[2] = PERIODICZ;

  periods_P[0] = PERIODICX_P;
  periods_P[1] = PERIODICY_P;
  periods_P[2] = PERIODICZ_P;

  cVERBOSE = false;             // communication verbose ?

}


/** Within CART_COMM, processes find about their new rank numbers, their cartesian coordinates,
  and their neighbors  */
void VCtopology3D::setup_vctopology(MPI_Comm old_comm) {
  // create a matrix with ranks, and neighbours for fields
  MPI_Cart_create(old_comm, 3, dims, periods, reorder, &CART_COMM);

  // create a matrix with ranks, and neighbours for Particles
  MPI_Cart_create(old_comm, 3, dims, periods_P, reorder, &CART_COMM_P);

  // field Communicator
  if (CART_COMM != MPI_COMM_NULL) {
    MPI_Comm_rank(CART_COMM, &cartesian_rank);
    MPI_Cart_coords(CART_COMM, cartesian_rank, 3, coordinates);

    MPI_Cart_shift(CART_COMM, XDIR, RIGHT, &xleft_neighbor, &xright_neighbor);
    MPI_Cart_shift(CART_COMM, YDIR, RIGHT, &yleft_neighbor, &yright_neighbor);
    MPI_Cart_shift(CART_COMM, ZDIR, RIGHT, &zleft_neighbor, &zright_neighbor);

  }
  else {
    // previous check that nprocs = XLEN*YLEN*ZLEN should prevent reaching this line.
    eprintf("A process is thrown away from the new topology for fields.");
  }
  // Particles Communicator
  //if (CART_COMM_P != MPI_COMM_NULL) {
  //  int pcl_coordinates[3];
  //  int pcl_cartesian_rank;
  //  MPI_Comm_rank(CART_COMM_P, &pcl_cartesian_rank);
  //  MPI_Cart_coords(CART_COMM_P, pcl_cartesian_rank, 3, pcl_coordinates);
  //  
  // This seems to be assumed elsewhere in the code.
  // We need to eliminate this assumption.
  // This becomes important if we want
  // a running program to have more than one
  // mesh and processor topology.
  // The MPI rank is for system-level code, e.g.
  // to identify the process from which debug is coming.
  // The cartesian rank is for application-level code.
  assert_eq(cartesian_rank, MPIdata::get_rank());


  if (CART_COMM_P != MPI_COMM_NULL) {
      MPI_Cart_shift(CART_COMM_P, XDIR, RIGHT, &xleft_neighbor_P, &xright_neighbor_P);
      MPI_Cart_shift(CART_COMM_P, YDIR, RIGHT, &yleft_neighbor_P, &yright_neighbor_P);
      MPI_Cart_shift(CART_COMM_P, ZDIR, RIGHT, &zleft_neighbor_P, &zright_neighbor_P);
  }

  _isPeriodicXlower_P = PERIODICX_P && (coordinates[0]==0);
  _isPeriodicXupper_P = PERIODICX_P && (coordinates[0]==dims[0]-1);
  _isPeriodicYlower_P = PERIODICY_P && (coordinates[1]==0);
  _isPeriodicYupper_P = PERIODICY_P && (coordinates[1]==dims[1]-1);
  _isPeriodicZlower_P = PERIODICZ_P && (coordinates[2]==0);
  _isPeriodicZupper_P = PERIODICZ_P && (coordinates[2]==dims[2]-1);

  _noXleftNeighbor_P = (xleft_neighbor_P == MPI_PROC_NULL);
  _noXrghtNeighbor_P = (xright_neighbor_P == MPI_PROC_NULL);
  _noYleftNeighbor_P = (yleft_neighbor_P == MPI_PROC_NULL);
  _noYrghtNeighbor_P = (yright_neighbor_P == MPI_PROC_NULL);
  _noZleftNeighbor_P = (zleft_neighbor_P == MPI_PROC_NULL);
  _noZrghtNeighbor_P = (zright_neighbor_P == MPI_PROC_NULL);

  _isBoundaryProcess_P =
    _noXleftNeighbor_P ||
    _noXrghtNeighbor_P ||
    _noYleftNeighbor_P ||
    _noYrghtNeighbor_P ||
    _noZleftNeighbor_P ||
    _noZrghtNeighbor_P;

}
/** destructor */
VCtopology3D::~VCtopology3D() {

}
/** print topology info */
void VCtopology3D::Print() {
  cout << endl;
  cout << "Virtual Cartesian Processors Topology" << endl;
  cout << "-------------------------------------" << endl;
  cout << "Processors grid: " << XLEN << "x" << YLEN << "x" << ZLEN << endl;
  cout << "Periodicity X: " << periods[0] << endl;
  cout << "Periodicity Y: " << periods[1] << endl;
  cout << "Periodicity Z: " << periods[2] << endl;
  cout << endl;
}
/** print cartesian rank of neighbors and coordinate of process */
void VCtopology3D::PrintMapping() {
  cout << endl;
  cout << "Mapping of process " << cartesian_rank << endl;
  cout << "----------------------" << endl;
  cout << "Coordinates: X = " << coordinates[0] << "; Y = " << coordinates[1] << "; Z = " << coordinates[2] << endl;
  cout << "Neighbors: xLeft = " << xleft_neighbor << "; xRight = " << xright_neighbor << "; yLeft = " << yleft_neighbor << "; yRight = " << yright_neighbor << "; zLeft = " << zleft_neighbor << "; zRight = " << zright_neighbor << endl;
  cout << endl;
}
