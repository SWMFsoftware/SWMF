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

#ifndef Connection_h
#define Connection_h

#include <mpi.h> // is there a way to forward declare mpi types?

// The combination of group (comm), tag, and neighbor
// should be unique for each connection.
//
// In a Cartesian topology, the challenges to this uniqueness
// occur for periodic boundary conditions when there is a
// dimension only 1 or 2 processes thick.
//
// In the 1-process-thick case, we can forgo MPI communication in
// that direction altogether.  But in the 2-process-thick case,
// we need to make a distinction between upward and downward
// channels of communication.  Distinguishing upward and downward
// directions is sufficient for any topology based on convex MPI
// subdomains, since in this case two subdomains can share at
// most two faces.  Any mechanism for distinguishing upward and
// downward directions for a pair of subdomains is sufficient;
// one way is to compare the positions of the centroids, first in
// the x direction, then in the y direction, and then in the z
// direction.  One might use tag=1 for downward communication and
// tag=2 for upward communication, and use MPI_COMM_MYSIM for the
// group.
//
namespace Direction
{
  enum Enum
  {
    DEFAULT = 0,
    PARTICLE_DN, // downward communication of particles
    PARTICLE_UP, // upward communication of particles
    XDN,
    XUP,
    YDN,
    YUP,
    ZDN,
    ZUP
  };
}
class Connection
{
 private: // data
  // In MPI a message envelope includes the following
  // information plus the rank of this process.
  int _rank; // rank within group of neighbors we're connecting to
  int _tag; // tag to attach to messages
  MPI_Comm _comm; // communicator group
 public: // init
  // construct a connection that creates self-communication
  // in place of null connection
  static Connection null2self(int rank_, int tag_, int self_tag, MPI_Comm comm_)
  {
    Connection c(rank_,tag_,comm_);
    if(c._rank==MPI_PROC_NULL)
    {
      c._rank = MPIdata::get_rank();
      c._tag = self_tag;
    }
    return c;
  }
  Connection():
    _rank(0),
    _tag(0),
    _comm(MPI_COMM_MYSIM)
  {}
  Connection(int rank_, int tag_, MPI_Comm comm_ = MPI_COMM_MYSIM)
  {
    init(rank_,tag_,comm_);
  }
 private:
  void init(int rank_, int tag_, MPI_Comm comm_ = MPI_COMM_MYSIM)
  {
    _rank = rank_;
    _tag = tag_;
    _comm = comm_;
  }
 public: // accessors
  int rank()const{return _rank;}
  int tag()const{return _tag;}
  MPI_Comm comm()const{return _comm;}
  static const char* tag_name(int tag)
  {
    switch(tag)
    {
      default: unsupported_value_error(tag);
      case Direction::XDN: return "XDN";
      case Direction::XUP: return "XUP";
      case Direction::YDN: return "YDN";
      case Direction::YUP: return "YUP";
      case Direction::ZDN: return "ZDN";
      case Direction::ZUP: return "ZUP";
    }
  }
  const char* tag_name()const
  {
    return tag_name(_tag);
  }
};

#endif
