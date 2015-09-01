#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include "BlockCommunicator.h"
#include "IDgenerator.h"
#include "../main/Parameters.cpp"
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"
#include "../utility/MPIdata.cpp"
#include "../utility/IDgenerator.cpp"
// not sure on what systems this will be available
#include <ext/stdio_filebuf.h>
using namespace std;

#if 0
  Changes to make:
  * use persistent communication to reduce MPI overhead?
  * use buffered mode?  Consider arguments at:
    https://blogs.cisco.com/performance/top-10-reasons-why-buffered-sends-are-evil/
    https://www.cac.cornell.edu/VW/mpip2p/buffsend.aspx
#endif

//int MPIrank;
//int MPInumprocs;

#define debugout cout << "(" << MPIdata::get_rank() << ")" \
  << "DEBUG " << __func__ << ", " << __FILE__ << ":" << __LINE__ << ": "

// hacked mechanisms to print output one process at a time
//
// for some reason I still get threads overwriting one another,
// although not as badly as without this.
//
void output_barrier()
{
  fflush(stdout);
  cout << flush;
  MPI_Barrier(MPI_COMM_WORLD);
}
int barrier_ret_nprocs()
{
  output_barrier();
  return MPIdata::get_nprocs();
}
int barrier_ret_1()
{
  output_barrier();
  return 1;
}
//#define criticalout \
//  for(int rank=0;rank<barrier_ret_nprocs();rank++) \
//      if(rank==MPIdata::get_rank()) 
#define criticalout
#define masterout \
  for(int i=0;i<barrier_ret_1();i++) \
    if(!MPIdata::get_rank())
// Creating MPI barriers to atomize output is a poor solution;
#define dout criticalout debugout
// a better solution is to write output to separate files:
//ofstream outfile;
//#define dout outfile;

// particle ("element")
struct Particle
{
  double u[8];
};

std::ostream& operator<<(std::ostream& os, const Particle& pcl)
{
  //os << showpos;
  os << "[";
  for(int i=0; i<7; i++)
    os << pcl.u[i] << " ";
  os << pcl.u[7] << "]"; // << "\n";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const Block<Particle>& block_)
{
  const aligned_vector(Particle)& block = block_.get_block();
  for(int k=0; k<block.size();k++)
  {
    os << "\n  block[" << k << "] = " << block[k];
  }
  return os;
}

void test_particle_communication(
    BlockCommunicator<Particle>& lowXrecv,
    BlockCommunicator<Particle>& hghXrecv,
    BlockCommunicator<Particle>& lowXsend,
    BlockCommunicator<Particle>& hghXsend)
{
    // activate receiving
    //
    lowXrecv.recv_start();
    hghXrecv.recv_start();

    // make sure that the current block in each sender is ready
    lowXsend.send_start();
    hghXsend.send_start();

    // create, send, and receive particles
    for(int i=0; i<24;i++)
    {
      Particle lowPcl, hghPcl;
      for(int j=0;j<8;j++)
      {
        const double num = 100+i+.1*j;
        lowPcl.u[j] = num+.09;
        hghPcl.u[j] = num+.01;
      }
      hghPcl.u[0] = MPIdata::get_rank();
      lowPcl.u[0] = MPIdata::get_rank();
      bool lowXblockSent = lowXsend.send(lowPcl);
      bool hghXblockSent = hghXsend.send(hghPcl);
      if(0) // if(lowXblockSent || hghXblockSent)
      {
        // check if particles have arrived,
        // and if so deal with them
        MPI_Status status;
        if(lowXrecv.test_recv_curr_block(status))
        {
          Block<Particle>& recv_block = lowXrecv.fetch_received_block(status);
          dout << ": upon sending particle " << i
            << " received low particle block " << recv_block.get_id()
            // << ":" << recv_block
            << endl;
          lowXrecv.release_received_block();
        }
        if(hghXrecv.test_recv_curr_block(status))
        {
          Block<Particle>& recv_block = hghXrecv.fetch_received_block(status);
          dout << ": upon sending particle " << i
            << " received hgh particle block " << recv_block.get_id()
            //<< ":" << recv_block
            << endl;
          hghXrecv.release_received_block();
        }
      }
    }

    // send any remaining unsent particles
    //
    lowXsend.send_complete();
    hghXsend.send_complete();

    dprintf("--- finished sending ---");

    const int incount=2;
    MPI_Request recv_requests[incount] = 
    {
      lowXrecv.get_curr_request(),
      hghXrecv.get_curr_request()
    };
    // wait on and deal with remaining incoming blocks
    while(!(lowXrecv.comm_finished() && hghXrecv.comm_finished()))
    {
      // wait for some blocks to be received
      //
      //int recv_indices[incount]; // which requests completed
      //MPI_Status recv_statuses[incount]; // status of completed requests
      //int outcount; // number of requests that returned true
      //MPI_Waitsome(incount, recv_requests, &outcount, recv_indices, recv_statuses);
      //dprintf("recv.comm_finished() = [%d, %d]",
      //  lowXrecv.comm_finished(),
      //  hghXrecv.comm_finished());

      //#if 0
      int recv_index;
      MPI_Status recv_status;
      MPI_Waitany(incount, recv_requests, &recv_index, &recv_status);
      switch(recv_index)
      {
        default:
          invalid_value_error(recv_index);
        case MPI_UNDEFINED:
          eprintf("recv_requests contains no active handles");
          break;
        case 0: // lowXrecv
         {
          Block<Particle>& recv_block = lowXrecv.fetch_received_block(recv_status);
          dprintf("received lowXrecv.%d",recv_block.get_id());
          // dout << ": received lowXrecv." << recv_block.get_id() << recv_block << endl;
          lowXrecv.release_received_block();
          recv_requests[0] = lowXrecv.get_curr_request();
         }
          break;
        case 1: // hghXrecv
         {
          Block<Particle>& recv_block = hghXrecv.fetch_received_block(recv_status);
          dprintf("received hghXrecv.%d",recv_block.get_id());
          // dout << ": received hghXrecv." << recv_block.get_id() << recv_block << endl;
          hghXrecv.release_received_block();
          recv_requests[1] = hghXrecv.get_curr_request();
         }
          break;
      }
      //#endif

      //MPI_Status status;
      //if(!lowXrecv.at_end() && lowXrecv.test_recv_curr_block(status))
      //{
      //  Block<Particle>& recv_block = lowXrecv.fetch_received_block(status);
      //  dout << "(" << MPIdata::get_rank() << ") line "
      //    << __LINE__ << ": received lowXrecv."
      //    << recv_block.get_id() << ":"
      //    << recv_block << endl;
      //  lowXrecv.release_received_block();
      //}
      //if(!hghXrecv.at_end() && hghXrecv.test_recv_curr_block(status))
      //{
      //  Block<Particle>& recv_block = hghXrecv.fetch_received_block(status);
      //  dout << "(" << MPIdata::get_rank() << ") line "
      //    << __LINE__ << ": received hghXrecv."
      //    << recv_block.get_id() << ":"
      //    << recv_block << endl;
      //  hghXrecv.release_received_block();
      //}
    }
}

void test_particle_communication()
{
  //MPI_Init(&argc, &argv);
  //MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  //MPI_Comm_size(MPI_COMM_WORLD, &MPInumprocs);

  // define separate communicators for opposite
  // directions of information flow
  MPI_Comm up_comm;
  MPI_Comm dn_comm;
  const int reorder = 1;
  const int ndims = 1;
  /*const*/ int dims[ndims] = {2};
  /*const*/ int periods[ndims] = {1};
  const int XDIR = 0;
  const int UPWARD = 1;
  const int DNWARD = -1;
  // create a communicator for upward communication
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &up_comm);
  // create a communicator for downward communication
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &dn_comm);

  int up_cart_rank;
  int dn_cart_rank;
  MPI_Comm_rank(up_comm, &up_cart_rank);
  MPI_Comm_rank(dn_comm, &dn_cart_rank);
  criticalout dprint(up_cart_rank);
  criticalout dprint(dn_cart_rank);

  int up_coords[ndims];
  int dn_coords[ndims];
  MPI_Cart_coords(up_comm, up_cart_rank, ndims, up_coords);
  MPI_Cart_coords(dn_comm, dn_cart_rank, ndims, dn_coords);
  criticalout dprint(up_coords[0]);
  criticalout dprint(dn_coords[0]);

  //MPI_Comm_dup(up_comm, &dn_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &up_comm);
  //MPI_Comm_dup(MPI_COMM_WORLD, &dn_comm);

  // This does not actually do a shift; rather, it returns the
  // ranks of the neighbors that would be used to do this shift.
  // We use this as a mechanism to identify neighbors.
  // Shifting upward means that rank_source will be
  // lower_neighbor and rank_dest will be upper_neighbor.
  int up_src; // from lower
  int up_dst; // to upper
  int dn_src; // from upper
  int dn_dst; // to lower
  MPI_Cart_shift(up_comm, XDIR, UPWARD, &up_src, &up_dst);
  //MPI_Cart_shift(dn_comm, XDIR, DNWARD, &lower_neighbor, &upper_neighbor);
  MPI_Cart_shift(dn_comm, XDIR, DNWARD, &dn_src, &dn_dst);
  criticalout dprint(up_src); // from lower
  criticalout dprint(up_dst); // to upper
  criticalout dprint(dn_src); // from upper
  criticalout dprint(dn_dst); // to lower
  //Connection lowXrecvConn(up_src,2,up_comm);
  //Connection hghXrecvConn(dn_src,1,dn_comm);
  //Connection lowXsendConn(dn_dst,1,dn_comm);
  //Connection hghXsendConn(up_dst,2,up_comm);
  assert(up_src==dn_dst);
  assert(dn_src==up_dst);
  Connection lowXrecvConn(up_src,Direction::XUP);
  Connection hghXrecvConn(dn_src,Direction::XDN);
  Connection lowXsendConn(dn_dst,Direction::XDN);
  Connection hghXsendConn(up_dst,Direction::XUP);

  // showing that we can propagate a message upward
  if(0)
  {
    const int count=16;
    char recvbuff[count];
    char sendbuff[count];
    const int sendcount = 1+ // because of terminating '\0'
      snprintf(sendbuff,count,"hello from %d", MPIdata::get_rank());
    criticalout dprint(sendbuff);
    // receive message from lower
    MPI_Request recv_request;
    MPI_Irecv(recvbuff, count, MPI_CHAR, up_src, 0, up_comm, &recv_request);
    // send a message to upper
    MPI_Request send_request;
    MPI_Isend(sendbuff, sendcount, MPI_CHAR, up_dst, 0, up_comm, &send_request);
    // wait for message to arrive
    MPI_Status recv_status;
    MPI_Wait(&recv_request, &recv_status);
    // print the message
    criticalout dprint(recvbuff);
  }

  // showing that upward and downward messages are kept straight
  if(0)
  {
    const int count=16;
    char up_recvbuff[count];
    char up_sendbuff[count];
    char dn_recvbuff[count];
    char dn_sendbuff[count];
    const int up_sendcount = 1+ // because of terminating '\0'
      snprintf(up_sendbuff,count,"up_hello from %d", MPIdata::get_rank());
    const int dn_sendcount = 1+ // because of terminating '\0'
      snprintf(dn_sendbuff,count,"dn_hello from %d", MPIdata::get_rank());
    criticalout dprint(up_sendbuff);
    criticalout dprint(dn_sendbuff);
    // receive message from lower
    MPI_Request up_recv_request;
    MPI_Request dn_recv_request;
    MPI_Irecv(up_recvbuff, count, MPI_CHAR, up_src, 0, up_comm, &up_recv_request);
    MPI_Irecv(dn_recvbuff, count, MPI_CHAR, dn_src, 0, dn_comm, &dn_recv_request);
    // send a message to upper
    MPI_Request up_send_request;
    MPI_Request dn_send_request;
    MPI_Isend(up_sendbuff, up_sendcount, MPI_CHAR, up_dst, 0, up_comm, &up_send_request);
    MPI_Isend(dn_sendbuff, dn_sendcount, MPI_CHAR, dn_dst, 0, dn_comm, &dn_send_request);
    // wait for message to arrive
    MPI_Status up_recv_status;
    MPI_Status dn_recv_status;
    MPI_Wait(&up_recv_request, &up_recv_status);
    MPI_Wait(&dn_recv_request, &dn_recv_status);
    // print the message
    criticalout dprint(up_recvbuff);
    criticalout dprint(dn_recvbuff);
  }

  // showing that we can propagate a particle upward
  if(0)
  {
    const int count=8;
    Particle send_pcl;
    for(int i=0;i<count;i++)
      send_pcl.u[i]=i+.1*MPIdata::get_rank()+.09;
    criticalout debugout << "send_pcl = " << send_pcl << endl;
    // receive message from lower
    Particle recv_pcl;
    MPI_Request recv_request;
    MPI_Irecv(recv_pcl.u, count, MPI_DOUBLE, up_src, 0, up_comm, &recv_request);
    // send a message to upper
    MPI_Request send_request;
    MPI_Isend(send_pcl.u, count, MPI_DOUBLE, up_dst, 0, up_comm, &send_request);
    // wait for message to arrive
    MPI_Status recv_status;
    MPI_Wait(&recv_request, &recv_status);
    // print the message
    criticalout debugout << "recv_pcl = " << recv_pcl << endl;
  }

  // showing that we can propagate a block of particles upward
  if(0)
  {
    masterout debugout << "=== propagating particles upward ===" << endl;

    const int blocksize=8; //
    const int count=8*blocksize;
    //
    // receive message from lower
    //
    Block<Particle> recv_pcls(blocksize, MPIdata::get_rank());
    recv_pcls.recv(lowXrecvConn);
    //
    // send particles
    //
    Block<Particle> send_pcls(blocksize, MPIdata::get_rank());
    aligned_vector(Particle)& send_block = send_pcls.fetch_block();
    // initialize particles
    for(int p=0;p<blocksize;p++)
    {
      Particle pcl;
      for(int i=0;i<8;i++)
        pcl.u[i]=p + .1*i+.01*MPIdata::get_rank()+.009;
      send_block.push_back(pcl);
    }
    criticalout debugout << "send_pcls = " << send_pcls << endl;
    // send a message to upper
    send_pcls.send(hghXsendConn);
    // wait for message to arrive
    recv_pcls.waitfor_recv();
    // print the message one process at a time
    criticalout debugout << "recv_pcls = " << recv_pcls << endl;
  }

  // showing that we can propagate a list of blocks of particles upward
  if(0)
  {
    masterout debugout << "=== propagating particles upward ===" << endl;

    const int blocksize=4;
    const int numblocks = 2;
    const int numpcls = blocksize*numblocks-1;
    //
    // receive message from lower
    //
    BlockCommunicator<Particle> recv_pcls(lowXrecvConn, blocksize, numblocks);
    recv_pcls.recv_start();
    //
    // send particles
    //
    BlockCommunicator<Particle> send_pcls(hghXsendConn, blocksize, numblocks);
    // create and send particles
    int num_blocks_sent = 0;
    for(int p=0;p<numpcls;p++)
    {
      Particle pcl;
      for(int i=0;i<8;i++)
        pcl.u[i]=p + .1*i+.01*MPIdata::get_rank()+.009;
      num_blocks_sent += send_pcls.send(pcl);
    }
    send_pcls.send_complete();
    //assert_eq(num_blocks_sent, numblocks);
    // print all the blocks that were sent.
    //for(send_pcls.rewind();!send_pcls.at_end();send_pcls.advance_block())
    //{
    //  Block<Particle>& curr_block = send_pcls.fetch_curr_block();
    //  dout << curr_block << endl;
    //}

    // read and print blocks as they arrive
    {
      //recv_pcls.rewind();
      while(!recv_pcls.comm_finished())
      {
        MPI_Status status;
        MPI_Wait(&recv_pcls.fetch_curr_block().fetch_request(),&status);
        //recv_pcls.waitfor_recv_curr_block(status);
        Block<Particle>& curr_block = recv_pcls.fetch_received_block(status);
        dout << curr_block << endl;
        recv_pcls.release_received_block();
      }
    }
  }

  // communicating particles in blocks
  if(1)
  {
    masterout debugout << "=== propagating particles both ways ===" << endl;
    const int blocksize = 4;

    // for each neighbor, create a receive communicator
    //
    BlockCommunicator<Particle> lowXrecv(lowXrecvConn, blocksize, 1);
    BlockCommunicator<Particle> hghXrecv(hghXrecvConn, blocksize, 1);

    // for each neighbor, create a send communicator
    //
    BlockCommunicator<Particle> lowXsend(lowXsendConn, blocksize, 1);
    BlockCommunicator<Particle> hghXsend(hghXsendConn, blocksize, 1);

    lowXrecv.post_recvs();
    hghXrecv.post_recvs();

    test_particle_communication(lowXrecv,hghXrecv,lowXsend,hghXsend);

    dprintf("=== hey, that was great, let's try it again! ===");

    test_particle_communication(lowXrecv,hghXrecv,lowXsend,hghXsend);
  }

  // probably MPI_Finalize takes care of this anyway...
  MPI_Comm_free(&up_comm);
  MPI_Comm_free(&dn_comm);
}

void test_pcl_id_generator()
{
  const int nop=100;
  doubleIDgenerator pclIDgenerator;
  pclIDgenerator.reserve_particles_in_range(0,nop);
  dprint(pclIDgenerator.generateID());
  dprint(pclIDgenerator.generateID());
  dprint(pclIDgenerator.generateID());
}

int main(int argc, char **argv)
{
  MPIdata::init(&argc, &argv);

  test_particle_communication();
  test_pcl_id_generator();

  MPIdata::finalize_mpi();
  //MPI_Finalize();
}
