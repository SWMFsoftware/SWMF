/*******************************************************************************************
  Particles3Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef Part3DCOMM_H
#define Part3DCOMM_H

//#include "CollectiveIO.h"
#include "ipicfwd.h"
#include "Alloc.h"
#include "Particle.h" // for ParticleType
// unfortunately this includes mpi.h, which includes 35000 lines:
#include "BlockCommunicator.h"
#include "aligned_vector.h"
#include "Larray.h"
#include "IDgenerator.h"
/**
 * 
 * class for particles of the same species with communications methods
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3Dcomm // :public Particles
{
public:
  /** constructor */
  Particles3Dcomm(int species, CollectiveIO * col,
    VirtualTopology3D * vct, Grid * grid);
  /** destructor */
  ~Particles3Dcomm();

  /** interpolation method GRID->PARTICLE order 1: CIC */
  // This does not belong in this class and is no longer in use.
  void interpP2G(Field * EMf);

 public: // handle boundary conditions
  // apply boundary conditions to all particles at the
  // end of a list of particles starting with index start
  // 
  // these are virtual so user can override these
  // to provide arbitrary custom boundary conditions
  virtual void apply_Xleft_BC(vector_SpeciesParticle& pcls, int start=0);
  virtual void apply_Yleft_BC(vector_SpeciesParticle& pcls, int start=0);
  virtual void apply_Zleft_BC(vector_SpeciesParticle& pcls, int start=0);
  virtual void apply_Xrght_BC(vector_SpeciesParticle& pcls, int start=0);
  virtual void apply_Yrght_BC(vector_SpeciesParticle& pcls, int start=0);
  virtual void apply_Zrght_BC(vector_SpeciesParticle& pcls, int start=0);
 private: // handle boundary conditions
  void apply_periodic_BC_global(vector_SpeciesParticle& pcl_list, int pstart);
  bool test_pcls_are_in_nonperiodic_domain(const vector_SpeciesParticle& pcls)const;
  bool test_pcls_are_in_domain(const vector_SpeciesParticle& pcls)const;
  bool test_outside_domain(const SpeciesParticle& pcl)const;
  bool test_outside_nonperiodic_domain(const SpeciesParticle& pcl)const;
  bool test_Xleft_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_x() < 0.; }
  bool test_Xrght_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_x() > Lx; }
  bool test_Yleft_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_y() < 0.; }
  bool test_Yrght_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_y() > Ly; }
  bool test_Zleft_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_z() < 0.; }
  bool test_Zrght_of_domain(const SpeciesParticle& pcl)
  { return pcl.get_z() > Lz; }
  void apply_nonperiodic_BCs_global(vector_SpeciesParticle&, int pstart);
  bool test_all_pcls_are_in_subdomain();
  void apply_BCs_globally(vector_SpeciesParticle& pcl_list);
  void apply_BCs_locally(vector_SpeciesParticle& pcl_list,
    int direction, bool apply_shift, bool do_apply_BCs);
 private: // communicate particles between processes
  void flush_send();
  bool send_pcl_to_appropriate_buffer(SpeciesParticle& pcl, int count[6]);
  int handle_received_particles(int pclCommMode=0);
 public:
  int separate_and_send_particles();
  void recommunicate_particles_until_done(int min_num_iterations=3);
  void communicate_particles();
  void pad_capacities();
 private:
  void resize_AoS(int nop);
  void resize_SoA(int nop);
  void copyParticlesToAoS();
  void copyParticlesToSoA();

 public:
  void convertParticlesToSynched();
  void convertParticlesToAoS();
  void convertParticlesToSoA();
  bool particlesAreSoA()const;

  /*! sort particles for vectorized push (needs to be parallelized) */
  //void sort_particles_serial_SoA_by_xavg();
  void sort_particles_serial();
  void sort_particles_serial_AoS();
  //void sort_particles_serial_SoA();

  // get accessors for optional arrays
  //
  //Larray<SpeciesParticle>& fetch_pcls(){ return _pcls; }
  //Larray<SpeciesParticle>& fetch_pclstmp(){ return _pclstmp; }

  // particle creation methods
  //
  void reserve_remaining_particle_IDs()
  {
    // reserve remaining particle IDs starting from getNOP()
    pclIDgenerator.reserve_particles_in_range(getNOP());
  }
  // create new particle
  void create_new_particle(
    double u, double v, double w, double q,
    double x, double y, double z)
  {
    const double t = pclIDgenerator.generateID();
    _pcls.push_back(SpeciesParticle(u,v,w,q,x,y,z,t));
  }
  // add particle to the list
  void add_new_particle(
    double u, double v, double w, double q,
    double x, double y, double z, double t)
  {
    _pcls.push_back(SpeciesParticle(u,v,w,q,x,y,z,t));
  }

  void delete_particle(int pidx)
  {
    _pcls[pidx]=_pcls.back();
    _pcls.pop_back();
    //_pcls.delete_element(pidx);
  }

  // inline get accessors
  //
  double get_dx(){return dx;}
  double get_dy(){return dy;}
  double get_dz(){return dz;}
  double get_invdx(){return inv_dx;}
  double get_invdy(){return inv_dy;}
  double get_invdz(){return inv_dz;}
  double get_xstart(){return xstart;}
  double get_ystart(){return ystart;}
  double get_zstart(){return zstart;}
  ParticleType::Type get_particleType()const { return particleType; }
  const SpeciesParticle& get_pcl(int pidx)const{ return _pcls[pidx]; }
  const vector_SpeciesParticle& get_pcl_list()const{ return _pcls; }
  const SpeciesParticle* get_pclptr(int id)const{ return &(_pcls[id]); }
  const double *getUall()  const { assert(particlesAreSoA()); return &u[0]; }
  const double *getVall()  const { assert(particlesAreSoA()); return &v[0]; }
  const double *getWall()  const { assert(particlesAreSoA()); return &w[0]; }
  const double *getQall()  const { assert(particlesAreSoA()); return &q[0]; }
  const double *getXall()  const { assert(particlesAreSoA()); return &x[0]; }
  const double *getYall()  const { assert(particlesAreSoA()); return &y[0]; }
  const double *getZall()  const { assert(particlesAreSoA()); return &z[0]; }
  static ID_field::Enum id_field()
  {
    return ID_field::Q; // q
    return ID_field::T; // t
  }
  const longid*getParticleIDall() const
  {
    longid* ret;
    switch(id_field())
    {
      default:
        invalid_value_error(id_field());
      case ID_field::Q:
        ret = (longid*) &q[0];
      case ID_field::T:
        ret = (longid*) &t[0];
    }
    return ret;
  }
  // accessors for particle with index indexPart
  //
  int getNOP()  const { return _pcls.size(); }
  // set particle components
  void setU(int i, double in){_pcls[i].set_u(in);}
  void setV(int i, double in){_pcls[i].set_v(in);}
  void setW(int i, double in){_pcls[i].set_w(in);}
  void setQ(int i, double in){_pcls[i].set_q(in);}
  void setX(int i, double in){_pcls[i].set_x(in);}
  void setY(int i, double in){_pcls[i].set_y(in);}
  void setZ(int i, double in){_pcls[i].set_z(in);}
  void setT(int i, double in){_pcls[i].set_t(in);}
  // fetch particle components
  double& fetchU(int i){return _pcls[i].fetch_u();}
  double& fetchV(int i){return _pcls[i].fetch_v();}
  double& fetchW(int i){return _pcls[i].fetch_w();}
  double& fetchQ(int i){return _pcls[i].fetch_q();}
  double& fetchX(int i){return _pcls[i].fetch_x();}
  double& fetchY(int i){return _pcls[i].fetch_y();}
  double& fetchZ(int i){return _pcls[i].fetch_z();}
  double& fetchT(int i){return _pcls[i].fetch_t();}
  // get particle components
  double getU(int i)const{return _pcls[i].get_u();}
  double getV(int i)const{return _pcls[i].get_v();}
  double getW(int i)const{return _pcls[i].get_w();}
  double getQ(int i)const{return _pcls[i].get_q();}
  double getX(int i)const{return _pcls[i].get_x();}
  double getY(int i)const{return _pcls[i].get_y();}
  double getZ(int i)const{return _pcls[i].get_z();}
  double getT(int i)const{return _pcls[i].get_t();}
  //int get_npmax() const {return npmax;}

  // computed get access
  //
  /** return the Kinetic energy */
  double getKe();
  /** return the maximum kinetic energy */
  double getMaxVelocity();
  /** return energy distribution */
  long long *getVelocityDistribution(int nBins, double maxVel);
  /** return the momentum */
  double getP();
  /** Print particles info: positions, velocities */
  void Print() const;
  /** Print the number of particles of this subdomain */
  void PrintNp() const;

public:
  // accessors
  //int get_ns()const{return ns;}
  // return number of this species
  int get_species_num()const{return ns;}
  int get_numpcls_in_bucket(int cx, int cy, int cz)const
  { return (*numpcls_in_bucket)[cx][cy][cz]; }
  int get_bucket_offset(int cx, int cy, int cz)const
  { return (*bucket_offset)[cx][cy][cz]; }

protected:
  // pointers to topology and grid information
  // (should be const)
  const Collective * col;
  const VirtualTopology3D * vct;
  const Grid * grid;
  //
  /** number of this species */
  int ns;
  /** maximum number of particles of this species on this domain. used for memory allocation */
  //int npmax;
  /** number of particles of this species on this domain */
  //int nop; // see getNOP();
  /** total number of particles */
  //long long np_tot;
  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;
  /** charge to mass ratio */
  double qom;
  /** recon thick */
  double delta;
  /** thermal velocity  - Direction X*/
  double uth;
  /** thermal velocity  - Direction Y*/
  double vth;
  /** thermal velocity  - Direction Z*/
  double wth;
  /** u0 Drift velocity - Direction X */
  double u0;
  /** v0 Drift velocity - Direction Y */
  double v0;
  /** w0 Drift velocity - Direction Z */
  double w0;
  // used to generate unique particle IDs
  doubleIDgenerator pclIDgenerator;

  ParticleType::Type particleType;
  //
  // AoS representation
  //
  //Larray<SpeciesParticle> _pcls;
  vector_SpeciesParticle _pcls;
  //
  // particles data
  //
  // SoA representation
  //
  // velocity components
  vector_double u;
  vector_double v;
  vector_double w;
  // charge
  vector_double q;
  // position
  vector_double x;
  vector_double y;
  vector_double z;
  // subcycle time
  vector_double t;
  // indicates whether this class is for tracking particles
  bool TrackParticleID;
  bool isTestParticle;
  double pitch_angle;
  double energy;

  // structures for sorting particles
  //
  /** Average position data (used during particle push) **/
  //
  //Larray<double>& _xavg;
  //Larray<double>& _yavg;
  //Larray<double>& _zavg;
  //
  // alternate temporary storage for sorting particles
  //
  //Larray<SpeciesParticle> _pclstmp;
  vector_SpeciesParticle _pclstmp;
  //
  // references for buckets for serial sort.
  //
  array3_int* numpcls_in_bucket;
  array3_int* numpcls_in_bucket_now; // accumulator used during sorting
  //array3_int* bucket_size; // maximum number of particles in bucket
  array3_int* bucket_offset;

  /** rank of processor in which particle is created (for ID) */
  int BirthRank[2];
  /** number of variables to be stored in buffer for communication for each particle  */
  int nVar;
  /** time step */
  double dt;
  //
  // Copies of grid data (should just put pointer to Grid in this class)
  //
  /** Simulation domain lengths */
  double xstart, xend, ystart, yend, zstart, zend, invVOL;
  /** Lx = simulation box length - x direction   */
  double Lx;
  /** Ly = simulation box length - y direction   */
  double Ly;
  /** Lz = simulation box length - z direction   */
  double Lz;
  /** grid spacings */
  double dx, dy, dz;
  /** number of grid nodes */
  int nxn, nyn, nzn;
  /** number of grid cells */
  int nxc, nyc, nzc;
  // convenience values from grid
  double inv_dx;
  double inv_dy;
  double inv_dz;
  //
  // Communication variables
  //
  /** buffers for communication */
  //
  // communicator for this species (duplicated from MPI_COMM_WORLD)
  MPI_Comm mpi_comm;
  // send buffers
  //
  BlockCommunicator<SpeciesParticle> sendXleft;
  BlockCommunicator<SpeciesParticle> sendXrght;
  BlockCommunicator<SpeciesParticle> sendYleft;
  BlockCommunicator<SpeciesParticle> sendYrght;
  BlockCommunicator<SpeciesParticle> sendZleft;
  BlockCommunicator<SpeciesParticle> sendZrght;
  //
  // recv buffers
  //
  BlockCommunicator<SpeciesParticle> recvXleft;
  BlockCommunicator<SpeciesParticle> recvXrght;
  BlockCommunicator<SpeciesParticle> recvYleft;
  BlockCommunicator<SpeciesParticle> recvYrght;
  BlockCommunicator<SpeciesParticle> recvZleft;
  BlockCommunicator<SpeciesParticle> recvZrght;

  /** bool for communication verbose */
  bool cVERBOSE;
  /** Boundary condition on particles:
          <ul>
          <li>0 = exit</li>
          <li>1 = perfect mirror</li>
          <li>2 = riemission</li>
          <li>3 = periodic condition </li>
          </ul>
          */
  /** Boundary Condition Particles: FaceXright */
  int bcPfaceXright;
  /** Boundary Condition Particles: FaceXleft */
  int bcPfaceXleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceYright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceYleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceZright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceZleft;
  //
  // Other variables
  //
  /** speed of light in vacuum */
  double c;
  /** restart variable for loading particles from restart file */
  int restart;
  /** Number of iteration of the mover*/
  int NiterMover;
  /** velocity of the injection of the particles */
  double Vinj;
  /** removed charge from species */
  double Q_removed;
  /** density of the injection of the particles */
  double Ninj;

 protected:

  // limits to apply to particle velocity
  //
  double umax;
  double vmax;
  double wmax;
  double umin;
  double vmin;
  double wmin;
};

// find the particles with particular IDs and print them
void print_pcls(vector_SpeciesParticle& pcls, int ns, longid* id_list, int num_ids);

typedef Particles3Dcomm Particles;

#endif
