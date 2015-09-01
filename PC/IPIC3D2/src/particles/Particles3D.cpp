/*******************************************************************************************
  Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/


#include <mpi.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include "asserts.h"

#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "CollectiveIO.h"
#include "Collective.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "Field.h"
#include "MPIdata.h"
#include "ipicdefs.h"
#include "TimeTasks.h"
#include "parallel.h"
#include "Particles3D.h"

#include "mic_particles.h"
#include "debug.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
// particles processed together
#define P_SAME_TIME 2

// if true then mover cannot move particle 
// more than one processor subdomain
//
static bool cap_velocity(){return false;}

/**
 * 
 * Class for particles of the same species
 * @date Fri Jun 4 2009
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** particles are uniformly distributed with zero velocity   */
void Particles3D::uniform_background(Field * EMf)
{
  for (int i = 1; i < grid->getNXC() - 1; i++)
  for (int j = 1; j < grid->getNYC() - 1; j++)
  for (int k = 1; k < grid->getNZC() - 1; k++)
  {
    for (int ii = 0; ii < npcelx; ii++)
    for (int jj = 0; jj < npcely; jj++)
    for (int kk = 0; kk < npcelz; kk++)
    {
      double x = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
      double y = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
      double z = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
      double u = 0.0;
      double v = 0.0;
      double w = 0.0;
      double q = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
      _pcls.push_back(SpeciesParticle(u,v,w,q,x,y,z,0));
    }
  }
  cout << "Velocity Maxwellian Distribution " << endl;
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
  <ul>
  <li> dim = 0 --> constant velocity on X direction </li>
  <li> dim = 1 --> constant velocity on Y direction </li>
  <li> dim = 2 --> constant velocity on Z direction </li>
  </ul>

*/
void Particles3D::constantVelocity(double vel, int dim, Field * EMf) {
  switch (dim) {
    case 0:
      for (int i = 0; i < getNOP(); i++)
      {
        setU(i,vel);
        setV(i,0.);
        setW(i,0.);
      }
      break;
    case 1:
      for (int i = 0; i < getNOP(); i++)
      {
        setU(i,0.);
        setV(i,vel);
        setW(i,0.);
      }
      break;
    case 2:
      for (int i = 0; i < getNOP(); i++)
      {
        setU(i,0.);
        setV(i,0.);
        setW(i,vel);
      }
      break;

  }

}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles3D::alt_maxwellian(Field * EMf) {
  eprintf("unimplemented");
}

#ifdef BATSRUS
/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::MaxwellianFromFluid(Field* EMf,Collective *col, int is){

  /*
   * Constuctiong the distrebution function from a Fluid model
   */

  // loop over grid cells and set position, velociy and charge of all particles indexed by counter
  // there are multiple (27 or so) particles per grid cell.
  int i,j,k,counter=0;
  for (i=1; i< grid->getNXC()-1;i++)
    for (j=1; j< grid->getNYC()-1;j++)
      for (k=1; k< grid->getNZC()-1;k++)
        MaxwellianFromFluidCell(col,is, i,j,k,counter,x,y,z,q,u,v,w,ParticleID);
}

void Particles3D::MaxwellianFromFluidCell(Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, longid* ParticleID)
{
  /*
   * grid           : local grid object (in)
   * col            : collective (global) object (in)
   * is             : species index (in)
   * i,j,k          : grid cell index on proc (in)
   * ip             : particle number counter (inout)
   * x,y,z          : particle position (out)
   * q              : particle charge (out)
   * vx,vy,vz       : particle velocity (out)
   * ParticleID     : particle tracking ID (out)
   */

  // loop over particles inside grid cell i,j,k
  for (int ii=0; ii < npcelx; ii++)
    for (int jj=0; jj < npcely; jj++)
      for (int kk=0; kk < npcelz; kk++){
        // Assign particle positions: uniformly spaced. x_cellnode + dx_particle*(0.5+index_particle)
        fetchX(ip) = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);
        fetchY(ip) = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
        fetchZ(ip) = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
        // q = charge
        fetchQ(ip) =  (qom/fabs(qom))*(col->getFluidRhoCenter(i,j,k,is)/npcel)*(1.0/grid->getInvVOL());
        // u = X velocity
        sample_maxwellian(
          fetchU(ip),fetchV(ip),fetchW(ip),
          col->getFluidUthx(i,j,k,is),
          col->getFluidVthx(i,j,k,is),
          col->getFluidWthx(i,j,k,is),
          col->getFluidUx(i,j,k,is),
          col->getFluidVx(i,j,k,is),
          col->getFluidWx(i,j,k,is));
        ip++ ;
      }
}
#endif

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Field * EMf)
{
  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  assert_eq(_pcls.size(),0);

  const double q_sgn = (qom / fabs(qom));
  // multipled by charge density gives charge per particle
  const double q_factor =  q_sgn * grid->getVOL() / npcel;

  for (int i = 1; i < grid->getNXC() - 1; i++)
  {
  for (int j = 1; j < grid->getNYC() - 1; j++)
  for (int k = 1; k < grid->getNZC() - 1; k++)
  {
    const double q = q_factor * EMf->getRHOcs(i, j, k, ns);
    for (int ii = 0; ii < npcelx; ii++)
    for (int jj = 0; jj < npcely; jj++)
    for (int kk = 0; kk < npcelz; kk++)
    {
      double u,v,w;
      sample_maxwellian(
        u,v,w,
        uth, vth, wth,
        u0, v0, w0);
      // could also sample positions randomly as in repopulate_particles();
      const double x = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
      const double y = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
      const double z = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
      create_new_particle(u,v,w,q,x,y,z);
    }
  }
  //dprintf("capacity=%d, size=%d", _pcls.capacity(), getNOP());
  }
  if(0)
  {
    dprintf("number of particles of species %d: %d", ns, getNOP());
    const int num_ids = 1;
    longid id_list[num_ids] = {0};
    print_pcls(_pcls,ns,id_list, num_ids);
  }
}


/** pitch_angle_energy initialization (Assume B on z only) for test particles */
void Particles3D::pitch_angle_energy(Field * EMf) {

    /* initialize random generator with different seed on different processor */
    srand(vct->getCartesian_rank() + 3 + ns);
    assert_eq(_pcls.size(),0);

    double p0, pperp0, gyro_phase;

    const double q_factor =  (qom / fabs(qom)) * grid->getVOL() / npcel;

    long long counter=0;

    for (int i=1; i< grid->getNXC()-1;i++)
        for (int j=1; j< grid->getNYC()-1;j++)
            for (int k=1; k< grid->getNZC()-1;k++){

            	// q = charge following electron (species 0)
            	const double q = q_factor * EMf->getRHOcs(i, j, k, 0);

                for (int ii=0; ii < npcelx; ii++)
                    for (int jj=0; jj < npcely; jj++)
                        for (int kk=0; kk < npcelz; kk++){
                        	const double x= (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);
                        	const double y= (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
                        	const double z= (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);

                            // velocity - assumes B is along z
                            p0=sqrt((energy+1)*(energy+1)-1);
                            const double w =p0*cos(pitch_angle);
                            pperp0=p0*sin(pitch_angle);
                            gyro_phase = 2*M_PI* rand()/(double)RAND_MAX;
                            const double u=pperp0*cos(gyro_phase);
                            const double v=pperp0*sin(gyro_phase);
                            counter++ ;

                            create_new_particle(u,v,w,q,x,y,z);
                        }
            }
    const int num_ids = 1;
    longid id_list[num_ids] = {0};
    if (vct->getCartesian_rank() == 0){
    	cout << "------------------------------------------" << endl;
        cout << "Initialize Test Particle "<< ns << " with pitch angle "<< pitch_angle << ", energy " << energy << ", qom " << qom << ", npcel "<< counter<< endl;
        cout << "------------------------------------------" << endl;
    }
}


/** Force Free initialization (JxB=0) for particles */
void Particles3D::force_free(Field * EMf)
{
  eprintf("this function was not properly implemented and needs to be revised.");
#if 0
  /* initialize random generator */
  srand(vct->getCartesian_rank() + 1 + ns);
  for (int i = 1; i < grid->getNXC() - 1; i++)
  for (int j = 1; j < grid->getNYC() - 1; j++)
  for (int k = 1; k < grid->getNZC() - 1; k++)
  {
    for (int ii = 0; ii < npcelx; ii++)
    for (int jj = 0; jj < npcely; jj++)
    for (int kk = 0; kk < npcelz; kk++)
    {
      double x = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
      double y = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
      double z = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
      // q = charge
      double q = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / invVOL);
      double shaperx = tanh((y - Ly / 2) / delta) / cosh((y - Ly / 2) / delta) / delta;
      double shaperz = 1.0 / (cosh((y - Ly / 2) / delta) * cosh((y - Ly / 2) / delta)) / delta;
      eprintf("shapery needs to be initialized.");
      eprintf("flvx etc. need to be initialized.");
      double shapery;
      // new drift velocity to satisfy JxB=0
      const double flvx = u0 * flvx * shaperx;
      const double flvz = w0 * flvz * shaperz;
      const double flvy = v0 * flvy * shapery;
      double u = c;
      double v = c;
      double w = c;
      while ((fabs(u) >= c) || (fabs(v) >= c) || (fabs(w) >= c))
      {
        sample_maxwellian(
          u, v, w,
          uth, vth, wth,
          flvx, flvy, flvz);
      }
      create_new_particle(u,v,w,q,x,y,z);
    }
  }
#endif
}

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0) {

  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);
  jx_mod *= alpha;
  jy_mod *= alpha;
  jz_mod *= alpha;
  for (int i = 0; i < getNOP(); i++) {
    fetchU(i) += jx_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jx_phase);
    fetchV(i) += jy_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jy_phase);
    fetchW(i) += jz_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jz_phase);
  }
}

/** explicit mover */
void Particles3D::mover_explicit(Field * EMf) {
  eprintf("unimplemented");
}
//
// Create a vectorized version of this mover as follows.
//
// Let N be the number of doubles that fit in the vector unit.
//
// Vectorization can be done in one of two ways:
// 1. trivially by writing loops over stride-1 data,
//    particularly if data is in SoA format, and
// 2. with intrinsics, e.g. by loading the components
//    of physical vectors in the vector unit and
//    performing physical vector operations.
//
// Process N:=sizeof(vector_unit)/sizeof(double) particles at a time:
//
// A. Initialize average position with old position.
// B. repeat the following steps:
//    1. Use average position to generate 8 node weights and indices:
//       for pcl=1:N: (3) position -> (8) weights and (8) node_indices
//    2. Use weights and fields at node indices to calculate sampled field:
//       for pcl=1:N: (8) weights, field at nodes -> (6 or 8) sampled field
//    3. Use sampled field to advance velocity:
//       for pcl=1:N: sampled field, velocity -> (3) new velocity
//    4. Use velocity to advance position:
//       for pcl=1:N: velocity, old position -> (3) average position
// C. Use old position and average position to update the position.
//
// Since the fields used by the pusher are in AoS format,
// step 2 can be vectorized trivially, independent of whether
// the particle representation is AoS or SoA, and is expected
// to dominate in the absence of vectorization.
//
// If the particle representation is AoS, then
// vectorizing steps 1, 3, and 4 requires intrinsics.
//
// If the particle representation is SoA then
// vectorizing step 3 requires intrinsics unless:
// 1. the Nx(6 or 8) array of sampled fields is transposed
//    from AoS format to SoA format or
// 2. particles are sorted by mesh cell and are subcycled
//    or resorted with each iteration to keep the average
//    position within the mesh cell;
//    in this case the cell nodes are the same for all
//    N particles, so the sampled fields can be directly
//    calculated in SoA format.
//
// Compare the vectorization notes at the top of sumMoments()
//
/** mover with a Predictor-Corrector scheme */
void Particles3D::mover_PC(Field * EMf) {
#pragma omp parallel
{
  convertParticlesToSoA();
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "***PC MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_double fieldForPcls = EMf->get_fieldForPcls();

  
  const double dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
  #pragma omp for schedule(static)
  // why does single precision make no difference in execution speed?
  //#pragma simd vectorlength(VECTOR_WIDTH)
  for (int pidx = 0; pidx < getNOP(); pidx++) {
    // copy the particle
    const double xorig = getX(pidx);
    const double yorig = getY(pidx);
    const double zorig = getZ(pidx);
    const double uorig = getU(pidx);
    const double vorig = getV(pidx);
    const double worig = getW(pidx);
    double xavg = xorig;
    double yavg = yorig;
    double zavg = zorig;
    double uavg;
    double vavg;
    double wavg;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {
      // interpolation G-->P
      const double ixd = floor((xavg - xstart) * inv_dx);
      const double iyd = floor((yavg - ystart) * inv_dy);
      const double izd = floor((zavg - zstart) * inv_dz);
      // interface of index to right of cell
      int ix = 2 + int(ixd);
      int iy = 2 + int(iyd);
      int iz = 2 + int(izd);

      // use field data of closest cell in domain
      //
      if (ix < 1) ix = 1;
      if (iy < 1) iy = 1;
      if (iz < 1) iz = 1;
      if (ix > nxc) ix = nxc;
      if (iy > nyc) iy = nyc;
      if (iz > nzc) iz = nzc;
      // index of cell of particle;
      const int cx = ix - 1;
      const int cy = iy - 1;
      const int cz = iz - 1;

      const double xi0   = xavg - grid->getXN(ix-1);
      const double eta0  = yavg - grid->getYN(iy-1);
      const double zeta0 = zavg - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - xavg;
      const double eta1  = grid->getYN(iy) - yavg;
      const double zeta1 = grid->getZN(iz) - zavg;

      double weights[8] ALLOC_ALIGNED;
      const pfloat weight0 = invVOL*xi0;
      const pfloat weight1 = invVOL*xi1;
      const pfloat weight00 = weight0*eta0;
      const pfloat weight01 = weight0*eta1;
      const pfloat weight10 = weight1*eta0;
      const pfloat weight11 = weight1*eta1;
      weights[0] = weight00*zeta0; // weight000
      weights[1] = weight00*zeta1; // weight001
      weights[2] = weight01*zeta0; // weight010
      weights[3] = weight01*zeta1; // weight011
      weights[4] = weight10*zeta0; // weight100
      weights[5] = weight10*zeta1; // weight101
      weights[6] = weight11*zeta0; // weight110
      weights[7] = weight11*zeta1; // weight111
      //weights[0] = xi0 * eta0 * zeta0 * qi * invVOL; // weight000
      //weights[1] = xi0 * eta0 * zeta1 * qi * invVOL; // weight001
      //weights[2] = xi0 * eta1 * zeta0 * qi * invVOL; // weight010
      //weights[3] = xi0 * eta1 * zeta1 * qi * invVOL; // weight011
      //weights[4] = xi1 * eta0 * zeta0 * qi * invVOL; // weight100
      //weights[5] = xi1 * eta0 * zeta1 * qi * invVOL; // weight101
      //weights[6] = xi1 * eta1 * zeta0 * qi * invVOL; // weight110
      //weights[7] = xi1 * eta1 * zeta1 * qi * invVOL; // weight111

      // creating these aliases seems to accelerate this method by about 30%
      // on the Xeon host, processor, suggesting deficiency in the optimizer.
      //
      const double* field_components[8];
      get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);

      //double Exl=0,Exl=0,Ezl=0,Bxl=0,Byl=0,Bzl=0;
      //for(int c=0; c<8; c++)
      //{
      //  Bxl += weights[c] * field_components[c][0];
      //  Byl += weights[c] * field_components[c][1];
      //  Bzl += weights[c] * field_components[c][2];
      //  Exl += weights[c] * field_components[c][0+DFIELD_3or4];
      //  Eyl += weights[c] * field_components[c][1+DFIELD_3or4];
      //  Ezl += weights[c] * field_components[c][2+DFIELD_3or4];
      //}
      // causes compile error on icpc
      //double sampled_field[8]={0,0,0,0,0,0,0,0} ALLOC_ALIGNED;
      double sampled_field[8] ALLOC_ALIGNED;
      for(int i=0;i<8;i++) sampled_field[i]=0;
      double& Bxl=sampled_field[0];
      double& Byl=sampled_field[1];
      double& Bzl=sampled_field[2];
      double& Exl=sampled_field[0+DFIELD_3or4];
      double& Eyl=sampled_field[1+DFIELD_3or4];
      double& Ezl=sampled_field[2+DFIELD_3or4];
      const int num_field_components=2*DFIELD_3or4;
      for(int c=0; c<8; c++)
      {
        const double* field_components_c=field_components[c];
        ASSUME_ALIGNED(field_components_c);
        const double weights_c = weights[c];
        #pragma simd
        for(int i=0; i<num_field_components; i++)
        {
          sampled_field[i] += weights_c*field_components_c[i];
        }
      }
      const double Omx = qdto2mc*Bxl;
      const double Omy = qdto2mc*Byl;
      const double Omz = qdto2mc*Bzl;

      // end interpolation
      const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const pfloat denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const pfloat ut = uorig + qdto2mc * Exl;
      const pfloat vt = vorig + qdto2mc * Eyl;
      const pfloat wt = worig + qdto2mc * Ezl;
      //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation 
      uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
      // update average position
      xavg = xorig + uavg * dto2;
      yavg = yorig + vavg * dto2;
      zavg = zorig + wavg * dto2;
    }                           // end of iteration
    // update the final position and velocity
    fetchX(pidx) = xorig + uavg * dt;
    fetchY(pidx) = yorig + vavg * dt;
    fetchZ(pidx) = zorig + wavg * dt;
    fetchU(pidx) = 2.0 * uavg - uorig;
    fetchV(pidx) = 2.0 * vavg - vorig;
    fetchW(pidx) = 2.0 * wavg - worig;
  }                             // END OF ALL THE PARTICLES
}
}

void Particles3D::mover_PC_AoS(Field * EMf)
{
#pragma omp parallel
{
  convertParticlesToAoS();
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "***AoS MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  
  const double dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
  #pragma omp for schedule(static)
  for (int pidx = 0; pidx < getNOP(); pidx++) {
    // copy the particle
    SpeciesParticle* pcl = &_pcls[pidx];
    ALIGNED(pcl);
    const double xorig = pcl->get_x();
    const double yorig = pcl->get_y();
    const double zorig = pcl->get_z();
    const double uorig = pcl->get_u();
    const double vorig = pcl->get_v();
    const double worig = pcl->get_w();
    double xavg = xorig;
    double yavg = yorig;
    double zavg = zorig;
    double uavg;
    double vavg;
    double wavg;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {

      // compute weights for field components
      //
      double weights[8] ALLOC_ALIGNED;
      int cx,cy,cz;
      grid->get_safe_cell_and_weights(xavg,yavg,zavg,cx,cy,cz,weights);

      const double* field_components[8] ALLOC_ALIGNED;
      get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);

      //double Exl=0,Exl=0,Ezl=0,Bxl=0,Byl=0,Bzl=0;
      //for(int c=0; c<8; c++)
      //{
      //  Bxl += weights[c] * field_components[c][0];
      //  Byl += weights[c] * field_components[c][1];
      //  Bzl += weights[c] * field_components[c][2];
      //  Exl += weights[c] * field_components[c][0+DFIELD_3or4];
      //  Eyl += weights[c] * field_components[c][1+DFIELD_3or4];
      //  Ezl += weights[c] * field_components[c][2+DFIELD_3or4];
      //}
      // causes compile error on icpc
      //double sampled_field[8]={0,0,0,0,0,0,0,0} ALLOC_ALIGNED;
      double sampled_field[8] ALLOC_ALIGNED;
      for(int i=0;i<8;i++) sampled_field[i]=0;
      double& Bxl=sampled_field[0];
      double& Byl=sampled_field[1];
      double& Bzl=sampled_field[2];
      double& Exl=sampled_field[0+DFIELD_3or4];
      double& Eyl=sampled_field[1+DFIELD_3or4];
      double& Ezl=sampled_field[2+DFIELD_3or4];
      const int num_field_components=2*DFIELD_3or4;
      for(int c=0; c<8; c++)
      {
        const double* field_components_c=field_components[c];
        ASSUME_ALIGNED(field_components_c);
        const double weights_c = weights[c];
        #pragma simd
        for(int i=0; i<num_field_components; i++)
        {
          sampled_field[i] += weights_c*field_components_c[i];
        }
      }
      const double Omx = qdto2mc*Bxl;
      const double Omy = qdto2mc*Byl;
      const double Omz = qdto2mc*Bzl;

      // end interpolation
      const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const pfloat denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const pfloat ut = uorig + qdto2mc * Exl;
      const pfloat vt = vorig + qdto2mc * Eyl;
      const pfloat wt = worig + qdto2mc * Ezl;
      //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
      const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
      // solve the velocity equation 
      uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
      vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
      wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
      // update average position
      xavg = xorig + uavg * dto2;
      yavg = yorig + vavg * dto2;
      zavg = zorig + wavg * dto2;
    }                           // end of iteration
    // update the final position and velocity
    //
    if(cap_velocity())
    {
      bool cap = (abs(uavg)>umax || abs(vavg)>vmax || abs(wavg)>wmax)
        ? true : false;
      // we could do something more smooth or sophisticated
      if(builtin_expect(cap,false))
      {
        if(true)
        {
          dprintf("capping velocity: abs(%g,%g,%g)>(%g,%g,%g)",
            uavg,vavg,wavg,umax,vmax,wmax);
        }
        if(uavg>umax) uavg=umax; else if(uavg<umin) uavg=umin;
        if(vavg>vmax) vavg=vmax; else if(vavg<vmin) vavg=vmin;
        if(wavg>wmax) wavg=wmax; else if(wavg<wmin) wavg=wmin;
      }
    }
    //
    pcl->set_x(xorig + uavg * dt);
    pcl->set_y(yorig + vavg * dt);
    pcl->set_z(zorig + wavg * dt);
    pcl->set_u(2.0 * uavg - uorig);
    pcl->set_v(2.0 * vavg - vorig);
    pcl->set_w(2.0 * wavg - worig);
  }                             // END OF ALL THE PARTICLES
}
}


void Particles3D::mover_PC_AoS_Relativistic(Field * EMf)
{
    // The local average number of PC iterations
	// The sum over all processes of the avg. numb. of PC iter
	#ifdef PRINTPCL
	  double innter_sum = 0.0,subcycle_sum=0.0;
	#endif
#pragma omp parallel
{
  convertParticlesToAoS();
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  #pragma omp master
  { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }


  const double PC_err = 1E-6;
  const int    nop=getNOP(); //for OpenMP
  int 		   pidx = 0;      //for OpenMP

#ifdef PRINTPCL
	#pragma omp for schedule(static) reduction(+:innter_sum,subcycle_sum)
#else
	#pragma omp for schedule(static)
#endif
  for (pidx = 0; pidx < nop; pidx++) {

	  // copy the particle
	  SpeciesParticle* pcl = &_pcls[pidx];
	  ALIGNED(pcl);

	  //determine # of subcyles
	  double weights[8] ALLOC_ALIGNED;
	  int cx,cy,cz;
	  grid->get_safe_cell_and_weights(pcl->get_x(),pcl->get_y(),pcl->get_z(),cx,cy,cz,weights);
	  const double* field_components[8] ALLOC_ALIGNED;
	  get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);
	  double Bxl = 0.0;
	  double Byl = 0.0;
	  double Bzl = 0.0;
	  for(int c=0; c<8; c++)
	  {
		Bxl += weights[c] * field_components[c][0];
		Byl += weights[c] * field_components[c][1];
		Bzl += weights[c] * field_components[c][2];
	  }
	  const double B_mag= sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
	  double dt_sub = M_PI*c/(4*abs(qom)*B_mag);
	  const int sub_cycles = int(dt/dt_sub)+1;
	  dt_sub = dt/double(sub_cycles);
	  double dto2_sub = .5 * dt_sub, qdto2mc_sub = qom*dto2_sub/c;

#ifdef PRINTPCL
	  subcycle_sum += sub_cycles;
	  double innter_avg = 0.0;
#endif

	  //start subcycling
	  for(int cyc_cnt=0;cyc_cnt<sub_cycles;cyc_cnt++)
	  {

		    const double xorig = pcl->get_x();
       	    const double yorig = pcl->get_y();
       	    const double zorig = pcl->get_z();
       	    const double uorig = pcl->get_u();
       	    const double vorig = pcl->get_v();
       	    const double worig = pcl->get_w();
       	    double xavg = xorig;
       	    double yavg = yorig;
       	    double zavg = zorig;
       	    double uavg_old, uavg = uorig;
       	    double vavg_old, vavg = vorig;
       	    double wavg_old, wavg = worig;
       	    const double gamma0 = 1.0/(sqrt(1.0 - uorig*uorig - vorig*vorig - worig*worig));
       	    double gamma1;
       	    pfloat ut,vt,wt;

            int innter = 0;
            double currErr = PC_err+1;//initialize to a larger value

            // calculate the average velocity iteratively
            while(currErr> PC_err*PC_err &&  innter < NiterMover){

				  // Save old v_avg
				  uavg_old = uavg;
				  vavg_old = vavg;
				  wavg_old = wavg;

				  // compute weights for field components
				  double weights[8] ALLOC_ALIGNED;
				  int cx,cy,cz;
				  grid->get_safe_cell_and_weights(xavg,yavg,zavg,cx,cy,cz,weights);

				  const double* field_components[8] ALLOC_ALIGNED;
				  get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);

				  double Exl = 0.0;
				  double Eyl = 0.0;
				  double Ezl = 0.0;
				  double Bxl = 0.0;
				  double Byl = 0.0;
				  double Bzl = 0.0;
				  for(int c=0; c<8; c++)
				  {
					Bxl += weights[c] * field_components[c][0];
					Byl += weights[c] * field_components[c][1];
					Bzl += weights[c] * field_components[c][2];
					Exl += weights[c] * field_components[c][0+DFIELD_3or4];
					Eyl += weights[c] * field_components[c][1+DFIELD_3or4];
					Ezl += weights[c] * field_components[c][2+DFIELD_3or4];
				  }
				  const double Omx = qdto2mc_sub*Bxl;
				  const double Omy = qdto2mc_sub*Byl;
				  const double Omz = qdto2mc_sub*Bzl;

				  // end interpolation
				  const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
				  pfloat denom = 1.0 / (1.0 + omsq);


				  //relativistic part

				  // solve the position equation
				  const pfloat ut = uorig*gamma0 + qdto2mc_sub * Exl;
				  const pfloat vt = vorig*gamma0 + qdto2mc_sub * Eyl;
				  const pfloat wt = worig*gamma0 + qdto2mc_sub * Ezl;

				  gamma1 = sqrt(1.0 + ut*ut + vt*vt + wt*wt);
				  Bxl /=gamma1;
				  Byl /=gamma1;
				  Bzl /=gamma1;
				  denom /=gamma1;

				  const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;

				  // solve the velocity equation
				  uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
				  vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
				  wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;

				  // update average position
				  xavg = xorig + uavg * dto2_sub;
				  yavg = yorig + vavg * dto2_sub;
				  zavg = zorig + wavg * dto2_sub;

				  currErr = ((uavg_old-uavg)*(uavg_old-uavg)+(vavg_old-vavg)*(vavg_old-vavg)+(wavg_old-wavg)*(wavg_old-wavg)) /
					     (uavg_old*uavg_old+vavg_old*vavg_old+wavg_old*wavg_old);
				  innter++;
            }  // end of iteration
#ifdef PRINTPCL
            innter_avg += innter;
#endif

            // relativistic velocity update
		    ut = uorig*gamma0;
		    vt = vorig*gamma0;
		    wt = worig*gamma0;
		    double vt_sq = ut*ut + vt*vt + wt*wt;
		    double vavg_sq = uavg*uavg + vavg*vavg + wavg*wavg;
		    double vt_vavg = ut*uavg + vt*vavg + wt*wavg;

		    double cfa = 1.0 - vavg_sq;
		    double cfb =-2.0*(-vt_vavg+gamma0*vavg_sq);
		    double cfc =-1.0-gamma0*gamma0*vavg_sq+2.0*gamma0*vt_vavg - vt_sq;

		    double delta_rel= cfb*cfb -4.0*cfa*cfc;
		    // update velocity
		    if (delta_rel < 0.0){
		        cout << "Relativity violated: gamma0=" << gamma0 << ",  vavg_sq=" << vavg_sq;
		        pcl->set_u((2.0*gamma1)*uavg - uorig*gamma0);
		        pcl->set_v((2.0*gamma1)*vavg - vorig*gamma0);
		        pcl->set_w((2.0*gamma1)*wavg - worig*gamma0);
		    } else {
		        gamma1 = (-cfb+sqrt(delta_rel))/2.0/cfa;
		        pcl->set_u((1.0 + gamma0/gamma1)*uavg - ut/gamma1);
		        pcl->set_v((1.0 + gamma0/gamma1)*vavg - vt/gamma1);
		        pcl->set_w((1.0 + gamma0/gamma1)*wavg - wt/gamma1);
		    }

		    // update the final position and velocity
		    pcl->set_x(xorig + uavg * dt_sub);
		    pcl->set_y(yorig + vavg * dt_sub);
		    pcl->set_z(zorig + wavg * dt_sub);

	  }//END OF subcycling
#ifdef PRINTPCL
	  innter_sum += innter_avg/sub_cycles;
#endif
  }// END OF ALL THE PARTICLES


  #pragma omp master 
  {
	  timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING);

#ifdef PRINTPCL
	  double local_subcycle = subcycle_sum/nop;
	  double local_innter   = innter_sum/nop;
	  double localAvgArr[2];
	  localAvgArr[0]=local_subcycle;
	  localAvgArr[1]=local_innter;
	  double globalAvgArr[2];
	  MPI_Reduce(&localAvgArr, &globalAvgArr, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  if (MPIdata::get_rank() == 0)
		  cout << "*** Relativistic AoS MOVER with Subcycling  species " << ns << " ***" 
		  << globalAvgArr[0]/MPIdata::get_nprocs()  << " subcyles ***" << globalAvgArr[1]/MPIdata::get_nprocs()<< " ITERATIONS   ****" << endl;
#endif
  }
}
}







// move the particle using MIC vector intrinsics
void Particles3D::mover_PC_AoS_vec_intr(Field * EMf)
{
 #ifndef __MIC__
	   if (MPIdata::get_rank() == 0)
			  cout << "*** MIC not enabled, switch to AoS mover " <<endl;
   	   mover_PC_AoS(EMf);
 #else
 #pragma omp parallel
 {
  convertParticlesToAoS();
  // Here and below x stands for all 3 physical position coordinates
  // and u stands for all 3 velocity coordinates.
  const F64vec8 dx_inv = make_F64vec8(get_invdx(), get_invdy(), get_invdz());
  // starting physical position of proper subdomain ("pdom", without ghosts)
  const F64vec8 pdom_xlow = make_F64vec8(get_xstart(),get_ystart(), get_zstart());
  F64vec8 umaximum = make_F64vec8(umax,vmax,wmax);
  F64vec8 uminimum = make_F64vec8(umin,vmin,wmin);
  //
  // compute canonical coordinates of subdomain (including ghosts)
  // relative to global coordinates.
  // x = physical position, X = canonical coordinates.
  //
  // starting position of cell in lower corner
  // of proper subdomain (without ghosts);
  // probably this is an integer value, but we won't rely on it.
  const F64vec8 pdom_Xlow = dx_inv*pdom_xlow;
  // g = including ghosts
  // starting position of cell in low corner
  const F64vec8 gdom_Xlow = pdom_Xlow - F64vec8(1.);
  // starting position of cell in high corner of ghost domain
  // in canonical coordinates
  const F64vec8 nXc = make_F64vec8(nxc,nyc,nzc);
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "***AoS_vec_intr MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  SpeciesParticle * pcls = &_pcls[0];
  ALIGNED(pcls);
  #pragma omp master
  { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
  const double dto2_d = .5 * dt;
  const double qdto2mc_d = qom * dto2_d / c;
  const F64vec8 dto2 = F64vec8(dto2_d);
  const F64vec8 qdto2mc = F64vec8(qdto2mc_d);
  #pragma omp for schedule(static)
  for (int pidx = 0; pidx < getNOP(); pidx+=2)
  {
    // copy the particle
    SpeciesParticle* pcl = &(pcls[pidx]);

    // gather position and velocity data from particles
    //
    F64vec8 pcl0 = *(F64vec8*)&(pcl[0]);
    F64vec8 pcl1 = *(F64vec8*)&(pcl[1]);
    const F64vec8 xorig = cat_hgh_halves(pcl0,pcl1);
    F64vec8 xavg = xorig;
    const F64vec8 uorig = cat_low_halves(pcl0,pcl1);

    // calculate the average velocity iteratively
    //
    // (could stop iteration when it is determined that
    // both particles are converged, e.g. if change in
    // xavg is sufficiently small)
    F64vec8 uavg;
    for (int iter = 0; iter < NiterMover; iter++)
    {
      // convert to canonical coordinates relative to subdomain with ghosts
      const F64vec8 gX = dx_inv*xavg - gdom_Xlow;
      F64vec8 cellXstart = floor(gX);
      // map to cell within the process subdomain (including ghosts);
      // this is triggered if xavg is outside the ghost subdomain
      // and results in extrapolation from the nearest ghost cell
      // rather than interpolation as in the usual case.
      cellXstart = maximum(cellXstart,F64vec8(0.));
      cellXstart = minimum(cellXstart,nXc);
      // get cell coordinates.
      const I32vec16 cell = round_to_nearest(cellXstart);
      // get field_components for each particle
      F64vec8 field_components0[8]; // first pcl
      F64vec8 field_components1[8]; // second pcl
      ::get_field_components_for_cell(
        field_components0,field_components1,fieldForPcls,cell);

      // get weights for field_components based on particle position
      //
      F64vec8 weights[2];
      const F64vec8 X = gX - cellXstart;
      construct_weights_for_2pcls(weights, X);

      // interpolate field to get fields
      F64vec8 fields[2];
      // sample fields for first particle
      fields[0] = sample_field_mic(weights[0],field_components0);
      // sample fields for second particle
      fields[1] = sample_field_mic(weights[1],field_components1);
      const F64vec8 B = cat_low_halves(fields[0],fields[1]);
      const F64vec8 E = cat_hgh_halves(fields[0],fields[1]);

      // use sampled field to push particle block
      //
      uavg = compute_uvg_for_2pcls(uorig, B, E, qdto2mc);
      // update average position
      xavg = xorig + uavg*dto2;
    } // end of iterative particle advance
    // update the final position and velocity
    if(cap_velocity())
    {
      uavg = minimum(uavg,umaximum);
      uavg = maximum(uavg,uminimum);
    }
    const F64vec8 xnew = xavg+(xavg - xorig);
    const F64vec8 unew = uavg+(uavg - uorig);
    const F64vec8 pcl0new = cat_low_halves(unew, xnew);
    const F64vec8 pcl1new = cat_hgh_halves(unew, xnew);
    copy012and456(pcl0,pcl0new);
    copy012and456(pcl1,pcl1new);

    // could save using no-read stores( _mm512_storenr_pd),
    // but we just read this, so presumably it is still in cache.
    _mm512_store_pd(&pcl[0], pcl0);
    _mm512_store_pd(&pcl[1], pcl1);
  }
  #pragma omp master
  { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
 }
 #endif
}

void Particles3D::mover_PC_AoS_vec(Field * EMf)
{
#pragma omp parallel
{
  convertParticlesToAoS();
  #pragma omp master
  if (vct->getCartesian_rank() == 0) {
    cout << "***AoS_vec MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();

  const int NUM_PCLS_MOVED_AT_A_TIME = 8;
  // make sure that we won't overrun memory
  int needed_capacity = roundup_to_multiple(getNOP(),NUM_PCLS_MOVED_AT_A_TIME);
  assert_le(needed_capacity,_pcls.capacity());

  #pragma omp master
  { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
  const double dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
  #pragma omp for schedule(static)
  for (int pidx = 0; pidx < getNOP(); pidx+=NUM_PCLS_MOVED_AT_A_TIME)
  {
    // copy the particles
    SpeciesParticle* pcl[NUM_PCLS_MOVED_AT_A_TIME];
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    {
      pcl[i] = &_pcls[pidx+i];
    }
    // actually, all the particles are aligned,
    // but the compiler should be able to see that.
    ALIGNED(pcl[0]);
    double xorig[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
    double uorig[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
    double  xavg[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
    double  uavg[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
    // gather data into vectors
    // #pragma simd collapse(2)
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    for(int j=0;j<3;j++)
    {
      xavg[i][j] = xorig[i][j] = pcl[i]->get_x(j);
      uorig[i][j] = pcl[i]->get_u(j);
    }
    // calculate the average velocity iteratively
    for (int innter = 0; innter < NiterMover; innter++) {

      // compute weights for field components
      //
      double weights[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(64)));
      int cx[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        grid->get_safe_cell_and_weights(xavg[i],cx[i],weights[i]);
      }

      const double* field_components[NUM_PCLS_MOVED_AT_A_TIME][8] __attribute__((aligned(64)));
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        get_field_components_for_cell(field_components[i],fieldForPcls,
          cx[i][0],cx[i][1],cx[i][2]);
      }

      double E[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
      double B[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
      // could do this with memset
      // #pragma simd collapse(2)
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        E[i][j]=0;
        B[i][j]=0;
      }
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      for(int c=0; c<8; c++)
      {
        B[i][j] += weights[i][c] * field_components[i][c][j];
        E[i][j] += weights[i][c] * field_components[i][c][j+DFIELD_3or4];
      }
      double Om[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        Om[i][j] = qdto2mc*B[i][j];
      }

      // can these dot products vectorize if
      // NUM_PCLS_MOVED_AT_A_TIME is large enough?
      double omsq[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(64)));
      double denom[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(64)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        omsq[i] = Om[i][0] * Om[i][0]
                + Om[i][1] * Om[i][1]
                + Om[i][2] * Om[i][2];
        denom[i] = 1.0 / (1.0 + omsq[i]);
      }
      // solve the position equation
      double ut[NUM_PCLS_MOVED_AT_A_TIME][3] __attribute__((aligned(64)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        ut[i][j] = uorig[i][j] + qdto2mc * E[i][j];
      }
      double udotOm[NUM_PCLS_MOVED_AT_A_TIME] __attribute__((aligned(64)));
      for(int i=0; i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        udotOm[i] = ut[i][0] * Om[i][0]
                  + ut[i][1] * Om[i][1]
                  + ut[i][2] * Om[i][2];
      }
      // solve the velocity equation 
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      {
        // these cross-products might not vectorize so well...
        uavg[i][0] = (ut[i][0] + (ut[i][1] * Om[i][2] - ut[i][2] * Om[i][1] + udotOm[i] * Om[i][0])) * denom[i];
        uavg[i][1] = (ut[i][1] + (ut[i][2] * Om[i][0] - ut[i][0] * Om[i][2] + udotOm[i] * Om[i][1])) * denom[i];
        uavg[i][2] = (ut[i][2] + (ut[i][0] * Om[i][1] - ut[i][1] * Om[i][0] + udotOm[i] * Om[i][2])) * denom[i];
      }
      // update average position
      // #pragma simd collapse(2)
      for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
      for(int j=0;j<3;j++)
      {
        xavg[i][j] = xorig[i][j] + uavg[i][j] * dto2;
      }
    } // end of iteration
    // update the final position and velocity (scatter)
    for(int i=0;i<NUM_PCLS_MOVED_AT_A_TIME;i++)
    for(int j=0;j<3;j++)
    {
      pcl[i]->set_x(j, xorig[i][j] + uavg[i][j] * dt);
      pcl[i]->set_u(j, 2.*uavg[i][j] - uorig[i][j]);
    }
  }
  #pragma omp master
  { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
}
}

// This currently computes extrapolated values based on field in
// original mesh cell (unstable?), but execution time suggests
// bound on performance.  For correct execution would need to
// sort by xavg with each iteration like in mover_PC_vectorized.
// But in fact this does not run any faster than mover_PC_AoS
//
//void Particles3D::mover_PC_AoS_vec_onesort(Field * EMf)
//{
//  convertParticlesToAoS();
//  #pragma omp master
//  if (vct->getCartesian_rank() == 0) {
//    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
//  }
//  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();
//
//  SpeciesParticle * pcls = fetch_pcls();
//  #pragma omp master
//  { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
//  const double dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
//
//  #pragma omp for collapse(2) // schedule(static)
//  for(int cx=0;cx<nxc;cx++)
//  for(int cy=0;cy<nyc;cy++)
//  for(int cz=0;cz<nzc;cz++)
//  //for(int cell=0; cell<ncells; cell++)
//  {
//    // Idea of this function is that we only need
//    // to do this once for each group of particles.
//    //
//    const double* field_components[8];
//    get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);
//
//    // push all particles in mesh cell
//    //
//    //const int numpcls_in_cell = numpcls_in_bucket_1d[cell];
//    const int numpcls_in_cell = get_numpcls_in_bucket(cx,cy,cz);
//    const int bucket_offset = get_bucket_offset(cx,cy,cz);
//    const int bucket_end = bucket_offset+numpcls_in_cell;
//    for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
//    {
//      SpeciesParticle* pcl = &pcls[pidx];
//      ALIGNED(pcl);
//      // copy the particle
//      const pfloat xorig = pcl->get_x();
//      const pfloat yorig = pcl->get_y();
//      const pfloat zorig = pcl->get_z();
//      const pfloat uorig = pcl->get_u();
//      const pfloat vorig = pcl->get_v();
//      const pfloat worig = pcl->get_w();
//      double xavg = xorig;
//      double yavg = yorig;
//      double zavg = zorig;
//      double uavg;
//      double vavg;
//      double wavg;
//      // calculate the average velocity iteratively
//      for (int innter = 0; innter < NiterMover; innter++) {
//
//        // compute weights for field components
//        //
//        double weights[8];
//        // xstart marks start of domain excluding ghosts
//        const double rel_xpos = xavg - xstart;
//        const double rel_ypos = yavg - ystart;
//        const double rel_zpos = zavg - zstart;
//        // cell position minus 1 (due to ghost cells)
//        const double cxm1_pos = rel_xpos * inv_dx;
//        const double cym1_pos = rel_ypos * inv_dy;
//        const double czm1_pos = rel_zpos * inv_dz;
//
//        // fraction of the distance from the right of the cell
//        const double w1x = cx - cxm1_pos;
//        const double w1y = cy - cym1_pos;
//        const double w1z = cz - czm1_pos;
//        // fraction of distance from the left
//        const double w0x = 1-w1x;
//        const double w0y = 1-w1y;
//        const double w0z = 1-w1z;
//        //
//        Grid::get_weights(weights, w0x, w0y, w0z, w1x, w1y, w1z);
//
//        //if(false) // this would fail
//        //{
//        //   int cx_,cy_,cz_;
//        //   grid->get_safe_cell_coordinates(xavg,yavg,zavg,cx_,cy_,cz_);
//        //   assert_eq(cx,cx_);
//        //   assert_eq(cy,cy_);
//        //   assert_eq(cz,cz_);
//        //}
//
//        pfloat Exl = 0.0;
//        pfloat Eyl = 0.0;
//        pfloat Ezl = 0.0;
//        pfloat Bxl = 0.0;
//        pfloat Byl = 0.0;
//        pfloat Bzl = 0.0;
//        for(int c=0; c<8; c++)
//        {
//          Bxl += weights[c] * field_components[c][0];
//          Byl += weights[c] * field_components[c][1];
//          Bzl += weights[c] * field_components[c][2];
//          Exl += weights[c] * field_components[c][0+DFIELD_3or4];
//          Eyl += weights[c] * field_components[c][1+DFIELD_3or4];
//          Ezl += weights[c] * field_components[c][2+DFIELD_3or4];
//        }
//        const double Omx = qdto2mc*Bxl;
//        const double Omy = qdto2mc*Byl;
//        const double Omz = qdto2mc*Bzl;
//
//        // end interpolation
//        const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
//        const pfloat denom = 1.0 / (1.0 + omsq);
//        // solve the position equation
//        const pfloat ut = uorig + qdto2mc * Exl;
//        const pfloat vt = vorig + qdto2mc * Eyl;
//        const pfloat wt = worig + qdto2mc * Ezl;
//        //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
//        const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
//        // solve the velocity equation 
//        uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
//        vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
//        wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
//        // update average position
//        xavg = xorig + uavg * dto2;
//        yavg = yorig + vavg * dto2;
//        zavg = zorig + wavg * dto2;
//      }
//      // update the final position and velocity
//      pcl->set_x(xorig + uavg * dt);
//      pcl->set_y(yorig + vavg * dt);
//      pcl->set_z(zorig + wavg * dt);
//      pcl->set_u(2.0 * uavg - uorig);
//      pcl->set_v(2.0 * vavg - vorig);
//      pcl->set_w(2.0 * wavg - worig);
//    }
//  }
//  #pragma omp master
//  { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
//}

///** mover with a Predictor-Corrector scheme */
//void Particles3D::mover_PC_vectorized(Field * EMf)
//{
//  convertParticlesToSoA();
//  assert_eq(nxc,nxn-1);
//  assert_eq(nyc,nyn-1);
//  assert_eq(nzc,nzn-1);
//  #pragma omp master
//  if (vct->getCartesian_rank() == 0) {
//    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
//  }
//  const_arr4_pfloat fieldForPcls = EMf->get_fieldForPcls();
//
//  // initialize average positions
//  #pragma omp for schedule(static)
//  for(int pidx = 0; pidx < getNOP(); pidx++)
//  {
//    _xavg[pidx] = x[pidx];
//    _yavg[pidx] = y[pidx];
//    _zavg[pidx] = z[pidx];
//  }
//
//  const pfloat dto2 = .5 * dt, qdto2mc = qom * dto2 / c;
//  for(int niter=1; niter<=NiterMover; niter++)
//  {
//    // sort particles based on the time-averaged position
//    if(niter>1) // on first iteration already was sorted to sum moments
//    {
//      #pragma omp master
//      {
//        timeTasks_begin_task(TimeTasks::MOVER_PCL_SORTING);
//        // this changes the definitions of x,y,z,u,v,w,_xavg,_yavg,_zavg,etc.
//        sort_particles_serial_SoA_by_xavg(grid,vct);
//        timeTasks_end_task(TimeTasks::MOVER_PCL_SORTING);
//      }
//      #pragma omp barrier
//    }
//
//    #pragma omp master
//    { timeTasks_begin_task(TimeTasks::MOVER_PCL_MOVING); }
//    // move particles in parallel
//    //
//    // iterate over mesh cells
//    //const int ncells=nxc*nyc*nzc;
//    //int *numpcls_in_bucket_1d = &numpcls_in_bucket[0][0][0];
//    //int *bucket_offset_1d = &bucket_offset[0][0][0];
//    ALIGNED(x);
//    ALIGNED(y);
//    ALIGNED(z);
//    ALIGNED(u);
//    ALIGNED(v);
//    ALIGNED(w);
//    ALIGNED(_xavg);
//    ALIGNED(_yavg);
//    ALIGNED(_zavg);
//    int serial_pidx = 0;
//    #pragma omp for collapse(2) // schedule(static)
//    for(int cx=0;cx<nxc;cx++)
//    for(int cy=0;cy<nyc;cy++)
//    for(int cz=0;cz<nzc;cz++)
//    //for(int cell=0; cell<ncells; cell++)
//    {
//      // interface to the right of cell
//      const int ix = cx+1;
//      const int iy = cy+1;
//      const int iz = cz+1;
//
//      const double* field_components[8];
//      field_components[0] = fieldForPcls[ix][iy][iz]; // field000
//      field_components[1] = fieldForPcls[ix][iy][cz]; // field001
//      field_components[2] = fieldForPcls[ix][cy][iz]; // field010
//      field_components[3] = fieldForPcls[ix][cy][cz]; // field011
//      field_components[4] = fieldForPcls[cx][iy][iz]; // field100
//      field_components[5] = fieldForPcls[cx][iy][cz]; // field101
//      field_components[6] = fieldForPcls[cx][cy][iz]; // field110
//      field_components[7] = fieldForPcls[cx][cy][cz]; // field111
//
//      // push all particles in mesh cell
//      //
//      //const int numpcls_in_cell = numpcls_in_bucket_1d[cell];
//      const int numpcls_in_cell = get_numpcls_in_bucket(cx,cy,cz);
//      const int bucket_offset = get_bucket_offset(cx,cy,cz);
//      const int bucket_end = bucket_offset+numpcls_in_cell;
//      // This pragma helps on Xeon but hurts on Xeon Phi.
//      // On the Phi we could accelerate by processing two particles at a time.
//      // there should be no function calls in this loop (except inlined calls)
//      #pragma simd
//      for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
//      {
//        // serial case: check that pidx is correct
//        //assert_eq(pidx,serial_pidx);
//        //serial_pidx++;
//        // confirm that particle is in correct cell
//        //if(true)
//        //{
//        //  int cx_,cy_,cz_;
//        //  get_safe_cell_for_pos(cx_,cy_,cz_,_xavg[pidx],_yavg[pidx],_zavg[pidx]);
//        //  if((cx_!=cx)
//        //   ||(cy_!=cy)
//        //   ||(cz_!=cz))
//        //  {
//        //    dprintf("\n\t cx =%d, cy =%d, cz =%d"
//        //            "\n\t cx_=%d, cy_=%d, cz_=%d"
//        //            "\n\t cxf=%g, cyf=%g, czf=%g",
//        //            cx,cy,cz,
//        //            cx_,cy_,cz_,
//        //            1+(_xavg[pidx]-xstart)*inv_dx,
//        //            1+(_yavg[pidx]-ystart)*inv_dy,
//        //            1+(_zavg[pidx]-zstart)*inv_dz);
//        //  }
//        //  assert_eq(cx_,cx);
//        //  assert_eq(cy_,cy);
//        //  assert_eq(cz_,cz);
//        //}
//
//        // copy the particle
//        const pfloat xorig = x[pidx];
//        const pfloat yorig = y[pidx];
//        const pfloat zorig = z[pidx];
//        const pfloat uorig = u[pidx];
//        const pfloat vorig = v[pidx];
//        const pfloat worig = w[pidx];
//
//        // compute weights for field components
//        //
//        double weights[8];
//        const double abs_xpos = _xavg[pidx];
//        const double abs_ypos = _yavg[pidx];
//        const double abs_zpos = _zavg[pidx];
//        // xstart marks start of domain excluding ghosts
//        const double rel_xpos = abs_xpos - xstart;
//        const double rel_ypos = abs_ypos - ystart;
//        const double rel_zpos = abs_zpos - zstart;
//        // cell position minus 1 (due to ghost cells)
//        const double cxm1_pos = rel_xpos * inv_dx;
//        const double cym1_pos = rel_ypos * inv_dy;
//        const double czm1_pos = rel_zpos * inv_dz;
//        // index of interface to right of cell
//        const int ix = cx + 1;
//        const int iy = cy + 1;
//        const int iz = cz + 1;
//        // fraction of the distance from the right of the cell
//        const double w1x = cx - cxm1_pos;
//        const double w1y = cy - cym1_pos;
//        const double w1z = cz - czm1_pos;
//        // fraction of distance from the left
//        const double w0x = 1-w1x;
//        const double w0y = 1-w1y;
//        const double w0z = 1-w1z;
//        const double weight00 = w0x*w0y;
//        const double weight01 = w0x*w1y;
//        const double weight10 = w1x*w0y;
//        const double weight11 = w1x*w1y;
//        weights[0] = weight00*w0z; // weight000
//        weights[1] = weight00*w1z; // weight001
//        weights[2] = weight01*w0z; // weight010
//        weights[3] = weight01*w1z; // weight011
//        weights[4] = weight10*w0z; // weight100
//        weights[5] = weight10*w1z; // weight101
//        weights[6] = weight11*w0z; // weight110
//        weights[7] = weight11*w1z; // weight111
//
//        pfloat Exl = 0.0;
//        pfloat Eyl = 0.0;
//        pfloat Ezl = 0.0;
//        pfloat Bxl = 0.0;
//        pfloat Byl = 0.0;
//        pfloat Bzl = 0.0;
//
//        // would expanding this out help to vectorize?
//        for(int c=0; c<8; c++)
//        {
//          Bxl += weights[c] * field_components[c][0];
//          Byl += weights[c] * field_components[c][1];
//          Bzl += weights[c] * field_components[c][2];
//          Exl += weights[c] * field_components[c][0+DFIELD_3or4];
//          Eyl += weights[c] * field_components[c][1+DFIELD_3or4];
//          Ezl += weights[c] * field_components[c][2+DFIELD_3or4];
//        }
//        const double Omx = qdto2mc*Bxl;
//        const double Omy = qdto2mc*Byl;
//        const double Omz = qdto2mc*Bzl;
//
//        // end interpolation
//        const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
//        const pfloat denom = 1.0 / (1.0 + omsq);
//        // solve the position equation
//        const pfloat ut = uorig + qdto2mc * Exl;
//        const pfloat vt = vorig + qdto2mc * Eyl;
//        const pfloat wt = worig + qdto2mc * Ezl;
//        //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
//        const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
//        // solve the velocity equation 
//        const pfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
//        const pfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
//        const pfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
//        // update average position
//        _xavg[pidx] = xorig + uavg * dto2;
//        _yavg[pidx] = yorig + vavg * dto2;
//        _zavg[pidx] = zorig + wavg * dto2;
//
//        // if it is the last iteration, update the position and velocity
//        // (hopefully this will not compromise vectorization...)
//        if(niter==NiterMover)
//        {
//          x[pidx] = xorig + uavg * dt;
//          y[pidx] = yorig + vavg * dt;
//          z[pidx] = zorig + wavg * dt;
//          u[pidx] = 2.0 * uavg - uorig;
//          v[pidx] = 2.0 * vavg - vorig;
//          w[pidx] = 2.0 * wavg - worig;
//        }
//      }
//    }
//    #pragma omp master
//    { timeTasks_end_task(TimeTasks::MOVER_PCL_MOVING); }
//  }
//}

/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_relativistic(Field * EMf)
{
  eprintf("Mover_relativistic not implemented");
  return (0);
}

inline void Particles3D::populate_cell_with_particles(
  int i, int j, int k, double q_per_particle,
  double dx_per_pcl, double dy_per_pcl, double dz_per_pcl)
{
  const double cell_low_x = grid->getXN(i,j,k);
  const double cell_low_y = grid->getYN(i,j,k);
  const double cell_low_z = grid->getZN(i,j,k);
  for (int ii=0; ii < npcelx; ii++)
  for (int jj=0; jj < npcely; jj++)
  for (int kk=0; kk < npcelz; kk++)
  {
    double u,v,w,q,x,y,z;
    sample_maxwellian(
      u,v,w,
      uth, vth, wth,
      u0, v0, w0);
    x = (ii + sample_u_double())*dx_per_pcl + cell_low_x;
    y = (jj + sample_u_double())*dy_per_pcl + cell_low_y;
    z = (kk + sample_u_double())*dz_per_pcl + cell_low_z;
    create_new_particle(u,v,w,q_per_particle,x,y,z);
  }
}

// This could be generalized to use fluid moments
// to generate particles.
//
void Particles3D::repopulate_particles()
{
  using namespace BCparticles;

  // if this is not a boundary process then there is nothing to do
  if(!vct->isBoundaryProcess_P()) return;

  // if there are no reemission boundaries then no one has anything to do
  const bool repop_bndry_in_X = !vct->getPERIODICX_P() &&
        (bcPfaceXleft == REEMISSION || bcPfaceXright == REEMISSION);
  const bool repop_bndry_in_Y = !vct->getPERIODICY_P() &&
        (bcPfaceYleft == REEMISSION || bcPfaceYright == REEMISSION);
  const bool repop_bndry_in_Z = !vct->getPERIODICZ_P() &&
        (bcPfaceZleft == REEMISSION || bcPfaceZright == REEMISSION);
  const bool repopulation_boundary_exists =
        repop_bndry_in_X || repop_bndry_in_Y || repop_bndry_in_Z;

  if(!repopulation_boundary_exists) return;


  // boundaries to repopulate
  //
  const bool repopulateXleft = (vct->noXleftNeighbor_P() && bcPfaceXleft == REEMISSION);
  const bool repopulateYleft = (vct->noYleftNeighbor_P() && bcPfaceYleft == REEMISSION);
  const bool repopulateZleft = (vct->noZleftNeighbor_P() && bcPfaceZleft == REEMISSION);
  const bool repopulateXrght = (vct->noXrghtNeighbor_P() && bcPfaceXright == REEMISSION);
  const bool repopulateYrght = (vct->noYrghtNeighbor_P() && bcPfaceYright == REEMISSION);
  const bool repopulateZrght = (vct->noZrghtNeighbor_P() && bcPfaceZright == REEMISSION);
  const bool do_repopulate = 
       repopulateXleft || repopulateYleft || repopulateZleft
    || repopulateXrght || repopulateYrght || repopulateZrght;
  // if this process has no reemission boundaries then there is nothing to do
  if(!do_repopulate)
    return;

  // there are better ways to obtain these values...
  //
  double  FourPI =16*atan(1.0);
  const double q_per_particle
    = (qom/fabs(qom))*(Ninj/FourPI/npcel)*(1.0/grid->getInvVOL());

  const int nxc = grid->getNXC();
  const int nyc = grid->getNYC(); const int nzc = grid->getNZC();
  // number of cell layers to repopulate at boundary
  const int num_layers = 3;
  const double xLow = num_layers*dx;
  const double yLow = num_layers*dy;
  const double zLow = num_layers*dz;
  const double xHgh = Lx-xLow;
  const double yHgh = Ly-yLow;
  const double zHgh = Lz-zLow;
  if(repopulateXleft || repopulateXrght) assert_gt(nxc, 2*num_layers);
  if(repopulateYleft || repopulateYrght) assert_gt(nyc, 2*num_layers);
  if(repopulateZleft || repopulateZrght) assert_gt(nzc, 2*num_layers);

  // delete particles in repopulation layers
  //
  const int nop_orig = getNOP();
  int pidx = 0;
  while(pidx < getNOP())
  {
    SpeciesParticle& pcl = _pcls[pidx];
    // determine whether to delete the particle
    const bool delete_pcl =
      (repopulateXleft && pcl.get_x() < xLow) ||
      (repopulateYleft && pcl.get_y() < yLow) ||
      (repopulateZleft && pcl.get_z() < zLow) ||
      (repopulateXrght && pcl.get_x() > xHgh) ||
      (repopulateYrght && pcl.get_y() > yHgh) ||
      (repopulateZrght && pcl.get_z() > zHgh);
    if(delete_pcl)
      delete_particle(pidx);
    else
      pidx++;
  }
  const int nop_remaining = getNOP();

  const double dx_per_pcl = dx/npcelx;
  const double dy_per_pcl = dy/npcely;
  const double dz_per_pcl = dz/npcelz;

  // starting coordinate of upper layer
  const int upXstart = nxc-1-num_layers;
  const int upYstart = nyc-1-num_layers;
  const int upZstart = nzc-1-num_layers;

  // inject new particles.
  //
  {
    // we shrink the imagined boundaries of the array as we go along to ensure
    // that we never inject particles twice in a single mesh cell.
    //
    // initialize imagined boundaries to full subdomain excluding ghost cells.
    //
    int xbeg = 1;
    int xend = nxc-2;
    int ybeg = 1;
    int yend = nyc-2;
    int zbeg = 1;
    int zend = nzc-2;
    if (repopulateXleft)
    {
      //cout << "*** Repopulate Xleft species " << ns << " ***" << endl;
      for (int i=1; i<= num_layers; i++)
      for (int j=ybeg; j<=yend; j++)
      for (int k=zbeg; k<=zend; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
      // these have all been filled, so never touch them again.
      xbeg += num_layers;
    }
    if (repopulateXrght)
    {      
      cout << "*** Repopulate Xright species " << ns << " ***" << endl;
      for (int i=upXstart; i<=xend; i++)
      for (int j=ybeg; j<=yend; j++)
      for (int k=zbeg; k<=zend; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
      // these have all been filled, so never touch them again.
      xend -= num_layers;
    }
    if (repopulateYleft)
    {     
      // cout << "*** Repopulate Yleft species " << ns << " ***" << endl;
      for (int i=xbeg; i<=xend; i++)
      for (int j=1; j<=num_layers; j++)
      for (int k=zbeg; k<=zend; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
      // these have all been filled, so never touch them again.
      ybeg += num_layers;
    }
    if (repopulateYrght)
    {     
      // cout << "*** Repopulate Yright species " << ns << " ***" << endl;
      for (int i=xbeg; i<=xend; i++)
      for (int j=upYstart; j<=yend; j++)
      for (int k=zbeg; k<=zend; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
      // these have all been filled, so never touch them again.
      yend -= num_layers;
    }
    if (repopulateZleft)
    {   
      //   cout << "*** Repopulate Zleft species " << ns << " ***" << endl;
      for (int i=xbeg; i<=xend; i++)
      for (int j=ybeg; j<=yend; j++)
      for (int k=1; k<=num_layers; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
    }
    if (repopulateZrght)
    {   
      //   cout << "*** Repopulate Zright species " << ns << " ***" << endl;
      for (int i=xbeg; i<=xend; i++)
      for (int j=ybeg; j<=yend; j++)
      for (int k=upZstart; k<=zend; k++)
      {
        populate_cell_with_particles(i,j,k,q_per_particle,
          dx_per_pcl, dy_per_pcl, dz_per_pcl);
      }
    }
  }
  const int nop_final = getNOP();
  const int nop_deleted = nop_orig - nop_remaining;
  const int nop_created = nop_final - nop_remaining;

  //dprintf("change in # particles: %d - %d + %d = %d",nop_orig, nop_deleted, nop_created, nop_final);

  //if (vct->getCartesian_rank()==0){
  //  cout << "*** number of particles " << getNOP() << " ***" << endl;
  //}
}

// Open BC for particles: duplicate particles on the boundary,.
// shift outside the box and update location to test if inside box
// if so, add to particle list
void Particles3D::openbc_particles_inflow()
{
  eprintf("openbc_particles_inflow has not yet been implemented");

/*
  if(!vct->isBoundaryProcess()) return;

  using namespace BCparticles;
  //This is a hack here to test injection on XLeft
  if(!vct->getPERIODICX() && vct->noXleftNeighbor() && bcPfaceXleft == OPENBCIn){

	dprintf( "******  Inflow OpenBC ******");
    int delpcl = 0;
    double delq=0.0;

    //delete those exiting particles
    int pidx = 0;
    double pclX, pclY, pclZ;
    while(pidx < getNOP()){
      SpeciesParticle& pcl = _pcls[pidx];
      pclX= pcl.get_x();
      pclY= pcl.get_y();
      pclZ= pcl.get_z();
      if( pclX < 0 || pclX > Lx  || pclY < 0 || pclY > Ly || pclZ < 0 || pclZ > Lz){
	delpcl ++;
	delq += pcl.get_q();
	delete_particle(pidx);
      }else
	pidx ++;
    }

    dprintf("delete %d pcl, %f charges, now NOP=%d",delpcl,delq,getNOP());
    const int newpartStartID = getNOP();

    //Create Maxwellian Uniformly distributed particles
    srand(vct->getCartesian_rank() + time(0) + ns);dprintf("Seed =%d",vct->getCartesian_rank() + time(0) + ns);
    const double FourPI = 16*atan(1.0);
    const double q_sgn  = (qom / fabs(qom));
    const double InjCellx = Vinj * dt; dprintf("Vinj=%f, dt=%f, InjCellx=%f",Vinj,dt,InjCellx);
    //const double InjVol   = InjCellx*(dy*(grid->getNYC()-2))*(dz*(grid->getNZC()-2));
    const double InjVol   = grid->getVOL()*InjCellx/dx;
    const double q_tot_particle =  q_sgn*col->getRHOinit(ns)/FourPI*InjVol;
    //const double q_tot_particle =  delq;
    const double q_per_particle =  q_sgn*col->getRHOinit(ns)/FourPI*InjVol/npcel;
    dprintf("q_tot_particle=%f, q_per_particle=%f",q_tot_particle,q_per_particle);

    int nop_inject=0;
    double q_inject = 0.0;
    const double InjCellStartX = 0.0-InjCellx;
    const double invRandMax = 1.0/double(RAND_MAX);


    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++){

    	/*
    	//The below is uniformly distributed in location
		for (int ii=0; ii < npcelx; ii++)
		  for (int jj=0; jj < npcely; jj++)
			for (int kk=0; kk < npcelz; kk++){

			  double u,v,w,x,y,z;
			  sample_maxwellian(u,v,w,uth, vth, wth,u0, v0, w0);
			  //dprintf("sample_maxwellian %f, %f, %f, %f, %f, %f,%f,%f,%f)",u,v,w,uth, vth, wth,u0, v0, w0);

			  //Uniform distribution
			  x = (ii + .5) * (InjCellx / npcelx) + InjCellStartX;
			  y = (jj + .5) * (dy / npcely) + grid->getYN(1, j, k);
			  z = (kk + .5) * (dz / npcelz) + grid->getZN(1, j, k);

			  //check location after one time step
			  x += u*dt;
			  y += v*dt;
			  z += w*dt;

			  //Add particle if it enters the domain box
			  //here may need communicating those exiting particles
			  if( x>0 && x<Lx && y>0 && y<Ly && z>0 && z<Lz){
				create_new_particle(u,v,w,q_per_particle,x,y,z);
				nop_inject ++;
				q_inject += q_per_particle;
			  }
		}

		//The below is randomly distributed in location
		for (int pclid=0; pclid < npcel; pclid++){

		  double u,v,w,x,y,z;
		  sample_maxwellian(u,v,w,uth, vth, wth,u0, v0, w0);

		  x = InjCellStartX		   + (rand()*invRandMax)*InjCellx;
		  y = grid->getYN(1, j, k) + (rand()*invRandMax)*dy;
		  z = grid->getZN(1, j, k) + (rand()*invRandMax)*dz;

		  //check location after one time step
		  x += u*dt;
		  y += v*dt;
		  z += w*dt;

		  //Add particle if it enters the domain box
		  //here may need communicating those exiting particles
		  if( x>0 && x<Lx && y>0 && y<Ly && z>0 && z<Lz){
			create_new_particle(u,v,w,q_tot_particle,x,y,z);
			nop_inject ++;
			q_inject += q_per_particle;
		  }
		}

     }
      
     //the below part is for dividing over the number of entering particles
      const int newnop = getNOP();
      q_inject=0.0;
      for(int startid=newpartStartID;startid<newnop;startid++){
          SpeciesParticle& pcl = _pcls[startid];
          pcl.fetch_q() /= nop_inject;
          q_inject += pcl.get_q();
	  	  if(startid==newpartStartID) dprintf("x, y,z, u, v, w, q=%f,%f,%f,%f,%f,%f,%f",pcl.get_x(),pcl.get_y(),pcl.get_z(),pcl.get_u(),pcl.get_v(),pcl.get_w(),pcl.get_q());
      }

      dprintf("create %d species %d  pcl,  now NOP= %d , inject Q=%f" ,nop_inject, ns, getNOP(),q_inject);
   
  }//end of injecting from Xleft
  */
}

// Open BC for particles: duplicate particles on the boundary,.
// shift outside the box and update location to test if inside box
// if so, add to particle list
void Particles3D::openbc_particles_outflow()
{
  // if this is not a boundary process then there is nothing to do
  if(!vct->isBoundaryProcess_P()) return;

  //The below is OpenBC outflow for all other boundaries
  using namespace BCparticles;

  const bool openXleft = !vct->getPERIODICX_P() && vct->noXleftNeighbor_P() &&  bcPfaceXleft == OPENBCOut;
  const bool openYleft = !vct->getPERIODICY_P() && vct->noYleftNeighbor_P() &&  bcPfaceYleft == OPENBCOut;
  const bool openZleft = !vct->getPERIODICZ_P() && vct->noZleftNeighbor_P() &&  bcPfaceZleft == OPENBCOut;

  const bool openXright = !vct->getPERIODICX_P() && vct->noXrghtNeighbor_P() && bcPfaceXright == OPENBCOut;
  const bool openYright = !vct->getPERIODICY_P() && vct->noYrghtNeighbor_P() && bcPfaceYright == OPENBCOut;
  const bool openZright = !vct->getPERIODICZ_P() && vct->noZrghtNeighbor_P() && bcPfaceZright == OPENBCOut;

  if(!openXleft && !openYleft && !openZleft && !openXright && !openYright && !openZright)  return;

  const int num_layers = 3;
  assert_gt(nxc-2, (openXleft+openXright)*num_layers); //excluding 2 ghost cells, #of cells should be larger than total # of openBC layers
  assert_gt(nyc-2, (openYleft+openYright)*num_layers);
  assert_gt(nzc-2, (openZleft+openZright)*num_layers);

  const double xLow = num_layers*dx;
  const double yLow = num_layers*dy;
  const double zLow = num_layers*dz;
  const double xHgh = Lx-xLow;
  const double yHgh = Ly-yLow;
  const double zHgh = Lz-zLow;

  const bool   apply_openBC[6]    = {openXleft, openXright,openYleft, openYright,openZleft, openZright};
  const double delete_boundary[6] = {0, Lx,0, Ly,0, Lz};
  const double open_boundary[6]   = {xLow, xHgh,yLow, yHgh,zLow, zHgh};

  const int nop_orig = getNOP();
  const int capacity_out = roundup_to_multiple(nop_orig*0.1,DVECWIDTH);
  vector_SpeciesParticle injpcls(capacity_out);


  for(int dir_cnt=0;dir_cnt<6;dir_cnt++){

    if(apply_openBC[dir_cnt]){
    	  //dprintf( "*** OpenBC for Direction %d on particle species %d",dir_cnt, ns);

		  int pidx = 0;
		  int direction = dir_cnt/2;
		  double delbry  = delete_boundary[dir_cnt];
		  double openbry = open_boundary[dir_cnt];
		  double location;
		  while(pidx < getNOP())
		  {
		     SpeciesParticle& pcl = _pcls[pidx];
		     location = pcl.get_x(direction);

		     // delete the exiting particle if out of box on the direction of OpenBC
		     if((dir_cnt%2==0 && location<delbry) ||(dir_cnt%2==1 && location>delbry))
		       delete_particle(pidx);
		     else{
		       pidx++;

		       //copy the particle within open boundary to inject particle list if their shifted location after 1 time step is within simulation box
		       if ((dir_cnt%2==0 && location<openbry) ||(dir_cnt%2==1 && location>openbry)){
		    	   double injx=pcl.get_x(0), injy=pcl.get_x(1), injz=pcl.get_x(2);
		    	   double inju=pcl.get_u(0), injv=pcl.get_u(1), injw=pcl.get_u(2);
		    	   double injq=pcl.get_q();

		    	   //shift 3 layers out, not mirror
		    	   if(direction == 0) injx = (dir_cnt%2==0) ?(injx-xLow):(injx+xLow);
		    	   if(direction == 1) injy = (dir_cnt%2==0) ?(injy-yLow):(injy+yLow);
		    	   if(direction == 2) injz = (dir_cnt%2==0) ?(injz-zLow):(injz+zLow);

		    	   injx = injx + inju*dt;
		    	   injy = injy + injv*dt;
		    	   injz = injz + injw*dt;

		    	   //Add particle if it enter that sub-domain or the domain box?
		    	   //assume create particle as long as it enters the domain box
		    	   if(injx>0 && injx<Lx && injy>0 && injy<Ly && injz>0 && injz<Lz){
		    		    injpcls.push_back(SpeciesParticle(inju,injv,injw,injq,injx,injy,injz,pclIDgenerator.generateID()));
		    	   }
		       }
		     }
		   }
	  }
  }

  //const int nop_remaining = getNOP();
  //const int nop_deleted = nop_orig - nop_remaining;
  const int nop_created = injpcls.size();

  //dprintf("change in # particles: %d - %d + %d = %d",nop_orig, nop_deleted, nop_created, nop_remaining);

  for(int outId=0;outId<nop_created;outId++)
	  _pcls.push_back(injpcls[outId]);
}

//Simply delete exiting test particles if openBC
void Particles3D::openbc_delete_testparticles()
{
  // if this is not a boundary process then there is nothing to do
  if(!vct->isBoundaryProcess_P()) return;

  //The below is OpenBC outflow for all other boundaries
  using namespace BCparticles;

  const bool openXleft = !vct->getPERIODICX_P() && vct->noXleftNeighbor_P() &&  (bcPfaceXleft == OPENBCOut||bcPfaceXleft == REEMISSION);
  const bool openYleft = !vct->getPERIODICY_P() && vct->noYleftNeighbor_P() &&  (bcPfaceYleft == OPENBCOut||bcPfaceYleft == REEMISSION);
  const bool openZleft = !vct->getPERIODICZ_P() && vct->noZleftNeighbor_P() &&  (bcPfaceZleft == OPENBCOut||bcPfaceZleft == REEMISSION);

  const bool openXright = !vct->getPERIODICX_P() && vct->noXrghtNeighbor_P() && (bcPfaceXright == OPENBCOut||bcPfaceXright == REEMISSION);
  const bool openYright = !vct->getPERIODICY_P() && vct->noYrghtNeighbor_P() && (bcPfaceYright == OPENBCOut||bcPfaceYright == REEMISSION);
  const bool openZright = !vct->getPERIODICZ_P() && vct->noZrghtNeighbor_P() && (bcPfaceZright == OPENBCOut||bcPfaceZright == REEMISSION);
  //dprintf( "Entered openbc_delete_testparticles %d %d %d %d %d %d", openXleft,openXright,openYleft,openYright,openZleft,openZright);
  if(!openXleft && !openYleft && !openZleft && !openXright && !openYright && !openZright)  return;

  const bool   apply_openBC[6]    = {openXleft, openXright,openYleft, openYright,openZleft, openZright};
  const double delete_boundary[6] = {0, Lx,0, Ly,0, Lz};

  for(int dir_cnt=0;dir_cnt<6;dir_cnt++){

    if(apply_openBC[dir_cnt]){

		  int pidx = 0;
		  int delnop = 0;
		  int direction = dir_cnt/2;
		  double delbry  = delete_boundary[dir_cnt];
		  double location;
		  while(pidx < getNOP())
		  {
		     SpeciesParticle& pcl = _pcls[pidx];
		     location = pcl.get_x(direction);

		     // delete the exiting particle if out of box on the direction of OpenBC
		     if((dir_cnt%2==0 && location<delbry) ||(dir_cnt%2==1 && location>delbry)){
		    	 delete_particle(pidx);
		    	 delnop++;
		     }
		     else {
		    	 pidx ++;
		     }
		   }
		  dprintf( "*** Delete %d Test Particle for OpenBC for Direction %d on particle species %d",delnop, dir_cnt, ns);
	  }

  }

}


/** apply a linear perturbation to particle distribution */
//
// This was not being used, so I took it out. -eaj
//
//void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Field * EMf) {
//
//  double value1 = 0.0, value2 = 0.0, max_value = 0.0, min_value = 0.0, phi, n;
//  // rescaling of amplitudes according to deltaBoB //
//
//  const double alpha = deltaBoB * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0) + EMf->getBz(1, 1, 0) * EMf->getBz(1, 1, 0)) / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);
//
//  Ex_mod *= alpha;
//  Ey_mod *= alpha;
//  Ez_mod *= alpha;
//  Bx_mod *= alpha;
//  By_mod *= alpha;
//  Bz_mod *= alpha;
//
//
//
//  // find the maximum value of f=1+delta_f/f0
//  for (register double vpar = -2 * uth; vpar <= 2 * uth; vpar += 0.0005)
//    for (register double vperp = 1e-10; vperp <= 2 * vth; vperp += 0.0005)
//      for (register double X = xstart; X <= xend; X += 2 * grid->getDX())
//        for (register double Y = ystart; Y <= yend; Y += 2 * grid->getDY()) {
//          value1 = 1 + delta_f(vpar, vperp, 0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(vpar, vperp);
//
//          if (value1 > max_value)
//            max_value = value1;
//
//
//        }
//
//
//
//  max_value *= 3.2;
//  phi = 1.48409;
//  n = 2.948687;                 // security factor...
//  if (ns == 1) {
//    max_value *= 3.0;
//    phi = -1.65858;
//    n = 2.917946;
//  }                             // security factor...
//  cout << "max-value=" << max_value << " min-value=" << min_value << endl;
//
//  /* initialize random generator */
//  srand(vct->getCartesian_rank() + 2);
//
//  int counter = 0;
//  int total_generated = 0;
//  const double q_value = (qom / fabs(qom)) * ((0.19635) / npcel) * (1.0 / invVOL);
//  for (int i = 1; i < grid->getNXC() - 1; i++)
//  for (int j = 1; j < grid->getNYC() - 1; j++)
//  {
//    const double x_factor = (dx / (npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))))) + grid->getXN(i, j, 0);
//    const double y_factor = (dy / npcely) + grid->getYN(i, j, 0);
//    for (int ii = 0; ii < npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))); ii++)
//    for (int jj = 0; jj < npcely; jj++)
//    {
//      const double x = (ii + .5) * x_factor;
//      const double y = (jj + .5) * y_factor;
//      eprintf("this was not being initialized.");
//      const double z = 0.;
//
//      // apply rejection method in velocity space
//      double u,v,w;
//      bool rejected = true;
//      while (rejected) {
//        total_generated++;
//        sample_maxwellian(
//          u,v,w,
//          uth,vth,wth,
//          u0,v0,w0);
//
//        // test: if rand < (1+delta_f/f0)/max_value --> accepted
//        if (rand() / (double) RAND_MAX <= (1 + delta_f(u, v, w, x, y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(u, sqrt(v * v + w * w))) / max_value)
//          rejected = false;
//      }
//      counter++;
//      add_new_particle(u,v,w,q_value,x,y,z,0.);
//    }
//  }
//  assert_eq(getNOP(),counter+1);
//  // if (vct->getCartesian_rank()==0)
//  cout << "Rejection method: " << getNOP() / double (total_generated) * 100 << " % of particles are accepted for species " << ns << " counter=" << counter << endl;
//}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double theta, Field * EMf) {
  const complex < double >I(0.0, 1.0);
  const double vperp = sqrt(v * v + w * w);
  const double vpar = u;
  const double kpar = kx;
  double kperp;
  if (ky == 0.0)                // because this formula is not valid for exactly parallel
    kperp = 1e-9;
  else
    kperp = ky;
  const double om_c = qom / c * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0)) / 2 / M_PI;
  const double phi = atan2(w, v);
  const double lambda = kperp * vperp / om_c;
  const complex < double >omega(omega_re, omega_i);

  const int lmax = 5;           // sum from -lmax to lmax

  double bessel_Jn_array[lmax + 2];
  double bessel_Jn_prime_array[lmax + 1];
  complex < double >a1[2 * lmax + 1], a2[2 * lmax + 1], a3[2 * lmax + 1];
  complex < double >factor, deltaf;

  // rotation of x,y
  double temp;
  temp = x;
  x = x * cos(theta) - y * sin(theta);
  y = temp * sin(theta) + y * cos(theta);


  /** for compilation issues comment this part: PUT in the math stuff */
  // calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
  factor = (kpar * vperp / omega * df0_dvpar(vpar, vperp) + (1.0 - (kpar * vpar / omega)) * df0_dvperp(vpar, vperp));
  for (register int l = -lmax; l < 0; l++) {  // negative index
    a1[l + lmax] = factor / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * 0.5 * pow(-1.0, -l) * (bessel_Jn_array[-l - 1] - bessel_Jn_array[-l + 1]);
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * pow(-1.0, -l) * bessel_Jn_array[-l];
  }

  for (register int l = 0; l < lmax + 1; l++) { // positive index
    a1[l + lmax] = factor / lambda * bessel_Jn_array[l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * bessel_Jn_prime_array[l];
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * bessel_Jn_array[l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * bessel_Jn_array[l];
  }

  deltaf = (0.0, 0.0);
  for (register int l = -lmax; l < lmax + 1; l++) {
    deltaf += (a3[l + lmax] * Ex_mod * exp(I * Ex_phase) + a1[l + lmax] * Ey_mod * exp(I * Ey_phase) + a2[l + lmax] * Ez_mod * exp(I * Ez_phase)) / (kpar * vpar + l * om_c - omega) * exp(-I * phi * (double) l);
  }
  deltaf *= I * qom * exp(I * lambda * sin(phi)) * exp(I * (2 * M_PI * kx * x + 2 * M_PI * ky * y));

  return (real(deltaf));
}

double Particles3D::df0_dvpar(double vpar, double vperp) {
  double result;
  result = -2 * (vpar - u0) / uth / uth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::df0_dvperp(double vpar, double vperp) {
  double result;
  result = -2 * (vperp) / vth / vth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::f0(double vpar, double vperp) {
  double result;
  result = exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

void Particles3D::RotatePlaneXY(double theta) {
  double temp, temp2;
  for (register int s = 0; s < getNOP(); s++) {
    temp = u[s];
    temp2 = v[s];
    u[s] = temp * cos(theta) + v[s] * sin(theta);
    v[s] = -temp * sin(theta) + temp2 * cos(theta);
  }
}

/*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
double Particles3D::deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center)
{
  int pidx = 0;
  double Q_removed=0.;
  while (pidx < _pcls.size())
  {
    SpeciesParticle& pcl = _pcls[pidx];
    double xd = pcl.get_x() - x_center;
    double yd = pcl.get_y() - y_center;
    double zd = pcl.get_z() - z_center;

    if ( (xd*xd+yd*yd+zd*zd) < R*R ){
      Q_removed += pcl.get_q();
      delete_particle(pidx);
    } else {
      pidx++;
    }
  }
  return(Q_removed);
}

double Particles3D::deleteParticlesInsideSphere2DPlaneXZ(double R, double x_center, double z_center)
{
  int pidx = 0;
  double Q_removed=0.;
  while (pidx < _pcls.size())
  {
    SpeciesParticle& pcl = _pcls[pidx];
    double xd = pcl.get_x() - x_center;
    double zd = pcl.get_z() - z_center;

    if ( (xd*xd+zd*zd) < R*R ){
      Q_removed += pcl.get_q();
      delete_particle(pidx);
    } else {
      pidx++;
    }
  }
  return(Q_removed);
}
