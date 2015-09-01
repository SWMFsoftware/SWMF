/*******************************************************************************************
  Particles3D.h  -  Class for particles of the same species, in a 2D space and 3 component velocity
  -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
 ********************************************************************************************/

#ifndef Part2D_H
#define Part2D_H

#include "Particles3Dcomm.h"
//#include "TimeTasks.h"

/**
 * 
 * Class for particles of the same species, in a 2D space and 3 component velocity
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3D:public Particles3Dcomm {

  public:
    /** constructor */
    //Particles3D();
    Particles3D(int species, CollectiveIO *col, VirtualTopology3D *vct, Grid * grid):
      Particles3Dcomm(species, col, vct, grid)
    {}
    /** destructor */
    ~Particles3D(){}
    /** Initial condition: uniform in space and motionless */
    void uniform_background(Field * EMf);
    /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
      <ul>
      <li> dim = 0 --> constant velocity on X direction </li>
      <li> dim = 1 --> constant velocity on Y direction </li>
      </ul>
      */
    void constantVelocity(double vel, int dim, Field * EMf);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void maxwellian(Field * EMf);
    /** pitch_angle_energy initialization (Assume B on z only) for test particles */
    void pitch_angle_energy(Field * EMf);
    /** Force Free initialization (JxB=0) for particles */
    void force_free(Field * EMf);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void alt_maxwellian(Field * EMf);
    /** Linear_perturbation */
    //void linear_perturbation(double deltaBX, double kx, double ky, double theta, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Field * EMf);
    /**Add a periodic perturbation in velocity exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
    void AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0);
    /** Linear delta f for bi-maxwellian plasma */
    double delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_ampl, double Ex_phase, double Ey_ampl, double Ey_phase, double Ez_ampl, double Ez_phase, double theta, Field * EMf);
    /** Derivative of f0 wrt vpar */
    double df0_dvpar(double vpar, double vperp);
    /** Derivative of f0 wrt vperp */
    double df0_dvperp(double vpar, double vperp);
    /** Equilibrium bi-maxwellian f0 */
    double f0(double vpar, double vperp);
    /** Rotate velocities in plane XY of angle theta */
    void RotatePlaneXY(double theta);
    /** mover with the esplicit non relativistic scheme */
    void mover_explicit(Field * EMf);
    /** mover with a Predictor-Corrector Scheme */
    void mover_PC(Field * EMf);
    /** array-of-structs version of mover_PC */
    void mover_PC_AoS(Field * EMf);
    /** Relativistic array-of-structs version of mover_PC with adaptive Subcycling and PC*/
    void mover_PC_AoS_Relativistic(Field * EMf);
    /* vectorized version of previous */
    void mover_PC_AoS_vec(Field * EMf);
    /* mic particle mover */
    void mover_PC_AoS_vec_intr(Field * EMf);
    /* this computes garbage */
    void mover_PC_AoS_vec_onesort(Field * EMf);
    /** vectorized version of mover_PC **/
    void mover_PC_vectorized(Field * EMf);
    /** relativistic mover with a Predictor-Corrector scheme */
    int mover_relativistic(Field * EMf);
   private:
    /** repopulate particles in a single cell */
    void populate_cell_with_particles(int i, int j, int k, double q,
      double dx_per_pcl, double dy_per_pcl, double dz_per_pcl);
   public:
    /** repopulate particles in boundary layer */
    void repopulate_particles();
    /*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
    double deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center);
    double deleteParticlesInsideSphere2DPlaneXZ(double R, double x_center, double z_center);
    /**Particles Open Boundary */
    void openbc_particles_outflow();
    void openbc_delete_testparticles();
    void openbc_particles_inflow();

#ifdef BATSRUS
    /*! Initial condition: given a fluid model (BATSRUS) */
    void MaxwellianFromFluid(Field* EMf,Collective *col, int is);
    /*! Initiate dist. func. for a single cell form a fluid model (BATSRUS) */
    void MaxwellianFromFluidCell(Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, longid* ParticleID);
#endif

};

#endif
