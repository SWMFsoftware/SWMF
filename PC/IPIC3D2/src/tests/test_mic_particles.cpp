#define MICVEC_DEFINE_OUTPUT_OPERATORS
#include <iostream>
#include <omp.h>
#include <stdio.h>
//#include "time.h" // for clock_gettime()
#include <stdint.h> // for uint64_t
#include <stdlib.h> // rand()
#include <sys/time.h>
#include <micvec.h>
#include <assert.h>
#include "Alloc.h"
#include "mic_particles.h"
#include "Particle.h"
#include "Grid.h"
#include "../utility/debug.cpp"
#include "../utility/asserts.cpp"
#include "../utility/errors.cpp"

#define printexpr(var) std::cout << "line " << __LINE__ << ": " \
  << #var << " = " << var << std::endl;

#define printarr(var,numel) std::cout << "line " << __LINE__ << ": " \
  << #var << " = ("; \
  for(int i=0;i<(numel-1);i++) std::cout << var[i] << ", "; \
  std::cout << var[numel-1] << ")" << std::endl;

inline void get_field_components_for_cell(
  arr1_double_get field_components[8],
  const_arr4_double fieldForPcls,
  int cx,int cy,int cz)
{
  // interface to the right of cell
  const int ix = cx+1;
  const int iy = cy+1;
  const int iz = cz+1;

  // is this faster?
  //
  //field_components[0] = fieldForPcls[ix][iy][iz]; // field000
  //field_components[1] = fieldForPcls[ix][iy][cz]; // field001
  //field_components[2] = fieldForPcls[ix][cy][iz]; // field010
  //field_components[3] = fieldForPcls[ix][cy][cz]; // field011
  //field_components[4] = fieldForPcls[cx][iy][iz]; // field100
  //field_components[5] = fieldForPcls[cx][iy][cz]; // field101
  //field_components[6] = fieldForPcls[cx][cy][iz]; // field110
  //field_components[7] = fieldForPcls[cx][cy][cz]; // field111
  //
  // or is this?
  //
  // creating these aliases seems to accelerate this method (by about 30%?)
  // on the Xeon host processor, suggesting deficiency in the optimizer.
  //
  arr3_double_fetch field0 = fieldForPcls[ix];
  arr3_double_fetch field1 = fieldForPcls[cx];
  arr2_double_fetch field00 = field0[iy];
  arr2_double_fetch field01 = field0[cy];
  arr2_double_fetch field10 = field1[iy];
  arr2_double_fetch field11 = field1[cy];
  field_components[0] = field00[iz]; // field000 
  field_components[1] = field00[cz]; // field001 
  field_components[2] = field01[iz]; // field010 
  field_components[3] = field01[cz]; // field011 
  field_components[4] = field10[iz]; // field100 
  field_components[5] = field10[cz]; // field101 
  field_components[6] = field11[iz]; // field110 
  field_components[7] = field11[cz]; // field111 
}

using namespace iPic3D;
using namespace std;

const double dx=1.;
const double dy=1.;
const double dz=1.;
const double invdx=1./dx;
const double invdy=1./dy;
const double invdz=1./dz;
const double xstart=0.;
const double ystart=0.;
const double zstart=0.;
const double dt=.1;
const double c =1.;
const double qom =1.;
const double dto2 = .5 * dt;
const double qdto2mc = qom * dto2 / c;
const int NiterMover=4;
const int nxc = 10;
const int nyc = 12;
const int nzc = 6;
const int nxn = nxc+1;
const int nyn = nyc+1;
const int nzn = nzc+1;
const int DFIELD_3or4 = 4;
array4<pfloat> fieldForPcls(nxn, nyn, nzn, 2*DFIELD_3or4);
const int nop=2;
SpeciesParticle pclsA[nop] ALLOC_ALIGNED;
SpeciesParticle pclsB[nop] ALLOC_ALIGNED;
SpeciesParticle* pcls = pclsA;

Grid grid(nxc,nyc,nzc,dx,dy,dz,xstart,ystart,zstart);

void init_fieldForPcls()
{
  const int nxn = fieldForPcls.dim1();
  const int nyn = fieldForPcls.dim2();
  const int nzn = fieldForPcls.dim3();
  const int nnn = fieldForPcls.dim4();
  dprint(fieldForPcls.get_size());
  dprint(nxn);
  dprint(nyn);
  dprint(nzn);
  dprint(nnn);
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    fieldForPcls[i][j][k][0] = i+.2;
    fieldForPcls[i][j][k][1] = j+.3;
    fieldForPcls[i][j][k][2] = k+.4;
    fieldForPcls[i][j][k][0+DFIELD_3or4] = i+.5;
    fieldForPcls[i][j][k][1+DFIELD_3or4] = j+.6;
    fieldForPcls[i][j][k][2+DFIELD_3or4] = k+.7;
  }
}

void init_pcls()
{
  for(int i=0;i<nop;i++)
  {
    double u=10;
    double v=20;
    double w=30;
    double q=10;
    double x=1;
    double y=2;
    double z=3;
    double t=0;
    pcls[i].set(
      u+.1*i, v+.1*i, w+.1*i, q+.1*i,
      x+.1*i, y+.1*i, z+.1*i, t+.1*i);
    pclsB[i] = pclsA[i];
  }
}

void init_state()
{
  init_fieldForPcls();
  init_pcls();
}

void test_get_field_components_for_cell()
{
  F64vec8 field_components0[8];
  F64vec8 field_components1[8];
  I32vec16 cx(0,0,0,0,0,0,0,0,0,3,2,1,0,4,3,2);
  get_field_components_for_cell(
    field_components0,
    field_components1,
    fieldForPcls,
    cx);
  std::cout << "cx=" << cx << endl;
  for(int i=0;i<8;i++)
  {
    std::cout << "field_components0[" << i << "]=" << field_components0[i] << endl;
  }
  for(int i=0;i<8;i++)
  {
    std::cout << "field_components1[" << i << "]=" << field_components1[i] << endl;
  }
}

// import numpy;
// X0 = numpy.array([.1, .2, .5, 0, .5, .5, .5, 0]);
// #X0 = X0r[::-1];
// X1 = 1 - X0;
// weights0 = range(8);
// weights0[0] = X0[0]*X0[1]*X0[2]; # weight000
// weights0[1] = X0[0]*X0[1]*X1[2]; # weight001
// weights0[2] = X0[0]*X1[1]*X0[2]; # weight010
// weights0[3] = X0[0]*X1[1]*X1[2]; # weight011
// weights0[4] = X1[0]*X0[1]*X0[2]; # weight100
// weights0[5] = X1[0]*X0[1]*X1[2]; # weight101
// weights0[6] = X1[0]*X1[1]*X0[2]; # weight110
// weights0[7] = X1[0]*X1[1]*X1[2]; # weight111   
// 
void test_construct_weights_for_2pcls()
{
  F64vec8 X0(0,.1,.2,.5,0,.5,.2,.1);
  F64vec8 weights[2];
  construct_weights_for_2pcls(weights, X0);
  for(int i=0;i<8;i++)
  {
    std::cout << "weights0[" << i << "]=" << weights[0][i] << endl;
  }
  for(int i=0;i<8;i++)
  {
    std::cout << "weights1[" << i << "]=" << weights[1][i] << endl;
  }
}

void test_maximum()
{
  F64vec8 u(1,-2,4,5,-6,7,-8,-9);
  F64vec8 v(1,2,-4,5,-6,7,8,-9);
  F64vec8 mx = maximum(u,v);
  F64vec8 mn = minimum(u,v);
  cout << "u=" << u << endl;
  cout << "v=" << v << endl;
  cout << "mx=" << mx << endl;
  cout << "mn=" << mn << endl;
}

// move the particle using MIC vector intrinsics
void test_version_of_mover_PC_AoS_vec_intr()
{
  dprintf("=== entering ===");
  // Here and below x stands for all 3 physical position coordinates
  // and u stands for all 3 velocity coordinates.
  const F64vec8 dx_inv = make_F64vec8(invdx,invdy,invdz);
  const F64vec8 pdom_xlow = make_F64vec8(xstart,ystart,zstart);
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

  ASSUME_ALIGNED(pcls);
  const F64vec8 dto2 = F64vec8(::dto2);
  const F64vec8 qdto2mc = F64vec8(::qdto2mc);
  for (int pidx = 0; pidx < nop; pidx+=2)
  {
    dprintf("processing particles %d and %d",pidx,pidx+1);
    // copy the particle
    SpeciesParticle* pcl = &(pcls[pidx]);
    //ALIGNED(pcl);

    // gather position and velocity data from particles
    //
    F64vec8 pcl0 = *(F64vec8*)&(pcl[0]);
    F64vec8 pcl1 = *(F64vec8*)&(pcl[1]);
    //F64vec8 pcl0 = *(F64vec8*)&(pcls[pidx]);
    //F64vec8 pcl1 = *(F64vec8*)&(pcls[pidx+1]);
    const F64vec8 xorig = cat_hgh_halves(pcl0,pcl1);
    F64vec8 xavg = xorig;
    const F64vec8 uorig = cat_low_halves(pcl0,pcl1);

    printexpr(uorig);
    printexpr(xorig);

    // calculate the average velocity iteratively
    //
    // (could stop iteration when it is determined that
    // both particles are converged, e.g. if change in
    // xavg is sufficiently small)
    F64vec8 uavg;
    for (int iter = 0; iter < NiterMover; iter++)
    {
      dprintf("starting iteration %d",iter);
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
      const F64vec8 X = gX - cellXstart;
      //printexpr(X);
      F64vec8 weights[2];
      construct_weights_for_2pcls(weights, X);
      //printexpr(weights[0]);
      //printexpr(weights[1]);

      // interpolate field to get fields
      F64vec8 fields[2];
      // sample fields for first particle
      fields[0] = sample_field_mic(weights[0],field_components0);
      // sample fields for second particle
      fields[1] = sample_field_mic(weights[1],field_components1);
      const F64vec8 B = cat_low_halves(fields[0],fields[1]);
      const F64vec8 E = cat_hgh_halves(fields[0],fields[1]);
      //printexpr(B);
      //printexpr(E);

      // use sampled field to push particle block
      //
      uavg = compute_uvg_for_2pcls(uorig, B, E, qdto2mc);
      // update average position
      xavg = xorig + uavg*dto2;
      printexpr(uavg);
      printexpr(xavg);
    } // end of iterative particle advance
    // update the final position and velocity
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
}

void test_version_of_mover_PC_AoS()
{
  dprintf("=== entering ===");
  for (int pidx = 0; pidx < nop; pidx++)
  {
    dprintf("processing particle %d",pidx);
    // copy the particle
    SpeciesParticle* pcl = &pclsB[pidx];
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
    dprintf("uorig: %g, %g, %g",uorig,vorig,worig);
    dprintf("xorig: %g, %g, %g",xorig,yorig,zorig);
    Grid* grid = &::grid;
    // calculate the average velocity iteratively
    for (int iter = 0; iter < NiterMover; iter++)
    {
      dprintf("starting iteration %d",iter);
      // compute weights for field components
      //
      double weights[8];
      int cx,cy,cz;
      dprintf("xavg: %g, %g, %g",xavg,yavg,zavg);
      grid->get_safe_cell_and_weights(xavg,yavg,zavg,cx,cy,cz,weights);
      //printarr(weights,8);

      arr1_double_get field_components[8];
      get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);

      double Ex = 0.0; double Ey = 0.0; double Ez = 0.0;
      double Bx = 0.0; double By = 0.0; double Bz = 0.0;
      for(int c=0; c<8; c++)
      {
        Bx += weights[c] * field_components[c][0];
        By += weights[c] * field_components[c][1];
        Bz += weights[c] * field_components[c][2];
        Ex += weights[c] * field_components[c][0+DFIELD_3or4];
        Ey += weights[c] * field_components[c][1+DFIELD_3or4];
        Ez += weights[c] * field_components[c][2+DFIELD_3or4];
      }
      const double Omx = qdto2mc*Bx;
      const double Omy = qdto2mc*By;
      const double Omz = qdto2mc*Bz;
      //dprintf("B: %g, %g, %g",Bx,By,Bz);
      //dprintf("E: %g, %g, %g",Bx,By,Bz);

      // end interpolation
      const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
      const pfloat denom = 1.0 / (1.0 + omsq);
      // solve the position equation
      const pfloat ut = uorig + qdto2mc * Ex;
      const pfloat vt = vorig + qdto2mc * Ey;
      const pfloat wt = worig + qdto2mc * Ez;
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
      dprintf("uavg: %g, %g, %g",uavg,vavg,wavg);
      dprintf("xavg: %g, %g, %g",xavg,yavg,zavg);
    }                           // end of iteration
    // update the final position and velocity
    pcl->set_x(xorig + uavg * dt);
    pcl->set_y(yorig + vavg * dt);
    pcl->set_z(zorig + wavg * dt);
    pcl->set_u(2.0 * uavg - uorig);
    pcl->set_v(2.0 * vavg - vorig);
    pcl->set_w(2.0 * wavg - worig);
  }
}

Grid3DCU::Grid3DCU(
  int nxc_, int nyc_, int nzc_,
  double dx_, double dy_, double dz_,
  double xStart_, double yStart_, double zStart_)
: nxc(nxc_), nyc(nyc_), nzc(nzc_),
  dx(dx_), dy(dy_), dz(dz_),
  xStart(xStart_), yStart(yStart_), zStart(zStart_)
{
  const double xWidth = dx*(nxc-2);
  const double yWidth = dy*(nyc-2);
  const double zWidth = dz*(nzc-2);

  xEnd = xStart + xWidth;
  yEnd = yStart + yWidth;
  zEnd = zStart + zWidth;

  init_derived_parameters();
}

// set derived convenience
void Grid3DCU::init_derived_parameters()
{
  // calculation conveniences
  //
  invVOL = 1.0 / (dx * dy * dz);
  invdx = 1.0 / dx;
  invdy = 1.0 / dy;
  invdz = 1.0 / dz;
  //
  nxn = nxc + 1;
  nyn = nyc + 1;
  nzn = nzc + 1;
  cxlast = nxc-1;
  cylast = nyc-1;
  czlast = nzc-1;

  // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
  pfloat_node_xcoord = new pfloat[nxn];
  pfloat_node_ycoord = new pfloat[nyn];
  pfloat_node_zcoord = new pfloat[nzn];
  node_xcoord = new double[nxn];
  node_ycoord = new double[nyn];
  node_zcoord = new double[nzn];
  for (int i=0; i<nxn; i++) node_xcoord[i] = xStart + (i - 1) * dx;
  for (int j=0; j<nyn; j++) node_ycoord[j] = yStart + (j - 1) * dy;
  for (int k=0; k<nzn; k++) node_zcoord[k] = zStart + (k - 1) * dz;
  for (int i=0; i<nxn; i++) pfloat_node_xcoord[i] = node_xcoord[i];
  for (int j=0; j<nyn; j++) pfloat_node_ycoord[j] = node_ycoord[j];
  for (int k=0; k<nzn; k++) pfloat_node_zcoord[k] = node_zcoord[k];
  // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
  center_xcoord = new double[nxc];
  center_ycoord = new double[nyc];
  center_zcoord = new double[nzc];
  for(int i=0; i<nxc; i++) center_xcoord[i] = .5*(node_xcoord[i]+node_xcoord[i+1]);
  for(int j=0; j<nyc; j++) center_ycoord[j] = .5*(node_ycoord[j]+node_ycoord[j+1]);
  for(int k=0; k<nzc; k++) center_zcoord[k] = .5*(node_zcoord[k]+node_zcoord[k+1]);
}

/** deallocate the local grid */
Grid3DCU::~Grid3DCU() {
  delete [] node_xcoord;
  delete [] node_ycoord;
  delete [] node_zcoord;
  delete [] center_xcoord;
  delete [] center_ycoord;
  delete [] center_zcoord;
}

int main()
{
  init_state();
  // printf("#\n");
  // printf("# calling test_get_field_components_for_cell:\n");
  // printf("#\n");
  // test_get_field_components_for_cell();
  // printf("#\n");
  // printf("# calling test_construct_weights_for_2pcls:\n");
  // printf("#\n");
  // test_construct_weights_for_2pcls();
  // printf("#\n");
  // printf("# calling test_maximum:\n");
  // printf("#\n");
  // test_maximum();
  test_version_of_mover_PC_AoS_vec_intr();
  test_version_of_mover_PC_AoS();
}
