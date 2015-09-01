/*******************************************************************************************************
  Grid3DCUCU.h  -  uniform cartesian 3D local grid for each process, including che guard cells
  -------------------

 *******************************************************************************************************/

#ifndef GRID3DCU_H
#define GRID3DCU_H

#include "arraysfwd.h"
#include "ipicfwd.h"
#include "math.h" // for floor
#include "assert.h"

/**
 * Uniform cartesian local grid 3D
 *
 * @date Fri Jun 4 2008
 * @par Copyright:
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 3.0
 *
 */
class Grid3DCU                  // :public Grid
{
public:
  /** constructor */
  Grid3DCU(CollectiveIO * col, VirtualTopology3D * vct);
  Grid3DCU(
    int nxc_, int nyc_, int nzc_,
    double dx_, double dy_, double dz_,
    double xStart_, double yStart_, double zStart_);
  /** destructor */
  ~Grid3DCU();
  /** allocate grid arrays for this domain */
  //void allocate(CollectiveIO * ptC, VirtualTopology3D * ptVCT);
  /** deallocate grid arrays for this domain */
  //void deallocate();
  /** print grid info */
  void print()const;
  /** calculate a derivative along a direction on nodes */
  //void derivN(arr3_double derN,
  //  const_arr4_double scFieldC, int ns, int dir)const;
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradC2N(arr3_double gradXN, arr3_double gradYN, arr3_double gradZN,
    const_arr3_double scFieldC)const;
  /** calculate gradient on nodes, given a scalar field defined on central points  */
  void gradN2C(arr3_double gradXC, arr3_double gradYC, arr3_double gradZC,
    const_arr3_double scFieldN)const;
  /** calculate divergence on central points, given a vector field defined on nodes  */
  void divN2C(arr3_double divC,
    const_arr3_double vecFieldXN,
    const_arr3_double vecFieldYN,
    const_arr3_double vecFieldZN)const;
  /** calculate divergence on nodes, given a vector field defined on central points  */
  void divC2N(arr3_double divN,
    const_arr3_double vecFieldXC,
    const_arr3_double vecFieldYC,
    const_arr3_double vecFieldZC)const;
  /** calculate curl on nodes, given a vector field defined on central points  */
  void curlC2N(arr3_double curlXN, arr3_double curlYN,
    arr3_double curlZN,
    const_arr3_double vecFieldXC,
    const_arr3_double vecFieldYC,
    const_arr3_double vecFieldZC)const;
  /** calculate curl on central points, given a vector field defined on nodes  */
  void curlN2C(arr3_double curlXC, arr3_double curlYC, arr3_double curlZC,
    const_arr3_double vecFieldXN,
    const_arr3_double vecFieldYN,
    const_arr3_double vecFieldZN)const;

  /** calculate divergence on central points, given a Tensor field defined on nodes  */
  void divSymmTensorN2C(arr3_double divCX, arr3_double divCY, arr3_double divCZ,
    const_arr4_double pXX,
    const_arr4_double pXY,
    const_arr4_double pXZ,
    const_arr4_double pYY,
    const_arr4_double pYZ,
    const_arr4_double pZZ, int ns)const;

  /** calculate laplacian on nodes, given a scalar field defined on nodes */
  void lapN2N(arr3_double lapN, const_arr3_double scFieldN, EMfields3D *EMf)const;
  /** calculate laplacian on central points, given a scalar field defined on central points for Poisson */
  //void lapC2Cpoisson(arr3_double lapC, arr3_double scFieldC)const;
  void lapC2Cpoisson(arr3_double lapC, arr3_double scFieldC, EMfields3D *EMf)const;
  /** calculate laplacian on central points, given a scalar field defined on central points */
  void lapC2C(arr3_double lapC, const_arr3_double scFieldC)const;

  /** calculate divergence on boundaries */
  void divBCleft(arr3_double divBC,
    const_arr3_double vectorX,
    const_arr3_double vectorY,
    const_arr3_double vectorZ, int leftActiveNode, int dirDER)const;
  /** calculate divergence on boundaries */
  void divBCright(arr3_double divBC,
    const_arr3_double vectorX,
    const_arr3_double vectorY,
    const_arr3_double vectorZ, int rightActiveNode, int dirDER)const;
  /** calculate derivative on boundaries */
  void derBC(arr3_double derBC,
    const_arr3_double vector, int leftActiveNode, int dirDER)const;


  /** interpolate on nodes from central points */
  void interpC2N(arr3_double vecFieldN, const_arr3_double vecFieldC)const;
  /** interpolate on central points from nodes */
  void interpN2C(arr3_double vecFieldC, const_arr3_double vecFieldN)const;
  /** interpolate on central points from nodes */
  void interpN2C(arr4_double vecFieldC, int ns, const_arr4_double vecFieldN)const;

private:
  void init_derived_parameters();

  // /////////// PRIVATE VARIABLES //////////////
  //
  // Grid conventions:
  //
  //   [The purpose of this note is not so much to teach as
  //   to establish a consistent vocabulary and set of naming
  //   conventions for grid parameters.]
  //
  //   The purpose of iPic3D is to solve a Vlasov-Maxwell
  //   plasma problem.  Each problem is solved on a box-shaped
  //   domain and is defined by initial conditions and boundary
  //   conditions.
  //
  //   Conceptually, the full problem domain is solved on
  //   a "guarded" global grid.  The guarded global grid
  //   consists of proper grid cells (that partition the domain)
  //   surrounded by a layer of "guard" or "ghost" cells that are
  //   used to handle boundary conditions.  At the beginning of
  //   each cycle of the algorithm, the guard cells are populated
  //   with data based on the choice of boundary conditions of
  //   the problem.  Refer to the layer of proper cells just
  //   inside the guard cells as "(proper) boundary cells".
  //   Guard cell data is generally populated with data from
  //   boundary cells.  For periodic boundary conditions, guard
  //   cell data is simply copied from the appropriate boundary
  //   cell.  For conducting wall and outflow boundaries, each
  //   guard cell is populated based on the data from the
  //   neighboring boundary cell, probably with reflection of
  //   location and often with negation of some values, depending
  //   on the type of boundaries involved.
  //
  //   The corners of the mesh cells are referred to as
  //   "nodes".  In iPic3D, electromagnetic field quantities are
  //   canonically represented as sample values at the nodes, and
  //   the value of the electromagnetic field in the interior of
  //   a grid cell is determined by bilinear interpolation of
  //   the field values at the nodes.  Bilinear interpolation
  //   is also used when converting between a representation as
  //   cell-average quantities and a representation as nodal
  //   samples.  This maintains second-order accuracy in space.
  //
  //   To parallelize with MPI, the full domain of the problem
  //   is partitioned into "(proper) subdomains", one for
  //   each process, such that the global grid partitions the
  //   subdomains.  Refer to the grid cells within a subdomain
  //   as its "(proper) subgrid".  Refer to the outer layer of
  //   cells within each subgrid as "boundary cells".  The grid
  //   defined in each process is a "guarded subgrid", which
  //   consists of the proper subgrid surrounded by a layer of
  //   guard cells corresponding either to the boundary cells of
  //   the neighboring subdomain or to ghost cells of the guarded
  //   global grid.
  //
  //   The guard cells thus lie outside the proper physical
  //   subdomain for which the process is responsible.  Guard
  //   cells are used to handle boundary conditions and
  //   interprocess communication.  At the beginning of each
  //   cycle of the algorithm, the guard cells are repopulated.
  //   Guard cells outside the problem domain are populated using
  //   boundary conditions.  Any guard cell inside the problem
  //   domain is populated from the appropriate boundary cell of
  //   the appropriate process.
  //
  // Conventions for naming of grid quantities:
  //
  //   n = number
  //   x = x-direction
  //   y = y-direction
  //   z = z-direction
  //   c = cell of grid
  //   n = node of grid
  //   
  //   nxc = number of x cells (in guarded grid)
  //   nxn = number of x nodes (in guarded grid)
  //
  //   nxc_g = number of x cells in guarded/"greater"/padded subgrid
  //   nxc_r = number of x cells in "real"/"restricted"/proper subgrid
  //   
private:
  const static bool suppress_runaway_particle_instability = true;
  const VirtualTopology3D& _vct;
  /** number of cells - X direction, including + 2 (guard cells) */
  int nxc;
  /** number of nodes - X direction, including + 2 extra nodes for guard cells */
  int nxn;
  /** number of cell - Y direction, including + 2 (guard cells) */
  int nyc;
  /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
  int nyn;
  /** number of cell - Z direction, including + 2 (guard cells) */
  int nzc;
  /** number of nodes - Z direction, including + 2 extra nodes for guard cells */
  int nzn;
  /** dx = space step - X direction */
  double dx;
  /** dy = space step - Y direction */
  double dy;
  /** dz = space step - Z direction */
  double dz;
  /** invdx = 1/dx */
  double invdx;
  /** invdy = 1/dy */
  double invdy;
  /** invdz = 1/dz */
  double invdz;
  /** volume of mesh cell*/
  double VOL;
  /** invol = inverse of volume*/
  double invVOL;
  /** index of last cell including ghost cells */
  // (precomputed for speed)
  int cxlast; // nxc-1;
  int cylast; // nyc-1;
  int czlast; // nzc-1;
  // number of cells excluding guard cells
  // (i.e. restricted to proper subdomain)
  int nxc_r;
  int nyc_r;
  int nzc_r;
  // number of subdomain cells in a regular (untruncated) subdomain grid
  int num_cells_rr;
  /** node coordinate */
  pfloat *pfloat_node_xcoord;
  pfloat *pfloat_node_ycoord;
  pfloat *pfloat_node_zcoord;
  double *node_xcoord;
  double *node_ycoord;
  double *node_zcoord;
  /** center coordinate */
  double *center_xcoord;
  double *center_ycoord;
  double *center_zcoord;
  /** local grid boundaries coordinate of proper subdomain */
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  double xStart_g, yStart_g, zStart_g;
  double epsilon;
  double nxc_minus_epsilon;
  double nyc_minus_epsilon;
  double nzc_minus_epsilon;

public: // accessors (inline)
  const VirtualTopology3D& get_vct()const{return _vct;}
  int getNXC()const{ return (nxc); }
  int getNXN()const{ return (nxn); }
  int getNYC()const{ return (nyc); }
  int getNYN()const{ return (nyn); }
  int getNZC()const{ return (nzc); }
  int getNZN()const{ return (nzn); }
  int get_nxc_r()const{return nxc_r;}
  int get_nyc_r()const{return nyc_r;}
  int get_nzc_r()const{return nzc_r;}
  int get_num_cells_rr()const{return num_cells_rr;}
  double getDX()const{ return (dx); }
  double getDY()const{ return (dy); }
  double getDZ()const{ return (dz); }
  double get_invdx()const{ return (invdx); }
  double get_invdy()const{ return (invdy); }
  double get_invdz()const{ return (invdz); }
  //
  // coordinate accessors
  //
  // calculated equivalents (preferred for accelerator?):
  double calcXN(int X)const{ return xStart+(X-1)*dx;}
  double calcYN(int Y)const{ return yStart+(Y-1)*dy;}
  double calcZN(int Z)const{ return zStart+(Z-1)*dz;}
  //const pfloat &get_pfloat_XN(int X)const{ return pfloat_node_xcoord[X];}
  //const pfloat &get_pfloat_YN(int Y)const{ return pfloat_node_ycoord[Y];}
  //const pfloat &get_pfloat_ZN(int Z)const{ return pfloat_node_zcoord[Z];}
  const double &getXN(int X)const{ return node_xcoord[X];}
  const double &getYN(int Y)const{ return node_ycoord[Y];}
  const double &getZN(int Z)const{ return node_zcoord[Z];}
  const double &getXC(int X)const{ return center_xcoord[X];}
  const double &getYC(int Y)const{ return center_ycoord[Y];}
  const double &getZC(int Z)const{ return center_zcoord[Z];}
  //
  // The following could be eliminated in favor of the previous
  // unless we truly anticipate generalizing to a deformed
  // logically cartesian mesh.  See issue #40.
  //
  const double &getXN(int X, int Y, int Z)const{ return node_xcoord[X];}
  const double &getYN(int X, int Y, int Z)const{ return node_ycoord[Y];}
  const double &getZN(int X, int Y, int Z)const{ return node_zcoord[Z];}
  const double &getXC(int X, int Y, int Z)const{ return center_xcoord[X];}
  const double &getYC(int X, int Y, int Z)const{ return center_ycoord[Y];}
  const double &getZC(int X, int Y, int Z)const{ return center_zcoord[Z];}
  //
  double getXstart()const{ return (xStart); }
  double getXend()const{ return (xEnd); }
  double getYstart()const{ return (yStart); }
  double getYend()const{ return (yEnd); } 
  double getZstart()const{ return (zStart); }
  double getZend()const{ return (zEnd); }
  double getVOL()const{ return (VOL); }
  double getInvVOL()const{ return (invVOL); }

  // inline methods to calculate mesh cell and weights.
  static void get_weights(double weights[8],
    double w0x, double w0y, double w0z,
    double w1x, double w1y, double w1z)
  {
    // which of the following is faster?
    //
    // this:
    //
    //const double weight00 = w0x*w0y;
    //const double weight01 = w0x*w1y;
    //const double weight10 = w1x*w0y;
    //const double weight11 = w1x*w1y;
    //weights[0] = weight00*w0z; // weight000
    //weights[1] = weight00*w1z; // weight001
    //weights[2] = weight01*w0z; // weight010
    //weights[3] = weight01*w1z; // weight011
    //weights[4] = weight10*w0z; // weight100
    //weights[5] = weight10*w1z; // weight101
    //weights[6] = weight11*w0z; // weight110
    //weights[7] = weight11*w1z; // weight111
    //
    // or this:
    //
    weights[0] = w0x*w0y*w0z; // weight000
    weights[1] = w0x*w0y*w1z; // weight001
    weights[2] = w0x*w1y*w0z; // weight010
    weights[3] = w0x*w1y*w1z; // weight011
    weights[4] = w1x*w0y*w0z; // weight100
    weights[5] = w1x*w0y*w1z; // weight101
    weights[6] = w1x*w1y*w0z; // weight110
    weights[7] = w1x*w1y*w1z; // weight111
  }
  void get_cell_coordinates(
    int& cx, int& cy, int& cz,
    double xpos, double ypos, double zpos)const
  {
      // xStart marks start of domain excluding ghosts
      const double rel_xpos = xpos - xStart;
      const double rel_ypos = ypos - yStart;
      const double rel_zpos = zpos - zStart;
      // cell position minus 1 (due to ghost cells)
      const double cxm1_pos = rel_xpos * invdx;
      const double cym1_pos = rel_ypos * invdy;
      const double czm1_pos = rel_zpos * invdz;
      cx = 1 + int(floor(cxm1_pos));
      cy = 1 + int(floor(cym1_pos));
      cz = 1 + int(floor(czm1_pos));
  }
  void make_grid_position_safe(double& cx_pos, double& cy_pos, double& cz_pos)const
  {
    // if the position is outside the domain, then map
    // it to the edge of the guarded subdomain
    //
    if (cx_pos < epsilon) cx_pos = epsilon;
    if (cy_pos < epsilon) cy_pos = epsilon;
    if (cz_pos < epsilon) cz_pos = epsilon;
    if (cx_pos > nxc_minus_epsilon) cx_pos = nxc_minus_epsilon;
    if (cy_pos > nyc_minus_epsilon) cy_pos = nyc_minus_epsilon;
    if (cz_pos > nzc_minus_epsilon) cz_pos = nzc_minus_epsilon;
  }
  void assert_cell_coordinates_safe(int cx, int cy, int cz)const
  {
    assert(cx >= 0);
    assert(cy >= 0);
    assert(cz >= 0);
    assert(cx < nxc);
    assert(cy < nyc);
    assert(cz < nzc);
  }
  void make_cell_coordinates_safe(int& cx, int& cy, int& cz)const
  {
    // if the cell is outside the domain, then treat it as
    // in the nearest ghost cell.
    //
    if (cx < 0) cx = 0;
    if (cy < 0) cy = 0;
    if (cz < 0) cz = 0;
    if (cx > cxlast) cx = cxlast; //nxc-1;
    if (cy > cylast) cy = cylast; //nyc-1;
    if (cz > czlast) cz = czlast; //nzc-1;
  }
  void get_safe_cell_coordinates(
    int& cx, int& cy, int& cz,
    double x, double y, double z)const
  {
    get_cell_coordinates(cx,cy,cz,x,y,z);
    make_cell_coordinates_safe(cx,cy,cz);
  }
  void get_safe_cell_and_weights(
    double xpos, double ypos, double zpos,
    int &cx, int& cy, int& cz,
    double weights[8])const
  {
    //convert_xpos_to_cxpos(xpos,ypos,zpos,cx_pos,cy_pos,cz_pos);
    // gxStart marks start of guarded domain (including ghosts)
    const double rel_xpos = xpos - xStart_g;
    const double rel_ypos = ypos - yStart_g;
    const double rel_zpos = zpos - zStart_g;
    // cell position (in guarded array)
    double cx_pos = rel_xpos * invdx;
    double cy_pos = rel_ypos * invdy;
    double cz_pos = rel_zpos * invdz;
    //
    if(suppress_runaway_particle_instability)
      make_grid_position_safe(cx_pos,cy_pos,cz_pos);
    //
    cx = int(floor(cx_pos));
    cy = int(floor(cy_pos));
    cz = int(floor(cz_pos));
    // this was the old algorithm.
    if(!suppress_runaway_particle_instability)
      make_cell_coordinates_safe(cx,cy,cz);
    assert_cell_coordinates_safe(cx,cy,cz);
  
    // fraction of distance from the left
    const double w0x = cx_pos - cx;
    const double w0y = cy_pos - cy;
    const double w0z = cz_pos - cz;
    // fraction of the distance from the right of the cell
    const double w1x = 1.-w0x;
    const double w1y = 1.-w0y;
    const double w1z = 1.-w0z;

    get_weights(weights, w0x, w0y, w0z, w1x, w1y, w1z);
  }
  void get_safe_cell_and_weights(double xpos[3], int cx[3], double weights[8])const
  {
    get_safe_cell_and_weights(xpos[0],xpos[1],xpos[2],cx[0],cx[1],cx[2],weights);
  }
};

typedef Grid3DCU Grid;

#endif
