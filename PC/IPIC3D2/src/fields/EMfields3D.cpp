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

#include <mpi.h>
#include "ipichdf5.h"
#include "EMfields3D.h"
#include "Collective.h"
#include "Basic.h"
#include "Com3DNonblk.h"
#include "VCtopology3D.h"
#include "Grid3DCU.h"
#include "CG.h"
#include "GMRES.h"
#include "Particles3Dcomm.h"
#include "TimeTasks.h"
#include "Moments.h"
#include "Parameters.h"
#include "ompdefs.h"
#include "debug.h"
#include "string.h"
#include "mic_particles.h"
#include "ipicmath.h" // for roundup_to_multiple
#include "Alloc.h"
#include "asserts.h"
#ifdef BATSRUS
#include "LinearSolver.h"
#endif
#ifndef NO_HDF5
#endif

#include <iostream>
#include <stdio.h>
//#include <sstream>

#define NG_F 1

using std::cout;
using std::endl;
using namespace iPic3D;

/*! constructor */
//
// We rely on the following rule from the C++ standard, section 12.6.2.5:
//
//   nonstatic data members shall be initialized in the order
//   they were declared in the class definition
//
// in particular, nxc, nyc, nzc and nxn, nyn, nzn are assumed
// initialized when subsequently used.
//
EMfields3D::EMfields3D(Collective * col, Grid * grid, VirtualTopology3D *vct) : 
  _col(*col),
  _grid(*grid),
  _vct(*vct),
  nxc(grid->getNXC()),
  nxn(grid->getNXN()),
  nyc(grid->getNYC()),
  nyn(grid->getNYN()),
  nzc(grid->getNZC()),
  nzn(grid->getNZN()),
  dx(grid->getDX()),
  dy(grid->getDY()),
  dz(grid->getDZ()),
  invVOL(grid->getInvVOL()),
  xStart(grid->getXstart()),
  xEnd(grid->getXend()),
  yStart(grid->getYstart()),
  yEnd(grid->getYend()),
  zStart(grid->getZstart()),
  zEnd(grid->getZend()),
  Lx(col->getLx()),
  Ly(col->getLy()),
  Lz(col->getLz()),
  ns(col->getNs()),
  c(col->getC()),
  dt(col->getDt()),
  th(col->getTh()),
  ue0(col->getU0(0)),
  ve0(col->getV0(0)),
  we0(col->getW0(0)),
  x_center(col->getx_center()),
  y_center(col->gety_center()),
  z_center(col->getz_center()),
  L_square(col->getL_square()),
  delt (c*th*dt), // declared after these
  //
  // array allocation: nodes
  //
  fieldForPcls  (nxn, nyn, nzn, 2*DFIELD_3or4),
  Ex   (nxn, nyn, nzn),
  Ey   (nxn, nyn, nzn),
  Ez   (nxn, nyn, nzn),
  Exth (nxn, nyn, nzn),
  Eyth (nxn, nyn, nzn),
  Ezth (nxn, nyn, nzn),
  Exthp(nxn, nyn, nzn),
  Eythp(nxn, nyn, nzn),
  Ezthp(nxn, nyn, nzn),
  Bxn  (nxn, nyn, nzn),
  Byn  (nxn, nyn, nzn),
  Bzn  (nxn, nyn, nzn),
  rhon (nxn, nyn, nzn),
  Jx   (nxn, nyn, nzn),
  Jy   (nxn, nyn, nzn),
  Jz   (nxn, nyn, nzn),
  Jxh  (nxn, nyn, nzn),
  Jyh  (nxn, nyn, nzn),
  Jzh  (nxn, nyn, nzn),
  divEn(nxn, nyn, nzn),  
  energyError_G(nxn,nyn,nzn),
  //
  // species-specific quantities
  //
  rhons (ns, nxn, nyn, nzn),
  rhocs (ns, nxc, nyc, nzc),
  Jxs   (ns, nxn, nyn, nzn),
  Jys   (ns, nxn, nyn, nzn),
  Jzs   (ns, nxn, nyn, nzn),
  Jxsh  (ns, nxn, nyn, nzn),
  Jysh  (ns, nxn, nyn, nzn),
  Jzsh  (ns, nxn, nyn, nzn),
  pXXsn (ns, nxn, nyn, nzn),
  pXYsn (ns, nxn, nyn, nzn),
  pXZsn (ns, nxn, nyn, nzn),
  pYYsn (ns, nxn, nyn, nzn),
  pYZsn (ns, nxn, nyn, nzn),
  pZZsn (ns, nxn, nyn, nzn),
  // array allocation: central points 
  //
  PHI  (nxc, nyc, nzc),
  Bxc  (nxc, nyc, nzc),
  Byc  (nxc, nyc, nzc),
  Bzc  (nxc, nyc, nzc),
  rhoc (nxc, nyc, nzc),
  rhoh (nxc, nyc, nzc),
  divEc(nxc,nyc,nzc),
  
  // temporary arrays
  //
  tempXC (nxc, nyc, nzc),
  tempYC (nxc, nyc, nzc),
  tempZC (nxc, nyc, nzc),
  //
  tempXN (nxn, nyn, nzn),
  tempYN (nxn, nyn, nzn),
  tempZN (nxn, nyn, nzn),
  tempC  (nxc, nyc, nzc),
  tempX  (nxn, nyn, nzn),
  tempY  (nxn, nyn, nzn),
  tempZ  (nxn, nyn, nzn),
  temp2X (nxn, nyn, nzn),
  temp2Y (nxn, nyn, nzn),
  temp2Z (nxn, nyn, nzn),
  imageX (nxn, nyn, nzn),
  imageY (nxn, nyn, nzn),
  imageZ (nxn, nyn, nzn),
  Dx (nxn, nyn, nzn),
  Dy (nxn, nyn, nzn),
  Dz (nxn, nyn, nzn),
  vectX (nxn, nyn, nzn),
  vectY (nxn, nyn, nzn),
  vectZ (nxn, nyn, nzn),
  divC  (nxc, nyc, nzc),
  //arr (nxc-2,nyc-2,nzc-2),
  // B_ext and J_ext should not be allocated unless used.
#ifdef BATSRUS
  fluidBxc (nxc,nyc,nzc),
  fluidByc (nxc,nyc,nzc),
  fluidBzc (nxc,nyc,nzc),
  fluidBxn (nxn,nyn,nzn),
  fluidByn (nxn,nyn,nzn),
  fluidBzn (nxn,nyn,nzn),
  fluidExc (nxc,nyc,nzc),
  fluidEyc (nxc,nyc,nzc),
  fluidEzc (nxc,nyc,nzc),
  fluidEx  (nxn,nyn,nzn),
  fluidEy  (nxn,nyn,nzn),
  fluidEz  (nxn,nyn,nzn),
  
  p0XXsn (ns, nxn, nyn, nzn),
  p0YYsn (ns, nxn, nyn, nzn),
  p0ZZsn (ns, nxn, nyn, nzn),
  p0XYsn (ns, nxn, nyn, nzn),
  p0XZsn (ns, nxn, nyn, nzn),
  p0YZsn (ns, nxn, nyn, nzn),

#endif
  Bx_ext(nxn,nyn,nzn),
  By_ext(nxn,nyn,nzn),
  Bz_ext(nxn,nyn,nzn),
  Bx_tot(nxn,nyn,nzn),
  By_tot(nxn,nyn,nzn),
  Bz_tot(nxn,nyn,nzn),
  Jx_ext(nxn,nyn,nzn),
  Jy_ext(nxn,nyn,nzn),
  Jz_ext(nxn,nyn,nzn)
{
  // Define mass matrix.
  M_GII = newArr5(double, nxn, nyn, nzn, ngp, n9);
  M_CI  = newArr4(double, nxc, nyc, nzc, ngp);

  // External imposed fields
  //
  B1x = col->getB1x();
  B1y = col->getB1y();
  B1z = col->getB1z();
  //if(B1x!=0. || B1y !=0. || B1z!=0.)
  //{
  //  eprintf("This functionality has not yet been implemented");
  //}
  Bx_ext.setall(0.);
  By_ext.setall(0.);
  Bz_ext.setall(0.);
  Bx_tot.setall(0.);
  By_tot.setall(0.);
  Bz_tot.setall(0.);
  //
  PoissonCorrection = false;
  if (col->getPoissonCorrection()=="yes"){
	  PoissonCorrection = true;
	  PoissonCorrectionCycle = col->getPoissonCorrectionCycle();
  }
  CGtol = col->getCGtol();
  GMREStol = col->getGMREStol();
  nGMRESRestart = col->get_nGMRESRestart();
  qom = new double[ns];
  for (int i = 0; i < ns; i++)
    qom[i] = col->getQOM(i);
  // boundary conditions: PHI and EM fields
  bcPHIfaceXright = col->getBcPHIfaceXright();
  bcPHIfaceXleft  = col->getBcPHIfaceXleft();
  bcPHIfaceYright = col->getBcPHIfaceYright();
  bcPHIfaceYleft  = col->getBcPHIfaceYleft();
  bcPHIfaceZright = col->getBcPHIfaceZright();
  bcPHIfaceZleft  = col->getBcPHIfaceZleft();

  bcEMfaceXright = col->getBcEMfaceXright();
  bcEMfaceXleft = col->getBcEMfaceXleft();
  bcEMfaceYright = col->getBcEMfaceYright();
  bcEMfaceYleft = col->getBcEMfaceYleft();
  bcEMfaceZright = col->getBcEMfaceZright();
  bcEMfaceZleft = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  delta = col->getDelta();
  Smooth = col->getSmooth();
  SmoothNiter = col->getSmoothNiter();
  // get the density background for the gem Challange
  rhoINIT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i = 0; i < ns; i++) {
    rhoINIT[i] = col->getRHOinit(i);
    if ((fabs(col->getW0(i)) != 0) || (fabs(col->getU0(i)) != 0)) // GEM and LHDI
      DriftSpecies[i] = true;
    else
      DriftSpecies[i] = false;
  }
  /*! parameters for GEM challenge */
  FourPI = 16 * atan(1.0);
  /*! Restart */
  restart1 = col->getRestart_status();



  /** index of nodes containing unknows in implicit solver*/
  inminsolve = 1;
  inmaxsolve = nxn-2;
  jnminsolve = 1;
  jnmaxsolve = nyn-2;
  knminsolve = 1;
  knmaxsolve = nzn-2;
  n3SolveNode = 3
    *       (inmaxsolve - inminsolve + 1)
    *       (jnmaxsolve - jnminsolve + 1)
    *       (knmaxsolve - knminsolve + 1);


  /* Index of cells solved by field solver. */
  icMinSolve = 1; icMaxSolve = nxc - 2;
  jcMinSolve = 1; jcMaxSolve = nyc - 2;
  kcMinSolve = 1; kcMaxSolve = nzc - 2;   
  nSolveCell  = (icMaxSolve - icMinSolve + 1)
    *           (jcMaxSolve - jcMinSolve + 1)
    *           (kcMaxSolve - kcMinSolve + 1);
  
  if(Parameters::get_VECTORIZE_MOMENTS())
  {
    // In this case particles are sorted
    // and there is no need for each thread
    // to sum moments in a separate array.
    sizeMomentsArray = 1;
  }
  else
  {
    sizeMomentsArray = omp_get_max_threads();
  }
  moments10Array = new Moments10*[sizeMomentsArray];
  moments13Array = new Moments13*[sizeMomentsArray];
  for(int i=0;i<sizeMomentsArray;i++)
  {
    moments10Array[i] = new Moments10(nxn,nyn,nzn);
    moments13Array[i] = new Moments13(nxn,nyn,nzn);
  }
        
  //----Define MPI Derived Data types for Halo Exchange.   Begin------------------
  // Three Grids:
  //  1) center array (nxc*nyc*nzc)
  //  2) node array (nxc+1)*(nyc+1)*(nzc+1)
  //  3) the auxiliary array for ghost cell swaping (nxc+2)*(nyc+2)*(nzc+2)
  nGridMpiData = 3; 
  int nSend;
  nSend = 1;      mpiDataIdx[nSend] = 0;
  nSend = ngp;    mpiDataIdx[nSend] = 1;
  nSend = ngp*n9; mpiDataIdx[nSend] = 2; 
  
  const int nKey = mpiDataIdx.size()*nGridMpiData;
  
  yzFacetype =new MPI_Datatype[nKey];  
  xzFacetype =new MPI_Datatype[nKey]; 
  xyFacetype =new MPI_Datatype[nKey]; 
  xEdgetype  =new MPI_Datatype[nKey]; 
  yEdgetype  =new MPI_Datatype[nKey]; 
  zEdgetype  =new MPI_Datatype[nKey]; 
  xEdgetype2 =new MPI_Datatype[nKey]; 
  yEdgetype2 =new MPI_Datatype[nKey]; 
  zEdgetype2 =new MPI_Datatype[nKey]; 
  cornertype =new MPI_Datatype[nKey]; 

  int nDouble, idx;
  for(int iLayAdd = 0; iLayAdd < nGridMpiData; iLayAdd++ ){
    const int nx = nxc + iLayAdd;
    const int ny = nyc + iLayAdd;
    const int nz = nzc + iLayAdd;
    map<int, int>::const_iterator mpiIter = mpiDataIdx.begin();
    while(mpiIter !=mpiDataIdx.end()){
      nDouble = mpiIter->first; 
      idx = (mpiIter->second) + iLayAdd*nGridMpiData;
      mpiIter++;

      //For face exchange on X dir
      MPI_Type_vector((ny-2),(nz-2)*nDouble,nz*nDouble, MPI_DOUBLE, &yzFacetype[idx]);
      MPI_Type_commit(&yzFacetype[idx]);
    
      //For face exchange on Y dir
      MPI_Type_hvector((nx-2),(nz-2)*nDouble,(nDouble*nz*ny*sizeof(double)), MPI_DOUBLE, &xzFacetype[idx]);
      MPI_Type_commit(&xzFacetype[idx]);

      MPI_Type_vector((ny-2), nDouble, nz*nDouble, MPI_DOUBLE, &yEdgetype[idx]);
      MPI_Type_commit(&yEdgetype[idx]);
    
      //For face exchangeg on Z dir
      MPI_Type_hvector((nx-2), 1, (nDouble*nz*ny*sizeof(double)), yEdgetype[idx], &xyFacetype[idx]);
      MPI_Type_commit(&xyFacetype[idx]);
    
      //2 yEdgeType can be merged into one message
      MPI_Type_hvector(2, 1, nDouble*(nz-1)*sizeof(double), yEdgetype[idx], &yEdgetype2[idx]);
      MPI_Type_commit(&yEdgetype2[idx]);
    
      MPI_Type_contiguous((nz-2)*nDouble,MPI_DOUBLE, &zEdgetype[idx]);
      MPI_Type_commit(&zEdgetype[idx]);
    
      MPI_Type_hvector(2, (nz-2)*nDouble,nDouble*(nx-1)*(ny*nz)*sizeof(double), MPI_DOUBLE, &zEdgetype2[idx]);
      MPI_Type_commit(&zEdgetype2[idx]);
    
      MPI_Type_vector((nx-2), nDouble, nDouble*ny*nz, MPI_DOUBLE, &xEdgetype[idx]);
      MPI_Type_commit(&xEdgetype[idx]);
      MPI_Type_hvector(2, 1, nDouble*(ny-1)*nz*sizeof(double), xEdgetype[idx], &xEdgetype2[idx]);
      MPI_Type_commit(&xEdgetype2[idx]);
    
      //corner used to communicate in x direction
      int blocklength[]={nDouble,nDouble,nDouble,nDouble};
      int displacements[]={0,(nz-1)*nDouble,(ny-1)*nz*nDouble,(ny*nz-1)*nDouble};
      MPI_Type_indexed(4, blocklength, displacements, MPI_DOUBLE, &cornertype[idx]);
      MPI_Type_commit(&cornertype[idx]);
   
    }//mpiIter

  }
  //----Define MPI Derived Data types for Halo Exchange   End------------------

    
    if (col->getWriteMethod() == "pvtk" || col->getWriteMethod() == "nbcvtk"){
    	//test Endian
    	int TestEndian = 1;
    	lEndFlag =*(char*)&TestEndian;

    	 //create process file view
		  int  size[3], subsize[3], start[3];

		  //3D subarray - reverse X, Z
		  subsize[0] = nzc-2; subsize[1] = nyc-2; subsize[2] = nxc-2;
		  size[0] = (nzc-2)*vct->getZLEN();size[1] = (nyc-2)*vct->getYLEN();size[2] = (nxc-2)*vct->getXLEN();
		  start[0]= vct->getCoordinates(2)*subsize[0];
		  start[1]= vct->getCoordinates(1)*subsize[1];
		  start[2]= vct->getCoordinates(0)*subsize[2];

		  MPI_Type_contiguous(3,MPI_FLOAT, &xyzcomp);
		  MPI_Type_commit(&xyzcomp);

		  MPI_Type_create_subarray(3, size, subsize, start,MPI_ORDER_C, xyzcomp, &procviewXYZ);
		  MPI_Type_commit(&procviewXYZ);

		  MPI_Type_create_subarray(3, size, subsize, start,MPI_ORDER_C, MPI_FLOAT, &procview);
		  MPI_Type_commit(&procview);

		  subsize[0] = nxc-2; subsize[1] =nyc-2; subsize[2] = nzc-2;
		  size[0]    = nxc;	  size[1] 	 =nyc;	 size[2] 	= nzc;
		  start[0]	 = 1;	  start[1]	 =1;	 start[2]	= 1;
		  MPI_Type_create_subarray(3, size, subsize, start,MPI_ORDER_C, MPI_FLOAT, &ghosttype);
		  MPI_Type_commit(&ghosttype);

    }

}

void EMfields3D::freeDataType(){
  const int nKey = mpiDataIdx.size()*nGridMpiData;
  for(int iCount = 0; iCount<nKey; iCount++){
    MPI_Type_free(&yzFacetype[iCount]);
    MPI_Type_free(&xzFacetype[iCount]);
    MPI_Type_free(&xyFacetype[iCount]);
    MPI_Type_free(&xEdgetype[iCount]);
    MPI_Type_free(&yEdgetype[iCount]);
    MPI_Type_free(&zEdgetype[iCount]);
    MPI_Type_free(&xEdgetype2[iCount]);
    MPI_Type_free(&yEdgetype2[iCount]);
    MPI_Type_free(&zEdgetype2[iCount]);
    MPI_Type_free(&cornertype[iCount]);    
  }

    delete [] yzFacetype;
    delete [] xzFacetype;
    delete [] xyFacetype;
    delete [] xEdgetype;
    delete [] yEdgetype;
    delete [] zEdgetype;
    delete [] xEdgetype2;
    delete [] yEdgetype2;
    delete [] zEdgetype2;
    delete [] cornertype;    
}


// This was Particles3Dcomm::interpP2G()
void EMfields3D::sumMomentsOld(const Particles3Dcomm& pcls)
{
  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  double const*const x = pcls.getXall();
  double const*const y = pcls.getYall();
  double const*const z = pcls.getZall();
  double const*const u = pcls.getUall();
  double const*const v = pcls.getVall();
  double const*const w = pcls.getWall();
  double const*const q = pcls.getQall();
  //
  const int is = pcls.get_species_num();

  const int nop = pcls.getNOP();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP
  TimeTasks timeTasksAcc;
  #pragma omp parallel private(timeTasks)
  {
    int thread_num = omp_get_thread_num();
    Moments10& speciesMoments10 = fetch_moments10Array(thread_num);
    speciesMoments10.set_to_zero();
    arr4_double moments = speciesMoments10.fetch_arr();
    // The following loop is expensive, so it is wise to assume that the
    // compiler is stupid.  Therefore we should on the one hand
    // expand things out and on the other hand avoid repeating computations.
    #pragma omp for
    for (int i = 0; i < nop; i++)
    {
      // compute the quadratic moments of velocity
      //
      const double ui=u[i];
      const double vi=v[i];
      const double wi=w[i];
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[10];
      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;

      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      const double xi0   = x[i] - grid->getXN(ix-1);
      const double eta0  = y[i] - grid->getYN(iy-1);
      const double zeta0 = z[i] - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - x[i];
      const double eta1  = grid->getYN(iy) - y[i];
      const double zeta1 = grid->getZN(iz) - z[i];
      const double qi = q[i];
      const double weight000 = qi * xi0 * eta0 * zeta0 * invVOL;
      const double weight001 = qi * xi0 * eta0 * zeta1 * invVOL;
      const double weight010 = qi * xi0 * eta1 * zeta0 * invVOL;
      const double weight011 = qi * xi0 * eta1 * zeta1 * invVOL;
      const double weight100 = qi * xi1 * eta0 * zeta0 * invVOL;
      const double weight101 = qi * xi1 * eta0 * zeta1 * invVOL;
      const double weight110 = qi * xi1 * eta1 * zeta0 * invVOL;
      const double weight111 = qi * xi1 * eta1 * zeta1 * invVOL;
      double weights[8];
      weights[0] = weight000;
      weights[1] = weight001;
      weights[2] = weight010;
      weights[3] = weight011;
      weights[4] = weight100;
      weights[5] = weight101;
      weights[6] = weight110;
      weights[7] = weight111;

      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        momentsArray[0] = moments[ix  ][iy  ][iz  ]; // moments000 
        momentsArray[1] = moments[ix  ][iy  ][iz-1]; // moments001 
        momentsArray[2] = moments[ix  ][iy-1][iz  ]; // moments010 
        momentsArray[3] = moments[ix  ][iy-1][iz-1]; // moments011 
        momentsArray[4] = moments[ix-1][iy  ][iz  ]; // moments100 
        momentsArray[5] = moments[ix-1][iy  ][iz-1]; // moments101 
        momentsArray[6] = moments[ix-1][iy-1][iz  ]; // moments110 
        momentsArray[7] = moments[ix-1][iy-1][iz-1]; // moments111 

        for(int m=0; m<10; m++)
        for(int c=0; c<8; c++)
        {
          momentsArray[c][m] += velmoments[m]*weights[c];
        }
      }
    }
    

    // reduction
    

    // reduce arrays
    {
      #pragma omp critical (reduceMoment0)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { rhons[is][i][j][k] += invVOL*moments[i][j][k][0]; }}
      #pragma omp critical (reduceMoment1)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1]; }}
      #pragma omp critical (reduceMoment2)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jys  [is][i][j][k] += invVOL*moments[i][j][k][2]; }}
      #pragma omp critical (reduceMoment3)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3]; }}
      #pragma omp critical (reduceMoment4)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4]; }}
      #pragma omp critical (reduceMoment5)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5]; }}
      #pragma omp critical (reduceMoment6)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6]; }}
      #pragma omp critical (reduceMoment7)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7]; }}
      #pragma omp critical (reduceMoment8)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8]; }}
      #pragma omp critical (reduceMoment9)
      for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
        { pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9]; }}
    }
    
    #pragma omp critical
    timeTasksAcc += timeTasks;
  }
  // reset timeTasks to be its average value for all threads
  timeTasksAcc /= omp_get_max_threads();
  timeTasks = timeTasksAcc;
  communicateGhostP2G(is);
}
//
// Create a vectorized version of this moment accumulator as follows.
//
// A. Moment accumulation
//
// Case 1: Assuming AoS particle layout and using intrinsics vectorization:
//   Process P:=N/4 particles at a time:
//   1. gather position coordinates from P particles and
//      generate Px8 array of weights and P cell indices.
//   2. for each particle, add 10x8 array of moment-weight
//      products to appropriate cell accumulator.
//   Each cell now has a 10x8 array of node-destined moments.
//   (See sumMoments_AoS_intr().)
// Case 2: Assuming SoA particle layout and using trivial vectorization:
//   Process N:=sizeof(vector_unit)/sizeof(double) particles at a time:
//   1. for pcl=1:N: (3) positions -> (8) weights, cell_index
//   2. For each of 10 moments:
//      a. for pcl=1:N: (<=2 of 3) charge velocities -> (1) moment
//      b. for pcl=1:N: (1) moment, (8) weights -> (8) node-destined moments
//      c. transpose 8xN array of node-destined moments to Nx8 array 
//      d. foreach pcl: add node-distined moments to cell of cell_index
//   Each cell now has a 10x8 array of node-destined moments.
//
//   If particles are sorted by mesh cell, then all moments are destined
//   for the same node; in this case, we can simply accumulate an 8xN
//   array of node-destined moments in each mesh cell and at the end
//   gather these moments at the nodes; to help performance and
//   code reuse, we will in each cell first transpose the 8xN array
//   of node-destined moments to an Nx8 array.
//
// B: Moment reduction
//
//   Gather the information from cells to nodes:
//   3. [foreach cell transpose node-destined moments:
//      10x8 -> 8x10, or rather (8+2)x8 -> 8x8 + 8x2]
//   4. at each node gather moments from cells.
//   5. [transpose moments at nodes if step 3 was done.]
//
//   We will likely omit steps 3 and 5; they could help to optimize,
//   but even without these steps, step 4 is not expected to dominate.
//
// Compare the vectorization notes at the top of mover_PC().
//
// This was Particles3Dcomm::interpP2G()
void EMfields3D::sumMoments(const Particles3Dcomm* part)
{
  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP
  #pragma omp parallel
  {
  for (int i = 0; i < ns; i++)
  {
    const Particles3Dcomm& pcls = part[i];
    assert_eq(pcls.get_particleType(), ParticleType::SoA);
    const int is = pcls.get_species_num();
    assert_eq(i,is);

    double const*const x = pcls.getXall();
    double const*const y = pcls.getYall();
    double const*const z = pcls.getZall();
    double const*const u = pcls.getUall();
    double const*const v = pcls.getVall();
    double const*const w = pcls.getWall();
    double const*const q = pcls.getQall();

    const int nop = pcls.getNOP();

    int thread_num = omp_get_thread_num();
    
    Moments10& speciesMoments10 = fetch_moments10Array(thread_num);
    arr4_double moments = speciesMoments10.fetch_arr();
    //
    // moments.setmode(ompmode::mine);
    // moments.setall(0.);
    // 
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    //
    // This barrier is not needed
    #pragma omp barrier
    // The following loop is expensive, so it is wise to assume that the
    // compiler is stupid.  Therefore we should on the one hand
    // expand things out and on the other hand avoid repeating computations.
    #pragma omp for // used nowait with the old way
    for (int i = 0; i < nop; i++)
    {
      // compute the quadratic moments of velocity
      //
      const double ui=u[i];
      const double vi=v[i];
      const double wi=w[i];
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[10];
      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;

      //
      // compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      const double xi0   = x[i] - grid->getXN(ix-1);
      const double eta0  = y[i] - grid->getYN(iy-1);
      const double zeta0 = z[i] - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - x[i];
      const double eta1  = grid->getYN(iy) - y[i];
      const double zeta1 = grid->getZN(iz) - z[i];
      const double qi = q[i];
      const double invVOLqi = invVOL*qi;
      const double weight0 = invVOLqi * xi0;
      const double weight1 = invVOLqi * xi1;
      const double weight00 = weight0*eta0;
      const double weight01 = weight0*eta1;
      const double weight10 = weight1*eta0;
      const double weight11 = weight1*eta1;
      double weights[8];
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

      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        arr2_double_fetch moments00 = moments[ix  ][iy  ];
        arr2_double_fetch moments01 = moments[ix  ][iy-1];
        arr2_double_fetch moments10 = moments[ix-1][iy  ];
        arr2_double_fetch moments11 = moments[ix-1][iy-1];
        momentsArray[0] = moments00[iz  ]; // moments000 
        momentsArray[1] = moments00[iz-1]; // moments001 
        momentsArray[2] = moments01[iz  ]; // moments010 
        momentsArray[3] = moments01[iz-1]; // moments011 
        momentsArray[4] = moments10[iz  ]; // moments100 
        momentsArray[5] = moments10[iz-1]; // moments101 
        momentsArray[6] = moments11[iz  ]; // moments110 
        momentsArray[7] = moments11[iz-1]; // moments111 

        for(int m=0; m<10; m++)
        for(int c=0; c<8; c++)
        {
          momentsArray[c][m] += velmoments[m]*weights[c];
        }
      }
    }
    

    // reduction
    

    // reduce moments in parallel
    //
    for(int thread_num=0;thread_num<get_sizeMomentsArray();thread_num++)
    {
      arr4_double moments = fetch_moments10Array(thread_num).fetch_arr();
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++)
      for(int j=0;j<nyn;j++)
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] += invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] += invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9];
      }
    }
    //
    // This was the old way of reducing;
    // did not scale well to large number of threads
    //{
    //  #pragma omp critical (reduceMoment0)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { rhons[is][i][j][k] += invVOL*moments[i][j][k][0]; }}
    //  #pragma omp critical (reduceMoment1)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1]; }}
    //  #pragma omp critical (reduceMoment2)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jys  [is][i][j][k] += invVOL*moments[i][j][k][2]; }}
    //  #pragma omp critical (reduceMoment3)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3]; }}
    //  #pragma omp critical (reduceMoment4)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4]; }}
    //  #pragma omp critical (reduceMoment5)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5]; }}
    //  #pragma omp critical (reduceMoment6)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6]; }}
    //  #pragma omp critical (reduceMoment7)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7]; }}
    //  #pragma omp critical (reduceMoment8)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8]; }}
    //  #pragma omp critical (reduceMoment9)
    //  for(int i=0;i<nxn;i++){for(int j=0;j<nyn;j++) for(int k=0;k<nzn;k++)
    //    { pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9]; }}
    //}
    
    // uncomment this and remove the loop below
    // when we change to use asynchronous communication.
    // communicateGhostP2G(is, vct);
  }
  }
  for (int i = 0; i < ns; i++)
  {
    communicateGhostP2G(i);
  }
}

void EMfields3D::sumMoments_AoS(const Particles3Dcomm* part, bool doCalcMomentsOnly)
{
  const Collective *col = &get_col();
  bool useAccurateJ = col->getuseAccurateJ() && (!doCalcMomentsOnly);
  const VirtualTopology3D *vct = &get_vct();

  bool doTestEMWave = false;

  #ifdef BATSRUS
  doTestEMWave = col->getdoTestEMWave();
  #endif
  
  if(doTestEMWave){
    for (int species_idx = 0; species_idx < ns; species_idx++)
      {
	const Particles3Dcomm& pcls = part[species_idx];
	assert_eq(pcls.get_particleType(), ParticleType::AoS);
	const int is = pcls.get_species_num();
	// Use eqValue() to set these arrays to zero!! -- Yuxi
	for(int i=0;i<nxn;i++)
	  for(int j=0;j<nyn;j++)
	    for(int k=0;k<nzn;k++)
	      {
		rhons[is][i][j][k] = 0; 
		Jxs  [is][i][j][k] = 0; 
		Jys  [is][i][j][k] = 0; 
		Jzs  [is][i][j][k] = 0; 
		pXXsn[is][i][j][k] = 0; 
		pXYsn[is][i][j][k] = 0; 
		pXZsn[is][i][j][k] = 0; 
		pYYsn[is][i][j][k] = 0; 
		pYZsn[is][i][j][k] = 0; 
		pZZsn[is][i][j][k] = 0;
		Jxh[i][j][k] = 0;
		Jyh[i][j][k] = 0;
		Jzh[i][j][k] = 0;
	      }
      }
    return;
  }


  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // To make memory use scale to a large number of threads, we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  //#ifdef _OPENMP

#ifdef BATSRUS
  string nameFunc="sumMoments_AoS";
  bool doTestFunc, doTestCell;
  ofstream * outfile;
  doTestFunc = do_test_func(outfile,nameFunc);
  if(doTestFunc){
    int ispecies=0;
    const Particles3Dcomm& pcls = part[ispecies];
    const int nop = pcls.getNOP();

    for(int pidx=0; pidx < nop; pidx++){
      const SpeciesParticle& pcl = pcls.get_pcl(pidx);
      const double ui=pcl.get_u();
      const double vi=pcl.get_v();
      const double wi=pcl.get_w();
      const double xi=pcl.get_x();
      const double yi=pcl.get_y();
      const double zi=pcl.get_z();

      const int ix = 2 + int (floor((pcl.get_x() - xstart) * inv_dx));
      const int iy = 2 + int (floor((pcl.get_y() - ystart) * inv_dy));
      const int iz = 2 + int (floor((pcl.get_z() - zstart) * inv_dz));

      doTestCell = do_test_cell(outfile,ix-1,iy-1,iz-1);
      if(doTestCell){
	(*outfile)<<xi<<"\t"<<yi<<"\t"<<zi<<"\t"
		  <<ui<<"\t"<<vi<<"\t"<<wi<<"\t"<<"\n";
      }
    }
  }
#endif


  // Used when 'useAccurateJ' is true.
  double weights_III[2][2][2];
  double alpha_DD[nDimMax][nDimMax];
  set_fieldForPcls();
  const_arr4_pfloat fieldForPcls = get_fieldForPcls();    
  if(useAccurateJ){
  for(int i1 = 0; i1 < nxn; i1++ )
    for(int j1 = 0; j1 < nyn; j1++)
      for(int k1 = 0; k1 < nzn; k1++)
	for(int gp = 0; gp < ngp; gp++)
	  for(int iR = 0; iR < n9; iR++)
	    M_GII[i1][j1][k1][gp][iR] = 0;
  }


  
#pragma omp parallel
  {
  for (int species_idx = 0; species_idx < ns; species_idx++)
  {
    const Particles3Dcomm& pcls = part[species_idx];
    assert_eq(pcls.get_particleType(), ParticleType::AoS);
    const int is = pcls.get_species_num();
    assert_eq(species_idx,is);

    const int nop = pcls.getNOP();

    int thread_num = omp_get_thread_num();
    Moments13& speciesMoments13 = fetch_moments13Array(thread_num);
    arr4_double moments = speciesMoments13.fetch_arr();
    //
    // moments.setmode(ompmode::mine);
    // moments.setall(0.);
    // 
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    //

    int dChunk = nxc-2; 
    for(int iChunk = 1; iChunk <= nxc - 2; iChunk+=dChunk){
#pragma omp barrier
#pragma omp for
    for (int pidx = 0; pidx < nop; pidx++)
    {
      const SpeciesParticle& pcl = pcls.get_pcl(pidx);
      // compute the quadratic moments of velocity
      //
      const double xp = pcl.get_x();
      const int iCell = 1 + int (floor((xp - xstart) * inv_dx));
      if(iCell>= iChunk && iCell<iChunk+dChunk){
      
      const double yp = pcl.get_y();
      const double zp = pcl.get_z();
      const double ui=pcl.get_u();
      const double vi=pcl.get_v();
      const double wi=pcl.get_w();
      const double uui=ui*ui;
      const double uvi=ui*vi;
      const double uwi=ui*wi;
      const double vvi=vi*vi;
      const double vwi=vi*wi;
      const double wwi=wi*wi;
      double velmoments[13];

      velmoments[0] = 1.;
      velmoments[1] = ui;
      velmoments[2] = vi;
      velmoments[3] = wi;
      velmoments[4] = uui;
      velmoments[5] = uvi;
      velmoments[6] = uwi;
      velmoments[7] = vvi;
      velmoments[8] = vwi;
      velmoments[9] = wwi;
      velmoments[10]= 0 ;
      velmoments[11]= 0 ;
      velmoments[12]= 0 ;

      
      //
      // Compute the weights to distribute the moments
      //
      const int ix = 2 + int (floor((xp - xstart) * inv_dx));
      const int iy = 2 + int (floor((yp - ystart) * inv_dy));
      const int iz = 2 + int (floor((zp - zstart) * inv_dz));
      const double xi0   = xp - grid->getXN(ix-1);
      const double eta0  = yp - grid->getYN(iy-1);
      const double zeta0 = zp - grid->getZN(iz-1);
      const double xi1   = grid->getXN(ix) - xp;
      const double eta1  = grid->getYN(iy) - yp;
      const double zeta1 = grid->getZN(iz) - zp;
      const double qi = pcl.get_q();
      const double invVOLqi = invVOL*qi;
      const double weight0 = invVOLqi * xi0;
      const double weight1 = invVOLqi * xi1;
      const double weight00 = weight0*eta0;
      const double weight01 = weight0*eta1;
      const double weight10 = weight1*eta0;
      const double weight11 = weight1*eta1;
      double weights[8];
      weights[0] = weight00*zeta0; // weight000
      weights[1] = weight00*zeta1; // weight001
      weights[2] = weight01*zeta0; // weight010
      weights[3] = weight01*zeta1; // weight011
      weights[4] = weight10*zeta0; // weight100
      weights[5] = weight10*zeta1; // weight101
      weights[6] = weight11*zeta0; // weight110
      weights[7] = weight11*zeta1; // weight111


      if(useAccurateJ){
	// The following calculations should be replaced by 'weights' -- Yuxi
	weights_III[0][0][0] = xi1 * eta1 * zeta1 * invVOL;
	weights_III[0][0][1] = xi1 * eta1 * zeta0 * invVOL;
	weights_III[0][1][0] = xi1 * eta0 * zeta1 * invVOL; 
	weights_III[0][1][1] = xi1 * eta0 * zeta0 * invVOL;
	weights_III[1][0][0] = xi0 * eta1 * zeta1 * invVOL;
	weights_III[1][0][1] = xi0 * eta1 * zeta0 * invVOL;
	weights_III[1][1][0] = xi0 * eta0 * zeta1 * invVOL;
	weights_III[1][1][1] = xi0 * eta0 * zeta0 * invVOL;	
	
	const double dto2 = .5 * dt, qdto2mc = qom[is] * dto2 / c;
	// Change the name 'weights'!!! It is also used outside the 'useAccurateJ' scope!!! -- Yuxi
	double weights[8] ALLOC_ALIGNED;
	int cx,cy,cz;

	// It seems the particle is always 'safe' at current stage? -Yuxi
	grid->get_safe_cell_and_weights(xp,yp,zp, cx,cy,cz,weights);

	const double* field_components[8] ALLOC_ALIGNED;
	get_field_components_for_cell(field_components,fieldForPcls,cx,cy,cz);

	double sampled_field[8] ALLOC_ALIGNED;
	for(int i=0;i<8;i++) sampled_field[i]=0;
	double& Bxl=sampled_field[0];
	double& Byl=sampled_field[1];
	double& Bzl=sampled_field[2];
	
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

	const pfloat udotOm = ui * Omx + vi * Omy + wi * Omz;
	// solve the velocity equation
	velmoments[10] = (ui + (vi * Omz - wi * Omy + udotOm * Omx)) * denom;
	velmoments[11] = (vi + (wi * Omx - ui * Omz + udotOm * Omy)) * denom;
	velmoments[12] = (wi + (ui * Omy - vi * Omx + udotOm * Omz)) * denom;

	// Assemble rotation matrix alpha_DD = Rot_DD*beta*q

	const double c0 = denom*invVOLqi*qdto2mc;
	alpha_DD[x_][x_] = (    1 + Omx*Omx)*c0;
	alpha_DD[x_][y_] = (  Omz + Omx*Omy)*c0;
	alpha_DD[x_][z_] = ( -Omy + Omx*Omz)*c0;
	alpha_DD[y_][x_] = ( -Omz + Omx*Omy)*c0;
	alpha_DD[y_][y_] = (    1 + Omy*Omy)*c0;
	alpha_DD[y_][z_] = (  Omx + Omy*Omz)*c0;
	alpha_DD[z_][x_] = (  Omy + Omx*Omz)*c0;
	alpha_DD[z_][y_] = ( -Omx + Omy*Omz)*c0;
	alpha_DD[z_][z_] = (    1 + Omz*Omz)*c0;

	int count = 0; 	
	int ip, jp, kp; // Indexes for g'
	double wg, wgp, weight; //W_pg, W_pg'

	const int iMin = ix - 1;
	const int jMin = iy - 1;
	const int kMin = iz - 1;
	const int iMax = ix + 1;
	const int jMax = iy + 1;
	const int kMax = iz + 1;

	int gp = -1; 
	for(int i1 = iMin; i1 < iMax; i1++ )
	  for(int j1 = jMin; j1 < jMax; j1++)
	    for(int k1 = kMin; k1 < kMax; k1++){
	      wg = weights_III[i1-ix+1][j1-iy+1][k1-iz+1];
	      for(int i2 = iMin; i2 < iMax; i2++){
		ip = i2 - i1 + 1;
		if(ip > 0){		
		  for(int j2 = jMin; j2 < jMax; j2++){
		    jp = j2 - j1 + 1;
		     for(int k2 = kMin; k2 < kMax; k2++){
		      kp = k2 - k1 + 1;
		      weight = wg*weights_III[i2-ix+1][j2-iy+1][k2-iz+1];
		      gp = ip*9 + jp*3 + kp;

		      double *M_I = M_GII[i1][j1][k1][gp];
		      
		      M_I[0] += alpha_DD[0][0]*weight;
		      M_I[1] += alpha_DD[0][1]*weight;
		      M_I[2] += alpha_DD[0][2]*weight;                                    
		      M_I[3] += alpha_DD[1][0]*weight;
		      M_I[4] += alpha_DD[1][1]*weight;
		      M_I[5] += alpha_DD[1][2]*weight;                                    
		      M_I[6] += alpha_DD[2][0]*weight;
		      M_I[7] += alpha_DD[2][1]*weight;
		      M_I[8] += alpha_DD[2][2]*weight;                                    
		     } // k2		  	      		  
		  } // j2                            
		}// if (ip > 0)
	      }// i2
	    } //k1

      } // useAccurateJ
      
      
      // add particle to moments
      {
        arr1_double_fetch momentsArray[8];
        arr2_double_fetch moments00 = moments[ix  ][iy  ];
        arr2_double_fetch moments01 = moments[ix  ][iy-1];
        arr2_double_fetch moments10 = moments[ix-1][iy  ];
        arr2_double_fetch moments11 = moments[ix-1][iy-1];
        momentsArray[0] = moments00[iz  ]; // moments000 
        momentsArray[1] = moments00[iz-1]; // moments001 
        momentsArray[2] = moments01[iz  ]; // moments010 
        momentsArray[3] = moments01[iz-1]; // moments011 
        momentsArray[4] = moments10[iz  ]; // moments100 
        momentsArray[5] = moments10[iz-1]; // moments101 
        momentsArray[6] = moments11[iz  ]; // moments110 
        momentsArray[7] = moments11[iz-1]; // moments111 

        for(int m=0; m<13; m++)
	  for(int c=0; c<8; c++)
	    {
	      momentsArray[c][m] += velmoments[m]*weights[c];
	    }
      }
      }
    } // pidx 
    } // iChunk
   

    if(useAccurateJ){
      int gps, gpr; // gp_send, gp_receive
      for(int i1 = 1; i1 < nxn-1; i1++)
	for(int j1 = 1; j1 < nyn-1; j1++)
	  for(int k1 = 1; k1 < nzn-1; k1++){
	    const int ip = 2, ipr = 0; 	    
	    const int ir = i1 + ip - 1;
	      for(int jp = 0; jp < 3; jp++){
		const int jr = j1 + jp - 1;
		const int jpr = 2 - jp;
		for(int kp = 0; kp < 3; kp++){
		    const int kr = k1 + kp - 1;
		    const int kpr = 2 - kp;
		    gpr = jpr*3 + kpr; // gpr = ipr*9 + jpr*3 + kpr
		    gps = 18 + jp*3 + kp; // gps = ip*9+jp*3+kp
		    
		    M_GII[ir][jr][kr][gpr][0] = M_GII[i1][j1][k1][gps][0];
		    M_GII[ir][jr][kr][gpr][1] = M_GII[i1][j1][k1][gps][1];
		    M_GII[ir][jr][kr][gpr][2] = M_GII[i1][j1][k1][gps][2];
		    M_GII[ir][jr][kr][gpr][3] = M_GII[i1][j1][k1][gps][3];
		    M_GII[ir][jr][kr][gpr][4] = M_GII[i1][j1][k1][gps][4];
		    M_GII[ir][jr][kr][gpr][5] = M_GII[i1][j1][k1][gps][5];
		    M_GII[ir][jr][kr][gpr][6] = M_GII[i1][j1][k1][gps][6];
		    M_GII[ir][jr][kr][gpr][7] = M_GII[i1][j1][k1][gps][7];
		    M_GII[ir][jr][kr][gpr][8] = M_GII[i1][j1][k1][gps][8];		    
		}
	      }
	      } 
    }
    
    // reduction
   

    // reduce moments in parallel
    //
    for(int thread_num=0;thread_num<get_sizeMomentsArray();thread_num++)
    {
      arr4_double moments = fetch_moments13Array(thread_num).fetch_arr();
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++)
      for(int j=0;j<nyn;j++)
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] += invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] += invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] += invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] += invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] += invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] += invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] += invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] += invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] += invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] += invVOL*moments[i][j][k][9];
        Jxsh [is][i][j][k] += invVOL*moments[i][j][k][10];
        Jysh [is][i][j][k] += invVOL*moments[i][j][k][11];
        Jzsh [is][i][j][k] += invVOL*moments[i][j][k][12];
      }
    }
    
  }
  }

  for (int i = 0; i < ns; i++)
  {
    communicateGhostP2G(i);
  }
  
  if(col->getuseAccurateJ()){
    communicateInterp(nxn, nyn, nzn, ngp, n9, M_GII, vct, this);
    communicateNode_P(nxn, nyn, nzn, ngp, n9, M_GII, vct, this);
  }
}



void EMfields3D::calc_cell_center_density(const Particles3Dcomm* part,
					  bool doCalcDensityOnly){
  const Collective *col = &get_col();
  if(col->get_DoCalcRhocDirectly()){
    // Calculate cell center densities rhocs from particles
    sum_cell_center_density(part, doCalcDensityOnly);
  }else{
    // Interpolate rhocs from rhons
    interpDensitiesN2C();	   
  }
}

void EMfields3D::sum_cell_center_density(const Particles3Dcomm* part, bool doCalcDensityOnly){
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();

  const Grid *grid = &get_grid();

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;

  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();

  
  const int iSpeciesCorrect = col->get_iSpeciesLightest();
  const string divECleanType = col->get_divECleanType();
  bool doCalcMatrix;
  bool doCorrectWeight = (divECleanType=="weight");

  for(int iSpecies = 0; iSpecies<ns; iSpecies++)
    for(int i = 0; i<nxc; i++)
      for(int j = 0; j<nyc; j++)
    	for(int k = 0; k<nzc; k++){
    	  rhocs[iSpecies][i][j][k] = 0; 
    	}

  if(!doCalcDensityOnly)
    for(int i = 0; i<nxc; i++)
      for(int j = 0; j<nyc; j++)
	for(int k = 0; k<nzc; k++)
	  for(int i4 = 0; i4<ngp; i4++)
	    M_CI[i][j][k][i4] = 0;   

#pragma omp parallel
  {
    double weights_III[2][2][2];  
    double weights_IIID[2][2][2][nDimMax];
    for (int species_idx = 0; species_idx < ns; species_idx++){
      const Particles3Dcomm& pcls = part[species_idx];
      assert_eq(pcls.get_particleType(), ParticleType::AoS);
      const int is = pcls.get_species_num();
      assert_eq(species_idx,is);

      const int nop = pcls.getNOP();

      doCalcMatrix = (!doCalcDensityOnly && 
		      (is==iSpeciesCorrect || 
		       divECleanType=="position_all") &&
		      divECleanType.substr(0,15) != "weight_estimate");
      
#pragma omp barrier
#pragma omp for
      for (int pidx = 0; pidx < nop; pidx++)
	{
	  const SpeciesParticle& pcl = pcls.get_pcl(pidx);
	  const double xp = pcl.get_x();      
	  const double yp = pcl.get_y();
	  const double zp = pcl.get_z();

	  const int ix = int (floor((xp - xstart) * inv_dx + 0.5));
	  const int iy = int (floor((yp - ystart) * inv_dy + 0.5));
	  const int iz = int (floor((zp - zstart) * inv_dz + 0.5));	    

	  const double xi0   = xp - grid->getXC(ix);
	  const double eta0  = yp - grid->getYC(iy);
	  const double zeta0 = zp - grid->getZC(iz);
	  const double xi1   = grid->getXC(ix+1) - xp;
	  const double eta1  = grid->getYC(iy+1) - yp;
	  const double zeta1 = grid->getZC(iz+1) - zp;

	  const double qi = pcl.get_q();
	  const double invVOLqi = invVOL*qi;

	  const double weight0 = invVOL*xi0;
	  const double weight1 = invVOL*xi1;
	  const double weight00 = weight0*eta0;
	  const double weight01 = weight0*eta1;
	  const double weight10 = weight1*eta0;
	  const double weight11 = weight1*eta1;	  

	  
          weights_III[1][1][1]= weight00 * zeta0; 
          weights_III[1][1][0]= weight00 * zeta1; 
          weights_III[1][0][1]= weight01 * zeta0; 
          weights_III[1][0][0]= weight01 * zeta1;
          weights_III[0][1][1]= weight10 * zeta0;
          weights_III[0][1][0]= weight10 * zeta1;
          weights_III[0][0][1]= weight11 * zeta0;
          weights_III[0][0][0]= weight11 * zeta1;

	  for(int ii = 0; ii<2; ii++)
	    for(int jj = 0; jj<2; jj++)
	      for(int kk=0; kk<2; kk++)
		rhocs[is][ix+ii][iy+jj][iz+kk] += weights_III[ii][jj][kk]*qi;

	  if(doCalcMatrix){
	    const int nPower = col->get_nPowerWeight();
	    
	    const int iMin = ix;
	    const int jMin = iy;
	    const int kMin = iz;
	    const int iMax = ix + 1;
	    const int jMax = iy + 1;
	    const int kMax = iz + 1;
	    

	    if(doCorrectWeight){
	      // The weights_III calculated for cell center density 
	      // interpolation is w_pc, which is the same for the weight
	      // correction with minimizing sum(0.5*(r_p-1)^2*q^nPower), where
	      // nPower is 2; 
	      // When nPower is 0, the interpolation function is w_pc*q.
	      // qi*w_pc. 

	      double coef=1.0;
	      if(nPower==1) coef = fabs(qi);
	      if(nPower==0) coef = qi*qi;
		
	      double weight, wg;
	      int gp, ip, jp, kp;
	      for (int i1 = iMin; i1 <= iMax; i1++)
		for (int j1 = jMin; j1 <= jMax; j1++)
		  for (int k1 = kMin; k1 <= kMax; k1++) {
		    wg = invVOL*weights_III[i1 - iMin][j1 - jMin][k1 - kMin];
		    for (int i2 = iMin; i2 <= iMax; i2++)
		      for (int j2 = jMin; j2 <= jMax; j2++)
			for (int k2 = kMin; k2 <= kMax; k2++) {
			  weight =
			    wg * weights_III[i2 - iMin][j2 - jMin][k2 - kMin];
			  ip = i2 - i1 + 1;
			  jp = j2 - j1 + 1;
			  kp = k2 - k1 + 1;
			  gp = ip * 9 + jp * 3 + kp;			
			  M_CI[i1][j1][k1][gp] += weight*coef;			
			}// k2
		  }// k1

	    }else{
	      // Particle position correction. 

	      weights_IIID[1][1][1][x_] = eta0*zeta0*invVOL;
	      weights_IIID[1][1][1][y_] = xi0*zeta0*invVOL;
	      weights_IIID[1][1][1][z_] = xi0*eta0*invVOL;
	      
	      // xi0*eta0*zeta1*invVOL;
	      weights_IIID[1][1][0][x_] = eta0*zeta1*invVOL; 
	      weights_IIID[1][1][0][y_] = xi0*zeta1*invVOL; 
	      weights_IIID[1][1][0][z_] = -xi0*eta0*invVOL; 

	      // xi0*eta1*zeta0*invVOL;
	      weights_IIID[1][0][1][x_] = eta1*zeta0*invVOL;
	      weights_IIID[1][0][1][y_] = -xi0*zeta0*invVOL;
	      weights_IIID[1][0][1][z_] = xi0*eta1*invVOL;

	      // xi0*eta1*zeta1*invVOL;
	      weights_IIID[1][0][0][x_] = eta1*zeta1*invVOL;
	      weights_IIID[1][0][0][y_] = -xi0*zeta1*invVOL;
	      weights_IIID[1][0][0][z_] = -xi0*eta1*invVOL;

	      // xi1*eta0*zeta0*invVOL;
	      weights_IIID[0][1][1][x_] = -eta0*zeta0*invVOL;
	      weights_IIID[0][1][1][y_] = xi1*zeta0*invVOL;
	      weights_IIID[0][1][1][z_] = xi1*eta0*invVOL;


	      // xi1*eta0*zeta1*invVOL;
	      weights_IIID[0][1][0][x_] = -eta0*zeta1*invVOL;
	      weights_IIID[0][1][0][y_] = xi1*zeta1*invVOL;
	      weights_IIID[0][1][0][z_] = -xi1*eta0*invVOL;
	      

	      // xi1*eta1*zeta0*invVOL;
	      weights_IIID[0][0][1][x_] = -eta1*zeta0*invVOL;
	      weights_IIID[0][0][1][y_] = -xi1*zeta0*invVOL;
	      weights_IIID[0][0][1][z_] = xi1*eta1*invVOL;

	      // xi1*eta1*zeta1*invVOL;
	      weights_IIID[0][0][0][x_] = -eta1*zeta1*invVOL;
	      weights_IIID[0][0][0][y_] = -xi1*zeta1*invVOL;
	      weights_IIID[0][0][0][z_] = -xi1*eta1*invVOL;	     
	      

	      double coef=1;
	      if(nPower==1) coef = fabs(qi);
	      if(nPower==0) coef = qi*qi;	    	    	 	    

	      double weight, wg_D[3];
	      int gp, ip, jp, kp;
	      for (int i1 = iMin; i1 <= iMax; i1++)
		for (int j1 = jMin; j1 <= jMax; j1++)
		  for (int k1 = kMin; k1 <= kMax; k1++) {
		    double *M_I = M_CI[i1][j1][k1];
		    for(int iDim = 0; iDim < nDimMax; iDim++)
		      wg_D[iDim] = invVOL*
			weights_IIID[i1 - iMin][j1 - jMin][k1 - kMin][iDim];
		    
		    for (int i2 = iMin; i2 <= iMax; i2++)
		      for (int j2 = jMin; j2 <= jMax; j2++)
			for (int k2 = kMin; k2 <= kMax; k2++) {			 

			  const double tmp0 = wg_D[0] *
			    weights_IIID[i2-iMin][j2-jMin][k2-kMin][0];
			  const double tmp1 = wg_D[1] *
			    weights_IIID[i2-iMin][j2-jMin][k2-kMin][1];
			  const double tmp2 = wg_D[2] *
			    weights_IIID[i2-iMin][j2-jMin][k2-kMin][2];
			  weight = tmp0 + tmp1 + tmp2; 
			  

			  // for(int iDim = 0; iDim<nDimMax; iDim++){
			  //   weight += wg_D[iDim] * 
			  //     weights_IIID[i2-iMin][j2-jMin][k2-kMin][iDim];
			  // }


			  ip = i2 - i1 + 1;
			  jp = j2 - j1 + 1;
			  kp = k2 - k1 + 1;
			  gp = ip * 9 + jp * 3 + kp;			
			  M_I[gp] += weight*coef;			
			  
			  // cout<<"i1 = "<<i1<<" j1 = "<<j1<<" k1 = "<<k1<<" gp = "<<gp
			  //     <<" weight = "<<weight<<" M = "<<M_CI[i1][j1][k1][gp]
			  //     <<endl;

			}// k2
		  }// k1

	    }

	  } 

        }// pidx
      


      
      
      
      for(int i=0;i<nxc;i++)
	for(int j=0;j<nyc;j++)
	  for(int k=0;k<nzc;k++){
	    rhocs[is][i][j][k] *=invVOL;
	  }



    }// iSpecies    
  }


  // for(int i =1; i<nxc-1; i++)
  //   for(int j = 1; j<nyc-1; j++)
  //     for(int k = 1; k<nzc-1; k++){
  // 	for(int gp=0; gp<27; gp++){
  // 	  double m = M_CI[i][j][k][gp];
	  
  // 	  if(m !=0) cout<<"0 i = "<<i<<" j = "<<j<<" k = "<<k<<" gp = "<<gp
  // 			<<" m = "<<m<<endl;
	  
  // 	}

	
  //     }




  if(!doCalcDensityOnly){
    communicateSwapAddGhostCell(nxc, nyc, nzc, ngp, M_CI, vct, this);
    communicateCenter_P(nxc, nyc, nzc, ngp, M_CI, vct, this);
  }



  // for(int i =1; i<nxc-1; i++)
  //   for(int j = 1; j<nyc-1; j++)
  //     for(int k = 1; k<nzc-1; k++){
  // 	for(int gp=0; gp<27; gp++){
  // 	  double m = M_CI[i][j][k][gp];
	  
  // 	  if(m !=0) cout<<"1 i = "<<i<<" j = "<<j<<" k = "<<k<<" gp = "<<gp
  // 			<<" m = "<<m<<endl;
	  
  // 	}

	
  //     }





  for(int is = 0; is < ns; is++){
    double ***moment = convert_to_arr3(rhocs[is]);
    communicateSwapAddGhostCell(nxc, nyc, nzc, moment, vct, this);
    communicateCenterBC_P(nxc, nyc, nzc, moment, 2, 2, 2, 2, 2, 2, vct, this);

    const double signq = qom[is]/(fabs(qom[is]));

#ifdef BATSRUS
    // Apply 1st-order accurate boundary conditions. 
    if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) 
      for (int i = 0; i < nyc; i++)
	for (int k = 0; k < nzc; k++){
	  // The right hand side returns the values at nodes. Only 1st-order accurate. 
	  rhocs[is][0][i][k] = signq*col->getPICRhoNum(0,i,k,is);
	  rhocs[is][1][i][k] = signq*col->getPICRhoNum(1,i,k,is);
	}

    
    if (vct->getYleft_neighbor_P() == MPI_PROC_NULL)
      for (int i = 0; i < nxc; i++)
	for (int k = 0; k < nzc; k++){
	  rhocs[is][i][0][k] = signq*col->getPICRhoNum(i,0,k,is);
	  rhocs[is][i][1][k] = signq*col->getPICRhoNum(i,1,k,is);
	}
    
    if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) 
      for (int i = 0; i < nxc; i++)
	for (int j = 0; j < nyc; j++){
	  rhocs[is][i][j][0] = signq*col->getPICRhoNum(i,j,0,is);
	  rhocs[is][i][j][1] = signq*col->getPICRhoNum(i,j,1,is);
	}


    if (vct->getXright_neighbor_P() == MPI_PROC_NULL) 
      for (int i = 0; i < nyc; i++)
	for (int k = 0; k < nzc; k++){
	  rhocs[is][nxc - 1][i][k] = signq*col->getPICRhoNum(nxc-1,i,k,is);
	  rhocs[is][nxc - 2][i][k] = signq*col->getPICRhoNum(nxc-2,i,k,is);
	}


    if (vct->getYright_neighbor_P() == MPI_PROC_NULL) 
      for (int i = 0; i < nxc; i++)
	for (int k = 0; k < nzc; k++){
	  rhocs[is][i][nyc - 1][k] = signq*col->getPICRhoNum(i,nyc-1,k,is);
	  rhocs[is][i][nyc - 2][k] = signq*col->getPICRhoNum(i,nyc-2,k,is);
	}


    if (vct->getZright_neighbor_P() == MPI_PROC_NULL) 
      for (int i = 0; i < nxc; i++)
	for (int j = 0; j < nyc; j++){
	  rhocs[is][i][j][nzc - 1] = signq*col->getPICRhoNum(i,j,nzc-1,is);
	  rhocs[is][i][j][nzc - 2] = signq*col->getPICRhoNum(i,j,nzc-2,is);
	}    
#endif
  }


  eqValue(0,rhon,nxc,nyc,nzc);
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxc; i++)
      for (register int j = 0; j < nyc; j++)
        for (register int k = 0; k < nzc; k++){
          rhoc[i][j][k] += rhocs[is][i][j][k];
	}

}


#ifdef __MIC__
// add moment weights to all ten moments for the cell of the particle
// (assumes that particle data is aligned with cache boundary and
// begins with the velocity components)
inline void addto_cell_moments(
  F64vec8* cell_moments,
  F64vec8 weights,
  F64vec8 vel)
{
  // broadcast particle velocities
  const F64vec8 u = F64vec8(vel[0]);
  const F64vec8 v = F64vec8(vel[1]);
  const F64vec8 w = F64vec8(vel[2]);
  // construct kronecker product of moments and weights
  const F64vec8 u_weights = u*weights;
  const F64vec8 v_weights = v*weights;
  const F64vec8 w_weights = w*weights;
  const F64vec8 uu_weights = u*u_weights;
  const F64vec8 uv_weights = u*v_weights;
  const F64vec8 uw_weights = u*w_weights;
  const F64vec8 vv_weights = v*v_weights;
  const F64vec8 vw_weights = v*w_weights;
  const F64vec8 ww_weights = w*w_weights;
  // add moment weights to accumulated moment weights in mesh mesh
  cell_moments[0] += weights;
  cell_moments[1] += u_weights;
  cell_moments[2] += v_weights;
  cell_moments[3] += w_weights;
  cell_moments[4] += uu_weights;
  cell_moments[5] += uv_weights;
  cell_moments[6] += uw_weights;
  cell_moments[7] += vv_weights;
  cell_moments[8] += vw_weights;
  cell_moments[9] += ww_weights;
}
#endif // __MIC__

// sum moments of AoS using MIC intrinsics
// 
// We could rewrite this without intrinsics also.  The core idea
// of this algorithm is that instead of scattering the data of
// each particle to its nodes, in each cell we accumulate the
// data that would be scattered and then scatter it at the end.
// By waiting to scatter, with each particle we work with an
// aligned 10x8 matrix rather than a 8x10 matrix, which means
// that for each particle we make 10 vector stores rather than
// 8*2=16 or 8*3=24 vector stores (for unaligned data).  This
// also avoids the expense of computing node indices for each
// particle.
//
// 1. compute vector of 8 weights using position
// 2. form kronecker product of weights with moments
//    by scaling the weights by each velocity moment;
//    add each to accumulated weights for this cell
// 3. after sum is complete, transpose weight-moment
//    product in each cell and distribute to its 8 nodes.
//    An optimized way:
//    A. transpose the first 8 weighted moments with fast 8x8
//       matrix transpose.
//    B. transpose 2x8 matrix of the last two weighted moments
//       and then use 8 masked vector adds to accumulate
//       to weights at nodes.
//    But the optimized way might be overkill since distributing
//    the sums from the cells to the nodes should not dominate
//    if the number of particles per mesh cell is large;
//    if the number of particles per mesh cell is small,
//    then a fully vectorized moment sum is hard to justify anyway.
//
// See notes at the top of sumMoments().
//
void EMfields3D::sumMoments_AoS_intr(const Particles3Dcomm* part)
{
#ifndef __MIC__
  eprintf("not implemented");
#else
  const Grid *grid = &get_grid();

  // define global parameters
  //
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  // Here and below x stands for all 3 physical position coordinates
  const F64vec8 dx_inv = make_F64vec8(inv_dx, inv_dy, inv_dz);
  // starting physical position of proper subdomain ("pdom", without ghosts)
  const F64vec8 pdom_xlow = make_F64vec8(xstart,ystart, zstart);
  //
  // X = canonical coordinates.
  //
  // starting position of cell in lower corner
  // of proper subdomain (without ghosts);
  // probably this is an integer value, but we won't rely on it.
  const F64vec8 pdom_Xlow = dx_inv*pdom_xlow;
  // g = including ghosts
  // starting position of cell in low corner
  const F64vec8 gdom_Xlow = pdom_Xlow - F64vec8(1.);
  // starting position of cell in high corner of physical domain
  // in canonical coordinates of ghost domain
  const F64vec8 nXcm1 = make_F64vec8(nxc-1,nyc-1,nzc-1);

  // allocate memory per mesh cell for accumulating moments
  //
  const int num_threads = omp_get_max_threads();
  array4<F64vec8>* cell_moments_per_thr
    = (array4<F64vec8>*) malloc(num_threads*sizeof(array4<F64vec8>));
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // use placement new to allocate array to accumulate moments for thread
    new(&cell_moments_per_thr[thread_num]) array4<F64vec8>(nxc,nyc,nzc,10);
  }
  //
  // allocate memory per mesh node for accumulating moments
  //
  array3<F64vec8>* node_moments_first8_per_thr
    = (array3<F64vec8>*) malloc(num_threads*sizeof(array3<F64vec8>));
  array4<double>* node_moments_last2_per_thr
    = (array4<double>*) malloc(num_threads*sizeof(array4<double>));
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // use placement new to allocate array to accumulate moments for thread
    new(&node_moments_first8_per_thr[thread_num]) array3<F64vec8>(nxn,nyn,nzn);
    new(&node_moments_last2_per_thr[thread_num]) array4<double>(nxn,nyn,nzn,2);
  }

  // The moments of a particle must be distributed to the 8 nodes of the cell
  // in proportion to the weight of each node.
  //
  // Refer to the kronecker product of weights and moments as
  // "weighted moments" or "moment weights".
  //
  // Each thread accumulates moment weights in cells.
  //
  // Because particles are not assumed to be sorted by mesh cell,
  // we have to wait until all particles have been processed
  // before we transpose moment weights to weighted moments;
  // the memory that we must allocate to sum moments is thus
  // num_thread*8 times as much as if particles were pre-sorted
  // by mesh cell (and num_threads times as much as if particles
  // were sorted by thread subdomain).
  //
  #pragma omp parallel
  {
    // array4<F64vec8> cell_moments(nxc,nyc,nzc,10);
    const int this_thread = omp_get_thread_num();
    assert_lt(this_thread,num_threads);
    array4<F64vec8>& cell_moments = cell_moments_per_thr[this_thread];

    for (int species_idx = 0; species_idx < ns; species_idx++)
    {
      const Particles3Dcomm& pcls = part[species_idx];
      assert_eq(pcls.get_particleType(), ParticleType::AoS);
      const int is = pcls.get_species_num();
      assert_eq(species_idx,is);

      // moments.setmode(ompmode::mine);
      // moments.setall(0.);
      // 
      F64vec8 *cell_moments1d = &cell_moments[0][0][0][0];
      int moments1dsize = cell_moments.get_size();
      for(int i=0; i<moments1dsize; i++) cell_moments1d[i]=F64vec8(0.);
      //
      // number or particles processed at a time
      const int num_pcls_per_loop = 2;
      const vector_SpeciesParticle& pcl_list = pcls.get_pcl_list();
      const int nop = pcl_list.size();
      // if the number of particles is odd, then make
      // sure that the data after the last particle
      // will not contribute to the moments.
      #pragma omp single // the implied omp barrier is needed
      {
        // make sure that we will not overrun the array
        assert_divides(num_pcls_per_loop,pcl_list.capacity());
        // round up number of particles
        int nop_rounded_up = roundup_to_multiple(nop,num_pcls_per_loop);
        for(int pidx=nop; pidx<nop_rounded_up; pidx++)
        {
          // (This is a benign violation of particle
          // encapsulation and requires a cast).
          SpeciesParticle& pcl = (SpeciesParticle&) pcl_list[pidx];
          pcl.set_to_zero();
        }
      }
      #pragma omp for
      for (int pidx = 0; pidx < nop; pidx+=2)
      {
        // cast particles as vectors
        // (assumes each particle exactly fits a cache line)
        const F64vec8& pcl0 = (const F64vec8&)pcl_list[pidx];
        const F64vec8& pcl1 = (const F64vec8&)pcl_list[pidx+1];
        // gather position data from particles
        // (assumes position vectors are in upper half)
        const F64vec8 xpos = cat_hgh_halves(pcl0,pcl1);

        // convert to canonical coordinates relative to subdomain with ghosts
        const F64vec8 gX = dx_inv*xpos - gdom_Xlow;
        F64vec8 cellXstart = floor(gX);
        // all particles at this point should be inside the
        // proper subdomain of this process, but maybe we
        // will need to enforce this because of inconsistency
        // of floating point arithmetic?
        //cellXstart = maximum(cellXstart,F64vec8(1.));
        //cellXstart = minimum(cellXstart,nXcm1);
        assert(!test_lt(cellXstart,F64vec8(1.)));
        assert(!test_gt(cellXstart,nXcm1));

        // get weights for field_components based on particle position
        //
        F64vec8 weights[2];
        const F64vec8 X = gX - cellXstart;
        construct_weights_for_2pcls(weights, X);

        // add scaled weights to all ten moments for the cell of each particle
        //
        // the cell that we will write to
        const I32vec16 cell = round_to_nearest(cellXstart);
        const int* c=(int*)&cell;
        F64vec8* cell_moments0 = &cell_moments[c[0]][c[1]][c[2]][0];
        F64vec8* cell_moments1 = &cell_moments[c[4]][c[5]][c[6]][0];
        addto_cell_moments(cell_moments0, weights[0], pcl0);
        addto_cell_moments(cell_moments1, weights[1], pcl1);
      }
      if(!this_thread) timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION);

      // reduction
      if(!this_thread) timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION);

      // reduce moments in parallel
      //
      // this code currently makes no sense for multiple threads.
      assert_eq(num_threads,1);
      {
        // For each thread, distribute moments from cells to nodes
        // and then sum moments at each node over all threads.
        //
        // (Alternatively we could sum over all threads and then
        // distribute to nodes; this alternative would be preferable
        // for vectorization efficiency but more difficult to parallelize
        // across threads).

        // initialize moment accumulators
        //
        memset(&node_moments_first8_per_thr[this_thread][0][0][0],
          0, sizeof(F64vec8)*node_moments_first8_per_thr[0].get_size());
        memset(&node_moments_last2_per_thr[this_thread][0][0][0][0],
          0, sizeof(double)*node_moments_last2_per_thr[0].get_size());

        // distribute moments from cells to nodes
        //
        #pragma omp for collapse(2)
        for(int cx=1;cx<nxc;cx++)
        for(int cy=1;cy<nyc;cy++)
        for(int cz=1;cz<nzc;cz++)
        {
          const int ix=cx+1;
          const int iy=cy+1;
          const int iz=cz+1;
          F64vec8* cell_mom = &cell_moments[cx][cy][cz][0];

          // scatter the cell's first 8 moments to its nodes
          // for each thread
          {
            F64vec8* cell_mom_first8 = cell_mom;
            // regard cell_mom_first8 as a pointer to 8x8 data and transpose
            transpose_8x8_double((double(*)[8]) cell_mom_first8);
            // scatter the moment vectors to the nodes
            array3<F64vec8>& node_moments_first8 = node_moments_first8_per_thr[this_thread];
            arr_fetch2(F64vec8) node_moments0 = node_moments_first8[ix];
            arr_fetch2(F64vec8) node_moments1 = node_moments_first8[cx];
            arr_fetch1(F64vec8) node_moments00 = node_moments0[iy];
            arr_fetch1(F64vec8) node_moments01 = node_moments0[cy];
            arr_fetch1(F64vec8) node_moments10 = node_moments1[iy];
            arr_fetch1(F64vec8) node_moments11 = node_moments1[cy];
            node_moments00[iz] += cell_mom_first8[0]; // node_moments_first8[ix][iy][iz]
            node_moments00[cz] += cell_mom_first8[1]; // node_moments_first8[ix][iy][cz]
            node_moments01[iz] += cell_mom_first8[2]; // node_moments_first8[ix][cy][iz]
            node_moments01[cz] += cell_mom_first8[3]; // node_moments_first8[ix][cy][cz]
            node_moments10[iz] += cell_mom_first8[4]; // node_moments_first8[cx][iy][iz]
            node_moments10[cz] += cell_mom_first8[5]; // node_moments_first8[cx][iy][cz]
            node_moments11[iz] += cell_mom_first8[6]; // node_moments_first8[cx][cy][iz]
            node_moments11[cz] += cell_mom_first8[7]; // node_moments_first8[cx][cy][cz]
          }

          // scatter the cell's last 2 moments to its nodes
          {
            array4<double>& node_moments_last2 = node_moments_last2_per_thr[this_thread];
            arr3_double_fetch node_moments0 = node_moments_last2[ix];
            arr3_double_fetch node_moments1 = node_moments_last2[cx];
            arr2_double_fetch node_moments00 = node_moments0[iy];
            arr2_double_fetch node_moments01 = node_moments0[cy];
            arr2_double_fetch node_moments10 = node_moments1[iy];
            arr2_double_fetch node_moments11 = node_moments1[cy];
            double* node_moments000 = node_moments00[iz];
            double* node_moments001 = node_moments00[cz];
            double* node_moments010 = node_moments01[iz];
            double* node_moments011 = node_moments01[cz];
            double* node_moments100 = node_moments10[iz];
            double* node_moments101 = node_moments10[cz];
            double* node_moments110 = node_moments11[iz];
            double* node_moments111 = node_moments11[cz];

            const F64vec8 mom8 = cell_mom[8];
            const F64vec8 mom9 = cell_mom[9];

            bool naive_last2 = true;
            if(naive_last2)
            {
              node_moments000[0] += mom8[0]; node_moments000[1] += mom9[0];
              node_moments001[0] += mom8[1]; node_moments001[1] += mom9[1];
              node_moments010[0] += mom8[2]; node_moments010[1] += mom9[2];
              node_moments011[0] += mom8[3]; node_moments011[1] += mom9[3];
              node_moments100[0] += mom8[4]; node_moments100[1] += mom9[4];
              node_moments101[0] += mom8[5]; node_moments101[1] += mom9[5];
              node_moments110[0] += mom8[6]; node_moments110[1] += mom9[6];
              node_moments111[0] += mom8[7]; node_moments111[1] += mom9[7];
            }
            else
            {
              // Let a=moment#8 and b=moment#9.
              // Number the nodes 0 through 7.
              //
              // This transpose changes data from the form
              //   [a0 a1 a2 a3 a4 a5 a6 a7]=mom8
              //   [b0 b1 b2 b3 b4 b5 b6 b7]=mom9
              // into the form
              //   [a0 b0 a2 b2 a4 b4 a6 b6]=out8
              //   [a1 b1 a3 b3 a5 b5 a7 b7]=out9
              F64vec8 out8, out9;
              trans2x2(out8, out9, mom8, mom9);

              // probably the compiler is not smart enough to recognize that
              // each line can be done with a single vector instruction:
              node_moments000[0] += out8[0]; node_moments000[1] += out8[1];
              node_moments001[0] += out9[0]; node_moments001[1] += out9[1];
              node_moments010[0] += out8[2]; node_moments010[1] += out8[3];
              node_moments011[0] += out9[2]; node_moments011[1] += out9[3];
              node_moments100[0] += out8[4]; node_moments100[1] += out8[5];
              node_moments101[0] += out9[4]; node_moments101[1] += out9[5];
              node_moments110[0] += out8[6]; node_moments110[1] += out8[7];
              node_moments111[0] += out9[6]; node_moments111[1] += out9[7];
            }
          }
        }

        // at each node add moments to moments of first thread
        //
        #pragma omp for collapse(2)
        for(int nx=1;nx<nxn;nx++)
        for(int ny=1;ny<nyn;ny++)
        {
          arr_fetch1(F64vec8) node_moments8_for_master
            = node_moments_first8_per_thr[0][nx][ny];
          arr_fetch2(double) node_moments2_for_master
            = node_moments_last2_per_thr[0][nx][ny];
          for(int thread_num=1;thread_num<num_threads;thread_num++)
          {
            arr_fetch1(F64vec8) node_moments8_for_thr
              = node_moments_first8_per_thr[thread_num][nx][ny];
            arr_fetch2(double) node_moments2_for_thr
              = node_moments_last2_per_thr[thread_num][nx][ny];
            for(int nz=1;nz<nzn;nz++)
            {
              node_moments8_for_master[nz] += node_moments8_for_thr[nz];
              node_moments2_for_master[nz][0] += node_moments2_for_thr[nz][0];
              node_moments2_for_master[nz][1] += node_moments2_for_thr[nz][1];
            }
          }
        }

        // transpose moments for field solver
        //
        #pragma omp for collapse(2)
        for(int nx=1;nx<nxn;nx++)
        for(int ny=1;ny<nyn;ny++)
        {
          arr_fetch1(F64vec8) node_moments8_for_master
            = node_moments_first8_per_thr[0][nx][ny];
          arr_fetch2(double) node_moments2_for_master
            = node_moments_last2_per_thr[0][nx][ny];
          arr_fetch1(double) rho_sxy = rhons[is][nx][ny];
          arr_fetch1(double) Jx__sxy = Jxs  [is][nx][ny];
          arr_fetch1(double) Jy__sxy = Jys  [is][nx][ny];
          arr_fetch1(double) Jz__sxy = Jzs  [is][nx][ny];
          arr_fetch1(double) pXX_sxy = pXXsn[is][nx][ny];
          arr_fetch1(double) pXY_sxy = pXYsn[is][nx][ny];
          arr_fetch1(double) pXZ_sxy = pXZsn[is][nx][ny];
          arr_fetch1(double) pYY_sxy = pYYsn[is][nx][ny];
          arr_fetch1(double) pYZ_sxy = pYZsn[is][nx][ny];
          arr_fetch1(double) pZZ_sxy = pZZsn[is][nx][ny];
          for(int nz=0;nz<nzn;nz++)
          {
            rho_sxy[nz] = invVOL*node_moments8_for_master[nz][0];
            Jx__sxy[nz] = invVOL*node_moments8_for_master[nz][1];
            Jy__sxy[nz] = invVOL*node_moments8_for_master[nz][2];
            Jz__sxy[nz] = invVOL*node_moments8_for_master[nz][3];
            pXX_sxy[nz] = invVOL*node_moments8_for_master[nz][4];
            pXY_sxy[nz] = invVOL*node_moments8_for_master[nz][5];
            pXZ_sxy[nz] = invVOL*node_moments8_for_master[nz][6];
            pYY_sxy[nz] = invVOL*node_moments8_for_master[nz][7];
            pYZ_sxy[nz] = invVOL*node_moments2_for_master[nz][0];
            pZZ_sxy[nz] = invVOL*node_moments2_for_master[nz][1];
          }
        }
      }
      if(!this_thread) timeTasks_end_task(TimeTasks::MOMENT_REDUCTION);
    }
  }

  // deallocate memory per mesh node for accumulating moments
  //
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // call destructor to deallocate arrays
    node_moments_first8_per_thr[thread_num].~array3<F64vec8>();
    node_moments_last2_per_thr[thread_num].~array4<double>();
  }
  free(node_moments_first8_per_thr);
  free(node_moments_last2_per_thr);

  // deallocate memory for accumulating moments
  //
  for(int thread_num=0;thread_num<num_threads;thread_num++)
  {
    // deallocate array to accumulate moments for thread
    cell_moments_per_thr[thread_num].~array4<F64vec8>();
  }
  free(cell_moments_per_thr);

  for (int i = 0; i < ns; i++)
  {
    communicateGhostP2G(i);
  }
#endif // __MIC__
}

inline void compute_moments(double velmoments[10], double weights[8],
  int i,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  ALIGNED(x);
  ALIGNED(y);
  ALIGNED(z);
  ALIGNED(u);
  ALIGNED(v);
  ALIGNED(w);
  ALIGNED(q);
  // compute the quadratic moments of velocity
  //
  const double ui=u[i];
  const double vi=v[i];
  const double wi=w[i];
  const double uui=ui*ui;
  const double uvi=ui*vi;
  const double uwi=ui*wi;
  const double vvi=vi*vi;
  const double vwi=vi*wi;
  const double wwi=wi*wi;
  //double velmoments[10];
  velmoments[0] = 1.;
  velmoments[1] = ui;
  velmoments[2] = vi;
  velmoments[3] = wi;
  velmoments[4] = uui;
  velmoments[5] = uvi;
  velmoments[6] = uwi;
  velmoments[7] = vvi;
  velmoments[8] = vwi;
  velmoments[9] = wwi;

  // compute the weights to distribute the moments
  //
  //double weights[8];
  const double abs_xpos = x[i];
  const double abs_ypos = y[i];
  const double abs_zpos = z[i];
  const double rel_xpos = abs_xpos - xstart;
  const double rel_ypos = abs_ypos - ystart;
  const double rel_zpos = abs_zpos - zstart;
  const double cxm1_pos = rel_xpos * inv_dx;
  const double cym1_pos = rel_ypos * inv_dy;
  const double czm1_pos = rel_zpos * inv_dz;
  //if(true)
  //{
  //  const int cx_inf = int(floor(cxm1_pos));
  //  const int cy_inf = int(floor(cym1_pos));
  //  const int cz_inf = int(floor(czm1_pos));
  //  assert_eq(cx-1,cx_inf);
  //  assert_eq(cy-1,cy_inf);
  //  assert_eq(cz-1,cz_inf);
  //}
  // fraction of the distance from the right of the cell
  const double w1x = cx - cxm1_pos;
  const double w1y = cy - cym1_pos;
  const double w1z = cz - czm1_pos;
  // fraction of distance from the left
  const double w0x = 1-w1x;
  const double w0y = 1-w1y;
  const double w0z = 1-w1z;
  // we are calculating a charge moment.
  const double qi=q[i];
  const double weight0 = qi*w0x;
  const double weight1 = qi*w1x;
  const double weight00 = weight0*w0y;
  const double weight01 = weight0*w1y;
  const double weight10 = weight1*w0y;
  const double weight11 = weight1*w1y;
  weights[0] = weight00*w0z; // weight000
  weights[1] = weight00*w1z; // weight001
  weights[2] = weight01*w0z; // weight010
  weights[3] = weight01*w1z; // weight011
  weights[4] = weight10*w0z; // weight100
  weights[5] = weight10*w1z; // weight101
  weights[6] = weight11*w0z; // weight110
  weights[7] = weight11*w1z; // weight111
}

// add particle to moments
inline void add_moments_for_pcl(double momentsAcc[8][10],
  int i,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  double velmoments[10];
  double weights[8];
  compute_moments(velmoments,weights,
    i, x, y, z, u, v, w, q,
    xstart, ystart, zstart,
    inv_dx, inv_dy, inv_dz,
    cx, cy, cz);

  // add moments for this particle
  {
    // which is the superior order for the following loop?
    for(int c=0; c<8; c++)
    for(int m=0; m<10; m++)
    {
      momentsAcc[c][m] += velmoments[m]*weights[c];
    }
  }
}


// vectorized version of previous method
// 
inline void add_moments_for_pcl_vec(double momentsAccVec[8][10][8],
  double velmoments[10][8], double weights[8][8],
  int i,
  int imod,
  double const * const x,
  double const * const y,
  double const * const z,
  double const * const u,
  double const * const v,
  double const * const w,
  double const * const q,
  double xstart,
  double ystart,
  double zstart,
  double inv_dx,
  double inv_dy,
  double inv_dz,
  int cx,
  int cy,
  int cz)
{
  ALIGNED(x);
  ALIGNED(y);
  ALIGNED(z);
  ALIGNED(u);
  ALIGNED(v);
  ALIGNED(w);
  ALIGNED(q);
  // compute the quadratic moments of velocity
  //
  const double ui=u[i];
  const double vi=v[i];
  const double wi=w[i];
  const double uui=ui*ui;
  const double uvi=ui*vi;
  const double uwi=ui*wi;
  const double vvi=vi*vi;
  const double vwi=vi*wi;
  const double wwi=wi*wi;
  //double velmoments[10];
  velmoments[0][imod] = 1.;
  velmoments[1][imod] = ui;
  velmoments[2][imod] = vi;
  velmoments[3][imod] = wi;
  velmoments[4][imod] = uui;
  velmoments[5][imod] = uvi;
  velmoments[6][imod] = uwi;
  velmoments[7][imod] = vvi;
  velmoments[8][imod] = vwi;
  velmoments[9][imod] = wwi;

  // compute the weights to distribute the moments
  //
  //double weights[8];
  const double abs_xpos = x[i];
  const double abs_ypos = y[i];
  const double abs_zpos = z[i];
  const double rel_xpos = abs_xpos - xstart;
  const double rel_ypos = abs_ypos - ystart;
  const double rel_zpos = abs_zpos - zstart;
  const double cxm1_pos = rel_xpos * inv_dx;
  const double cym1_pos = rel_ypos * inv_dy;
  const double czm1_pos = rel_zpos * inv_dz;
  //if(true)
  //{
  //  const int cx_inf = int(floor(cxm1_pos));
  //  const int cy_inf = int(floor(cym1_pos));
  //  const int cz_inf = int(floor(czm1_pos));
  //  assert_eq(cx-1,cx_inf);
  //  assert_eq(cy-1,cy_inf);
  //  assert_eq(cz-1,cz_inf);
  //}
  // fraction of the distance from the right of the cell
  const double w1x = cx - cxm1_pos;
  const double w1y = cy - cym1_pos;
  const double w1z = cz - czm1_pos;
  // fraction of distance from the left
  const double w0x = 1-w1x;
  const double w0y = 1-w1y;
  const double w0z = 1-w1z;
  // we are calculating a charge moment.
  const double qi=q[i];
  const double weight0 = qi*w0x;
  const double weight1 = qi*w1x;
  const double weight00 = weight0*w0y;
  const double weight01 = weight0*w1y;
  const double weight10 = weight1*w0y;
  const double weight11 = weight1*w1y;
  weights[0][imod] = weight00*w0z; // weight000
  weights[1][imod] = weight00*w1z; // weight001
  weights[2][imod] = weight01*w0z; // weight010
  weights[3][imod] = weight01*w1z; // weight011
  weights[4][imod] = weight10*w0z; // weight100
  weights[5][imod] = weight10*w1z; // weight101
  weights[6][imod] = weight11*w0z; // weight110
  weights[7][imod] = weight11*w1z; // weight111

  // add moments for this particle
  {
    for(int c=0; c<8; c++)
    for(int m=0; m<10; m++)
    {
      momentsAccVec[c][m][imod] += velmoments[m][imod]*weights[c][imod];
    }
  }
}

void EMfields3D::sumMoments_vectorized(const Particles3Dcomm* part)
{
  const Grid *grid = &get_grid();

  const double inv_dx = grid->get_invdx();
  const double inv_dy = grid->get_invdy();
  const double inv_dz = grid->get_invdz();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  #pragma omp parallel
  {
  for (int species_idx = 0; species_idx < ns; species_idx++)
  {
    const Particles3Dcomm& pcls = part[species_idx];
    assert_eq(pcls.get_particleType(), ParticleType::SoA);
    const int is = pcls.get_species_num();
    assert_eq(species_idx,is);

    double const*const x = pcls.getXall();
    double const*const y = pcls.getYall();
    double const*const z = pcls.getZall();
    double const*const u = pcls.getUall();
    double const*const v = pcls.getVall();
    double const*const w = pcls.getWall();
    double const*const q = pcls.getQall();

    const int nop = pcls.getNOP();
    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
    Moments10& speciesMoments10 = fetch_moments10Array(0);
    arr4_double moments = speciesMoments10.fetch_arr();
    //
    // moments.setmode(ompmode::ompfor);
    //moments.setall(0.);
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    #pragma omp for // because shared
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    
    // prevent threads from writing to the same location
    for(int cxmod2=0; cxmod2<2; cxmod2++)
    for(int cymod2=0; cymod2<2; cymod2++)
    // each mesh cell is handled by its own thread
    #pragma omp for collapse(2)
    for(int cx=cxmod2;cx<nxc;cx+=2)
    for(int cy=cymod2;cy<nyc;cy+=2)
    for(int cz=0;cz<nzc;cz++)
    {
     //dprint(cz);
     // index of interface to right of cell
     const int ix = cx + 1;
     const int iy = cy + 1;
     const int iz = cz + 1;
     {
      // reference the 8 nodes to which we will
      // write moment data for particles in this mesh cell.
      //
      arr1_double_fetch momentsArray[8];
      arr2_double_fetch moments00 = moments[ix][iy];
      arr2_double_fetch moments01 = moments[ix][cy];
      arr2_double_fetch moments10 = moments[cx][iy];
      arr2_double_fetch moments11 = moments[cx][cy];
      momentsArray[0] = moments00[iz]; // moments000 
      momentsArray[1] = moments00[cz]; // moments001 
      momentsArray[2] = moments01[iz]; // moments010 
      momentsArray[3] = moments01[cz]; // moments011 
      momentsArray[4] = moments10[iz]; // moments100 
      momentsArray[5] = moments10[cz]; // moments101 
      momentsArray[6] = moments11[iz]; // moments110 
      momentsArray[7] = moments11[cz]; // moments111 

      const int numpcls_in_cell = pcls.get_numpcls_in_bucket(cx,cy,cz);
      const int bucket_offset = pcls.get_bucket_offset(cx,cy,cz);
      const int bucket_end = bucket_offset+numpcls_in_cell;

      bool vectorized=false;
      if(!vectorized)
      {
        // accumulators for moments per each of 8 threads
        double momentsAcc[8][10];
        memset(momentsAcc,0,sizeof(double)*8*10);
        for(int i=bucket_offset; i<bucket_end; i++)
        {
          add_moments_for_pcl(momentsAcc, i,
            x, y, z, u, v, w, q,
            xstart, ystart, zstart,
            inv_dx, inv_dy, inv_dz,
            cx, cy, cz);
        }
        for(int c=0; c<8; c++)
        for(int m=0; m<10; m++)
        {
          momentsArray[c][m] += momentsAcc[c][m];
        }
      }
      if(vectorized)
      {
        double velmoments[10][8];
        double weights[8][8];
        double momentsAccVec[8][10][8];
        memset(momentsAccVec,0,sizeof(double)*8*10*8);
        #pragma simd
        for(int i=bucket_offset; i<bucket_end; i++)
        {
          add_moments_for_pcl_vec(momentsAccVec, velmoments, weights,
            i, i%8,
            x, y, z, u, v, w, q,
            xstart, ystart, zstart,
            inv_dx, inv_dy, inv_dz,
            cx, cy, cz);
        }
        for(int c=0; c<8; c++)
        for(int m=0; m<10; m++)
        for(int i=0; i<8; i++)
        {
          momentsArray[c][m] += momentsAccVec[c][m][i];
        }
      }
     }
    }
    #pragma omp master
    { timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION); }

    // reduction
    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION); }
    {
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++){
      for(int j=0;j<nyn;j++){
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] = invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] = invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] = invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] = invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] = invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] = invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] = invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] = invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] = invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] = invVOL*moments[i][j][k][9];
      }}}
    }
    #pragma omp master
    { timeTasks_end_task(TimeTasks::MOMENT_REDUCTION); }
    // uncomment this and remove the loop below
    // when we change to use asynchronous communication.
    // communicateGhostP2G(is);
  }
  }
  for (int i = 0; i < ns; i++)
  {
    communicateGhostP2G(i);
  }
}

void EMfields3D::sumMoments_vectorized_AoS(const Particles3Dcomm* part)
{
  const Grid *grid = &get_grid();

  const double inv_dx = grid->get_invdx();
  const double inv_dy = grid->get_invdy();
  const double inv_dz = grid->get_invdz();
  const int nxn = grid->getNXN();
  const int nyn = grid->getNYN();
  const int nzn = grid->getNZN();
  const double xstart = grid->getXstart();
  const double ystart = grid->getYstart();
  const double zstart = grid->getZstart();
  #pragma omp parallel
  {
  for (int species_idx = 0; species_idx < ns; species_idx++)
  {
    const Particles3Dcomm& pcls = part[species_idx];
    assert_eq(pcls.get_particleType(), ParticleType::AoS);
    const int is = pcls.get_species_num();
    assert_eq(species_idx,is);

    const int nop = pcls.getNOP();
    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOMENT_ACCUMULATION); }
    Moments10& speciesMoments10 = fetch_moments10Array(0);
    arr4_double moments = speciesMoments10.fetch_arr();
    //
    // moments.setmode(ompmode::ompfor);
    //moments.setall(0.);
    double *moments1d = &moments[0][0][0][0];
    int moments1dsize = moments.get_size();
    #pragma omp for // because shared
    for(int i=0; i<moments1dsize; i++) moments1d[i]=0;
    
    // prevent threads from writing to the same location
    for(int cxmod2=0; cxmod2<2; cxmod2++)
    for(int cymod2=0; cymod2<2; cymod2++)
    // each mesh cell is handled by its own thread
    #pragma omp for collapse(2)
    for(int cx=cxmod2;cx<nxc;cx+=2)
    for(int cy=cymod2;cy<nyc;cy+=2)
    for(int cz=0;cz<nzc;cz++)
    {
     //dprint(cz);
     // index of interface to right of cell
     const int ix = cx + 1;
     const int iy = cy + 1;
     const int iz = cz + 1;
     {
      // reference the 8 nodes to which we will
      // write moment data for particles in this mesh cell.
      //
      arr1_double_fetch momentsArray[8];
      arr2_double_fetch moments00 = moments[ix][iy];
      arr2_double_fetch moments01 = moments[ix][cy];
      arr2_double_fetch moments10 = moments[cx][iy];
      arr2_double_fetch moments11 = moments[cx][cy];
      momentsArray[0] = moments00[iz]; // moments000 
      momentsArray[1] = moments00[cz]; // moments001 
      momentsArray[2] = moments01[iz]; // moments010 
      momentsArray[3] = moments01[cz]; // moments011 
      momentsArray[4] = moments10[iz]; // moments100 
      momentsArray[5] = moments10[cz]; // moments101 
      momentsArray[6] = moments11[iz]; // moments110 
      momentsArray[7] = moments11[cz]; // moments111 

      // accumulator for moments per each of 8 threads
      double momentsAcc[8][10][8];
      const int numpcls_in_cell = pcls.get_numpcls_in_bucket(cx,cy,cz);
      const int bucket_offset = pcls.get_bucket_offset(cx,cy,cz);
      const int bucket_end = bucket_offset+numpcls_in_cell;

      // data is not stride-1, so we do *not* use
      // #pragma simd
      {
        // accumulators for moments per each of 8 threads
        double momentsAcc[8][10];
        memset(momentsAcc,0,sizeof(double)*8*10);
        for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
        {
          const SpeciesParticle* pcl = &pcls.get_pcl(pidx);
          // This depends on the fact that the memory
          // occupied by a particle coincides with
          // the alignment interval (64 bytes)
          ALIGNED(pcl);
          double velmoments[10];
          double weights[8];
          // compute the quadratic moments of velocity
          //
          const double ui=pcl->get_u();
          const double vi=pcl->get_v();
          const double wi=pcl->get_w();
          const double uui=ui*ui;
          const double uvi=ui*vi;
          const double uwi=ui*wi;
          const double vvi=vi*vi;
          const double vwi=vi*wi;
          const double wwi=wi*wi;
          //double velmoments[10];
          velmoments[0] = 1.;
          velmoments[1] = ui;
          velmoments[2] = vi;
          velmoments[3] = wi;
          velmoments[4] = uui;
          velmoments[5] = uvi;
          velmoments[6] = uwi;
          velmoments[7] = vvi;
          velmoments[8] = vwi;
          velmoments[9] = wwi;
        
          // compute the weights to distribute the moments
          //
          //double weights[8];
          const double abs_xpos = pcl->get_x();
          const double abs_ypos = pcl->get_y();
          const double abs_zpos = pcl->get_z();
          const double rel_xpos = abs_xpos - xstart;
          const double rel_ypos = abs_ypos - ystart;
          const double rel_zpos = abs_zpos - zstart;
          const double cxm1_pos = rel_xpos * inv_dx;
          const double cym1_pos = rel_ypos * inv_dy;
          const double czm1_pos = rel_zpos * inv_dz;
          //if(true)
          //{
          //  const int cx_inf = int(floor(cxm1_pos));
          //  const int cy_inf = int(floor(cym1_pos));
          //  const int cz_inf = int(floor(czm1_pos));
          //  assert_eq(cx-1,cx_inf);
          //  assert_eq(cy-1,cy_inf);
          //  assert_eq(cz-1,cz_inf);
          //}
          // fraction of the distance from the right of the cell
          const double w1x = cx - cxm1_pos;
          const double w1y = cy - cym1_pos;
          const double w1z = cz - czm1_pos;
          // fraction of distance from the left
          const double w0x = 1-w1x;
          const double w0y = 1-w1y;
          const double w0z = 1-w1z;
          // we are calculating a charge moment.
          const double qi=pcl->get_q();
          const double weight0 = qi*w0x;
          const double weight1 = qi*w1x;
          const double weight00 = weight0*w0y;
          const double weight01 = weight0*w1y;
          const double weight10 = weight1*w0y;
          const double weight11 = weight1*w1y;
          weights[0] = weight00*w0z; // weight000
          weights[1] = weight00*w1z; // weight001
          weights[2] = weight01*w0z; // weight010
          weights[3] = weight01*w1z; // weight011
          weights[4] = weight10*w0z; // weight100
          weights[5] = weight10*w1z; // weight101
          weights[6] = weight11*w0z; // weight110
          weights[7] = weight11*w1z; // weight111
        
          // add moments for this particle
          {
            // which is the superior order for the following loop?
            for(int c=0; c<8; c++)
            for(int m=0; m<10; m++)
            {
              momentsAcc[c][m] += velmoments[m]*weights[c];
            }
          }
        }
        for(int c=0; c<8; c++)
        for(int m=0; m<10; m++)
        {
          momentsArray[c][m] += momentsAcc[c][m];
        }
      }
     }
    }
    #pragma omp master
    { timeTasks_end_task(TimeTasks::MOMENT_ACCUMULATION); }

    // reduction
    #pragma omp master
    { timeTasks_begin_task(TimeTasks::MOMENT_REDUCTION); }
    {
      #pragma omp for collapse(2)
      for(int i=0;i<nxn;i++){
      for(int j=0;j<nyn;j++){
      for(int k=0;k<nzn;k++)
      {
        rhons[is][i][j][k] = invVOL*moments[i][j][k][0];
        Jxs  [is][i][j][k] = invVOL*moments[i][j][k][1];
        Jys  [is][i][j][k] = invVOL*moments[i][j][k][2];
        Jzs  [is][i][j][k] = invVOL*moments[i][j][k][3];
        pXXsn[is][i][j][k] = invVOL*moments[i][j][k][4];
        pXYsn[is][i][j][k] = invVOL*moments[i][j][k][5];
        pXZsn[is][i][j][k] = invVOL*moments[i][j][k][6];
        pYYsn[is][i][j][k] = invVOL*moments[i][j][k][7];
        pYZsn[is][i][j][k] = invVOL*moments[i][j][k][8];
        pZZsn[is][i][j][k] = invVOL*moments[i][j][k][9];
      }}}
    }
    #pragma omp master
    { timeTasks_end_task(TimeTasks::MOMENT_REDUCTION); }
    // uncomment this and remove the loop below
    // when we change to use asynchronous communication.
    // communicateGhostP2G(is);
  }
  }
  for (int i = 0; i < ns; i++)
  {
    communicateGhostP2G(i);
  }
}

/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(arr3_double vectPhys, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++)
        vectPhys[i][j][k] = *vectSolver++;

}
/** method to convert a 1D field in a 3D field not considering guard cells*/
void solver2phys(arr3_double vectPhys1, arr3_double vectPhys2, arr3_double vectPhys3, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vectPhys1[i][j][k] = *vectSolver++;
        vectPhys2[i][j][k] = *vectSolver++;
	vectPhys3[i][j][k] = *vectSolver++;
	}
}
/** method to convert a 3D field into a 1D field whithin given index range */
void solver2phys(arr3_double vectPhys1, arr3_double vectPhys2, arr3_double vectPhys3, double *vectSolver,		 
		 int iMin, int iMax , int jMin, int jMax, int kMin, int kMax){
  for (register int i=iMin; i <= iMax; i++)
    for (register int j=jMin; j <= jMax; j++)
      for (register int k=kMin; k <= kMax; k++){
	vectPhys1[i][j][k] = *vectSolver++;
	vectPhys2[i][j][k] = *vectSolver++;
	vectPhys3[i][j][k] = *vectSolver++;
      }
}
void solver2phys(arr3_double vectPhys, double *vectSolver,		 
		 int iMin, int iMax , int jMin, int jMax, int kMin, int kMax){
  for (register int i=iMin; i <= iMax; i++)
    for (register int j=jMin; j <= jMax; j++)
      for (register int k=kMin; k <= kMax; k++){
	vectPhys[i][j][k] = *vectSolver++;
      }
}

/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(double *vectSolver, const arr3_double vectPhys, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++)
        *vectSolver++ = vectPhys.get(i,j,k);
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
void phys2solver(double *vectSolver, const arr3_double vectPhys1, const arr3_double vectPhys2, const arr3_double vectPhys3, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        *vectSolver++ = vectPhys1.get(i,j,k);
        *vectSolver++ = vectPhys2.get(i,j,k);
        *vectSolver++ = vectPhys3.get(i,j,k);
      }
}
/** method to convert a 3D field into a 1D field whithin given index range */
/** Necessary for MHD-IPIC coupling */
void phys2solver(double* vectSolver, const arr3_double vectPhys1,
		 const arr3_double vectPhys2,const arr3_double vectPhys3,
		 int iMin, int iMax , int jMin, int jMax, int kMin, int kMax){
  for (register int i=iMin; i <= iMax; i++)
    for (register int j=jMin; j <= jMax; j++)
      for (register int k=kMin; k <= kMax; k++){
	*vectSolver++ =  vectPhys1.get(i,j,k);
	*vectSolver++ =  vectPhys2.get(i,j,k);
	*vectSolver++ =  vectPhys3.get(i,j,k);
      }
}
void phys2solver(double* vectSolver, const arr3_double vectPhys,
		 int iMin, int iMax , int jMin, int jMax, int kMin, int kMax){
  for (register int i=iMin; i <= iMax; i++)
    for (register int j=jMin; j <= jMax; j++)
      for (register int k=kMin; k <= kMax; k++){
	*vectSolver++ =  vectPhys.get(i,j,k);
      }
}


/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::calculateE(int cycle)
{
  const Collective *col = &get_col();
  const VirtualTopology3D * vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;

  // Solve dE = E^{n+1} - E^{n} instead E^{n+1}
  bool doSolveForChange = col->getCase()=="BATSRUS";  
  //doSolveForChange = false;

  bool useIPIC3DSolver = true; 
#ifdef BATSRUS
  useIPIC3DSolver = col->get_useIPIC3DSolver();
#endif
  
  array3_double gradPHIX (nxn, nyn, nzn);
  array3_double gradPHIY (nxn, nyn, nzn);
  array3_double gradPHIZ (nxn, nyn, nzn);

  double *xkrylov = new double[n3SolveNode];  // 3 E components
  double *bkrylov = new double[n3SolveNode];  // 3 components  
  double *dxkrylov;
  if(doSolveForChange)
    dxkrylov = new double[n3SolveNode]; // 3 components  
  // set to zero all the stuff 
  eqValue(0.0, xkrylov, n3SolveNode);
  eqValue(0.0, bkrylov, n3SolveNode);
  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, gradPHIX, nxn, nyn, nzn);
  eqValue(0.0, gradPHIY, nxn, nyn, nzn);
  eqValue(0.0, gradPHIZ, nxn, nyn, nzn);
  
  if(doSolveForChange)
    eqValue(0.0, dxkrylov, n3SolveNode);


  if (PoissonCorrection &&  cycle%PoissonCorrectionCycle == 0) {
#ifdef BATSRUS
    calculate_PHI(iPIC3D_PoissonImage,
		  col->get_PoissonTol(),
		  col->get_PoissonIter(),true);
#else
    calculate_PHI(NULL, 0.1, 20, true);
#endif
    
    // calculate the gradient
    grid->gradC2N(gradPHIX, gradPHIY, gradPHIZ, PHI);
    sub(Ex, gradPHIX, nxn, nyn, nzn);
    sub(Ey, gradPHIY, nxn, nyn, nzn);
    sub(Ez, gradPHIZ, nxn, nyn, nzn);
  }                             // end of divergence cleaning

  if (vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER ***" << endl;

  // prepare the source 
  MaxwellSource(bkrylov);
  phys2solver(xkrylov,Ex,Ey,Ez,inminsolve,inmaxsolve,jnminsolve,jnmaxsolve,knminsolve,knmaxsolve);
  
  if(doSolveForChange){
    MaxwellImage(dxkrylov,xkrylov,false);   // LHS at t_n  moved to RHS as we solve for the difference
    addscale(-1.0,1.0,bkrylov,dxkrylov, n3SolveNode);    // Add to get the new b (RHS)
    eqValue (0.0, xkrylov, n3SolveNode);              // x=0 as we solve for the change
  }
  
#ifdef BATSRUS
  if(! useIPIC3DSolver){
    linear_solver_matvec_c = iPIC3D_MaxwellImage;

    int nVarSolve = 3;
    int nIter = col->get_EFieldIter();
    double EFieldTol = col->get_EFieldTol();
    int nDimIn = nDimMax;
    int nI = inmaxsolve - inminsolve + 1;
    int nJ = jnmaxsolve - jnminsolve + 1;
    int nK = knmaxsolve - knminsolve + 1;
    int nBlock = 1;
    MPI_Fint iComm = MPI_Comm_c2f(MPI_COMM_MYSIM);
    double precond_matrix_II[1][1];
    precond_matrix_II[0][0] = 0; 
    // parameter to choose preconditioner types
    //0:No precondition; 1: BILU; 2:DILU;
    //[-1,0): MBILU;   
    double PrecondParam=0;
    int lTest = vct->getCartesian_rank() == 0;

    linear_solver_wrapper("GMRES", &EFieldTol, &nIter, &nVarSolve, &nDimIn,
			  &nI, &nJ, &nK, &nBlock, &iComm, bkrylov,
			  xkrylov, &PrecondParam, precond_matrix_II[0], 
			  &lTest);
  }
#endif
  if(useIPIC3DSolver){
  // solve A.xKrylob = bKrylov for x being the change in E field (dE)
  GMRES(&Field::MaxwellImage, xkrylov, n3SolveNode,
	bkrylov, nGMRESRestart, 200, GMREStol,doSolveForChange, this);
  }

 
  
  // move from krylov space to physical space: Eth = dE
  solver2phys(Exth,Eyth,Ezth,xkrylov,inminsolve,inmaxsolve,jnminsolve,jnmaxsolve,knminsolve,knmaxsolve);


  if(doSolveForChange){
    // Caclulate the fully implicit solution E* = E_n + dE
    // Eth = Eth + 1.0*E, where Eth = dE originally, and Eth=E* after the update
    addscale(1.0,Exth,Ex,nxn,nyn,nzn);
    addscale(1.0,Eyth,Ey,nxn,nyn,nzn);
    addscale(1.0,Ezth,Ez,nxn,nyn,nzn);
  }

  #ifdef BATSRUS
  // Apply BATSRUS boundary condition to Eth. 
  if(col->getCase()=="BATSRUS") fixE_BATSRUS(Exth,Eyth,Ezth,false);
  #endif

  // apply smoothing to Eth
  smoothE();

  // E_n+1 = 1/theta*E* -(1.0-theta)/theta*E_n = [ E* + (theta - 1)E ] / theta  
  // E = 1/theta*Eth -(1.0-theta)/theta*E = [ Eth + (theta - 1)E ] / theta
  // for theta=1 (usual setting): E = Eth
  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);
  
  #ifdef BATSRUS
  // Apply BATSRUS boundary condition to E_n+1.
  // It is not necessary for theta == 1.
  if(col->getCase()=="BATSRUS") fixE_BATSRUS(Ex,Ey,Ez,false);
  #endif

  
  // communicate so the interpolation can have good values
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);

  if(col->getCase()!="BATSRUS"){
    // OpenBC Inflow: this needs to be integrate to Halo Exchange BC
    OpenBoundaryInflowE(Exth, Eyth, Ezth, nxn, nyn, nzn);
    OpenBoundaryInflowE(Ex, Ey, Ez, nxn, nyn, nzn);
  }

  // Calculate and store div(E)
  grid->divN2C(divEc, Ex, Ey, Ez);
  communicateCenterBC(nxc, nyc, nzc, divEc, 2, 2, 2, 2, 2, 2, vct, this);
  grid->interpC2N(divEn, divEc);
  
  const double particleTheta = col->get_particleTheta();
  
  const double c0 = (th==1? 0:(1-particleTheta)/(1-th));
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for(int k = 0; k < nzn; k++){
	Exthp[i][j][k] = (1-c0)*Ex[i][j][k] + c0*Exth[i][j][k];
	Eythp[i][j][k] = (1-c0)*Ey[i][j][k] + c0*Eyth[i][j][k];
	Ezthp[i][j][k] = (1-c0)*Ez[i][j][k] + c0*Ezth[i][j][k];	
      }

  //calc_energy_error(); 
  //fix_energy();

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;
  if(doSolveForChange) delete[] dxkrylov;

}


//------------------------------------------------------
void EMfields3D::calculate_PHI(MATVEC FuncImage, double krylovTol,
			       int nIter, bool useFloatPHI, bool doCalcError){
  /*
    Input:
    1) FuncImage: the function to calculate A*x for linear solver. 
    2) krylovTol: the tolerance of the iterative solver. Default is 1e-6.
    3) nIter: the maximum iteration number. Default is 50. 
    4) useFloatPHI: if useFloatPHI is true, then use float boundary 
       condition for PHI, otherwise, PHI is fixed to 0 at boundaries. 
       Default is 'true'.
   */
  const Collective *col = &get_col();
  const VirtualTopology3D * vct = &get_vct();
  const Grid *grid = &get_grid();

  bool useIPIC3DSolver = true; 
#ifdef BATSRUS
  useIPIC3DSolver = col->get_useIPIC3DSolver();
#endif

  double xkrylovPoisson[nSolveCell];
  double bkrylovPoisson[nSolveCell];
  array3_double divE     (nxc, nyc, nzc);

  eqValue(0.0, divE, nxc, nyc, nzc);
  eqValue(0.0, xkrylovPoisson, nSolveCell);
  eqValue(0.0, PHI, nxc, nyc, nzc);
  
  if(doCalcError){
    grid->divN2C(divE, Ex, Ey, Ez);
  }

  scale(tempC, rhoc, -FourPI, nxc, nyc, nzc);
  
  sum(divE, tempC, nxc, nyc, nzc);
  // move to krylov space

  phys2solver(bkrylovPoisson, divE, icMinSolve,icMaxSolve,
	      jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);
  if (vct->getCartesian_rank() == 0) cout << "*** DIVERGENCE CLEANING ***" << endl;
#ifdef BATSRUS
  if(! useIPIC3DSolver){
    linear_solver_matvec_c = FuncImage;

    int nVarSolve = 1;
    int nDimIn = nDimMax;
    int nI = icMaxSolve - icMinSolve + 1;
    int nJ = jcMaxSolve - jcMinSolve + 1;
    int nK = kcMaxSolve - kcMinSolve + 1;
    int nBlock = 1;
    MPI_Fint iComm = MPI_Comm_c2f(MPI_COMM_MYSIM);

    double **precond_matrix_II;
    
    int nVector = nVarSolve*nVarSolve*nI*nJ*nK;
    precond_matrix_II = newArr2(double, 2*nDimIn+1, nVector);

    int idx; 
    int nDimTrue = nDimIn; 
    if(nI==1) nDimTrue--;
    if(nJ==1) nDimTrue--;
    if(nK==1) nDimTrue--;
    
    for(int i = 0; i<nI; i++)
      for(int j = 0; j<nJ; j++)
	for(int k = 0; k<nK; k++)
	  {
	    idx = k + j*nK + i*nK*nJ;
	    precond_matrix_II[0][idx] =  -2*nDimTrue;	    
	    precond_matrix_II[1][idx] = (i !=    0)? 1:0; // iSub
	    precond_matrix_II[2][idx] = (i != nI-1)? 1:0; // iSup
	    precond_matrix_II[3][idx] = (j !=    0)? 1:0; // jSub
	    precond_matrix_II[4][idx] = (j != nJ-1)? 1:0; // jSup	    
	    precond_matrix_II[5][idx] = (k !=    0)? 1:0; // kSub
	    precond_matrix_II[6][idx] = (k != nK-1)? 1:0; // kSup      
    }

    // parameter to choose preconditioner types
    //0:No precondition; 1: BILU; 2:DILU;
    //[-1,0): MBILU;   
    double PrecondParam=2;
    int lTest = vct->getCartesian_rank() == 0;

    linear_solver_wrapper("GMRES", &krylovTol, &nIter, &nVarSolve, &nDimIn,
			  &nI, &nJ, &nK, &nBlock, &iComm, bkrylovPoisson,
			  xkrylovPoisson, &PrecondParam, precond_matrix_II[0], 
			  &lTest);
    delArr2(precond_matrix_II, 2*nDimIn+1);
  }
#endif

    
  if(useIPIC3DSolver){
  GMRES(&Field::PoissonImage, xkrylovPoisson, nSolveCell, bkrylovPoisson, nGMRESRestart, 200, krylovTol,false, this);
}

  solver2phys(PHI, xkrylovPoisson, icMinSolve,icMaxSolve,
    jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);

#ifdef BATSRUS
  if(useFloatPHI) fixPHI_BATSRUS();
#endif  
  communicateCenterBC(nxc, nyc, nzc, PHI, 2, 2, 2, 2, 2, 2, vct,this);

}


void EMfields3D::matvec_weight_correction(double *image, double *vector){

  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  // allocate 2 three dimensional service vectors
  array3_double temp(nxc, nyc, nzc);
  array3_double im(nxc, nyc, nzc);
  eqValue(0.0, image, nSolveCell);
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, icMinSolve,icMaxSolve,
	      jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);

  communicateCenterBC_P(nxc, nyc, nzc, temp, 2, 2, 2, 2, 2, 2, vct, this);

  int gp; // g' 
  for(int i =1; i<nxc-1; i++)
    for(int j = 1; j<nyc-1; j++)
      for(int k = 1; k<nzc-1; k++)
	  for(int i2 = i-1; i2 <= i+1; i2++)
	    for(int j2 = j-1; j2 <= j+1; j2++)
	      for(int k2 = k-1; k2 <= k+1; k2++){
		gp = (i2-i+1)*9+(j2-j+1)*3+k2-k+1;
		im[i][j][k] += temp[i2][j2][k2]*M_CI[i][j][k][gp];	       
		
	      }	

  scale(im,FourPI*FourPI,nxc,nyc,nzc);      
  
  // move from physical space to krylov space
  phys2solver(image, im, icMinSolve, icMaxSolve,
	      jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);

}

/*! Calculate sorgent for Maxwell solver */
void EMfields3D::MaxwellSource(double *bkrylov)
{
  const Collective *col = &get_col();
  const VirtualTopology3D * vct = &get_vct();
  const Grid *grid = &get_grid();

  eqValue(0.0, tempC, nxc, nyc, nzc);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, tempXN, nxn, nyn, nzn);
  eqValue(0.0, tempYN, nxn, nyn, nzn);
  eqValue(0.0, tempZN, nxn, nyn, nzn);
  eqValue(0.0, temp2X, nxn, nyn, nzn);
  eqValue(0.0, temp2Y, nxn, nyn, nzn);
  eqValue(0.0, temp2Z, nxn, nyn, nzn);


  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct, this);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct, this);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct, this);

  if (get_col().getCase()=="ForceFree") 	fixBforcefree();
  if (get_col().getCase()=="GEM")       	fixBnGEM();
  if (get_col().getCase()=="GEMnoPert") 	fixBnGEM();
  if (get_col().getCase()=="GEMDoubleHarris") 	fixBnGEM();
#ifdef BATSRUS
  if (get_col().getCase()=="BATSRUS")           fixB_BATSRUS();
#endif
  
  // OpenBC:
  if (get_col().getCase()!="BATSRUS")
    OpenBoundaryInflowB(Bxc,Byc,Bzc,nxc,nyc,nzc);

  if (get_col().getCase()=="GEM")       		fixBcGEM();
  if (get_col().getCase()=="GEMnoPert") 		fixBcGEM();
  if (get_col().getCase()=="GEMDoubleHarris") 	fixBcGEM();

  // prepare curl of B for known term of Maxwell solver: for the source term
  grid->curlC2N(tempXN, tempYN, tempZN, Bxc, Byc, Bzc);
  scale(temp2X, Jxh, -FourPI / c, nxn, nyn, nzn);
  scale(temp2Y, Jyh, -FourPI / c, nxn, nyn, nzn);
  scale(temp2Z, Jzh, -FourPI / c, nxn, nyn, nzn);

  sum(temp2X, tempXN, nxn, nyn, nzn);
  sum(temp2Y, tempYN, nxn, nyn, nzn);
  sum(temp2Z, tempZN, nxn, nyn, nzn);
  scale(temp2X, delt, nxn, nyn, nzn);
  scale(temp2Y, delt, nxn, nyn, nzn);
  scale(temp2Z, delt, nxn, nyn, nzn);

  
  // -gradRhoRatio*delt*delt*grad(4pi*rhoh)---------
  double gradRhoRatio = get_col().get_gradRhoRatio(); 
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct, this);
  grid->gradC2N(tempX, tempY, tempZ, rhoh);      
  scale(tempX, -delt * delt * FourPI * gradRhoRatio, nxn, nyn, nzn);
  scale(tempY, -delt * delt * FourPI * gradRhoRatio, nxn, nyn, nzn);
  scale(tempZ, -delt * delt * FourPI * gradRhoRatio, nxn, nyn, nzn);
  //----------------------------

  // sum E, past values
  sum(tempX, Ex, nxn, nyn, nzn);
  sum(tempY, Ey, nxn, nyn, nzn);
  sum(tempZ, Ez, nxn, nyn, nzn);

  
  // The theta used for particles (pth) may be different from the theta for 
  // EM field (th). E^(n+pth) = (1-pth/th)*E^n + pth/th*E^(n+th); 
  double cn = -(1 - col->get_particleTheta()/th);
  if(fabs(cn)>1e-6){
    MUdot(Dx, Dy, Dz, Ex, Ey, Ez);
    scale(Dx, cn, nxn, nyn, nzn);
    scale(Dy, cn, nxn, nyn, nzn);
    scale(Dz, cn, nxn, nyn, nzn);
    sum(tempX, Dx, nxn, nyn, nzn);
    sum(tempY, Dy, nxn, nyn, nzn);
    sum(tempZ, Dz, nxn, nyn, nzn);  
  }

  // sum curl(B) + jhat part
  sum(tempX, temp2X, nxn, nyn, nzn);
  sum(tempY, temp2Y, nxn, nyn, nzn);
  sum(tempZ, temp2Z, nxn, nyn, nzn);
  
  if (get_col().getCase()!="BATSRUS"){
    // Boundary condition in the known term
    // boundary condition: Xleft
    if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
      perfectConductorLeftS(tempX, tempY, tempZ, 0);
    // boundary condition: Xright
    if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
      perfectConductorRightS(tempX, tempY, tempZ, 0);
    // boundary condition: Yleft
    if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
      perfectConductorLeftS(tempX, tempY, tempZ, 1);
    // boundary condition: Yright
    if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
      perfectConductorRightS(tempX, tempY, tempZ, 1);
    // boundary condition: Zleft
    if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
      perfectConductorLeftS(tempX, tempY, tempZ, 2);
    // boundary condition: Zright
    if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
      perfectConductorRightS(tempX, tempY, tempZ, 2);
  }

  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, inminsolve,inmaxsolve,jnminsolve,jnmaxsolve,knminsolve,knmaxsolve);  
}


/*! Mapping of Maxwell image to give to solver */
//
// In the field solver, there is one layer of ghost cells.  The
// nodes on the ghost cells define two outer layers of nodes: the
// outermost nodes are clearly in the interior of the neighboring
// subdomain and can naturally be referred to as "ghost nodes",
// but the second-outermost layer is on the boundary between
// subdomains and thus does not clearly belong to any one process.
// Refer to these shared nodes as "boundary nodes".
//
// To compute the laplacian, we first compute the gradient
// at the center of each cell by differencing the values at
// the corners of the cell.  We then compute the laplacian
// (i.e. the divergence of the gradient) at each node by
// differencing the cell-center values in the cells sharing
// the node.
//
// The laplacian is required to be defined on all boundary
// and interior nodes.  
// 
// In the krylov solver, we make no attempt to use or to
// update the (outer) ghost nodes, and we assume (presumably
// correctly) that the boundary nodes are updated identically
// by all processes that share them.  Therefore, we must
// communicate gradient values in the ghost cells.  The
// subsequent computation of the divergence requires that
// this boundary communication first complete.
//
// An alternative way would be to communicate outer ghost node
// values after each update of Eth.  In this case, there would
// be no need for the 10=3*3+1 boundary communications in the body
// of MaxwellImage() entailed in the calls to lapN2N plus the
// call needed prior to the call to gradC2N.  Of course,
// we would then need to communicate the 3 components of the
// electric field for the outer ghost nodes prior to each call
// to MaxwellImage().  This second alternative would thus reduce
// communication by over a factor of 3.  Essentially, we would
// replace the cost of communicating cell-centered differences
// for ghost cell values with the cost of directly computing them.
//
// Also, while this second method does not increase the potential
// to avoid exposing latency, it can make it easier to do so.
//
// Another change that I would propose: define:
//
//   array4_double physical_vector(3,nxn,nyn,nzn);
//   arr3_double vectX = physical_vector[0];
//   arr3_double vectY = physical_vector[1];
//   arr3_double vectZ = physical_vector[2];
//   vector = &physical_vector[0][0][0][0];
//
// It is currently the case that boundary nodes are
// duplicated in "vector" and therefore receive a weight
// that is twice, four times, or eight times as much as
// other nodes in the Krylov inner product.  The definitions
// above would imply that ghost nodes also appear in the
// inner product.  To avoid this issue, we could simply zero
// ghost nodes before returning to the Krylov solver.  With
// the definitions above, phys2solver() would simply zero
// the ghost nodes and solver2phys() would populate them via
// communication.  Note that it would also be possible, if
// desired, to give duplicated nodes equal weight by
// rescaling their values in these two methods.
//
void EMfields3D::MaxwellImage(double *im, double* vector, bool doSolveForChange)
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  eqValue(0.0, im, n3SolveNode);
  eqValue(0.0, imageX, nxn, nyn, nzn);
  eqValue(0.0, imageY, nxn, nyn, nzn);
  eqValue(0.0, imageZ, nxn, nyn, nzn);
  eqValue(0.0, tempX, nxn, nyn, nzn);
  eqValue(0.0, tempY, nxn, nyn, nzn);
  eqValue(0.0, tempZ, nxn, nyn, nzn);
  eqValue(0.0, Dx, nxn, nyn, nzn);
  eqValue(0.0, Dy, nxn, nyn, nzn);
  eqValue(0.0, Dz, nxn, nyn, nzn);
  // move from krylov space to physical space
  solver2phys(vectX,vectY,vectZ,vector,inminsolve,inmaxsolve,jnminsolve,jnmaxsolve,knminsolve,knmaxsolve);  
  
#ifdef BATSRUS
  if(col->getCase()=="BATSRUS") fixE_BATSRUS(vectX,vectY,vectZ,doSolveForChange);
#endif

  grid->lapN2N(imageX, vectX,this);
  grid->lapN2N(imageY, vectY,this);
  grid->lapN2N(imageZ, vectZ,this);

  // (-delt*delt - cDiff)*lap(E) --------------- (1)
  double cDiff = get_col().get_cDiff()*0.5*dx*dx;
  scale(imageX, -delt * delt - cDiff, nxn, nyn, nzn);
  scale(imageY, -delt * delt - cDiff, nxn, nyn, nzn);
  scale(imageZ, -delt * delt - cDiff, nxn, nyn, nzn);
  //------------------------------------------------


  // delt*delt*(1-gradRhoRatio)*grad(div(E))--------- (2)  
  double coefDiffusion = get_col().get_ratioDivC2C();
  calc_gradDivE(tempX, tempY, tempZ, vectX, vectY, vectZ, coefDiffusion);

  double gradRhoRatio = col->get_gradRhoRatio();
  double coef = delt*delt*(1-gradRhoRatio);  
  scale(tempX, coef, nxn, nyn, nzn);
  scale(tempY, coef, nxn, nyn, nzn);
  scale(tempZ, coef, nxn, nyn, nzn);
  sum(imageX, tempX, nxn, nyn, nzn);
  sum(imageY, tempY, nxn, nyn, nzn);
  sum(imageZ, tempZ, nxn, nyn, nzn);
  //----------------------------------------


  // -delt*delt*gradRhoRatio*grad(div(D))--------- (3)
  MUdot(Dx, Dy, Dz, vectX, vectY, vectZ);
  grid->divN2C(divC, Dx, Dy, Dz);  
  communicateCenterBC(nxc, nyc, nzc, divC, 2, 2, 2, 2, 2, 2, vct, this);
  grid->gradC2N(tempX, tempY, tempZ, divC);
  scale(tempX, -delt*delt*gradRhoRatio, nxn, nyn, nzn);
  scale(tempY, -delt*delt*gradRhoRatio, nxn, nyn, nzn);
  scale(tempZ, -delt*delt*gradRhoRatio, nxn, nyn, nzn);
  sum(imageX, tempX, nxn, nyn, nzn);
  sum(imageY, tempY, nxn, nyn, nzn);
  sum(imageZ, tempZ, nxn, nyn, nzn);
  //----------------------------------------------------------   

  // (1)+(2)+(3)+D
  // The theta used for particles (pth) may be different from the theta for 
  // EM field (th). E^(n+pth) = (1-pth/th)*E^n + pth/th*E^(n+th); 
  double cth = col->get_particleTheta()/th;
  scale(Dx, cth, nxn, nyn, nzn);
  scale(Dy, cth, nxn, nyn, nzn);
  scale(Dz, cth, nxn, nyn, nzn);
  sum(imageX, Dx, nxn, nyn, nzn);
  sum(imageY, Dy, nxn, nyn, nzn);
  sum(imageZ, Dz, nxn, nyn, nzn);

  // (1)+(2)+(3)+D+E
  sum(imageX, vectX, nxn, nyn, nzn);
  sum(imageY, vectY, nxn, nyn, nzn);
  sum(imageZ, vectZ, nxn, nyn, nzn);

  communicateNodeBC(nxn, nyn, nzn, imageX, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, imageY, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, imageZ, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);


  if(col->getCase() !="BATSRUS"){
    // boundary condition: Xleft
    if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
      perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 0);
    // boundary condition: Xright
    if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
      perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 0);
    // boundary condition: Yleft
    if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
      perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 1);
    // boundary condition: Yright
    if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
      perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 1);
    // boundary condition: Zleft
    if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
      perfectConductorLeft(imageX, imageY, imageZ, vectX, vectY, vectZ, 2);
    // boundary condition: Zright
    if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
      perfectConductorRight(imageX, imageY, imageZ, vectX, vectY, vectZ, 2);

    // OpenBC
    OpenBoundaryInflowEImage(imageX, imageY, imageZ, vectX, vectY, vectZ, nxn, nyn, nzn);
  }

  // move from physical space to krylov space
  phys2solver(im, imageX, imageY, imageZ,inminsolve,inmaxsolve,jnminsolve,jnmaxsolve,knminsolve,knmaxsolve);
  
}

void EMfields3D::calc_gradDivE(arr3_double gradDivEx_G, 
			       arr3_double gradDivEy_G, 
			       arr3_double gradDivEz_G, 
			       arr3_double Ex_G, 
			       arr3_double Ey_G, 
			       arr3_double Ez_G,
			       double coefDiff){
  /* 
     output = (1- coefDiff)*grad(div1(E)) + coefDiff*grad(div2(E)), 
     where div1 is the common compact operator, while div2 is a 
     diffusive operator.      
   */

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  grid->divN2C(divC, Ex_G, Ey_G, Ez_G);
  if(coefDiff>0){
    grid->interpN2C(tempXC,Ex_G);
    grid->interpN2C(tempYC,Ey_G);
    grid->interpN2C(tempZC,Ez_G);
    communicateCenterBC(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct, this);
    communicateCenterBC(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct, this);
    communicateCenterBC(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct, this);
    grid->divC2C(tempC, tempXC, tempYC, tempZC);  

    double coef1, coef2; 
    coef1 = coefDiff;
    coef2 = 1- coef1;
    for (int i = 1; i < nxc - 1; i++)
      for (int j = 1; j < nyc - 1; j++)
	for (int k = 1; k < nzc - 1; k++) {
	  divC[i][j][k] = tempC[i][j][k]*coef1 + divC[i][j][k]*coef2; 
	}
  }
  
  communicateCenterBC(nxc, nyc, nzc, divC, 2, 2, 2, 2, 2, 2, vct, this);
  grid->gradC2N(gradDivEx_G, gradDivEy_G, gradDivEz_G, divC);  

  communicateNodeBC(nxn, nyn, nzn, gradDivEx_G, 2, 2, 2, 2, 2, 2, vct, this);
  communicateNodeBC(nxn, nyn, nzn, gradDivEy_G, 2, 2, 2, 2, 2, 2, vct, this);
  communicateNodeBC(nxn, nyn, nzn, gradDivEz_G, 2, 2, 2, 2, 2, 2, vct, this);

}


/*! Calculate PI dot (vectX, vectY, vectZ) */
void EMfields3D::PIdot(arr3_double PIdotX, arr3_double PIdotY, arr3_double PIdotZ, const_arr3_double vectX, const_arr3_double vectY, const_arr3_double vectZ, int ns)
{
  const Grid *grid = &get_grid();
  double beta, edotb, omcx, omcy, omcz, denom;

  beta = .5 * qom[ns] * dt / c;
  for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++) {
        omcx = beta * (Bxn[i][j][k] + Bx_ext[i][j][k]);
        omcy = beta * (Byn[i][j][k] + By_ext[i][j][k]);
        omcz = beta * (Bzn[i][j][k] + Bz_ext[i][j][k]);

        edotb = vectX.get(i,j,k) * omcx + vectY.get(i,j,k) * omcy + vectZ.get(i,j,k) * omcz;
        denom = 1 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        PIdotX.fetch(i,j,k) += (vectX.get(i,j,k) + (vectY.get(i,j,k) * omcz - vectZ.get(i,j,k) * omcy + edotb * omcx)) * denom;
        PIdotY.fetch(i,j,k) += (vectY.get(i,j,k) + (vectZ.get(i,j,k) * omcx - vectX.get(i,j,k) * omcz + edotb * omcy)) * denom;
        PIdotZ.fetch(i,j,k) += (vectZ.get(i,j,k) + (vectX.get(i,j,k) * omcy - vectY.get(i,j,k) * omcx + edotb * omcz)) * denom;
      }
}
/*! Calculate MU dot (vectX, vectY, vectZ) */
void EMfields3D::MUdot(arr3_double MUdotX, arr3_double MUdotY, arr3_double MUdotZ,
  arr3_double vectX, arr3_double vectY, arr3_double vectZ)
{
  const Grid *grid = &get_grid();
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  double beta, edotb, omcx, omcy, omcz, denom;

  if(col->getuseAccurateJ()){

    communicateNodeBC(nxn, nyn, nzn, vectX, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, vectY, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, vectZ, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);

    double c0 = FourPI*delt/c;
    
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	for (int k = 1; k < nzn - 1; k++) {
	  MUdotX[i][j][k] = 0.0;
	  MUdotY[i][j][k] = 0.0;
	  MUdotZ[i][j][k] = 0.0;	  

	  int gp; // g' 
	  for(int i2 = i-1; i2 <= i+1; i2++)
	    for(int j2 = j-1; j2 <= j+1; j2++)
	      for(int k2 = k-1; k2 <= k+1; k2++){
		gp = (i2-i+1)*9+(j2-j+1)*3+k2-k+1;
		const double *M_I = M_GII[i][j][k][gp];
		const double& vctX = vectX[i2][j2][k2];
		const double& vctY = vectY[i2][j2][k2];
		const double& vctZ = vectZ[i2][j2][k2];
		MUdotX[i][j][k] += (vctX*M_I[0] + vctY*M_I[1] + vctZ*M_I[2])*c0;
		MUdotY[i][j][k] += (vctX*M_I[3] + vctY*M_I[4] + vctZ*M_I[5])*c0; 
		MUdotZ[i][j][k] += (vctX*M_I[6] + vctY*M_I[7] + vctZ*M_I[8])*c0;
	      }
	}
    return;
  }

  
  for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
      for (int k = 1; k < nzn - 1; k++) {
        MUdotX[i][j][k] = 0.0;
        MUdotY[i][j][k] = 0.0;
        MUdotZ[i][j][k] = 0.0;
      }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int k = 1; k < nzn - 1; k++) {
          omcx = beta * (Bxn[i][j][k] + Bx_ext[i][j][k]);
          omcy = beta * (Byn[i][j][k] + By_ext[i][j][k]);
          omcz = beta * (Bzn[i][j][k] + Bz_ext[i][j][k]);	  
          edotb = vectX.get(i,j,k) * omcx + vectY.get(i,j,k) * omcy + vectZ.get(i,j,k) * omcz;
          denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
          MUdotX.fetch(i,j,k) += (vectX.get(i,j,k) + (vectY.get(i,j,k) * omcz - vectZ.get(i,j,k) * omcy + edotb * omcx)) * denom;
          MUdotY.fetch(i,j,k) += (vectY.get(i,j,k) + (vectZ.get(i,j,k) * omcx - vectX.get(i,j,k) * omcz + edotb * omcy)) * denom;
          MUdotZ.fetch(i,j,k) += (vectZ.get(i,j,k) + (vectX.get(i,j,k) * omcy - vectY.get(i,j,k) * omcx + edotb * omcz)) * denom;
        }
  }  
}

void EMfields3D::boxSmooth(arr3_double vector, int nx, int ny, int nz, int nSmooth)
{
  const VirtualTopology3D *vct = &get_vct();
  double ***temp = newArr3(double, nx, ny, nz);
  communicateNodeBC(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct, this);   
  for(int iSmooth = 0; iSmooth<nSmooth; iSmooth++){
    for (int i = 1; i < nx - 1; i++)
      for (int j = 1; j < ny - 1; j++)
	for (int k = 1; k < nz - 1; k++){
	  temp[i][j][k] = 0; 
	  for(int ii = -1; ii<2; ii++)
	    for(int jj=-1; jj<2; jj++)
	      for(int kk=-1; kk<2; kk++)
		temp[i][j][k] += vector[i+ii][j+jj][k+kk];
	  temp[i][j][k] *= 1./27;
	}
    for (int i = 1; i < nx - 1; i++)
      for (int j = 1; j < ny - 1; j++)
	for (int k = 1; k < nz - 1; k++)
	  vector[i][j][k] = temp[i][j][k];

    communicateNodeBC(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct, this);   
  }
  delArr3(temp, nx, ny);
}


/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smooth(arr3_double vector, int type)
{
  if(Smooth==1.0) return;
  const Collective *col = &get_col();
#ifdef BATSRUS
  if(!col->getdoSmoothAll()) return;
#endif  
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double alpha   = Smooth;
  double beta3D  = (1-alpha)/6.0;
  double beta2D  = (1-alpha)/4.0;
  int nx, ny, nz;
  switch (type) {
  	  case (0):
			   nx = grid->getNXC();
			   ny = grid->getNYC();
			   nz = grid->getNZC();
			   break;
  	  case (1):
			   nx = grid->getNXN();
			   ny = grid->getNYN();
			   nz = grid->getNZN();
			   break;
  }
  double ***temp = newArr3(double, nx, ny, nz);
  for (int icount = 1; icount < SmoothNiter + 1; icount++) {
      switch (type) {
        case (0):
          communicateCenterBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct, this);
          break;
        case (1):
          communicateNodeBoxStencilBC_P(nx, ny, nz, vector, 2, 2, 2, 2, 2, 2, vct, this);
          break;
      }

      for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++)
          for (int k = 1; k < nz - 1; k++){
	    // getSmoothFactor is node based. but i,j,k can be cell indexes 
#ifdef BATSRUS
	    alpha = col->getSmoothFactor(i,j,k);
#else
	    alpha = Smooth;
#endif
	    beta3D  = (1-alpha)/6.0;
	    
            temp[i][j][k] = alpha * vector[i][j][k] + beta3D * (vector[i - 1][j][k] + vector[i + 1][j][k] + vector[i][j - 1][k] + vector[i][j + 1][k] + vector[i][j][k - 1] + vector[i][j][k + 1]);
          //temp[i][j][k] = alpha * vector[i][j][k] + beta2D * (vector[i - 1][j][k] + vector[i + 1][j][k] + vector[i][j][k - 1] + vector[i][j][k + 1]);
  }

      for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++)
          for (int k = 1; k < nz - 1; k++)
            vector[i][j][k] = temp[i][j][k];

  }
  delArr3(temp, nx, ny);
}

/* Interpolation smoothing: Smoothing (vector must already have ghost cells)
 * TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smoothE()
{
  if(Smooth==1.0) return;
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();

  double alpha   = Smooth;
  double beta3D  = (1-alpha)/6.0;
  double beta2D  = (1-alpha)/4.0;

  double ***temp = newArr3(double, nxn, nyn, nzn);

  for (int icount = 1; icount < SmoothNiter + 1; icount++) {
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);
      
      // Exth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++){
#ifdef BATSRUS
	    alpha = col->getSmoothFactor(i,j,k);
#else
	    alpha = Smooth;
#endif
	    beta3D  = (1-alpha)/6.0;
            temp[i][j][k] = alpha * Exth[i][j][k] + beta3D * (Exth[i - 1][j][k] + Exth[i + 1][j][k] + Exth[i][j - 1][k] + Exth[i][j + 1][k] + Exth[i][j][k - 1] + Exth[i][j][k + 1]);
      	    //temp[i][j][k] = alpha * Ex[i][j][k] + beta2D * (Ex[i - 1][j][k] + Ex[i + 1][j][k] + Ex[i][j][k - 1] + Ex[i][j][k + 1]);
	  }

      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Exth[i][j][k] = temp[i][j][k];
      // Eyth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++){
#ifdef BATSRUS
	    alpha = col->getSmoothFactor(i,j,k);
#else
	    alpha = Smooth;
#endif
	    beta3D  = (1-alpha)/6.0;

            temp[i][j][k] = alpha * Eyth[i][j][k] + beta3D * (Eyth[i - 1][j][k] + Eyth[i + 1][j][k] + Eyth[i][j - 1][k] + Eyth[i][j + 1][k] + Eyth[i][j][k - 1] + Eyth[i][j][k + 1]);
      //temp[i][j][k] = alpha * Ey[i][j][k] + beta2D * (Ey[i - 1][j][k] + Ey[i + 1][j][k] + Ey[i][j][k - 1] + Ey[i][j][k + 1]);
	  }

      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Eyth[i][j][k] = temp[i][j][k];
      // Ezth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++){
#ifdef BATSRUS
	    alpha = col->getSmoothFactor(i,j,k);
#else
	    alpha = Smooth;
#endif
	    beta3D  = (1-alpha)/6.0;

            temp[i][j][k] = alpha * Ezth[i][j][k] + beta3D * (Ezth[i - 1][j][k] + Ezth[i + 1][j][k] + Ezth[i][j - 1][k] + Ezth[i][j + 1][k] + Ezth[i][j][k - 1] + Ezth[i][j][k + 1]);
      //temp[i][j][k] = alpha * Ez[i][j][k] + beta2D * (Ez[i - 1][j][k] + Ez[i + 1][j][k] + Ez[i][j][k - 1] + Ez[i][j][k + 1]);
	  }

      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Ezth[i][j][k] = temp[i][j][k];
      
  }


  if(col->getCase()=="BATSRUS"){
    communicateNodeBoxStencilBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
    communicateNodeBoxStencilBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
    communicateNodeBoxStencilBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);
  }
  
  delArr3(temp, nxn, nyn);  
}

/* SPECIES: Interpolation smoothing TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector type = 1 --> node based vector */
void EMfields3D::smooth(double value, arr4_double vector, int is, int type)
{
  eprintf("Smoothing for Species not implemented in 3D");
}

/*! fix the B boundary when running gem , This assume non-periodic condition on Y dimension*/
void EMfields3D::fixBcGEM()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Bxc[i][nyc - 2][k] = Bxc[i][nyc - 1][k];
        Bxc[i][nyc - 3][k] = Bxc[i][nyc - 1][k];
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z;
        Bzc[i][nyc - 2][k] = B0z;
        Bzc[i][nyc - 3][k] = B0z;
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bxc[i][1][k] = Bxc[i][0][k];
        Bxc[i][2][k] = Bxc[i][0][k];
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z;
        Bzc[i][1][k] = B0z;
        Bzc[i][2][k] = B0z;
      }
  }
}

void EMfields3D::fixBnGEM()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxn[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Bxn[i][nyc - 2][k] = Bxc[i][nyc - 1][k];
        Bxn[i][nyc - 3][k] = Bxc[i][nyc - 1][k];
        Byn[i][nyc - 1][k] = B0y;
        Bzn[i][nyc - 1][k] = B0z;
        Bzn[i][nyc - 2][k] = B0z;
        Bzn[i][nyc - 3][k] = B0z;
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bxc[i][1][k] = Bxc[i][0][k];
        Bxc[i][2][k] = Bxc[i][0][k];
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z;
        Bzc[i][1][k] = B0z;
        Bzc[i][2][k] = B0z;
      }
  }
}

/*! fix the B boundary when running forcefree */
void EMfields3D::fixBforcefree()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][nyc - 1][k] = B0x * tanh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);
        Byc[i][nyc - 1][k] = B0y;
        Bzc[i][nyc - 1][k] = B0z / cosh((grid->getYC(i, nyc - 1, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 2][k] = B0z / cosh((grid->getYC(i, nyc - 2, k) - Ly / 2) / delta);;
        Bzc[i][nyc - 3][k] = B0z / cosh((grid->getYC(i, nyc - 3, k) - Ly / 2) / delta);
      }
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    for (int i = 0; i < nxc; i++)
      for (int k = 0; k < nzc; k++) {
        Bxc[i][0][k] = B0x * tanh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Byc[i][0][k] = B0y;
        Bzc[i][0][k] = B0z / cosh((grid->getYC(i, 0, k) - Ly / 2) / delta);
        Bzc[i][1][k] = B0z / cosh((grid->getYC(i, 1, k) - Ly / 2) / delta);
        Bzc[i][2][k] = B0z / cosh((grid->getYC(i, 2, k) - Ly / 2) / delta);
      }
  }
}

// This method assumes mirror boundary conditions;
// we therefore need to double the density on the boundary
// nodes to incorporate the mirror particles from the mirror
// cell just outside the domain.
//
/*! adjust densities on boundaries that are not periodic */
void EMfields3D::adjustNonPeriodicDensities(int is)
{
  const VirtualTopology3D *vct = &get_vct();
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][1][i][k] *= 2;
        Jxs  [is][1][i][k] *= 2;
        Jys  [is][1][i][k] *= 2;
        Jzs  [is][1][i][k] *= 2;
        pXXsn[is][1][i][k] *= 2;
        pXYsn[is][1][i][k] *= 2;
        pXZsn[is][1][i][k] *= 2;
        pYYsn[is][1][i][k] *= 2;
        pYZsn[is][1][i][k] *= 2;
        pZZsn[is][1][i][k] *= 2;
    }
  }
  
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][i][1][k] *= 2;
        Jxs  [is][i][1][k] *= 2;
        Jys  [is][i][1][k] *= 2;
        Jzs  [is][i][1][k] *= 2;
        pXXsn[is][i][1][k] *= 2;
        pXYsn[is][i][1][k] *= 2;
        pXZsn[is][i][1][k] *= 2;
        pYYsn[is][i][1][k] *= 2;
        pYZsn[is][i][1][k] *= 2;
        pZZsn[is][i][1][k] *= 2;
    }
  }
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
    {
        rhons[is][i][j][1] *= 2;
        Jxs  [is][i][j][1] *= 2;
        Jys  [is][i][j][1] *= 2;
        Jzs  [is][i][j][1] *= 2;
        pXXsn[is][i][j][1] *= 2;
        pXYsn[is][i][j][1] *= 2;
        pXZsn[is][i][j][1] *= 2;
        pYYsn[is][i][j][1] *= 2;
        pYZsn[is][i][j][1] *= 2;
        pZZsn[is][i][j][1] *= 2;
    }
  }
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][nxn - 2][i][k] *= 2;
        Jxs  [is][nxn - 2][i][k] *= 2;
        Jys  [is][nxn - 2][i][k] *= 2;
        Jzs  [is][nxn - 2][i][k] *= 2;
        pXXsn[is][nxn - 2][i][k] *= 2;
        pXYsn[is][nxn - 2][i][k] *= 2;
        pXZsn[is][nxn - 2][i][k] *= 2;
        pYYsn[is][nxn - 2][i][k] *= 2;
        pYZsn[is][nxn - 2][i][k] *= 2;
        pZZsn[is][nxn - 2][i][k] *= 2;
    }
  }
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int k = 1; k < nzn - 1; k++)
    {
        rhons[is][i][nyn - 2][k] *= 2;
        Jxs  [is][i][nyn - 2][k] *= 2;
        Jys  [is][i][nyn - 2][k] *= 2;
        Jzs  [is][i][nyn - 2][k] *= 2;
        pXXsn[is][i][nyn - 2][k] *= 2;
        pXYsn[is][i][nyn - 2][k] *= 2;
        pXZsn[is][i][nyn - 2][k] *= 2;
        pYYsn[is][i][nyn - 2][k] *= 2;
        pYZsn[is][i][nyn - 2][k] *= 2;
        pZZsn[is][i][nyn - 2][k] *= 2;
    }
  }
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
    for (int j = 1; j < nyn - 1; j++)
    {
        rhons[is][i][j][nzn - 2] *= 2;
        Jxs  [is][i][j][nzn - 2] *= 2;
        Jys  [is][i][j][nzn - 2] *= 2;
        Jzs  [is][i][j][nzn - 2] *= 2;
        pXXsn[is][i][j][nzn - 2] *= 2;
        pXYsn[is][i][j][nzn - 2] *= 2;
        pXZsn[is][i][j][nzn - 2] *= 2;
        pYYsn[is][i][j][nzn - 2] *= 2;
        pYZsn[is][i][j][nzn - 2] *= 2;
        pZZsn[is][i][j][nzn - 2] *= 2;
    }
  }

}

void EMfields3D::ConstantChargeOpenBCv2()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

	ff = qom[is]/fabs(qom[is]);

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = rhons[is][4][j][k];
          rhons[is][1][j][k] = rhons[is][4][j][k];
          rhons[is][2][j][k] = rhons[is][4][j][k];
          rhons[is][3][j][k] = rhons[is][4][j][k];
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-3][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-2][j][k] = rhons[is][nx-5][j][k];
          rhons[is][nx-1][j][k] = rhons[is][nx-5][j][k];
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = rhons[is][i][4][k];
          rhons[is][i][1][k] = rhons[is][i][4][k];
          rhons[is][i][2][k] = rhons[is][i][4][k];
          rhons[is][i][3][k] = rhons[is][i][4][k];
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-3][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-2][k] = rhons[is][i][ny-5][k];
          rhons[is][i][ny-1][k] = rhons[is][i][ny-5][k];
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = rhons[is][i][j][4];
          rhons[is][i][j][1] = rhons[is][i][j][4];
          rhons[is][i][j][2] = rhons[is][i][j][4];
          rhons[is][i][j][3] = rhons[is][i][j][4];
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-3] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-2] = rhons[is][i][j][nz-5];
          rhons[is][i][j][nz-1] = rhons[is][i][j][nz-5];
        }
    }
  }

}

void EMfields3D::ConstantChargeOpenBC()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double ff;

  int nx = grid->getNXN();
  int ny = grid->getNYN();
  int nz = grid->getNZN();

  for (int is = 0; is < ns; is++) {

    ff = qom[is]/fabs(qom[is]);

    if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcEMfaceXleft ==2)) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][0][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][1][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][3][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcEMfaceXright ==2)) {
      for (int j=0; j < ny;j++)
        for (int k=0; k < nz;k++){
          rhons[is][nx-4][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-3][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-2][j][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][nx-1][j][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcEMfaceYleft ==2))  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][0][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][1][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][3][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcEMfaceYright ==2))  {
      for (int i=0; i < nx;i++)
        for (int k=0; k < nz;k++){
          rhons[is][i][ny-4][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-3][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-2][k] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][ny-1][k] = ff * rhoINIT[is] / FourPI;
        }
    }

    if(vct->getZleft_neighbor()==MPI_PROC_NULL && (bcEMfaceZleft ==2))  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][0] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][1] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][3] = ff * rhoINIT[is] / FourPI;
        }
    }


    if(vct->getZright_neighbor()==MPI_PROC_NULL && (bcEMfaceZright ==2))  {
      for (int i=0; i < nx;i++)
        for (int j=0; j < ny;j++){
          rhons[is][i][j][nz-4] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-3] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-2] = ff * rhoINIT[is] / FourPI;
          rhons[is][i][j][nz-1] = ff * rhoINIT[is] / FourPI;
        }
    }
  }

}

void EMfields3D::ConstantChargePlanet(double R,
  double x_center, double y_center, double z_center)
{
  const Grid *grid = &get_grid();

  double xd;
  double yd;
  double zd;
  double ff;

  for (int is = 0; is < ns; is++) {

	ff = qom[is]/fabs(qom[is]);

    for (int i = 1; i < nxn; i++) {
      for (int j = 1; j < nyn; j++) {
        for (int k = 1; k < nzn; k++) {

          xd = grid->getXN(i,j,k) - x_center;
          yd = grid->getYN(i,j,k) - y_center;
          zd = grid->getZN(i,j,k) - z_center;

          if ((xd*xd+yd*yd+zd*zd) <= R*R) {
            rhons[is][i][j][k] = ff * rhoINIT[is] / FourPI;
          }

        }
      }
    }
  }

}

void EMfields3D::ConstantChargePlanet2DPlaneXZ(double R,  double x_center,double z_center)
{
  const Grid *grid = &get_grid();
  //if (get_vct().getCartesian_rank() == 0)
      //cout << "*** Constant Charge 2D Planet ***" << endl;

  assert_eq(nyn,4);
  double xd;
  double zd;

  for (int is = 0; is < ns; is++) {
    const double sign_q = qom[is]/(fabs(qom[is]));
    for (int i = 1; i < nxn; i++)
        for (int k = 1; k < nzn; k++) {

          xd = grid->getXN(i,1,k) - x_center;
          zd = grid->getZN(i,1,k) - z_center;

          if ((xd*xd+zd*zd) <= R*R) {
            rhons[is][i][1][k] = sign_q * rhoINIT[is] / FourPI;
            rhons[is][i][2][k] = sign_q * rhoINIT[is] / FourPI;
          }

	}
  }

}

/*! Populate the field data used to push particles */
// 
// 
//
void EMfields3D::set_fieldForPcls()
{
  #pragma omp parallel for collapse(3)
  for(int i=0;i<nxn;i++)
  for(int j=0;j<nyn;j++)
  for(int k=0;k<nzn;k++)
  {
    fieldForPcls[i][j][k][0] = (pfloat) (Bxn[i][j][k] + Bx_ext[i][j][k]);
    fieldForPcls[i][j][k][1] = (pfloat) (Byn[i][j][k] + By_ext[i][j][k]);
    fieldForPcls[i][j][k][2] = (pfloat) (Bzn[i][j][k] + Bz_ext[i][j][k]);
    fieldForPcls[i][j][k][0+DFIELD_3or4] = (pfloat) Exthp[i][j][k];
    fieldForPcls[i][j][k][1+DFIELD_3or4] = (pfloat) Eythp[i][j][k];
    fieldForPcls[i][j][k][2+DFIELD_3or4] = (pfloat) Ezthp[i][j][k];
  }
}

/*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
void EMfields3D::calculateB()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getCartesian_rank() == 0)
    cout << "*** B CALCULATION ***" << endl;

  // calculate the curl of Eth
  grid->curlN2C(tempXC, tempYC, tempZC, Exth, Eyth, Ezth);  
  
  // update the magnetic field
  addscale(-c * dt, 1, Bxc, tempXC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Byc, tempYC, nxc, nyc, nzc);
  addscale(-c * dt, 1, Bzc, tempZC, nxc, nyc, nzc);

  // communicate ghost 
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct,this);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct,this);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct,this);

  if (get_col().getCase()=="ForceFree") 	fixBforcefree();
  if (get_col().getCase()=="GEM")       	fixBnGEM();
  if (get_col().getCase()=="GEMnoPert") 	fixBnGEM();
  if (get_col().getCase()=="GEMDoubleHarris") 	fixBnGEM();
#ifdef BATSRUS
  if (get_col().getCase()=="BATSRUS")           fixB_BATSRUS();
  #endif
  
  // OpenBC:
  if (get_col().getCase()!="BATSRUS")
    OpenBoundaryInflowB(Bxc,Byc,Bzc,nxc,nyc,nzc);

  // interpolate C2N
  grid->interpC2N(Bxn, Bxc);
  grid->interpC2N(Byn, Byc);
  grid->interpC2N(Bzn, Bzc);

  if (get_col().getCase()=="GEM")       	fixBcGEM();
  if (get_col().getCase()=="GEMnoPert") 	fixBcGEM();
  if (get_col().getCase()=="GEMDoubleHarris") 	fixBcGEM();

  communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct, this);

#ifdef BATSRUS
  if (get_col().getCase()=="BATSRUS") fixB_BATSRUS();
#endif  

}


/*
void EMfields3D::calc_energy_error(){
  const Collective *col = &get_col();
  bool doFixEnergy = col->get_doFixEnergy();
  if(!doFixEnergy) return;   

  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  const double invFourPi = 1./FourPI; 
  
  eqValue(0.0, energyError_G, nxn, nyn,nzn); 

  if(false){
    MUdot(Dx, Dy, Dz, Exthp, Eythp, Ezthp);
    // Coeficient "4*pi*delt/c" is multiplied to M*E inside MUdot. Undo the multiplication.
    double c1 = 1/(FourPI*delt/c);
    scale(Dx, c1, nxn, nyn, nzn);
    scale(Dy, c1, nxn, nyn, nzn);
    scale(Dz, c1, nxn, nyn, nzn);
   
    // Calculate the energy error on each node. 
    for (int i = 1; i < nxn-1; i++)
      for (int j = 1; j < nyn-1; j++)
	for(int k = 1; k < nzn-1; k++){
	  energyError_G[i][j][k] =
	    -dt*((Jxh[i][j][k] + Dx[i][j][k])*(Exthp[i][j][k] - Exth[i][j][k]) + 
		 (Jyh[i][j][k] + Dy[i][j][k])*(Eythp[i][j][k] - Eyth[i][j][k])+ 
		 (Jzh[i][j][k] + Dz[i][j][k])*(Ezthp[i][j][k] - Ezth[i][j][k]));
	}     

  }

  double coefDiffusion = get_col().get_ratioDivC2C();
  
  if(coefDiffusion != 0){  
    calc_gradDivE(tempXN, tempYN, tempZN, Exth, Eyth, Ezth, coefDiffusion);
    calc_gradDivE(tempX, tempY, tempZ, Exth, Eyth, Ezth, 0);
    for (int i = 1; i < nxn-1; i++)
      for (int j = 1; j < nyn-1; j++)
	for(int k = 1; k < nzn-1; k++){
	  energyError_G[i][j][k] += 2*invFourPi*delt*delt*
	    ((tempXN[i][j][k]-tempX[i][j][k])*Exth[i][j][k] + 
	     (tempYN[i][j][k]-tempY[i][j][k])*Eyth[i][j][k] + 
	     (tempZN[i][j][k]-tempZ[i][j][k])*Ezth[i][j][k]); 
	}
  }

  // Smooth the energy errors. 
  int nSmooth = col->get_nSmoothEnergy();
  boxSmooth(energyError_G, nxn, nyn, nzn, nSmooth);
}

void EMfields3D::fix_energy(){
  const Collective *col = &get_col();
  bool doFixEnergy = col->get_doFixEnergy();

  if(!doFixEnergy) return;   
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  const double invFourPi = 1./FourPI; 
  
  // Compress or stretch the electric field so that the energy is conserved. 
  double energyE, energyError, energyENew; 
  double scaleRatio, scaleRatio1, scaleRatio2; 
  //double Estatx, Estaty, Estatz, Ewavex, Ewavey, Ewavez; 
  for (int i = 1; i < nxn-1; i++)
    for (int j = 1; j < nyn-1; j++)
      for(int k = 1; k < nzn-1; k++){
	energyE = 0.5*invFourPi*(Ex[i][j][k] * Ex[i][j][k] + 
				 Ey[i][j][k] * Ey[i][j][k] + 
				 Ez[i][j][k] * Ez[i][j][k]);
	 
	energyError = energyError_G[i][j][k]; 
	energyENew = energyError + energyE;

	if(energyENew >= 0){ 
	  scaleRatio = sqrt((energyENew)/energyE); 
	  Ex[i][j][k] *=scaleRatio; 
	  Ey[i][j][k] *=scaleRatio; 
	  Ez[i][j][k] *=scaleRatio; 
	}	
      }   

  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);
}

*/

/*! initialize EM field with transverse electric waves 1D and rotate anticlockwise (theta degrees) */
void EMfields3D::initEM_rotate(double B, double theta)
{
  const Grid *grid = &get_grid();

  // initialize E and rhos on nodes
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      Ex[i][j][0] = 0.0;
      Ey[i][j][0] = 0.0;
      Ez[i][j][0] = 0.0;
      Bxn[i][j][0] = B * cos(theta * M_PI / 180);
      Byn[i][j][0] = B * sin(theta * M_PI / 180);
      Bzn[i][j][0] = 0.0;
      rhons[0][i][j][0] = 0.07957747154595; // electrons: species is now first index
      rhons[1][i][j][0] = 0.07957747154595; // protons: species is now first index
    }
  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);


  for (int is = 0; is < ns; is++)
    grid->interpN2C(rhocs, is, rhons);
}

/*!Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid) {

  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  ne_mod *= alpha;
  ni_mod *= alpha;
  // cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      rhons[0][i][j][0] += ne_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ne_phase);
      rhons[1][i][j][0] += ni_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ni_phase);
    }

  for (int is = 0; is < ns; is++)
    grid->interpN2C(rhocs, is, rhons);
}


/*!Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid) {

  double alpha;

  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;

  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      Ex[i][j][0] += Ex_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ex_phase);
      Ey[i][j][0] += Ey_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ey_phase);
      Ez[i][j][0] += Ez_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ez_phase);
      Bxn[i][j][0] += Bx_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bx_phase);
      Byn[i][j][0] += By_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + By_phase);
      Bzn[i][j][0] += Bz_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bz_phase);

    }

  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);

}


/*! Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
void EMfields3D::calculateHatFunctions()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  const Collective *col = &get_col();
  // smoothing
  smooth(rhoc,0);
    
#ifdef BATSRUS
  // calculate j hat
  if(col->getdoTestEMWave()){
    eqValue(0, Jxh, nxn, nyn, nzn);
    eqValue(0, Jyh, nxn, nyn, nzn);
    eqValue(0, Jzh, nxn, nyn, nzn);
    eqValue(0, rhoh, nxc, nyc, nzc);
    return;
  }
  
  if(get_col().getCase()=="BATSRUS")
    for (int is=0; is < ns; is++){

      if(get_col().getuseAccurateJ() && col->getnIsotropic()>=0){
	eprintf("Error: Rotation of Jxsh/Jysh/Jzsh at boundary has not been implemented!!!");
      }
      
      // Pass Jxs[is] is better!! -Yuxi      
      fixVarBCnode(Jxs,&Collective::getPICJx<int>,col->getnIsotropic(),is);
      fixVarBCnode(Jys,&Collective::getPICJy<int>,col->getnIsotropic(),is);
      fixVarBCnode(Jzs,&Collective::getPICJz<int>,col->getnIsotropic(),is);
      fixVarBCnode(pXXsn,&Collective::getPICPxx<int>,col->getnIsotropic(),is);      
      fixVarBCnode(pXYsn,&Collective::getPICPxy<int>,col->getnIsotropic(),is);
      fixVarBCnode(pXZsn,&Collective::getPICPxz<int>,col->getnIsotropic(),is);
      fixVarBCnode(pYYsn,&Collective::getPICPyy<int>,col->getnIsotropic(),is);
      fixVarBCnode(pYZsn,&Collective::getPICPyz<int>,col->getnIsotropic(),is);
      fixVarBCnode(pZZsn,&Collective::getPICPzz<int>,col->getnIsotropic(),is);
    }
#endif

  if(col->getuseAccurateJ()){
    eqValue(0, Jxh, nxn, nyn, nzn);
    eqValue(0, Jyh, nxn, nyn, nzn);
    eqValue(0, Jzh, nxn, nyn, nzn);

    for(int is = 0; is < ns; is++){
      sum(Jxh, Jxsh, nxn, nyn, nzn, is);
      sum(Jyh, Jysh, nxn, nyn, nzn, is);
      sum(Jzh, Jzsh, nxn, nyn, nzn, is);
    }
  }else{
    for (int is=0; is < ns; is++){
      grid->divSymmTensorN2C(tempXC, tempYC, tempZC, pXXsn, pXYsn, pXZsn, pYYsn, pYZsn, pZZsn, is);

      scale(tempXC, -dt / 2.0, nxc, nyc, nzc);
      scale(tempYC, -dt / 2.0, nxc, nyc, nzc);
      scale(tempZC, -dt / 2.0, nxc, nyc, nzc);
      // communicate before interpolating
      communicateCenterBC_P(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct, this);
      communicateCenterBC_P(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct, this);
      communicateCenterBC_P(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct, this);

      grid->interpC2N(tempXN, tempXC);
      grid->interpC2N(tempYN, tempYC);
      grid->interpC2N(tempZN, tempZC);

      if(col->getuseExplicitMover()){
	eqValue(0.0, tempXN, nxn, nyn, nzn);
	eqValue(0.0, tempYN, nxn, nyn, nzn);
	eqValue(0.0, tempZN, nxn, nyn, nzn);
      }
      
      sum(tempXN, Jxs, nxn, nyn, nzn, is);
      sum(tempYN, Jys, nxn, nyn, nzn, is);
      sum(tempZN, Jzs, nxn, nyn, nzn, is);

      // PIDOT
      PIdot(Jxh, Jyh, Jzh, tempXN, tempYN, tempZN, is);

    }
  }
  
  // smooth j
  smooth(Jxh, 1);
  smooth(Jyh, 1);
  smooth(Jzh, 1);

  // calculate rho hat = rho - (dt*theta)div(jhat)
  grid->divN2C(tempXC, Jxh, Jyh, Jzh);
  scale(tempXC, -dt * th, nxc, nyc, nzc);
  sum(tempXC, rhoc, nxc, nyc, nzc);
  eq(rhoh, tempXC, nxc, nyc, nzc);

  // communicate rhoh
  communicateCenterBC_P(nxc, nyc, nzc, rhoh, 2, 2, 2, 2, 2, 2, vct, this);
  #ifdef BATSRUS
  if(get_col().getCase()=="BATSRUS")
    fixVarBCcell(rhoh,&Collective::getFluidZero<int>,col->getnCharge(),0);
  #endif
}
/*! Image of Poisson Solver */
void EMfields3D::PoissonImage(double *image, double *vector, bool doSolveForChange)
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  // allocate 2 three dimensional service vectors
  array3_double temp(nxc, nyc, nzc);
  array3_double im(nxc, nyc, nzc);
  eqValue(0.0, image, nSolveCell);
  eqValue(0.0, temp, nxc, nyc, nzc);
  eqValue(0.0, im, nxc, nyc, nzc);
  // move from krylov space to physical space and communicate ghost cells
  solver2phys(temp, vector, icMinSolve,icMaxSolve,
	      jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);
  // calculate the laplacian
  grid->lapC2Cpoisson(im, temp, this);
  // move from physical space to krylov space
  phys2solver(image, im, icMinSolve, icMaxSolve,
	      jcMinSolve,jcMaxSolve,kcMinSolve,kcMaxSolve);
}
/*! interpolate charge density and pressure density from node to center */
void EMfields3D::interpDensitiesN2C()
{
  const VirtualTopology3D *vct = &get_vct();

  get_grid().interpN2C(rhoc, rhon);
  communicateCenterBC_P(nxc, nyc, nzc, rhoc, 2, 2, 2, 2, 2, 2, vct, this);

  for (int is = 0; is < ns; is++){
    get_grid().interpN2C(rhocs, is, rhons);
    double ***moment = convert_to_arr3(rhocs[is]);
    communicateCenterBC_P(nxc, nyc, nzc, moment, 2, 2, 2, 2, 2, 2, vct, this);
  }
}
/*! communicate ghost for grid -> Particles interpolation */
void EMfields3D::communicateGhostP2G(int ns)
{
  // interpolate adding common nodes among processors
  timeTasks_set_communicating();

  const VirtualTopology3D *vct = &get_vct();

  double ***moment0 = convert_to_arr3(rhons[ns]);
  double ***moment1 = convert_to_arr3(Jxs  [ns]);
  double ***moment2 = convert_to_arr3(Jys  [ns]);
  double ***moment3 = convert_to_arr3(Jzs  [ns]);
  double ***moment4 = convert_to_arr3(pXXsn[ns]);
  double ***moment5 = convert_to_arr3(pXYsn[ns]);
  double ***moment6 = convert_to_arr3(pXZsn[ns]);
  double ***moment7 = convert_to_arr3(pYYsn[ns]);
  double ***moment8 = convert_to_arr3(pYZsn[ns]);
  double ***moment9 = convert_to_arr3(pZZsn[ns]);
  double ***moment10;
  double ***moment11;
  double ***moment12;

  // add the values for the shared nodes
  // Call NonBlocking Halo Exchange + Interpolation  
  communicateInterp(nxn, nyn, nzn, moment0, vct,this);
  communicateInterp(nxn, nyn, nzn, moment1, vct,this);
  communicateInterp(nxn, nyn, nzn, moment2, vct,this);
  communicateInterp(nxn, nyn, nzn, moment3, vct,this);
  communicateInterp(nxn, nyn, nzn, moment4, vct,this);
  communicateInterp(nxn, nyn, nzn, moment5, vct,this);
  communicateInterp(nxn, nyn, nzn, moment6, vct,this);
  communicateInterp(nxn, nyn, nzn, moment7, vct,this);
  communicateInterp(nxn, nyn, nzn, moment8, vct,this);
  communicateInterp(nxn, nyn, nzn, moment9, vct,this);

  if(get_col().getuseAccurateJ()){
    moment10 = convert_to_arr3(Jxsh[ns]);
    moment11 = convert_to_arr3(Jysh[ns]);
    moment12 = convert_to_arr3(Jzsh[ns]);
    
    communicateInterp(nxn, nyn, nzn, moment10, vct,this);
    communicateInterp(nxn, nyn, nzn, moment11, vct,this);
    communicateInterp(nxn, nyn, nzn, moment12, vct,this);
  }
  
  if(get_col().getCase()=="BATSRUS"){
    // setNodePlasmaBC is already called at the beginning of each iteration
    // to set BC for IPIC3D. Here, at the end of iteration, it is called
    // again to correcto the boundary nodes so that the output and the
    // hat values (see calculatehatfunctions()) can have the right values. 
#ifdef BATSRUS
    setNodePlasmaBC(ns);
#endif
  }else{
    // calculate the correct densities on the boundaries
    adjustNonPeriodicDensities(ns);
  }
  //Call Nonblocking Halo Exchange
  communicateNode_P(nxn, nyn, nzn, moment0, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment1, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment2, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment3, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment4, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment5, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment6, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment7, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment8, vct, this);
  communicateNode_P(nxn, nyn, nzn, moment9, vct, this);    

  if(get_col().getuseAccurateJ()){
    communicateNode_P(nxn, nyn, nzn, moment10, vct, this);
    communicateNode_P(nxn, nyn, nzn, moment11, vct, this);
    communicateNode_P(nxn, nyn, nzn, moment12, vct, this);    
  }
  
}

/*! communicate ghost for grid -> Particles interpolation */
//void EMfields3D::communicateGhostMomentsX()
//{
//  const VirtualTopology3D * vct = &get_vct();
//  timeTasks_set_communicating();
//  
//  // start communication of moments in the X direction
//  for(int is=0; is<ns; i++)
//  {
//    for(int im=0; im<10; im++)
//    {
//      // copy data from faces into buffers
//      // send buffer data right and left
//    }
//  }
//
//  // receive and parse communication
//}

void EMfields3D::setZeroDerivedMoments()
{
  for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
        Jx  [i][j][k] = 0.0;
        Jxh [i][j][k] = 0.0;
        Jy  [i][j][k] = 0.0;
        Jyh [i][j][k] = 0.0;
        Jz  [i][j][k] = 0.0;
        Jzh [i][j][k] = 0.0;
        rhon[i][j][k] = 0.0;
      }
  for (register int i = 0; i < nxc; i++)
    for (register int j = 0; j < nyc; j++)
      for (register int k = 0; k < nzc; k++) {
        rhoc[i][j][k] = 0.0;
        rhoh[i][j][k] = 0.0;
      }
}

void EMfields3D::setZeroPrimaryMoments() {

  // set primary moments to zero
  //
  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          rhons[kk][i][j][k] = 0.0;
          Jxs  [kk][i][j][k] = 0.0;
          Jys  [kk][i][j][k] = 0.0;
          Jzs  [kk][i][j][k] = 0.0;
          pXXsn[kk][i][j][k] = 0.0;
          pXYsn[kk][i][j][k] = 0.0;
          pXZsn[kk][i][j][k] = 0.0;
          pYYsn[kk][i][j][k] = 0.0;
          pYZsn[kk][i][j][k] = 0.0;
          pZZsn[kk][i][j][k] = 0.0;
	  Jxsh [kk][i][j][k] = 0.0;
          Jysh [kk][i][j][k] = 0.0;
          Jzsh [kk][i][j][k] = 0.0;
        }

}
/*! set to 0 all the densities fields */
void EMfields3D::setZeroDensities() {
  setZeroDerivedMoments();
  setZeroPrimaryMoments();
}

/*!SPECIES: Sum the charge density of different species on NODES */
void EMfields3D::sumOverSpecies()
{
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++)
          rhon[i][j][k] += rhons[is][i][j][k];
}

/*!SPECIES: Sum current density for different species */
void EMfields3D::sumOverSpeciesJ() {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          Jx[i][j][k] += Jxs[is][i][j][k];
          Jy[i][j][k] += Jys[is][i][j][k];
          Jz[i][j][k] += Jzs[is][i][j][k];
        }
}



/*! initialize Magnetic and Electric Field with initial configuration */
void EMfields3D::init()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++) {
      for (int j = 0; j < nyn; j++) {
        for (int k = 0; k < nzn; k++) {
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = B0x;
          Byn[i][j][k] = B0y;
          Bzn[i][j][k] = B0z;
        }
      }
    }

    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {                        // READING FROM RESTART
  #ifdef NO_HDF5
    eprintf("restart requires compiling with HDF5");
  #else
    col->read_field_restart(vct,grid,Bxn,Byn,Bzn,Bxc,Byc,Bzc,Ex,Ey,Ez,&rhons,ns);

    // communicate species densities to ghost nodes
    for (int is = 0; is < ns; is++) {
      double ***moment0 = convert_to_arr3(rhons[is]);
      communicateNode_P(nxn, nyn, nzn, moment0, vct, this);
    }

    if (col->getCase()=="Dipole") {
      ConstantChargePlanet(col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
    }else if(col->getCase()=="Dipole2D") {
    	ConstantChargePlanet2DPlaneXZ(col->getL_square(),col->getx_center(),col->getz_center());
      }

    ConstantChargeOpenBC();

    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct, this);

#ifndef BATSRUS
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
#endif

    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct,this);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct,this);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct,this);

    // communicate E
    communicateNodeBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct, this);

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  #endif // NO_HDF5
  }
}


#ifdef BATSRUS
/*! initiliaze EM for GEM challange */
void EMfields3D::initBATSRUS(VirtualTopology3D * vct, Grid * grid,
			     Collective * col)
{
  if (vct->getCartesian_rank() == 0) {
    cout << "------------------------------------------" << endl;
    cout << "         Initialize EMfield from BATSRUS          " << endl;
    cout << "------------------------------------------" << endl;
  }

  if (vct->getXleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    inminsolve = 1+NG_F;
    icMinSolve = 1+NG_F;
  }
  if (vct->getXright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    inmaxsolve = nxn-2-NG_F;
    icMaxSolve = nxc-2-NG_F;
  }
  if (vct->getYleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    jnminsolve = 1+NG_F;
    jcMinSolve = 1+NG_F;
  }
  if (vct->getYright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    jnmaxsolve = nyn-2-NG_F;
    jcMaxSolve = nyc-2-NG_F;
  }
  if (vct->getZleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    knminsolve = 1+NG_F;
    kcMinSolve = 1+NG_F;
  }
  if (vct->getZright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    knmaxsolve = nzn-2-NG_F;
    kcMaxSolve = nzc-2-NG_F;
  }

  /** number of unknowns solved for in implicit solver*/
  nSolveCell  = (icMaxSolve - icMinSolve + 1)
    *           (jcMaxSolve - jcMinSolve + 1)
    *           (kcMaxSolve - kcMinSolve + 1);

  n3SolveNode = 3
    *       (inmaxsolve - inminsolve + 1)
    *       (jnmaxsolve - jnminsolve + 1)
    *       (knmaxsolve - knminsolve + 1);

    // loop over species and cell centers: fill in charge density
    for (int is=0; is < ns; is++)
      for (int i=0; i < nxn; i++)
	for (int j=0; j < nyn; j++)
	  for (int k=0; k < nzn; k++)
	    {
	      // WARNING getFluidRho contains "case" statment
	      rhons[is][i][j][k] = col->getPICRhoNum(i,j,k,is);
	    }

    // interpolate from cell centers to nodes (corners of cells)
    for (int is=0 ; is<ns; is++)
      grid->interpN2Cfull(rhocs[is],rhons[is]);
      
    // copy over Electric and Magnetic fields
    SyncWithFluid(col,grid,vct);

    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for (int k=0; k < nzc; k++)
	  {
	    Bxc[i][j][k] = fluidBxc[i][j][k];
	    Byc[i][j][k] = fluidByc[i][j][k];
	    Bzc[i][j][k] = fluidBzc[i][j][k];
	  }
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++)
	  {
	    Bxn[i][j][k] = fluidBxn[i][j][k];
	    Byn[i][j][k] = fluidByn[i][j][k];
	    Bzn[i][j][k] = fluidBzn[i][j][k];
	    Ex[i][j][k]  = fluidEx[i][j][k];
	    Ey[i][j][k]  = fluidEy[i][j][k];
	    Ez[i][j][k]  = fluidEz[i][j][k];

	    Bx_ext[i][j][k] = 0.0;
	    By_ext[i][j][k] = 0.0;
	    Bz_ext[i][j][k] = 0.0;
	  }
    
}


 void EMfields3D::calculateFluidPressure(){
   const VirtualTopology3D *vct = &get_vct();
   string nameSub = "calculateFluidPressure";   
   double ux, uy, uz, rhoMass, pRatio;
   pRatio = 0.01;
   for (register int is = 0; is < ns; is++)
     for (register int i = 0; i < nxn; i++)
       for (register int j = 0; j < nyn; j++)
	 for (register int k = 0; k < nzn; k++) {
	   rhoMass = rhons[is][i][j][k]/qom[is];
	   if(rhoMass>1e-99){
	     ux = Jxs[is][i][j][k]/rhons[is][i][j][k];
	     uy = Jys[is][i][j][k]/rhons[is][i][j][k];
	     uz = Jzs[is][i][j][k]/rhons[is][i][j][k];

	     p0XXsn[is][i][j][k] = pXXsn[is][i][j][k]/qom[is] - rhoMass*ux*ux;
	     p0YYsn[is][i][j][k] = pYYsn[is][i][j][k]/qom[is] - rhoMass*uy*uy;
	     p0ZZsn[is][i][j][k] = pZZsn[is][i][j][k]/qom[is] - rhoMass*uz*uz;

	     p0XYsn[is][i][j][k] = pXYsn[is][i][j][k]/qom[is] - rhoMass*ux*uy;
	     p0XZsn[is][i][j][k] = pXZsn[is][i][j][k]/qom[is] - rhoMass*ux*uz;
	     p0YZsn[is][i][j][k] = pYZsn[is][i][j][k]/qom[is] - rhoMass*uy*uz;
	   }else{
	     p0XXsn[is][i][j][k] = 0;
	     p0YYsn[is][i][j][k] = 0;
	     p0ZZsn[is][i][j][k] = 0;

	     p0XYsn[is][i][j][k] = 0;
	     p0XZsn[is][i][j][k] = 0;
	     p0YZsn[is][i][j][k] = 0;
	   }
	 } // for 
 }

void EMfields3D::setNodePlasmaBC(int is)
{
  const VirtualTopology3D *vct = &get_vct();
  const Collective *col = &get_col();
  const double signq = qom[is]/(fabs(qom[is]));

  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL ) {
    for (int i = 1; i < nyn - 1; i++)
      for (int k = 1; k < nzn - 1; k++)
	{
	  rhons[is][1][i][k] = signq*col->getPICRhoNum(1,i,k,is);
	  Jxs  [is][1][i][k] = col->getPICJx(1,i,k,is);
	  Jys  [is][1][i][k] = col->getPICJy(1,i,k,is);
	  Jzs  [is][1][i][k] = col->getPICJz(1,i,k,is);
	  pXXsn[is][1][i][k] = col->getPICPxx(1,i,k,is);
	  pXYsn[is][1][i][k] = col->getPICPxy(1,i,k,is);
	  pXZsn[is][1][i][k] = col->getPICPxz(1,i,k,is);
	  pYYsn[is][1][i][k] = col->getPICPyy(1,i,k,is);
	  pYZsn[is][1][i][k] = col->getPICPyz(1,i,k,is);
	  pZZsn[is][1][i][k] = col->getPICPzz(1,i,k,is);
	}
  }
  
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int k = 1; k < nzn - 1; k++)
	{
	  rhons[is][i][1][k] = signq*col->getPICRhoNum(i,1,k,is);
	  Jxs  [is][i][1][k] = col->getPICJx(i,1,k,is);
	  Jys  [is][i][1][k] = col->getPICJy(i,1,k,is);
	  Jzs  [is][i][1][k] = col->getPICJz(i,1,k,is);
	  pXXsn[is][i][1][k] = col->getPICPxx(i,1,k,is);
	  pXYsn[is][i][1][k] = col->getPICPxy(i,1,k,is);
	  pXZsn[is][i][1][k] = col->getPICPxz(i,1,k,is);
	  pYYsn[is][i][1][k] = col->getPICPyy(i,1,k,is);
	  pYZsn[is][i][1][k] = col->getPICPyz(i,1,k,is);
	  pZZsn[is][i][1][k] = col->getPICPzz(i,1,k,is);
	}
  }
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	{
	  rhons[is][i][j][1] = signq*col->getPICRhoNum(i,j,1,is);
	  Jxs  [is][i][j][1] = col->getPICJx(i,j,1,is);
	  Jys  [is][i][j][1] = col->getPICJy(i,j,1,is);
	  Jzs  [is][i][j][1] = col->getPICJz(i,j,1,is);
	  pXXsn[is][i][j][1] = col->getPICPxx(i,j,1,is);
	  pXYsn[is][i][j][1] = col->getPICPxy(i,j,1,is);
	  pXZsn[is][i][j][1] = col->getPICPxz(i,j,1,is);
	  pYYsn[is][i][j][1] = col->getPICPyy(i,j,1,is);
	  pYZsn[is][i][j][1] = col->getPICPyz(i,j,1,is);
	  pZZsn[is][i][j][1] = col->getPICPzz(i,j,1,is);
	}
  }
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nyn - 1; i++)
      for (int k = 1; k < nzn - 1; k++)
	{
	  rhons[is][nxn - 2][i][k] = signq*col->getPICRhoNum(nxn-2,i,k,is);
	  Jxs  [is][nxn - 2][i][k] = col->getPICJx(nxn-2,i,k,is);
	  Jys  [is][nxn - 2][i][k] = col->getPICJy(nxn-2,i,k,is);
	  Jzs  [is][nxn - 2][i][k] = col->getPICJz(nxn-2,i,k,is);
	  pXXsn[is][nxn - 2][i][k] = col->getPICPxx(nxn-2,i,k,is);
	  pXYsn[is][nxn - 2][i][k] = col->getPICPxy(nxn-2,i,k,is);
	  pXZsn[is][nxn - 2][i][k] = col->getPICPxz(nxn-2,i,k,is);
	  pYYsn[is][nxn - 2][i][k] = col->getPICPyy(nxn-2,i,k,is);
	  pYZsn[is][nxn - 2][i][k] = col->getPICPyz(nxn-2,i,k,is);
	  pZZsn[is][nxn - 2][i][k] = col->getPICPzz(nxn-2,i,k,is);
	}
  }
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int k = 1; k < nzn - 1; k++)
	{
	  rhons[is][i][nyn - 2][k] = signq*col->getPICRhoNum(i,nyn-2,k,is);
	  Jxs  [is][i][nyn - 2][k] = col->getPICJx(i,nyn-2,k,is);
	  Jys  [is][i][nyn - 2][k] = col->getPICJy(i,nyn-2,k,is);
	  Jzs  [is][i][nyn - 2][k] = col->getPICJz(i,nyn-2,k,is);
	  pXXsn[is][i][nyn - 2][k] = col->getPICPxx(i,nyn-2,k,is);
	  pXYsn[is][i][nyn - 2][k] = col->getPICPxy(i,nyn-2,k,is);
	  pXZsn[is][i][nyn - 2][k] = col->getPICPxz(i,nyn-2,k,is);
	  pYYsn[is][i][nyn - 2][k] = col->getPICPyy(i,nyn-2,k,is);
	  pYZsn[is][i][nyn - 2][k] = col->getPICPyz(i,nyn-2,k,is);
	  pZZsn[is][i][nyn - 2][k] = col->getPICPzz(i,nyn-2,k,is);
	}
  }
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
	{
	  rhons[is][i][j][nzn - 2] = signq*col->getPICRhoNum(i,j,nzn-2,is);
	  Jxs  [is][i][j][nzn - 2] = col->getPICJx(i,j,nzn-2,is);
	  Jys  [is][i][j][nzn - 2] = col->getPICJy(i,j,nzn-2,is);
	  Jzs  [is][i][j][nzn - 2] = col->getPICJz(i,j,nzn-2,is);
	  pXXsn[is][i][j][nzn - 2] = col->getPICPxx(i,j,nzn-2,is);
	  pXYsn[is][i][j][nzn - 2] = col->getPICPxy(i,j,nzn-2,is);
	  pXZsn[is][i][j][nzn - 2] = col->getPICPxz(i,j,nzn-2,is);
	  pYYsn[is][i][j][nzn - 2] = col->getPICPyy(i,j,nzn-2,is);
	  pYZsn[is][i][j][nzn - 2] = col->getPICPyz(i,j,nzn-2,is);
	  pZZsn[is][i][j][nzn - 2] = col->getPICPzz(i,j,nzn-2,is);
	}
  }
  
}


/*! Copy field form the fluid solution for boundary condition*/
void EMfields3D::SyncWithFluid(CollectiveIO *col,Grid *grid,VirtualTopology3D *vct)
{
  int minDn;
  string funcName="EMfields3D::SyncWithFluid";

  // bool doTestPart = false;

  // if(doTestPart){
  //   double ip0 = 10, ip1=20;

  //   double q0 = 8e-3;

  //   double E0 = 0.5*q0/(dx*dy*dz)*dx*FourPI;
    
  //   for (int i=0; i < nxn; i++)
  //     for (int j=0; j < nyn; j++)
  // 	for (int k=0; k < nzn; k++){
  // 	  fluidEx[i][j][k] = 0;
  // 	  if(i < ip0 || i > ip1) fluidEx[i][j][k] = E0;
  // 	  if(i > ip0 && i < ip1) fluidEx[i][j][k] = -E0;

	  
  // 	  fluidEy[i][j][k] = 0;
  // 	  fluidEz[i][j][k] = 0;
  // 	  fluidBxn[i][j][k] = 0;
  // 	  fluidByn[i][j][k] = 0;
  // 	  fluidBzn[i][j][k] = 0;
  // 	}
    
  // }

  bool useUserIC;
  useUserIC = false;

  if(useUserIC){
    // Users can define EM initial conditions here, and then 
    // the initial conditions from MHD will be ignored. 
    double Ex0=0, Ey0=0, Ez0=0; 
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++){
	  double xn, yn, zn, kdotx;
	  if(i<nxn/3){
	    Ex0 = 0;
	  }else if(i<2*nxn/3){
	    Ex0 = 1e-1; 
	  }else{
	    Ex0 = 0; 
	  }
	  fluidEx[i][j][k] = Ex0;
	  fluidEy[i][j][k] = Ey0;
	  fluidEz[i][j][k] = Ez0;
	  fluidBxn[i][j][k] = 0;
	  fluidByn[i][j][k] = 0;
	  fluidBzn[i][j][k] = 0;
	}

  }else if(col->getdoTestEMWave()){
    double k_D[3], E0_D[3], phi0;
    int x_=0, y_=1, z_=2;
    
    for(int i=0; i<3; i++){
      k_D[i] = col->getwaveVec_D(i);
      E0_D[i] = col->getamplE_D(i);
    }

    phi0 = col->getphase0deg()*FourPI/360.0;
    
    double k0 = sqrt(k_D[x_]*k_D[x_] + k_D[y_]*k_D[y_] + k_D[z_]*k_D[z_]);
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++){
	  double xn, yn, zn, kdotx;
	  xn = grid->getXN(i,j,k);
	  yn = grid->getYN(i,j,k);
	  zn = grid->getZN(i,j,k);
	  kdotx = xn*k_D[x_] + yn*k_D[y_] + zn*k_D[z_];

	  fluidEx[i][j][k] = E0_D[x_]*sin(kdotx + phi0);
	  fluidEy[i][j][k] = E0_D[y_]*sin(kdotx + phi0);
	  fluidEz[i][j][k] = E0_D[z_]*sin(kdotx + phi0);
	  fluidBxn[i][j][k] = (k_D[y_]*fluidEz[i][j][k] - k_D[z_]*fluidEy[i][j][k])/k0;// + E0_D[x_];
	  fluidByn[i][j][k] = (k_D[z_]*fluidEx[i][j][k] - k_D[x_]*fluidEz[i][j][k])/k0;// + E0_D[y_];
	  fluidBzn[i][j][k] = (k_D[x_]*fluidEy[i][j][k] - k_D[y_]*fluidEx[i][j][k])/k0;// + E0_D[z_];	  	 
	}
    
  }else{
  
    // loop over cell centers and fill in magnetic and electric fields
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++)
	  {
	    col->setFluidFieldsNode(&fluidEx[i][j][k],&fluidEy[i][j][k],&fluidEz[i][j][k],
				    &fluidBxn[i][j][k],&fluidByn[i][j][k],&fluidBzn[i][j][k],i,j,k);
	  }
  }
  // interpolate from cell centers to nodes (corners of cells)
  grid->interpN2Cfull(fluidExc,fluidEx);
  grid->interpN2Cfull(fluidEyc,fluidEy);
  grid->interpN2Cfull(fluidEzc,fluidEz);
  grid->interpN2Cfull(fluidBxc,fluidBxn);
  grid->interpN2Cfull(fluidByc,fluidByn);
  grid->interpN2Cfull(fluidBzc,fluidBzn);

  // Get new timestep from fluid solver
  dt  = col->getFluidDt();
  delt = c * th * dt;

  for(int is=0; is<ns;is++) setNodePlasmaBC(is);
}

/** fix the E boundary when running with batsrus*/
void EMfields3D::fixE_BATSRUS(arr3_double VecX, arr3_double VecY, arr3_double VecZ, bool doSolveForChange){
  const VirtualTopology3D *vct = &get_vct();

  double coeffBcE;

  if(doSolveForChange)  coeffBcE=0.0;
  else                  coeffBcE=1.0;

  if (vct->getZright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    for (int i=0; i < nxn;i++)
      for (int j=0; j < nyn;j++)
	for (int k=0; k <= NG_F;k++){
	  VecX[i][j][nzn-1-k] = fluidEx[i][j][nzn-1-k]*coeffBcE;
	  VecY[i][j][nzn-1-k] = fluidEy[i][j][nzn-1-k]*coeffBcE;
	  VecZ[i][j][nzn-1-k] = fluidEz[i][j][nzn-1-k]*coeffBcE;
	}
  }
  if (vct->getZleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    for (int i=0; i < nxn;i++)
      for (int j=0; j< nyn;j++)
	for (int k=0; k <= NG_F;k++){
	  VecX[i][j][k] = fluidEx[i][j][k]*coeffBcE;
	  VecY[i][j][k] = fluidEy[i][j][k]*coeffBcE;
	  VecZ[i][j][k] = fluidEz[i][j][k]*coeffBcE;
	}
  }
  if (vct->getYright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    for (int i=0; i < nxn;i++)
      for (int j=0;j<=NG_F;j++)
	for (int k=0; k < nzn;k++){
	  VecX[i][nyn-1-j][k] = fluidEx[i][nyn-1-j][k]*coeffBcE;
	  VecY[i][nyn-1-j][k] = fluidEy[i][nyn-1-j][k]*coeffBcE;
	  VecZ[i][nyn-1-j][k] = fluidEz[i][nyn-1-j][k]*coeffBcE;
	}
  }
  if (vct->getYleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    for (int i=0; i < nxn;i++)
      for (int j=0;j<=NG_F;j++)
	for (int k=0; k < nzn;k++){
	  VecX[i][j][k] = fluidEx[i][j][k]*coeffBcE;
	  VecY[i][j][k] = fluidEy[i][j][k]*coeffBcE;
	  VecZ[i][j][k] = fluidEz[i][j][k]*coeffBcE;
	}
  }
  if (vct->getXright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    for (int i=0;i<=NG_F;i++)
      for (int j=0; j < nyn;j++)
	for (int k=0; k < nzn;k++){
	  VecX[nxn-1-i][j][k] = fluidEx[nxn-1-i][j][k]*coeffBcE;
	  VecY[nxn-1-i][j][k] = fluidEy[nxn-1-i][j][k]*coeffBcE;
	  VecZ[nxn-1-i][j][k] = fluidEz[nxn-1-i][j][k]*coeffBcE;
	}
  }
  if (vct->getXleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    for (int i=0;i<=NG_F;i++)
      for (int j=0; j < nyn;j++)
	for (int k=0; k < nzn;k++){
	  VecX[i][j][k] = fluidEx[i][j][k]*coeffBcE;
	  VecY[i][j][k] = fluidEy[i][j][k]*coeffBcE;
	  VecZ[i][j][k] = fluidEz[i][j][k]*coeffBcE;
	}
  }
}





/** fix the B boundary when running with BATSRUS*/
void EMfields3D::fixB_BATSRUS(){
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getZright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    // CELL
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for(int ng=1;ng<=NG_F+1; ng++){
	  Bxc[i][j][nzc-ng] = fluidBxc[i][j][nzc-ng];
	  Byc[i][j][nzc-ng] = fluidByc[i][j][nzc-ng];
	  Bzc[i][j][nzc-ng] = fluidBzc[i][j][nzc-ng];
	}
    // NODE
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for(int ng=1;ng<=NG_F+1; ng++){
	  Bxn[i][j][nzn-ng] = fluidBxn[i][j][nzn-ng];
	  Byn[i][j][nzn-ng] = fluidByn[i][j][nzn-ng];
	  Bzn[i][j][nzn-ng] = fluidBzn[i][j][nzn-ng];
	}
  }

  if (vct->getZleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    // CELL
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for(int ng=0;ng<NG_F+1; ng++){
	  Bxc[i][j][ng] = fluidBxc[i][j][ng];
	  Byc[i][j][ng] = fluidByc[i][j][ng];
	  Bzc[i][j][ng] = fluidBzc[i][j][ng];
	}
    // NODE
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for(int ng=0;ng<NG_F+1; ng++){
	  Bxn[i][j][ng] = fluidBxn[i][j][ng];
	  Byn[i][j][ng] = fluidByn[i][j][ng];
	  Bzn[i][j][ng] = fluidBzn[i][j][ng];
	}
  }


  if (vct->getYright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    // CELL
    for (int i=0; i < nxc;i++)
      for (int k=0; k < nzc;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  Bxc[i][nyc-ng][k] = fluidBxc[i][nyc-ng][k];
	  Byc[i][nyc-ng][k] = fluidByc[i][nyc-ng][k];
	  Bzc[i][nyc-ng][k] = fluidBzc[i][nyc-ng][k];
	}
    // NODE
    for (int i=0; i < nxn;i++)
      for (int k=0; k < nzn;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  Bxn[i][nyn-ng][k] = fluidBxn[i][nyn-ng][k];
	  Byn[i][nyn-ng][k] = fluidByn[i][nyn-ng][k];
	  Bzn[i][nyn-ng][k] = fluidBzn[i][nyn-ng][k];
	}
  }

  if (vct->getYleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    // CELL
    for (int i=0; i < nxc;i++)
      for (int k=0; k < nzc;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  Bxc[i][ng][k] = fluidBxc[i][ng][k];
	  Byc[i][ng][k] = fluidByc[i][ng][k];
	  Bzc[i][ng][k] = fluidBzc[i][ng][k];
	}
    // NODE
    for (int i=0; i < nxn;i++)
      for (int k=0; k < nzn;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  Bxn[i][ng][k] = fluidBxn[i][ng][k];
	  Byn[i][ng][k] = fluidByn[i][ng][k];
	  Bzn[i][ng][k] = fluidBzn[i][ng][k];
	}
  }

  if (vct->getXright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    // CELL
    for (int j=0; j < nyc;j++)
      for (int k=0; k < nzc;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  Bxc[nxc-ng][j][k] = fluidBxc[nxc-ng][j][k];
	  Byc[nxc-ng][j][k] = fluidByc[nxc-ng][j][k];
	  Bzc[nxc-ng][j][k] = fluidBzc[nxc-ng][j][k];
	}
    // NODE
    for (int j=0; j < nyn;j++)
      for (int k=0; k < nzn;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  Bxn[nxn-ng][j][k] = fluidBxn[nxn-ng][j][k];
	  Byn[nxn-ng][j][k] = fluidByn[nxn-ng][j][k];
	  Bzn[nxn-ng][j][k] = fluidBzn[nxn-ng][j][k];
	}
  }

  if (vct->getXleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    // CELL
    for (int j=0; j < nyc;j++)
      for (int k=0; k < nzc;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  Bxc[ng][j][k] = fluidBxc[ng][j][k];
	  Byc[ng][j][k] = fluidByc[ng][j][k];
	  Bzc[ng][j][k] = fluidBzc[ng][j][k];
	}
    // NODE
    for (int j=0; j < nyn;j++)
      for (int k=0; k < nzn;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  Bxn[ng][j][k] = fluidBxn[ng][j][k];
	  Byn[ng][j][k] = fluidByn[ng][j][k];
	  Bzn[ng][j][k] = fluidBzn[ng][j][k];
	}
  }
}




/** fix the PHI boundary when running with BATSRUS. 

           |-----|-------|-------|---------|------
global ic:    0     1        2        3        4 ..... nxc-2   nxc-1
where nxc include ghost cells. 

The index for PHI is from 0 to nxc-1, but for MHD-EPIC, cell 0 and nxc-1 is 
meaningless and the real ghost cells are 1 and nxc-2. The Poisson correction
takes cells 1 and nxc-2 as boundary condition, and the floating BC is set here. 
The cells from 2 and nxc-3 are solved by implicit solver. See the definition 
of icMinSolve ... kcMaxSolve. 

For convenience, the PHI value at cell 0 and nxc-1 is also proved based on 
floating BC. 
*/
void EMfields3D::fixPHI_BATSRUS(){
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (vct->getZright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    // CELL
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for(int ng=1;ng<=NG_F+1; ng++){
	  PHI[i][j][nzc-ng] = PHI[i][j][nzc-NG_F-2];
	}
  }

  if (vct->getZleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICZ()){
    // CELL
    for (int i=0; i < nxc; i++)
      for (int j=0; j < nyc; j++)
	for(int ng=0;ng<NG_F+1; ng++){
	  PHI[i][j][ng] = PHI[i][j][NG_F+1];
	}
  }


  if (vct->getYright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    // CELL
    for (int i=0; i < nxc;i++)
      for (int k=0; k < nzc;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  PHI[i][nyc-ng][k] = PHI[i][nyc-NG_F-2][k];
	}
  }

  if (vct->getYleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICY()){
    // CELL
    for (int i=0; i < nxc;i++)
      for (int k=0; k < nzc;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  PHI[i][ng][k] = PHI[i][NG_F+1][k];
	}
  }

  if (vct->getXright_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    // CELL
    for (int j=0; j < nyc;j++)
      for (int k=0; k < nzc;k++)
	for(int ng=1;ng<=NG_F+1;ng++){
	  PHI[nxc-ng][j][k] = PHI[nxc-NG_F-2][j][k];
	}
  }

  if (vct->getXleft_neighbor()==MPI_PROC_NULL && !vct->getPERIODICX()){
    // CELL
    for (int j=0; j < nyc;j++)
      for (int k=0; k < nzc;k++)
	for(int ng=0;ng<NG_F+1;ng++){
	  PHI[ng][j][k] = PHI[NG_F+1][j][k];
	}
  }
}

inline void EMfields3D::fixVarBCcell(arr3_double Var, 
				     double (CollectiveIO::*fluidVar)(const int, const int,const int, const int)const, int nOverlap, int is){

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  
  // total number of overlap+ghost cells
  int nBC;

  // Length of overlap+ghost region
  double LXbc, LYbc, LZbc;

  // distangs to left/right boundary
  double LeftX,RightX,LeftY,RightY,LeftZ,RightZ;

  // min distangs to boundary
  double dist;

  // grid positions
  double x, y, z;
  
  // there is nothing to do for a fully periodic system
  if(vct->getPERIODICX() && vct->getPERIODICY() && vct->getPERIODICZ()) return ;

  if(nOverlap < 0) 
    return; // do nothing
  else if(nOverlap == 0)
    nBC = NG_P; // only fix ghost cells
  else
    nBC = nOverlap +NG_P;
  
  LXbc = nBC*dx;
  LYbc = nBC*dy;
  LZbc = nBC*dz;

  // If Periodic
  LeftX = 2.0;
  RightX = LeftX;
  LeftY = 2.0;
  RightY = LeftY;
  LeftZ = 2.0;
  RightZ = LeftZ;

  for(int i=0;i<nxc;i++){
    x = grid->getXC(i,1,1);
    for (int j=0; j < nyc;j++){
      y = grid->getYC(1,j,1);
      for (int k=0; k < nzc;k++){
	z = grid->getZC(1,1,k);

	if(!vct->getPERIODICX()){
	  LeftX  = x/LXbc;
	  RightX = (Lx -x)/LXbc; 
	}
	
	if(!vct->getPERIODICY()){	
	  LeftY  = y/LYbc;
	  RightY = (Ly-y)/LYbc;
	}

	if(!vct->getPERIODICZ()){	
	  LeftZ  = z/LZbc;
	  RightZ = (Lz-z)/LZbc;
	}

	dist = minval6(LeftX,RightX,LeftY,RightY,LeftZ,RightZ);

	if(dist < 0.0) // we are in ghost cells
	  Var[i][j][k] =  (col->*fluidVar)(i,j,k,is);
	else if (dist < 1.0) // this is the overlap region
	  Var[i][j][k] = Var[i][j][k]*dist + (col->*fluidVar)(i,j,k,is)*(1.0-dist);

      }}}
}

/** fix the charge separation on the  boundary when running with batsrus*/
inline void EMfields3D::fixVarBCnode(arr4_double Var, 
				     double (CollectiveIO::*fluidVar)(int, int, int, int)const, int nOverlap, int is){

  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  
  // total number of overlap+ghost cells
  int nBC;

  // Length of overlap+ghost region
  double LXbc, LYbc, LZbc;

  // distangs to left/right boundary
  double LeftX,RightX,LeftY,RightY,LeftZ,RightZ;

  // min distangs to boundary
  double dist;

  // grid positions
  double x, y, z;
		
  // there is nothing to do for a fully periodic system
  if(vct->getPERIODICX() && vct->getPERIODICY() && vct->getPERIODICZ()) return ;

  if(nOverlap < 0) 
    return; // do nothing
  else if(nOverlap == 0)
    nBC = NG_P; // only fix ghost nodes
  else
    nBC = nOverlap +NG_P;

  LXbc = nBC*dx;
  LYbc = nBC*dy;
  LZbc = nBC*dz;

  // If Periodic
  LeftX = 2.0;
  RightX = LeftX;
  LeftY = 2.0;
  RightY = LeftY;
  LeftZ = 2.0;
  RightZ = LeftZ;

  for(int i=0;i<nxn;i++){
    x = grid->getXN(i,1,1);
    for (int j=0; j < nyn;j++){
      y = grid->getYN(1,j,1);
      for (int k=0; k < nzn;k++){
	z = grid->getZN(1,1,k);

	if(!vct->getPERIODICX()){
	  LeftX  = x/LXbc;
	  RightX = (Lx -x)/LXbc; 
	}

	if(!vct->getPERIODICY()){	
	  LeftY  = y/LYbc;
	  RightY = (Ly-y)/LYbc;
	}

	if(!vct->getPERIODICZ()){	
	  LeftZ  = z/LZbc;
	  RightZ = (Lz-z)/LZbc;
	}

	dist = minval6(LeftX,RightX,LeftY,RightY,LeftZ,RightZ);

	if(dist < 0.0) // we are in ghost cells
	  Var[is][i][j][k] =  (col->*fluidVar)(i,j,k,is);
	else if (dist < 1.0 && nOverlap>0) // this is the overlap region
	  Var[is][i][j][k] = Var[is][i][j][k]*dist + (col->*fluidVar)(i,j,k,is)*(1.0-dist);

      }}}
}

double EMfields3D::calDtMax(double dx, double dy, double dz, double dt, int &iError){
  /* Calculate max dt that satisfies the accuracy condition: max(uth*dt/dx) < 1.  */
  const Collective *col = &get_col();
  double uthLocal, p0, rho0;
  double dtMax=1e10;
  double uthLimit = col->get_maxUth();
  iError = 1; 
  for(int is=0; is<ns; is++){
    uthLocal = 0;
    for (int i=1; i<nxn-1; i++)
      for (int j=1; j<nyn-1; j++)
	for(int k=1; k<nzn-1; k++){
	  p0 = (pXXsn[is][i][j][k] + pYYsn[is][i][j][k] + pZZsn[is][i][j][k])/3;
	  rho0 = rhons[is][i][j][k];
	  if(rho0 !=0)
	    uthLocal = max(uthLocal, sqrt(p0/rho0));
	}

    double uth;          
    MPI_Allreduce(&uthLocal, &uth, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_MYSIM);
    
    uth = max(uth, 1e-99);
    dtMax = min(dtMax, min(dx/uth, min(dy/uth, dz/uth)));
    
    if(get_vct().getCartesian_rank() == 0){
      cout<<"is= "<<is<<" uth= "<<uth<<" dx/uth= "<<dx/uth
	  <<" dy/uth= "<<dy/uth<<" dz/uth= "<<dz/uth<<endl;      
    }

    if(uth > uthLimit) iError = -1;
    
  }// is

  return dtMax;
}

void EMfields3D:: write_plot_field(string filename, string *var_I, int nVar,
				   double** pointList_ID, long nPoint, bool isCoord,
				   const double No2OutL, const double No2OutV,
				   const double No2OutB, const double No2OutRho,
				   const double No2OutP, const double No2OutJ){

  string nameSub = "write_plot_field";
  bool doTestFunc;
  doTestFunc = do_test_func(nameSub);
  //---------------------
  if(doTestFunc){
    cout<<nameSub<<" start"<<endl;
  }
    
  const Collective *col = &get_col();
  const Grid *grid = &get_grid();
  const bool isBinary = col->getdoSaveBinary();

  string * varOut_I;
  varOut_I = new string[nVar+3]; // x/y/z + var 

  varOut_I[0] = 'X';
  varOut_I[1] = 'Y';
  varOut_I[2] = 'Z';
  for (int iVar=0; iVar<nVar; iVar++) varOut_I[iVar+3] = var_I[iVar];

  if(isBinary){
    FILE *outFile;
    int nRecord, nSizeDouble, nSizeInt;
    nSizeInt = sizeof(int);
    assert_eq(nSizeInt,4);
    nSizeDouble = sizeof(double);
    nRecord = (nVar+4)*nSizeDouble;
    if(doTestFunc){
      cout<<"size of int is "<<sizeof(int)
	  <<"size of double is "<<sizeof(double)
	  <<endl;
    }
    outFile = fopen(filename.c_str(),"wb");
    
    double data0;

    for(long iPoint = 0; iPoint<nPoint; iPoint++){
      // The PostIDL.f90 was originally designed for Fortran output. In order to
      // use PostIDL.f90, we should follow the format of Fortran
      // binary output. Each line is a record. Before and after
      // each record, use 4 byte to save the length of this record. 
      fwrite(&nRecord, nSizeInt, 1, outFile);
      data0 = dx*No2OutL;
      fwrite(&data0, nSizeDouble, 1, outFile);
      for(int iVar=0; iVar<nVar+3; iVar++){
	data0 = getVar(varOut_I[iVar],
		       pointList_ID[iPoint][0],
		       pointList_ID[iPoint][1],
		       pointList_ID[iPoint][2],
		       isCoord,
		       No2OutL, No2OutV, No2OutB,
		       No2OutRho, No2OutP, No2OutJ);
	fwrite(&data0, nSizeDouble, 1, outFile);
      }
      fwrite(&nRecord, nSizeInt, 1, outFile);
    }
    fclose(outFile);
  }else{
    ofstream outFile;
    outFile.open(filename.c_str(),fstream::out | fstream::trunc);
    outFile<<std::scientific;
    outFile.precision(7);
    for(long iPoint = 0; iPoint<nPoint; iPoint++){
      outFile<<dx*No2OutL;
      for(int iVar=0; iVar<nVar+3; iVar++){
	outFile<<"\t"<<getVar(varOut_I[iVar],
			      pointList_ID[iPoint][0],
			      pointList_ID[iPoint][1],
			      pointList_ID[iPoint][2],
			      isCoord,
			      No2OutL, No2OutV,
			      No2OutB, No2OutRho, No2OutP, No2OutJ);
      }	
      outFile<<"\n";
    }
    if(outFile.is_open()) outFile.close(); 
  }
  
  delete [] varOut_I;
}

double EMfields3D:: getVar(string var, double iIn, double jIn, double kIn, bool isCoord, 
			   const double No2OutL, const double No2OutV,
			   const double No2OutB, const double No2OutRho,
			   const double No2OutP, const double No2OutJ){
  /*
    if isCoord is true: iIn,jIn,kIn are the location (x/y/z) of the output point.
    else: iIn,jIn,kIn are the node index, eventhough they are double number.
   */
  
  const Collective *col = &get_col();
  const Grid *grid = &get_grid();
  int i, j, k;
  double value;

  if(isCoord){
    int ix, iy, iz;
    double weight_I[8];

    // The interpolate weights should be stored in an array so that it
    // only need to be calculated once. Yuxi 
    grid->getInterpolateWeight(iIn,jIn,kIn,ix,iy,iz,weight_I);

    const double w000 = weight_I[0];
    const double w001 = weight_I[1];
    const double w010 = weight_I[2];
    const double w011 = weight_I[3];
    const double w100 = weight_I[4];
    const double w101 = weight_I[5];
    const double w110 = weight_I[6];
    const double w111 = weight_I[7];

    value = 0;
    value += w000*getVar(var,ix,iy,iz,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w001*getVar(var,ix,iy,iz-1,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w010*getVar(var,ix,iy-1,iz,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w011*getVar(var,ix,iy-1,iz-1,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w100*getVar(var,ix-1,iy,iz,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w101*getVar(var,ix-1,iy,iz-1,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w110*getVar(var,ix-1,iy-1,iz,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);
    value += w111*getVar(var,ix-1,iy-1,iz-1,false,No2OutL, No2OutV,
			 No2OutB, No2OutRho, No2OutP, No2OutJ);

    return value;
  }else{
    i = iIn; j=jIn; k=kIn;
  }
  
  if(var.substr(0,2)=="qc"){
    // "qc", "qcS0", "qdS1"...
    if(var.size()==2){
      value = rhoc[i][j][k];
    }else{
      // "qcS0"...
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      value = rhocs[is][i][j][k];
    }
    value *= No2OutV*No2OutB/No2OutL;
  }else if(var.substr(0,1)=="q"){
    // "q", "qS0", "qS1"...
    if(var.size()==1){
      value = 0;
      for(int is=0; is<ns; is++){
	value += rhons[is][i][j][k];
      }      
    }else{
      // "qS0"...
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      value = rhons[is][i][j][k];
    }
    value *= No2OutV*No2OutB/No2OutL;
  }else if(var.substr(0,4)=="rhoc"){
    // "rho", "rhoS0", "rhoS1"...
    if(var.size()==4){
      value = 0;
      for(int is=0; is<ns; is++){
	value += rhocs[is][i][j][k]/qom[is];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = rhocs[is][i][j][k]/qom[is];
      }
    }
    value *= No2OutRho;
  }else if(var.substr(0,3)=="rho"){
    // "rho", "rhoS0", "rhoS1"...
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += rhons[is][i][j][k]/qom[is];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = rhons[is][i][j][k]/qom[is];
      }
    }
    value *= No2OutRho;
  }else if(var.substr(0,2)=="pS"){
    // pS0, pS1...
    bool isFluidP = true;
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    if(is>=ns){
      value=0;
    }else{
      value = (calPxx(is,i,j,k,isFluidP) +
	       calPyy(is,i,j,k,isFluidP) +
	       calPzz(is,i,j,k,isFluidP)) /3.0;
      value *= No2OutP;
    }
  }else if(var.substr(0,3)=="pXX" || var.substr(0,3) == "kXX"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPxx(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = calPxx(is,i,j,k,isFluidP);
      }
    }
    value *= No2OutP;
  }else if(var.substr(0,3)=="pYY" || var.substr(0,3) == "kYY"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPyy(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{	
	value = calPyy(is,i,j,k,isFluidP);
      }
    }
    value *= No2OutP; 
  }else if(var.substr(0,3)=="pZZ" || var.substr(0,3) == "kZZ"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPzz(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = calPzz(is,i,j,k,isFluidP);
      }
    }
    value *= No2OutP;
  }else if(var.substr(0,3)=="pXY" || var.substr(0,3) == "kXY"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPxy(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      value = calPxy(is,i,j,k,isFluidP);
    }
    value *= No2OutP;
  }else if(var.substr(0,3)=="pXZ" || var.substr(0,3) == "kXZ"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPxz(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = calPxz(is,i,j,k,isFluidP);
      }
    }
    value *= No2OutP;
  }else if(var.substr(0,3)=="pYZ" || var.substr(0,3) == "kYZ"){
    bool isFluidP;
    isFluidP = var.substr(0,1)=="p";
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += calPyz(is,i,j,k,isFluidP);
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = calPyz(is,i,j,k,isFluidP);
      }
    }
    value *= No2OutP;
  }else if(var.substr(0,3)=="jxh"){
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jxsh[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = Jxsh[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,3)=="jyh"){
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jysh[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = Jysh[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,3)=="jzh"){
    if(var.size()==3){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jzsh[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = Jzsh[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,3)=="uxh"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0; 
    if(is<ns){
      if(rhons[is][i][j][k] != 0)
	value = Jxsh[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,3)=="uyh"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0;
    if(is<ns){ 
      if(rhons[is][i][j][k] != 0)
	value = Jysh[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,3)=="uzh"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0; 
    if(is<ns){ 
      if(rhons[is][i][j][k] != 0)
	value = Jzsh[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,2)=="jx"){
    if(var.size()==2){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jxs[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = Jxs[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,2)=="jy"){
    if(var.size()==2){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jys[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
      value = Jys[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,2)=="jz"){
    if(var.size()==2){
      value = 0;
      for(int is=0; is<ns; is++){
	value += Jzs[is][i][j][k];
      }      
    }else{
      string::size_type pos;
      stringstream ss;
      int is;
      pos = var.find_first_of("0123456789");
      ss<<var.substr(pos);
      ss>>is;
      if(is>=ns){
	value=0;
      }else{
	value = Jzs[is][i][j][k];
      }
    }
    value *= No2OutJ;
  }else if(var.substr(0,2)=="ux"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0; 
    if(is<ns){
      if(rhons[is][i][j][k] != 0)
	value = Jxs[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,2)=="uy"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0;
    if(is<ns){ 
      if(rhons[is][i][j][k] != 0)
	value = Jys[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,2)=="uz"){
    string::size_type pos;
    stringstream ss;
    int is;
    pos = var.find_first_of("0123456789");
    ss<<var.substr(pos);
    ss>>is;
    value = 0; 
    if(is<ns){
      if(rhons[is][i][j][k] != 0)
	value = Jzs[is][i][j][k]/rhons[is][i][j][k];
      value *= No2OutV;
    }
  }else if(var.substr(0,5)=="divEc"){
    value = divEc[i][j][k];
    value *= No2OutV*No2OutB/No2OutL;
  }else if(var.substr(0,3)=="phi"){
    value = PHI[i][j][k];
  }else if(var.substr(0,10)=="energyDiff"){ 
    // The unit is [E]^2
    value = energyError_G[i][j][k];                                          
    value *= No2OutV*No2OutB*No2OutV*No2OutB; 
  }else if(var.substr(0,2)=="Ex"){
    value = Ex[i][j][k];
    value *= No2OutV*No2OutB;
  }else if(var.substr(0,2)=="Ey"){
    value = Ey[i][j][k];
    value *= No2OutV*No2OutB;
  }else if(var.substr(0,2)=="Ez"){
    value = Ez[i][j][k];
    value *= No2OutV*No2OutB;
  }else if(var.substr(0,2)=="Bx"){
    value  = Bxn[i][j][k];
    value *= No2OutB;
  }else if(var.substr(0,2)=="By"){
    value = Byn[i][j][k];
    value *= No2OutB;
  }else if(var.substr(0,2)=="Bz"){
    value = Bzn[i][j][k];
    value *= No2OutB;
  }else if(var.substr(0,6)=="smooth"){
    value = col->getSmoothFactor(i,j,k);
  }else if(var.substr(0,1)=="X"){
    if(col->getdoRotate()){
      // Save coord in PIC coordinates.
      value = (grid->getXN(i))*No2OutL;
    }else{
      // Save in MHD coordinates.
      value = (grid->getXN(i) + col->getFluidStartX())*No2OutL;
    }
  }else if(var.substr(0,1)=="Y"){
    if(col->getdoRotate()){
      value = (grid->getYN(j))*No2OutL;
    }else{
      value = (grid->getYN(j) + col->getFluidStartY())*No2OutL;
    }
  }else if(var.substr(0,1)=="Z"){
    double z0=0;
    if(!(col->getdoRotate())) z0 = col->getFluidStartZ();
    value = col->getnDim()==2?
      0:(grid->getZN(k) + z0)*No2OutL;	  
  }else{
    value = 0;
  }
  return value;
}

double EMfields3D:: calPxx(const int is, const int i, const int j, const int k, const bool isFluidP){
  // pXXsn also includes the energy from bulk motion. This function calculates
  // 'real' pXX.
  if(not isFluidP){
    return pXXsn[is][i][j][k];
  }else{
    return p0XXsn[is][i][j][k];
  }

}

double EMfields3D:: calPyy(const int is, const int i, const int j, const int k, const bool isFluidP){
  if(not isFluidP){
    return pYYsn[is][i][j][k];
  }else{
    return p0YYsn[is][i][j][k];
  }

}

double EMfields3D:: calPzz(const int is, const int i, const int j, const int k, const bool isFluidP){
  if(not isFluidP){
    return pZZsn[is][i][j][k];
  }else{
    return p0ZZsn[is][i][j][k];
  }

}

double EMfields3D:: calPxy(const int is, const int i, const int j, const int k, const bool isFluidP){
  if(not isFluidP){
    return pXYsn[is][i][j][k];
  }else{
    return p0XYsn[is][i][j][k];
  }
}

double EMfields3D:: calPxz(const int is, const int i, const int j, const int k, const bool isFluidP){
  if(not isFluidP){
    return pXZsn[is][i][j][k];
  }else{
    return p0XZsn[is][i][j][k];
  }
}

double EMfields3D:: calPyz(const int is, const int i, const int j, const int k, const bool isFluidP){
  if(not isFluidP){
    return pYZsn[is][i][j][k];
  }else{
    return p0YZsn[is][i][j][k];
  }
}

#endif


/*! initiliaze uniform EM */
void EMfields3D::initUniform()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  if (restart1 == 0) {
    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize Uniform EM field" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          for (int is = 0; is < ns; is++) {
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }	  
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = B0x;
          Byn[i][j][k] = B0y; 
          Bzn[i][j][k] = B0z;
        }
    
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x;
          Byc[i][j][k] = B0y; 
          Bzc[i][j][k] = B0z;
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}

/*! initiliaze EM for GEM challange */
void EMfields3D::initGEM()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  if (restart1 == 0) {
    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Perturbation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
          Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXN(i, j, k) - Lx / 2;
          ypert = grid->getYN(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly );
          Byc[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXC(i, j, k) - Lx / 2;
          ypert = grid->getYC(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxc[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byc[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}


void EMfields3D::initNullPoints()
{
	const VirtualTopology3D *vct = &get_vct();
	const Grid *grid = &get_grid();
	if (restart1 ==0){
		if (vct->getCartesian_rank() ==0){
			cout << "----------------------------------------" << endl;
			cout << "       Initialize 3D null point(s)" << endl;
			cout << "----------------------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i] << endl;
			}
			cout << "Smoothing Factor = " << Smooth << endl;
			cout << "-------------------------" << endl;
		}

        for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
		for (int k=0; k < nzn; k++){
		   // initialize the density for species
		   for (int is=0; is < ns; is++)
			   rhons[is][i][j][k] = rhoINIT[is]/FourPI;

			// electric field
			Ex[i][j][k] =  0.0;
			Ey[i][j][k] =  0.0;
			Ez[i][j][k] =  0.0;
			// Magnetic field
			Bxn[i][j][k] = -B0x*cos(2.*M_PI*grid->getXN(i,j,k)/Lx)*sin(2.*M_PI*grid->getYN(i,j,k)/Ly);
			Byn[i][j][k] = B0x*cos(2.*M_PI*grid->getYN(i,j,k)/Ly)*(-2.*sin(2.*M_PI*grid->getZN(i,j,k)/Lz) + sin(2.*M_PI*grid->getXN(i,j,k)/Lx));
			Bzn[i][j][k] = 2.*B0x*cos(2.*M_PI*grid->getZN(i,j,k)/Lz)*sin(2.*M_PI*grid->getYN(i,j,k)/Ly);
		}

		for (int i=0; i <nxc; i++)
		for (int j=0; j <nyc; j++)
		for (int k=0; k <nzc; k++) {
			Bxc[i][j][k] = -B0x*cos(2.*M_PI*grid->getXC(i,j,k)/Lx)*sin(2.*M_PI*grid->getYC(i,j,k)/Ly);
			Byc[i][j][k] = B0x*cos(2.*M_PI*grid->getYC(i,j,k)/Ly)*(-2.*sin(2.*M_PI*grid->getZC(i,j,k)/Lz) + sin(2.*M_PI*grid->getXC(i,j,k)/Lx));
			Bzc[i][j][k] = 2.*B0x*cos(2.*M_PI*grid->getZC(i,j,k)/Lz)*sin(2.*M_PI*grid->getYC(i,j,k)/Ly);
		}

	    // currents are used to calculate in the Maxwell's solver
	    // The ion current is equal to 0 (all current is on electrons)
	    for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
		for (int k=0; k < nzn; k++){
			Jxs[1][i][j][k] = 0.0; // ion species is species 1
			Jys[1][i][j][k] = 0.0; // ion species is species 1
			Jzs[1][i][j][k] = 0.0; // ion species is species 1
		}

	    // calculate the electron current from
        eqValue(0.0,tempXN,nxn,nyn,nzn);
        eqValue(0.0,tempYN,nxn,nyn,nzn);
        eqValue(0.0,tempZN,nxn,nyn,nzn);
        grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc); // here you calculate curl(B)
        // all current is on electrons, calculated from Ampere's law
        for (int i=0; i < nxn; i++)
        for (int j=0; j < nyn; j++)
        for (int k=0; k < nzn; k++){  // electrons are species 0
			Jxs[0][i][j][k] = c*tempXN[i][j][k]/FourPI; // ion species is species 1
			Jys[0][i][j][k] = c*tempYN[i][j][k]/FourPI; // ion species is species 1
			Jzs[0][i][j][k] = c*tempZN[i][j][k]/FourPI; // ion species is species 1
		}

		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init();  // use the fields from restart file
	}
}

void EMfields3D::initTaylorGreen()
{
	const VirtualTopology3D *vct = &get_vct();
	const Grid *grid = &get_grid();
	if (restart1 ==0){
		if (vct->getCartesian_rank() ==0){
			cout << "----------------------------------------" << endl;
			cout << "       Initialize Taylor-Green flow     " << endl;
			cout << "----------------------------------------" << endl;
			cout << "B0                               = " << B0x << endl;
            cout << "u0                               = " << ue0 << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i] << endl;
			}
			cout << "Smoothing Factor = " << Smooth << endl;
			cout << "-------------------------" << endl;
		}

        for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
		for (int k=0; k < nzn; k++){
		   // initialize the density for species
		   for (int is=0; is < ns; is++) {
             rhons[is][i][j][k] = rhoINIT[is]/FourPI;
             
             // The flow will be initialized from currents
             Jxs[is][i][j][k] = ue0 * rhons[is][i][j][k] * sin(2.*M_PI*grid->getXC(i,j,k)/Lx) * cos(2.*M_PI*grid->getYC(i,j,k)/Ly) * cos(2.*M_PI*grid->getZC(i,j,k)/Lz); 
             Jys[is][i][j][k] = -ue0 * rhons[is][i][j][k] * cos(2.*M_PI*grid->getXC(i,j,k)/Lx) * sin(2.*M_PI*grid->getYC(i,j,k)/Ly) * cos(2.*M_PI*grid->getZC(i,j,k)/Lz);
             Jzs[is][i][j][k] = 0.; // Z velocity is zero
           }

			// electric field
			Ex[i][j][k] =  0.0;
			Ey[i][j][k] =  0.0;
			Ez[i][j][k] =  0.0;
			// Magnetic field
			Bxn[i][j][k] = B0x * cos(2.*M_PI*grid->getXN(i,j,k)/Lx) * sin(2.*M_PI*grid->getYN(i,j,k)/Ly) * sin(2.*M_PI*grid->getZN(i,j,k)/Lz);
			Byn[i][j][k] = B0x * sin(2.*M_PI*grid->getXN(i,j,k)/Lx) * cos(2.*M_PI*grid->getYN(i,j,k)/Ly) * sin(2.*M_PI*grid->getZN(i,j,k)/Lz);
			Bzn[i][j][k] = -2. * B0x  * sin(2.*M_PI*grid->getXN(i,j,k)/Lx) * sin(2.*M_PI*grid->getYN(i,j,k)/Ly) * cos(2.*M_PI*grid->getZN(i,j,k)/Lz);
		}

		for (int i=0; i <nxc; i++)
		for (int j=0; j <nyc; j++)
		for (int k=0; k <nzc; k++) {
			Bxc[i][j][k] = B0x * cos(2.*M_PI*grid->getXC(i,j,k)/Lx) * sin(2.*M_PI*grid->getYC(i,j,k)/Ly) * sin(2.*M_PI*grid->getZC(i,j,k)/Lz);
			Byc[i][j][k] = B0x * sin(2.*M_PI*grid->getXC(i,j,k)/Lx) * cos(2.*M_PI*grid->getYC(i,j,k)/Ly) * sin(2.*M_PI*grid->getZC(i,j,k)/Lz);
			Bzc[i][j][k] = -2. * B0x * sin(2.*M_PI*grid->getXC(i,j,k)/Lx) * sin(2.*M_PI*grid->getYC(i,j,k)/Ly) * cos(2.*M_PI*grid->getZC(i,j,k)/Lz);
		}

		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init();  // use the fields from restart file
	}
}

void EMfields3D::initOriginalGEM()
{
  const Grid *grid = &get_grid();
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          const double yM = grid->getYN(i, j, k) - .5 * Ly;
          Bxn[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          Bxn[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byn[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          const double yM = grid->getYC(i, j, k) - .5 * Ly;
          Bxc[i][j][k] = B0x * tanh(yM / delta);
          // add the initial GEM perturbation
          const double xM = grid->getXC(i, j, k) - .5 * Lx;
          Bxc[i][j][k] -= (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * xM / Lx) * sin(M_PI * yM / Ly);
          Byc[i][j][k] = B0y + (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * xM / Lx) * cos(M_PI * yM / Ly);
          Bzc[i][j][k] = B0z;
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}

void EMfields3D::initGEMDoubleHarris()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  if (restart1 == 0){

    if (vct->getCartesian_rank() ==0){
      cout << "------------------------------------------" << endl;
      cout << "Initialize Double Harris Sheet with Perturbation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness)  = " << delta << endl;
      for (int i=0; i < ns; i++){
	cout << "rho species " << i <<" = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++){
	  const double xM = grid->getXN(i,j,k) - 0.5*Lx;
	  const double yB = grid->getYN(i,j,k) - 0.25*Ly;
	  const double yT = grid->getYN(i,j,k) - 0.75*Ly;
	  const double yBd = yB/delta;
	  const double yTd = yT/delta;
	  // initialize the density for species
	  for (int is=0; is < ns; is++){
	    if (DriftSpecies[is]){
	      const double sech_yBd = 1./cosh(yBd)+1e-5;
	      const double sech_yTd = 1./cosh(yTd)+1e-5;
	      rhons[is][i][j][k] = rhoINIT[is]*sech_yBd*sech_yBd/FourPI;
	      rhons[is][i][j][k]+= rhoINIT[is]*sech_yTd*sech_yTd/FourPI;
	    }
	    else
	      rhons[is][i][j][k] = rhoINIT[is]/FourPI;
	  }
	  // electric field
	  Ex[i][j][k] =  0.0;
	  Ey[i][j][k] =  0.0;
	  Ez[i][j][k] =  0.0;
	  // Magnetic field
	  Bxn[i][j][k] = B0x*(-1.0+tanh(yBd)-tanh(yTd));
	  Byn[i][j][k] = B0y;
	  Bzn[i][j][k] = B0z;
	  // add the initial X perturbation
	  xpert = grid->getXN(i, j, k) - Lx / 2;
	  ypert = yB;
	  exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
	  Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
	  Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);

	}

    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, 1, 1, 2, 2, 1, 1, vct, this);
    communicateNodeBC(nxn, nyn, nzn, Byn, 1, 1, 1, 1, 1, 1, vct, this);
    communicateNodeBC(nxn, nyn, nzn, Bzn, 1, 1, 2, 2, 1, 1, vct, this);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct, this);
    communicateCenterBC(nxc, nyc, nzc, Byc, 1, 1, 1, 1, 1, 1, vct, this);
    communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct, this);
    for (int is=0 ; is<ns; is++)
      grid->interpN2C(rhocs,is,rhons);
  } else {
    init();  // use the fields from restart file
  }
}

void EMfields3D::initDoublePeriodicHarrisWithGaussianHumpPerturbation()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  // perturbation localized in X
  const double pertX = 0.4;
  const double deltax = 8. * delta;
  const double deltay = 4. * delta;
  if (restart1 == 0) {
    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              const double sech_yBd = 1. / cosh(yBd);
              const double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxn[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byn[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxn[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byn[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);

          // guide field
          Bzn[i][j][k] = B0z;
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct, this);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct, this);

    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          const double xM = grid->getXN(i, j, k) - .5 * Lx;
          const double yB = grid->getYN(i, j, k) - .25 * Ly;
          const double yT = grid->getYN(i, j, k) - .75 * Ly;
          const double yBd = yB / delta;
          const double yTd = yT / delta;
          Bxc[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxc[i][j][k] += 0.;
          Byc[i][j][k] = B0y;
          // add the initial X perturbation
          const double xMdx = xM / deltax;
          const double yBdy = yB / deltay;
          const double yTdy = yT / deltay;
          const double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxc[i][j][k] -= (B0x * pertX) * humpB * (2.0 * yBdy);
          Byc[i][j][k] += (B0x * pertX) * humpB * (2.0 * xMdx);
          // add the second initial X perturbation
          const double humpT = exp(-xMdx * xMdx - yTdy * yTdy);
          Bxc[i][j][k] += (B0x * pertX) * humpT * (2.0 * yTdy);
          Byc[i][j][k] -= (B0x * pertX) * humpT * (2.0 * xMdx);
          // guide field
          Bzc[i][j][k] = B0z;
        }
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct,this);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct,this);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct,this);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}


/*! initialize GEM challenge with no Perturbation with dipole-like tail topology */
void EMfields3D::initGEMDipoleLikeTailNoPert()
{
  const Grid *grid = &get_grid();
  // parameters controling the field topology
  // e.g., x1=Lx/5,x2=Lx/4 give 'separated' fields, x1=Lx/4,x2=Lx/3 give 'reconnected' topology

  double x1 = Lx / 6.0;         // minimal position of the gaussian peak 
  double x2 = Lx / 4.0;         // maximal position of the gaussian peak (the one closer to the center)
  double sigma = Lx / 15;       // base sigma of the gaussian - later it changes with the grid
  double stretch_curve = 2.0;   // stretch the sin^2 function over the x dimension - also can regulate the number of 'knots/reconnecitons points' if less than 1
  double skew_parameter = 0.50; // skew of the shape of the gaussian
  double pi = 3.1415927;
  double r1, r2, delta_x1x2;

  if (restart1 == 0) {

    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYN(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYN(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXN(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field

          delta_x1x2 = x1 - x2 * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0 / stretch_curve) * (0.25 * FourPI) / 180.0));

          r1 = (grid->getYC(i, j, k) - (x1 + delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));
          r2 = (grid->getYC(i, j, k) - ((Lx - x1) - delta_x1x2)) * (1.0 - skew_parameter * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)) * (sin(((grid->getXC(i, j, k) - Lx / 2) / Lx * 180.0) * (0.25 * FourPI) / 180.0)));

          // tail-like field topology
          Bxn[i][j][k] = B0x * 0.5 * (-exp(-((r1) * (r1)) / (sigma * sigma)) + exp(-((r2) * (r2)) / (sigma * sigma)));

          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }

}

/*! initialize GEM challenge with no Perturbation */
void EMfields3D::initGEMnoPert()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  if (restart1 == 0) {

    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------------" << endl;
      cout << "Initialize GEM Challenge without Perturbation" << endl;
      cout << "----------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          Byn[i][j][k] = B0y;
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          Byc[i][j][k] = B0y;
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}

// new init, random problem
void EMfields3D::initRandomField()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  double **modes_seed = newArr2(double, 7, 7);
  if (restart1 ==0){
    // initialize
    if (get_vct().getCartesian_rank() ==0){
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i=0; i < ns; i++){
	cout << "rho species " << i <<" = " << rhoINIT[i];
	if (DriftSpecies[i])
	  cout << " DRIFTING " << endl;
	else
	  cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    double kx;
    double ky;
        
    /*       stringstream num_proc;
	     num_proc << vct->getCartesian_rank() ;
	     string cqsat = SaveDirName + "/RandomNumbers" + num_proc.str() + ".txt";
        ofstream my_file(cqsat.c_str(), fstream::binary);
	for (int m=-3; m < 4; m++)
            for (int n=-3; n < 4; n++){
            modes_seed[m+3][n+3] = rand() / (double) RAND_MAX;
            my_file <<"modes_seed["<< m+3<<"][" << "\t" << n+3 << "] = " << modes_seed[m+3][n+3] << endl;
            }
              my_file.close();
    */
    modes_seed[0][0] = 0.532767;
    modes_seed[0][1] = 0.218959;
    modes_seed[0][2] = 0.0470446;
    modes_seed[0][3] = 0.678865;
    modes_seed[0][4] = 0.679296;
    modes_seed[0][5] = 0.934693;
    modes_seed[0][6] = 0.383502;
    modes_seed[1][0] = 0.519416;
    modes_seed[1][1] = 0.830965;
    modes_seed[1][2] = 0.0345721;
    modes_seed[1][3] = 0.0534616;
    modes_seed[1][4] = 0.5297;
    modes_seed[1][5] = 0.671149;
    modes_seed[1][6] = 0.00769819;
    modes_seed[2][0] = 0.383416;
    modes_seed[2][1] = 0.0668422;
    modes_seed[2][2] = 0.417486;
    modes_seed[2][3] = 0.686773;
    modes_seed[2][4] = 0.588977;
    modes_seed[2][5] = 0.930436;
    modes_seed[2][6] = 0.846167;
    modes_seed[3][0] = 0.526929;
    modes_seed[3][1] = 0.0919649;
    modes_seed[3][2] = 0.653919;
    modes_seed[3][3] = 0.415999;
    modes_seed[3][4] = 0.701191;
    modes_seed[3][5] = 0.910321;
    modes_seed[3][6] = 0.762198;
    modes_seed[4][0] = 0.262453;
    modes_seed[4][1] = 0.0474645;
    modes_seed[4][2] = 0.736082;
    modes_seed[4][3] = 0.328234;
    modes_seed[4][4] = 0.632639;
    modes_seed[4][5] = 0.75641;
    modes_seed[4][6] = 0.991037;
    modes_seed[5][0] = 0.365339;
    modes_seed[5][1] = 0.247039;
    modes_seed[5][2] = 0.98255;
    modes_seed[5][3] = 0.72266;
    modes_seed[5][4] = 0.753356;
    modes_seed[5][5] = 0.651519;
    modes_seed[5][6] = 0.0726859;
    modes_seed[6][0] = 0.631635;
    modes_seed[6][1] = 0.884707;
    modes_seed[6][2] = 0.27271;
    modes_seed[6][3] = 0.436411;
    modes_seed[6][4] = 0.766495;
    modes_seed[6][5] = 0.477732;
    modes_seed[6][6] = 0.237774;

    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
	for (int k=0; k < nzn; k++){
	  // initialize the density for species
	  for (int is=0; is < ns; is++){
	    rhons[is][i][j][k] = rhoINIT[is]/FourPI;
	  }
	  // electric field
	  Ex[i][j][k] =  0.0;
	  Ey[i][j][k] =  0.0;
	  Ez[i][j][k] =  0.0;
	  // Magnetic field
	  Bxn[i][j][k] =  0.0;
	  Byn[i][j][k] =  0.0;
	  Bzn[i][j][k] =  B0z;
	  for (int m=-3; m < 4; m++)
	    for (int n=-3; n < 4; n++){

	      kx=2.0*M_PI*m/Lx;
	      ky=2.0*M_PI*n/Ly;
	      Bxn[i][j][k] += -B0x*ky*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	      Byn[i][j][k] += B0x*kx*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	      // Bzn[i][j][k] += B0x*cos(grid->getXN(i,j,k)*kx+grid->getYN(i,j,k)*ky+2.0*M_PI*modes_seed[m+3][n+3]);
	    }
	}
	  // communicate ghost
	  communicateNodeBC(nxn, nyn, nzn, Bxn, 1, 1, 2, 2, 1, 1, vct, this);
	  communicateNodeBC(nxn, nyn, nzn, Byn, 1, 1, 1, 1, 1, 1, vct, this);
	  communicateNodeBC(nxn, nyn, nzn, Bzn, 1, 1, 2, 2, 1, 1, vct, this);

	  // initialize B on centers
	  grid->interpN2C(Bxc, Bxn);
	  grid->interpN2C(Byc, Byn);
	  grid->interpN2C(Bzc, Bzn);
	  // communicate ghost
	  communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct,this);
	  communicateCenterBC(nxc, nyc, nzc, Byc, 1, 1, 1, 1, 1, 1, vct,this);
	  communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct,this);
	  for (int is=0 ; is<ns; is++)
            grid->interpN2C(rhocs,is,rhons);
	} else {
    init(); // use the fields from restart file
    }
  delArr2(modes_seed, 7);
}


/*! Init Force Free (JxB=0) */
void EMfields3D::initForceFree()
{
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();
  if (restart1 == 0) {

    // initialize
    if (get_vct().getCartesian_rank() == 0) {
      cout << "----------------------------------------" << endl;
      cout << "Initialize Force Free with Perturbation" << endl;
      cout << "----------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
      }
      cout << "Smoothing Factor = " << Smooth << endl;
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          Bxn[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXN(i, j, k) / Lx) * sin(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          Byn[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXN(i, j, k) / Lx) * cos(M_PI * (grid->getYN(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzn[i][j][k] = B0z / cosh((grid->getYN(i, j, k) - Ly / 2) / delta);
        }
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the perturbation
          Bxc[i][j][k] += (B0x / 10.0) * (M_PI / Ly) * cos(2 * M_PI * grid->getXC(i, j, k) / Lx) * sin(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          Byc[i][j][k] = B0y - (B0x / 10.0) * (2 * M_PI / Lx) * sin(2 * M_PI * grid->getXC(i, j, k) / Lx) * cos(M_PI * (grid->getYC(i, j, k) - Ly / 2) / Ly);
          // guide field
          Bzc[i][j][k] = B0z / cosh((grid->getYC(i, j, k) - Ly / 2) / delta);
        }

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}
/*! Initialize the EM field with constants values or from restart */
void EMfields3D::initBEAM(double x_center, double y_center, double z_center,
  double radius)
{
  const Grid *grid = &get_grid();

  double distance;
  // initialize E and rhos on nodes
  if (restart1 == 0) {
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = 0.0;
          Byn[i][j][k] = 0.0;
          Bzn[i][j][k] = 0.0;
          distance = (grid->getXN(i, j, k) - x_center) * (grid->getXN(i, j, k) - x_center) / (radius * radius) + (grid->getYN(i, j, k) - y_center) * (grid->getYN(i, j, k) - y_center) / (radius * radius) + (grid->getZN(i, j, k) - z_center) * (grid->getZN(i, j, k) - z_center) / (4 * radius * radius);
          // plasma
          rhons[0][i][j][k] = rhoINIT[0] / FourPI;  // initialize with constant density
          // electrons
          rhons[1][i][j][k] = rhoINIT[1] / FourPI;
          // beam
          if (distance < 1.0)
            rhons[2][i][j][k] = rhoINIT[2] / FourPI;
          else
            rhons[2][i][j][k] = 0.0;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = 0.0;
          Byc[i][j][k] = 0.0;
          Bzc[i][j][k] = 0.0;


        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(); // use the fields from restart file
  }
}

/*! Initialise a combination of magnetic dipoles */
void EMfields3D::initDipole()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  // initialize
  if (vct->getCartesian_rank() ==0){
      cout << "------------------------------------------" << endl;
      cout << "Initialise a Magnetic Dipole " << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "B1x   (external dipole field) - X  = " << B1x << endl;
      cout << "B1y                              = " << B1y << endl;
      cout << "B1z                              = " << B1z << endl;
      cout << "L_square - no magnetic field inside a sphere with radius L_square  = " << L_square << endl;
      cout << "Center dipole - X                = " << x_center << endl;
      cout << "Center dipole - Y                = " << y_center << endl;
      cout << "Center dipole - Z                = " << z_center << endl;
      cout << "Solar Wind drift velocity        = " << ue0 << endl;
  }


  double distance;
  double x_displ, y_displ, z_displ, fac1;

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){
        for (int is=0; is < ns; is++){
          rhons[is][i][j][k] = rhoINIT[is]/FourPI;
        }
        Ex[i][j][k] = ebc[0];
        Ey[i][j][k] = ebc[1];
        Ez[i][j][k] = ebc[2];

        double blp[3];
        // radius of the planet
        double a=L_square;

        double xc=x_center;
        double yc=y_center;
        double zc=z_center;

        double x = grid->getXN(i,j,k);
        double y = grid->getYN(i,j,k);
        double z = grid->getZN(i,j,k);

        double r2 = ((x-xc)*(x-xc)) + ((y-yc)*(y-yc)) + ((z-zc)*(z-zc));

        // Compute dipolar field B_ext

        if (r2 > a*a) {
            x_displ = x - xc;
            y_displ = y - yc;
            z_displ = z - zc;
            fac1 =  -B1z*a*a*a/pow(r2,2.5);
	    Bx_ext[i][j][k] = 3*x_displ*z_displ*fac1;
	    By_ext[i][j][k] = 3*y_displ*z_displ*fac1;
	    Bz_ext[i][j][k] = (2*z_displ*z_displ -x_displ*x_displ -y_displ*y_displ)*fac1;
        }
        else { // no field inside the planet
            Bx_ext[i][j][k]  = 0.0;
            By_ext[i][j][k]  = 0.0;
            Bz_ext[i][j][k]  = 0.0;
        }
        Bxn[i][j][k] = B0x;// + Bx_ext[i][j][k]
        Byn[i][j][k] = B0y;// + By_ext[i][j][k]
        Bzn[i][j][k] = B0z;// + Bz_ext[i][j][k]

      }
    }
  }

	grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);dprintf("1 Bzc[1][15][0]=%f",Bzc[1][15][0]);

	communicateCenterBC_P(nxc,nyc,nzc,Bxc,col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5],vct, this);
	communicateCenterBC_P(nxc,nyc,nzc,Byc,col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5],vct, this);
	communicateCenterBC_P(nxc,nyc,nzc,Bzc,col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5],vct, this);

	for (int is=0 ; is<ns; is++)
		grid->interpN2C(rhocs,is,rhons);

	if (restart1 != 0) { // EM initialization from RESTART
		init();  // use the fields from restart file
	}

}

/*! Initialise a 2D magnetic dipoles according to paper L.K.S Two-way coupling of a global Hall ....*/
void EMfields3D::initDipole2D()
{
  const Collective *col = &get_col();
  const VirtualTopology3D *vct = &get_vct();
  const Grid *grid = &get_grid();

  // initialize
  if (vct->getCartesian_rank() ==0){
      cout << "------------------------------------------" << endl;
      cout << "Initialise a 2D Magnetic Dipole on XY Plane" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "B1x   (external dipole field)    = " << B1x << endl;
      cout << "B1y                              = " << B1y << endl;
      cout << "B1z                              = " << B1z << endl;
      cout << "L_square - no magnetic field inside a sphere with radius L_square  = " << L_square << endl;
      cout << "Center dipole - X                = " << x_center << endl;
      cout << "Center dipole - Y                = " << y_center << endl;
      cout << "Center dipole - Z                = " << z_center << endl;
      cout << "Solar Wind drift velocity        = " << ue0 << endl;
      cout << "2D Smoothing Factor              = " << Smooth << endl;
      cout << "Smooth Iteration                 = " << SmoothNiter << endl;
  }


  double distance;
  double x_displ, z_displ, fac1;

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){
        for (int is=0; is < ns; is++){
          rhons[is][i][j][k] = rhoINIT[is]/FourPI;
        }
        Ex[i][j][k] = ebc[0];
        Ey[i][j][k] = ebc[1];
        Ez[i][j][k] = ebc[2];

        double blp[3];
        double a=L_square;

        double xc=x_center;
        double zc=z_center;

        double x = grid->getXN(i,j,k);
        double z = grid->getZN(i,j,k);

        double r2 = ((x-xc)*(x-xc)) + ((z-zc)*(z-zc));

        // Compute dipolar field B_ext

        if (r2 > a*a) {
            x_displ = x - xc;
            z_displ = z - zc;

            fac1 =  -B1z*a*a/(r2*r2);//fac1 = D/4?

			Bx_ext[i][j][k] = 2*x_displ*z_displ*fac1;
			By_ext[i][j][k] = 0.0;
			Bz_ext[i][j][k] = (z_displ*z_displ - x_displ*x_displ)*fac1;

        }
        else { // no field inside the planet
            Bx_ext[i][j][k]  = 0.0;
            By_ext[i][j][k]  = 0.0;
            Bz_ext[i][j][k]  = 0.0;
        }

        Bxn[i][j][k] = B0x;// + Bx_ext[i][j][k]
        Byn[i][j][k] = B0y;// + By_ext[i][j][k]
        Bzn[i][j][k] = B0z;// + Bz_ext[i][j][k]

      }
    }
  }


	grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);

	communicateCenterBC_P(nxc,nyc,nzc,Bxc,col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5],vct, this);
	communicateCenterBC_P(nxc,nyc,nzc,Byc,col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5],vct, this);
	communicateCenterBC_P(nxc,nyc,nzc,Bzc,col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5],vct, this);

	for (int is=0 ; is<ns; is++)
		grid->interpN2C(rhocs,is,rhons);



	if (restart1 != 0) { // EM initialization from RESTART
		init();  // use the fields from restart file
	}
}

/*! Calculate the susceptibility on the boundary leftX */
void EMfields3D::sustensorLeftX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[1][j][k]+Bx_ext[1][j][k]);
        omcy = beta * (Byn[1][j][k]+By_ext[1][j][k]);
        omcz = beta * (Bzn[1][j][k]+Bz_ext[1][j][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][1][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary rightX */
void EMfields3D::sustensorRightX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[nxn - 2][j][k]+Bx_ext[nxn - 2][j][k]);
        omcy = beta * (Byn[nxn - 2][j][k]+By_ext[nxn - 2][j][k]);
        omcz = beta * (Bzn[nxn - 2][j][k]+Bz_ext[nxn - 2][j][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][nxn - 2][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[i][1][k]+Bx_ext[i][1][k]);
        omcy = beta * (Byn[i][1][k]+By_ext[i][1][k]);
        omcz = beta * (Bzn[i][1][k]+Bz_ext[i][1][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][1][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[i][nyn - 2][k]+Bx_ext[i][nyn - 2][k]);
        omcy = beta * (Byn[i][nyn - 2][k]+By_ext[i][nyn - 2][k]);
        omcz = beta * (Bzn[i][nyn - 2][k]+Bz_ext[i][nyn - 2][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][nyn - 2][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * (Bxn[i][j][1]+Bx_ext[i][j][1]);
        omcy = beta * (Byn[i][j][1]+By_ext[i][j][1]);
        omcz = beta * (Bzn[i][j][1]+Bz_ext[i][j][1]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][1] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * (Bxn[i][j][nzn - 2]+Bx_ext[i][j][nzn - 2]);
        omcy = beta * (Byn[i][j][nzn - 2]+By_ext[i][j][nzn - 2]);
        omcz = beta * (Bzn[i][j][nzn - 2]+Bz_ext[i][j][nzn - 2]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][nyn - 2] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }
}

/*! Perfect conductor boundary conditions: LEFT wall */
void EMfields3D::perfectConductorLeft(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
  int dir)
{
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0:  // boundary condition on X-DIRECTION 
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorLeftX(susxx, susyx, suszx);
      for (int i=1; i <  nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[1][i][j] = vectorX.get(1,i,j) - (Ex[1][i][j] - susyx[i][j]*vectorY.get(1,i,j) - suszx[i][j]*vectorZ.get(1,i,j) - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[1][i][j] = vectorY.get(1,i,j) - 0.0*vectorY.get(2,i,j);
          imageZ[1][i][j] = vectorZ.get(1,i,j) - 0.0*vectorZ.get(2,i,j);
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorLeftY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[i][1][j] = vectorX.get(i,1,j) - 0.0*vectorX.get(i,2,j);
          imageY[i][1][j] = vectorY.get(i,1,j) - (Ey[i][1][j] - susxy[i][j]*vectorX.get(i,1,j) - suszy[i][j]*vectorZ.get(i,1,j) - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][1][j] = vectorZ.get(i,1,j) - 0.0*vectorZ.get(i,2,j);
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorLeftZ(susxz, susyz, suszz);
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          imageX[i][j][1] = vectorX.get(i,j,1);
          imageY[i][j][1] = vectorX.get(i,j,1);
          imageZ[i][j][1] = vectorZ.get(i,j,1) - (Ez[i][j][1] - susxz[i][j]*vectorX.get(i,j,1) - susyz[i][j]*vectorY.get(i,j,1) - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions: RIGHT wall */
void EMfields3D::perfectConductorRight(
  arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX,
  const_arr3_double vectorY,
  const_arr3_double vectorZ,
  int dir)
{
  double beta, omcx, omcy, omcz, denom;
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorRightX(susxx, susyx, suszx);
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[nxn-2][i][j] = vectorX.get(nxn-2,i,j) - (Ex[nxn-2][i][j] - susyx[i][j]*vectorY.get(nxn-2,i,j) - suszx[i][j]*vectorZ.get(nxn-2,i,j) - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[nxn-2][i][j] = vectorY.get(nxn-2,i,j) - 0.0 * vectorY.get(nxn-3,i,j);
          imageZ[nxn-2][i][j] = vectorZ.get(nxn-2,i,j) - 0.0 * vectorZ.get(nxn-3,i,j);
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);       
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorRightY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          imageX[i][nyn-2][j] = vectorX.get(i,nyn-2,j) - 0.0*vectorX.get(i,nyn-3,j);
          imageY[i][nyn-2][j] = vectorY.get(i,nyn-2,j) - (Ey[i][nyn-2][j] - susxy[i][j]*vectorX.get(i,nyn-2,j) - suszy[i][j]*vectorZ.get(i,nyn-2,j) - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][nyn-2][j] = vectorZ.get(i,nyn-2,j) - 0.0*vectorZ.get(i,nyn-3,j);
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION RIGHT
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorRightZ(susxz, susyz, suszz);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nyn-1;j++){
          imageX[i][j][nzn-2] = vectorX.get(i,j,nzn-2);
          imageY[i][j][nzn-2] = vectorY.get(i,j,nzn-2);
          imageZ[i][j][nzn-2] = vectorZ.get(i,j,nzn-2) - (Ez[i][j][nzn-2] - susxz[i][j]*vectorX.get(i,j,nzn-2) - susyz[i][j]*vectorY.get(i,j,nzn-2) - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);       
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::perfectConductorLeftS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION LEFT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[1][i][j] = 0.0;
          vectorY[1][i][j] = ebc[1];
          vectorZ[1][i][j] = ebc[2];
          //+//          vectorX[1][i][j] = 0.0;
          //+//          vectorY[1][i][j] = 0.0;
          //+//          vectorZ[1][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][1][j] = ebc[0];
          vectorY[i][1][j] = 0.0;
          vectorZ[i][1][j] = ebc[2];
          //+//          vectorX[i][1][j] = 0.0;
          //+//          vectorY[i][1][j] = 0.0;
          //+//          vectorZ[i][1][j] = 0.0;
        }
      break;
    case 2: // boundary condition on Z-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][1] = ebc[0];
          vectorY[i][j][1] = ebc[1];
          vectorZ[i][j][1] = 0.0;
          //+//          vectorX[i][j][1] = 0.0;
          //+//          vectorY[i][j][1] = 0.0;
          //+//          vectorZ[i][j][1] = 0.0;
        }
      break;
  }
}

/*! Perfect conductor boundary conditions for source: RIGHT WALL */
void EMfields3D::perfectConductorRightS(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[nxn-2][i][j] = 0.0;
          vectorY[nxn-2][i][j] = ebc[1];
          vectorZ[nxn-2][i][j] = ebc[2];
          //+//          vectorX[nxn-2][i][j] = 0.0;
          //+//          vectorY[nxn-2][i][j] = 0.0;
          //+//          vectorZ[nxn-2][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][nyn-2][j] = ebc[0];
          vectorY[i][nyn-2][j] = 0.0;
          vectorZ[i][nyn-2][j] = ebc[2];
          //+//          vectorX[i][nyn-2][j] = 0.0;
          //+//          vectorY[i][nyn-2][j] = 0.0;
          //+//          vectorZ[i][nyn-2][j] = 0.0;
        }
      break;
    case 2:
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][nzn-2] = ebc[0];
          vectorY[i][j][nzn-2] = ebc[1];
          vectorZ[i][j][nzn-2] = 0.0;
          //+//          vectorX[i][j][nzn-2] = 0.0;
          //+//          vectorY[i][j][nzn-2] = 0.0;
          //+//          vectorZ[i][j][nzn-2] = 0.0;
        }
      break;
  }
}


void EMfields3D::OpenBoundaryInflowEImage(arr3_double imageX, arr3_double imageY, arr3_double imageZ,
  const_arr3_double vectorX, const_arr3_double vectorY, const_arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();
  // Assuming E = - ve x B
  double injE[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,injE);
  scale(injE,-1.0,3);

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[0][j][k] = vectorX[0][j][k] - injE[0];
        imageY[0][j][k] = vectorY[0][j][k] - injE[1];
        imageZ[0][j][k] = vectorZ[0][j][k] - injE[2];
      }
  }
  /*
  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright == 2) {
    for (int j=1; j < ny-1;j++)
      for (int k=1; k < nz-1;k++){
        imageX[nx-1][j][k] = vectorX[nx-1][j][k]- injE[0];
        imageY[nx-1][j][k] = vectorY[nx-1][j][k]- injE[1];
        imageZ[nx-1][j][k] = vectorZ[nx-1][j][k]- injE[2];

      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][0][k] = vectorX[i][0][k]-injE[0];
        imageY[i][0][k] = vectorY[i][0][k]-injE[1];
        imageZ[i][0][k] = vectorZ[i][0][k]-injE[2];
      }

  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int k=1; k < nz-1;k++){
        imageX[i][ny-1][k] = vectorX[i][ny-1][k]-injE[0];
        imageY[i][ny-1][k] = vectorY[i][ny-1][k]-injE[1];
        imageZ[i][ny-1][k] = vectorZ[i][ny-1][k]-injE[2];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][0] = vectorX[i][j][0]-injE[0];
        imageY[i][j][0] = vectorY[i][j][0]-injE[1];
        imageZ[i][j][0] = vectorZ[i][j][0]-injE[2];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=1; i < nx-1;i++)
      for (int j=1; j < ny-1;j++){
        imageX[i][j][nz-1] = vectorX[i][j][nz-1]-injE[0];
        imageY[i][j][nz-1] = vectorY[i][j][nz-1]-injE[1];
        imageZ[i][j][nz-1] = vectorZ[i][j][nz-1]-injE[2];
      }
  }
  */
}

void EMfields3D::OpenBoundaryInflowB(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2 && nx>10) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
          
	vectorX[0][j][k] = B0x;
        vectorY[0][j][k] = B0y;
        vectorZ[0][j][k] = B0z;

	vectorX[1][j][k] = B0x;
	vectorY[1][j][k] = B0y;
	vectorZ[1][j][k] = B0z;
		
	vectorX[2][j][k] = B0x;
	vectorY[2][j][k] = B0y;
	vectorZ[2][j][k] = B0z;

	vectorX[3][j][k] = B0x;
	vectorY[3][j][k] = B0y;
	vectorZ[3][j][k] = B0z;

      }
  }

  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2 && nx>10 ) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){

        vectorX[nx-4][j][k] = vectorX[nx-5][j][k];
        vectorY[nx-4][j][k] = vectorY[nx-5][j][k];
        vectorZ[nx-4][j][k] = vectorZ[nx-5][j][k];

        vectorX[nx-3][j][k] = vectorX[nx-5][j][k];
        vectorY[nx-3][j][k] = vectorY[nx-5][j][k];
        vectorZ[nx-3][j][k] = vectorZ[nx-5][j][k];

        vectorX[nx-2][j][k] = vectorX[nx-5][j][k];
        vectorY[nx-2][j][k] = vectorY[nx-5][j][k];
        vectorZ[nx-2][j][k] = vectorZ[nx-5][j][k];

        vectorX[nx-1][j][k] = vectorX[nx-5][j][k];
        vectorY[nx-1][j][k] = vectorY[nx-5][j][k];
        vectorZ[nx-1][j][k] = vectorZ[nx-5][j][k];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2 && ny> 10)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){

    	  vectorX[i][0][k] = vectorX[i][4][k];
    	  vectorY[i][0][k] = vectorY[i][4][k];
    	  vectorZ[i][0][k] = vectorZ[i][4][k];

    	  vectorX[i][1][k] = vectorX[i][4][k];
    	  vectorY[i][1][k] = vectorY[i][4][k];
    	  vectorZ[i][1][k] = vectorZ[i][4][k];

    	  vectorX[i][2][k] = vectorX[i][4][k];
    	  vectorY[i][2][k] = vectorY[i][4][k];
    	  vectorZ[i][2][k] = vectorZ[i][4][k];

    	  vectorX[i][3][k] = vectorX[i][4][k];
    	  vectorY[i][3][k] = vectorY[i][4][k];
    	  vectorZ[i][3][k] = vectorZ[i][4][k];
      } 
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==2 && ny>10)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k< nz;k++){

    	vectorX[i][ny-4][k] = vectorX[i][ny-5][k];
        vectorY[i][ny-4][k] = vectorY[i][ny-5][k];
        vectorZ[i][ny-4][k] = vectorZ[i][ny-5][k];

        vectorX[i][ny-3][k] = vectorX[i][ny-5][k];
        vectorY[i][ny-3][k] = vectorY[i][ny-5][k];
        vectorZ[i][ny-3][k] = vectorZ[i][ny-5][k];

        vectorX[i][ny-2][k] = vectorX[i][ny-5][k];
        vectorY[i][ny-2][k] = vectorY[i][ny-5][k];
        vectorZ[i][ny-2][k] = vectorZ[i][ny-5][k];

        vectorX[i][ny-1][k] = vectorX[i][ny-5][k];
        vectorY[i][ny-1][k] = vectorY[i][ny-5][k];
        vectorZ[i][ny-1][k] = vectorZ[i][ny-5][k];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2 && nz > 10)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){

    	  vectorX[i][j][0] = vectorX[i][j][4];
    	  vectorY[i][j][0] = vectorY[i][j][4];
    	  vectorZ[i][j][0] = vectorZ[i][j][4];

    	  vectorX[i][j][1] = vectorX[i][j][4];
    	  vectorY[i][j][1] = vectorY[i][j][4];
    	  vectorZ[i][j][1] = vectorZ[i][j][4];

    	  vectorX[i][j][2] = vectorX[i][j][4];
    	  vectorY[i][j][2] = vectorY[i][j][4];
    	  vectorZ[i][j][2] = vectorZ[i][j][4];

    	  vectorX[i][j][3] = vectorX[i][j][4];
    	  vectorY[i][j][3] = vectorY[i][j][4];
    	  vectorZ[i][j][3] = vectorZ[i][j][4];

      } 
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2 && nz>10)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){

    	vectorX[i][j][nz-4] = vectorX[i][j][nz-5];
        vectorY[i][j][nz-4] = vectorY[i][j][nz-5];
        vectorZ[i][j][nz-4] = vectorZ[i][j][nz-5];

        vectorX[i][j][nz-3] = vectorX[i][j][nz-5];
        vectorY[i][j][nz-3] = vectorY[i][j][nz-5];
        vectorZ[i][j][nz-3] = vectorZ[i][j][nz-5];

        vectorX[i][j][nz-2] = vectorX[i][j][nz-5];
        vectorY[i][j][nz-2] = vectorY[i][j][nz-5];
        vectorZ[i][j][nz-2] = vectorZ[i][j][nz-5];

        vectorX[i][j][nz-1] = vectorX[i][j][nz-5];
        vectorY[i][j][nz-1] = vectorY[i][j][nz-5];
        vectorZ[i][j][nz-1] = vectorZ[i][j][nz-5];
      }
  }

  /*
  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
		vectorX[nx-2][j][k] = B0x;
		vectorY[nx-2][j][k] = B0y;
		vectorZ[nx-2][j][k] = B0z;

        vectorX[nx-1][j][k] = B0x;
        vectorY[nx-1][j][k] = B0y;
        vectorZ[nx-1][j][k] = B0z;
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
		vectorX[i][1][k] = B0x;
		vectorY[i][1][k] = B0y;
		vectorZ[i][1][k] = B0z;

        vectorX[i][0][k] = B0x;
        vectorY[i][0][k] = B0y;
        vectorZ[i][0][k] = B0z;
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2)  {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
		vectorX[i][ny-2][k] = B0x;
		vectorY[i][ny-2][k] = B0y;
		vectorZ[i][ny-2][k] = B0z;

        vectorX[i][ny-1][k] = B0x;
        vectorY[i][ny-1][k] = B0y;
        vectorZ[i][ny-1][k] = B0z;
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
		vectorX[i][j][1] = B0x;
		vectorY[i][j][1] = B0y;
		vectorZ[i][j][1] = B0z;

        vectorX[i][j][0] = B0x;
        vectorY[i][j][0] = B0y;
        vectorZ[i][j][0] = B0z;
      }
  }


  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2)  {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
		vectorX[i][j][nz-2] = B0x;
		vectorY[i][j][nz-2] = B0y;
		vectorZ[i][j][nz-2] = B0z;

        vectorX[i][j][nz-1] = B0x;
        vectorY[i][j][nz-1] = B0y;
        vectorZ[i][j][nz-1] = B0z;
      }
  }
  */
}

void EMfields3D::OpenBoundaryInflowE(arr3_double vectorX, arr3_double vectorY, arr3_double vectorZ,
  int nx, int ny, int nz)
{
  const VirtualTopology3D *vct = &get_vct();
  // Assuming E = - ve x B
  double injE[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,injE);
  scale(injE,-1.0,3);

  if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
	//vectorX[0][j][k] = injE[0];
	//vectorY[0][j][k] = injE[1];
	//vectorZ[0][j][k] = injE[2];

        vectorX[1][j][k] = injE[0];
        vectorY[1][j][k] = injE[1];
        vectorZ[1][j][k] = injE[2];


        vectorX[2][j][k] = injE[0];
        vectorY[2][j][k] = injE[1];
        vectorZ[2][j][k] = injE[2];


        vectorX[3][j][k] = injE[0];
        vectorY[3][j][k] = injE[1];
        vectorZ[3][j][k] = injE[2];
      } 
  }
  /*
  if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==2) {
    for (int j=0; j < ny;j++)
      for (int k=0; k < nz;k++){
	//vectorX[nx-2][j][k] = injE[0];
	//	vectorY[nx-2][j][k] = injE[1];
	//	vectorZ[nx-2][j][k] = injE[2];

        vectorX[nx-1][j][k] = injE[0];
        vectorY[nx-1][j][k] = injE[1];
        vectorZ[nx-1][j][k] = injE[2];
      }
  }

  if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//vectorX[i][1][k] = injE[0];
	//ectorY[i][1][k] = injE[1];
	//vectorZ[i][1][k] = injE[2];

        vectorX[i][0][k] = injE[0];
        vectorY[i][0][k] = injE[1];
        vectorZ[i][0][k] = injE[2];
      }
  }

  if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==2) {
    for (int i=0; i < nx;i++)
      for (int k=0; k < nz;k++){
	//vectorX[i][ny-2][k] = injE[0];
	//vectorY[i][ny-2][k] = injE[1];
	//vectorZ[i][ny-2][k] = injE[2];

        vectorX[i][ny-1][k] = injE[0];
        vectorY[i][ny-1][k] = injE[1];
        vectorZ[i][ny-1][k] = injE[2];
      }
  }

  if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//vectorX[i][j][1] = injE[0];
	//vectorY[i][j][1] = injE[1];
	//vectorZ[i][j][1] = injE[2];

        vectorX[i][j][0] = injE[0];
        vectorY[i][j][0] = injE[1];
        vectorZ[i][j][0] = injE[2];
      }
  }

  if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright ==2) {
    for (int i=0; i < nx;i++)
      for (int j=0; j < ny;j++){
	//vectorX[i][j][nz-2] = injE[0];
	//vectorY[i][j][nz-2] = injE[1];
	//vectorZ[i][j][nz-2] = injE[2];

        vectorX[i][j][nz-1] = injE[0];
        vectorY[i][j][nz-1] = injE[1];
        vectorZ[i][j][nz-1] = injE[2];
      }
   }*/
}

/*! get Electric Field component X array cell without the ghost cells */
//arr3_double EMfields3D::getExc()
//{
//  array3_double tmp(nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, Ex);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[i][j][k];
//  return arr;
//}
/*! get Electric Field component Y array cell without the ghost cells */
//arr3_double EMfields3D::getEyc()
//{
//  array3_double tmp(nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, Ey);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[i][j][k];
//  return arr;
//}
/*! get Electric Field component Z array cell without the ghost cells */
//arr3_double EMfields3D::getEzc()
//{
//  array3_double tmp(nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, Ez);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[i][j][k];
//  return arr;
//}
/*! get Magnetic Field component X array cell without the ghost cells */
//arr3_double EMfields3D::getBxc() {
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=Bxc[i][j][k];
//  return arr;
//}
/*! get Magnetic Field component Y array cell without the ghost cells */
//arr3_double EMfields3D::getByc() {
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=Byc[i][j][k];
//  return arr;
//}
/*! get Magnetic Field component Z array cell without the ghost cells */
//arr3_double EMfields3D::getBzc() {
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=Bzc[i][j][k];
//  return arr;
//}
/*! get species density component X array cell without the ghost cells */
//arr3_double EMfields3D::getRHOcs(int is)
//{
//  array4_double tmp(ns,nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, is, rhons);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[is][i][j][k];
//  return arr;
//}

/*! get Magnetic Field component X array species is cell without the ghost cells */
//arr3_double EMfields3D::getJxsc(int is)
//{
//  array4_double tmp(ns,nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, is, Jxs);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[is][i][j][k];
//  return arr;
//}

/*! get current component Y array species is cell without the ghost cells */
//arr3_double EMfields3D::getJysc(int is)
//{
//  array4_double tmp(ns,nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, is, Jys);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[is][i][j][k];
//  return arr;
//}
/*! get current component Z array species is cell without the ghost cells */
//arr3_double EMfields3D::getJzsc(int is)
//{
//  array4_double tmp(ns,nxc,nyc,nzc);
//  get_grid().interpN2C(tmp, is, Jzs);
//
//  for (int i = 1; i < nxc-1; i++)
//    for (int j = 1; j < nyc-1; j++)
//      for (int k = 1; k < nzc-1; k++)
//        arr[i-1][j-1][k-1]=tmp[is][i][j][k];
//  return arr;
//}
/*! get the electric field energy */
double EMfields3D::getEenergy(void) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++)
        localEenergy += .5 * dx * dy * dz * (Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI);

  MPI_Allreduce(&localEenergy, &totalEenergy, 1, MPI_DOUBLE, MPI_SUM, (&get_vct())->getFieldComm());
  return (totalEenergy);

}
/*! get the magnetic field energy */
double EMfields3D::getBenergy(void) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;
  for (int i = 1; i < nxc - 1; i++)
    for (int j = 1; j < nyc - 1; j++)
      for (int k = 1; k < nzc - 1; k++){
        Bxt = Bxc[i][j][k];
        Byt = Byc[i][j][k];
        Bzt = Bzc[i][j][k];
        localBenergy += .5*dx*dy*dz*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM,(&get_vct())->getFieldComm());
  return (totalBenergy);
}

/*! get bulk kinetic energy*/
double EMfields3D::getBulkEnergy(int is) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  for (int i = 1; i < nxn - 2; i++)
    for (int j = 1; j < nyn - 2; j++)
      for (int k = 1; k < nzn - 2; k++){
	if(rhons[is][i][j][k] !=0)
	  localBenergy += .5 * dx * dy * dz * (Jxs[is][i][j][k] * Jxs[is][i][j][k] + Jys[is][i][j][k] * Jys[is][i][j][k] + Jzs[is][i][j][k] * Jzs[is][i][j][k]) / (rhons[is][i][j][k]);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, (&get_vct())->getFieldComm());
  return (totalBenergy);
}



/*! Print info about electromagnetic field */
void EMfields3D::print(void) const {
}

/*! destructor*/
EMfields3D::~EMfields3D() {
  delete [] qom;
  delete [] rhoINIT;
  delete [] DriftSpecies;
  for(int i=0;i<sizeMomentsArray;i++) { delete moments10Array[i]; }
  for(int i=0;i<sizeMomentsArray;i++) { delete moments13Array[i]; }
  delete [] moments10Array;
  delete [] moments13Array;
  delArr5(M_GII, nxn, nyn, nzn, ngp);
  delArr4(M_CI,  nxc, nyc, nzc);
  freeDataType();
}
