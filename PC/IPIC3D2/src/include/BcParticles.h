/***************************************************************************
  BcParticles.h  -  Library to manage boundary conditions for particles
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef BcParticles_H
#define BcParticles_H

//#include <ipicmath.h>

//#include "VirtualTopology3D.h"
//#include "Basic.h"
//#include "errors.h"

namespace BCparticles
{
  enum Enum
  {
    EXIT = 0,
    PERFECT_MIRROR = 1,
    REEMISSION = 2,
    OPENBCOut = 3,
    OPENBCIn = 4
  };
}

//inline void BCpclLeft(double& x, double& u, double& v, double& w, double Lx, double ut, double vt, double wt, int bcFaceXleft)
//{
//  assert_lt(x, 0.);
//  switch (bcFaceXleft)
//  {
//    default:
//      unsupported_value_error(bcFaceXleft);
//    case BCparticles::PERFECT_MIRROR:
//      x = -x;
//      u = -u;
//      break;
//    case BCparticles::REEMISSION:
//      x = -x;
//      sample_maxwellian(u,v,w, ut,vt,wt);
//      u = fabs(u);
//      break;
//    case BCparticles::EXIT:
//      break;
//  }
//}
//inline void BCpclRight(double& x, double& u, double& v, double& w, double Lx, double ut, double vt, double wt, int bcFaceXright)
//{
//  assert_gt(x, Lx);
//  switch (bcFaceXright)
//  {
//    default:
//      unsupported_value_error(bcFaceXright);
//    case BCparticles::PERFECT_MIRROR:
//      x = 2 * Lx - x;
//      u = -u;
//      break;
//    case BCparticles::REEMISSION:
//      x = 2 * Lx - x;
//      sample_maxwellian(u,v,w, ut,vt,wt);
//      u = -fabs(u);
//      break;
//    case BCparticles::EXIT:
//      break;
//  }
//}

///** set the boundary condition for a particle */
//void BCpclLeft(double *x, double *u, double Lx, double ut, int bcFaceXleft)
//{
//  assert_lt(*x,0);
//  switch (bcFaceXleft)
//  {
//    case BCparticles::PERFECT_MIRROR:
//      MODULO(x, Lx);
//      *u = -*u;
//      break;
//    case BCparticles::REEMISSION:
//      MODULO(x, Lx);
//      sample_maxwellian(*u,ut);
//      *u = fabs(*u);
//      break;
//  }
//}
//void BCpclRght(double *x, double *u, double Lx, double ut, int bcFaceXright)
//{
//  assert_gt(*x, Lx);
//  switch (bcFaceXright)
//  {
//    case BCparticles::PERFECT_MIRROR:
//      *x = 2 * Lx - *x;
//      *u = -*u;
//      break;
//    case BCparticles::REEMISSION:
//      *x = 2 * Lx - *x;
//      sample_maxwellian(*u,ut);
//      *u = -fabs(*u);
//  }
//}
//void BCpart(double *x, double *u, double Lx, double ut, int bcFaceXright, int bcFaceXleft)
//{
//  if (*x > Lx) {
//    switch (bcFaceXright) {
//      case BCparticles::PERFECT_MIRROR:
//        *x = 2 * Lx - *x;
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        *x = 2 * Lx - *x;
//        sample_maxwellian(*u,ut);
//        *u = -fabs(*u);
//    }
//  }
//  else if (*x < 0) {
//    switch (bcFaceXleft) {
//      case BCparticles::PERFECT_MIRROR:
//        MODULO(x, Lx);
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        MODULO(x, Lx);
//        sample_maxwellian(*u,ut);
//        *u = fabs(*u);
//        break;
//    }
//  }
//}
//void BCpart(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt, int bcFaceXright, int bcFaceXleft)
//{
//  if (*x > Lx) {
//    switch (bcFaceXright) {
//      case BCparticles::PERFECT_MIRROR:
//        *x = 2 * Lx - *x;
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        *x = 2 * Lx - *x;
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *u = -fabs(*u);
//    }
//  }
//  else if (*x < 0) {
//    switch (bcFaceXleft) {
//      case BCparticles::PERFECT_MIRROR:
//        *x = -*x;
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        *x = -*x;
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *u = fabs(*u);
//        break;
//    }
//  }
//}
///** set the boundary condition on boundaries for particle in 3D*/
//void BCpart(double *x, double *y, double *z, double *u, double *v, double *w, double Lx, double Ly, double Lz, double ut, double vt, double wt, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct)
//{
//  if (*x > Lx && vct->getXright_neighbor() == MPI_PROC_NULL) {
//    switch (bcFaceXright) {
//      case BCparticles::PERFECT_MIRROR:
//        MODULO(x, Lx);
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        double harvest, prob, theta;
//        MODULO(x, Lx);
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *u = -fabs(*u);
//        break;
//    }
//  }
//  if (*x < 0 && vct->getXleft_neighbor() == MPI_PROC_NULL) {
//    switch (bcFaceXleft) {
//      case BCparticles::PERFECT_MIRROR:
//        MODULO(x, Lx);
//        *u = -*u;
//        break;
//      case BCparticles::REEMISSION:
//        MODULO(x, Lx);
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *u = fabs(*u);
//        break;
//    }
//  }
//  if (*y > Ly && vct->getYright_neighbor() == MPI_PROC_NULL) {
//    switch (bcFaceYright) {
//      case BCparticles::PERFECT_MIRROR:
//        MODULO(y, Ly);
//        *v = -*v;
//        break;
//      case BCparticles::REEMISSION:
//        MODULO(y, Ly);
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *v = -fabs(*v);
//        break;
//    }
//  }
//  if (*y < 0 && vct->getYleft_neighbor() == MPI_PROC_NULL) {
//    switch (bcFaceYleft) {
//      case BCparticles::PERFECT_MIRROR:
//        MODULO(y, Ly);
//        *v = -*v;
//        break;
//      case BCparticles::REEMISSION:
//        MODULO(y, Ly);
//        sample_maxwellian(*u,*v,*w, ut,vt,wt);
//        *v = fabs(*v);
//        break;
//    }
//  }
//}
#endif
