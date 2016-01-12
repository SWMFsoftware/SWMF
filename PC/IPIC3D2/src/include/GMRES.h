#ifndef GMRES_new2_H
#define GMRES_new2_H

#include "ipicfwd.h"

typedef void (EMfields3D::*FIELD_IMAGE) (double *, double *, bool);
typedef void (*GENERIC_IMAGE) (double *, double *);

void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen, const double *b, int m, int max_iter, double tol,bool doSolveForChange, EMfields3D * field);
void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);

#endif
