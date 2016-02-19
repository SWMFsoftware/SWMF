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

/*******************************************************************************************
  CG.h  -  Conjugate Gradient Solver
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#ifndef CG_H
#define CG_H

#include "ipicfwd.h"

// These declarations are currently needed because Field is not anymore a class
// CG needs a pointer to the function that solves the fields.
// 
// To avoid changing all the code we typedef Field as of type EMfields3D (which is
// not derived anymore from Field). This will be improved in future releases.

class EMfields3D;
typedef EMfields3D Field;
typedef void (Field::*FIELD_IMAGE) (double *, double *, bool);
typedef void (*GENERIC_IMAGE) (double *, double *);

bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, FIELD_IMAGE FunctionImage, bool doSolveForChange, Field * field);

#endif
