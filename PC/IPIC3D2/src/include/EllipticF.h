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
EllipticF.h  -  Elliptic Functions Tool
                            -------------------
developers: Found on the web
********************************************************************************************/

#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

double drf(double x, double y, double z, int *piErr);
double drd(double x, double y, double z, int *piErr);

#define EllipticF(k,ierr) drf(0.0,1.0-pow(k,2),1.0,&ierr)

#define EllipticE(k,ierr) (drf(0.0,1.0-pow(k,2),1.0,&ierr)-(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

#endif
