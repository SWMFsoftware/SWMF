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

/***************************************************************************
  BcFields3D.h  -  Library to manage boundary conditions for fields 
  -------------------
begin                : Fri Jan 2009
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef BcFields_H
#define BcFields_H

#include "VCtopology3D.h"

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, double ***vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);

/** set the boundary condition on boundaries */
void BCface_P(int nx, int ny, int nz, double ***vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);

/** set the boundary condition on boundaries Particles*/
void BCface_P(int nx, int ny, int nz, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, const VirtualTopology3D * vct);

#endif
