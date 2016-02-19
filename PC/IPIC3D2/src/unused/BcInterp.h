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
  BcInterp.h  -  Library to manage boundary conditions for fields 
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef BcInterp_H
#define BcInterp_H
/** set the boundary condition on boundaries */
inline void BCfaceInterp(int xStartActiveN, int yStartActiveN, int zStartActiveN, int xEndActiveN, int yEndActiveN, int zEndActiveN, double ***vector, double value, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, double dx, double dy, double dz, VirtualTopology3D * vct) {
  // XLEFT
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = yStartActiveN; i <= yEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[xStartActiveN][i][j] *= value;

        break;

      case 1:
        break;
      case 2:
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // multiply for value on boundaries
        for (int i = yStartActiveN; i <= yEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[xEndActiveN][i][j] *= value;
        break;
      case 1:
        break;
      case 2:
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[i][yStartActiveN][j] *= value;
        break;

      case 1:
        break;

      case 2:
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[i][yEndActiveN][j] *= value;
        break;

      case 1:
        break;

      case 2:
        break;

    }


  }
  // ZRIGHT
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = yStartActiveN; j <= yEndActiveN; j++)
            vector[i][j][zStartActiveN] *= value;
        break;

      case 1:
        break;

      case 2:
        break;
    }
  }
  // ZLEFT
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = yStartActiveN; j <= yEndActiveN; j++)
            vector[i][j][zEndActiveN] *= value;
        break;

      case 1:
        break;

      case 2:
        break;
    }
  }

}
/** set the boundary condition on boundaries */
inline void BCfaceInterp(int xStartActiveN, int yStartActiveN, int zStartActiveN, int xEndActiveN, int yEndActiveN, int zEndActiveN, int ns, double ****vector, double value, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, double dx, double dy, double dz, VirtualTopology3D * vct) {
  // XLEFT
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = yStartActiveN; i <= yEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[xStartActiveN][i][j][ns] *= value;
        break;

      case 1:
        break;
      case 2:
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // multiply for value on boundaries
        for (int i = yStartActiveN; i <= yEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[xEndActiveN][i][j][ns] *= value;
        break;
      case 1:
        break;
      case 2:
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[i][yStartActiveN][j][ns] *= value;
        break;

      case 1:
        break;

      case 2:
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = zStartActiveN; j <= zEndActiveN; j++)
            vector[i][yEndActiveN][j][ns] *= value;
        break;

      case 1:
        break;

      case 2:
        break;

    }


  }
  // ZRIGHT
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = yStartActiveN; j <= yEndActiveN; j++)
            vector[i][j][zStartActiveN][ns] *= value;
        break;

      case 1:
        break;

      case 2:
        break;
    }
  }
  // ZLEFT
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // multiply for value on boundaries
        for (int i = xStartActiveN; i <= xEndActiveN; i++)
          for (int j = yStartActiveN; j <= yEndActiveN; j++)
            vector[i][j][zEndActiveN][ns] *= value;
        break;

      case 1:
        break;

      case 2:
        break;
    }
  }


}
#endif
