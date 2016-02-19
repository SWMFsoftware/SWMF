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

#include "BcFields3D.h"

/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, double ***vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  // XLEFT
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = -vector[1][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = vector[1][i][j];
          }
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = -vector[nx - 2][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = vector[nx - 2][i][j];
          }
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = -vector[i][1][j];

          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = vector[i][1][j];
          }
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = -vector[i][ny - 2][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = vector[i][ny - 2][j];
          }
        break;

    }


  }
  // Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = -vector[i][j][1];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = 0.0;

          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = vector[i][j][1];
          }
        break;
    }
  }
  // Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = -vector[i][j][nz - 2];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = vector[i][j][nz - 2];
          }
        break;
    }
  }

}


// / particles
/** set the boundary condition on boundaries */
void BCface_P(int nx, int ny, int nz, double ***vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  // XLEFT
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = -vector[1][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[0][i][j] = vector[1][i][j];
          }
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = -vector[nx - 2][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[nx - 1][i][j] = vector[nx - 2][i][j];
          }
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = -vector[i][1][j];

          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][0][j] = vector[i][1][j];
          }
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = -vector[i][ny - 2][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[i][ny - 1][j] = vector[i][ny - 2][j];
          }
        break;

    }


  }
  // Zleft
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = -vector[i][j][1];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = 0.0;

          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][0] = vector[i][j][1];
          }
        break;
    }
  }
  // Zright
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = -vector[i][j][nz - 2];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[i][j][nz - 1] = vector[i][j][nz - 2];
          }
        break;
    }
  }

}

// SPECIES
/** set the boundary condition on boundaries */
void BCface(int nx, int ny, int nz, int ns, double ****vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  // XLEFT
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = -vector[ns][1][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = vector[ns][1][i][j];
          }
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = -vector[ns][nx - 2][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = vector[ns][nx - 2][i][j];
          }
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = -vector[ns][i][1][j];

          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = vector[ns][i][1][j];
          }
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = -vector[ns][i][ny - 2][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = vector[ns][i][ny - 2][j];
          }
        break;

    }


  }
  // Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = -vector[ns][i][j][1];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = 0.0;

          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = vector[ns][i][j][1];
          }
        break;
    }
  }
  // Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = -vector[ns][i][j][nz - 2];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = vector[ns][i][j][nz - 2];
          }
        break;
    }
  }

}

// SPECIES
/** set the boundary condition on boundaries Particles*/
void BCface_P(int nx, int ny, int nz, int ns, double ****vector,
  int bcFaceXright, int bcFaceXleft,
  int bcFaceYright, int bcFaceYleft,
  int bcFaceZright, int bcFaceZleft,
  const VirtualTopology3D * vct)
{
  // XLEFT
  if (vct->getXleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceXleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = -vector[ns][1][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][0][i][j] = vector[ns][1][i][j];
          }
        break;


    }

  }
  // XRIGHT
  if (vct->getXright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceXright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = -vector[ns][nx - 2][i][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < ny; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][nx - 1][i][j] = vector[ns][nx - 2][i][j];
          }
        break;

    }

  }
  // YLEFT
  if (vct->getYleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceYleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = -vector[ns][i][1][j];

          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][0][j] = vector[ns][i][1][j];
          }
        break;

    }
  }
  // YRIGHT
  if (vct->getYright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceYright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = -vector[ns][i][ny - 2][j];
          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < nz; j++) {
            vector[ns][i][ny - 1][j] = vector[ns][i][ny - 2][j];
          }
        break;

    }


  }
  // Zleft
  if (vct->getZleft_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceZleft) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = -vector[ns][i][j][1];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = 0.0;

          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][0] = vector[ns][i][j][1];
          }
        break;
    }
  }
  // Zright
  if (vct->getZright_neighbor_P() == MPI_PROC_NULL) {
    switch (bcFaceZright) {

      case 0:                  // Dirichilet = 0 Second Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = -vector[ns][i][j][nz - 2];


          }
        break;

      case 1:                  // Dirichilet = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = 0.0;
          }
        break;

      case 2:                  // Neumann = 0 First Order
        for (int i = 0; i < nx; i++)
          for (int j = 0; j < ny; j++) {
            vector[ns][i][j][nz - 1] = vector[ns][i][j][nz - 2];
          }
        break;
    }
  }

}
