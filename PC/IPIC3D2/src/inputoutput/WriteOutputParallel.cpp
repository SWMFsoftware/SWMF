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
#include "WriteOutputParallel.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "VCtopology3D.h"
#include "errors.h"
#include <string>
#include <sstream>
#include <iomanip>
using std::string;

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef PHDF5
  string       grpname;
  string       dtaname;

  stringstream filenmbr;
  string       filename;

  /* ------------------- */
  /* Setup the file name */
  /* ------------------- */

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  /* ---------------------------------------------------------------------------- */
  /* Define the number of cells in the globa and local mesh and set the mesh size */
  /* ---------------------------------------------------------------------------- */

  int nxc = grid->getNXC();
  int nyc = grid->getNYC();
  int nzc = grid->getNZC();

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { nxc-2,            nyc-2,            nzc-2 };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  /* --------------------------------------- */
  /* Declare and open the parallel HDF5 file */
  /* --------------------------------------- */

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getComm());

  outputfile.CreatePHDF5file(L, dglob, dlocl, false);

  /* ------------------------ */
  /* Write the Electric field */
  /* ------------------------ */

  outputfile.WritePHDF5dataset("Fields", "Ex", EMf->getExc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ey", EMf->getEyc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Ez", EMf->getEzc(), nxc-2, nyc-2, nzc-2);

  /* ------------------------ */
  /* Write the Magnetic field */
  /* ------------------------ */

  outputfile.WritePHDF5dataset("Fields", "Bx", EMf->getBxc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "By", EMf->getByc(), nxc-2, nyc-2, nzc-2);
  outputfile.WritePHDF5dataset("Fields", "Bz", EMf->getBzc(), nxc-2, nyc-2, nzc-2);

  /* ----------------------------------------------- */
  /* Write the moments for each species */
  /* ----------------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream snmbr;
    snmbr << is;
    const string num = snmbr.str();

    // Charge Density
    outputfile.WritePHDF5dataset("Fields", string("Rho_")+num , EMf->getRHOcs(is), nxc-2, nyc-2, nzc-2);
    // Current
    outputfile.WritePHDF5dataset("Fields", string("Jx_")+num, EMf->getJxsc(is), nxc-2, nyc-2, nzc-2);
    outputfile.WritePHDF5dataset("Fields", string("Jy_")+num, EMf->getJysc(is), nxc-2, nyc-2, nzc-2);
    outputfile.WritePHDF5dataset("Fields", string("Jz_")+num, EMf->getJzsc(is), nxc-2, nyc-2, nzc-2);
  }

  outputfile.ClosePHDF5file();

#else  
  eprintf(
    " The input file requests the use of the Parallel HDF5 functions,\n"
    " but the code has been compiled using the sequential HDF5 library.\n"
    " Recompile the code using the parallel HDF5 options\n"
    " or change the input file options. ");
#endif

}
