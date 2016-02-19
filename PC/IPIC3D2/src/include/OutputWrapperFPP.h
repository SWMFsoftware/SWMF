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

#ifndef OutputWrapperFPP_h
#define OutputWrapperFPP_h
// ===
// OutputWrapperFPP: output wrapper for file-per-process output
//
//   This class should provide a mechanism to avoid having
//   repeatedly opening and closing the same file.
// ===
#include "ipicfwd.h"
#include "PSKOutput.h"
#ifndef NO_HDF5
#include "PSKhdf5adaptor.h"
#endif

using namespace PSK;

class OutputWrapperFPP
{
 private:
  #ifndef NO_HDF5
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output
  #endif // NO_HDF5
  int cartesian_rank;
  string SaveDirName;
  string RestartDirName;
  string output_file;
  string restart_file;
 public:
  void init_output_files(
    Collective    *col,
    VCtopology3D  *vct,
    Grid3DCU      *grid,
    EMfields3D    *EMf,
    Particles3D   *part,
    int 		  ns,
    Particles3D   *testpart,
    int 		  nstestpart);
  void append_output(const char* tag, int cycle);
  void append_output(const char* tag, int cycle, int sample);
  void append_restart(int cycle);

#ifdef BATSRUS
  void getFluidState(int nDim, int nPoint, double *Xyz_I, double *data_I,
		     int nVar);
#endif
  
};

#endif // OutputWrapperFPP_h
