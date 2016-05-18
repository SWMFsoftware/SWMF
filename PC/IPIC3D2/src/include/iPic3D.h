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
  iPIC3D.cpp  -  Main file for 3D simulation
  -------------------
 ************************************************************************** */

#ifndef _IPIC3D_H_
#define _IPIC3D_H_

class Timing;

#ifndef NO_MPI
#include "mpi.h"
#endif
#include "ipicfwd.h"
#include "assert.h"
#include <iostream>
#include <string>
using std::string;
using std::stringstream;
using std::cout;
using std::endl;
#ifndef NO_HDF5
class OutputWrapperFPP;
#endif
namespace iPic3D {

  class c_Solver {

  public:
    ~c_Solver();
    c_Solver():
    col(0),
      vct(0),
      grid(0),
      EMf(0),
      part(0),
#ifndef NO_HDF5
      outputWrapperFPP(0),
#endif
      Ke(0),
      BulkEnergy(0),
      momentum(0),
      Qremoved(0),
      my_clock(0)
	{}
    int Init(int argc, char **argv, double inittime = 0.0,
	     stringstream *param = NULL, int iIPIC = 0, int *paramint = NULL,
	     double *griddim = NULL, double *paramreal = NULL,
	     stringstream *ss = NULL, bool doCoupling=false);
    void CalculateMoments();
    void CalculateField(int cycle);
    bool ParticlesMover();
    void CalculateB();
    //
    // output methods
    //
    void WriteRestart(int cycle);
    void WriteConserved(int cycle);
    void WriteVelocityDistribution(int cycle);
    void WriteVirtualSatelliteTraces();
    void WriteFields(int cycle);
    void WriteParticles(int cycle);
    void WriteTestParticles(int cycle);
    void WriteOutput(int cycle);
    void Finalize();

    int FirstCycle() { return (first_cycle); }
    int get_myrank() { return (myrank); }
    int LastCycle();
    
  private:
    void pad_particle_capacities();
    void convertParticlesToSoA();
    void convertParticlesToAoS();
    void convertParticlesToSynched();
    void sortParticles();

  private:
    //static MPIdata * mpi;
    Collective    *col;
    VCtopology3D  *vct;
    Grid3DCU      *grid;
    EMfields3D    *EMf;
    Particles3D   *part;
    Particles3D   *testpart;
    double        *Ke;
    double        *BulkEnergy;
    double        *momentum;
    double        *Qremoved;
    Timing        *my_clock;

#ifndef NO_HDF5
    OutputWrapperFPP& fetch_outputWrapperFPP(){
      assert(outputWrapperFPP);
      return *outputWrapperFPP;}
    OutputWrapperFPP *outputWrapperFPP;
#endif

    //bool verbose;
    string SaveDirName;
    string RestartDirName;
    string cqsat;
    string cq;
    string ds;
    string num_proc_str;
    int restart_cycle;
    int restart_status;
    int first_cycle;
    int ns;
    int nstestpart;
    int nprocs;
    int myrank;
    int nsat;
    int nDistributionBins;
    double Eenergy;
    double Benergy;
    double TOTenergy;
    double TOTmomentum;

    //the below used for IO
    MPI_Request *headerReq;
    MPI_Request *dataReq;
    MPI_Request *footReq;
    float *testpclPos;
    int    pclbuffersize;
    float *testpclVel;
    MPI_File fh;
    MPI_Status*  status;
    float**** fieldwritebuffer;
    MPI_Request fieldreqArr[4];//E+B+Je+Ji
    MPI_File    fieldfhArr[4];
    MPI_Status  fieldstsArr[4];
    int fieldreqcounter;
    
    float*** momentwritebuffer;
    MPI_Request momentreqArr[14];//rho+PXX+PXY+PXZ++PYY+PYZ+PZZ for species0,1
    MPI_File    momentfhArr[14];
    MPI_Status  momentstsArr[14];
    int momentreqcounter;

#ifdef BATSRUS
  public:
    void SetParam(int *paramint, double *griddim, double *paramreal, stringstream *ss);
    void SyncWithFluid(int cycle);

    void GetNgridPnt(int *nPoint);
    void GetGridPnt(double *Pos_I);
    void setStateVar(double *State_I, int *iPoint_I);
    void getStateVar(int nDim, int nPoint, double *Xyz_I, double *data_I, int nVar);
    void findProcForPoint(int nPoint, double *Xyz_I, int *iProc_I);    
    double getDt();
    double getSItime();
    void setSIDt(double SIDt, bool isSWMFDt);
    double calSIDt();
    void updateSItime();
    void SetCycle(int iCycle);

    // IDL output. Begin-----------------
    static const int nDimMax=3;
    int nPlotFile;
    double **plotRange_ID;
    int **plotIndexRange_ID; // Local index range.
    string *nameSnapshot_I;
    static const int nVarMax = 50; 
    string **Var_II;
    int *nVar_I;
    long *nCell_I;
    bool IsBinary;
    int nByte;
    long nGlobalNode;
    string *outputFormat_I;
    string *outputUnit_I;
    string *plotVar_I;
    bool *doOutputParticles_I;
    int *iSpeciesOutput_I;
    double No2OutL, No2OutV, No2OutB, No2OutRho, No2OutP, No2OutJ;
    
    void write_plot_idl(int cycle);
    void write_plot_init();
    void write_plot_header(int iPlot, int cycle);
    void write_plot_data(int iPlot, int cycle);
    void set_output_unit(int iPlot);
    
    // IDL output. End-----------------
#endif

    
  };

}

#endif
