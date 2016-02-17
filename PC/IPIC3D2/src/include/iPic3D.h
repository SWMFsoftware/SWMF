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
    void setSIDt(double SIDt);
    void SetCycle(int iCycle);
#endif

    
  };

}

#endif
