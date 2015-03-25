//$Id$

#ifndef _TEMPLATE_H_
#define _TEMPLATE_H_

// AMPS core code
#include "pic.h"

// Exosphere model
#include "Exosphere.h"

// self-explanatory
#include "constants.h"

/*****************************************************************************
 * Below are headers with physicals models to be used in simulation,
 * their location is src/models/<name-of-the-model>
 * UNCOMMENT required models below or add necessary ones
 *****************************************************************************/
// charge exchange physical model
// #include "ChargeExchange.h"
//----------------------------------
// surface sputtering model
// #include "Sputtering.h"
//----------------------------------
// electron impact ionization model
// #include "ElectronImpact.h"
//----------------------------------
// photolytic reactions model
// #include "PhotolyticReactions.h"

/*****************************************************************************
 * This preprocessor branch is added
 * in order for application to be able to work in conjunction with SPICE
 * or in stand-alone mode (for purpose of debugging, testing etc.)
 * may remove it if only one option is possible
 *****************************************************************************/
#if _EXOSPHERE__ORBIT_CALCULATION__MODE_ == _PIC_MODE_ON_
// SPICE IS used for this application: include SPICE header
#include "SpiceUsr.h"
#else
// SPICE is NOT used for this application: include empty substitutes
#include "SpiceEmptyDefinitions.h"
#endif

/*****************************************************************************
 * The following namespace contains variables, functions etc. 
 * which are specific to this particular application
 *****************************************************************************/
namespace Template {
  // it is based on a generic application named Exosphere
  using namespace Exosphere;
  //--------------------------------------------------------------------------
  // specific initialization procedures
  void Init_BeforeParser();
  //--------------------------------------------------------------------------
  // methods for sampling & printing physical parameters to a separate file
  namespace Sampling{
    // it is based on a generic Exosphere::Sampling
    using namespace Exosphere::Sampling;
    /*************************************************************************
     * User may define procedure for sampling any parameters;
     * below are declarations for sampling a generic parameter Gamma:
     * if needed - use as template, otherwise - remove / comment out
     *************************************************************************/
    namespace Gamma {
      const  int    nSampleIntervals = 100;
      const  double GammaMax = 1.0;
      const  double GammaMin = 0.0;
      // function implementations are in Template_Sampling_Gamma.cpp
      extern double SamplingBuffer[nSampleIntervals];
      void Init();
      void SampleData();
      // the following functions is called by the code's core;
      // prints data to a separate file (see implementation)
      void PrintSampledData(int DataOutputFileNumber);
    }
  }
  //--------------------------------------------------------------------------
  // methods for coupling with SWMF
  namespace Coupling {
    /*************************************************************************
     * Simulation can be performed within SWMF; if data is to be passed 
     * to framework from AMPS, a "Send" procedure is required
     *************************************************************************/
    // implementations are in Template_Coupling.cpp
    void Send(char *NameVar, int *nVarIn, int *nDimIn, 
	      int *nPoint, double *Xyz_DI, double *Data_VI);
  }
  //--------------------------------------------------------------------------
  // methods for printing physical parameters to the general AMPS output file
  namespace Output{
    /*************************************************************************
     * Additional data may be printed to the general output file
     * as well as carried by particles;
     * below are declarations for functions needed for the core code
     * to include a generic parameter Alpha in the printing queue;
     * if needed - use as template, otherwise - remove / comment out
     *************************************************************************/
    extern int TotalDataLength;
    extern int AlphaOffset; 
    void Init();
    // request for space in the particles' data buffer
    int RequestDataBuffer(int offset);
    // printing of variables' names
    void PrintVariableList(FILE* fout,int DataSetNumber);
    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
    // printing of variables' values
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
  }
  //--------------------------------------------------------------------------
  // methods for the reaction processing specific for the application
  namespace Physics {
    /*************************************************************************
     * Physical models' procedures are included in the beginning of the header;
     * below is the physics processor which uses procedures from these models;
     * physical processors are either continuous with characteristic life-time
     * between events, e.g. interaction with background plasma etc.,
     * or incident-based, e.g intercation with surface of planet etc.
     *************************************************************************/
    // Particle reactions, e.g. photoionization, electron impact etc.
    namespace Continuous {
      // effective net life-time    
      double LifeTime(double *x, int spec, long int ptr,
		      bool &ReactionAllowedFlag,
		      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      
      //----------------------------------------------------------------------
      // Processor for continuous physical processor:
      int ReactionProcessor(double *xInit,double *xFinal,double *vFinal,
			    long int ptr,int &spec,
			    PIC::ParticleBuffer::byte *ParticleData, 
			    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
    }
    //------------------------------------------------------------------------
    // Processes at key events, e.g. particle hits a surface of a planet
    namespace Incidential {
      // processes when particle hits a surface
      int ParticleSphereInteraction(int spec, long int ptr, 
				    double *x, double *v,
				    double &dtTotal,
				    void *NodeDataPonter,
				    void *SphereDataPointer);      
    }
    //------------------------------------------------------------------------
    // self-explanatory
    void inline TotalParticleAcceleration(double *accl, int spec, long int ptr,
					  double *x, double *v, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {

      accl[0]=0.0; accl[1]=0.0;  accl[2]=0.0; 

    }
  }
}

#endif
