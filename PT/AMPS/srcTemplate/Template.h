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
#endif//_EXOSPHERE__ORBIT_CALCULATION__MODE_ == _PIC_MODE_ON_

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
      const  int    nSampleInterval = 100;
      const  double GammaMax = 1.0;
      const  double GammaMin = 0.0;
      // function implementations are in Template_Sampling_Gamma.cpp
      extern double SampleBuffer[nSampleInterval];
      void Init();
      void SampleData();
      // the following functions is called by the code's core;
      // prints data to a separate file (see implementation)
      void PrintSampleData(int DataOutputFileNumber);
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
      //......................................................................
      // local values of position, velocity and acceleration
      double x_LOCAL[3], v_LOCAL[3], accl_LOCAL[3]={0.0};
      memcpy(x_LOCAL,x,3*sizeof(double));
      memcpy(v_LOCAL,v,3*sizeof(double));

      //......................................................................
      // calculate Lorentz force if needed
#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
      // find electro-magnetic field
      double E[3],B[3];
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 
      // if no coupler is used -> use characteristic values
      memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
      memcpy(B,Exosphere::swB_Typical,3*sizeof(double));
#else 
      // if coupler is used -> get values from it
      //......................................................................
      // find the cell based on particles' position x_LOCAL and block startNode
      // input: x_LOCAL, startNode; output: nd, i,j,k
      long int nd;  // cell's number in the block
      int i,j,k;    // cell's coordinates in the block

      // fail-safe check: if the block doesn't exist => exit
      if (startNode->block==NULL) 
	exit(__LINE__,__FILE__,"Error: the block is not initialized");

      // flag: true - exit if a point is not found in the block / false: don't
      nd = PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false);

      // fail-safe check: if a point isn't found, try seacrhing in other blocks
      if (nd==-1) {
	// try to found the block the point is in;
	// starting point for search is block startNode
	startNode=PIC::Mesh::mesh.findTreeNode(x_LOCAL,startNode);
	nd=PIC::Mesh::mesh.fingCellIndex(x_LOCAL,i,j,k,startNode,false);
	// if still not found => exit
	if (nd==-1) 
	  exit(__LINE__,__FILE__,"Error: the cell is not found");
      }

      //......................................................................
      // finally, get fields' values at the cell
      PIC::CPLR::GetBackgroundElectricField(E,x_LOCAL,nd,startNode);
      PIC::CPLR::GetBackgroundMagneticField(B,x_LOCAL,nd,startNode);
#endif//_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__OFF_ 

      //......................................................................
      // calculate acceleraton due to Lorentz force
      double ElCharge=PIC::MolecularData::GetElectricCharge(spec);
      double mass=PIC::MolecularData::GetMass(spec);
      accl_LOCAL[0]+=ElCharge*(E[0]+v_LOCAL[1]*B[2]-v_LOCAL[2]*B[1])/mass;
      accl_LOCAL[1]+=ElCharge*(E[1]-v_LOCAL[0]*B[2]+v_LOCAL[2]*B[0])/mass;
      accl_LOCAL[2]+=ElCharge*(E[2]+v_LOCAL[0]*B[1]-v_LOCAL[1]*B[0])/mass;
#endif//_FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_

      //......................................................................
      // calculate gravitational force if needed
#if _FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_
      // find the distance to the massive body's center, here it's {0, 0, 0} 
      double r, r2;
      r2=x_LOCAL[0]*x_LOCAL[0]+x_LOCAL[1]*x_LOCAL[1]+x_LOCAL[2]*x_LOCAL[2];
      r=sqrt(r2);
      double MassTemplateBody = 0.0;
      // calculate acceleration due to gravity
      for (int idim=0;idim<DIM;idim++) {
	accl_LOCAL[idim]-=GravityConstant*MassTemplateBody/r2*x_LOCAL[idim]/r;
      }
#endif//_FORCE_GRAVITY_MODE_ == _PIC_MODE_ON_

      //......................................................................
      //copy the local values of the acceleration to the global ones
      memcpy(accl,accl_LOCAL,3*sizeof(double));
    }
  }
}

#endif//_TEMPLATE_H_
