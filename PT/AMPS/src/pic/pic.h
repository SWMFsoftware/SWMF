//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//===============================================
//$Id$
//===============================================
//the header for the new particle model

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>



#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

#ifndef _PIC_
#define _PIC_

//the global model settings
#include "picGlobal.dfn"

//load the macros that defined symbolic references to the species
#include "picSpeciesMacro.dfn"

//load the user defined settings
/*
#if _PIC_USER_DEFINITION_MODE_ == _PIC_USER_DEFINITION_MODE__ENABLED_
#include "UserDefinition.PIC.h"
#endif
*/



#if _PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE_ == _PIC__USER_DEFINED__LOAD_SPECIES_MACRO__MODE__ON_
$MARKER:SPECIES-MACRO-DEFINIETION-USED-IN-SIMULATION$
#endif

#include "ifileopr.h"
#include "specfunc.h"
#include "constants.h"
#include "rnd.h"

//include the appropriate mesh header
#if DIM == 3
#include "meshAMR3d.h"
#elif DIM == 2
#include "meshAMR2d.h"
#else 
#include "meshAMR1d.h"
#endif

#include "meshAMRinternalSurface.h"




/*
#define Kbol 1.3806503E-23
#define ElectronCharge 1.602176565E-19
#define ElectronMass 9.10938291E-31
#define ProtonMass 1.67262158E-27
#define eV2J ElectronCharge
*/



namespace PIC {

  //Global constants of the PIC solver
  //extern int nTotalSpecies;
  static const int nTotalSpecies=1;

  //The currect and total number of processors used in the simulation
  extern int ThisThread,nTotalThreads;

  //the output and prefix for the diagnostic information
  extern FILE* DiagnospticMessageStream;
  extern char DiagnospticMessageStreamName[_MAX_STRING_LENGTH_PIC_];

  //the directory for output files
  extern char OutputDataFileDirectory[_MAX_STRING_LENGTH_PIC_];


  //the variables that controls the sampling procedure
  //LastSampleLength: the length of the last collected sample
  //CollectingSampleCounter: the number of iterations passed in the current sample loop
  //RequireSampleLength: the maximum number of iterations required in the current sample loop
  //DataOutputFileNumber: the number of the NEXT output file
  //SamplingMode: if _RESTART_SAMPLING_MODE_ -> after each output of the data file, the sampling buffer is flushed, _ACCUMULATE_SAMPLING_MODE_ -> the sampling data are saved and used for the next output of the flow file

  #define _RESTART_SAMPLING_MODE_    0
  #define _ACCUMULATE_SAMPLING_MODE_ 1
  extern long int LastSampleLength,CollectingSampleCounter,RequiredSampleLength,DataOutputFileNumber;
  extern int SamplingMode;

  //the tags for the data exchenge between processors
  #define _PIC_SUBDOMAIN_BOUNDARY_LAYER_SAMPLING_DATA_EXCHANGE_TAG_   0
  #define _PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_     1

  //handle run time signals and exeptions
  void SignalHandler(int);

  //perform one time step
  void TimeStep();
//  void Sampling();

  //init the particle solver
  void InitMPI();
  void Init_BeforeParser();
  void Init_AfterParser();

  //the list of functions used to exchenge the execution statiscics
  typedef void (*fExchangeExecutionStatistics) (CMPI_channel*,long int);
  extern vector<fExchangeExecutionStatistics> ExchangeExecutionStatisticsFunctions;



  //structural elements of the solver
  namespace Parser {
    void Run(char*);
    void readMain(CiFileOperations&);
    void readGeneral(CiFileOperations&);

  }

  //ray tracing and calculation of the shadow regions on the NASTRAN surfaces
  namespace RayTracing {
    extern unsigned int nCallsTestDirectAccess;

    bool GetBlockExitPoint(double *xBlockMin,double *xBlockMax,double *x0Ray,double *lRay,double *xBlockExit, double *xFaceExitLocal, int &nExitFace);
    bool TestDirectAccess(double *xStart,double *xTarget);

    void SetCutCellShadowAttribute(double *xLightSource, bool ParallelExecution=false);
  }

  //define the test-run parameters
  namespace ModelTestRun {
    extern bool mode;
    extern int nTotalIteraction;
  }

  //run tame calculation of the check sums
  namespace RunTimeSystemState {
    void GetParticleFieldCheckSum(char *msg=NULL);
    void GetParticleFieldCheckSum(long int nline,char *fname);

    void GetDomainDecompositionCheckSum(char *msg=NULL);
    void GetDomainDecompositionCheckSum(long int nline,char *fname);

    void GetMeanParticleMicroscopicParameters(FILE* fout,char *msg=NULL);
    void GetMeanParticleMicroscopicParameters(FILE* fout,long int nline,char *fname);
    void GetMeanParticleMicroscopicParameters(const char *fname);
  }


  namespace MolecularData {

    //Available molecular models
    #define _HS_MOL_MODEL_   0
    #define _VHS_MOL_MODEL_  1
    #define _VSS_MOL_MODEL_  2
    extern int MolModelCode;
    void SetMolType(int);
    int GetMolType();

    //modeling of external species
    #define _EXTERNAL_SPECIES_ON_  true
    #define _EXTERNAL_SPECIES_OFF_ false
    extern bool ExternalSpeciesModelingFlag;

    //modeling of unimolecular reactions
    #define _UNIMOLECULAR_REACTIONS_ON_   true
    #define _UNIMOLECULAR_REACTIONS_OFF_  false
    extern bool UnimolecularReactionFlag;

    //modeling of internal degrees of freedom
    #define _INTERNAL_DEGREES_OF_FREEDOM_ON_  true
    #define _INTERNAL_DEGRESS_OF_FREEDOM_OFF_ false
    extern bool InternalDegreesOfFreedomModelingFlag;

    //molecular models
    namespace MolecularModels {
      namespace HS {
        //the table of the constant collsion cross sections and reference diameter
        static const double ConstantCollisionCrossSectionTable[1][1]={{0.0}};
        static const double ConstantReferenceDiameter[1][1]={{0.0}};


        inline double GetTotalCrossSection(int s0,int s1) {return ConstantCollisionCrossSectionTable[s0][s1];}
        inline double GetDiam(int s0,int s1) {return ConstantReferenceDiameter[s0][s1];}
        inline double GetRefDiam(int s0,int s1) {return ConstantReferenceDiameter[s0][s1];}

      }

      namespace VHS {
      using namespace HS;

      }

      namespace VSS {
      using namespace VHS;

      }
    }

    inline double GetRefDiam(int s0,int s1) {
#if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__HS_
      return MolecularModels::HS::GetRefDiam(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }

    inline double GetDiam(int s0,int s1) {
#if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__HS_
      return MolecularModels::HS::GetRefDiam(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }

    inline double GetTotalCrossSect(double Vrel,int ptr0,int s0,int ptr1,int s1,long int ncell) {
#if _PIC__PARTICLE_COLLISION_MODEL_ == _PIC__PARTICLE_COLLISION_MODEL__HS_
      return MolecularModels::HS::GetTotalCrossSection(s0,s1);
#else
      exit(__LINE__,__FILE__,"not implemented");
      return 1.0;
#endif
    }


    //init the molecular data buffers
    void Init();

    //mass of particles
    static const double MolMass[]={0.0};
    static const double ElectricChargeTable[]={0.0};

//    extern double *MolMass;
//    void SetMass(double,int);
    inline double GetMass(int spec) {return MolMass[spec];}
    inline double GetElectricCharge(int spec) {return ElectricChargeTable[spec];}


    //get and set the value of the electric charge for a species
/*    extern double *ElectricCharge;
    int GetElectricCharge(int);
    void SetElectricCharge(double,int);*/

    //get and set the species numbers and chemical symbols
    static const char ChemTable[][_MAX_STRING_LENGTH_PIC_]={"nothin is defined"};


//    extern char **ChemTable; // <- The table of chemical symbols used in the simulation
    extern char **LoadingSpeciesList; // <- the list of species that CAN BE locased in the simualtion


    int GetSpecieNumber(char*);
    void SetChemSymbol(char*,int);
    void GetChemSymbol(char*,int);
    const char* GetChemSymbol(int);

    //set and get the specie type (gas, external, background)
    static const int SpcecieTypeTable[]={-1};
//    extern int *SpcecieTypeTable;
//    void SetSpecieType(int SpcecieType,int spec);
    int GetSpecieType(int spec);


    namespace Parser {
      void InitChemTable(CiFileOperations&);
      int GetSpeciesNumber(int&,CiFileOperations&);
      void SpeciesBlock(char*,int,CiFileOperations&);
      void run(CiFileOperations&);
    }

  }

  namespace ParticleBuffer {
    typedef unsigned char byte;

    //the set of parameters that determine the 'double' parameters of a particle
    #define _PIC_PARTICLE_DATA_SPECIEDID_OFFSET_ 0
    #define _PIC_PARTICLE_DATA_NEXT_OFFSET_ (_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_ + sizeof(unsigned char))
    #define _PIC_PARTICLE_DATA_PREV_OFFSET_ (_PIC_PARTICLE_DATA_NEXT_OFFSET_ + sizeof(long int))
    #define _PIC_PARTICLE_DATA_VELOCITY_OFFSET_  (_PIC_PARTICLE_DATA_PREV_OFFSET_ + sizeof(long int))
    #define _PIC_PARTICLE_DATA_POSITION_OFFSET_  (_PIC_PARTICLE_DATA_VELOCITY_OFFSET_+ 3*sizeof(double))

    #define _PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_ (_PIC_PARTICLE_DATA_POSITION_OFFSET_+DIM*sizeof(double))

    //the individual particle's weight corection
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
    #define _PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_ _PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_
    #define _PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_  sizeof(double)
#elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
    #define _PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_ -1
    #define _PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_  0
#else
   exit(__LINE__,__FILE__,"Error: unknown option");
#endif


    //electrically charged dust
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_    (_PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_)
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_  (_PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_+sizeof(double))

#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_  (_PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_+2*sizeof(double))
   #define _PIC_PARTICLE_DATA__DUST_GRAIN__DATA_LENGTH_   3*sizeof(double)
#else
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_  -1
   #define _PIC_PARTICLE_DATA__DUST_GRAIN__DATA_LENGTH_    2*sizeof(double)
#endif

#else
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_MASS_OFFSET_    -1
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_RADIUS_OFFSET_  -1
   #define _PIC_PARTICLE_DATA__DUST_GRAIN_CHARGE_OFFSET_  -1
   #define _PIC_PARTICLE_DATA__DUST_GRAIN__DATA_LENGTH_    0
#endif


    //the total length of the default particle data
    #define _PIC_PARTICLE_DATA__DEFAULT_DATA_LENGTH_  (_PIC_PARTICLE_DATA__BISIC_DATA_LENGTH_+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET__DATA_LENGTH_+_PIC_PARTICLE_DATA__DUST_GRAIN__DATA_LENGTH_)


    //the total length of a data allocated for a particle
    extern long int ParticleDataLength;

    //The particle buffer's internal data
    extern byte *ParticleDataBuffer;
    extern long int MaxNPart,NAllPart,FirstPBufferParticle;

    //Request additional data for a particle
    void RequestDataStorage(long int &offset,int TotalDataLength);

    //the basic data access functions for a particle
    byte *GetParticleDataPointer(long int);

    //check the total particles number
    void CheckParticleList();

    /*
    double *GetX(long int);
    void GetX(double*,long int);
    void SetX(double*,long int);

    double *GetX(byte*);
    void GetX(double*,byte*);
    void SetX(double*,byte*);

    double *GetV(long int);
    void GetV(double*,long int);
    void SetV(double*,long int);

    double *GetV(byte*);
    void GetV(double*,byte*);
    void SetV(double*,byte*);

    unsigned int GetI(byte*);
    void SetI(unsigned int,byte*);

    unsigned int GetI(long int);
    void SetI(unsigned int,long int);

    long int GetPrev(long int);
    long int GetNext(long int);
    void SetPrev(long int,long int);
    void SetNext(long int,long int);

    long int GetPrev(byte*);
    long int GetNext(byte*);
    void SetPrev(long int,byte*);
    void SetNext(long int,byte*);

    double GetIndividualStatWeightCorrection(long int);
    double GetIndividualStatWeightCorrection(byte*);
    void SetIndividualStatWeightCorrection(double,long int);
    void SetIndividualStatWeightCorrection(double,byte*);
    */

    //==========================================================
    //get the idividual particle weight correction
    inline double GetIndividualStatWeightCorrection(long int ptr) {
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      return *((double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_));
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      return 1;
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }

    inline double GetIndividualStatWeightCorrection(byte *ParticleDataStart) {
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      return *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_));
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      return 1;
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }

    inline void SetIndividualStatWeightCorrection(double WeightCorrectionFactor,long int ptr) {
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(WeightCorrectionFactor)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      //do nothing
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }

    inline void SetIndividualStatWeightCorrection(double WeightCorrectionFactor,byte *ParticleDataStart) {
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ == _PIC_DEBUGGER_MODE_ON_
      if (!isfinite(WeightCorrectionFactor)) exit(__LINE__,__FILE__,"Error: Floating Point Exeption");
#endif
#endif

    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      exit(__LINE__,__FILE__,"Error: SetIndividualStatWeightCorrection cannot be used with _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_");
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }
    //==========================================================
    //get the particle position
    inline double *GetX(long int ptr) {
      return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
    }

    inline double *GetX(byte *ParticleDataStart) {
      return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
    }

    inline void GetX(double* x,long int ptr) {
      /*
      register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
      */

      memcpy(x,ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_,DIM*sizeof(double));
    }

    inline void GetX(double* x,byte *ParticleDataStart) {
      /*
      register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
      */

      memcpy(x,ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_,DIM*sizeof(double));
    }

    inline void SetX(double* x,long int ptr) {
//      register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
//      register int idim;

//=====================  DEBUG ================
/*

      if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+10<pow(2439.0e3,2)) {
        double r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

        cout << pow(2439.0e3,2)- r*r << "   "  << 2439.0e3-r << __FILE__ << "@" << __LINE__ << endl;

      }
*/
//=================   END DEBUG =============

//      for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];

      memcpy(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_,x,DIM*sizeof(double));
    }

    inline void SetX(double* x,byte *ParticleDataStart) {
//      register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
//      register int idim;

      //=====================  DEBUG ================

/*
            if ((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])+10<pow(2439.0e3,2)) {
              double r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

              cout << pow(2439.0e3,2)- r*r << "   "  << 2439.0e3-r   << __FILE__ << "@" << __LINE__ << endl;

            }
*/
      //=================   END DEBUG =============

//      for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];

      memcpy(ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_,x,DIM*sizeof(double));
    }

    //==========================================================
    //get the particle velocity
    inline double *GetV(long int ptr) {
      return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
    }

    inline double *GetV(byte *ParticleDataStart) {
      return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
    }

    inline void GetV(double* v,long int ptr) {
      /*
      register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
      */

      memcpy(v,ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,3*sizeof(double));
    }

    inline void GetV(double* v,byte *ParticleDataStart) {
      /*
      register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
      */

      memcpy(v,ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,3*sizeof(double));
    }

    inline void SetV(double* v,long int ptr) {
      /*
      register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
      */

/*      if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>1.0e9) {
        exit(__LINE__,__FILE__,"the velocity is too large");
      }*/

      memcpy(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,v,3*sizeof(double));
    }

    inline void SetV(double* v,byte *ParticleDataStart) {
      /*
      register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
      */

/*      if (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>1.0e9) {
        exit(__LINE__,__FILE__,"the velocity is too large");
      }*/

      memcpy(ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,v,3*sizeof(double));
    }


    //==========================================================
    //get the particle's species ID
    //the first 7 bits will be used for specie ID, the last 8th bit will be used to control wether the particle is allocated or not

    inline unsigned int GetI(byte* ParticleDataStart) {
      return ((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
    }

    inline unsigned int GetI(long int ptr) {
      return ((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
    }

    inline void SetI(int spec,byte* ParticleDataStart) {
      unsigned char flag,t=spec;

      flag=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x80);
      t|=flag;

      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    inline void SetI(int spec,long int ptr) {
      unsigned char flag,t=spec;

      flag=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x80);
      t|=flag;

      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    inline bool IsParticleAllocated(byte* ParticleDataStart) {
      unsigned char flag;

      flag=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x80);
      return (flag==0) ? false : true;
    }

    inline bool IsParticleAllocated(long int ptr) {
      unsigned char flag;

      flag=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x80);
      return (flag==0) ? false : true;
    }

    inline void SetParticleDeleted(byte* ParticleDataStart) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    inline void SetParticleDeleted(long int ptr) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    inline void SetParticleAllocated(byte* ParticleDataStart) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
      t|=0x80;

      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    inline void SetParticleAllocated(long int ptr) {
      unsigned char t;

      t=((*((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))) & 0x7f);
      t|=0x80;

      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=t;
    }

    //==========================================================
    //get/set prev
    inline long int GetPrev(long int ptr) {
      return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_PREV_OFFSET_));
    }

    inline long int GetPrev(byte* ParticleDataStart) {
      return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_PREV_OFFSET_));
    }

    inline void SetPrev(long int prev,long int ptr) {
      *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_PREV_OFFSET_))=prev;
    }

    inline void SetPrev(long int prev,byte* ParticleDataStart) {
      *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_PREV_OFFSET_))=prev;
    }

    //get/set next
    inline long int GetNext(long int ptr) {
      return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_NEXT_OFFSET_));
    }

    inline long int GetNext(byte* ParticleDataStart) {
      return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_NEXT_OFFSET_));
    }

    inline void SetNext(long int next,long int ptr) {
      *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_NEXT_OFFSET_))=next;
    }

    inline void SetNext(long int next,byte* ParticleDataStart) {
      *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_NEXT_OFFSET_))=next;
    }
    //========================================================

    //the particle buffer procedure
    void Init(long int);
    long int GetMaxNPart();
    long int GetAllPartNum();
    long int GetParticleDataLength();

    long int GetNewParticle();
    long int GetNewParticle(long int&);

    void DeleteParticle(long int);
    void DeleteParticle(long int,long int&);

    void CloneParticle(long int,long int);

    void ExcludeParticleFromList(long int,long int&);

    void SaveImageFile(int);
    void LoadImageFile(int);

    void PackParticleData(char*,long int);
    void UnPackParticleData(char*,long int);

    unsigned long GetChecksum();
  }



  namespace Mesh {
    class cDataCenterNode;

	  //the limiting size of the domain and the function controlling the local mesh resolution
	  extern double xmin[3],xmax[3];

	  typedef double (*fLocalMeshResolution) (double*);
	  extern fLocalMeshResolution LocalMeshResolution;

	  //the offset of the sampled infomation that is stored in 'center nodes'
    extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;

	  //the data and order that the data are saved in the associated data buffer of 'center nodes'
    //3. The offset of the data buffer for 'completed sample'
	  //4. The offset of the data buffer for 'collecting sample'

	  //sampling data each specie
	  //a. for each specie
	    //1. particle weight
	    //2. particle number
	    //3. particle velocity[3]
	    //4. particle pow(velocity[3],2)
	  //b. sampling data requested for involved physical models and external species
    extern int sampledParticleWeghtRelativeOffset,sampledParticleNumberRelativeOffset,sampledParticleNumberDensityRelativeOffset;
    extern int sampledParticleVelocityRelativeOffset,sampledParticleVelocity2RelativeOffset,sampledParticleSpeedRelativeOffset;
    extern int sampledParticleNormalParallelVelocityRelativeOffset,sampledParticleNormalParallelVelocity2RelativeOffset;
    extern int sampledExternalDataRelativeOffset;
    extern int sampleSetDataLength;



    //user defiend functions for printing the 'center node' data into an output file
    typedef void (*fPrintVariableListCenterNode)(FILE* fout,int DataSetNumber);
    typedef void (*fPrintDataCenterNode)(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,cDataCenterNode *CenterNode);
    typedef void (*fInterpolateCenterNode)(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,cDataCenterNode *CenterNode);

    extern vector<fPrintVariableListCenterNode> PrintVariableListCenterNode;
    extern vector<fPrintDataCenterNode> PrintDataCenterNode;
    extern vector<fInterpolateCenterNode> InterpolateCenterNode;

    //the class defining the 'central node' that contains the sampling data
    class cDataCenterNode : public cBasicCenterNode {
    public:
	    //parameters that defines the parameters of the associated data used for sampling and code running
      static int totalAssociatedDataLength,LocalParticleVolumeInjectionRateOffset;

//      long int FirstCellParticle,tempParticleMovingList;

	    char *associatedDataPointer;

	    inline int AssociatedDataLength() {
              return totalAssociatedDataLength;
            }

	    void SetAssociatedDataBufferPointer(char* ptr) {
              associatedDataPointer=ptr;
            }

	    inline char* GetAssociatedDataBufferPointer() {
              return associatedDataPointer;
            }

	    //clean the sampling buffers
	    void cleanDataBuffer() {
	      cBasicCenterNode::cleanDataBuffer();

//	      FirstCellParticle=-1,tempParticleMovingList=-1;

	      int i,length=totalAssociatedDataLength/sizeof(double);
	      double *ptr;
	      for (i=0,ptr=(double*)associatedDataPointer;i<length;i++,ptr++) *ptr=0.0;

	      if (totalAssociatedDataLength%sizeof(double)) exit(__LINE__,__FILE__,"Error: the cell internal buffers contains data different from double");
 	    }

	    //init the buffers
	    cDataCenterNode() : cBasicCenterNode() {
	      associatedDataPointer=NULL;
	    }

	    //get the sampled macroscopic parameter of the flow
	    double GetParticleNumber(int s) {
	      double res=0.0;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          res=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNumberRelativeOffset)))/PIC::LastSampleLength;
        }

        return res;
	    }

      double GetRealParticleNumber(int s) {
        double res=0.0;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          res=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)))/PIC::LastSampleLength;
        }

        return res;
      }

      double GetNumberDensity(int s) {
        double res=0.0;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          res=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNumberDensityRelativeOffset)))/PIC::LastSampleLength;
        }

        return res;
      }

      double GetCompleteSampleCellParticleWeight(int s) {
        return *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset));
      }

      void GetBulkVelocity(double *v,int s) {
        int idim;
        double TotalWeight,*SampledData;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) {
          exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        }
        #endif

        if (PIC::LastSampleLength!=0) {
          TotalWeight=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)));
          SampledData=3*s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocityRelativeOffset);

          if (TotalWeight>0.0) for (idim=0;idim<DIM;idim++) v[idim]=SampledData[idim]/TotalWeight; /// /PIC::LastSampleLength;
          else for (idim=0;idim<DIM;idim++) v[idim]=0.0;
        }
        else for (idim=0;idim<DIM;idim++) v[idim]=0.0;
      }

      void GetBulkVelocitySquared(double *v2,int s) {
        int idim;
        double TotalWeight,*SampledData;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          TotalWeight=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)));
          SampledData=3*s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocity2RelativeOffset);

          if (TotalWeight>0.0) for (idim=0;idim<DIM;idim++) v2[idim]=SampledData[idim]/TotalWeight; // /PIC::LastSampleLength;
          else for (idim=0;idim<DIM;idim++) v2[idim]=0.0;
        }
        else for (idim=0;idim<DIM;idim++) v2[idim]=0.0;
      }

      double GetMeanParticleSpeed(int s) {
        double TotalWeight,*SampledData,res=0.0;



//==============================  DEBUGGER ===============
//        return Measure;

//============================== END DEBUGGER ============


        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          TotalWeight=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)));
          SampledData=s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleSpeedRelativeOffset);

          if (TotalWeight>0.0) res=(*SampledData)/TotalWeight; // /PIC::LastSampleLength;
        }

        return res;
      }

      double GetTranslationalTemperature(int s) {
        int idim;
        double res=0.0,v[3]={0.0,0.0,0.0},v2[3]={0.0,0.0,0.0};

        GetBulkVelocity(v,s);
        GetBulkVelocitySquared(v2,s);

        for (idim=0;idim<3;idim++) res+=v2[idim]-v[idim]*v[idim];

        return PIC::MolecularData::GetMass(s)*res/(3.0*Kbol);
      }

      double GetParallelTranslationalTemperature(int s) {
#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
        return GetTranslationalTemperature(s);
#else
        double v=0.0,v2=0.0,w=0.0,res=0.0;

        w=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)));

        v=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNormalParallelVelocityRelativeOffset)));
        v2=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNormalParallelVelocity2RelativeOffset)));

        if (w>0.0) {
          v/=w,v2/=w;
          res=PIC::MolecularData::GetMass(s)*(v2-v*v)/Kbol;
        }

        return res;
#endif
      }

      double GetTangentialTranslationalTemperature(int s) {
#if _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE_ == _PIC_SAMPLE__PARALLEL_TANGENTIAL_TEMPERATURE__MODE__OFF_
        return GetTranslationalTemperature(s);
#else
      return 0.5*(3.0*GetTranslationalTemperature(s)-GetParallelTranslationalTemperature(s));
#endif
      }

	    //print the sampled data into a file
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread);
      /*
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread) {
        int idim;

        struct cOutputData {
          double NumberDesnity,ParticleNumber,v[3],MeanParticleSpeed,TranslationalTemeprature;
        } OutputData;

        if (pipe->ThisThread==CenterNodeThread) {
          OutputData.NumberDesnity=GetNumberDensity(DataSetNumber);
          OutputData.ParticleNumber=GetParticleNumber(DataSetNumber);
          GetBulkVelocity(OutputData.v,DataSetNumber);
          OutputData.MeanParticleSpeed=GetMeanParticleSpeed(DataSetNumber);
          OutputData.TranslationalTemeprature=GetTranslationalTemperature(DataSetNumber);
        }


        if (pipe->ThisThread==0) {
          if (CenterNodeThread!=0) pipe->recv((char*)&OutputData,sizeof(OutputData),CenterNodeThread);

          fprintf(fout,"%e  %e ",OutputData.NumberDesnity,OutputData.ParticleNumber);
          for (idim=0;idim<DIM;idim++) fprintf(fout,"%e ",OutputData.v[idim]);
          fprintf(fout,"%e %e ",OutputData.MeanParticleSpeed,OutputData.TranslationalTemeprature);
        }
        else pipe->send((char*)&OutputData,sizeof(OutputData));

        //print the user defind 'center node' data
        list<fPrintDataCenterNode>::iterator fptr;

        for (fptr=PrintDataCenterNode.begin();fptr!=PrintDataCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber,pipe,CenterNodeThread,this);

        //print data sampled by the user defined sampling functions
        if (PIC::IndividualModelSampling::PrintSampledData.size()!=0) {
          for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintSampledData.size();i++) PIC::IndividualModelSampling::PrintSampledData[i](fout,DataSetNumber,pipe,CenterNodeThread,this);
        }

      }
      */

      void PrintFileDescriptior(FILE* fout,int DataSetNumber) {
        char sym[_MAX_STRING_LENGTH_PIC_];

        PIC::MolecularData::GetChemSymbol(sym,DataSetNumber);
        fprintf(fout,"TITLE=\"specie=%s\"",sym);
      }


      void PrintVariableList(FILE* fout,int DataSetNumber);
      /*
      void PrintVariableList(FILE* fout,int DataSetNumber) {
       fprintf(fout,", \"Number Density\", \"Particle Number\"");
       for (int idim=0;idim<DIM;idim++) fprintf(fout,", \"V%i\"",idim);
       fprintf(fout,", \"Speed\", \"Translational Temperature\"");

       //print the user defind 'center node' data
       list<fPrintVariableListCenterNode>::iterator fptr;
       for (fptr=PrintVariableListCenterNode.begin();fptr!=PrintVariableListCenterNode.end();fptr++) (*fptr)(fout,DataSetNumber);

       //print varialbes sampled by the user defined sampling procedures
       if (PIC::IndividualModelSampling::PrintVariableList.size()!=0) {
         for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintVariableList.size();i++) PIC::IndividualModelSampling::PrintVariableList[i](fout,DataSetNumber);
       }

      }
      */

      void Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients);
      /*
      void Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients) {
        int i,s,idim;
        double c;



        //==============================  DEBUGGER ===============
                 if (nInterpolationCoeficients!=0) Measure=InterpolationList[0]->Measure;

        //============================== END DEBUGGER ============


        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if (associatedDataPointer==NULL) exit(__LINE__,__FILE__,"Error: The associated data buffer is not initialized");
        #endif

        double InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedBulkVelocity[3]={0.0,0.0,0.0},InterpolatedBulk2Velocity[3]={0.0,0.0,0.0};
        double InterpolatedParticleSpeed=0.0;
        double pWeight;

        for (s=0;s<PIC::nTotalSpecies;s++) {
          InterpolatedParticleWeight=0.0,InterpolatedParticleNumber=0.0,InterpolatedParticleNumberDeinsity=0.0,InterpolatedParticleSpeed=0.0;
          for (idim=0;idim<3;idim++) InterpolatedBulkVelocity[idim]=0.0,InterpolatedBulk2Velocity[idim]=0.0;

          //interpolate the sampled data
          for (i=0;i<nInterpolationCoeficients;i++) {
            pWeight=InterpolationList[i]->GetRealParticleNumber(s);
            c=PIC::LastSampleLength*InterpolationCoeficients[i];

            InterpolatedParticleWeight+=c*pWeight;
            InterpolatedParticleNumber+=c*InterpolationList[i]->GetParticleNumber(s);
            InterpolatedParticleNumberDeinsity+=c*InterpolationList[i]->GetNumberDensity(s);


            for (idim=0;idim<3;idim++) {
              InterpolatedBulkVelocity[idim]+=c*(*(idim+3*s+(double*)(InterpolationList[i]->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocityRelativeOffset)));
              InterpolatedBulk2Velocity[idim]+=c*(*(idim+3*s+(double*)(InterpolationList[i]->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocity2RelativeOffset)));
            }

            InterpolatedParticleSpeed+=c*(*(s+(double*)(InterpolationList[i]->associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleSpeedRelativeOffset)));
          }

          //stored the interpolated data in the associated data buffer
          *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNumberRelativeOffset))=InterpolatedParticleNumber;
          *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset))=InterpolatedParticleWeight;
          *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))=InterpolatedParticleNumberDeinsity;
          *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleSpeedRelativeOffset))=InterpolatedParticleSpeed;

          for (i=0;i<3;i++) {
            *(i+3*s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocityRelativeOffset))=InterpolatedBulkVelocity[i];
            *(i+3*s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocity2RelativeOffset))=InterpolatedBulk2Velocity[i];
          }

          *(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleSpeedRelativeOffset))=InterpolatedParticleSpeed;

        }

        //print the user defind 'center node' data
        list<fInterpolateCenterNode>::iterator fptr;

        for (fptr=InterpolateCenterNode.begin();fptr!=InterpolateCenterNode.end();fptr++) (*fptr)(InterpolationList,InterpolationCoeficients,nInterpolationCoeficients,this);

        //interpolate data sampled by user defiend sampling procedures
        if (PIC::IndividualModelSampling::InterpolateCenterNodeData.size()!=0) {
          for (unsigned int i=0;i<PIC::IndividualModelSampling::PrintVariableList.size();i++) PIC::IndividualModelSampling::InterpolateCenterNodeData[i](fout,DataSetNumber);
        }
      }

      */
    };


    //the class that contains the run information for the cell's corners
    class cDataCornerNode : public cBasicCornerNode {
    public:
    };

    //the data stored in a block
    //1. Local Time Step [NS]: depending on the model mode there will be a 'global' time step for the simulation, 'global' time step for the cell or individual time step for each simulated species
    //2. Local particle weight [NS]: depending on the model mode there will be a 'global' weight for the simulation, 'global' weight for the cell or individual weigh for each simulated species

    class cDataBlockAMR : public cBasicBlockAMR<cDataCornerNode,cDataCenterNode> {
    public:
      static int LocalTimeStepOffset,LocalParticleWeightOffset,totalAssociatedDataLength;
      char *associatedDataPointer;

      long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],tempParticleMovingListTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

      int AssociatedDataLength() {
        return totalAssociatedDataLength;
      }

      void SetAssociatedDataBufferPointer(char* ptr) {
        associatedDataPointer=ptr;
      }


      char* GetAssociatedDataBufferPointer() {
        return associatedDataPointer;
      }


      cDataBlockAMR () : cBasicBlockAMR<cDataCornerNode,cDataCenterNode> () {
        associatedDataPointer=NULL;
      }


      //exchenge of the data between processors
      void sendBoundaryLayerBlockData(CMPI_channel *pipe);
      void recvBoundaryLayerBlockData(CMPI_channel *pipe,int From);

      //send the block to abother processor
      void sendMoveBlockAnotherProcessor(CMPI_channel *pipe);
      void recvMoveBlockAnotherProcessor(CMPI_channel *pipe,int From);

      //clean the sampling buffers
      void cleanDataBuffer() {
        int i,length=totalAssociatedDataLength/sizeof(double);
        double *ptr;

        //clean the associated data buffers
        for (i=0,ptr=(double*)associatedDataPointer;i<length;i++,ptr++) *ptr=0.0;

        //clean the base class' data
        cBasicBlockAMR<cDataCornerNode,cDataCenterNode>::cleanDataBuffer();

        //clean the Particle Tables
        length=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
        for (i=0;i<length;i++) FirstCellParticleTable[i]=-1,tempParticleMovingListTable[i]=-1;
      }


      //set and get the local time step
      void SetLocalTimeStep(double dt, int spec);
      double GetLocalTimeStep(int spec);
      void SetLocalParticleWeight(double weight, int spec);
      double GetLocalParticleWeight(int spec);

      //print into a output file the blocks' parameters: the local time step, the local weight
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int BlockThread) {

        struct cOutputData {
          double dtLocal,wLocal;
        } OutputData;



        if (pipe->ThisThread==BlockThread) {
          OutputData.dtLocal=GetLocalTimeStep(DataSetNumber);
          OutputData.wLocal=GetLocalParticleWeight(DataSetNumber);
        }

        if (pipe->ThisThread==0) {
           if (BlockThread!=0) pipe->recv((char*)&OutputData,sizeof(OutputData),BlockThread);

           fprintf(fout,"%e  %e  ",OutputData.dtLocal,OutputData.wLocal);
         }
         else pipe->send((char*)&OutputData,sizeof(OutputData));
      }


      void PrintVariableList(FILE* fout) {
        fprintf(fout,", \"Local Time Step\" ");
        fprintf(fout,", \"log10(Local Particle Weight)\" ");
      }


    };


    //init the sampling buffers of the cell's data
    void initCellSamplingDataBuffer();
    void SetCellSamplingDataRequest();

    //return time step and the particle's weights
    inline double GetLocalTimeStep(int spec,cDataBlockAMR* block) {
      return *(spec+(double*)(cDataBlockAMR::LocalTimeStepOffset+block->GetAssociatedDataBufferPointer()));
    }

    inline void SetLocalTimeStep(double dt,int spec,cDataBlockAMR* block) {
      *(spec+(double*)(cDataBlockAMR::LocalTimeStepOffset+block->GetAssociatedDataBufferPointer()))=dt;
    }

    double GetLocalParticleWeight(int,cDataBlockAMR*);
    void SetLocalParticleWeight(double,int,cDataBlockAMR*);

    void flushCompletedSamplingBuffer(cDataCenterNode*);
    void flushCollectingSamplingBuffer(cDataCenterNode*);
    void switchSamplingBuffers();






    //the computational mesh
    #if DIM == 3
    extern cMeshAMR3d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #elif DIM == 2
    extern cMeshAMR2d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #else
    extern cMeshAMR1d<cDataCornerNode,cDataCenterNode,cDataBlockAMR > mesh;
    #endif



    //init the computational mesh
    void Init(double*,double*,fLocalMeshResolution);
    void buildMesh();
    void loadMesh(char*);

  }

  //volume injection of model particles
#if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  namespace VolumeParticleInjection {

    //init the model
    void Init();

    //Generte new particle internal properties
    typedef void (*fGenerateInternalParticleProperties)(long int ptr,int spec,int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
    extern fGenerateInternalParticleProperties GenerateInternalParticleProperties;

    //generate a random position in a cell
    void GetRandomCellPosition(double *x,int iCell,int jCell,int kCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //The list of volume injection processes
    //the injection rate of particle due to a specific injection process
    typedef void (*fSpeciesInjectionRate)(bool *InjectionFlag,double *Rate, int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the function that process that injection reaction (generate the particles)
    typedef long int (*fInjectionProcessor)(int iCellIndex,int jCellIndex,int kCellIndex,PIC::Mesh::cDataCenterNode *cell, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the limit of the local time step due to the injection process
    typedef double (*fLocalTimeStepLimit)(int spec,bool& TimeStepLimitationImposed, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

    //the array that stores the total injection rate by the volume injection
    extern double *SourceRate;

    struct cVolumeInjectionDescriptor {
      fSpeciesInjectionRate SpeciesInjectionRate;
      fInjectionProcessor InjectionProcessor;
      fLocalTimeStepLimit LocalTimeStepLimit;
    };

    const int nMaxInjectionProcessEntries=128;
    extern int nRegistratedInjectionProcesses;

    extern cVolumeInjectionDescriptor VolumeInjectionDescriptor[nMaxInjectionProcessEntries];

    void inline RegisterVolumeInjectionProcess(fSpeciesInjectionRate f1, fInjectionProcessor f2, fLocalTimeStepLimit f3) {
      if (nRegistratedInjectionProcesses==nMaxInjectionProcessEntries-1) {
        exit(__LINE__,__FILE__,"Error: the volume injection processes buffer is overflow, increase the length of the buffer ('PIC::VolumeParticleInjection::nMaxInjectionProcessEntries'");
      }

      if ((f1==NULL)||(f2==NULL)||(f3==NULL)) exit(__LINE__,__FILE__,"Error: one of the functions is not defined");

      if (SourceRate==NULL) {
        if (PIC::nTotalSpecies==0) exit(__LINE__,__FILE__,"Error: neede to know the total number of species before initializing the sampling buffer");

        SourceRate=new double[PIC::nTotalSpecies];
        for (int s=0;s<PIC::nTotalSpecies;s++) SourceRate[s]=0.0;
      }

      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].SpeciesInjectionRate=f1;
      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].InjectionProcessor=f2;
      VolumeInjectionDescriptor[nRegistratedInjectionProcesses].LocalTimeStepLimit=f3;

      nRegistratedInjectionProcesses++;
    }


    long int  InjectParticle();

    //paericle weight injection rates
    void InitTotalInjectionRate();
    double GetTotalInjectionRate(int);
    double GetTotalTimeStepInjection(int spec);
    double GetBlockInjectionRate(int spec,PIC::Mesh::cDataBlockAMR *block);
    double GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell);
  }
#endif

  //sampling functions
  namespace Sampling {

    //the minimum number of iteration for output of the datafile
    extern int minIterationNumberForDataOutput;

    //'SaveOutputDataFile' determines weather the data file will be created for a particular species (output of data files can be suppress for particular species, such as external, dust....)
    extern bool *SaveOutputDataFile;

    //sample the normal and tangential kinetic temperatures: constant origin of the direction of the normal
    static const double constNormalDirection__SampleParallelTangentialTemperature[3]={0.0,0.0,0.0};

    namespace ExternalSamplingLocalVariables {

      //the external procedures for sampling particle data
      typedef void (*fSamplingProcessor) ();

      //the procedure that prints the sampled data into a file
      typedef void (*fPrintOutputFile)(int);

      extern const int nMaxSamplingRoutines;
      extern int SamplingRoutinesRegistrationCounter;

      extern fSamplingProcessor *SamplingProcessor;
      extern fPrintOutputFile *PrintOutputFile;

      inline void RegisterSamplingRoutine(fSamplingProcessor p,fPrintOutputFile f) {
        if (SamplingRoutinesRegistrationCounter==0) {
          SamplingProcessor=new fSamplingProcessor[nMaxSamplingRoutines];
          PrintOutputFile=new fPrintOutputFile[nMaxSamplingRoutines];

          for (int i=0;i<nMaxSamplingRoutines;i++) SamplingProcessor[i]=NULL,PrintOutputFile[i]=NULL;
        }
        else if (SamplingRoutinesRegistrationCounter==nMaxSamplingRoutines-1) {
          exit(__LINE__,__FILE__,"Error: SamplingRoutinesRegistrationCounter exeeds its maximum value: increse the value of PIC::Sampling::ExternalSamplingProcedures::nMaxSamplingRoutines");
        }

        SamplingProcessor[SamplingRoutinesRegistrationCounter]=p;
        PrintOutputFile[SamplingRoutinesRegistrationCounter]=f;

        SamplingRoutinesRegistrationCounter++;
      }


    }



    void Sampling();

    //sample the particle data
    inline void SampleParticleData(char* ParticleData,char *SamplingBuffer,PIC::Mesh::cDataBlockAMR *block,PIC::Mesh::cDataCenterNode *cell,double TimeStepFraction) {
      double Speed2,*v,LocalParticleWeight,v2;
      int s,idim;
      double *sampledVelocityOffset,*sampledVelocity2Offset;


      Speed2=0.0;

      s=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)ParticleData);
      v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)ParticleData);

      LocalParticleWeight=block->GetLocalParticleWeight(s);
      LocalParticleWeight*=TimeStepFraction*PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)ParticleData);

      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberRelativeOffset))+=1;
      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;


      sampledVelocityOffset=3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocityRelativeOffset);
      sampledVelocity2Offset=3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocity2RelativeOffset);

      for (idim=0;idim<3;idim++) {
        v2=v[idim]*v[idim];
        Speed2+=v2;

        *(idim+sampledVelocityOffset)+=v[idim]*LocalParticleWeight;
        *(idim+sampledVelocity2Offset)+=v2*LocalParticleWeight;
      }

      *(s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleSpeedRelativeOffset))+=sqrt(Speed2)*LocalParticleWeight;

    }
  }

  //colecular collisions
  namespace MolecularCollisions {

    //collisions between model particles
    namespace ParticleCollisionModel {
      //sample collision statistics
      extern int CollsionFrequentcySamplingOffset;

      int RequestSamplingData(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void Init();

      //model of the particle collisions
      void ntc();
    }

    //models for calculation of the relative velocity after a collision
    namespace VelocityScattering {
      namespace HS {
        inline void VelocityAfterCollision(double *vrel,int s0,int s1) {
          double Vrc,V[3];
          double CosKsi,SinKsi,CosEps,SinEps,D,c;

          CosKsi=2.0*rnd()-1.0;
          SinKsi=sqrt(1.0-CosKsi*CosKsi);

          c=2*Pi*rnd();
          SinEps=sin(c);
          CosEps=cos(c);

          D=sqrt(vrel[1]*vrel[1]+vrel[2]*vrel[2]);
          if (D>1.0E-6) {
            Vrc=sqrt(vrel[0]*vrel[0]+vrel[1]*vrel[1]+vrel[2]*vrel[2]);
            V[0]=CosKsi*vrel[0]+SinKsi*SinEps*D;
            V[1]=CosKsi*vrel[1]+SinKsi*(Vrc*vrel[2]*CosEps-vrel[0]*vrel[1]*SinEps)/D;
            V[2]=CosKsi*vrel[2]-SinKsi*(Vrc*vrel[1]*CosEps+vrel[0]*vrel[2]*SinEps)/D;
          }
          else {
            V[0]=CosKsi*vrel[0];
            V[1]=SinKsi*CosEps*vrel[0];
            V[2]=SinKsi*SinEps*vrel[0];
          }

          memcpy(vrel,V,3*sizeof(double));
        }
      }
    }

    //collisions with the background atmosphere
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
    namespace BackgroundAtmosphere {

      //the total number of the background species, the mass table, the table of cosntant collision cross sections with the model species
      static const int nTotalBackgroundSpecies=1;
      static const double BackgroundSpeciesMassTable[]={0.0};
      static const double BackgroundAtmosphereConstantCrossSectionTable[PIC::nTotalSpecies][nTotalBackgroundSpecies];
      static const int Background2ModelSpeciesConversionTable[]={-1};

      inline int GetTotalNumberBackgroundSpecies() {return nTotalBackgroundSpecies;}
      inline double GetBackgroundMolecularMass(int spec) {return BackgroundSpeciesMassTable[spec];}


/*
      //the default value of the user-defined function that calculates the collision cross section
      #define _PIC_BACKGROUND_ATMOSPHERE__COLLISION_CROSS_SECTION_FUNCTION_(spec,BackgroundSpecieNumber,modelParticleData,BackgroundAtmosphereParticleData,TranslationalEnergy,cr2) (0.0)

      //get local, cell mean, and cell maximum density of the background species
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_LOCAL_DENSITY_(x,BackgroundSpecieNumber,cell,node) (0.0)
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_CELL_MEAN_DENSITY_(BackgroundSpecieNumber,cell,node) (0.0)
      #define _PIC_BACKGROUND_ATMOSPHERE__BACKGROUND_SPECIES_CELL_MAXIMUM_DENSITY_(BackgroundSpecieNumber,cell,node) (0.0)

      //the user-defined function for calcualtion of the scattering angle
      #define _PIC_BACKGROUND_ATMOSPHERE__COLLISION_SCATTERING_ANGLE_(Vrel,TranslationalEnergy,spec,BackgroundSpecieNumber) (0.0)

      //the condition the remove the model particle from the simulation after the collision with the background particle
      #define _PIC_BACKGROUND_ATMOSPHERE__REMOVE_CONDITION_MODEL_PARTICLE_(modelParticleData) (true)

      //the conditions to inject the background atmosphere particle after a collision with the model particle
      #define _PIC_BACKGROUND_ATMOSPHERE__INJECT_CONDITION_BACKGROUND_PARTICLE_(BackgroundAtmosphereParticleData) (false)

      //Evaluate GetSigmaCrMax in a cell
      #define _PIC_BACKGROUND_ATMOSPHERE__GET_SIGMA_CR_MAX(spec,BackgroundSpecieNumber,modelParticleData) (0.0)

      //generate the background atmosphere particle
      #define _PIC_BACKGROUND_ATMOSPHERE__GENERATE_BACKGROUND_PARTICLE_(BackgroundAtmosphereParticleData,BackgroundSpecieNumber,cell,node) (exit(__LINE__,__FILE__,"Error: _PIC_BACKGROUND_ATMOSPHERE__GENERATE_BACKGROUND_PARTICLE_ was not set"))

      //define the mode for loading the definition file
      #define _PIC_BACKGROUND_ATMOSPHERE__LOAD_USER_DEFINITION__MODE_ _PIC_MODE_OFF_
      #define _PIC_BACKGROUND_ATMOSPHERE__UDER_DEFINITION_ "UserDefinition.PIC.BackgroundAtmosphere.h"
*/

      //define functions for calculation of the properties of the background atmosphere species
      double GetCollisionCrossSectionBackgoundAtmosphereParticle(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData,PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,double TranslationalEnergy,double cr2);
      double GetSigmaCrMax(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData);
      double GetCollisionScatteringAngle(double* Vrel,double TranslationalEnergy,int spec,int BackgroundSpecieNumber);

      void GenerateBackgoundAtmosphereParticle(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

      double GetBackgroundLocalNumberDensity(int BackgroundSpecieNumber,double *x);
      double GetCellMeanBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      double GetCellMaximumBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
      double GetCellLocalBackgroundNumberDensity(double x[3],int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

      //the conditions to keep the background particle or remove a model particle after a collision
      bool KeepConditionModelParticle(PIC::ParticleBuffer::byte *ModelParticleData);



      //include the user defined properties of the background atmosphere
      #if _PIC_BACKGROUND_ATMOSPHERE__LOAD_USER_DEFINITION__MODE_ ==  _PIC_MODE_ON_
      #include _PIC_BACKGROUND_ATMOSPHERE__UDER_DEFINITION_
      #endif


      //Sampling of the model data
      extern int LocalTotalCollisionFreqSamplingOffset;
      int RequestSamplingData(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

      //sample the global rate of termalization of exospehric species and the rate of injection of new exospheric particles from the background atmosphere
      extern double *ThermalizationRate,*CollisionSourceRate;

      //init the model
      void Init_BeforeParser();
      void Init_AfterParser();

      //output sampled model parameters
      void OutputSampledModelData(int);
      void SampleModelData();

      //processor of the collisions
      void CollisionProcessor();
      void RemoveThermalBackgroundParticles();
    }
#endif

  }

  namespace IndividualModelSampling {

    //reserve memory to store sampling data
    typedef int (*fRequestSamplingData)(int);
    extern vector<fRequestSamplingData> RequestSamplingData;

    //reserve memoty in a cell associated data buffer for non-sampling data
    typedef int (*fRequestStaticCellData)(int);
    extern vector<fRequestStaticCellData> RequestStaticCellData;

    //the list of user defined sampling procedures
    typedef void (*fSamplingProcedure)();
    extern vector<fSamplingProcedure> SamplingProcedure;

    //print the variable list
    typedef void (*fPrintVariableList)(FILE* fout,int nDataSet);
    extern vector<fPrintVariableList> PrintVariableList;

    //interpolate center node data
    typedef void (*fInterpolateCenterNodeData)(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode* cell);
    extern vector<fInterpolateCenterNodeData> InterpolateCenterNodeData;

    //print the sampled node data
    typedef void (*fPrintSampledData)(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode* cell);
    extern vector<fPrintSampledData> PrintSampledData;
  }

  //sample and output the particle's distribution function
  namespace DistributionFunctionSample {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    extern int v2SamplingMode,speedSamplingMode;

    //the range of the velocity scale and the number of nodes in the sample
    extern double vMax,vMin;
    extern long int nSampledFunctionPoints;
    extern double dV,dV2,dSpeed;

    //the sampling buffers
    extern double **SamplingBuffer;

    //sampling data offsets
    extern int Sample_Velocity_Offset,Sample_Speed_Offset,Sample_V2_Offset,SampleDataLength;

    //get the offset to the beginig of the sampling data for a particular samplePoint, spec,.....
    long int GetSampleDataOffset(int spec,int nInterval);

    //sampling  locations
    extern double **SamplingLocations;
    extern int nSamleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void Init(double ProbeLocations[][DIM],int nProbeLocations);

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
  }


  //sample and output the particle's distribution of particle's flux
  namespace ParticleFluxDistributionSample {
    //the init flag
    extern bool SamplingInitializedFlag;

    //the modes for sampling of the v^2 and the absolute value of velocity
    extern const int _LINEAR_SAMPLING_SCALE_,_LOGARITHMIC_SAMPLING_SCALE_;
    extern int v2SamplingMode,speedSamplingMode;

    //the range of the velocity scale and the number of nodes in the sample
    extern double vMax,vMin;
    extern long int nSampledFunctionPoints;
    extern double dV2,dSpeed;

    //the sampling buffers
    extern double **SamplingBuffer;
    extern double **SamplingFlux;

    //sampling data offsets
    extern int Sample_Speed_Offset,Sample_V2_Offset,SampleDataLength;

    //get the offset to the beginig of the sampling data for a particular samplePoint, spec,.....
    long int GetSampleDataOffset(int spec,int nInterval);

    //sampling  locations
    extern double **SamplingLocations;
    extern double **SamplingPointingDirections;
    extern double maxSamplingConeAngle,cosMaxSamplingConeAngle;

    extern int nSamleLocations;
    extern cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** SampleNodes;
    extern long int *SampleLocalCellNumber;

    void Init(double ProbeLocations[][DIM],double ProbeDirections[][DIM],double maxConeAngles,int nProbeLocations);

    void SampleDistributionFnction();
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
    void printMacroscopicParameters(char *fname,int spec);
  }


  //procedures for distribution of particle velocities
  namespace Distribution {

    //the macroses defined the modes for generation of the particle parameters
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_ -> the particle properties are generated with the requested distributuion, no weight correction are needed
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_ -> to recover the requasted distribution, a particle weght correction is needed

    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_         0
    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_   1


    void MaxwellianVelocityDistribution(double *v,const double *BulkFlowVelocity,const double Temp,const int spec);
    double InjectMaxwellianDistribution(double *v,const double *BulkFlowVelocity,const double Temp,double *ExternalNormal,int spec,int WeightCorrectionMode=_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_);
  }


  /*
  namespace ParticleBuffer {
    typedef unsigned char byte;

    //the set of parameters that determine the 'double' parameters of a particle
    #define _PIC_PARTICLE_DATA_SPECIEDID_OFFSET_ 0
    #define _PIC_PARTICLE_DATA_NEXT_OFFSET_ (_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_ + sizeof(unsigned char))
    #define _PIC_PARTICLE_DATA_PREV_OFFSET_ (_PIC_PARTICLE_DATA_NEXT_OFFSET_ + sizeof(long int))
    #define _PIC_PARTICLE_DATA_VELOCITY_OFFSET_  (_PIC_PARTICLE_DATA_PREV_OFFSET_ + sizeof(long int))
    #define _PIC_PARTICLE_DATA_POSITION_OFFSET_  (_PIC_PARTICLE_DATA_VELOCITY_OFFSET_+ 3*sizeof(double))

    //the offset for the variable that contains the 'particle weight correction'
    extern int _PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_;

    //the total length of a data allocated for a particle
    extern long int ParticleDataLength;

    //The particle buffer's internal data
    extern byte *ParticleDataBuffer;
    extern long int MaxNPart,NAllPart,FirstPBufferParticle;

    //Request additional data for a particle
    void RequestDataStorage(long int &offset,int TotalDataLength);

    //the basic data access functions for a particle
    byte *GetParticleDataPointer(long int);

    double *GetX(long int);
    void GetX(double*,long int);
    void SetX(double*,long int);

    double *GetX(byte*);
    void GetX(double*,byte*);
    void SetX(double*,byte*);

    double *GetV(long int);
    void GetV(double*,long int);
    void SetV(double*,long int);

    double *GetV(byte*);
    void GetV(double*,byte*);
    void SetV(double*,byte*);

    unsigned int GetI(byte*);
    void SetI(unsigned int,byte*);

    unsigned int GetI(long int);
    void SetI(unsigned int,long int);

    long int GetPrev(long int);
    long int GetNext(long int);
    void SetPrev(long int,long int);
    void SetNext(long int,long int);

    long int GetPrev(byte*);
    long int GetNext(byte*);
    void SetPrev(long int,byte*);
    void SetNext(long int,byte*);

    double GetIndividualStatWeightCorrection(long int);
    double GetIndividualStatWeightCorrection(byte*);
    void SetIndividualStatWeightCorrection(double,long int);
    void SetIndividualStatWeightCorrection(double,byte*);

    //the particle buffer procedure
    void Init(long int);
    long int GetMaxNPart();
    long int GetAllPartNum();
    long int GetParticleDataLength();

    long int GetNewParticle();
    long int GetNewParticle(long int&);

    void DeleteParticle(long int);
    void DeleteParticle(long int,long int&);

    void CloneParticle(long int,long int);

    void ExcludeParticleFromList(long int,long int&);

    void SaveImageFile(int);
    void LoadImageFile(int);

    void PackParticleData(char*,long int);
    void UnPackParticleData(char*,long int);

    unsigned long GetChecksum();
  }

*/


  //the mode of the internal degrees of freedom
  namespace IDF {

    static const int nTotalVibtationalModes[]={-1};
    static const int nTotalRotationalModes[]={-1};
    static const int nSpeciesMaxVibrationalModes=0;
    static const double CharacteristicVibrationalTemperature[]={0.0};    //the rule of access CharacteristicVibrationalTemperature[nmode+s*nSpeciesMaxVibrationalModes]
    static const double RotationZnumber[]={0.0};

    extern int _ROTATIONAL_ENERGY_SAMPLE_DATA_OFFSET_;
    extern int _VIBRATIONAL_ENERGY_SAMPLE_DATA_OFFSET_[PIC::nTotalSpecies];
    extern int _TOTAL_SAMPLE_PARTICLE_WEIGHT_SAMPLE_DATA_OFFSET_;

    namespace LB {
      extern int _ROTATIONAL_ENERGY_OFFSET_,_VIBRATIONAL_ENERGY_OFFSET_;


      inline double GetRotE(PIC::ParticleBuffer::byte *ParticleDataStart) {
        return *((double*)(ParticleDataStart+_ROTATIONAL_ENERGY_OFFSET_));
      }

      inline void SetRotE(double e,PIC::ParticleBuffer::byte *ParticleDataStart) {
         *((double*)(ParticleDataStart+_ROTATIONAL_ENERGY_OFFSET_))=e;
      }

      inline double GetVibE(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
        if (nmode>=0) return *(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_));

        double res=0.0;
        int n,nVibModes,s;

        s=PIC::ParticleBuffer::GetI(ParticleDataStart);
        nVibModes=nTotalVibtationalModes[s];

        for (n=0;n<nVibModes;n++) res+=*(n+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_));
        return res;
      }

      inline void SetVibE(double e,int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
        *(nmode+(double*)(ParticleDataStart+_VIBRATIONAL_ENERGY_OFFSET_))=e;
      }

      void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart);
      void InitRotTemp(double RotTemp,PIC::ParticleBuffer::byte *ParticleDataStart);

      double GetCellRotTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      double GetCellMeanRotE(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellMeanVibE(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell);

      //request data for the model
      int RequestSamplingData(int offset);

      //init
      void Init();

      //output the model data
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintVariableList(FILE* fout,int DataSetNumber);

      //calcualte the temperature index
      //get the temperature index
      inline double GetTempIndex(int s0,int s1) {
        static const double TemepratureIndex[1][1]={0.0};

        return TemepratureIndex[s0][s1];
      }
    }

    namespace qLB {
    using namespace LB;

//      extern int NmaxVibModes; //the maximum number of the vibrational modes for the set of used species
      extern int _VIBRATIONAL_GROUND_LEVEL_SAMPLE_DATA_OFFSET_,_VIBRATIONAL_FIRST_EXITED_LEVEL_SAMPLE_DATA_OFFSET_;

      int RequestSamplingData(int offset);
      void Init();

      double GetVibE(int nmode_in,PIC::ParticleBuffer::byte *ParticleDataStart);
      void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart);
      void SetVibLevel(double VibQuantumNumber,int nmode,PIC::ParticleBuffer::byte *ParticleDataStart);
      int GetVibLevel(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart);

      //the ground and the first exited vibrational states
      double GetGroundVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);
      double GetFirstExitedVibLevelPopulationFraction(int nmode,int s,PIC::Mesh::cDataCenterNode* cell);

      double GetCellVibTemp(int s,PIC::Mesh::cDataCenterNode* cell);
      double GetCellVibTemp(int nmode_in,int s,PIC::Mesh::cDataCenterNode* cell);

      void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell);
    }

    inline double GetRotE(PIC::ParticleBuffer::byte *ParticleDataStart) {
      return LB::GetRotE(ParticleDataStart);
    }

    inline double GetVibE(int nmode,PIC::ParticleBuffer::byte *ParticleDataStart) {
      return LB::GetVibE(nmode,ParticleDataStart);
    }

    inline void RedistributeEnergy(PIC::ParticleBuffer::byte *ptr0,PIC::ParticleBuffer::byte *ptr1,double& vrel,bool* ChangeParticlePropertiesFlag,PIC::Mesh::cDataCenterNode* cell) {
      LB::RedistributeEnergy(ptr0,ptr1,vrel,ChangeParticlePropertiesFlag,cell);
    }

    inline void InitRotTemp(double RotTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
      LB::InitRotTemp(RotTemp,ParticleDataStart);
    }

    inline void InitVibTemp(double VibTemp,PIC::ParticleBuffer::byte *ParticleDataStart) {
      LB::InitVibTemp(VibTemp,ParticleDataStart);
    }

    inline void Init() {
      LB::Init();
    }
  }


  namespace Mover {

//    #include "UserDefinition.PIC.Mover.h"

    //the return codes of the moving procedures
    #define _PARTICLE_REJECTED_ON_THE_FACE_ -1
    #define _PARTICLE_DELETED_ON_THE_FACE_   0
    #define _PARTICLE_CROSS_THE_FACE_        1
    #define _PARTICLE_LEFT_THE_DOMAIN_       2
    #define _PARTICLE_MOTION_FINISHED_       3

//    typedef void (*fTotalParticleAcceleration)(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
    typedef int (*fSpeciesDependentParticleMover) (long int,double,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    typedef int (*fSpeciesDependentParticleMover_BoundaryInjection) (long int,double,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,bool);


    //the vector containing the species specific particle moving procedures
//    extern fSpeciesDependentParticleMover *MoveParticleTimeStep;
//    extern fTotalParticleAcceleration TotalParticleAcceleration;
//    extern fSpeciesDependentParticleMover_BoundaryInjection *MoveParticleBoundaryInjection;

    //process a particle when it leaves the boundary of the computational domain
    typedef int (*fProcessOutsideDomainParticles) (long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
    extern fProcessOutsideDomainParticles ProcessOutsideDomainParticles;

    typedef int (*fProcessTriangleCutFaceIntersection) (long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace);
    extern fProcessTriangleCutFaceIntersection ProcessTriangleCutFaceIntersection;


    void Init();
    void MoveParticles();

    int UniformWeight_UniformTimeStep_noForce(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag);

    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,bool FirstBoundaryFlag);


    void TotalParticleAcceleration_default(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
  }


  namespace ParticleWeightTimeStep {
    typedef double (*fSetFunction) (int,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);

    extern fSetFunction LocalParticleWeight,LocalTimeStep,LocalBlockInjectionRate;
    extern double maxReferenceInjectedParticleNumber;


    //when the global particle weight/time step are used, the following are the buffers where these parameters are stored
    extern double *GlobalParticleWeight,*GlobalTimeStep;

    double GetMaximumBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);
    double GetTotalBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    void initParticleWeight(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    void initParticleWeight_ConstantWeight();
    void initParticleWeight_ConstantWeight(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);


    void initTimeStep(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

    void copyLocalParticleWeightDistribution(int specTarget,int specSource,double ProportionaltyCoefficient=1.0);
    void copyLocalTimeStepDistribution(int specTarger,int specSource,double ProportionaltyCoefficient=1.0);

  }

  namespace Parallel {


     //count the number of particles that were send and recieve by the thread
     extern long int sendParticleCounter,recvParticleCounter,IterationNumberAfterRebalancing;
     extern double RebalancingTime,CumulativeLatency;

     //the factor the trrigeres the emergency load rebalancing. The condition for the rebalancing:
     //(PIC::Parallel::CumulativeLatency>PIC::Parallel::EmergencyLoadRebalancingFactor*PIC::Parallel::RebalancingTime)
     extern double EmergencyLoadRebalancingFactor;

     //exchenge paricles between iterations
     void ExchangeParticleData();

     //Latency of the run
     extern double Latency;
  }

  namespace Debugger {
    //contains functions that are used for debugging the code

    //InfiniteLoop==false ==> no problem found; InfiniteLoop==true => the actual number of particles does not consider the that in teh particle buffer
    bool InfiniteLoop(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL);
    void FindDoubleReferencedParticle(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=NULL);
  }

  namespace Alarm {
    //finish the execution of the code at a particular value of the walltime
    extern bool AlarmInitialized,WallTimeExeedsLimit;
    extern double StartTime;
    extern double RequestedExecutionWallTime;

    inline void SetAlarm(double requestedWallTime) {
      AlarmInitialized=true;
      StartTime=MPI_Wtime();
      RequestedExecutionWallTime=requestedWallTime;
    }

    inline void FinishExecution() {
      MPI_Finalize();
      printf("$PREFIX:!!!!! Execution is finished by the alarm (PIC::Alarm) !!!!!!\n");
      exit(__LINE__,__FILE__,"!!!!! Execution is finished by the alarm !!!!!!");
    }
  }

  namespace ColumnIntegration {

    //define 3 nodes on the surface of a bounding plane; index value: 0 -> xmin component of the coordinate, 1 -> xmax component of the coordinate
    static const int x0PlaneNodeIndex[6][3]={ {0,0,0},{1,0,0},       {0,0,0},{0,1,0},           {0,0,0},{0,0,1}};
    static const int x1PlaneNodeIndex[6][3]={ {0,1,0},{1,1,0},       {1,0,0},{1,1,0},           {1,0,0},{1,0,1}};
    static const int x2PlaneNodeIndex[6][3]={ {0,0,1},{1,0,1},       {0,0,1},{0,1,1},           {0,1,0},{0,1,1}};
    static const int PlaneNormal[6][3]=     { {1,0,0},{1,0,0},       {0,1,0},{0,1,0},           {0,0,1},{0,0,1}};

    struct cBoundingBoxFace {
      double x[3];
      double e0[3];
      double e1[3];
      double Normal[3];
      double e0Length;
      double e1Length;
    };

    extern cBoundingBoxFace BoundingBoxFace[6];
    extern bool InitializedFlag;

    //control initialization of the model
    extern bool ModelInitFlag;

    void Init();
    bool FindIntegrationLimits(double *x0,double *l,double& IntegrationPathLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xStartNode,double *xFinish,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &xFinishNode);

//    //get a single integrated value
//    double GetCoulumnIntegral(double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xStartNode,double *l,double IntegrationPathLength,double (*Integrand)(double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
//    double GetCoulumnIntegral(double *x0,double *l,double (*Integrand)(double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));

    //get values for multiple integrals
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *xStart,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xStartNode,double *l,double IntegrationPathLength,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
    void GetCoulumnIntegral(double *ResultVector,int ResultVectorLength,double *x0,double *l,void (*Integrand)(double*,int,double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*));
  }


  //namespace CPLR contains definitions of all couplers used in AMPS
  namespace CPLR {

    //coupling with SWMF
    namespace SWMF {
      extern int MagneticFieldOffset,TotalDataLength,PlasmaDensityOffset,BulkVelocityOffset,PlasmaPressureOffset;

      //init the coupler
      void init();
      void ConvertMpiCommunicatorFortran2C(signed int* iComm,signed int* iProc,signed int* nProc);

      //output the interpolated data into a file
      int RequestDataBuffer(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

      //prepare the list of the coordinates for the interpolation
      void ResetCenterPointProcessingFlag();
      void GetCenterPointNumber(int *nCenterPoints);
      void GetCenterPointCoordinates(double *x);
      void RecieveCenterPointData(char* ValiableList, int nVarialbes,double *data,int *index);

      //calcualte physical parameters
      inline void GetBackgroundMagneticField(double *B,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
        register int idim;
        register double *offset=(double*)(MagneticFieldOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

        for (idim=0;idim<3;idim++) B[idim]=offset[idim];
      }

      inline void GetBackgroundPlasmaVelocity(double *v,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
        register int idim;
        register double *offset=(double*)(BulkVelocityOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

        for (idim=0;idim<3;idim++) v[idim]=offset[idim];
      }

      inline void GetBackgroundElectricField(double *E,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
        double B[3],v[3];

        GetBackgroundMagneticField(B,x,nd,node);
        GetBackgroundPlasmaVelocity(v,x,nd,node);

        E[0]=-(v[1]*B[2]-B[1]*v[2]);
        E[1]=-(-v[0]*B[2]+B[0]*v[2]);
        E[2]=-(v[0]*B[1]-B[0]*v[1]);
      }

    }


    //coupling of AMPS through the ICES tool
    namespace ICES {
       extern char locationICES[_MAX_STRING_LENGTH_PIC_]; //location of the data and the dace cases

       void Init();
       void SetLocationICES(const char*);

       //the total number of bytes used to store the ICES data vector; the offset of the data associated with the ICES data vector
       extern int TotalAssociatedDataLength,AssociatedDataOffset;

       //the offsets for the plasma parameters loaded with ICES
       extern int ElectricFieldOffset,MagneticFieldOffset,PlasmaPressureOffset,PlasmaNumberDensityOffset,PlasmaTemperatureOffset,PlasmaBulkVelocityOffset,DataStatusOffsetSWMF;

       //the offsets for parameters loaded from the DSMC model
       extern int NeutralBullVelocityOffset,NeutralNumberDensityOffset,NeutralTemperatureOffset,DataStatusOffsetDSMC;

       //calcualte the total number of cells in the mesh
       long int getTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

       //create the trajectory file
       void createCellCenterCoordinateList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);


       //retrive the SWMF data file
       #define _PIC_ICES__STATUS_OK_   0

       class cDataNodeSWMF {
       public:
         double swNumberDensity,swTemperature,swPressure,E[3],B[3],swVel[3];
         int status;

         void flush() {
           swNumberDensity=0.0,swTemperature=0.0,swPressure=0.0;
           for (int i=0;i<3;i++) E[i]=0.0,B[i]=0.0,swVel[i]=0.0;
         }
       };

       class cDataNodeDSMC {
       public:
         double neutralNumberDensity,neutralTemperature,neutralVel[3];
         int status;

         void flush() {
           neutralNumberDensity=0.0,neutralTemperature=0.0;

           for (int i=0;i<3;i++) neutralVel[i]=0.0;
         }
       };

       void retriveSWMFdata(const char *DataFile);
       void retriveDSMCdata(const char *Case,const char *DataFile,const char *MeshFile);


       void readSWMFdata(const double MeanIonMass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree); //MeanIonMass -> the mean ion mass of the plasma flow in [amu]
       void readDSMCdata(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);

       //user defined pre-processor of the data that is readed by ICES
       typedef void (*fDSMCdataPreProcessor)(double *x,cDataNodeDSMC& data);
       typedef void (*fSWMFdataPreProcessor)(double *x,cDataNodeSWMF& data);

       extern fDSMCdataPreProcessor DSMCdataPreProcessor;
       extern fSWMFdataPreProcessor SWMFdataPreProcessor;

       void PrintVariableList(FILE* fout,int DataSetNumber);
       void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
       void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);

       //print the ion flux at a sphere
       void PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius);

       //calculate the values of the located parameters
       inline void GetBackgroundElectricField(double *E,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         register int idim;
         register double *offset=(double*)(ElectricFieldOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

         for (idim=0;idim<3;idim++) E[idim]=offset[idim];
       }

       inline void GetBackgroundMagneticField(double *B,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         register int idim;
         register double *offset=(double*)(MagneticFieldOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

         for (idim=0;idim<3;idim++) B[idim]=offset[idim];
       }

       inline void GetBackgroundPlasmaVelocity(double *vel,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         register int idim;
         register double *offset=(double*)(PlasmaBulkVelocityOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer());

         for (idim=0;idim<3;idim++) vel[idim]=offset[idim];
       }

       inline double GetBackgroundPlasmaPressure(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         return *((double*)(PlasmaPressureOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer()));
       }

       inline double GetBackgroundPlasmaNumberDensity(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         return *((double*)(PlasmaNumberDensityOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer()));
       }

       inline double GetBackgroundPlasmaTemperature(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         return *((double*)(PlasmaTemperatureOffset+node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer()));
       }

       inline void GetBackgroundFieldsVector(double *E,double *B,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
         register int idim;
         register char *offset=node->block->GetCenterNode(nd)->GetAssociatedDataBufferPointer();

         double *e=(double*)(offset+ElectricFieldOffset);
         double *b=(double*)(offset+MagneticFieldOffset);

         for (idim=0;idim<3;idim++) B[idim]=b[idim],E[idim]=e[idim];
       }
    }


    inline void GetBackgroundElectricField(double *E,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
      #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
      ICES::GetBackgroundElectricField(E,x,nd,node);
      #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
      SWMF::GetBackgroundElectricField(E,x,nd,node);
      #else
      exit(__LINE__,__FILE__,"not implemented");
      #endif
     }

     inline void GetBackgroundMagneticField(double *B,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       ICES::GetBackgroundMagneticField(B,x,nd,node);
       #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
       SWMF::GetBackgroundMagneticField(B,x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif
     }


     inline void GetBackgroundPlasmaVelocity(double *vel,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       ICES::GetBackgroundPlasmaVelocity(vel,x,nd,node);
       #elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
       SWMF::GetBackgroundPlasmaVelocity(vel,x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif
     }

     inline double GetBackgroundPlasmaPressure(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       double res=0.0;

       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       res=ICES::GetBackgroundPlasmaPressure(x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif

       return res;
     }

     inline double GetBackgroundPlasmaNumberDensity(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       double res=0.0;

       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       res=ICES::GetBackgroundPlasmaNumberDensity(x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif

       return res;
     }

     inline double GetBackgroundPlasmaTemperature(double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       double res=0.0;

       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       res=ICES::GetBackgroundPlasmaTemperature(x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif

       return res;
     }

     inline void GetBackgroundFieldsVector(double *E,double *B,double *x,long int nd,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
       #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__ICES_
       GetBackgroundFieldsVector(E,B,x,nd,node);
       #else
       exit(__LINE__,__FILE__,"not implemented");
       #endif
     }


  }

  namespace ChemicalReactions {

    namespace PhotolyticReactions {
      //the photolytic reactions: the model is initialized if _PIC_PHOTOLYTIC_REACTIONS_MODE_ = _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
      //the total photolytic lifetime of the species;
      extern double *ConstantTotalLifeTime;

/*
      typedef double (*fTotalLifeTime)(double *x,int spec,long int ptr,bool &ReactionAllowedFlag);
      extern fTotalLifeTime *TotalLifeTime;

      typedef int (*fReactionProcessor)(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData);
      extern fReactionProcessor *ReactionProcessorTable;
*/

      void Init();
//      void SetReactionProcessor(fReactionProcessor f,int spec);
//      void SetSpeciesTotalPhotolyticLifeTime(fTotalLifeTime f,int spec);

      //The return codes of the photolytic model
      #define _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_      0
      #define _PHOTOLYTIC_REACTION_OCCURES_                  1

      #define _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_        2
      #define _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_ 3

      //the default function that returns the constant life time value
      double TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag);


      //multiply the lifetime by the following constant (use it for example for adjectment to variation of a heliocentric distance)
      static const double ReactionLifetimeMultiplier=1.0;

      //the total number of the reactions that are modeled in the current model run
      static const int nTotalUnimolecularReactions=1;

      //the maximum number of the reaction products in the simulated reactions
      static const int nMaxUnimolecularReactionProducts=1;

      //the descriptor of the Unimolecular reactions
      struct cUnimoleculecularReactionDescriptor {
//      public:
//        cUnimoleculecularReactionDescriptor() {}

        int ReactionId;
        double LifeTime;
        double ReactionRate;
        int SourceSpecie;
        int nProducts;
        int ProductSpecies[nMaxUnimolecularReactionProducts];
      };

      //the array of descriptors of the reactions
      static const cUnimoleculecularReactionDescriptor UnimoleculecularReactionDescriptor[nTotalUnimolecularReactions]={0,0.0,0.0,0,0,{0}};

      //the maximum number of reactions in which a species can particiepate
      static const int nMaxSpeciesUnimolecularReactionNumber=1;

      //the list of the reactions in wich a space can particcipate
      static const int SpeciesUnimolecularReactionList[PIC::nTotalSpecies][nMaxSpeciesUnimolecularReactionNumber]={{0}};

      //the total reaction rate for individual species
      static const double TotalSpecieUnimolecularReactionRate[PIC::nTotalSpecies]={0.0};


      void InitProductStatWeight();

      //the default function for processing the photolytic reactions -> the particle is removed if the reaction occured
      inline int PhotolyticReactionProcessor_default(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {
        return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
      }

      //the manager of the photolytic reaction module
      //if the reaction has occured-> spes returns the species number of the transformed particles, TimeInterval return the time when the reaction had occured

       int PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval);
/*      inline int PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval) {
        int code=_PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
        register double p,lifetime,c;
        bool flag;

        lifetime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(x,spec,ptr,flag);
        if (flag==false) return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;


        c=exp(-TimeInterval/lifetime);
        p=1.0-c; //the probability for reaction to occur

        if (rnd()<p) {
          //determine the time offset when the reaction had occured
          TimeInterval=-lifetime*log(1.0+rnd()*(c-1.0));

          return _PHOTOLYTIC_REACTION_OCCURES_;
        }

        return code;
      }*/

    }

    namespace GenericParticleTranformation {
      //contains functions that are used to describe transformations (changing of internal parameters of a particle) that are not due to chemical reactions
      //dt <- can be limited by the function
      typedef int (*fTransformationIndicator)(double *x,double *v,int spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double &dt,bool &TransformationTimeStepLimitFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
//      extern fTransformationIndicator *TransformationIndicator;

      typedef int (*fTransformationProcessor)(double *xInit,double *xFinal,double *v,int& spec,long int ptr,PIC::ParticleBuffer::byte *ParticleData,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
//      extern fTransformationProcessor *TransformationProcessor;

      inline void Init() {
        /*
        TransformationIndicator=new fTransformationIndicator[PIC::nTotalSpecies];
        TransformationProcessor=new fTransformationProcessor[PIC::nTotalSpecies];

        //only one transformation model can be used!!!
        if (_PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_) exit(__LINE__,__FILE__,"Error: only one particle transformation model can be used");

        for (int s=0;s<PIC::nTotalSpecies;s++) TransformationIndicator[s]=NULL,TransformationProcessor[s]=NULL;
        */
      }

      inline void SetSpeciesModel(fTransformationIndicator Indicator, fTransformationProcessor Processor,int spec) {
        /*
         if (TransformationIndicator==NULL) Init();

         TransformationIndicator[spec]=Indicator;
         TransformationProcessor[spec]=Processor;
         */
      }

/*
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_       0
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_  1
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_        2
*/
    }
  }

  namespace BC {

    //the list of blocks that are connected to the bounding box, where the injection boundary conditions are applied
    extern list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > boundingBoxInjectionBlocksList;

    typedef bool (*fBlockInjectionIndicator)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionIndicator BlockInjectionBCindicatior;

    typedef long int (*fBlockInjectionBC)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionBC userDefinedBoundingBlockInjectionFunction;

    //the number of injected particles and injection rate
    extern long int nTotalInjectedParticles;
    extern long int *nInjectedParticles;
    extern double *ParticleProductionRate;

    void InitBoundingBoxInjectionBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=PIC::Mesh::mesh.rootTree);

    //model the particle injection for the current time step
    void InjectionBoundaryConditions();

    //calculate of the injection rate of particles distributed with Maxwellian distribution
    double CalculateInjectionRate_MaxwellianDistribution(const double NumberDesnity,const double Temp,const double *BulkVelocity,double *ExternalNormal,const int spec);

    //user-defined generic particle-injection function
    typedef long int (*fUserDefinedParticleInjectionFunction)();
    extern fUserDefinedParticleInjectionFunction UserDefinedParticleInjectionFunction;

    namespace InternalBoundary {

      namespace Sphere {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalSphericalData> InternalSpheres;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        //the spherical body as a internal body
        /*
        class cSurfaceDataSphere  {
        public:
          unsigned int faceat;
          double *SamplingBuffer;

          cSurfaceDataSphere() {
            faceat=-1,SamplingBuffer=NULL;
          }
        };
        */




        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalSphere();
//        cSurfaceDataSphere* GetSphereSurfaceData(cInternalBoundaryConditionsDescriptor);



        inline double* GetCompletedSamplingBuffer(cInternalSphericalData* Sphere) {return Sphere->SamplingBuffer+completedCellSampleDataPointerOffset;}
        inline double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {return GetCompletedSamplingBuffer((cInternalSphericalData*)Descriptor.BoundaryElement);}

        inline double* GetCollectingSamplingBuffer(cInternalSphericalData* Sphere) {return Sphere->SamplingBuffer+collectingCellSampleDataPointerOffset;}
        inline double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor Descriptor) {return GetCollectingSamplingBuffer((cInternalSphericalData*)Descriptor.BoundaryElement);}


        //====================================================
        //the offset of sampling data for a particular specie
        inline int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
          return 2*SurfaceElement*TotalSampleSetLength+completedCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
        }

        inline int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement) {
          return 2*SurfaceElement*TotalSampleSetLength+collectingCellSampleDataPointerOffset+SpeciesSampleDataOffset[spec];
        }




//        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
//        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
//        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
//        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data (3D)
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);


        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }

      namespace RotationBody {
      using namespace Sphere;
        extern cAMRheap<cInternalRotationBodyData> InternalRotationBody;

        cInternalBoundaryConditionsDescriptor RegisterInternalRotationBody();
        void PrintDefaultDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalRotationBodyData *RotationBody,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      }

      namespace NastranSurface {
      using namespace Sphere;
        extern cAMRheap<cInternalNastranSurfaceData> InternalNastranSurface;

        cInternalBoundaryConditionsDescriptor RegisterInternalNastranSurface();
        void PrintDefaultDataStateVector(FILE* fout,long int nElement,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalNastranSurfaceData *NastranSurface,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
      }

      namespace Circle {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalCircleData> InternalCircles;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalCircle();
        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data (3D)
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,long int nPolarPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalCircleData *Circle,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalCircleData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleCircleInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleCircleInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);


        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }

      namespace Sphere_1D {
        extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;
        extern int sampledFluxDownRelativeOffset,sampledFluxUpRelativeOffset;
        extern int sampledMeanVelocityDownRelativeOffset,sampledMeanVelocityUpRelativeOffset;
        extern int sampledMeanEnergyDownRelativeOffset,sampledMeanEnergyUpRelativeOffset;
        extern int sampledSurfaceNumberDensityRelativeOffset;

        extern long int TotalSampleSetLength;
        extern long int *SpeciesSampleDataOffset;
        extern long int *SpeciesSampleUserDefinedDataOffset;
        extern long int TotalSurfaceElementNumber;

        extern bool UserDefinedSamplingProcedureFlag;


        extern cAMRheap<cInternalSphere1DData> InternalSpheres;

        //default sampling:
        //1. particle flux Down/Up
        //2. mean particles' velocity Down/Up
        //3. mean particles' energy Down/Up
        //4. surface number density

        //the spherical body as a internal body
        /*
        class cSurfaceDataSphere  {
        public:
          unsigned int faceat;
          double *SamplingBuffer;

          cSurfaceDataSphere() {
            faceat=-1,SamplingBuffer=NULL;
          }
        };
        */




        void Init();
        void Init(long int *RequestedSamplingSetDataLength,long int *UserDefinedSampleDataRelativeOffset);

        cInternalBoundaryConditionsDescriptor RegisterInternalSphere();
//        cSurfaceDataSphere* GetSphereSurfaceData(cInternalBoundaryConditionsDescriptor);
        double* GetCompletedSamplingBuffer(cInternalBoundaryConditionsDescriptor);
        double* GetCollectingSamplingBuffer(cInternalBoundaryConditionsDescriptor);

        //the offset of the sampling data for a particular specie
        int completeSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int collectingSpecieSamplingDataOffset(int spec,long int SurfaceElement);
        int SurfaceElementSamplingSetLength();
        void switchSamplingBuffers();

        //print the 'USER DEFINED' surface data
        typedef void (*fPrintVariableList)(FILE*);
        extern fPrintVariableList PrintUserDefinedVariableList;

        typedef void (*fPrintTitle)(FILE*);
        extern fPrintTitle PrintUserDefinedTitle;

        typedef void (*fPrintDataStateVector)(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
        extern fPrintDataStateVector PrintUserDefinedDataStateVector;

        //print default surface data
        void PrintDefaultVariableList(FILE*);
        void PrintDefaultTitle(FILE*);
        void PrintDefaultDataStateVector(FILE* fout,cInternalSphere1DData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

        //clear the sampling buffers
        void flushCollectingSamplingBuffer(cInternalSphericalData* Sphere);

        //particle-spherical surface interaction
        typedef int (*fParticleSphereInteraction)(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);
        int ParticleSphereInteraction_SpecularReflection(int spec,long int ptr,double *x,double *v,double &dtTotal, void* NodeDataPointer,void *SphereDataPointer);


        //Sampling of the particles data
        typedef void (*fSampleParticleData)(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        void SampleDefaultParticleData(long int ptr,double *x,double *v,double &dtTotal, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode,cInternalBoundaryConditionsDescriptor* sphereDescriptor);
        extern fSampleParticleData SampleParticleData;
      }


    }



  }


}


#endif

//include headers for individual physical models
#if _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE_ == _PIC_MODEL__DUST__ELECTRIC_CHARGE_MODE__ON_
#include "pic__model__electrically_charged_dust.h"
#endif

/* the headers are added by the configuration scripts (exosphere.pl) during compiling of the code
//inlude headers for the user defined models
#if _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ == _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__ON_
#include "UserDefinition.PIC.PhysicalModelHeaderList.h"
#elif _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE_ == _PIC__USER_DEFINED__USER_PHYSICAL_MODEL__HEADER_LIST_MODE__OFF_
//do nothing
#else
!!!! ERROR:  The option is unknown !!!!!!
#endif
*/


