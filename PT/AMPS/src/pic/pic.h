//===============================================
//$Id$
//===============================================
//the header for the new particle model


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

//load the user defined settings
#if _PIC_USER_DEFINITION_MODE_ == _PIC_USER_DEFINITION_MODE__ENABLED_
#include "UserDefinition.PIC.h"
#endif

#include "ifileopr.h"
#include "specfunc.h"
#include "constants.h"

//include the appropriate mesh header
#if DIM == 3
#include "meshAMR3d.h"
#elif DIM == 2
#include "meshAMR2d.h"
#else 
#include "meshAMR1d.h"
#endif

#include "meshAMRinternalSurface.h"

//the maximum length of the strings
#define _MAX_STRING_LENGTH_PIC_  500


/*
#define Kbol 1.3806503E-23
#define ElectronCharge 1.602176565E-19
#define ElectronMass 9.10938291E-31
#define ProtonMass 1.67262158E-27
#define eV2J ElectronCharge
*/



namespace PIC {

  //Global constants of the PIC solver
  extern int nTotalSpecies;

  //The currect and total number of processors used in the simulation
  extern int ThisThread,nTotalThreads;


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

    //init the molecular data buffers
    void Init();

    //mass of particles
    extern double *MolMass;
    void SetMass(double,int);
    double GetMass(int);


    //get and set the value of the electric charge for a species
    extern double *ElectricCharge;
    int GetElectricCharge(int);
    void SetElectricCharge(double,int);

    //get and set the species numbers and chemical symbols
    extern char **ChemTable; // <- The table of chemical symbols used in the simulation
    extern char **LoadingSpeciesList; // <- the list of species that CAN BE locased in the simualtion


    int GetSpecieNumber(char*);
    void SetChemSymbol(char*,int);
    void GetChemSymbol(char*,int);

    //set and get the specie type (gas, external, background)
    extern int *SpcecieTypeTable;
    void SetSpecieType(int SpcecieType,int spec);
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
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      //do nothing
    #else
      exit(__LINE__,__FILE__,"Error: unknown option");
    #endif
    }

    inline void SetIndividualStatWeightCorrection(double WeightCorrectionFactor,byte *ParticleDataStart) {
    #if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
      *((double*) (ParticleDataStart+_PIC_PARTICLE_DATA_WEIGHT_CORRECTION_OFFSET_)) =WeightCorrectionFactor;
    #elif _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_OFF_
      //do nothing
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

      memcpy(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,v,3*sizeof(double));
    }

    inline void SetV(double* v,byte *ParticleDataStart) {
      /*
      register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
      */

      memcpy(ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_,v,3*sizeof(double));
    }


    //==========================================================
    //get the particle's species ID

    inline unsigned int GetI(byte* ParticleDataStart) {
      return *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
    }

    inline unsigned int GetI(long int ptr) {
      return *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
    }

    inline void SetI(int spec,byte* ParticleDataStart) {
      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=(unsigned char)spec;
    }

    inline void SetI(int spec,long int ptr) {
      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=(unsigned char)spec;
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

      long int FirstCellParticle,tempParticleMovingList;

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

	      FirstCellParticle=-1,tempParticleMovingList=-1;

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

      void GetBulkVelocity(double *v,int s) {
        int idim;
        double TotalWeight,*SampledData;

        #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
        if ((s<0)||(s>=PIC::nTotalSpecies)) exit(__LINE__,__FILE__,"Error: 's' is out of the range");
        #endif

        if (PIC::LastSampleLength!=0) {
          TotalWeight=(*(s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleWeghtRelativeOffset)));
          SampledData=3*s+(double*)(associatedDataPointer+PIC::Mesh::completedCellSampleDataPointerOffset+PIC::Mesh::sampledParticleVelocityRelativeOffset);

          if (TotalWeight>0.0) for (idim=0;idim<DIM;idim++) v[idim]=SampledData[idim]/TotalWeight/PIC::LastSampleLength;
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

          if (TotalWeight>0.0) for (idim=0;idim<DIM;idim++) v2[idim]=SampledData[idim]/TotalWeight/PIC::LastSampleLength;
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

          if (TotalWeight>0.0) res=(*SampledData)/TotalWeight/PIC::LastSampleLength;
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
        cBasicBlockAMR<cDataCornerNode,cDataCenterNode>cleanDataBuffer();
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
    double GetBlockInjectionRate(int spec,PIC::Mesh::cDataBlockAMR *block);
    double GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell);
  }
#endif

  //sampling functions
  namespace Sampling {

    //the minimum number of iteration for output of the datafile
    extern int minIterationNumberForDataOutput;

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

    //collisions with the background atmosphere
#if _PIC_BACKGROUND_ATMOSPHERE_MODE_ == _PIC_BACKGROUND_ATMOSPHERE_MODE__ON_
    namespace BackgroundAtmosphere {

      //include the user defined properties of the background atmosphere
      #include "UserDefinition.BackgroundAtmosphere.h"


      //Sampling of the model data
      extern int LocalTotalCollisionFreqSamplingOffset;
      int RequestSamplingData(int offset);
      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);


      //init the model
      void Init_BeforeParser();



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

    void SampleDistributionFnction(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* sampleNode);
    void flushSamplingBuffers();

    void printDistributionFunction(char *fname,int spec);
  }

  //procedures for distribution of particle velocities
  namespace Distribution {

    //the macroses defined the modes for generation of the particle parameters
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_ -> the particle properties are generated with the requested distributuion, no weight correction are needed
    //_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_ -> to recover the requasted distribution, a particle weght correction is needed

    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_         0
    #define _PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_   1


    void MaxwellianVelocityDistribution(double *v,double *BulkFlowVelocity,double Temp,int spec);
    double InjectMaxwellianDistribution(double *v,double *BulkFlowVelocity,double Temp,double *ExternalNormal,int spec,int WeightCorrectionMode=_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__NO_WEIGHT_CORRECTION_);
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

  namespace Mover {

    #include "UserDefinition.PIC.Mover.h"

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
    typedef int (*fProcessOutsideDomainParticles) (long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
    extern fProcessOutsideDomainParticles ProcessOutsideDomainParticles;


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
      printf("!!!!! Execution is finished by the alarm (PIC::Alarm) !!!!!!\n");
      exit(__LINE__,__FILE__,"!!!!! Execution is finished by the alarm !!!!!!");
    }
  }

  namespace ICES {

     extern char locationICES[_MAX_STRING_LENGTH_PIC_]; //location of the data and the dace cases

     void Init();
     void SetLocationICES(const char*);

     //the total number of bytes used to store the ICES data vector; the offset of the data associated with the ICES data vector
     extern int TotalAssociatedDataLength,AssociatedDataOffset;

     //the offsets for the plasma parameters loaded with ICES
     extern int ElectricFieldOffset,MagneticFieldOffset,PlasmaPressureOffset,PlasmaNumberDensityOffset,PlasmaTemperatureOffset,PlasmaBulkVelocityOffset;

     //the offsets for parameters loaded from the DSMC model
     extern int NeutralBullVelocityOffset,NeutralNumberDensityOffset,NeutralTemperatureOffset;

     //calcualte the total number of cells in the mesh
     long int getTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

     //create the trajectory file
     void createCellCenterCoordinateList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=PIC::Mesh::mesh.rootTree);


     //retrive the SWMF data file
     struct cDataNodeSWMF {
       double swNumberDensity,swTemperature,swPressure,E[3],B[3],swVel[3];
     };

     struct cDataNodeDSMC {
       double neutralNumberDensity,neutralTemperature,neutralVel[3];
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

  namespace ChemicalReactions {

    namespace PhotolyticReactions {
      //the photolytic reactions: the model is initialized if _PIC_PHOTOLYTIC_REACTIONS_MODE_ = _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
      //the total photolytic lifetime of the species;
      extern double *ConstantTotalLifeTime;

      typedef double (*fTotalLifeTime)(double *x,int spec,long int ptr,bool &ReactionAllowedFlag);
      extern fTotalLifeTime *TotalLifeTime;

      typedef int (*fReactionProcessor)(double *xInit,double *xFinal,long int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData);
      extern fReactionProcessor *ReactionProcessorTable;

      void Init();
      void SetReactionProcessor(fReactionProcessor f,int spec);
      void SetSpeciesTotalPhotolyticLifeTime(fTotalLifeTime f,int spec);

      //The return codes of the photolytic model
      #define _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_      0
      #define _PHOTOLYTIC_REACTION_OCCURES_                  1

      #define _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_        2
      #define _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_ 3

      //the default function that returns the constant life time value
      double TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag);

      //the manager of the photolytic reaction module
      //if the reaction has occured-> spes returns the species number of the transformed particles, TimeInterval return the time when the reaction had occured

      inline int PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval) {
        int code=_PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
        register double p,lifetime,c;
        bool flag;

        lifetime=TotalLifeTime[spec](x,spec,ptr,flag);
        if (flag==false) return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;


        c=exp(-TimeInterval/lifetime);
        p=1.0-c; //the probability for reaction to occur

        if (rnd()<p) {
          //determine the time offset when the reaction had occured
          TimeInterval=-lifetime*log(1.0+rnd()*(c-1.0));

          return _PHOTOLYTIC_REACTION_OCCURES_;
        }

        return code;
      }

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

      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_       0
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__TRANSFORMATION_OCCURED_  1
      #define _GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_        2
    }
  }

  namespace BC {

    //the list of blocks that are connected to the bounding box, where the injection boundary conditions are applied
    extern list<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* > boundingBoxInjectionBlocksList;

    typedef bool (*fBlockInjectionIndicator)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionIndicator BlockInjectionBCindicatior;

    typedef long int (*fBlockInjectionBC)(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
    extern fBlockInjectionBC userDefinedBoundingBlockInjectionFunction;

    //the numbe of injected particles
    extern long int nInjectedParticles;

    void InitBoundingBoxInjectionBlockList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode=PIC::Mesh::mesh.rootTree);

    //model the particle injection for the current time step
    void InjectionBoundaryConditions();

    //calculate of the injection rate of particles distributed with Maxwellian distribution
    double CalculateInjectionRate_MaxwellianDistribution(double NumberDesnity,double Temp,double *BulkVelocity,double *ExternalNormal,int spec);


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




