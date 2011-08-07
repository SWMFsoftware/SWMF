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

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

#ifndef _PIC_
#define _PIC_

#include "picGlobal.dfn"
#include "ifileopr.h"
#include "specfunc.h"

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
#define Kbol 1.3806503E-23




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


  //perform one time step
  void TimeStep();
  void Sampling();

  //init the particle solver
  void Init();

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
    extern char **ChemTable;
    int GetSpecieNumber(char*);
    void SetChemSymbol(char*,int);
    void GetChemSymbol(char*,int);

    namespace Parser {
      void InitChemTable(CiFileOperations&);
      int GetSpeciesNumber(int&,CiFileOperations&);
      void SpeciesBlock(int,CiFileOperations&);
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
      register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
    }

    inline void GetX(double* x,byte *ParticleDataStart) {
      register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
    }

    inline void SetX(double* x,long int ptr) {
      register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];
    }

    inline void SetX(double* x,byte *ParticleDataStart) {
      register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];
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
      register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
    }

    inline void GetV(double* v,byte *ParticleDataStart) {
      register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
    }

    inline void SetV(double* v,long int ptr) {
      register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
    }

    inline void SetV(double* v,byte *ParticleDataStart) {
      register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
      register int idim;

      for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
    }


    //==========================================================
    //get the particle's species ID

    inline unsigned int GetI(byte* ParticleDataStart) {
      return *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
    }

    inline unsigned int GetI(long int ptr) {
      return *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
    }

    inline void SetI(unsigned int spec,byte* ParticleDataStart) {
      *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=spec;
    }

    inline void SetI(unsigned int spec,long int ptr) {
      *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=spec;
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
	  //the limiting size of the domain and the function controlling the local mesh resolution
	  extern double xmin[3],xmax[3];

	  typedef double (*fLocalMeshResolution) (double*);
	  extern fLocalMeshResolution LocalMeshResolution;

	  //the offset of the sampled infomation that is stored in 'center nodes'
    extern int completedCellSampleDataPointerOffset,collectingCellSampleDataPointerOffset;

    //the flag determines if 'external' sampling procedures are used
    extern bool ExternalSamplingProcedureDefinedFlag;

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




    //the class defining the 'central node' that contains the sampling data
    class cDataCenterNode : public cBasicCenterNode {
    public:
	    //parameters that defines the parameters of the associated data used for sampling and code running
      static int totalAssociatedDataLength;

      long int FirstCellParticle,tempParticleMovingList;

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

	    //clean the sampling buffers
	    void cleanDataBuffer() {
	      cBasicCenterNode::cleanDataBuffer();

	      FirstCellParticle=-1,tempParticleMovingList=-1;

	      int i,length=totalAssociatedDataLength/sizeof(double);
	      double *ptr;
	      for (i=0,ptr=(double*)associatedDataPointer;i<length;i++,ptr++) *ptr=0.0;
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
      }

      void PrintFileDescriptior(FILE* fout,int DataSetNumber) {
        char sym[_MAX_STRING_LENGTH_PIC_];

        PIC::MolecularData::GetChemSymbol(sym,DataSetNumber);
        fprintf(fout,"TITLE=\"specie=%s\"",sym);
      }

      void PrintVariableList(FILE* fout,int DataSetNumber) {
       fprintf(fout,", \"Number Density\", \"Particle Number\"");
       for (int idim=0;idim<DIM;idim++) fprintf(fout,", \"V%i\"",idim);
       fprintf(fout,", \"Speed\", \"Translational Temperature\"");
      }

      void Interpolate(cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients) {
        int i,s,idim;
        double c;

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
      }
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


//      void SendNodeData(CMPI_channel *pipe,int DataSetTag);
      /*
      void SendNodeData(CMPI_channel *pipe,int DataSetTag) {
        int iCell,jCell,kCell;
        long int LocalCellNumber;
        PIC::Mesh::cDataCenterNode *cell=NULL;

        #if DIM == 3
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
        #elif DIM == 2
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
        #elif DIM == 1
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
        #else
        exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
        #endif

        if ((DataSetTag==_PIC_SUBDOMAIN_BOUNDARY_LAYER_SAMPLING_DATA_EXCHANGE_TAG_)||(DataSetTag==_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_)) {
          for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
            LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
            cell=GetCenterNode(LocalCellNumber);

            pipe->send(cell->associatedDataPointer,cell->totalAssociatedDataLength);
          }

        }
        else exit(__LINE__,__FILE__,"Error: unknown option");

        //send all blocks' data when the blocks is moved to another processor
        if (DataSetTag==_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_) {
          int Signal;

          const int _CENTRAL_NODE_NUMBER_SIGNAL_=1;
          const int _NEW_PARTICLE_SIGNAL_=       2;
          const int _END_COMMUNICATION_SIGNAL_=  3;

          for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(iCell,jCell,kCell);
            Particle=sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle;

            if  (Particle!=-1) {
              pipe.send(_CENTRAL_NODE_NUMBER_SIGNAL_);
              pipe.send(LocalCellNumber);

              while (Particle!=-1) {
                PIC::ParticleBuffer::PackParticleData(buffer,Particle);
                pipe.send(_NEW_PARTICLE_SIGNAL_);
                pipe.send(buffer,PIC::ParticleBuffer::ParticleDataLength);

                NextParticle=PIC::ParticleBuffer::GetNext(Particle);
                PIC::ParticleBuffer::DeleteParticle(Particle);
                Particle=NextParticle;
              }

              sendNode->block->GetCenterNode(LocalCellNumber)->FirstCellParticle=-1;
            }
          }

          pipe.send(_END_COMMUNICATION_SIGNAL_);
        }


      }
      */

//      void RecvNodeData(CMPI_channel *pipe,int DataSetTag,int From);
      /*
      void RecvNodeData(CMPI_channel *pipe,int DataSetTag,int From) {
        int iCell,jCell,kCell;
        long int LocalCellNumber;
        PIC::Mesh::cDataCenterNode *cell=NULL;

        #if DIM == 3
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=_BLOCK_CELLS_Z_;
        #elif DIM == 2
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=_BLOCK_CELLS_Y_,kCellMax=1;
        #elif DIM == 1
        static const int iCellMax=_BLOCK_CELLS_X_,jCellMax=1,kCellMax=1;
        #else
        exit(__LINE__,__FILE__,"Error: the value of the parameter is not recognized");
        #endif

        if (DataSetTag==_PIC_SUBDOMAIN_BOUNDARY_LAYER_SAMPLING_DATA_EXCHANGE_TAG_) {
          for (kCell=0;kCell<kCellMax;kCell++) for (jCell=0;jCell<jCellMax;jCell++) for (iCell=0;iCell<iCellMax;iCell++) {
            LocalCellNumber=getCenterNodeLocalNumber(iCell,jCell,kCell);
            cell=GetCenterNode(LocalCellNumber);

            pipe->recv(cell->associatedDataPointer,cell->totalAssociatedDataLength,From);
          }

        }
        else exit(__LINE__,__FILE__,"Error: unknown option");
      }
      */

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
    double GetLocalTimeStep(int,cDataBlockAMR*);
    void SetLocalTimeStep(double,int,cDataBlockAMR*);

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
    exit(__LINE__,__FILE__,"Not implemented yet");
    #endif



    //init the computational mesh
    void Init(double*,double*,fLocalMeshResolution);
    void buildMesh();
    void loadMesh(char*);

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
    void MaxwellianVelocityDistribution(double *v,double *BulkFlowVelocity,double Temp,int spec);
    void InjectMaxwellianDistribution(double *v,double *BulkFlowVelocity,double Temp,double *ExternalNormal,int spec);
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
    //the return codes of the moving procedures
    #define _PARTICLE_REJECTED_ON_THE_FACE_ -1
    #define _PARTICLE_DELETED_ON_THE_FACE_   0
    #define _PARTICLE_CROSS_THE_FACE_        1
    #define _PARTICLE_LEFT_THE_DOMAIN_       2
    #define _PARTICLE_MOTION_FINISHED_       3

    #include "ParticleAcceleration.PIC.h"

    typedef int (*fSpeciesDependentParticleMover) (long int,double,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);

    //the vector containing the species specific particle moving procedures
    extern fSpeciesDependentParticleMover *MoveParticleTimeStep;

    void Init();
    void MoveParticles();

    int UniformWeight_UniformTimeStep_noForce(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
    int UniformWeight_UniformTimeStep_noForce_TraceTrajectory(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

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

  }

  namespace Parallel {


     //count the number of particles that were send and recieve by the thread
     extern long int sendParticleCounter,recvParticleCounter;

     //exchenge paricles between iterations
     void ExchangeParticleData();
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

    }



  }


}



#endif
