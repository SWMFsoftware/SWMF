//$Id$

#ifndef _PIC_EARTH_H_
#define _PIC_EARTH_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


#include "Exosphere.h"

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif


#include "constants.h"
#include "Earth.dfn"

#include "GCR_Badavi2011ASR.h"

//class that is used for keeping information of the injected faces
class cBoundaryFaceDescriptor {
public:
  int nface;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;

  cBoundaryFaceDescriptor() {
    nface=-1,node=NULL;
  }
};

class cCompositionGroupTable {
public:
  int FistGroupSpeciesNumber;
  int nModelSpeciesGroup;
  double minVelocity,maxVelocity; //the velocity range of particles from a given species group that corresponds to the energy range from Earth::BoundingBoxInjection
  double GroupVelocityStep;   //the velocity threhold after which the species number in the group is switched
  double maxEnergySpectrumValue;

  double inline GetMaxVelocity(int spec) {
    int nGroup=spec-FistGroupSpeciesNumber;

    if ((nGroup<0)||(spec>=FistGroupSpeciesNumber+nModelSpeciesGroup)) exit(__LINE__,__FILE__,"Error: cannot recogniza the velocit group");
    return minVelocity+(nGroup+1)*GroupVelocityStep;
  }

  double inline GetMinVelocity(int spec) {
    int nGroup=spec-FistGroupSpeciesNumber;

    if ((nGroup<0)||(spec>=FistGroupSpeciesNumber+nModelSpeciesGroup)) exit(__LINE__,__FILE__,"Error: cannot recogniza the velocit group");
    return minVelocity+nGroup*GroupVelocityStep;
  }

  cCompositionGroupTable() {
    FistGroupSpeciesNumber=-1,nModelSpeciesGroup=-1;
    minVelocity=-1.0,maxVelocity=-1.0,GroupVelocityStep=-1.0;
  }
};

namespace Earth {
  using namespace Exosphere;
  
  //that fucntion that created the SampledDataRecoveryTable
  void DataRecoveryManager(list<pair<string,list<int> > >&,int,int);

  //composition table of the GCR composition
  extern cCompositionGroupTable *CompositionGroupTable;
  extern int *CompositionGroupTableIndex;
  extern int nCompositionGroups;

  //the mesh parameters
  namespace Mesh {
    extern char sign[_MAX_STRING_LENGTH_PIC_];
  }
  
  //sampling of the energetic particle fluency
  namespace Sampling {

    //the flag determine whether sampling will be performed
    extern bool SamplingMode;

    //the number of the spherical layers of the sampling spheres, and shell's radii
    const int nSphericalShells=1;
    extern double SampleSphereRadii[nSphericalShells];


    //sampling of the particle fluency
    namespace Fluency {
      extern double minSampledEnergy;
      extern double maxSampledEnergy;

      extern int nSampledLevels;
      extern double dLogEnergy;
    }

    //Sphere object used for sampling
    extern cInternalSphericalData SamplingSphericlaShell[nSphericalShells];

    //output sampled data
    void PrintVariableList(FILE*);
    void PrintTitle(FILE*);
    void PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);

    //the function that will be called to print sampled data
    void PrintManager(int);
    void SamplingManager();

   //Init sampling procedure
    void Init();

    //Sample GyroFrequency and GyroRadius of the model particles
    namespace ParticleData {
      extern int _GYRO_RADIUS_SAMPLING_OFFSET_;
      extern int _GYRO_FRECUENCY_SAMPLING_OFFSET_;
      extern bool SamplingMode;

      int RequestSamplingData(int offset);
      void SampleParticleData(char* tempParticleData,double LocalParticleWeight,char *SamplingData,int s,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);

      void PrintVariableList(FILE* fout,int DataSetNumber);
      void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
      void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);


      void Init();
    }
  }

  //particle mover: call Boris-relativistic, and sample particles that intersects the sampling shells
  int ParticleMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);



  //calculation of the cutoff rigidity
  namespace CutoffRigidity {
    const double RigidityTestRadiusVector=400.0E3+_RADIUS_(_TARGET_);
    const double RigidityTestMinEnergy=1.0*MeV2J;
    const double RigidityTestMaxEnergy=1.0E4*MeV2J;

    namespace OutputDataFile {
      void PrintVariableList(FILE* fout);
      void PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,
          long int SurfaceElementsInterpolationListLength,
          cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
    }

    namespace ParticleInjector {
      double GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer);
      bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
          double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
    }

    //save the location of the particle origin, and the particle rigidity
    extern long int InitialRigidityOffset,InitialLocationOffset;

    //sample the rigidity mode
    extern bool SampleRigidityMode;

    //sphere for sampling of the cutoff regidity
    extern double ***CutoffRigidityTable;

    //process particles that leaves that computational domain
    int ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);

    //reversed time particle integration procedure
    int ReversedTimeRelativisticBoris(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

    //init the cutoff sampling model
    void AllocateCutoffRigidityTable();
    void Init_BeforeParser();
  }

  //impulse source of the energetic particles
  namespace ImpulseSource {
    extern bool Mode;  //the model of using the model
    extern double TimeCounter;

    struct cImpulseSourceData {
      int spec;
      double time;
      bool ProcessedFlag;
      double x[3];
      double Source;
      double NumberInjectedParticles;
    };

    namespace EnergySpectrum {
      const int Mode_Constatant=0;

      extern int Mode;

      namespace Constant {
        extern double e;
      }
    }

    extern cImpulseSourceData ImpulseSourceData[];
    extern int nTotalSourceLocations;

    long int InjectParticles();
    void InitParticleWeight();
  }

  //injection of new particles into the system
  namespace BoundingBoxInjection {
    //Energy limis of the injected particles
    const double minEnergy=1.0*MeV2J;
    const double maxEnergy=1.0E4*MeV2J;

    //the list of the faces located on the domain boundary through which particles can be injected
    extern cBoundaryFaceDescriptor *BoundaryFaceDescriptor;
    extern int nTotalBoundaryInjectionFaces;

/*
    //init and populate the tables used by the particle injectino procedure
    void InitBoundingBoxInjectionTable(double** &BoundaryFaceLocalInjectionRate,double* &maxBoundaryFaceLocalInjectionRate,
        double* &BoundaryFaceTotalInjectionRate,double (*InjectionRate)(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode));
*/

    //model that specify injectino of the gakactic cosmic rays
    namespace GCR {
/*      //data buffers used for injection of GCR thought the boundary of the computational domain
      extern double **BoundaryFaceLocalInjectionRate;
      extern double *maxBoundaryFaceLocalInjectionRate;
      extern double *BoundaryFaceTotalInjectionRate;*/


      extern double *InjectionRateTable;  //injection rate for a given species/composition group component
      extern double *maxEnergySpectrumValue;  //the maximum value of the energy epectra for a given species/composition group component

      //init the model
      void Init();

      //init and populate the tables used by the particle injectino procedure
      void InitBoundingBoxInjectionTable();

      //source rate and geenration of GCRs
      double InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
      void GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double *x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal);

      //init the model
      void Init();

    }

    //model that specifies injectino of SEP
    namespace SEP {

      //source rate and generation of SEP
      double InjectionRate(int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
      void GetNewParticle(PIC::ParticleBuffer::byte *ParticleData,double *x,int spec,int nface,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode,double *ExternalNormal);

      //init the model
      void Init();
    }

    //general injection functions
    bool InjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    long int InjectionProcessor(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);

    long int InjectionProcessor(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    long int InjectionProcessor(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    double InjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);


  }

  //interaction of the particles with the surface of the Earth
  namespace BC {
    int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v, double &dtTotal, void *NodeDataPonter,void *SphereDataPointer);
    double sphereInjectionRate(int spec,void *SphereDataPointer);
  }

  //Init the model
  void Init();

}

#endif /* EARTH_H_ */
