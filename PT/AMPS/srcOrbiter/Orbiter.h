/*
 * Orbiter.h
 *
 *  Created on: May 19, 2014
 *      Author: vtenishe
 */
//$Id$


#ifndef _ORBITER_H_
#define _ORBITER_H_

#include "Exosphere.h"
#include "Orbiter.dfn"
#include "SingleVariableDiscreteDistribution.h"

namespace Orbiter {
using namespace Exosphere;

  //Exchange model data
  void ExchangeModelData();

  //init the Orbiter model
  void Init_BeforeParser();

  //injection models
  namespace InjectionModel {

    //source from a group of faces that shares the same face attribute
    namespace FaceEjection {

      struct cFaceTableEntry {
        double TotalArea;
        int *FaceTable;
        int FaceTableLength;
        int faceat;
        cSingleVariableDiscreteDistribution<int> InjectionFaceGenerator;
      };

      struct cInjectionData {
        int faceat;
        double SourceRate;
        double Temperature;
        int Species;
        cFaceTableEntry *FaceTable;
      };

      //this source is considered active only when InjectionDataTableLength is not zero, which is a default value
      extern cFaceTableEntry *InjectionFaceTable;
      extern int InjectionFaceTableLength;

      extern cInjectionData InjectionDataTable[];
      extern int InjectionDataTableLength;

      //call the init function when needed
      extern bool SourceInitFlag;
      void Init();

      //inject model particles
      long int InjectParticles();
    }

    //point source of the model particles
    namespace PointSource {

      //source rate of the injeted particles
      struct cInjectionData {
        double SourceRate;
        double Temperature;
        double Location[3];
        int Species;
      };

      //this source is considered active only when InjectionDataTableLength is not zero, which is a default value
      extern cInjectionData InjectionDataTable[];
      extern int InjectionDataTableLength;

      //inject particles functions
      long int  InjectParticles();
    }

    //particle injection from a ring
    namespace Ring {
      //source rate, and temperature of the injected species
      extern double SourceRate[PIC::nTotalSpecies];
      extern double SourceTemperature[PIC::nTotalSpecies];

      //location of the center of the ring, and the coordinate frame related to the ring
      extern double e0[3],e1[3],e2[3],x0[3],Radius;

      //inject model particles
      long int  InjectParticles();
    }

    //particle desorption
    namespace Desorption {

      //inject model particles
      long int  InjectParticles();
    }

    //diffusion of particles from the surface
    namespace Diffusion {

      //inject model particles
      long int InjectParticles();
    }

    //the function controls injecting of the model particles from all sources
    long int  InjectParticles();
    double GetTotalInjectionRate(int spec);
  }

  //adsorption model
  namespace Adsorption {

    double GetStickingCoefficient();
  }

  //upstream boundary conditions
  namespace UpstreamBC {
    extern double Velocity[3];
    extern double NumberDensity[PIC::nTotalSpecies];
    extern double Temperature;
    extern bool UpstreamSourceMode;

    bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
    double BoundingBoxInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode);
  }

  //mesh data
  namespace Mesh {
    extern char sign[_MAX_STRING_LENGTH_PIC_];
  }

  //domain size limits
  namespace DomainSize {
    extern double xMinOffset[3],xMaxOffset[3];
  }

  //the set of function needed for calculating of the projection area of the spacecraft
  extern double ProjectionOrbiterSurfaceArea;
  double CalculateProjectionArea();

  //the name of the surace mesh file
  namespace SurfaceModel {
    //constant defining the mesh file formats
    const int MeshFileFormat_NASTRAN=0;
    const int MeshFileFormat_CEA=1;

    extern char MeshFileName[_MAX_STRING_LENGTH_PIC_];
    extern int MeshFileFormat;

    //scale the size of the surface mesh
    extern double ScalingFactor;
  }

  //sample calulated aerodynamcis properties
  namespace Sampling {

     //sampling of the drag coefficient
     namespace DragCoefficient {
       extern bool SamplingMode;
       extern double *dpX,*dpY,*dpZ;

       void Init();

       void SamplingProcessor();
       void PrintOutputFile(int);
     }

  }

  double GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement);
  double GetTotalProduction(int spec,void *BoundaryElement);

  bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
  //(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, char *tempParticleData,int BoundaryElementType,void *BoundaryElement);
}



#endif /* _ORBITER_H_ */


