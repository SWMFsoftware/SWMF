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

namespace Orbiter {
using namespace Exosphere;

  //init the Orbiter model
  void Init_BeforeParser();

  //upstream boundary conditions
  namespace UpstreamBC {
    extern double Velocity[3];
    extern double NumberDensity[PIC::nTotalSpecies];
    extern double Temperature;

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


