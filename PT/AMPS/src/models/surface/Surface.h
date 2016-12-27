/*
 * Surface.h
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 */

#ifndef _MODELS_SURFACE_SURFACE_H_
#define _MODELS_SURFACE_SURFACE_H_

#include "Surface.dfn"
#include "pic.h"

//the model of the gas-surface interaction
namespace Surface {

   //the models of the surface temeprature
   namespace Temeprature {

     //isothermal
     namespace Isothremal {
       extern double Temp;
     }

     //Radiative Equilibrium
     namespace RadiativeEquilibrium {

     }
   }

   double GetSurfaceTemeprature(CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

   //diffuse refrelction
   namespace DiffuseReflection {
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //specular refrelction
   namespace SpecularReflection {
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //quasi-specular reflection
   namespace QuasiSpecularReflection {
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //Cercignani-Lampis-Lord model
   namespace CLL {
     extern double TangentialMomentumAccommodationCoefficient[PIC::nTotalSpecies];
     extern double NormalPartKineticEnergyAccommodationCoefficient[PIC::nTotalSpecies];

     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //Maxwell model
   namespace MaxwellReflection {
     extern double AccommodationCoefficient[PIC::nTotalSpecies];

     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //process gas/surface interaction event
   int ParticleInteractionProcessor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
}



#endif /* _MODELS_SURFACE_SURFACE_H_ */
