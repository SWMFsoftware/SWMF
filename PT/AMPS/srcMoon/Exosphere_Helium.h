/*
 * Mercury.h
 *
 *  Created on: Feb 10, 2012
 *      Author: vtenishe
 */

//$Id$


#ifndef _EXOSPHERE__HELIUM_
#define _EXOSPHERE__HELIUM_

#include "He.h"
#include "Exosphere.h"
#include "UserDefinition.Exosphere.h"
#include "Moon.h"





      namespace HeliumDesorption {
      using namespace Moon;
        //typical solar wind conditions far from the planet

        extern double SourceRate[PIC::nTotalSpecies],maxLocalSourceRate[PIC::nTotalSpecies];
        //5% of solar wind ions are assumed Helium e.g., LeBlanc 2011 Bochsler, 2001
        const static double sw_Helium_Ion_Flux_1AU = -Exosphere_swVelocity_Typical[1]*Exosphere_swNumberDensity_Typical*0.05*Pi*_MOON__RADIUS_*_MOON__RADIUS_;

        //the object for distribution of injection positino on the planet's surface
        extern cSingleVariableDiscreteDistribution<int> SurfaceInjectionDistribution[PIC::nTotalSpecies];
        double GetSurfaceElementProductionRate(int nElement,int *spec);

        inline double GetTotalProductionRate(int spec,void *SphereDataPointer) {return SourceRate[spec];}

        double GetSurfaceElementProductionRate(int spec,int SurfaceElement,void *SphereDataPointer);


               bool GenerateParticleProperties(int spec,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double* v_IAU_OBJECT, double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, cInternalSphericalData* Sphere);


      }

#endif /* OBJECT_H_ */
