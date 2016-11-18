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

  double GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement);
  double GetTotalProduction(int spec,void *BoundaryElement);

  bool GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, int BoundaryElementType,void *BoundaryElement);
  //(int spec, double *x_SO_OBJECT,double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0, double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode, char *tempParticleData,int BoundaryElementType,void *BoundaryElement);
}



#endif /* _ORBITER_H_ */


