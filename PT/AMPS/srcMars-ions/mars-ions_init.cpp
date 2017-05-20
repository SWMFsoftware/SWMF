
//$Id$
//functions for calculation of the time step and mesh resolution

/*
 * mars-ion__TimeStep_MeshResolution.cpp
 *
 *  Created on: May 26, 2015
 *      Author: vtenishe
 */


#include "pic.h"
#include "mars-ions.h"


//the mesh resolution at the lower boundary
double localSphericalSurfaceResolution(double *x) {
  double res;

  switch (_PIC_NIGHTLY_TEST_MODE_) {
  case _PIC_MODE_ON_: 
    res=_RADIUS_(_TARGET_)/20.0;
    break;
  default:
     res=_RADIUS_(_TARGET_)/150.0;
  }

  return res;
}

//the mesh resulution within the domain
double localResolution(double *x) {
  return 5.0*_RADIUS_(_TARGET_)/20.0;
}

//the local time step
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize=startNode->GetCharacteristicCellSize();
  double CharacteristicSpeed=1.0E5;

  switch (spec) {
  case _H_PLUS_SPEC_:
    CharacteristicSpeed=1.0E6;
    break;
  default:
    CharacteristicSpeed=1.0E5;
  }

  return 0.3*CellSize/CharacteristicSpeed;
}
