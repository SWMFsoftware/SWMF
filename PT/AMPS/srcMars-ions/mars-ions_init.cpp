
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
  return _RADIUS_(_TARGET_)/20.0;
}

//the mesh resulution within the domain
double localResolution(double *x) {
  return 5.0*_RADIUS_(_TARGET_)/20.0;
}

//the local time step
double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize=startNode->GetCharacteristicCellSize();
  double CharacteristicSpeed=1.0E5;

  return 0.3*CellSize/CharacteristicSpeed;
}
