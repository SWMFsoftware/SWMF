/*
 * ParticleTracker.h
 *
 *  Created on: Nov 21, 2014
 *      Author: vtenishe
 */

//$Id$

#ifndef PARTICLETRACKER_H_
#define PARTICLETRACKER_H_

#include "pic.h"

void TotalParticleAcceleration(double *accl,int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode);
bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData);

#endif /* PARTICLETRACKER_H_ */
