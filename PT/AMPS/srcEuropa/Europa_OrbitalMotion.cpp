/*
 * Europa_OrbitalMotion.cpp
 *
 *  Created on: Mar 24, 2012
 *      Author: vtenishe
 */

//Orbital motion of Europa
//$Id$


#include "pic.h"

double Europa::OrbitalMotion::AccumulatedPlanetRotation=0.0,Europa::OrbitalMotion::TotalSimulationTime=0.0,Europa::OrbitalMotion::TAA=0.0;
double Europa::OrbitalMotion::CoordinateFrameRotationRate=0.0;

//direction to the Sun and the angle of the rotation between planetary axes and the direction to the Sun on the Z-plane
double Europa::OrbitalMotion::SunDirection_IAU_EUROPA[3]={0.0,0.0,0.0},Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=0.0;

// Transformation matrix from IAU to local coordinate system with respect to the Sun
double Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

//transformation matrix: IAU -> GALL_ESOM (used for calcualtion of the surface temeprature)
SpiceDouble Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[6][6];

//The locations of the Sun nad Jupiter in the frame wherethe simulation is performed (cenetered on Europa)
double Europa::OrbitalMotion::xJupiter_SO[3]={0.0,0.0,0.0}; ///,Europa::OrbitalMotion::xSun_SO[3]={0.0,0.0,0.0};


//update location of the Sun and Jupiter
void Europa::OrbitalMotion::UpdateSunJupiterLocation() {
  SpiceDouble StateSun[6],StateJupiter[6],ltlocal;
  int idim;

  spkezr_c("Sun",et,Europa::SO_FRAME,"none","EUROPA",StateSun,&ltlocal);
  spkezr_c("Jupiter",et,Europa::SO_FRAME,"none","EUROPA",StateJupiter,&ltlocal);

  for (idim=0;idim<3;idim++) xJupiter_SO[idim]=1.0E3*StateJupiter[idim],Europa::xSun_SO[idim]=1.0E3*StateSun[idim];
}
