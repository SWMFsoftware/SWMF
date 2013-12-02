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

//SPICE ephemeris time
SpiceDouble Europa::OrbitalMotion::et=0.0,Europa::OrbitalMotion::lt=0.0;

//direction to the Sun and the angle of the rotation between planetary axes and the direction to the Sun on the Z-plane
double Europa::OrbitalMotion::SunDirection_IAU_EUROPA[3]={0.0,0.0,0.0},Europa::OrbitalMotion::PlanetAxisToSunRotationAngle=0.0;

//matrixes for transformation GALL_EPHIOD->IAU and IAU->GALL_EPHIOD coordinate frames
SpiceDouble Europa::OrbitalMotion::GALL_EPHIOD_to_IAU_TransformationMartix[6][6],Europa::OrbitalMotion::IAU_to_GALL_EPHIOD_TransformationMartix[6][6];

// Transformation matrix from IAU to local coordinate system with respect to the Sun
double Europa::OrbitalMotion::IAU_to_LS_TransformationMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

//transformation matrix: IAU -> GALL_ESOM (used for calcualtion of the surface temeprature)
SpiceDouble Europa::OrbitalMotion::IAU_to_GALL_ESOM_TransformationMatrix[6][6];
