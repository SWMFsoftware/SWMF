//$Id$
//Evaluate density and pressure along the spacecraft trajectory using the Liouville theoreme

/*
 * RosinaMeasurements_Liouville.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "RosinaMeasurements.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif


void RosinaSample::Liouville::EvaluateLocation(int iPoint) {
  int iTest,iSurfaceElement;
  double c,l[3],xLocation[3],rLocation,xIntersection[3];


  const int nTotalTests=10000000;
  double cosAngleThrouhold=cos(25.0/180.0*Pi);

  const double CometocentricDistanceThrehold=5.0E3/tan(25.0/180.0*Pi);

  //set the location of the Sun, and shadowing flags according top the time of the observation


  rLocation=Vector3D::Length(xLocation);
  cosAngleThrouhold=(rLocation>3.0E3) ? cos(atan(5.0E3/rLocation)) : 0.0;

  //perform the integration loop
  for (iTest=0;iTest<nTotalTests;iTest++) {
    //generate the line direction
    if (rLocation<CometocentricDistanceThrehold) {
      //close to the nucleus limit
      Vector3D::Distribution::Uniform(l);
    }
    else {
      //far from the nucleus limit
      do {
        Vector3D::Distribution::Uniform(l);
        c=-Vector3D::DotProduct(xLocation,l)/rLocation;
      }
      while ((c<cosAngleThrouhold)||(c<0.0));
    }

    //find the first face of intersection with the nucleus
    iSurfaceElement=PIC::RayTracing::FindFistIntersectedFace(xLocation,l,xIntersection);

    //parameters of the primary species source at that surface element

    //pointing direcion of the surface element normal vector






  }


}



