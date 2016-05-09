//$Id$

/*
 * EnceladusMultiPlume_SourceLocation.cpp
 *
 *  Created on: May 3, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "Exosphere.h"
#include "EnceladusMultiPlume.h"

//the total number of the Tiger Stripes and Individual Plumes
//int EnceladusMultiPlume::nTotalTigerStripes=1;
//int EnceladusMultiPlume::nTotalIndividualPlumes=1;

//TigerStripes
EnceladusMultiPlume::cTigerStripe EnceladusMultiPlume::TigerStripeTable[EnceladusMultiPlume::nTotalTigerStripes];

//IndividualPlumes
EnceladusMultiPlume::cIndividualPlume EnceladusMultiPlume::IndividualPlumeTable[EnceladusMultiPlume::nTotalIndividualPlumes];

//convert the (Lat, wLon) location of the plumes into Cartesian coordinates
void EnceladusMultiPlume::RecalculateSourceLocations() {
  int idim;
  SpiceDouble x[3];
  SpiceInt      id;
  SpiceBoolean  found;

  //recalculate locations of the individual plumes
  bodn2c_c ( "ENCELADUS", &id, &found );

  for (int np=0;np<nTotalIndividualPlumes;np++) {
    //srfrec_c  converts planetocentric latitude and longitude of a surface
    //point on a specified body to rectangular coordinates.
    //planetocentic Longitude == east Longitude -> convert WLongitude to east Longitude before calling srfrec_c
    srfrec_c (id, (360.0-IndividualPlumeTable[np].wLon)*rpd_c(), IndividualPlumeTable[np].Lat*rpd_c(),x);
//    std::cout << "plume=" << np << ", x=" << x[0]*1.0E3 << "  " << x[1]*1.0E3 << "  " << x[2]*1.0E3 << std::endl;

    for (idim=0;idim<3;idim++) IndividualPlumeTable[np].xLoacation[idim]=x[idim]*1.0E3;
  }

  //determine the plume direction
  for (int np=0;np<nTotalIndividualPlumes;np++) {
    //get the plume direction vector
    //construc the coordinate system for the plume

    //radial direction -> z-axys
    //direction to the local north -> x-axys
    //y-axys -> the "right" coordinate system

    double l,l1,coord[3][3];

    //get the direction to the local north
    const double north[]={0.0,0.0,1.0};

    for (idim=0,l=0.0;idim<3;idim++) l+=pow(IndividualPlumeTable[np].xLoacation[idim],2);
    for (idim=0,l=sqrt(l);idim<3;idim++) coord[2][idim]=IndividualPlumeTable[np].xLoacation[idim]/l;

    for (idim=0,l=0.0;idim<3;idim++) l+=north[idim]*coord[2][idim];

    for (idim=0,l1=0.0;idim<3;idim++) {
      coord[0][idim]=north[idim]-l*coord[2][idim];
      l1+=pow(coord[0][idim],2);
    }

    for (idim=0,l1=sqrt(l1);idim<3;idim++) coord[0][idim]/=l1;

    coord[1][0]=-(coord[2][1]*coord[0][2]-coord[0][1]*coord[2][2]);
    coord[1][1]=+(coord[2][0]*coord[0][2]-coord[0][0]*coord[2][2]);
    coord[1][2]=-(coord[2][0]*coord[0][1]-coord[0][0]*coord[2][1]);

    for (idim=0;idim<3;idim++) IndividualPlumeTable[np].lPointing[idim]=cos(IndividualPlumeTable[np].TiltAngle*rpd_c())*coord[2][idim];   //vertical direction
    for (idim=0;idim<3;idim++) IndividualPlumeTable[np].lPointing[idim]+=sin(IndividualPlumeTable[np].TiltAngle*rpd_c())*sin(IndividualPlumeTable[np].AzimuthAngle*rpd_c())*coord[1][idim]; //direction to the east
    for (idim=0;idim<3;idim++) IndividualPlumeTable[np].lPointing[idim]+=sin(IndividualPlumeTable[np].TiltAngle*rpd_c())*cos(IndividualPlumeTable[np].AzimuthAngle*rpd_c())*coord[0][idim]; //direction to the north

    Vector3D::Normalize(IndividualPlumeTable[np].lPointing);
  }
}



