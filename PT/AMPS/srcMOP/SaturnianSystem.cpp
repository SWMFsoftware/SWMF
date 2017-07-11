//$Id$
//function related to the Saturnian system of the outer plamet magnetosphere model

/*
 * SaturnianSystem.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#include "MOP.h"

double MOP::SaturninanSystem::GetLocalMeshResolution(double *x) {
  double Res=1*_SATURN__RADIUS_;

  const double rOrbitEnceladus=238000.0E3;

  //find the location of the point of the center of the Enceladus' orbit closest to the 'x'
  double r,r2=0.0;
  double l_xy[3]={x[0],x[1],x[2]};


  //detemine the plane defined by the or Enceladus
  static bool initflag=false;
  static double e0[3],e1[3],e2[3];

  if (initflag==false) {
    initflag=true;

    SpiceDouble lt,et,StateEnceladus[6];
    int idim;

    utc2et_c(Exosphere::SimulationStartTimeString,&et);
    spkezr_c("Enceladus",et,"SSO","none","Saturn",StateEnceladus,&lt);
    for (idim=0;idim<3;idim++) e0[idim]=StateEnceladus[idim];

    et+=5.0*3600.0;

    spkezr_c("Enceladus",et,"SSO","none","Saturn",StateEnceladus,&lt);
    for (idim=0;idim<3;idim++) e1[idim]=StateEnceladus[idim];

    //create the frame of reference related to the orbital plane of Enceladus
    Vector3D::Normalize(e0);
    Vector3D::Orthogonalize(e0,e1);
    Vector3D::Normalize(e1);
    Vector3D::CrossProduct(e2,e0,e1);
  }

  Vector3D::Orthogonalize(e2,l_xy);
  if (l_xy[0]*l_xy[0]+l_xy[1]*l_xy[1]+l_xy[2]*l_xy[2]<1.0E-10) return Res;

  Vector3D::Normalize(l_xy);

  for (int idim=0;idim<3;idim++) {
    l_xy[idim]*=rOrbitEnceladus;
    r2+=pow(x[idim]-l_xy[idim],2);
  }

  //parameters of the mesh resolution enhancement when a point close to the orbit of Enceladus
  const double R0=1.0*_ENCELADUS__RADIUS_,Res0=0.3*_ENCELADUS__RADIUS_;
  const double R1=15.0*_ENCELADUS__RADIUS_,Res1=1.0*_ENCELADUS__RADIUS_;

  r=sqrt(r2);

  if (r<R0) {
    Res=10*Res0;
  }
  else if (r<R1) {
    Res=10*(Res0+(r-R0)/(R1-R0)*(Res1-Res0));
  }

  return Res;
}



