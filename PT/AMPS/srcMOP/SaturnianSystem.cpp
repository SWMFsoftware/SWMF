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

  //find the location of the point of the center of the Enceladus' obbit closest to the 'x'
  int idim;
  double r,r2=0.0,lz[3]={0.0,0.0,1.0};
  double l_xy[3]={x[0],x[1],x[2]};

  Vector3D::Orthogonalize(lz,l_xy);
  if (l_xy[0]*l_xy[0]+l_xy[1]*l_xy[1]+l_xy[2]*l_xy[2]<1.0E-10) return Res;

  Vector3D::Normalize(l_xy);

  for (idim=0;idim<3;idim++) {
    l_xy[idim]*=((idim!=2) ? rOrbitEnceladus : 0.0);
    r2+=pow(x[idim]-l_xy[idim],2);
  }

  //parameters of the mesh resolution enhancement when a point close to the orbit of Enceladus
  const double R0=1.0*_ENCELADUS__RADIUS_,Res0=0.3*_ENCELADUS__RADIUS_;
  const double R1=5.0*_ENCELADUS__RADIUS_,Res1=1.0*_ENCELADUS__RADIUS_;

  r=sqrt(r2);

  if (r<R0) Res=10*Res0;
  else if (r<R1) Res=10*(Res0+(r-R0)/(R1-R0)*(Res1-Res0));

  return Res;
}



