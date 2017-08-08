//$Id$
//function related to the Saturnian system of the outer plamet magnetosphere model

/*
 * SaturnianSystem.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#include "MOP.h"

//rotation vector of Saturn
double MOP::SaturninanSystem::Saturn::RotationAxis[3]={0.0,0.0,0.0};

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

//get Saturn's magnetic dipole field
void MOP::SaturninanSystem::Saturn::GetMagneticFieldDipole(double *B,double *x) {

}

//calculate electric and magnetic fields in Saturn's magnetosphere
void MOP::SaturninanSystem::Magnetosphere::GetMagneticField(double *B,double *x) {
  KMAG::GetMagneticField(B,x);
}

double MOP::SaturninanSystem::Magnetosphere::GetCorotationSpeed(double *x) {
  static const double rmin=0.0;
  static const double rmax=20.0;
  static const int nPoints=51;
  static const double dr=(rmax-rmin)/(nPoints-1);

  //the corrotatoin velocity profile is digitized from
  static const double CorotationVelocity[nPoints]={108.571,4110.533,9059.282,12514.883,16904.430,19956.104,23724.578,27711.359,31308.219,35543.676,39451.062,42035.762,45033.734,
                                            48125.098,50812.531,52759.797,53344.691,56891.340,60098.277,61305.414,64553.215,73669.703,78616.109,79434.461,78681.492,78319.594,
                                            77280.594,76308.094,75764.078,76295.250,81022.195,87594.836,94218.867,99334.547,1.027e5,1.056e5,1.080e5,1.096e5,1.115e5,1.129e5,
                                            1.140e5,1.152e5,1.157e5,1.164e5,1.168e5,1.172e5,1.176e5,1.179e5,1.179e5,1.180e5,1.180e5};


  //calculate the corrotation speed
  double r,U;

  r=fabs(Vector3D::DotProduct(x,MOP::SaturninanSystem::Saturn::RotationAxis))/_SATURN__RADIUS_;

  if (r<rmin) U=CorotationVelocity[0];
  else if (r>rmax) U=CorotationVelocity[nPoints-1];
  else {
    int i;

    i=(int)((r-rmin)/dr);

    if (i<nPoints) {
      double t;

      t=(r-rmin)/dr-i;
      U=(1.0-t)*CorotationVelocity[i]+t*CorotationVelocity[i+1];
    }
    else U=CorotationVelocity[nPoints-1];
  }

  return U;
}

void MOP::SaturninanSystem::Magnetosphere::GetElectricField(double *E,double *x) {
  double U,LocalCorotationVelocity[3],B[3],c;
  int idim;

  U=GetCorotationSpeed(x);
  Vector3D::CrossProduct(LocalCorotationVelocity,x,MOP::SaturninanSystem::Saturn::RotationAxis);

  for (c=0.0,idim=0;idim<3;idim++) c+=pow(LocalCorotationVelocity[idim],2);
  if (c>0.0) for (c=U/sqrt(c),idim=0;idim<3;idim++) LocalCorotationVelocity[idim]*=c;

  //E=-v\timesB
  GetMagneticField(B,x);
  Vector3D::CrossProduct(E,B,LocalCorotationVelocity);
}








