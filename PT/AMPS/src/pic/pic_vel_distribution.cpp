//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the function for generation of particles velocty distribution

#include "pic.h"

//get particle velocity distribution with Maxwellian
void PIC::Distribution::MaxwellianVelocityDistribution(double *v,double *BulkFlowVelocity,double Temp,int spec) {
  double beta;
  int idim;

  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*Temp));

  for (idim=0;idim<DIM;idim++) v[idim]=sqrt(-log(rnd()))/beta*sin(2.0*Pi*rnd())+BulkFlowVelocity[idim];
}


//====================================================
//Get the particle velocity that is injected with Maxwellian distribution
double PIC::Distribution::InjectMaxwellianDistribution(double *v,double *BulkFlowVelocity,double Temp,double *ExternalNormal,int spec,int WeightCorrectionMode) {
  int idim;
  double sc,beta,u,c,a,c1;
  register double dotProduct=0.0;

  double ParticleWeightCorrection=1.0;

#if DIM != 1
  double v1[3],v2[3],v3[3]={0.0,0.0,0.0};
#endif


  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*Temp));

  do {
    dotProduct=0.0;
    for (idim=0;idim<DIM;idim++) dotProduct+=BulkFlowVelocity[idim]*ExternalNormal[idim];
    sc=beta*fabs(dotProduct);

    do {
      if (dotProduct<=0.0) {
        do u=-10.0+20.0*rnd(); while(u+sc<=0.0);
        c=sqrt(2.0+sc*sc);
        a=2.0*(u+sc)/(sc+c)*exp(0.5+0.5*sc*(sc-c)-u*u);
      }
      else {
        do u=-10.0+20.0*rnd(); while(u+sc>=0.0);
        c=0.5*(-sqrt(sc*sc+2.0)-sc);
        a=(u+sc)*exp(-u*u)/(c+sc)/exp(-c*c);
      }

      if (WeightCorrectionMode==_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_) {
        ParticleWeightCorrection=a;
        break;
      }
    }
    while (a<rnd());

#if DIM == 1
    v[0]=-(fabs(u+sc)/beta)*ExternalNormal[0];

    c=sqrt(-log(rnd()))/beta;
    c1=rnd();

    v[1]=c*sin(2.0*Pi*c1);
    v[2]=c*cos(2.0*Pi*c1);
#elif DIM == 2
    exit(__LINE__,__FILE__,"not tested yet");

    register double t;
    t=-(fabs(u+sc)/beta);
    v3[0]=t*ExternalNormal[0],v3[1]=t*ExternalNormal[1],v3[2]=0.0;

    v2[0]=ExternalNormal[1],v2[1]=-ExternalNormal[0],v2[2]=0.0;
    v1[0]=0.0,v1[1]=0.0,v1[2]=1.0;

    c=sqrt(-log(rnd()))/beta;
    c1=rnd();

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v1[idim];
    for (idim=0;idim<3;idim++) v1[idim]=(c*sin(2.0*Pi*c1)+dotProduct)*v1[idim];

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v2[idim];
    for (idim=0;idim<3;idim++) v2[idim]=(c*cos(2.0*Pi*c1)+dotProduct)*v2[idim];

    for (idim=0;idim<3;idim++) v[idim]=v1[idim]+v2[idim]+v3[idim];
#elif DIM == 3
    register double t;

    for (t=-(fabs(u+sc)/beta),idim=0;idim<DIM;idim++) v3[idim]=t*ExternalNormal[idim];

    if (fabs(ExternalNormal[0])>1.0E-3) {
      v2[0]=ExternalNormal[1],v2[1]=-ExternalNormal[0],v2[2]=0.0;
    }
    else {
      v2[0]=0.0,v2[1]=ExternalNormal[2],v2[2]=-ExternalNormal[1];
    }

    for (dotProduct=0.0,idim=0;idim<DIM;idim++) dotProduct+=v2[idim]*v2[idim];
    for (dotProduct=sqrt(dotProduct),idim=0;idim<DIM;idim++) v2[idim]/=dotProduct;

    v1[0]=v2[1]*ExternalNormal[2]-v2[2]*ExternalNormal[1];
    v1[1]=-(v2[0]*ExternalNormal[2]-v2[2]*ExternalNormal[0]);
    v1[2]=v2[0]*ExternalNormal[1]-v2[1]*ExternalNormal[0];

    for (dotProduct=0.0,idim=0;idim<DIM;idim++) dotProduct+=v1[idim]*v1[idim];
    for (dotProduct=sqrt(dotProduct),idim=0;idim<DIM;idim++) v1[idim]/=dotProduct;

    c=sqrt(-log(rnd()))/beta;
    c1=rnd();

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v1[idim];
    for (idim=0;idim<3;idim++) v1[idim]=(c*sin(2.0*Pi*c1)+dotProduct)*v1[idim];

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v2[idim];
    for (idim=0;idim<3;idim++) v2[idim]=(c*cos(2.0*Pi*c1)+dotProduct)*v2[idim];

    for (idim=0;idim<3;idim++) v[idim]=v1[idim]+v2[idim]+v3[idim];
#else
    exit(__LINE__,__FILE__,"Error: unknown option");
#endif


    for (dotProduct=0.0,idim=0;idim<DIM;idim++) dotProduct+=v[idim]*ExternalNormal[idim];
  }
  while (dotProduct>=0.0);

  return ParticleWeightCorrection;
}

