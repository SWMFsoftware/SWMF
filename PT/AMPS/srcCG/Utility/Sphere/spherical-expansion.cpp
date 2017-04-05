//$Id$


#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <vector>

#define _DO_NOT_LOAD_GLOBAL_H_
#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_

#include "../../../src/general/rnd.h"
#include "../../../src/general/specfunc.h"

const double Flux=1.0;
const double MaxAltitude=20.0E3;
const double R0=1.73E3;
const double Temp=200.0;
const double MolMass=_H2O__MASS_;

const double SamplingTime=3600.0;
const int nTotalTests=10000;
const int nRadialIntervals=200;

//====================================================
//Get the particle velocity that is injected with Maxwellian distribution
double InjectMaxwellianDistribution(double *v,double *ExternalNormal) {
  int idim;
  double sc,beta,u,c,a,c1;
  register double dotProduct=0.0;

  double ParticleWeightCorrection=1.0;

  double v1[3],v2[3],v3[3]={0.0,0.0,0.0};


  beta=sqrt(MolMass/(2.0*Kbol*Temp));

  do {
    dotProduct=0.0;
    sc=0.0; //beta*fabs(dotProduct);

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

    }
    while (a<rnd());

    register double t;

    for (t=-(fabs(u+sc)/beta),idim=0;idim<3;idim++) v3[idim]=t*ExternalNormal[idim];

    if (fabs(ExternalNormal[0])>1.0E-3) {
      v2[0]=ExternalNormal[1],v2[1]=-ExternalNormal[0],v2[2]=0.0;
    }
    else {
      v2[0]=0.0,v2[1]=ExternalNormal[2],v2[2]=-ExternalNormal[1];
    }

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=v2[idim]*v2[idim];
    for (dotProduct=sqrt(dotProduct),idim=0;idim<3;idim++) v2[idim]/=dotProduct;

    v1[0]=v2[1]*ExternalNormal[2]-v2[2]*ExternalNormal[1];
    v1[1]=-(v2[0]*ExternalNormal[2]-v2[2]*ExternalNormal[0]);
    v1[2]=v2[0]*ExternalNormal[1]-v2[1]*ExternalNormal[0];

    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=v1[idim]*v1[idim];
    for (dotProduct=sqrt(dotProduct),idim=0;idim<3;idim++) v1[idim]/=dotProduct;

    c=sqrt(-log(rnd()))/beta;
    c1=rnd();

//    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v1[idim];
    for (idim=0;idim<3;idim++) v1[idim]=(c*sin(2.0*Pi*c1)+dotProduct)*v1[idim];

//    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=BulkFlowVelocity[idim]*v2[idim];
    for (idim=0;idim<3;idim++) v2[idim]=(c*cos(2.0*Pi*c1)+dotProduct)*v2[idim];

    for (idim=0;idim<3;idim++) v[idim]=v1[idim]+v2[idim]+v3[idim];


    for (dotProduct=0.0,idim=0;idim<3;idim++) dotProduct+=v[idim]*ExternalNormal[idim];
  }
  while (dotProduct>=0.0);

  return 0.0;
}


int main() {
  int i,iTest,idim;
  double ParticleWeight,TimeStep,x[3],v[3],ExternalNormal[3],r;

  //init the sampling mesh
  double *SamplingBufferDensity=new double[nRadialIntervals];
  double *SamplingBufferFlux=new double[nRadialIntervals];
  double dr=MaxAltitude/nRadialIntervals;

  for (i=0;i<nRadialIntervals;i++) SamplingBufferDensity[i]=0.0,SamplingBufferFlux[i]=0.0;
  
  //get the particle weight, and the time step 
  ParticleWeight=Flux*4.0*Pi*pow(R0,2) * SamplingTime/nTotalTests;
  TimeStep=dr/1.0E3/8.0;

  //perform testing
  for (iTest=0;iTest<nTotalTests;iTest++) {
    //generate a new particle 
    Vector3D::Distribution::Uniform(x);

    for (idim=0;idim<3;idim++) {
      ExternalNormal[idim]=-x[idim];
      x[idim]*=R0;
    }

    InjectMaxwellianDistribution(v,ExternalNormal); 

    //trace the particle until it left the domain 
    do {
      for (idim=0;idim<3;idim++) x[idim]+=v[idim]*TimeStep;

      r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      i=(r-R0)/dr;

      if (i<nRadialIntervals) {
        SamplingBufferDensity[i]+=ParticleWeight*TimeStep;
        SamplingBufferFlux[i]+=ParticleWeight*sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*TimeStep;
      }
    }
    while (i<nRadialIntervals); 
  }

  //output sampled results 
  printf("VAIABLES=\"r\", \"Altitude\", \"Density\", \"Total Flux\"\n");

  for (i=0;i<nRadialIntervals;i++) {
    double Volume=4.0/3.0*Pi*(pow(R0+(i+1)*dr,3)-pow(R0+i*dr,3));

    printf("%e %e %e %e\n",R0+(i+0.5)*dr,(i+0.5)*dr,SamplingBufferDensity[i]/SamplingTime/Volume,SamplingBufferFlux[i]/SamplingTime/Volume);
  }
 
  return 1;
}
