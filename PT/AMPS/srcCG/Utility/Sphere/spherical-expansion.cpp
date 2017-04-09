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
const int nTotalTests=100000;
const int nRadialIntervals=200;

const int VelocityDistributionMode__Maxwellian=0;
const int VelocityDistributionMode__Random=1;

const int VelocityDistributionMode=VelocityDistributionMode__Maxwellian;

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
  double speed,f,fmax,beta,ParticleWeight,TimeStep,x[3],v[3],ExternalNormal[3],r;

  //init the sampling mesh
  double *SamplingBufferDensity=new double[nRadialIntervals];
  double *SamplingBufferFlux=new double[nRadialIntervals];
  double dr=MaxAltitude/nRadialIntervals;
  double MeanGeneratedSpeed=0.0;

  for (i=0;i<nRadialIntervals;i++) SamplingBufferDensity[i]=0.0,SamplingBufferFlux[i]=0.0;
  
  //get the particle weight, and the time step 
  ParticleWeight=Flux*4.0*Pi*pow(R0,2) * SamplingTime/nTotalTests;
  TimeStep=dr/1.0E3/8.0;
  
  beta=sqrt(MolMass/(2.0*Kbol*Temp)); 
  fmax=2.0*pow(1.0/beta,2)*exp(-1.0);
  

  //perform testing
  for (iTest=0;iTest<nTotalTests;iTest++) {
    //generate a new particle 
    Vector3D::Distribution::Uniform(x);

    for (idim=0;idim<3;idim++) {
      ExternalNormal[idim]=-x[idim];
      x[idim]*=R0;
    }

    switch (VelocityDistributionMode) {
    case VelocityDistributionMode__Maxwellian: 
      InjectMaxwellianDistribution(v,ExternalNormal); 
      break;
    case VelocityDistributionMode__Random:
      //distribute the direction of the vector
      do {
        Vector3D::Distribution::Uniform(v);
      }
      while(Vector3D::DotProduct(v,x)<=0.0); 

      //distribute the speed 
      do {
        speed=5.0/beta*rnd();
      }
      while (pow(speed,2)*exp(-pow(speed*beta,2))/fmax<rnd()); 
  
//speed=4.848213e+02;

      for (idim=0;idim<3;idim++) v[idim]*=speed;
      MeanGeneratedSpeed+=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    }

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
  printf("Mean Generated Speed=%e\n",MeanGeneratedSpeed/nTotalTests); 
  printf("VAIABLES=\"r\", \"Altitude\", \"Density\", \"Total Flux\", \"Mean Speed\"\n");

  for (i=0;i<nRadialIntervals;i++) {
    double Volume=4.0/3.0*Pi*(pow(R0+(i+1)*dr,3)-pow(R0+i*dr,3));

    printf("%e %e %e %e %e\n",R0+(i+0.5)*dr,(i+0.5)*dr,SamplingBufferDensity[i]/SamplingTime/Volume,SamplingBufferFlux[i]/SamplingTime/Volume,SamplingBufferFlux[i]/SamplingBufferDensity[i]);
  }

  //evaluate the density, velocity, and flux with a Monte Carlo technique 
  double A,r2,cosTheta,SpeedSample=0.0,DensitySample,FluxSample,ll[3];  
  double xSample[3]={0.0,0.0,18.6E3}; //IMPORTANT: the pont of sampling is on the POSITIVE Z-AXIS 

  DensitySample=0.0,FluxSample=0.0;
  A=Flux*2.0*pow(beta,4)/Pi;

  for (iTest=0;iTest<nTotalTests;iTest++) { 
    Vector3D::Distribution::Uniform(x);
    r2=0.0,cosTheta=0.0; 

    if (x[2]<0.0) continue; //the sampling point cannot be seen from the point 'x'

    for (idim=0;idim<3;idim++) {
      x[idim]*=R0;

      ll[idim]=xSample[idim]-x[idim];
      r2+=pow(ll[idim],2);
      cosTheta+=x[idim]*ll[idim];
    }

    cosTheta/=R0*sqrt(r2);

    //increment the sampled parameters 
    DensitySample+=cosTheta*sqrt(Pi)/(4.0*r2*pow(beta,3));
    FluxSample+=cosTheta/(2.0*r2*pow(beta,4));
  }

  DensitySample*=A*(4.0*Pi*pow(R0,2))/nTotalTests;
  FluxSample*=A*(4.0*Pi*pow(R0,2))/nTotalTests;
 
  printf("Monte Carlo Test at xSample=%e, %e,%e\nDensity=%e\nFlux=%e\nSpeed=%e\n",xSample[0],xSample[1],xSample[2],DensitySample,FluxSample,FluxSample/DensitySample);

  std::cout << "2.0/(beta*sqrt(Pi))=" << 2.0/(beta*sqrt(Pi)) << std::endl;

  //macrospcopic paramters obtaind with the parmicle model at the same altitude
   i=(sqrt(xSample[0]*xSample[0]+-xSample[1]*xSample[1]+xSample[2]*xSample[2])-R0)/dr;

  printf("VAIABLES=\"r\", \"Altitude\", \"Density\", \"Total Flux\", \"Mean Speed\"\n");
  double Volume=4.0/3.0*Pi*(pow(R0+(i+1)*dr,3)-pow(R0+i*dr,3));
  printf("%e %e %e %e %e\n",R0+(i+0.5)*dr,(i+0.5)*dr,SamplingBufferDensity[i]/SamplingTime/Volume,SamplingBufferFlux[i]/SamplingTime/Volume,SamplingBufferFlux[i]/SamplingBufferDensity[i]);
          
  return 1;
}
