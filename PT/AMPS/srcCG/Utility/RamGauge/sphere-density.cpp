

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <vector>

#define _DO_NOT_LOAD_GLOBAL_H_
#define _COMPILATION_MODE_ _COMPILATION_MODE__MPI_

#include "../../../src/general/rnd.h"
#include "../../../src/general/specfunc.h"

namespace Gauge {
  double SphereRadius=1.0;
  double OpeningRadius=0.1;
  double Temp=293.0;
  double MolMass=18*_AMU_;
}

int nTotalTests=100000;
double AveragingPhysicalTime=100.0;

double BackgroundDensity=100.0;
double BackgroundVelocity[3]={0.0,0.0,-100.0};

void EvaluateCase(double& Density) {
  int idim,iTest,nface,nd0,nd1,nd2;
  double ParticleWeight,x[3],v[3],norm[3],x0[3],e0[3],e1[3],ParticleTimeCounter=0.0;


  //get the numerical flux through the face 
  double NumericalFlux=Pi*pow(Gauge::OpeningRadius,2)*Vector3D::Length(BackgroundVelocity);

  //if the numerical flux is positive -> determine the contribution of the flux to the density inside the gauge 
  double t1,t2,tmin,tmax,tStart;

  ParticleWeight=NumericalFlux*AveragingPhysicalTime/nTotalTests;
    
  for (iTest=0;iTest<nTotalTests;iTest++) {
    //generate a new particle 
    double a=0.0,b=0.0,c=0.0,d,t1,t2;

    Vector3D::Distribution::Circle::Unifirm(x,Gauge::OpeningRadius);

    d=pow(Gauge::OpeningRadius,2)-x[0]*x[0]-x[1]*x[1];
    x[2]=(d>0.0) ? sqrt(d) : Gauge::OpeningRadius;   

    for (idim=0;idim<3;idim++) {
      v[idim]=BackgroundVelocity[idim];
    }

    do {
      //calculate the flight time to the next interaction with the walls of the sphere
      a=0.0,b=0.0,c=0.0,d,tmin; 
    
      for (idim=0;idim<3;idim++) {
        a+=pow(v[idim],2);
        b+=2.0*v[idim]*x[idim];
        c+=pow(x[idim],2);
      }

      c-=pow(Gauge::SphereRadius,2);
      d=pow(b,2)-4.0*a*c;

      d=sqrt(d);
      t1=(-b+d)/(2.0*a);
      t2=(-b-d)/(2.0*a);

      tmax=std::max(t1,t2);

      //increment the time counter
      ParticleTimeCounter+=tmax*ParticleWeight;
      
     //determine the new direction of the particle velocity
     for (idim=0;idim<3;idim++) x[idim]+=v[idim]*tmax;
 
      do {
        Vector3D::Distribution::Uniform(v);

        double beta=sqrt(Gauge::MolMass/(2.0*Kbol*Gauge::Temp));

        for (idim=0;idim<3;idim++) v[idim]=sqrt(-log(rnd()))/beta*sin(2.0*Pi*rnd());
      }
      while (Vector3D::DotProduct(x,v)>0.0);
    }
    while ((x[0]*x[0]+x[1]*x[1]>pow(Gauge::OpeningRadius,2))||(x[2]<0.0));
 }


//  std::cout << "Density=" << ParticleTimeCounter/AveragingPhysicalTime/(Gauge::Height*Pi*pow(Gauge::Radius,2)) << ", Volume=" << Gauge::Height*Pi*pow(Gauge::Radius,2) << std::endl;

  Density=ParticleTimeCounter/AveragingPhysicalTime/(4.0/3.0*Pi*pow(Gauge::SphereRadius,3));
}


int main() {
  int i;
  double phi,Density;

  const double Speed=100.0;

  std::cout << "Balance of fluxes: " <<  4.0*Speed/sqrt(8.0*Kbol*Gauge::Temp/(Pi*Gauge::MolMass)) << "\n";

  for (i=1;i<10;i++) {
    Gauge::OpeningRadius=0.5*Gauge::SphereRadius/double(i);

    BackgroundVelocity[0]=-Speed*cos(phi);
    BackgroundVelocity[1]=0.0;
    BackgroundVelocity[2]=-Speed*sin(phi);

    EvaluateCase(Density);
    std::cout << Gauge::OpeningRadius << "\t  " <<  Density << "\t  " << 1.0+sin(phi) << std::endl;
  }
  
  return 1;
}





 

