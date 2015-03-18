//$Id$
//the simeple test particle model

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

const double Pi=2.0*acos(0.0);
const double PiTimes2=2.0*Pi;
const double Kbol=1.3806503E-23;
const double GravityConstant=6.67300E-11;
const double eV2J=1.602176565E-19;
const double _AMU_=1.66053886E-27;

#define _GEOMETRY_MODE__LINEAR_   0
#define _GEOMETRY_MODE__SPHERICAL_ 1

#define _BOUNDARY_MODE__LEFT_DOMAIN_     0
#define _BOUNDARY_MODE__SCATTERED_BACK_  1

#define _DENSITY_SAMPLING__TOTAL_                      0
#define _DENSITY_SAMPLING__INITAL_INJECTION_           1
#define _DENSITY_SAMPLING__LEFT_BOUNDARY_REFLECTION_   2
#define _DENSITY_SAMPLING__RIGHT_BOUNDARY_REFLECTION_  3

#define _VELOCITY_DISTRIBUTION_SAMPLING__INITIAL_INJECTION_          0
#define _VELOCITY_DISTRIBUTION_SAMPLING__LEFT_BOUNDARY_REFLECTION_   1


unsigned long int rndLastSeed;

double rnd() {
  rndLastSeed*=48828125;
  rndLastSeed&=2147483647; // pow(2,31) - 1
  if (rndLastSeed==0) rndLastSeed=1;

  return double(rndLastSeed/2147483648.0); //(pow(2,31) - 1) + 1
}


//---------------------------  The model specific data ---------------------------------
#define _GEOMETRY_MODE_ _GEOMETRY_MODE__SPHERICAL_

const int nDensitySampleIntervals=10000;
const int nVelocityDistributionSampledIntervals=500;
const double MaxVelocity=10.0E3;
const double dVelocityDistributionInterval=MaxVelocity/nVelocityDistributionSampledIntervals;

const int nTotalTestParticles=2000;

const double xMin=1.569E6;
const double xMax=3.0*xMin;

const double ParticleMass=32*_AMU_;
const double SourceRate=1.0E27;

double VelocityDistributionSampling[2][nVelocityDistributionSampledIntervals];
double DensitySampling[4][nDensitySampleIntervals];
int ModelParticleState;

double SetSamplingTime() {  //set up the time interval for sampling. must be more that the time needed for a particle to pass through the domain
  return (xMax-xMin)*10.0/500.0;
}

double SetTimeStep() { //evaluate the time set used for integration of the particle trajectory
  return (xMax-xMin)/nDensitySampleIntervals/MaxVelocity/2.0;
}

void GetParticleAcceleration(double *a,double *x) { //calcualte the acceleration of a particle
  double r2=0.0,c;
  int idim;

  static const double MassEuropa=4.8E22;

  for (idim=0;idim<3;idim++) r2+=pow(x[idim],2);
  c=GravityConstant*MassEuropa/pow(r2,1.5);

  for (idim=0;idim<3;idim++) a[idim]=-c*x[idim];


  a[0]=0.0,a[1]=0.0,a[2]=0.0;
}

void InitParticleVelocity(double *v) { //init the particle velocity
  double r,E,Speed,lVel[3];

  static const double _O2__MASS_=32*_AMU_;

  static const double U_o2=0.015*eV2J;  //the community accepted parameter of the energy distribution      //Burger 2010-SSR
  static const double Emax_o2=_O2__MASS_*MaxVelocity*MaxVelocity/2.0; //the maximum energy of the ejected particle

  do {
    r=rnd();
    E=r*U_o2/(1.0-r);    //Burger 2010-SSR
  }
  while (E>Emax_o2);

  Speed=sqrt(E*2.0/_O2__MASS_);

  double Phi      = PiTimes2*rnd();
  double SinTheta = pow(rnd(),0.5);
  double ExternalNormal[3]={1.0,0.0,0.0};


  lVel[0] = cos(Phi) * SinTheta;
  lVel[1] = sin(Phi) * SinTheta;
  lVel[2] = pow(1.0 - SinTheta*SinTheta, 0.5);
  // actual normal is obtained from {0,0,1} by 2 rotations:
  //   1: around initial z axis by angle A
  //   2: around new     y axis by angle B
  //
  // rotation matrix:  / CosB*CosA -SinA SinB*CosA \
  //                  |  CosB*SinA  CosA SinB*SinA  |
  //                   \-SinB         0  CosB      /
  //
  //         => normal = {SinB*CosA, SinB*SinA, CosB}
  double CosB = ExternalNormal[2];

  if(CosB < 1.0 && CosB > -1.0) {
    double SinB = pow(1.0 - CosB*CosB, 0.5);
    double CosA, SinA;
    CosA = ExternalNormal[0] / SinB;
    SinA = ExternalNormal[1] / SinB;
    v[0]=Speed*( lVel[0]*CosB*CosA - lVel[1]*SinA + lVel[2]*SinB*CosA);
    v[1]=Speed*( lVel[0]*CosB*SinA + lVel[1]*CosA + lVel[2]*SinB*SinA);
    v[2]=Speed*(-lVel[0]*SinB                     + lVel[2]*CosB     );
  }
  else    // if abs(CosB)==1 no rotation is needed
    for (int idim=0;idim<3;idim++) v[idim]=lVel[idim]*Speed;



  v[0]=2000.0,v[1]=0.0,v[2]=0.0;
}


int ProcessLeftBoundaryIntersection(double *x,double *v) {  //process boundary intersection at x==xmin
  double r2,c;
  int idim;

  static const double SurfaceTemeprature=90.0;

  //move the particle to the surface of the sphere
  for (r2=0.0,idim=0;idim<3;idim++) r2+=pow(x[idim],2);
  c=xMin/sqrt(r2);
  for (idim=0;idim<3;idim++) x[idim]*=c;


  //init velocity of the particle
  double c1,beta=sqrt(ParticleMass/(2.0*Kbol*SurfaceTemeprature));

  v[0]=sqrt(-log(rnd()))/beta;  //radial velocity

  c=sqrt(-log(rnd()))/beta;
  c1=rnd();

  v[1]=c*sin(2.0*Pi*c1);  //tangential velocity
  v[2]=c*cos(2.0*Pi*c1);  //tangential velocity

  return _BOUNDARY_MODE__SCATTERED_BACK_;
}

int ProcessRightBoundaryIntersection(double *x,double *v) { //process boundary intersection at x==xmax
  return _BOUNDARY_MODE__LEFT_DOMAIN_;
}

//---------------------------- The service functions -----------------------------------
void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);

  if (seed==-1) seed=thread;

  rndLastSeed=seed;
}

int main(int argc, char **argv) {
  int thread,i;
  double SamplingTime,TimeStep;
  double const dx=(xMax-xMin)/nDensitySampleIntervals;
  double ParticleWeight,TimeCounter;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);

  //init the random number generator
  rnd_seed(-1);
  SamplingTime=SetSamplingTime();
  TimeStep=SetTimeStep();

  //init the sampling buffers
  for (int mode=0;mode<2;mode++) for (i=0;i<nVelocityDistributionSampledIntervals;i++) VelocityDistributionSampling[mode][i]=0.0;
  for (int mode=0;mode<4;mode++) for (i=0;i<nDensitySampleIntervals;i++) DensitySampling[mode][i]=0.0;


  //calcualte the particle weight
  ParticleWeight=SourceRate*SamplingTime/nTotalTestParticles;

  //loop through the particles
  double halfTimeStep=0.5*TimeStep;
  int nTestParticles,idim;
  double v[3],a[3],x[3],vMiddle[3],aMiddle[3],xMiddle[3];


  for (nTestParticles=0;nTestParticles<2*nTotalTestParticles;nTestParticles++) {
    TimeCounter=(2.0*rnd()-1.0)*SamplingTime;
    InitParticleVelocity(v);

    ModelParticleState=_DENSITY_SAMPLING__INITAL_INJECTION_;
    x[0]=xMin,x[1]=0.0,x[2]=0.0;

    //sample the initial velocity of the injected particles
    double Speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    int iSpeedInterval=Speed/dVelocityDistributionInterval;
    if (iSpeedInterval<nVelocityDistributionSampledIntervals) VelocityDistributionSampling[_VELOCITY_DISTRIBUTION_SAMPLING__INITIAL_INJECTION_][iSpeedInterval]++;


    while (TimeCounter<SamplingTime) {
      //predictor
      GetParticleAcceleration(a,x);

      for (idim=0;idim<3;idim++) {
        xMiddle[idim]=x[idim]+halfTimeStep*v[idim];
        vMiddle[idim]=v[idim]+halfTimeStep*a[idim];
      }

      //corrector
      GetParticleAcceleration(aMiddle,xMiddle);

      for (idim=0;idim<3;idim++) {
        x[idim]+=TimeStep*vMiddle[idim];
        v[idim]+=TimeStep*aMiddle[idim];
      }

      TimeCounter+=TimeStep;

      //in case of the spherical symmetry rotate the particle's location to the x-axis
      if (_GEOMETRY_MODE_==_GEOMETRY_MODE__SPHERICAL_) {
        double r=0.0,vRadial=0.0,vTangential=0.0,l[3];

        for (idim=0;idim<3;idim++) r+=pow(x[idim],2);
        r=sqrt(r);

        for (idim=0;idim<3;idim++) {
          l[idim]=x[idim]/r;
          vRadial+=v[idim]*l[idim];
        }

        for (idim=0;idim<3;idim++) vTangential+=pow(v[idim]-vRadial*l[idim],2);

        vTangential=sqrt(vTangential);

        //set up the new velocity and location vectors
        x[0]=r,x[1]=0.0,x[2]=0.0;
        v[0]=vRadial,v[1]=vTangential,v[2]=0.0;
      }

      //sample the particle data
      int ncell,code;

      if (x[0]<=xMin) {
        code=ProcessLeftBoundaryIntersection(x,v);
        ModelParticleState=_DENSITY_SAMPLING__LEFT_BOUNDARY_REFLECTION_;
        if (code==_BOUNDARY_MODE__LEFT_DOMAIN_) break; // the particle left domain and tracking of its trajectory is terminated

        //sample the initial velocity of the injected particles
        double Speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        int iSpeedInterval=Speed/dVelocityDistributionInterval;
        if (iSpeedInterval<nVelocityDistributionSampledIntervals) VelocityDistributionSampling[_VELOCITY_DISTRIBUTION_SAMPLING__LEFT_BOUNDARY_REFLECTION_][iSpeedInterval]++;

      }
      else if (x[0]>=xMax) {
        code=ProcessRightBoundaryIntersection(x,v);
        ModelParticleState=_DENSITY_SAMPLING__RIGHT_BOUNDARY_REFLECTION_;
        if (code==_BOUNDARY_MODE__LEFT_DOMAIN_) break; // the particle left domain and tracking of its trajectory is terminated
      }
      else if (TimeCounter>0.0) {
        ncell=(int)((x[0]-xMin)/dx);
        DensitySampling[0][ncell]+=TimeStep*ParticleWeight;
        DensitySampling[ModelParticleState][ncell]+=TimeStep*ParticleWeight;
      }

    }
  }




  //collect sampled data from all processors and output density into the data file
  double DensitySamplingBuffer[4][nDensitySampleIntervals];
  double VelocityDistributionBuffer[2][nVelocityDistributionSampledIntervals];

  MPI_Reduce(&DensitySampling[0][0],&DensitySamplingBuffer[0][0],4*nDensitySampleIntervals,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&VelocityDistributionSampling[0][0],&VelocityDistributionBuffer[0][0],2*nVelocityDistributionSampledIntervals,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (thread==0) {
    FILE *fout;

    //calculate density
    for (int ncell=0;ncell<nDensitySampleIntervals;ncell++) {
      double c=1.0/SamplingTime/(4.0/3.0*Pi*(pow(xMin+(1+ncell)*dx,3)-pow(xMin+ncell*dx,3)));

      for (int mode=0;mode<4;mode++) DensitySamplingBuffer[mode][ncell]*=c;
    }

    //save the density profile
    fout=fopen("data.dat","w");
    fprintf(fout,"VARIABLES=\"x\", \"Altitude\", \"Total Density\", \"Density_INITAL_INJECTION_\", \"Density_LEFT_BOUNDARY_REFLECTION_\", \"Density_RIGHT_BOUNDARY_REFLECTION_\"\n");


    for (int ncell=0;ncell<nDensitySampleIntervals;ncell++) {
      fprintf(fout,"%e  %e  %e %e  %e  %e \n",xMin+(0.5+ncell)*dx,(0.5+ncell)*dx,
          DensitySamplingBuffer[0][ncell],DensitySamplingBuffer[1][ncell],
          DensitySamplingBuffer[2][ncell],DensitySamplingBuffer[3][ncell]
      );

    }

    fclose(fout);

    //save the velocity distribution
    fout=fopen("VelocityDistributionFunction.dat","w");
    fprintf(fout,"VARIABLES=\"v\", \"fInit(v)\", \"fLeftBoundaryReflection(v)\"\n");

    double fmax[2]={-1.0,-1.0};

    for (int mode=0;mode<2;mode++) {
      for (i=0;i<nVelocityDistributionSampledIntervals;i++) if (fmax[mode]<VelocityDistributionBuffer[mode][i]) fmax[mode]=VelocityDistributionBuffer[mode][i];
      for (i=0;i<nVelocityDistributionSampledIntervals;i++) if (fmax[mode]>0.0) VelocityDistributionBuffer[mode][i]/=fmax[mode];
    }

    for (i=0;i<nVelocityDistributionSampledIntervals;i++) {
      fprintf(fout,"%e  %e  %e\n",(0.5+i)*dVelocityDistributionInterval,VelocityDistributionBuffer[0][i],VelocityDistributionBuffer[1][i]);
    }

    fclose(fout);

    //calculate the column integral
    fout=fopen("ColumnDentiyIntegral.dat","w");
    fprintf(fout,"VARIABLES=\"x\", \"Altitude\", \"Column Density\"\n");

    for (int ncell=0;ncell<nDensitySampleIntervals;ncell++) {
      double x[3],l[3],res=0.0,dens,r,IntegratedLength=0.0;
      int idim,i;

      //initial value
      x[0]=xMin+(0.5+ncell)*dx,x[1]=-xMax,x[2]=0.0;
      l[0]=0.0,l[1]=1.0,l[2]=0.0;

      while (IntegratedLength<2.0*xMax) {
        r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        i=(int)((r-xMin)/dx);

        if ((xMin<=r)&&(r<=xMax)) res+=DensitySamplingBuffer[0][i]*dx;

        for (idim=0;idim<3;idim++) x[idim]+=l[idim]*dx;
        IntegratedLength+=dx;
      }

      fprintf(fout,"%e  %e  %e\n",xMin+(0.5+ncell)*dx,(0.5+ncell)*dx,res);
    }

    fclose(fout);
  }



   MPI_Finalize();
   return 1;
 }



