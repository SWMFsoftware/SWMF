//$Id$
//a simple 2d-axysymmertic model of the distribution of the neutrals and ions in the vicinity of the Rosetta spacecraft 

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

const int nZcells=200;
const int nRcells=200;
const double rMin=1.0,rMax=50.0,rKill=500.0;
const double Temp=200.0;
const int nTimeFrames=1;

const double InjectionTimeInterval=30.0*60.0;
const double maxSystemTime=2.0*3600.0;
const double dtMax=1.0E-4;
const double TimeFrameInterval=maxSystemTime/nTimeFrames;

const double maxPhiIncrement=2.0*Pi/100.0; 
const double maxRrelativeIncrement=0.01;

//double ****SampleCounter=NULL,****buffer=NULL; //////[TimeFrame][nZcells][nRcells][2];
double SampleCounter[nZcells][nRcells][2][nTimeFrames];
double SampleCounterClose[nZcells][nRcells][2][nTimeFrames];
static const double dR=rMax/nRcells,dZ=2.0*rMax/nZcells;

const double LifeTime=1.0E4;
const int nTotalParticles=100000;

const double xMinPannel=1.0;
const double xMaxPannel=20.0;

const double mh2o=2.99151E-26;
const double Kbol=1.3806503E-23;
const double beta=sqrt(mh2o/(2.0*Kbol*Temp));

const double U=1.0*1.602176565E-19;

unsigned long int rndLastSeed;

double rnd() {
  rndLastSeed*=48828125;
  rndLastSeed&=2147483647; // pow(2,31) - 1
  if (rndLastSeed==0) rndLastSeed=1;

  return double(rndLastSeed/2147483648.0); //(pow(2,31) - 1) + 1
}


void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);

  if (seed==-1) seed=thread;

  rndLastSeed=seed;
}



int main(int argc, char **argv) {
  int thread;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);

   rnd_seed(-1);

/*
  //init the counter [TimeFrame][nZcells][nRcells][2]; 
  SampleCounter=new double*** [nTimeFrames];
  SampleCounter[0]=new double** [nTimeFrames*nZcells];
  SampleCounter[0][0]=new double* [nTimeFrames*nZcells*nRcells];
  SampleCounter[0][0][0]=new double [nTimeFrames*nZcells*nRcell 
*/




  double *t,*t1;
  int i;
  for (i=0,t=&SampleCounter[0][0][0][0];i<nZcells*nRcells*2*nTimeFrames;i++,t++) *t=0.0; 
  for (i=0,t=&SampleCounterClose[0][0][0][0];i<nZcells*nRcells*2*nTimeFrames;i++,t++) *t=0.0;
 
  //loop over particles
  double x[3],v[3],SystemTime,IonizationTime;
  double r,phi,dr,dphi;

  int npart;
  double c,c1,M,M0[3],E,dt;
  double xOrbitPlane,yOrbitPlane;

  double e0[3],e1[3],l;
  int idim,iR,iZ,iTimeFrame;

  double RadialMotionMultiplier=0.0;
  double AngularVelocityPrev,AngularVelocityNow;

  for (npart=0;npart<nTotalParticles;npart++) {
   //the first particle is injected at the moment when evaporation from the panels started, the last particle - when the evaporation is finished
   SystemTime=double(npart)/double(nTotalParticles)*InjectionTimeInterval; 


if (thread==0) {
  if (npart%100==0) std::cout << npart << std::endl;
}



   //generate particle velocity 
do {
   do {
     //normal component
     v[2]=sqrt(-log(rnd()))/beta;

     //tangential component 
     c=sqrt(-log(rnd()))/beta;
     c1=rnd();

     v[0]=c*sin(2.0*Pi*c1);
     v[1]=c*cos(2.0*Pi*c1);

     //location of the particle
     r=sqrt(rnd()*(xMaxPannel*xMaxPannel-xMinPannel*xMinPannel)+xMinPannel*xMinPannel);      
   
     //the start time of the paticle trajectory integration
     do {
       IonizationTime=-LifeTime*log(rnd());
    }
     while (IonizationTime+SystemTime>maxSystemTime);
   }
  while ((v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*IonizationTime*IonizationTime>rMax*rMax);

  //initial position of the particle
  SystemTime+=IonizationTime;
  r=pow((x[0]=r+v[0]*IonizationTime),2);
  r+=pow((x[1]=v[1]*IonizationTime),2)+pow((x[2]=v[2]*IonizationTime),2);
  r=sqrt(r),phi=0.0;

  //angular momentum
  M0[0]=x[1]*v[2]-v[1]*x[2];
  M0[1]=-(x[0]*v[2]-v[0]*x[2]);
  M0[2]=x[0]*v[1]-v[0]*x[1];

  M=mh2o*sqrt(pow(M0[0],2)+pow(M0[1],2)+pow(M0[2],2)); 
  E=0.5*mh2o*(pow(v[0],2)+pow(v[1],2)+pow(v[2],2))-U*rMin/r;  

  //integrate the particle equation of motion. The partice will stay in the same plane. get the coordinate frame related to that frame; 
  for (idim=0,l=0.0;idim<3;idim++) {
    e0[idim]=x[idim];
    l+=pow(e0[idim],2);
  }
  
  for (idim=0,c=0.0,l=sqrt(l);idim<3;idim++) {
    e0[idim]/=l;
    c+=e0[idim]*v[idim];
  }
}
while (c==0.0); //generate only thoth particles that has a non-zero radial velocity components 
RadialMotionMultiplier=c/fabs(c);


  for (idim=0,l=0.0;idim<3;idim++) {
    e1[idim]=v[idim]-c*e0[idim];
    l+=pow(e1[idim],2);
  }

  for (idim=0,l=sqrt(l);idim<3;idim++) e1[idim]/=l;

  dt=dtMax; 
  AngularVelocityPrev=M/(mh2o*pow(r,2));

int nLoopCounter=0;

  while ((rMin<r)&&(r<rKill)&&(SystemTime<maxSystemTime)) {
    nLoopCounter++;

    do {
      AngularVelocityNow=M/(mh2o*pow(r,2));
      dphi=AngularVelocityNow*dt;     

 
double a0,a1,a2;

a0=2.0/mh2o*(E+U*rMin/r);
a1=pow(M/(mh2o*r),2);
a2=a0-a1;

      dr=sqrt(2.0/mh2o*(E+U*rMin/r)-pow(M/(mh2o*r),2))*dt;

      if ((fabs(dphi)>maxPhiIncrement)||(dr>maxRrelativeIncrement*r)) dt/=2.0; 
      else break;
    } while (true);


    //determine the sign of dr 
    if (nLoopCounter!=1) {
      //with increase of the distance the radial velocity is decreasing: the first -> particle is approaching the s/c, the latter -> the particle is going away

      if (AngularVelocityNow>AngularVelocityPrev) {
        RadialMotionMultiplier=-1.0;
      }
      else {
        RadialMotionMultiplier=1.0;
      }
    }


    AngularVelocityPrev=AngularVelocityNow; 

    //update the location of the particle
    dr*=RadialMotionMultiplier;
    r+=dr;
    phi+=dphi;

    //convert the new location into the original frame of reference
    xOrbitPlane=r*sin(phi);
    yOrbitPlane=r*cos(phi);

    for (idim=0;idim<3;idim++) {
      x[idim]=xOrbitPlane*e0[idim]+yOrbitPlane*e1[idim];
    }

    //determine location of the particle on the sampling mesh
    iZ=(x[2]+rMax)/dZ;
    iR=sqrt(x[0]*x[0]+x[1]*x[1])/dR;
    iTimeFrame=SystemTime/TimeFrameInterval; 
    if (iTimeFrame>=nTimeFrames) iTimeFrame=nTimeFrames-1; 

    if ((r<rMin)||(r>rKill)) {
      //go to the next particle
      break;
    }

    if ((iZ>=0)&&(iZ<nZcells)&&(iR<nRcells)) {
      //the point is with in domain
      SampleCounter[iZ][iR][0][iTimeFrame]+=dt/dtMax;
      if (dr<0.0) SampleCounter[iZ][iR][1][iTimeFrame]+=dt/dtMax; 
    }

    //sample into the 'close' samepl counter
    iZ=(x[2]+0.1*rMax)/(0.1*dZ);
    iR=sqrt(x[0]*x[0]+x[1]*x[1])/(0.1*dR); 

    if ((iZ>=0)&&(iZ<nZcells)&&(iR<nRcells)) {
      //the point is with in domain
      SampleCounterClose[iZ][iR][0][iTimeFrame]+=dt/dtMax;
      if (dr<0.0) SampleCounterClose[iZ][iR][1][iTimeFrame]+=dt/dtMax;
    }

    //increment the time counter
    SystemTime+=dt;
  }
  }

  //summ all counters
  double buffer[nZcells][nRcells][2][nTimeFrames];

  MPI_Reduce(&SampleCounter[0][0][0][0],&buffer[0][0][0][0],2*nZcells*nRcells*nTimeFrames,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for (i=0,t=&SampleCounter[0][0][0][0],t1=&buffer[0][0][0][0];i<nZcells*nRcells*2*nTimeFrames;i++,t++,t1++) *t=*t1;    

  MPI_Reduce(&SampleCounterClose[0][0][0][0],&buffer[0][0][0][0],2*nZcells*nRcells*nTimeFrames,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for (i=0,t=&SampleCounterClose[0][0][0][0],t1=&buffer[0][0][0][0];i<nZcells*nRcells*2*nTimeFrames;i++,t++,t1++) *t=*t1;

  //output the data file 
  FILE *fout;
  
  if (thread==0) { 
//    char fname[100];
//    sprintf(fname,"data.frame=%03d.dat",iTimeFrame);

    fout=fopen("data.dat","w");

    for (iTimeFrame=0;iTimeFrame<nTimeFrames;iTimeFrame++) {
      fprintf(fout,"VARIABLES=\"r\", \"z\", \"all\", \"to the s/c\"\n");
      fprintf(fout,"ZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nRcells,nZcells); 

      for (iZ=0;iZ<nZcells;iZ++) for (iR=0;iR<nRcells;iR++) {
        r=(0.5+iR)*dR;

        fprintf(fout,"%e %e %e %e\n",r,(0.5+iZ)*dZ-rMax,SampleCounter[iZ][iR][0][iTimeFrame]/r,SampleCounter[iZ][iR][1][iTimeFrame]/r);   
      }
    }

    fclose(fout); 

 //   sprintf(fname,"data-close.frame=%03d.dat",iTimeFrame);
    fout=fopen("data-close.dat","w");
    
    for (iTimeFrame=0;iTimeFrame<nTimeFrames;iTimeFrame++) {
      fprintf(fout,"VARIABLES=\"r\", \"z\", \"all\", \"to the s/c\"\n");
      fprintf(fout,"ZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nRcells,nZcells);

      for (iZ=0;iZ<nZcells;iZ++) for (iR=0;iR<nRcells;iR++) {
        r=(0.5+iR)*0.1*dR;

        fprintf(fout,"%e %e %e %e\n",r,(0.5+iZ)*0.1*dZ-0.1*rMax,SampleCounterClose[iZ][iR][0][iTimeFrame]/r,SampleCounterClose[iZ][iR][1][iTimeFrame]/r);
      }
    }

    fclose(fout);
  }

  MPI_Finalize();
  return 1.0;
}  
