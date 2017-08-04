//$Id$
//Interface to Khurana model of magnetic field in Saturn's magnetosphere

/*
 * KMAGInterface.cpp
 *
 *  Created on: Aug 1, 2017
 *      Author: vtenishe
 */


#include "KMAGInterface.h"

extern "C"{
  /*R,  THETA and PHI are radial distance, colatitude (in radians) and longitude (in radians)
  of the spacecraft. The program returns, local time (LT, in hours) and the model
  magnetic field in system 3 shperical coordinates.*/
  void kmag_(double* TIME,char* EPOCH,double* R,double* THETA,double* PHI,double* RLT,double* BY_IMF,double* BZ_IMF,double* BRM,double* BTM,double* BPM,double* Dp);

  /* INPUTS: From, a chracter*3 variable denoting the incoming coordinate system
     To:   a chracter*3 variable denoting the outgoing coordinate system
     Vecin a variable of dimension 3 containing the incoming vector
     Vecout: a variable of dimension 3 containing the outgoing vector
     time  A double precision variable denoting number of seconds from an epoch
     epoch   A five letter character which can be either 'ctime', 'etime' or 'J2000'/upper or lower case*/
  void krot_(char* From,char* To,double* Vecin,double* Vecout,double* time,char* epoch);
}


double KMAG::et=536500800.000000;  //the J2000 epoch - the number of seconds since 01/01/2000; The default number corresponds to 2017-01-01 00:00:00.000
char KMAG::Epoch[100]="J2000";
char KMAG::FRAME[100]="KSO";
double KMAG::SolarWind::IMF[3]={0.0,0.2,0.2};
double KMAG::SolarWind::DynamicPressure=0.02;

void KMAG::SetEpochTime(double etIn,const char* EpochIn) {
  et=etIn;
  sprintf(Epoch,EpochIn);
}

void KMAG::GetMagneticField(double *B,double *x) {
  double vecin[3],vecout[3]={0.0,0.0,0.0},r2=0.0;
  char SYS3FRAME[100]="S3C";
  int idim;

  //convert 'x' from meters to Saturn's radii
  for (idim=0;idim<3;idim++) {
    vecin[idim]=x[idim]/_SATURN__RADIUS_;
    r2+=pow(vecin[idim],2);
  }

  //if the point is inside the planet -> return B=0
  if (r2<=1.0) {
    for (idim=0;idim<3;idim++) B[idim]=0.0;
    return;
  }

  //convert the location from the solar orbiter to System III coordinate frame
  krot_(FRAME,SYS3FRAME,vecin,vecout,&et,Epoch);

  //calculate the magnetic fiels in System III coordinate fram
  double er,th,ph,rlt=0.0,brm=0.0,btm=0.0,bpm=0.0;

  er=sqrt(vecout[0]*vecout[0]+vecout[1]*vecout[1]+vecout[2]*vecout[2]);
  th=acos(vecout[2]/er);
  ph=atan2(vecout[1],vecout[0]);

  kmag_(&et,Epoch,&er,&th,&ph,&rlt,&SolarWind::IMF[1],&SolarWind::IMF[2],&brm,&btm,&bpm,&SolarWind::DynamicPressure);

  //convert the calculated magnetis field into System III frame
  double costh,sinth,cosph,sinph;

  costh=cos(th);
  sinth=sin(th);
  cosph=cos(ph);
  sinph=sin(ph);

  vecin[0]=brm*sinth*cosph+btm*costh*cosph-bpm*sinph;
  vecin[1]=brm*sinth*sinph+btm*costh*sinph+bpm*cosph;
  vecin[2]=brm*costh-btm*sinth;

  //rotate magnetic field vector from System III in KSO (solar orbiter)
  krot_(SYS3FRAME,FRAME,vecin,vecout,&et,Epoch);

  for (idim=0;idim<3;idim++) B[idim]=vecout[idim]*1e-9;
}

