//$Id$
//calculate trajecotry of Cassisni and orientation of UVIS instrument


//Important: description of the UVIS ports is in 'cas_uvis_v06.ti'

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


#include <sys/time.h>
#include <sys/resource.h>


#include "SpiceUsr.h"

#define      WDSIZE            320
#define Pi 3.141592654

const double TrajectoryExtractTimeIncrement=1;

//define the instrument code
#define _INMS_  0
#define _UVIS_  1

//define the UVIS observation ports
#define _UVIS_FUV_OCC_   -82842
#define _UVIS_SOLAR_OCC_ -82848

//the number of trajectories, definition of trajectories
struct cObservation {
  SpiceChar Tag[WDSIZE];
  SpiceChar Time[WDSIZE];
  double StartTimeIncrement,FinishTimeIncrement;
  int Instrument;
  int ObservationPort;
};


const int nTotalObservations=7;
cObservation ObservationList[nTotalObservations]={
	{"E2","2005-07-14T19:55:22",-1700,300,_INMS_,_INMS_},   //Waite-2006-science, Dong-2011-JGR
  {"E3","2008-03-12T19:06:12",-1600,1600,_INMS_,_INMS_},  //Dong-2011-JGR
  {"E5","2008-10-9T19:06:43",-1600,1600,_INMS_,_INMS_},  //Dong-2011-JGR
  {"E7","2009-11-2T7:42:00",-1600,1600,_INMS_,_INMS_},   //Dong-2011-JGR

  {"E2","2005-07-14T19:55:22",-400,600,_UVIS_,_UVIS_FUV_OCC_},  //Hansen-2006-science.pdf  (HSP and FUV data), Tian-2007-icarus
  {"24Oct2007","2007-10-24T16:59:50.3",-400,600,_UVIS_,_UVIS_FUV_OCC_},    //Hansen-2011-GRL.pdf,Hansen-2008-nature, Closest point of ray to star to limb (km) 15.64
  {"18May2010","2010-05-18T5:51:44.45",-400,900,_UVIS_,_UVIS_SOLAR_OCC_}};   //Hansen-2011-GRL.pdf (solar occultation)



struct cPlume {
  double Lat;
  double WLon;

  double x[3];
};

const int nTotPlumes=8;
const double Latitude[nTotPlumes]={-81.5,-79.2,-81.2,-73.2,-78.7,-87.1,-74.7,-82.1};
const double WLongitude[nTotPlumes]={31.2,313.2,294.2,148.4,72.6,237.0,28.9,115.5};
cPlume Plumes[nTotPlumes];


/*----------------------------------------------------------------------------------*/
void OrbitVariation() {

  //calculate the variation of Enceladus' orbit parameters
  SpiceDouble et,lt,etStart;
  SpiceDouble state0[6],state1[6];
  int nt;

  double v,r,vRadial,dt,v1,a;
  FILE *fout=fopen("Enceladus.Orbit.dat","w");

  fprintf(fout,"VARIABLES=\"et\", \"r\", \"v\", \"vRadial\", \"a\"\n");

  utc2et_c(ObservationList[0].Time,&et);
  etStart=et;

  for (nt=0;nt<500;nt++) {
    spkezr_c("Enceladus",et,"J2000","none","Saturn",state0,&lt);

    v=sqrt(pow(state0[3],2)+pow(state0[4],2)+pow(state0[5],2));
    dt=1500.0/v;

    spkezr_c("Enceladus",et+dt,"J2000","none","Saturn",state1,&lt);

    r=sqrt(pow(state0[0],2)+pow(state0[1],2)+pow(state0[2],2));
    vRadial=(state0[0]*state0[3]+state0[1]*state0[4]+state0[2]*state0[5])/r;

    v1=sqrt(pow(state1[3],2)+pow(state1[4],2)+pow(state1[5],2));
    a=fabs(v1-v)/dt;

    et+=100.0*3600.0/500.0;

    fprintf(fout,"%e  %e  %e  %e  %e  \n", et-etStart,r,v,vRadial,a);
  }

  fclose(fout);


  printf("VARIABLES=\"n\", \"r\", \"v\", \"vRadial\", \"a\"\n");

  for (int n=0;n<nTotalObservations;n++) {
    utc2et_c(ObservationList[n].Time,&et);

    spkezr_c("Enceladus",et,"J2000","none","Saturn",state0,&lt);

    v=sqrt(pow(state0[3],2)+pow(state0[4],2)+pow(state0[5],2));
    dt=1500.0/v;

    spkezr_c("Enceladus",et+dt,"J2000","none","Saturn",state1,&lt);

    r=sqrt(pow(state0[0],2)+pow(state0[1],2)+pow(state0[2],2));
    vRadial=(state0[0]*state0[3]+state0[1]*state0[4]+state0[2]*state0[5])/r;

    v1=sqrt(pow(state1[3],2)+pow(state1[4],2)+pow(state1[5],2));
    a=fabs(v1-v)/dt;

    et+=100.0*3600.0/500.0;

    printf("%i  %e  %e  %e  %e  \n", n,r,v,vRadial,a);
  }

}

/*----------------------------------------------------------------------------------*/
//set up individual plumes; get the positinos of the plumes
void InitIndividualPlume(int np,FILE *fout) {
  SpiceDouble x[3];
  SpiceInt      id,n;
  SpiceBoolean  found;

  bodn2c_c ( "ENCELADUS", &id, &found );

  //srfrec_c  converts planetocentric latitude and longitude of a surface
  //point on a specified body to rectangular coordinates.
  //planetocentic Longitude == east Longitude -> convert WLongitude to east Longitude before calling srfrec_c
  srfrec_c ( id, (360.0-WLongitude[np])*rpd_c(), Latitude[np]*rpd_c(),x);

  Plumes[np].Lat=Latitude[np];
  Plumes[np].WLon=WLongitude[np];
  Plumes[np].x[0]=x[0]*1.0E3,Plumes[np].x[1]=x[1]*1.0E3,Plumes[np].x[2]=x[2]*1.0E3;

  std::cout << "plume=" << np << ", x=" << Plumes[np].x[0] << "  " << Plumes[np].x[1] << "  " << Plumes[np].x[2] << std::endl;

  SpiceDouble radii[3],re,rp,f,spglon,spglat,spgalt;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html

  recpgr_c ( "ENCELADUS",  Plumes[np].x,  re,     f,&spglon, &spglat, &spgalt   );

  fprintf(fout,"ZONE T=\"Vent=%i\", I=1, J=1, F=POINT \n",np);
  fprintf(fout,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);
}

/*----------------------------------------------------------------------------------*/
//get trajectory and ground track of an flyby

void GetTrajectoryGroundTrack(int nObservation,FILE *fGroundTrack) {
  char fname[WDSIZE];
  FILE *fTrajectory;

  SpiceDouble et,lt,t,dEt,etStart;
  SpiceDouble state[6];
  SpiceDouble spoint[3];
  SpiceDouble trgepc;
  SpiceDouble srfvec[3];

  SpiceInt n;
  SpiceDouble radii[3],re,rp,f,spglon,spglat,spgalt,obspos[3],opglon,opglat,opgalt,r,rmin=-1.0;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html

  sprintf(fname,"Trajectory.%s.set=%i.dat",ObservationList[nObservation].Tag,nObservation);
  fTrajectory=fopen(fname,"w");
  fprintf(fTrajectory," \"UTC\", \"Time\", \t \"x\", \"y\", \"z\", \"r\" \t \"vx\", \"vy\", \"vz\", \"v\",\t \"l0\", \"l1\", \"l2\", \"ray of sight height\", \"Distance to the limb\"\n");


  int nt,nTimePoints=(ObservationList[nObservation].FinishTimeIncrement-ObservationList[nObservation].StartTimeIncrement)/TrajectoryExtractTimeIncrement;
  fprintf(fGroundTrack,"ZONE T=\"INMS %s: Ground Track\", I=%i, J=1, F=POINT \n",ObservationList[nObservation].Tag,nTimePoints);

  dEt=(ObservationList[nObservation].FinishTimeIncrement-ObservationList[nObservation].StartTimeIncrement)/(nTimePoints-1);
  utc2et_c(ObservationList[nObservation].Time,&etStart);

  for (nt=0;nt<nTimePoints;nt++) {
    et=etStart+ObservationList[nObservation].StartTimeIncrement+nt*dEt;
    spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);

    //get the sub-observer point
    subpnt_c("Intercept: ellipsoid","ENCELADUS",et,"IAU_ENCELADUS","none","Cassini",spoint,&trgepc,srfvec);
    recpgr_c ( "ENCELADUS",  spoint,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html
    vsub_c ( spoint, srfvec, obspos );
    recpgr_c ( "ENCELADUS",  obspos,  re,    f,&opglon, &opglat, &opgalt   );


    //fprintf(fGroundTrack,"%e  %e\n",360.0-spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);  ////sqrt(spoint[0]*spoint[0]+spoint[1]*spoint[1]));
    fprintf(fGroundTrack,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);

//      r=sqrt(pow(state[0],2)+pow(state[1],2)+pow(state[2],2))-252.1;
    r=opgalt;

    if ((rmin<0.0)||(r<rmin)) rmin=r;

    const int lenout=35;
    SpiceChar utcstr[100];

    et2utc_c(et,"C",6,lenout,utcstr);
    t=ObservationList[nObservation].StartTimeIncrement+nt*dEt;

    fprintf(fTrajectory,"%s ",utcstr);
    fprintf(fTrajectory,"%e \t %e  %e  %e  %e", t,state[0]*1.0E3,state[1]*1.0E3,state[2]*1.0E3,r*1.0E3);
    fprintf(fTrajectory,"\t %e  %e  %e  %e\n",state[3]*1.0E3,state[4]*1.0E3,state[5]*1.0E3,1.0E3*sqrt(state[3]*state[3]+state[4]*state[4]+state[5]*state[5]));
  }

  std::cout << std::endl << "INMS ground track: " << ObservationList[nObservation].Tag << std::endl;
  std::cout << "Minimum Altitude = " << rmin << std::endl;

  fclose(fTrajectory);
}



/*----------------------------------------------------------------------------------*/
void GetRayMinAltitude(SpiceDouble &spglon,SpiceDouble &spglat,SpiceDouble &spgalt,SpiceDouble *l,SpiceDouble* xmin,SpiceDouble et,int ObservationPortID) {
  SpiceDouble state[6],lt,rotate[3][3];
  int idim;

  #define  MAXBND 4
  SpiceChar    shape  [WDSIZE];
  SpiceChar    frame  [WDSIZE];
  SpiceDouble  bsight [3];
  SpiceInt     nbounds;
  SpiceDouble  bounds [MAXBND][3];

  SpiceInt n;
  SpiceDouble radii[3],re,rp,f;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html

  spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);

  getfov_c (ObservationPortID,MAXBND,WDSIZE, WDSIZE,shape, frame, bsight, &nbounds, bounds );
  pxform_c (frame,"IAU_ENCELADUS",et,rotate);
  mxv_c(rotate,bsight,l);

  //calculate the closest point between the line of sight and Enceladus
  SpiceDouble t=0.0;

  for (idim=0;idim<3;idim++) t-=state[idim]*l[idim]; //it shound be -=!!!
  for (idim=0;idim<3;idim++) xmin[idim]=state[idim]+t*l[idim];

  recpgr_c ( "ENCELADUS",  xmin,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html
}

double GetRayMinAltitude(SpiceDouble et,int ObservationPortID) {
  SpiceDouble spglon,spglat,spgalt,l[3],xmin[3];

  GetRayMinAltitude(spglon,spglat,spgalt,l,xmin,et,ObservationPortID);
  return spgalt;
}


double GetRayMinLimbDistance(SpiceDouble et,int ObservationPortID) {
  SpiceDouble state[6],lt,rotate[3][3];

  #define  MAXBND 4
  SpiceChar    shape  [WDSIZE];
  SpiceChar    frame  [WDSIZE];
  SpiceDouble  bsight [3],l[3];
  SpiceInt     nbounds;
  SpiceDouble  bounds [MAXBND][3];

  SpiceInt n;
  SpiceDouble radii[3],re,rp,f;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html

  spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);

  getfov_c (ObservationPortID,MAXBND,WDSIZE, WDSIZE,shape, frame, bsight, &nbounds, bounds );
  pxform_c (frame,"IAU_ENCELADUS",et,rotate);
  mxv_c(rotate,bsight,l);

  //calculate the closest point between the line of sight and the limb
  SpiceEllipse        limb;
  SpicePlane        plane;
  SpiceInt nxpts;
  SpiceDouble          xpt1[3],xpt2[3],vec1[3],vec2[3],sep1,sep2,limbClosestPoint;

  edlimb_c (radii[0],radii[1],radii[2], state, &limb );
  psv2pl_c ( state, state, l, &plane );
  inelpl_c ( &limb, &plane, &nxpts, xpt1, xpt2 );

  vsub_c   ( xpt1, state, vec1 );
  vsub_c   ( xpt2, state, vec2 );

  sep1 = vsep_c ( vec1, l );
  sep2 = vsep_c ( vec2, l );

  if (sep1<sep2) {
    limbClosestPoint=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])*sin(sep1);
  }
  else {
    limbClosestPoint=sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2])*sin(sep2);
  }
  return limbClosestPoint;
}
/*----------------------------------------------------------------------------------*/


//get the altitude of the line of sight
void GetLineOfSightAltitude(char* Time) {
  SpiceDouble l[3],et,spglon,spglat,spgalt,xmin[3];
  utc2et_c(Time,&et);

  GetRayMinAltitude(spglon,spglat,spgalt,l,xmin,et,_UVIS_FUV_OCC_);
  std::cout << "Minimal altitude of the line of sight = " << spgalt << " at UTC = " << Time << std::endl;
}

/*----------------------------------------------------------------------------------*/
//get pointing direction of UVIS, and the ground track of the line of sight

void GetLineOfSightGroundTrack(int nObservation,SpiceDouble et, FILE *fGroundTrack,double &minAltitude, double &minLimbDistance,int ObservationPortID) {
  SpiceDouble spglon,spglat,spgalt,l[3],xmin[3];

  const int nGroundTrackPoints=1000;
  const double LineOfSightLength=1.0E3;

  SpiceInt n;
  SpiceDouble radii[3],re,rp,f;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html


  SpiceDouble t,dt,x[3];
  int np;

  GetRayMinAltitude(spglon,spglat,spgalt,l,xmin,et,ObservationPortID);
  if ((minAltitude<0.0)||(minAltitude>spgalt)) minAltitude=spgalt;

  t=GetRayMinLimbDistance(et,ObservationPortID);
  if ((minLimbDistance<0.0)||(minLimbDistance>t)) minLimbDistance=t;

  dt=(ObservationList[nObservation].FinishTimeIncrement-ObservationList[nObservation].StartTimeIncrement)/(nGroundTrackPoints-1);

  //get the zone name
  const int lenout=35;
  SpiceChar utcstr[100];

  et2utc_c(et,"C",6,lenout,utcstr);
  fprintf(fGroundTrack,"ZONE T=\"UVIS %s: Ground Track, UTC=%s\", I=%i, J=1, F=POINT \n",ObservationList[nObservation].Tag,utcstr,nGroundTrackPoints);

  for (np=0;np<nGroundTrackPoints;np++) {
    for (int idim=0;idim<3;idim++) x[idim]=xmin[idim]+LineOfSightLength*l[idim]*(ObservationList[nObservation].StartTimeIncrement+np*dt);
    recpgr_c ( "ENCELADUS",  x,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html

    //fprintf(fGroundTrack,"%e  %e\n",/*360.0-*/spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);  ////sqrt(spoint[0]*spoint[0]+spoint[1]*spoint[1]));
    fprintf(fGroundTrack,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);
  }

}

/*----------------------------------------------------------------------------------*/
void TigerStripes() {

  struct cTigerStripe {
    char Tag[WDSIZE];
    double cStart[2];
    double cFinish[2];

    double xStart[3];
    double xFinish[3];

    double latStart;
    double wlonStart;
    double latFinish;
    double wlonFinish;

    int TigerStripeType;
  };

  const int nTotalTigerStripes=8;

  cTigerStripe TigerStripe[nTotalTigerStripes]={
      {"Alexandria",{18.28442,1.98199},{3.92776, 18.13063},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,0},
      {"Bagdad-Strong",{-2.64108,  2.65766},{5.28217, -8.96396},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,1},
      {"Bagdad",{8.60045,  -12.34234},{-13.61174, 14.27928},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,0},
      {"Cairo-Strong",{5.28217,  5.56307},{11.78330,  -6.73423},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,1},
      {"Cairo",{0.13544, 13.80631},{11.78330, -6.93693},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,0},
      {"Damascus-Strong",{-11.98646, 1.77928},{-4.33409,  -11.46396},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,1},
      {"Damascus",{-16.18510,  6.84685},{-2.09932,  -14.16666},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,0},
      {"Unknown",{4.15350,19.68469},{-2.48307,14.36937},{0.0,0.0,0.0},{0.0,0.0,0.0},0.0,0.0,0.0,0.0,0}
  };



  //calculate position of the start and end point of each stripe
  int nstripe,n;
  double lat,wlon,c;
  SpiceInt      id;
  SpiceBoolean  found;

  SpiceDouble radii[3],re,rp,f,spglon,spglat,spgalt;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html

  bodn2c_c ( "ENCELADUS", &id, &found );

  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    //start point
    c=sqrt(pow(TigerStripe[nstripe].cStart[0],2)+pow(TigerStripe[nstripe].cStart[1],2));
    lat=(-90.0+c)/180.0*Pi;

    wlon=acos(TigerStripe[nstripe].cStart[0]/c);
    if (TigerStripe[nstripe].cStart[1]<0.0) wlon=2.0*Pi-wlon;
    wlon+=Pi/2.0;
    if (wlon>2.0*Pi) wlon-=2.0*Pi;

    TigerStripe[nstripe].latStart=lat,TigerStripe[nstripe].wlonStart=wlon;
    srfrec_c ( id, 2.0*3.141592654-wlon,lat,TigerStripe[nstripe].xStart);


    //finish point
    c=sqrt(pow(TigerStripe[nstripe].cFinish[0],2)+pow(TigerStripe[nstripe].cFinish[1],2));
    lat=(-90.0+c)/180.0*Pi;

    wlon=acos(TigerStripe[nstripe].cFinish[0]/c);
    if (TigerStripe[nstripe].cFinish[1]<0.0) wlon=2.0*Pi-wlon;
    wlon+=Pi/2.0;
    if (wlon>2.0*Pi) wlon-=2.0*Pi;

    TigerStripe[nstripe].latFinish=lat,TigerStripe[nstripe].wlonFinish=wlon;
    //georec_c ( wlon, lat, 0.0, re, f, TigerStripe[nstripe].xFinish ); //Geodetic longitude==
    srfrec_c ( id, 2.0*3.141592654-wlon,lat,TigerStripe[nstripe].xFinish);
  }


  //output the stripes into a file
  FILE *fout;
  double x[3],l[3];
  int npoint,idim;

  const int nTotalPoints=200;

  fout=fopen("StripeGroundTrack.dat","w");
  fprintf(fout," VARIABLES=\"W Lon\", \"R Projection\"\n");

  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    fprintf(fout,"ZONE T=\"Stripe: %s\", I=%i, J=1, F=POINT \n",TigerStripe[nstripe].Tag,nTotalPoints);


    for (idim=0;idim<3;idim++) l[idim]=TigerStripe[nstripe].xFinish[idim]-TigerStripe[nstripe].xStart[idim];

    for (npoint=0;npoint<nTotalPoints;npoint++) {
      for (idim=0;idim<3;idim++) x[idim]=TigerStripe[nstripe].xStart[idim]+double(npoint)/double(nTotalPoints-1)*l[idim];

/*      reclat_c ( x, &spgalt, &spglon, &spglat );
      if (spglon<0.0) spglon+=2.0*Pi;
      spglon=2.0*Pi-spglon;*/

      recpgr_c ( "ENCELADUS",  x,  re,     f,&spglon, &spglat, &spgalt   );
      fprintf(fout,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);
    }
  }

  fclose(fout);


  //prepare the list of points (3D) along the stripes
  fout=fopen("TigerStripe3D.h","w");


  double dx=1.0; //in km

  fprintf(fout,"const int nTotalTigerStripes=%i;\n",nTotalTigerStripes);
  fprintf(fout,"const double dxTigerStripe=%e;\n",dx*1.0E3);

  fprintf(fout,"const int nTypesTigerStripe=2;\n");
  fprintf(fout,"double TigerStripeProductionRate[2]={0.0,0.0};\n");
  fprintf(fout,"double TigerStripeBulkVelocity=0.0,TigerStripeTemeprature=0.0;\n");

  //output the array of types of the stripes
  fprintf(fout,"const int TypeTigerStripe[%i]={",nTotalTigerStripes);
  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    fprintf(fout,"%i",TigerStripe[nstripe].TigerStripeType);

    if (nstripe!=nTotalTigerStripes-1) fprintf(fout,",");
    else fprintf(fout,"};\n");
  }


  //create array with the number of points for each stripe
  int maxStripePointNumber=0;

  fprintf(fout,"\nconst int StripePointListLength[%i]={",nTotalTigerStripes);

  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    for (idim=0;idim<3;idim++) l[idim]=TigerStripe[nstripe].xFinish[idim]-TigerStripe[nstripe].xStart[idim];
    int npoints=(int)(sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2])/dx)+1;

    if (maxStripePointNumber<npoints) maxStripePointNumber=npoints;

    fprintf(fout,"%i",npoints);

    if (nstripe!=nTotalTigerStripes-1) {
      fprintf(fout,",");
    }
    else {
      fprintf(fout,"};");
    }
  }

  //output points for each stripe
  fprintf(fout,"\nconst double StripePointList[%i][%i][3]={\n",nTotalTigerStripes,maxStripePointNumber);

  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    int np,npoints;

     //evaluate the number of points along the stripe
    for (idim=0;idim<3;idim++) l[idim]=TigerStripe[nstripe].xFinish[idim]-TigerStripe[nstripe].xStart[idim];
    npoints=(int)(sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2])/dx)+1;

    fprintf(fout,"\n{\n");

    for (np=0;np<maxStripePointNumber;np++) {
      fprintf(fout,"{");

      for (idim=0;idim<3;idim++) {
        if (np<npoints) fprintf(fout,"%e",1.0e3*(TigerStripe[nstripe].xStart[idim]+l[idim]*double(np)/double(npoints-1)));
        else fprintf(fout,"0.0");

        if (idim==2) {
          if (np==maxStripePointNumber-1) fprintf(fout,"}\n");
          else fprintf(fout,"},\n");
        }
        else fprintf(fout,",");
      }
    }

    if (nstripe!=nTotalTigerStripes-1) fprintf(fout,"},\n");
    else fprintf(fout,"}\n");
  }

  fprintf(fout,"};\n");
  fclose(fout);

  //output positions of the plumes

  fout=fopen("PointSource3D.h","w");
  fprintf(fout,"\nconst int nTotalPointSources=8;\n\nconst double xPointSources[8][3]={\n");

  for (int np=0;np<8;np++) {
    fprintf(fout,"{%e,%e,%e}",Plumes[np].x[0],Plumes[np].x[1],Plumes[np].x[2]);
    if (np!=7) fprintf(fout,",");
    fprintf(fout,"\n");
  }

  fprintf(fout,"};\n");
  fclose(fout);

}



/*----------------------------------------------------------------------------------*/
int main() {
  //furnish the kernels
  const char LocationSPICE[]="/Users/vtenishe/SPICE/Kernels/CASSINI/";

  const int nKernelsLoad=7+4+5+7+6+5;
  const char KernelList[nKernelsLoad][WDSIZE]={
     //???
   "100514BP_SCPSE_10134_10145.bsp", "10134_10139ra.bc", "cas_uvis_v06.ti", "cpck20May2010.tpc", "cpck20May2010_Nav.tpc",
   "cas_v40.tf", "cas00155.tsc",

   //E2: UVIS
   "050825R_SCPSE_05186_05205.bsp", "05192_05197ra.bc", "cpck05Jul2005.tpc", "cpck05Jul2005_Nav.tpc",

   //E3: 2008-03-12T19:06:12
   "08072_08077ra.bc", "080306AP_SE_08050_08076.bsp", "080228AP_SCPSE_08050_08076.bsp", "cpck10Mar2008.tpc", "cpck10Mar2008_Nav.tpc",

   //E5: 2008-10-9T19:06:43
   "08292_08297ra.bc", "081010AP_SCPSE_08282_08315.bsp", "081010AP_SE_08282_08315.bsp", "cpck31Oct2008.tpc", "cpck31Oct2008_Nav.tpc",
   "cpck16Sep2008.tpc", "cpck16Sep2008_Nav.tpc",

   //E7: 2009-11-2T7:42:00
   "09305_09307ra.bc", "091104AP_SCPSE_09305_09329.bsp", "cpck18Nov2009.tpc", "cpck18Nov2009_Nav.tpc", "cpck07Oct2009_Nav.tpc",
   "cpck07Oct2009.tpc",

   //UVIS occultation 24 Oct 2007
   "07297_07297ra.bc", "cpck18Oct2007.tpc", "cpck18Oct2007_Nav.tpc", "071015AP_SCPSE_07288_07328.bsp","07297_07302ra.bc"
  };



  furnsh_c("/Users/vtenishe/SPICE/Kernels/NAIF/naif0010.tls");

  for (int iKernel=0;iKernel<nKernelsLoad;iKernel++) {
    char FullName[WDSIZE];

    sprintf(FullName,"%s%s",LocationSPICE,KernelList[iKernel]);
    furnsh_c(FullName);
  }

  //init positions of individual plumes
  int i;
  FILE *fout;

  fout=fopen("VentLocation.dat","w");
  fprintf(fout," VARIABLES=\"W Lon\", \"R Projection\"\n");

  for (i=0;i<nTotPlumes;i++) InitIndividualPlume(i,fout);
  fclose(fout);

  //get the ground track of the observations
  FILE *fGroundTrack=fopen("TrajectoryGroundTrack.dat","w");
  fprintf(fGroundTrack," VARIABLES=\"W Lon\", \"R Projection\"\n");


  FILE *fPointingUVIS=fopen("PointingUVIS.dat","w");
  fprintf(fPointingUVIS,"VARIABLES=\"x\",\"y\",\"z\",  \"l0\",\"l1\",\"l2\"\n");

  for (i=0;i<nTotalObservations;i++) {
    if (ObservationList[i].Instrument==_INMS_) {
      GetTrajectoryGroundTrack(i,fGroundTrack);
    }

    if (ObservationList[i].Instrument==_UVIS_) {
      FILE *fGroundTrackUVIS;
      char fname[WDSIZE];
      SpiceDouble et;

      std::cout << "\nUVIS Ground Track: " << ObservationList[i].Tag << std::endl;

      sprintf(fname,"GroundTrackUVIS.%s.dat",ObservationList[i].Tag);
      fGroundTrackUVIS=fopen(fname,"w");
      fprintf(fGroundTrackUVIS,"VARIABLES=\"W Lon\", \"R Projection\"\n");

      //calcualte the ground track of the UVIS observations
      double dt,t,minAltitude=-1.0,minLimbDistance=-1.0,minAltitudeTime,minLimbDistanceTime;

      const int nObservationTracks=500;

      dt=(ObservationList[i].FinishTimeIncrement-ObservationList[i].StartTimeIncrement)/(nObservationTracks-1);
      utc2et_c(ObservationList[i].Time,&et);

      //calculate the pointing direction and position of the s/c
      SpiceDouble spglon,spglat,spgalt,l[3],xmin[3],state[6],lt;

      GetRayMinAltitude(spglon,spglat,spgalt,l,xmin,et,ObservationList[i].ObservationPort);
      spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);
      fprintf(fPointingUVIS,"%e %e %e        %e %e %e\n",state[0],state[1],state[2],l[0],l[1],l[2]);


      //calculate the ground tracks during the observation
      for (int j=0;j<nObservationTracks;j++) {
        GetLineOfSightGroundTrack(i,et+dt*j+ObservationList[i].StartTimeIncrement,fGroundTrackUVIS,minAltitude,minLimbDistance,ObservationList[i].ObservationPort);
      }

      //calculate the minimum value of the altitude and distance from the limb
      dt=(ObservationList[i].FinishTimeIncrement-ObservationList[i].StartTimeIncrement)/(nObservationTracks-1);
      minAltitude=-1.0,minLimbDistance=-1.0;

      for (int j=0;j<nObservationTracks;j++) {
        t=GetRayMinLimbDistance(et+dt*j+ObservationList[i].StartTimeIncrement,ObservationList[i].ObservationPort);
        if ((minLimbDistance<0.0)||(minLimbDistance>t)) minLimbDistance=t,minLimbDistanceTime=dt*j+ObservationList[i].StartTimeIncrement;

        t=GetRayMinAltitude(et+dt*j+ObservationList[i].StartTimeIncrement,ObservationList[i].ObservationPort);
        if ((minAltitude<0.0)||(minAltitude>t)) minAltitude=t,minAltitudeTime=dt*j+ObservationList[i].StartTimeIncrement;
      }

      std::cout << "Min Altitude = " << minAltitude << ", at " << minAltitudeTime << " sec" << std::endl;
      std::cout << "Min Distance from the limb = " << minLimbDistance << ", at " << minLimbDistanceTime << " sec" << std::endl;

      fclose(fGroundTrackUVIS);
    }
  }

  fclose(fGroundTrack);
  fclose(fPointingUVIS);

  //get position of the s/c at the following locations
  const int nTestTrajectoryPoints=5;

  char TestTrajectoryPoints[nTestTrajectoryPoints][WDSIZE]={
      "2005-07-14T19:53:35","2005-07-14T19:53:50","2005-07-14T19:54:10","2005-07-14T19:54:35","2005-07-14T19:54:50"
  };

  for (int i=0;i<nTestTrajectoryPoints;i++) {
    GetLineOfSightAltitude(TestTrajectoryPoints[i]);
  }


  //get the tracj of the stripes
  TigerStripes();

  //variation of Enceladus' orbit parameters
  OrbitVariation();

  return 1;
}

