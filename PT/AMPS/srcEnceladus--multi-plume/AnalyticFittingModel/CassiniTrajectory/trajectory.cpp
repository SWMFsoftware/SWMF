//$Id$
//calculate trajecotry of Cassisni and orientation of UVIS instrument


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


const double TrajectoryExtractTimeRange=600,TrajectoryExtractTimeIncrement=1;

#define _INMS_  0
#define _UVIS_  1

//the number of trajectories, definition of trajectories
struct cObservation {
  SpiceChar Tag[WDSIZE];
  SpiceChar Time[WDSIZE];
  int StartTimeIncrement,FinishTimeIncrement;
  int Instrument;
};


const int nTotalObservations=7;
cObservation ObservationList[nTotalObservations]={
	{"E2","2005-07-14T19:55:22",-1700,300,_INMS_},   //Waite-2006-science, Dong-2011-JGR
    {"E3","2008-03-12T19:06:12",-300,600,_INMS_},  //Dong-2011-JGR
    {"E5","2008-10-9T19:06:43",-300,600,_INMS_},  //Dong-2011-JGR
    {"E7","2009-11-2T7:42:00",-20,40,_INMS_},   //Dong-2011-JGR

    {"E2","2005-07-14T19:55:22",-400e3,600e3,_UVIS_},
    {"24Oct2007","2007-10-24T16:59:50.3",-400e3,600e3,_UVIS_},
    {"18May2010","2010-05-18T5:51:44.45",-400e3,900e3,_UVIS_}};



struct cPlume {
  double Lat;
  double WLon;

  double x[3];
};

const int nTotPlumes=8;
const double Latitude[nTotPlumes]={-81.5,-79.2,-81.2,-73.2,-78.7,-87.1,-74.7,-82.1};
const double WLongitude[nTotPlumes]={31.2,313.2,294.2,148.4,72.6,237.0,28.9,115.5};
cPlume Plumes[nTotPlumes];

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
//    georec_c ( wlon, lat, 0.0, re, f, TigerStripe[nstripe].xStart );
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

/*
  double c=1.3;

  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    TigerStripe[nstripe].xStart[0]=c*TigerStripe[nstripe].xStartProjection[0];
    TigerStripe[nstripe].xStart[1]=c*TigerStripe[nstripe].xStartProjection[1];
    TigerStripe[nstripe].xStart[2]=-sqrt(pow(c*250.0,2)-pow(TigerStripe[nstripe].xStart[0],2)-pow(TigerStripe[nstripe].xStart[1],2));

    TigerStripe[nstripe].xFinish[0]=c*TigerStripe[nstripe].xFinishProjection[0];
    TigerStripe[nstripe].xFinish[1]=c*TigerStripe[nstripe].xFinishProjection[1];
    TigerStripe[nstripe].xFinish[2]=-sqrt(pow(c*250.0,2)-pow(TigerStripe[nstripe].xFinish[0],2)-pow(TigerStripe[nstripe].xFinish[1],2));
  }

*/


  //output the stripes into a file
  FILE *fout;
  double x[3],l[3];
  int npoint,idim;
  char fname[200];


  const int nTotalPoints=50;




  for (nstripe=0;nstripe<nTotalTigerStripes;nstripe++) {
    sprintf(fname,"%s.GroundTrack.dat",TigerStripe[nstripe].Tag);
    fout=fopen(fname,"w");
    fprintf(fout," VARIABLES=\"W Lon\", \"R Projection\"\n");


    //rotate initial and final points of the stripe
    //srfrec_c  converts planetocentric latitude and longitude of a surface
    //point on a specified body to rectangular coordinates.
    //planetocentic Longitude == east Longitude -> convert WLongitude to east Longitude before calling srfrec_c

/*    recpgr_c ( "ENCELADUS",  TigerStripe[nstripe].xStart,  re,     f,&spglon, &spglat, &spgalt   );
    spglon=3.141592654/2.0-spglon;
    srfrec_c ( id, 2.0*3.141592654-spglon,spglat,TigerStripe[nstripe].xStart);

    recpgr_c ( "ENCELADUS",  TigerStripe[nstripe].xFinish,  re,     f,&spglon, &spglat, &spgalt   );
    spglon=3.141592654/2.0-spglon;
    srfrec_c ( id, 2.0*3.141592654-spglon,spglat,TigerStripe[nstripe].xFinish);*/


    for (idim=0;idim<3;idim++) l[idim]=TigerStripe[nstripe].xFinish[idim]-TigerStripe[nstripe].xStart[idim];

    for (npoint=0;npoint<nTotalPoints;npoint++) {
      for (idim=0;idim<3;idim++) x[idim]=TigerStripe[nstripe].xStart[idim]+double(npoint)/double(nTotalPoints-1)*l[idim];

      //get latitude and longitude
//      recpgr_c ( "ENCELADUS",  x,  re,     f,&spglon, &spglat, &spgalt   );



      reclat_c ( x, &spgalt, &spglon, &spglat );
      if (spglon<0.0) spglon+=2.0*Pi;
      spglon=2.0*Pi-spglon;




      fprintf(fout,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);
    }

    fclose(fout);
  }

  //output locations of the plumes
  for (npoint=0;npoint<nTotPlumes;npoint++) {
    sprintf(fname,"Source=%i.GroundTrack.dat",npoint+1);
    fout=fopen(fname,"w");
    fprintf(fout," VARIABLES=\"W Lon\", \"R Projection\"\n");

    recpgr_c ( "ENCELADUS",  Plumes[npoint].x,  re,     f,&spglon, &spglat, &spgalt   );
    fprintf(fout,"%e  %e\n",spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);

    fclose(fout);
  }


  //prepare the list of points (3D) along the stripes
  fout=fopen("TigerStripe3D.h","w");


  double dx=10.0; //in km

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
    fprintf(fout,"{%e,%e,%e}",1.0e3*Plumes[np].x[0],1.0e3*Plumes[np].x[1],1.0e3*Plumes[np].x[2]);
    if (np!=7) fprintf(fout,",");
    fprintf(fout,"\n");
  }

  fprintf(fout,"};\n");
  fclose(fout);

}

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


  //process the observation's settings
  SpiceInt n,nObservation;

  SpiceDouble et,lt,t;
  SpiceDouble state[6];
  SpiceDouble spoint[3];
  SpiceDouble trgepc;
  SpiceDouble srfvec[3];
  SpiceDouble rotate[3][3];

  #define  MAXBND 4
  SpiceChar    shape  [WDSIZE];
  SpiceChar    frame  [WDSIZE];
  SpiceDouble  bsight [3],uvisPointing[3],rsight_hight,rsight_limb;
  SpiceInt     nbounds;
  SpiceDouble  bounds [MAXBND][3];


  SpiceDouble radii[3],re,rp,f,spglon,spglat,spgalt,obspos[3],opglon,opglat,opgalt;
  bodvrd_c ("ENCELADUS","RADII",3,&n,radii);
  re=radii[0];
  rp=radii[2];
  f=(re-rp)/re;  //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/subpnt_c.html


  //calcualte the positions of the plums
  for (int np=0;np<nTotPlumes;np++) {
    SpiceDouble x[3];
    SpiceInt      id;
    SpiceBoolean  found;

    bodn2c_c ( "ENCELADUS", &id, &found );

    //srfrec_c  converts planetocentric latitude and longitude of a surface
    //point on a specified body to rectangular coordinates.
    //planetocentic Longitude == east Longitude -> convert WLongitude to east Longitude before calling srfrec_c
    srfrec_c ( id, (360.0-WLongitude[np])*rpd_c(), Latitude[np]*rpd_c(),x);

    std::cout << "plume=" << np << ", x=" << x[0]*1.0E3 << "  " << x[1]*1.0E3 << "  " << x[2]*1.0E3 << std::endl;

    Plumes[np].Lat=Latitude[np];
    Plumes[np].WLon=WLongitude[np];
    Plumes[np].x[0]=x[0],Plumes[np].x[1]=x[1],Plumes[np].x[2]=x[2];
  }

  std::cout << "\n\n";

  //tiger stripes
  TigerStripes();

  //calcualte the position of the s/c during E2 at -70 and -40 sec from CA
  utc2et_c(ObservationList[0].Time,&et);

  subpnt_c("Intercept: ellipsoid","ENCELADUS",et-70.0,"IAU_ENCELADUS","none","Cassini",spoint,&trgepc,srfvec);
  recpgr_c ( "ENCELADUS",  spoint,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html
  std::cout << "E2 (t=-70 sec): lat=" << spglat*180.0/3.141592654 << ", lon=" << /*360.0-*/spglon*180.0/3.141592654 << std::endl;

  subpnt_c("Intercept: ellipsoid","ENCELADUS",et-40.0,"IAU_ENCELADUS","none","Cassini",spoint,&trgepc,srfvec);
  recpgr_c ( "ENCELADUS",  spoint,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html
  std::cout << "E2 (t=-40 sec): lat=" << spglat*180.0/3.141592654 << ", lon=" << /*360.0-*/spglon*180.0/3.141592654 << std::endl;


  //calculate the ground track of the line of sight of UVIS
    for (int nUvisGrountTrack=0;nUvisGrountTrack<nTotalGrountTrackUVIS;nUvisGrountTrack++) {
      FILE *fUvisGroundTrack;
      char fname[200];

      for (int nTimeIncrement=0;nTimeIncrement<nGroundTrackUvisTimeIncrements;nTimeIncrement++) {
        sprintf(fname,"%s.UvisGroundTrack.TimeIncrement=%i.dat",GroundTrackUVIS[nUvisGrountTrack].Tag,nTimeIncrement);
        fUvisGroundTrack=fopen(fname,"w");
        fprintf(fUvisGroundTrack," VARIABLES=\"W Lon\", \"R Projection\", \"Altitude\"\n");

        utc2et_c(GroundTrackUVIS[nUvisGrountTrack].Time,&et);
        et+=GroundTrackUVIS[nUvisGrountTrack].TimeIncrementPonts[nTimeIncrement];

        spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);

        //get pointing of UVIS
        getfov_c (-82845,MAXBND,WDSIZE, WDSIZE,shape, frame, bsight, &nbounds, bounds );

        pxform_c (frame,"CASSINI_SC_COORD",et,rotate);
        pxform_c ("CASSINI_SC_COORD","IAU_ENCELADUS",et,rotate);

    //        pxform_c (frame,"IAU_ENCELADUS",et,rotate);

        pxform_c ("CASSINI_UVIS_FUV","IAU_ENCELADUS",et,rotate);


        mxv_c(rotate,bsight,uvisPointing);

        //get the distance between the position of the s/c to the point of the closest approach of the line of sight
        double d,minr=-1.0,minr_first=-1.0,r;

        d=-(state[0]*uvisPointing[0]+state[1]*uvisPointing[1]+state[2]*uvisPointing[2]);

        //go through the points on the line of sight
        for (double dd=-1000.0;dd<1000.0;dd+=1.0) {
          double x[3];
          int idim;

          for (idim=0,r=0.0;idim<3;idim++) {
            x[idim]=state[idim]+(d+dd)*uvisPointing[idim];
            r+=pow(x[idim],2);
          }

          r=sqrt(r)-250.0;

          if (minr_first<0.0) minr_first=r;
          if ((minr<0.0)||(minr>r)) minr=r;

          //find the latitude and longitude of the point
          recpgr_c ( "ENCELADUS",  x,  re,     f,&spglon, &spglat, &spgalt   );
          fprintf(fUvisGroundTrack,"%e  %e  %e\n",/*360.0-*/spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654,r);
        }

        fclose(fUvisGroundTrack);

        if ((minr>=minr_first)||(minr>=r)) {
          std::cout << "WARNINIG: for UVIS ground track " << fname << ", Time interval=" << nTimeIncrement << " the minimum altitude point is at the boundary of the track" << std::endl;
        }

        std::cout << "UVIS ground track "<< fname << ", Time interval=" << nTimeIncrement << ": the minimum altitude of the line of sight is " << minr << std::endl;
      }
    }


  //calculate parameters of Cassini's trajectory
  for (nFlyby=0;nFlyby<nTotalFlybyes;nFlyby++) {
    FILE *fout, *fCassiniTrack;
    char fname[200];
    double r,rmin=-1.0;
    double lmin_uvis=-1.0,lmin_uvis_time=-1.0;

    sprintf(fname,"%s.ExtractedOrbit.dat",FlybyList[nFlyby].Tag);
    fout=fopen(fname,"w");
    fprintf(fout," \"UTC\", \"Time\", \t \"x\", \"y\", \"z\", \"r\" \t \"vx\", \"vy\", \"vz\", \"v\",\t \"l0\", \"l1\", \"l2\", \"ray of sight height\", \"Distance to the limb\"\n");


    sprintf(fname,"%s.CassiniTrack.dat",FlybyList[nFlyby].Tag);
    fCassiniTrack=fopen(fname,"w");
    fprintf(fCassiniTrack," VARIABLES=\"W Lon\", \"R Projection\"\n");

    utc2et_c(FlybyList[nFlyby].Time,&et);

    for (t=FlybyList[nFlyby].StartTimeIncrement,et+=FlybyList[nFlyby].StartTimeIncrement;t<=FlybyList[nFlyby].FinishTimeIncrement;t+=TrajectoryExtractTimeIncrement,et+=TrajectoryExtractTimeIncrement) {
      spkezr_c("Cassini",et,"IAU_ENCELADUS","none","ENCELADUS",state,&lt);

      //get the sub-observer point
      subpnt_c("Intercept: ellipsoid","ENCELADUS",et,"IAU_ENCELADUS","none","Cassini",spoint,&trgepc,srfvec);
      recpgr_c ( "ENCELADUS",  spoint,  re,     f,&spglon, &spglat, &spgalt   ); //http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/recpgr_c.html
      vsub_c ( spoint, srfvec, obspos );
      recpgr_c ( "ENCELADUS",  obspos,  re,    f,&opglon, &opglat, &opgalt   );


      fprintf(fCassiniTrack,"%e  %e\n",/*360.0-*/spglon*180.0/3.141592654,90.0+spglat*180.0/3.141592654);  ////sqrt(spoint[0]*spoint[0]+spoint[1]*spoint[1]));

//      r=sqrt(pow(state[0],2)+pow(state[1],2)+pow(state[2],2))-252.1;
      r=opgalt;

      if ((rmin<0.0)||(r<rmin)) rmin=r;

      const int lenout=35;
      SpiceChar utcstr[100];

      et2utc_c(et,"C",6,lenout,utcstr);

      fprintf(fout,"%s ",utcstr);
      fprintf(fout,"%e \t %e  %e  %e  %e", t,state[0]*1.0E3,state[1]*1.0E3,state[2]*1.0E3,r*1.0E3);
      fprintf(fout,"\t %e  %e  %e  %e",state[3]*1.0E3,state[4]*1.0E3,state[5]*1.0E3,1.0E3*sqrt(state[3]*state[3]+state[4]*state[4]+state[5]*state[5]));

      if (FlybyList[nFlyby].uvis==true) {
        //get the field of view of the instrument
        getfov_c (-82845,MAXBND,WDSIZE, WDSIZE,shape, frame, bsight, &nbounds, bounds );

        pxform_c (frame,"CASSINI_SC_COORD",et,rotate);
        pxform_c ("CASSINI_SC_COORD","IAU_ENCELADUS",et,rotate);

//        pxform_c (frame,"IAU_ENCELADUS",et,rotate);

        pxform_c ("CASSINI_UVIS_FUV","IAU_ENCELADUS",et,rotate);


        mxv_c(rotate,bsight,uvisPointing);

        //calculate the height the line of sight
        double d;


        d=state[0]*uvisPointing[0]+state[1]*uvisPointing[1]+state[2]*uvisPointing[2];
        rsight_hight=sqrt(pow(state[0]-d*uvisPointing[0],2)+pow(state[1]-d*uvisPointing[1],2)+pow(state[2]-d*uvisPointing[2],2))-252.1;

        SpiceDouble pnear[3];
        npedln_c(radii[0],radii[1],radii[2],state, uvisPointing,pnear,&rsight_hight);
        if ((lmin_uvis<0.0)||(lmin_uvis>rsight_hight)) lmin_uvis=rsight_hight,lmin_uvis_time=t;


        //calculate the minimum distance between the line of sight from the limb
        SpiceEllipse        limb;
        SpicePlane        plane;
        SpiceInt nxpts;
        SpiceDouble          xpt1[3],xpt2[3],vec1[3],vec2[3],sep1,sep2;

        edlimb_c (radii[0],radii[1],radii[2], state, &limb );
        psv2pl_c ( state, state, uvisPointing, &plane );
        inelpl_c ( &limb, &plane, &nxpts, xpt1, xpt2 );

        vsub_c   ( xpt1, state, vec1 );
        vsub_c   ( xpt2, state, vec2 );

        sep1 = vsep_c ( vec1, uvisPointing );
        sep2 = vsep_c ( vec2, uvisPointing );

        if (sep1<sep2) {
          rsight_limb=sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]+vec1[2]*vec1[2])*sin(sep1);
        }
        else {
          rsight_limb=sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]+vec2[2]*vec2[2])*sin(sep2);
        }




      }
      else uvisPointing[0]=0.0,uvisPointing[1]=0.0,uvisPointing[2]=0.0,rsight_hight=0.0,rsight_limb=0.0;


      fprintf(fout,"\t %e  %e  %e  %e  %e\n",uvisPointing[0],uvisPointing[1],uvisPointing[2],rsight_hight,rsight_limb);
    }

    std::cout << "Flyby:" << FlybyList[nFlyby].Tag <<  std::endl;
    std::cout << "\t Min Altitude=" << rmin << " km" << std::endl;
    if (FlybyList[nFlyby].uvis==true) std::cout << "\t Min distance of the line of sight from the limb=" << lmin_uvis << " km, time=" << lmin_uvis_time << std::endl;

    fclose(fout);
    fclose(fCassiniTrack);
  }


  return 1;
}





