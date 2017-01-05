//interface to Geopack 2008
//$Id$

#include <stdio.h>
#include <time.h>
#include <strings.h>


#include "GeopackInterface.h"

#include "constants.h"
#include "constants.PlanetaryData.h"
#include "ifileopr.h"

extern "C"{
  void recalc_08_(int*,int*,int*,int*,int*,double*,double*,double*);
  void sphcar_08_(double*,double*,double*,double*,double*,double*,int*);
  void bspcar_08_(double*,double*,double*,double*,double*,double*,double*,double*);
  void igrf_geo_08_(double*,double*,double*,double*,double*,double*);
}


void Geopack::Init(const char* Epoch,double *SolarWindVelocity) {
  CiFileOperations Parser;

  //conver the epoch string to the Geopack format
  int cnt=0;
  char tEpoch[300],str[100],*endptr;
  int Year,Month,Day,Hour,Minute,Second,DayOfYear;

  do {
    tEpoch[cnt]=Epoch[cnt];

    if ((tEpoch[cnt]=='T')||(tEpoch[cnt]=='-')||(tEpoch[cnt]==':')) tEpoch[cnt]=' ';
    cnt++;
  }
  while ((Epoch[cnt]!='\n')&&(Epoch[cnt]!=0));

  Parser.CutInputStr(str,tEpoch);
  Year=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Month=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Day=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Hour=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Minute=(int)strtol(str,&endptr,10);

  Parser.CutInputStr(str,tEpoch);
  Second=(int)strtol(str,&endptr,10);

  struct tm calendar = {0};
  calendar.tm_year = Year - 1900;
  calendar.tm_mon  = Month - 1;
  calendar.tm_mday = Day;
  time_t date = mktime ( &calendar );
  DayOfYear=calendar.tm_yday + 1;

  //init Geopack
  double VGSE[3];
  double defaultSolarWindVelocity[3]={-400.0,0.0,0.0};

  if (SolarWindVelocity!=NULL) {
    for (int idim=0;idim<3;idim++) VGSE[idim]=SolarWindVelocity[idim]/1.0E3;
  }
  else {
    for (int idim=0;idim<3;idim++) VGSE[idim]=defaultSolarWindVelocity[idim];
  }

  recalc_08_(&Year,&DayOfYear,&Hour,&Minute,&Second,VGSE+0,VGSE+1,VGSE+2);
}


void Geopack::IGRF::GetMagneticField(double *B,double *x) {
  /*
   *       SUBROUTINE IGRF_GEO_08 (R,THETA,PHI,BR,BTHETA,BPHI)
c
C  CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN THE SPHERICAL GEOGRAPHIC
C  (GEOCENTRIC) COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE MODEL
C  COEFFICIENTS  (e.g., http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html, revised: 22 March, 2005)
C
C  BEFORE THE FIRST CALL OF THIS SUBROUTINE, OR IF THE DATE (IYEAR AND IDAY) WAS CHANGED,
C  THE MODEL COEFFICIENTS SHOULD BE UPDATED BY CALLING THE SUBROUTINE RECALC_08
C
C-----INPUT PARAMETERS:
C
C   R, THETA, PHI - SPHERICAL GEOGRAPHIC (GEOCENTRIC) COORDINATES:
C   RADIAL DISTANCE R IN UNITS RE=6371.2 KM, COLATITUDE THETA AND LONGITUDE PHI IN RADIANS
C
C-----OUTPUT PARAMETERS:
C
C     BR, BTHETA, BPHI - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
C      (POSITIVE BR OUTWARD, BTHETA SOUTHWARD, BPHI EASTWARD)


    SUBROUTINE SPHCAR_08 (R,THETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICE VERSA
C    (THETA AND PHI IN RADIANS).
C
C                  J>0            J<0
C-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,THETA,PHI
C
C  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 WHEN CONVERTING
C        FROM CARTESIAN TO SPHERICAL COORDS (I.E., FOR J<0)
C
C   LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES AND MORE
C                         COMMENTS ADDED)
C
C   AUTHOR:  N. A. TSYGANENKO
C
   */

  int idim,J;
  double R,THETA,PHI,xLocal[3],BR,BTHETA,BPHI;

  for (idim=0;idim<3;idim++) xLocal[idim]=x[idim]/_EARTH__RADIUS_;

  //convert cartesian vector into the spherical coordinate system
  J=-1;
  sphcar_08_(&R,&THETA,&PHI,xLocal+0,xLocal+1,xLocal+3,&J);

  //extract the magnetic field vector
  igrf_geo_08_(&R,&THETA,&PHI,&BR,&BTHETA,&BPHI);

  //conver the mgnetic field vector to the cartesian coordinate frame
  /*       SUBROUTINE BSPCAR_08 (THETA,PHI,BR,BTHETA,BPHI,BX,BY,BZ)
C
C   CALCULATES CARTESIAN FIELD COMPONENTS FROM LOCAL SPHERICAL ONES
C
C-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
C              BR,BTHETA,BPHI -  LOCAL SPHERICAL COMPONENTS OF THE FIELD
C-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
C
C   LAST MOFIFICATION:  APRIL 1, 2003 (ONLY SOME NOTATION CHANGES)
C
C   WRITTEN BY:  N. A. TSYGANENKO
*/

  bspcar_08_(&THETA,&PHI,&BR,&BTHETA,&BPHI,B+0,B+1,B+2);

  for (idim=0;idim<3;idim++) B[idim]*=_NANO_;
}

//=======================================================================
//calculate the total magnetic field
void Geopack::GetMagneticField(double *B,double *x) {

}






















