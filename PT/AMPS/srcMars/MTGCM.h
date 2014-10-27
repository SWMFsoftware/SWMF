//==========================================================================
//$Id$
//==========================================================================
//the class provides reading and interpolation of data produced by MTGCM model

#ifndef _MTGCM_
#define _MTGCM_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include "specfunc.h"
#include "ifileopr.h"


#define _MTGCM_READER_STRING_LENGTH_ 2000

//the exeption handling modes
#define _MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__RETURN_ZERO_  0
#define _MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__ERROR_EXIT_   1

//interpolation modes at locations where holes in the data file are detected
#define _MTGCM_INTERPOLATION_MODE_UNDEFINED_                             -1
#define _MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT_                   0
#define _MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_                      1
#define _MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_   2

#define _TGCM_VALID_VALUE_LIMIT_ 1.0E35

#define Pi 3.14159265358979323846264338327950288419716939937510582

class cDataSetMTGCM {
public:
  char Label[_MTGCM_READER_STRING_LENGTH_],DataFileName[_MTGCM_READER_STRING_LENGTH_];
  long int nAltitudePoints,nLongitudePoints,nLatitudePoints;
  double minAltitude,maxAltitude,dAltitude;
  double *LongitudeCoordinates,*LatitudeCoordinates;
  double dLon,dLat,lonMin,lonMax,latMin,latMax;
  double ***DataBuffer;
  int nLonFileLines,nLatFileLines;

  double PlanetRadius;  //the radius of the planer
  int OutsideDomainInterpolationMode;

  //the exeption handler mode
  int ExeptionHandlerMode;

  //Warning flags
  bool NegativeScaleHightFound;
  bool PrintWarningMessage__DataSetNotLoaded;


  cDataSetMTGCM() {
    sprintf(Label,"\n");
    nAltitudePoints=0,nLongitudePoints=0,nLatitudePoints=0;
    minAltitude=0,maxAltitude=0,dAltitude=0;
    LongitudeCoordinates=NULL,LatitudeCoordinates=NULL;
    dLon=0.0,dLat=0.0,lonMin=0.0,lonMax=0.0,latMin=0.0,latMax=0.0;
    DataBuffer=NULL;
    nLonFileLines=0,nLatFileLines=0;
    PlanetRadius=-1.0;
    OutsideDomainInterpolationMode=_MTGCM_INTERPOLATION_MODE_UNDEFINED_;
    ExeptionHandlerMode=_MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__RETURN_ZERO_;

    NegativeScaleHightFound=false;
    PrintWarningMessage__DataSetNotLoaded=false;
  }

  void ReadDataFile(const char *fname) {
    CiFileOperations ifile;
    char *endptr,str[_MTGCM_READER_STRING_LENGTH_],str1[_MTGCM_READER_STRING_LENGTH_];
    double h,hLast=0.0,c,lonLast=0.0,latLast=0.0;
    int nx,ny,i,j,offset;

    sprintf(DataFileName,"%s",fname);

    //check initialization of the basic parameters:
    if (PlanetRadius<0.0) exit(__LINE__,__FILE__,"Error: the planet radiaus is not initialized");
    if (OutsideDomainInterpolationMode==_MTGCM_INTERPOLATION_MODE_UNDEFINED_) exit(__LINE__,__FILE__,"Error: the outside domain interpolation mode is not defined");

    ifile.openfile(fname);
    ifile.GetInputStr(Label,_MTGCM_READER_STRING_LENGTH_);

    minAltitude=-1,maxAltitude=-1,dAltitude=-1;

    //first pass through the file: load LongitudeCoordinates,LatitudeCoordinates; deternime nAltitudePoints
    while (ifile.eof()==false) {
      if (ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_)==false) break;
      ifile.CutInputStr(str1,str);

      if (strcmp("NX",str1)==0) {
        ifile.CutInputStr(str1,str);
        nx=strtol(str1,&endptr,10);

        if (nLongitudePoints==0) {
          nLongitudePoints=nx;
          ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
          nLonFileLines=1;

          LongitudeCoordinates=new double[nx];

          for (i=0;i<nLongitudePoints;i++) {
            if (strcmp("",str)==0) {
              ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
              nLonFileLines++;
            }

            ifile.CutInputStr(str1,str);
            c=strtod(str1,&endptr)/180.0*Pi;
            LongitudeCoordinates[i]=c;

            if (i==0) {
              lonMin=c;
            }
            else {
              if (i==1) dLon=c-lonLast;

              if (fabs(dLon-c+lonLast)>1.0E-5*dLon) exit(__LINE__,__FILE__,"Error: non uniform distribution of longitude points");
              if (i==nLongitudePoints-1) lonMax=c;
            }

            lonLast=c;
          }
        }
        else if (nLongitudePoints!=nx) exit(__LINE__,__FILE__,"Error: nx is not consistent");
      }
      else if (strcmp("NY",str1)==0) {
        ifile.CutInputStr(str1,str);
        ny=strtol(str1,&endptr,10);

        if (nLatitudePoints==0) {
          nLatitudePoints=ny;
          ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
          nLatFileLines=1;

          LatitudeCoordinates=new double[ny];

          for (i=0;i<nLatitudePoints;i++) {
            if (strcmp("",str)==0) {
              ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
              nLatFileLines++;
            }

            ifile.CutInputStr(str1,str);
            c=strtod(str1,&endptr)/180.0*Pi;
            LatitudeCoordinates[i]=c;

            if (i==0) {
              latMin=c;
            }
            else {
              if (i==1) dLat=c-latLast;

              if (fabs(dLat-c+latLast)>1.0E-5*dLat) exit(__LINE__,__FILE__,"Error: non uniform distribution of longitude points");
              if (i==nLatitudePoints-1) latMax=c;
            }

            latLast=c;
          }

        }
        else if (nLatitudePoints!=ny) exit(__LINE__,__FILE__,"Error: ny is not consistent");
      }
      else if (strcmp("UT",str1)==0) {
        ifile.CutInputStr(str1,str);
        ifile.CutInputStr(str1,str);
        ifile.CutInputStr(str1,str);

        h=strtod(str1,&endptr)*1.0E3;

        if ((h<minAltitude)||(minAltitude<0.0)) minAltitude=h;
        if ((h>maxAltitude)||(maxAltitude<0.0)) maxAltitude=h;

        if (nAltitudePoints==0) dAltitude=h,hLast=h;
        else if (nAltitudePoints==1) dAltitude=h-dAltitude,hLast=h;
        else {
          if (fabs(h-hLast-dAltitude)>1.0E-5*dAltitude) exit(__LINE__,__FILE__,"Error: dAlt is ont uniform");
          hLast=h;
        }

        nAltitudePoints++;
      }
    }


    //allocate the data buffer: DataBuffer[nAltitude][nLatitude][nLongitude];
    DataBuffer=new double** [nAltitudePoints];

    DataBuffer[0]=new double* [nAltitudePoints*nLatitudePoints];
    for (i=1;i<nAltitudePoints;i++) DataBuffer[i]=i*nLatitudePoints+DataBuffer[0];

    DataBuffer[0][0]=new double [nAltitudePoints*nLatitudePoints*nLongitudePoints];

    for (offset=0,i=0;i<nAltitudePoints;i++) for (j=0;j<nLatitudePoints;j++) {
      DataBuffer[i][j]=offset+DataBuffer[0][0];
      offset+=nLongitudePoints;
    }

    //rewind the file and read the data set
    int iAlt=0,iLat=0,iLon=0;

    ifile.rewindfile();

    while (ifile.eof()==false) {
      if (ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_)==false) break;
      ifile.CutInputStr(str1,str);

      if (strcmp("NY",str1)==0) {
        for (i=0;i<nLatFileLines;i++) ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

        ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

        for (iLat=0;iLat<nLatitudePoints;iLat++) for (iLon=0;iLon<nLongitudePoints;iLon++) {
          if (strcmp("",str)==0) ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

          ifile.CutInputStr(str1,str);
          DataBuffer[iAlt][iLat][iLon]=strtod(str1,&endptr);
        }

        ++iAlt;
      }
    }

    ifile.closefile();
  }


  //==========================================================================
  //determine wether the value of the point can be determined
  bool DataValueDefined(double *x) {
    double r;
    bool res;

    if ((nAltitudePoints==0)||(nLatitudePoints==0)||(nLongitudePoints==0)) {
      if (PrintWarningMessage__DataSetNotLoaded==false) {
        printf("WARNINIG: MTGCM interpolation module is accesssed before reading the dataset[file=%s, line=%i]\n",__FILE__,__LINE__);
        PrintWarningMessage__DataSetNotLoaded=true;
      }

      return false;
    }

    r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    res=(r-PlanetRadius-minAltitude>0.0) ? true : false;

    return res;
  }


  //==========================================================================
  //interpolate the stored value
  double Interpolate(double *x) {
    double res,Lon,Lat,r,h;
    int iAlt,iLon,iLat,idim;
    bool AboveDomainLimit=false;
    double cAlt,cLon,cLat; //interpolation coefficient

    for (r=0.0,idim=0;idim<3;idim++) r+=x[idim]*x[idim];
    r=sqrt(r);

    Lat=asin(x[2]/r);

    if (x[0]*x[0]+x[1]*x[1]<1.0E-20*x[2]*x[2]) Lon=0.0;
    else {
      Lon=acos(x[0]/sqrt(x[0]*x[0]+x[1]*x[1]));
      if (x[1]<0.0) Lon*=-1.0;
    }

    iAlt=(int)((r-PlanetRadius-minAltitude)/dAltitude);
    iLon=(int)((Lon-lonMin)/dLon);
    iLat=(int)((Lat-latMin)/dLat);

    cAlt=(r-PlanetRadius-minAltitude)/dAltitude-iAlt;
    cLon=(Lon-lonMin)/dLon-iLon;
    cLat=(Lat-latMin)/dLat-iLat;

    if (iLat<0) iLat=0,cLat=0.0;
    if (iLat>=nLatitudePoints-1) iLat=nLatitudePoints-2,cLat=1.0;

    if (iLon<0) iLon=0,cLon=0.0;
    if (iLon>=nLongitudePoints-1) iLon=nLongitudePoints-2,cLon=1.0;

    if (iAlt<0) {
      if (ExeptionHandlerMode==_MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__RETURN_ZERO_) return 0.0;
      exit(__LINE__,__FILE__,"Error: the point is below the lower boundary");
    }

    if (iAlt>=nAltitudePoints-1) AboveDomainLimit=true,iAlt=nAltitudePoints-1;
    else {
      for (int i=0;i<2;i++) for (int j=0;j<2;j++) if (DataBuffer[iAlt+1][iLat+i][iLon+j]>_TGCM_VALID_VALUE_LIMIT_) AboveDomainLimit=true;
    }


    if (AboveDomainLimit==false) {
      //!!!!! Linear interpolation with the altitude, Latitude and Longitude


cAlt=0.0,cLon=0.0,cLat=0.0;


      res=(1.0-cAlt)*((1.0-cLon)*(1.0-cLat)*DataBuffer[iAlt][iLat][iLon]+cLon*(1.0-cLat)*DataBuffer[iAlt][iLat][iLon+1]+(1.0-cLon)*cLat*DataBuffer[iAlt][iLat+1][iLon]+cLon*cLat*DataBuffer[iAlt][iLat+1][iLon+1]) +
          cAlt*((1.0-cLon)*(1.0-cLat)*DataBuffer[iAlt+1][iLat][iLon]+cLon*(1.0-cLat)*DataBuffer[iAlt+1][iLat][iLon+1]+(1.0-cLon)*cLat*DataBuffer[iAlt+1][iLat+1][iLon]+cLon*cLat*DataBuffer[iAlt+1][iLat+1][iLon+1]);

    }
    else { //the point is above the domain - use case dependent vertical interpolation
      int i,j,iAltInit=iAlt;
      double c=0.0,c1=0.0;

      for (res=0.0,i=0;i<2;i++) for (j=0;j<2;j++) {
        c=(i==0) ? (1.0-cLon) : cLon;
        c*=((j==0) ? (1.0-cLat) : cLat);

        iAlt=iAltInit;

        switch (OutsideDomainInterpolationMode) {
        case _MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT_: case _MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_:
          do {
            //estimate the scale hight
            for (;iAlt>0;iAlt--) if (DataBuffer[iAlt][iLat+j][iLon+i]<_TGCM_VALID_VALUE_LIMIT_) break;
            for (;iAlt>0;iAlt--) if ((DataBuffer[iAlt][iLat+j][iLon+i]>0.0)&&(DataBuffer[iAlt-1][iLat+j][iLon+i]>0.0)) break;

            h=(iAlt!=0) ? -1.0/dAltitude*log(DataBuffer[iAlt][iLat+j][iLon+i]/DataBuffer[iAlt-1][iLat+j][iLon+i]) : 0.0;

            if ((OutsideDomainInterpolationMode==_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_)&&(h<0.0)) {
              if (NegativeScaleHightFound==false) {
                NegativeScaleHightFound=true;
                printf("WARNING:  Negative value of the scale height is found by MTGCM interpolation procedure (filename is %s)\n",DataFileName);
              }

              if (--iAlt<=0) {
                if (ExeptionHandlerMode==_MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__RETURN_ZERO_) return 0.0;
                exit(__LINE__,__FILE__,"Error: iAlt is out of the range");
              }
            }
          }
          while ((OutsideDomainInterpolationMode==_MTGCM_INTERPOLATION_MODE_VERTICAL_SCALE_HIGHT__FORCE_POSITIVE_)&&(h<0.0));

          c1=DataBuffer[iAlt][iLat+j][iLon+i]*exp(-h*(r-PlanetRadius-minAltitude-iAlt*dAltitude));

          break;
        case _MTGCM_INTERPOLATION_MODE_VERTICAL_CONSTANT_:
          for (;iAlt>=0;iAlt--) if ((c1=DataBuffer[iAlt][iLat+j][iLon+i])<_TGCM_VALID_VALUE_LIMIT_) break;

          if (iAlt<0) {
            if (ExeptionHandlerMode==_MTGCM_INTERPOLATION_EXEPTION_HANDLING_MODE__RETURN_ZERO_) return 0.0;
            exit(__LINE__,__FILE__,"Error: iAlt is out of the range");
          }
          break;
        default:
          exit(__LINE__,__FILE__,"Error: option is not recognized");
        }

        res+=c*c1;
      }

    }

    return res;
  }

};


#endif
