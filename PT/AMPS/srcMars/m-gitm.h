/*
 * m-gitm.h
 *
 *  Created on: May 10, 2013
 *      Author: vtenishe
 */
//$Id$


#ifndef _M_GITM_H_
#define _M_GITM_H_

#include "MTGCM.h"

//the order of data variables in the M-GITM output file
#define _Tn_MGITM_    3
#define _Ti_MGITM_    4
#define _Te_MGITM_    5
#define _nCO2_MGITM_  6
#define _nO_MGITM_    7
#define _nN2_MGITM_   8
#define _nCO_MGITM_   9
#define _nO2_MGITM_   10
#define _nO2P_MGITM_  11
#define _nOP_MGITM_   12
#define _nCO2P_MGITM_ 13
#define _Ne_MGITM_    14
#define _UN_MGITM_    15
#define _VN_MGITM_    16
#define _WN_MGITM_    17

class cDataSetMGITM : public cDataSetMTGCM {
public:
    int VariableID;

    cDataSetMGITM() : cDataSetMTGCM() {
		VariableID=-1;
	}

	void ReadDataFile(int nvar,const char* fname) {
	  CiFileOperations ifile;
	  char *endptr,str[_MTGCM_READER_STRING_LENGTH_],str1[_MTGCM_READER_STRING_LENGTH_];

	  sprintf(DataFileName,"%s",fname);
	  VariableID=nvar;

	  //check initialization of the basic parameters:
	  if (PlanetRadius<0.0) exit(__LINE__,__FILE__,"Error: the planet radiaus is not initialized");
	  if (OutsideDomainInterpolationMode==_MTGCM_INTERPOLATION_MODE_UNDEFINED_) exit(__LINE__,__FILE__,"Error: the outside domain interpolation mode is not defined");

	  ifile.openfile(fname);

	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  nLongitudePoints=strtol(str1,&endptr,10);

	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  ifile.CutInputStr(str1,str);
	  nAltitudePoints=strtol(str1,&endptr,10);

	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);
	  ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

	  //first pass thought the file: determine the range of variation of the Lon, Lat and altitude
	  minAltitude=-1,maxAltitude=-1,dAltitude=-1;
	  lonMin=-100.0,lonMax=-1.0,latMin=-100.0,latMax=-1.0;
	  nLatitudePoints=0;

      double t;

	  while (ifile.eof()==false) {
		  nLatitudePoints++;
		  if (ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_)==false) break;

		  ifile.CutInputStr(str1,str);
		  t=strtod(str1,&endptr)/180.0*Pi;
//		  if (t>Pi) t-=2.0*Pi;


		  if (t>lonMax) lonMax=t;
		  if ((lonMin<-10.0)||(lonMin>t)) lonMin=t;

		  ifile.CutInputStr(str1,str);
		  t=strtod(str1,&endptr)/180.0*Pi;
		  if (t>latMax) latMax=t;
		  if ((latMin<-10.0)||(latMin>t)) latMin=t;

		 // ifile.CutInputStr(str1,str);
		  ifile.CutInputStr(str1,str);
		  t=strtod(str1,&endptr)*1.0E3;
		  if (t>maxAltitude) maxAltitude=t;
		  if ((minAltitude<0.0)||(minAltitude>t)) minAltitude=t;
	  }

	  nLatitudePoints/=nLongitudePoints*nAltitudePoints;

	  dLon=(lonMax-lonMin)/(nLongitudePoints-1);
	  dLat=(latMax-latMin)/(nLatitudePoints-1);
	  dAltitude=(maxAltitude-minAltitude)/(nAltitudePoints-1);

	  //allocate the data buffer: DataBuffer[nAltitude][nLatitude][nLongitude];
	  int i,j,offset;
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
	  for (i=0;i<8;i++) ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);


	  for (iAlt=0;iAlt<nAltitudePoints;iAlt++) for (iLat=0;iLat<nLatitudePoints;iLat++) for (iLon=0;iLon<nLongitudePoints;iLon++) {
		ifile.GetInputStr(str,_MTGCM_READER_STRING_LENGTH_);

		for (i=0;i<=VariableID;i++) ifile.CutInputStr(str1,str);
		DataBuffer[iAlt][iLat][iLon]=strtod(str1,&endptr);
	  }

	  ifile.closefile();
	}
};


#endif /* M_GITM_H_ */
