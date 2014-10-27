/*
 * UserDefinition.Exosphere.h
 *
 *  Created on: Jun 23, 2012
 *      Author: vtenishe
 */

//$Id$
//containes case specific user definition for the exospheric model

#ifndef USERDEFINITION_EXOSPHERE_H_
#define USERDEFINITION_EXOSPHERE_H_

//path to the SPICE Kernels directory
const char SPICE_Kernels_PATH[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/SPICE/Kernels";

//SPICE Kernels to be loaded
const int nFurnishedSPICEkernels=6;
const char SPICE_Kernels[nFurnishedSPICEkernels][_MAX_STRING_LENGTH_PIC_]={"/MESSENGER/kernels/spk/msgr_de405_de423s.bsp","NAIF/naif0010.tls","OTHER/Moon.LSO.tf","OTHER/GSE.tf",
    "MESSENGER/kernels/fk/msgr_dyn_v600.tf","/MESSENGER/kernels/pck/pck00009_MSGR_v10.tpc"};


//time stams for sampling remote column density observations
const int nReferenceGroundBasedObservations=8; //25;
const char ReferenceGroundBasedObservationTime[nReferenceGroundBasedObservations][_MAX_STRING_LENGTH_PIC_]={
    "1988-05-28T12:00:00","1988-10-02T12:00:00","1994-04-22T12:00:00","1993-11-29T12:00:00","1996-04-03T12:00:00",
    "1996-09-27T12:00:00","1997-03-24T12:00:00","2002-07-16T12:00:00"};


/*
   "2001-05-25T00:00:00","2006-01-12T00:00:00","2005-12-18T00:00:00","2006-06-17T00:00:00","2003-10-04T00:00:00",
    "2006-10-21T00:00:00","2003-05-07T00:00:00","1999-04-27T00:00:00","1998-05-28T00:00:00","1990-12-10T00:00:00",
    "1990-12-04T00:00:00","1988-01-11T00:00:00","2000-06-05T00:00:00","2003-02-06T00:00:00","2002-12-16T00:00:00",
    "2002-08-21T00:00:00","2003-12-13T00:00:00","2006-06-12T00:00:00","2008-07-13T00:00:00","2006-11-09T00:00:00"};
    */




#undef  _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_
#define _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_            _EXOSPHERE_SOURCE__OFF_


#undef _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_
#define _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_     _EXOSPHERE_SOURCE__ON_

#undef _EXOSPHERE_SOURCE__THERMAL_DESORPTION_
#define _EXOSPHERE_SOURCE__THERMAL_DESORPTION_               _EXOSPHERE_SOURCE__ON_


//redefine the value of the constant surface density
//the constant surface sodium density
//===================  Begin: Set Constant Surface Content (Killen-2012-JGR, Yakshinskiy&Madey-1999-?  ===============
#undef _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_
#define _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(el) (3.0E16)

#undef _EXOSPHERE__SURFACE_CONTENT_
//#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__BALANCE_FLUXES_
#define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_

#endif /* USERDEFINITION_EXOSPHERE_H_ */
