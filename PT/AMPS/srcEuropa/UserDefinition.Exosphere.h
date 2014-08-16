/*
 * UserDefinition.Exosphere.h
 *
 *  Created on: Aug 1, 2012
 *      Author: vtenishe
 */

//$Id$
//containes case specific user definition for the exospheric model

#ifndef USERDEFINITION_EXOSPHERE_H_
#define USERDEFINITION_EXOSPHERE_H_


//path to the SPICE Kernels directory

//SPICE Kernels to be loaded
//const int nFurnishedSPICEkernels=5+12+1;
//const char SPICE_Kernels[nFurnishedSPICEkernels][_MAX_STRING_LENGTH_PIC_]={"spk/msgr_de405_de423s.bsp","fk/msgr_dyn_v600.tf","../../NAIF/naif0010.tls","pck/pck00009_MSGR_v10.tpc","fk/msgr_v210.tf",
//    "ik/msgr_epps_v100.ti","ck/msgr20110413.bc","ck/msgr20110414.bc","ck/msgr20110415.bc","ck/msgr20110416.bc","ck/msgr20110417.bc","ck/msgr20110418.bc","ck/msgr20110419.bc","ck/msgr20110420.bc","ck/msgr20110421.bc",
//    "sclk/messenger_1486.tsc","spk/msgr_20040803_20140823_od266sc_0.bsp","../../OTHER/GSE.tf"};



//time stams for sampling remote column density observations
/*
const int nReferenceGroundBasedObservations=5; //25;
const char ReferenceGroundBasedObservationTime[nReferenceGroundBasedObservations][_MAX_STRING_LENGTH_PIC_]={
    "2008-05-18T00:00:00","2008-07-06T00:00:00","2008-11-07T00:00:00","2007-11-12T00:00:00","2007-06-03T00:00:00"};
*/


/*    "2001-05-25T00:00:00","2006-01-12T00:00:00","2005-12-18T00:00:00","2006-06-17T00:00:00","2003-10-04T00:00:00",
    "2006-10-21T00:00:00","2003-05-07T00:00:00","1999-04-27T00:00:00","1998-05-28T00:00:00","1990-12-10T00:00:00",
    "1990-12-04T00:00:00","1988-01-11T00:00:00","2000-06-05T00:00:00","2003-02-06T00:00:00","2002-12-16T00:00:00",
    "2002-08-21T00:00:00","2003-12-13T00:00:00","2006-06-12T00:00:00","2008-07-13T00:00:00","2006-11-09T00:00:00"};*/



#endif /* USERDEFINITION_EXOSPHERE_H_ */
