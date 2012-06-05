//=======================================
//$Id$
//=======================================
//the file contains information for differerent objects in the Solar System

#ifndef CONSTANTS_PLANETARYDATA_H_
#define CONSTANTS_PLANETARYDATA_H_

#define _AU_ 149598000.0E3

/*--------------------------   MARS  -----------------------------*/
#define _MARS__RADIUS_ 3.396E6
#define _MARS__MASS_ 6.4185E23
/*--------------------------   END MARS  -----------------------------*/


/*--------------------------   ENCELADUS  -----------------------------*/
#define _ENCELADUS__RADIUS_  252.1E3
#define _ENCELADUS__MASS_    1.08022E20

/*--------------------------   END ENCELADUS  -----------------------------*/


/*--------------------------   MERCURY  -----------------------------*/
#define _MERCURY__MASS_ 3.3022E23
#define _MERCURY__RADIUS_ 2439.7E3
#define _MERCUTY__APHELION_ (0.466697*_AU_)
#define _MERCURY__PERIHELION_ (0.307499*_AU_)
#define _MERCURY__SIDEREAL_ORBITAL_PERIOD_ (87.969*24.0*3600.0)
#define _MERCURY__SIDERIAL_ROTATIONAL_PERIOD_ (1407.6*24.0*3600.0)

/*--------------------------   END MERCURY  -----------------------------*/


/*--------------------------   SUN --------------------------------------*/
#define _SUN__MASS_ 1.9891E30

/*--------------------------   END SUN --------------------------------------*/

/*-------------------------- COMET CHURYUMOV-GERASIMENCO ----------------*/
#define _CG__MASS_ 1.01E13

//#define _CG__RADIUS_ 1.98E3
#define _CG__RADIUS_ 2.0E3

#endif /* CONSTANTS_PLANETARYDATA_H_ */
