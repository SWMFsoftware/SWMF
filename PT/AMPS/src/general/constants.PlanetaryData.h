//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================
//$Id$
//=======================================
//the file contains information for differerent objects in the Solar System

#ifndef CONSTANTS_PLANETARYDATA_H_
#define CONSTANTS_PLANETARYDATA_H_

#define _AU_ 149598000.0E3

/*--------------------------   MARS  -----------------------------*/
#define _MARS__RADIUS_ 3.396E6

//#define _MARS__RADIUS_ 3402000.0

#define _MARS__MASS_ 6.4185E23
/*--------------------------   END MARS  -----------------------------*/


/*--------------------------     TITAN   -----------------------------*/
#define _TITAN__RADIUS_ 2.575E6
#define _TITAN__MASS_ 1.35E23
/*--------------------------  END TITAN  -----------------------------*/

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

/*--------------------------   MOON  -----------------------------*/
#define _MOON__RADIUS_    1737.10E3
#define _MOON__MASS_      7.3477E22
/*--------------------------  END MOON  -----------------------------*/


/*--------------------------   EARTH  -----------------------------*/
#define _EARTH__RADIUS_    6371.0E3
#define _EARTH__MASS_      5.9736E24
/*--------------------------  END EARTH  -----------------------------*/


/*--------------------------   EUROPA  -----------------------------*/
#define _EUROPA__RADIUS_    1.569E6
#define _EUROPA__MASS_      4.8E22
/*--------------------------  END EUROPA  -----------------------------*/

/*--------------------------   SUN --------------------------------------*/
#define _SUN__MASS_ 1.9891E30

/*--------------------------   END SUN --------------------------------------*/

/*-------------------------- COMET CHURYUMOV-GERASIMENCO ----------------*/
#define _CG__MASS_ 1.01E13

//#define _CG__RADIUS_ 1.98E3
#define _CG__RADIUS_ 2.0E3

#endif /* CONSTANTS_PLANETARYDATA_H_ */
