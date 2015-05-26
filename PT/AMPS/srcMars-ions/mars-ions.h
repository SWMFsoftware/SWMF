
//$Id$
//header of the model of the minor ion distribution aroud Mars

/*
 * mars-ion.h
 *
 *  Created on: May 26, 2015
 *      Author: vtenishe
 */

#ifndef _MARS_ION_H_
#define _MARS_ION_H_

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
#include <dirent.h>

#include "pic.h"
#include "Exosphere.h"

#include "specfunc.h"
#include "ifileopr.h"

#include "mars-ions.dfn"

namespace MarsIon {
  using namespace Exosphere;


  //init the model
  inline void Init_BeforeParser() {}
  inline void Init_AfterParser() {}

  //process interaction of the particles with the boundaries of the domain and the surface of the planet
  int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);

}




#endif /* _MARS_ION_H_ */
