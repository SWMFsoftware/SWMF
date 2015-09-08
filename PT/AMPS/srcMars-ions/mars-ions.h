
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
#include "constants.h"

#include "specfunc.h"
#include "ifileopr.h"

#include "mars-ions.dfn"

namespace MarsIon {
  using namespace Exosphere;

  /* _MARS_ION_H_ */
  // user defined global time step
  extern double UserGlobalTimeStep;

  //  injection boundary condition
  extern double InjectionVelocity[3];
  extern double InjectionNDensity;
  extern double InjectionTemperature;

  // computational domain size
  extern double DomainXMin[3];
  extern double DomainXMax[3];
  extern double DomainDXMin;
  extern double DomainDXMax;
  /* _MARS_ION_H_ */


  //init the model
  void Init_BeforeParser();
  inline void Init_AfterParser() {}

  //process interaction of the particles with the boundaries of the domain and the surface of the planet
  int ParticleSphereInteraction(int spec,long int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer);

  namespace SourceProcesses {
    double GetCellInjectionRate(int spec,double *xMiddle);
    double GetCellInjectionRate(int spec,PIC::Mesh::cDataCenterNode *cell);
    double GetBlockInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node);
  }

  namespace Output{

    extern int TotalDataLength;

    namespace OplusSource {
      extern int RateOffset;
      extern int BulkVelocityOffset;
      extern int TemperatureOffset;
    }

    void Init();

    void PrintVariableList(FILE* fout,int DataSetNumber);
    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);

    int RequestDataBuffer(int offset);
  }


  //init the background data
  void InitBackgroundData();
}




#endif /* _MARS_ION_H_ */
