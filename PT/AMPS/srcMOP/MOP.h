//$Id$
//heander for modelingof the Jovian system radiation environemt



/*
 * MOP.h
 *
 *  Created on: Jun 21, 2017
 *      Author: vtenishe
 */

#ifndef _SRC_MOP_H_
#define _SRC_MOP_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "pic.h"
#include "Dust.h"
#include "Exosphere.h"
#include "constants.h"

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
#include "SpiceUsr.h"
#else
#include "SpiceEmptyDefinitions.h"
#endif

#include "MOP.dfn"
#include "KMAGInterface.h"


namespace MOP {

  //Saturnian system
  namespace SaturninanSystem {




    //parameters that characterize the planet
    namespace Saturn {
      extern double RotationAxis[3];

      void GetMagneticFieldDipole(double *B,double *x);
    }

    //parameters that characterize the magnetosphere
    namespace Magnetosphere {
      //get magnetic files in the Saturn's magnetosphere
      void GetMagneticField(double *B,double *x);
      void GetElectricField(double *E,double*x);

      //get the corrotation speed of the magnetosphere
      double GetCorotationSpeed(double *x);

      //get the ionization lifetime of the neutral species
      namespace Lifetime {
        double Get(int spec,double *x);
      }
    }

    //Enceladus is the main source of water in the Jovian system
    namespace Enceladus {
      //location and velocity of Encelsdus
      extern double xEnceladus[3],vEnceladus[3];

      //the total souce rate
      extern double TotalSourceRate,SourceTemperature;

      //model of the volatiles injection from Enceladus
      //PIC::ParticleWeightTimeStep::fUserDefinedExtraSourceRate PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate
      //PIC::BC::UserDefinedParticleInjectionFunction=
      long int InjectParticles();
      double SourceRate(int spec);
    }

    //the requested local mesh resolution accounted for the topology of the planet system
    double GetLocalMeshResolution(double *x);
  }

  //Jupiter system
  namespace JovianSystem {

    namespace Jupiter {

    }

    namespace Europa {

    }
  }

  //sampling and output of the model data
  namespace Sampling {
    namespace Output {
     //function that are used for sample data output
    //print out of the output file
    void PrintVariableList(FILE* fout,int DataSetNumber);
    void Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode);
    void PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode);
    }
  }

  //interface for MOP
  double SourceRate(int spec);
  long int InjectParticles();

  //the requested local mesh resolution accounted for the topology of the planet system
  double GetLocalMeshResolution(double *x);

  //init the model
  void Init();
}



#endif /* _SRC_MOP_H_ */
