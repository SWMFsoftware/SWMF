/*
 * pic_kameleon.cpp
 *
 *  Created on: Jan 1, 2015
 *      Author: vtenishe
 */
//$Id$
//AMPS' coupler functions that call CCMC's data reader Kameleon for interpolatino of the pre-calculated plasma model results stodered in the .cdf format

#include "pic.h"

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__CCMC_
#undef Pi //Kameleon uses its own cpnstant 'Pi'; The definition of the "global" 'Pi' is remove to eliminated the naming conflict with Kameleon library
#undef DIM

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <ccmc/LFM.h>
#include <ccmc/LFMInterpolator.h>
#include <ccmc/Interpolator.h>
#include <ccmc/Point3f.h>
#include <ccmc/Kameleon.h>
#include <ccmc/Tracer.h>
#include <string>
#include <fstream>
#include <iomanip>
#include "time.h"
#include <vector>

double PIC::CPLR::CCMC::PlasmaSpeciesAtomicMass=1.0*_AMU_;


void PIC::CPLR::CCMC::LFM::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static ccmc::LFM *lfm=NULL;
  static ccmc::LFMInterpolator *interpolator=NULL;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //allocate the interpolator and read the data file
    lfm=new ccmc::LFM;
    lfm->open(fname);

    /*
     * Load variables (pressure,density and e-field already loaded)
     * It is necessary to load the variables prior to interpolating
     */
    lfm->loadVariable("bx");
    lfm->loadVariable("by");
    lfm->loadVariable("bz");
    lfm->loadVariable("ux");
    lfm->loadVariable("uy");
    lfm->loadVariable("uz");

    interpolator=new ccmc::LFMInterpolator(lfm);
  }

  //symbolic names of the interpolated variables
  static const std::string bVariables[] = {"bx", "by", "bz"};
  static const std::string eVariables[] = {"ex", "ey", "ez"};
  static const std::string vVariables[] = {"ux", "uy", "uz"};

  //perform the interpolation loop
  int i,j,k,nd,idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;



  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double x[3],T,n,p;

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
      x[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
      x[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);

      //locate the cell
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
      offset=CenterNode->GetAssociatedDataBufferPointer();

      //save the interpolated values
      for (idim=0;idim<3;idim++) {
        *(idim+(double*)(offset+MagneticFieldOffset))=interpolator->interpolate(bVariables[idim],x[0],x[1],x[2]);
        *(idim+(double*)(offset+ElectricFieldOffset))=interpolator->interpolate(eVariables[idim],x[0],x[1],x[2]);
        *(idim+(double*)(offset+PlasmaBulkVelocityOffset))=interpolator->interpolate(vVariables[idim],x[0],x[1],x[2]);
      }

      p=interpolator->interpolate("p",x[0],x[1],x[2]);
      n=interpolator->interpolate("rho",x[0],x[1],x[2])/PlasmaSpeciesAtomicMass;
      T=p/(n*Kbol);

      *(double*)(offset+PlasmaPressureOffset)=p;
      *(double*)(offset+PlasmaNumberDensityOffset)=n;
      *(double*)(offset+PlasmaTemperatureOffset)=T;
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(fname,startNode->downNode[nDownNode]);
  }


  if (startNode==PIC::Mesh::mesh.rootTree) {
    //de-allocate the interpolation procedure
    delete lfm;
    delete interpolator;

    lfm=NULL,interpolator=NULL;
  }
}

#endif

