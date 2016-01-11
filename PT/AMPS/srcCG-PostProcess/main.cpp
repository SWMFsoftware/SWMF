//$Id$

/*
 * main.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include "PostProcess3D.h"
#include "SpiceUsr.h"


#define _DUST_MODE_  0
#define _GAS_MODE_   1

#define _MODE_ _GAS_MODE_
//#define _MODE_  _DUST_MODE_

cPostProcess3D amps;

//process the dust file
struct cColumnIntegrationSet {
  int (*IntegrantVectorLength)();
  void (*IntegrantVector)(double* Data,double *Location);
  void (*PrintVariableList)(FILE* fout);
};




//=================================================
//processing of the dust output file
namespace GAS {
  int IntegrantVectorLength() {return 3;}
  void PrintVariableList(FILE* fout) {fprintf(fout," \"n\", \"Background Spes 1\", \"Background Spec 2\"");}

  void IntegrantVector(double *data,double *x) {
    cPostProcess3D::cStencil Stencil;
    int i;

    amps.GetInterpolationStencil(x,&Stencil);

    for (i=0;i<IntegrantVectorLength();i++) data[i]=0.0;

    for (i=0;i<8;i++) {
      data[0]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][6];

      data[1]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][38];
      data[2]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][42];
    }
  }
}

namespace DUST {
  int nRadii=10;

  struct cVariablePair {
    int nVar;
    double GrainRadius;
  };

  cVariablePair VariablePair[]={
      {51,0.5*(1.000000E-07+1.995262E-07)},
      {56,0.5*(1.995262E-07+3.981072E-07)},
      {61,0.5*(3.981072E-07+7.943282E-07)},
      {66,0.5*(7.943282E-07+1.584893E-06)},
      {71,0.5*(1.584893E-06+3.162278E-06)},
      {76,0.5*(3.162278E-06+6.309573E-06)},
      {81,0.5*(6.309573E-06+1.258925E-05)},
      {86,0.5*(1.258925E-05+2.511886E-05)},
      {91,0.5*(2.511886E-05+5.011872E-05)},
      {96,0.5*(5.011872E-05+1.000000E-04)}
  };

  int IntegrantVectorLength() {return 1+nRadii;}

  void PrintVariableList(FILE* fout) {
    fprintf(fout," \"Column Densty\"");

    for (int i=0;i<nRadii;i++) fprintf(fout,", \"Column Density(%e)\"",VariablePair[i].GrainRadius);
  }

  void IntegrantVector(double *data,double *x) {
    cPostProcess3D::cStencil Stencil;
    int i,iRadius;

    amps.GetInterpolationStencil(x,&Stencil);

    for (i=0;i<IntegrantVectorLength();i++) data[i]=0.0;

    for (i=0;i<8;i++) data[0]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][46];

    for (iRadius=0;iRadius<nRadii;iRadius++) for (i=0;i<8;i++) {
      data[1+iRadius]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
    }
  }

}


int main(int argc,char **argv) {

  const char SimulationStartTimeString[]="2015-04-12T07:14:00";

  //init SPICE
  int i;
  SpiceDouble StateRosetta[6],StateSun[6],et,lt;
  double xObservation[3]={1.0E6,0,0},xPrimary[3]={0,0,0},xSecondary[3]={0,1.0E6,0};

  furnsh_c("kernels.tm");
  utc2et_c(SimulationStartTimeString,&et);
  spkezr_c("ROSETTA",et,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateRosetta,&lt);
  spkezr_c("SUN",et,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateSun,&lt);


  for (i=0;i<3;i++) {
    xObservation[i]=1.0E3*StateRosetta[i];
    xPrimary[i]=0.0;
    xSecondary[i]=1.0E3*StateSun[i];
  }


  amps.InitMPI();
  amps.SetBlockSize(5,5,5);

  if (_MODE_==_GAS_MODE_) {
    amps.LoadDataFile("pic.H2O.s=0.out=10.dat",".");
  }
  else if (_MODE_==_DUST_MODE_) {
    amps.LoadDataFile("pic.DUST%0.s=2.out=10.dat",".");
  }


  amps.PrintVariableList();

  //load the nucleus mesh
  CutCell::ReadNastranSurfaceMeshLongFormat("SHAP5_stefano.bdf","/Volumes/data/AMPS_DATA/ROSETTA");


//  for (int n=0;n<amps.data.nNodes;n++) amps.data.data[n][38]=sqrt(pow(amps.data.data[n][0],2)+pow(amps.data.data[n][1],2)+pow(amps.data.data[n][2],2));




  //process the gas output file
  cPostProcess3D::cColumnIntegral::cColumnIntegrationSet IntegrationSet;

  if (_MODE_==_GAS_MODE_) {
    IntegrationSet.IntegrantVector=GAS::IntegrantVector;
    IntegrationSet.IntegrantVectorLength=GAS::IntegrantVectorLength;
    IntegrationSet.PrintVariableList=GAS::PrintVariableList;
  }
  else {
    IntegrationSet.IntegrantVector=DUST::IntegrantVector;
    IntegrationSet.IntegrantVectorLength=DUST::IntegrantVectorLength;
    IntegrationSet.PrintVariableList=DUST::PrintVariableList;
  }

  amps.ColumnIntegral.Map.Circular(xObservation,xPrimary,xSecondary,1.5,100,100,"map.dat",&IntegrationSet);
  amps.SaveDataFile("Res.dat", amps.data);


  amps.FinalizeMPI();
  return 1;
}


