/*
 * main.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include "PostProcess3D.h"
#include "SpiceUsr.h"

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

void ProcessGasOutputFile() {

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



  amps.SetBlockSize(5,5,5);
  amps.LoadDataFile("pic.H2O.s=0.out=10.dat",".");
  amps.PrintVariableList();

  //load the nucleus mesh
  CutCell::ReadNastranSurfaceMeshLongFormat("SHAP5_stefano.bdf","/Volumes/data/AMPS_DATA/ROSETTA");


//  for (int n=0;n<amps.data.nNodes;n++) amps.data.data[n][38]=sqrt(pow(amps.data.data[n][0],2)+pow(amps.data.data[n][1],2)+pow(amps.data.data[n][2],2));




  //process the gas output file
  cPostProcess3D::cColumnIntegral::cColumnIntegrationSet Gas;

  Gas.IntegrantVector=GAS::IntegrantVector;
  Gas.IntegrantVectorLength=GAS::IntegrantVectorLength;
  Gas.PrintVariableList=GAS::PrintVariableList;

  amps.ColumnIntegral.Map.Circular(xObservation,xPrimary,xSecondary,1.5,100,100,"map.dat",&Gas);
  amps.SaveDataFile("Res.dat", amps.data);


  return 1;
}


