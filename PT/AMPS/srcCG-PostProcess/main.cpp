//$Id$

/*
 * main.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include "PostProcess3D.h"
#include "SpiceUsr.h"
#include "DustScatteringEfficientcy.h"


#define _DUST_MODE_  0
#define _GAS_MODE_   1

//#define _MODE_ _GAS_MODE_
#define _MODE_  _DUST_MODE_


//header of the functions for calculating of the dust scattering efficientcy




//post-processing object
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

  int IntegrantVectorLength() {return 2+2*nRadii;}

  void PrintVariableList(FILE* fout) {
    fprintf(fout," \"Column Densty\", \"Total Brightness\"");

    for (int i=0;i<nRadii;i++) fprintf(fout,", \"Column Density(%e)\", \"Brightness(%e)\"",VariablePair[i].GrainRadius,VariablePair[i].GrainRadius);
  }

  void IntegrantVector(double *data,double *x) {
    cPostProcess3D::cStencil Stencil;
    int i,iRadius;
    double GrainRadius,ScatteringEfficentcy,TotalBrightness=0.0;

    amps.GetInterpolationStencil(x,&Stencil);

    for (i=0;i<IntegrantVectorLength();i++) data[i]=0.0;

    for (i=0;i<8;i++) data[0]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][46];

    for (iRadius=0;iRadius<nRadii;iRadius++) {
      GrainRadius=VariablePair[iRadius].GrainRadius;

      ScatteringEfficentcy=
        LK::GetScatteringEfficeintcy(GrainRadius,LK::Ice2Dust0_899999976__Porosity0_649122834::Data,LK::Ice2Dust0_899999976__Porosity0_649122834::nDataPoints);

      for (i=0;i<8;i++) {
        data[2+2*iRadius]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
      }

      data[2+2*iRadius+1]=ScatteringEfficentcy*data[2+2*iRadius];
      TotalBrightness+=data[2+2*iRadius+1];
    }

    data[1]=TotalBrightness;
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
    amps.PrintVariableList();
  }
  else if (_MODE_==_DUST_MODE_) {
    amps.LoadDataFile("pic.DUST%0.s=2.out=3.dat",".");
    amps.PrintVariableList();

    //load the trajectory data file
    const int nTrajectoryFiles=4;
    char TrajectoryFileName[nTrajectoryFiles][100]={"amps.TrajectoryTracking.out=3.s=2.DUST%0.dat","amps.TrajectoryTracking.out=3.s=3.DUST%1.dat",
        "amps.TrajectoryTracking.out=3.s=4.DUST%2.dat","amps.TrajectoryTracking.out=3.s=5.DUST%3.dat"
    };

    for (int i=0;i<nTrajectoryFiles;i++) amps.ParticleTrajectory.LoadDataFile(TrajectoryFileName[i],".");
    amps.ParticleTrajectory.PrintVariableList();

    //load the surface data
    amps.SurfaceData.LoadDataFile("amps.cut-cell.surface-data.out=3.dat",".");
    amps.SurfaceData.PrintVariableList();

    //get the list of the faces that production rate above 85%
    double MinSourceRate=-1.0,MaxSouceRate=-1.0,LimitSourceRate;
    vector<int> FaceList;
    int n,ncell;

    for (n=0;n<amps.SurfaceData.nNodes;n++) {
      if ((MinSourceRate<0.0)||(MinSourceRate>amps.SurfaceData.data[n][3])) MinSourceRate=amps.SurfaceData.data[n][3];
      if ((MaxSouceRate<0.0)||(MaxSouceRate<amps.SurfaceData.data[n][3])) MaxSouceRate=amps.SurfaceData.data[n][3];
    }

    LimitSourceRate=MinSourceRate+0.85*(MaxSouceRate-MinSourceRate);

    for (n=0;n<amps.SurfaceData.nNodes;n++) if (amps.SurfaceData.data[n][3]>LimitSourceRate) {
      for (int i=0;i<amps.SurfaceData.NodeBall[n].size();i++) {
        bool found=false;

        ncell=amps.SurfaceData.NodeBall[n][i];
        for (int ii=0;ii<FaceList.size();ii++) if (FaceList[ii]==ncell) {
          found=true;
          break;
        }

        if (found==false) FaceList.push_back(ncell);
      }
    }

    printf("Face List begin:\n");
    for (int i=0;i<FaceList.size();i++) printf("%i ",FaceList[i]);
    printf("Face List end:\n");

    //output the set of the trajectories emerged from the maximum source rate faces
    int cnt=0;
    amps.ParticleTrajectory.PrintDataFileHeader("maxSourceRateTrajectories.dat");

    for (int n=0;n<amps.ParticleTrajectory.nTotalTrajectories;n++) {
      int nface=amps.ParticleTrajectory.IndividualTrajectories[n].Data[7+0*amps.ParticleTrajectory.nTrajectoryVariables];

      for (int iface=0;iface<FaceList.size();iface++) if (FaceList[iface]==nface) {
        amps.ParticleTrajectory.AddTrajectoryDataFile(&amps.ParticleTrajectory.IndividualTrajectories[n],cnt,"maxSourceRateTrajectories.dat");
        cnt++;

        break;
      }

    }





  }




  //load the nucleus mesh
  CutCell::ReadNastranSurfaceMeshLongFormat("SHAP5_stefano.bdf","/Volumes/data/AMPS_DATA/ROSETTA");
  amps.ParticleTrajectory.PrintSurfaceData("surface.dat");


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

  amps.ColumnIntegral.Map.Circular(xObservation,xPrimary,xSecondary,2,200,200,"map.dat",&IntegrationSet);
  amps.SaveDataFile("Res.dat", amps.data);


  amps.FinalizeMPI();
  return 1;
}


