//$Id$

/*
 * main.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>

#include "PostProcess3D.h"
#include "VIRTIS-M.h"
#include "SpiceUsr.h"
#include "DustScatteringEfficientcy.h"


#define _DUST_MODE_  0
#define _GAS_MODE_   1

#define _DUST_CASE__1SPEC_1GROUP_   0
#define _DUST_CASE__4SPEC_1GROUP_   1
#define _DUST_CASE__1SPEC_10GROUP_  2
#define _DUST_CASE__4SPEC_10GROUP_  4

//#define _MODE_ _GAS_MODE_
#define _MODE_  _DUST_MODE_
#define _DUST_CASE_ _DUST_CASE__1SPEC_1GROUP_

std::string OutputNumber="6";

//header of the functions for calculating of the dust scattering efficientcy




//post-processing object
cPostProcess3D amps;

//acceptance probability of a particle trajectory when printed into a file
double AcceptParticleTrajectory(int nTrajectory) { 
  return amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][8]/amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][4]; 
}


//process the dust file
struct cColumnIntegrationSet {
  int (*IntegrantVectorLength)();
  void (*IntegrantVector)(double* Data,double *Location);
  void (*PrintVariableList)(FILE* fout);
  void (*PostProcessColumnIntegralVector)(double*);
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

  void PostProcessColumnIntegralVector(double* data) {}
}

namespace DUST {
  struct cVariablePair {
    int nVar;
    double GrainRadius;
  };

#if _DUST_CASE_ == _DUST_CASE__1SPEC_1GROUP_
  int nRadii=1;
  int MeanDensityOffset=43;

  cVariablePair VariablePair[]={
      {48,0.5*(1.000000E-07+1.000000E-04)}
  };
#elif _DUST_CASE_ == _DUST_CASE__4SPEC_10GROUP_
  int nRadii=10;

  int MeanDensityOffset=46;

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
#else
#error "Dont know what is that"
#endif



  int IntegrantVectorLength() {return 5+3*nRadii;}

  void PrintVariableList(FILE* fout) {
    fprintf(fout," \"Column Densty\", \"Total Brightness\", \"Number of Trajectory points\", \"Number of Independent Trajectories\", \"Mean Velocity\",");

    for (int i=0;i<nRadii;i++) fprintf(fout,", \"Column Density(%e)\", \"Brightness(%e)\", \"Velocity(%e)\"",VariablePair[i].GrainRadius,VariablePair[i].GrainRadius,VariablePair[i].GrainRadius);
  }

  void PostProcessColumnIntegralVector(double* data) {

    //total speed
    if (data[0]>0.0) data[4]/=data[0];

    //speed for individual dust groupls
    for (int iRadius=0;iRadius<nRadii;iRadius++) if (data[5+3*iRadius]>0.0) data[5+3*iRadius+2]/=data[5+3*iRadius];
  }

  void IntegrantVector(double *data,double *x) {
    cPostProcess3D::cStencil Stencil;
    int i,iRadius;
    double GrainRadius,ScatteringEfficentcy,TotalBrightness=0.0;

    amps.GetInterpolationStencil(x,&Stencil);

    for (i=0;i<IntegrantVectorLength();i++) data[i]=0.0;

    for (i=0;i<8;i++) {
      data[0]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][MeanDensityOffset];
      data[4]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][MeanDensityOffset]*amps.data.data[Stencil.Node[i]][4+MeanDensityOffset];
    }

    //get the number of the trajectories
    cPostProcess3D::cCell*  cl=amps.GetCell(x);
    cPostProcess3D::cBlock* bl=amps.GetBlock(x);
    data[2]=cl->TrajectoryPoints.size()/(bl->xmax[0]-bl->xmin[0])/(bl->xmax[1]-bl->xmin[1])/(bl->xmax[2]-bl->xmin[2]);
    data[3]=cl->IndividualTrajectories.size()/(bl->xmax[0]-bl->xmin[0])/(bl->xmax[1]-bl->xmin[1])/(bl->xmax[2]-bl->xmin[2]);

    for (iRadius=0;iRadius<nRadii;iRadius++) {
      GrainRadius=VariablePair[iRadius].GrainRadius;

      ScatteringEfficentcy=
        LK::GetScatteringEfficeintcy(GrainRadius,LK::Ice2Dust0_899999976__Porosity0_649122834::Data,LK::Ice2Dust0_899999976__Porosity0_649122834::nDataPoints);

      for (i=0;i<8;i++) {
        data[5+3*iRadius]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
        data[5+3*iRadius+2]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][4+VariablePair[iRadius].nVar]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
      }

      data[5+3*iRadius+1]=ScatteringEfficentcy*data[5+3*iRadius];
      TotalBrightness+=data[5+3*iRadius+1];
    }

    data[1]=TotalBrightness;
  }

}


int main(int argc,char **argv) {


  amps.InitMPI();
  amps.SetBlockSize(5,5,5);

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


  std::cout << "Heliocentric Distance " << sqrt(StateSun[0]*StateSun[0]+StateSun[1]*StateSun[1]+StateSun[2]*StateSun[2])*1.0E3/_AU_ << std::endl;



  if (_MODE_==_GAS_MODE_) {
    amps.LoadDataFile("pic.H2O.s=0.out=10.dat",".");
    amps.PrintVariableList();
  }
  else if (_MODE_==_DUST_MODE_) {

    const int nTrajectoryFiles=1;
    std::string out="6";

    const int nTrajectoryFileMax=10;
    std::string TrajectoryFileName[nTrajectoryFileMax];

    for (int i=0;i<nTrajectoryFileMax;i++) {
      std::string s;
      std::stringstream t;

      t << 2+i;
      s=t.str();
      TrajectoryFileName[i]="amps.TrajectoryTracking.out="+out+".s="+s+".DUST%0.dat";
    }

    //load the trajectory data file
    std::string fDataFile="pic.DUST%0.s=2.out="+out+".dat";
    amps.LoadDataFile(fDataFile.c_str(),".");
    amps.PrintVariableList();

/*    char TrajectoryFileName[nTrajectoryFiles][100]={"amps.TrajectoryTracking.out=2.s=2.DUST%0.dat","amps.TrajectoryTracking.out=2.s=3.DUST%1.dat",
        "amps.TrajectoryTracking.out=2.s=4.DUST%2.dat","amps.TrajectoryTracking.out=2.s=5.DUST%3.dat"
    };*/

/*    const int nTrajectoryFiles=1;
    char TrajectoryFileName[nTrajectoryFiles][100]={"amps.TrajectoryTracking.out=3.s=0.H2O.dat"
    };*/

    for (int i=0;i<nTrajectoryFiles;i++) amps.ParticleTrajectory.LoadDataFile(TrajectoryFileName[i].c_str(),".");
    amps.ParticleTrajectory.PrintVariableList();
    amps.AssignParticleTrajectoriesToCells();

    //load the surface data
    std::string SurfaceDataFileName="amps.cut-cell.surface-data.out="+out+".dat";
    amps.SurfaceData.LoadDataFile(SurfaceDataFileName.c_str(),".");
    amps.SurfaceData.PrintVariableList();

    //get the list of the faces that production rate above 65%
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
      int nface=amps.ParticleTrajectory.IndividualTrajectories[n].Data[0][7];

      for (int iface=0;iface<FaceList.size();iface++) if (FaceList[iface]==nface) {
        amps.ParticleTrajectory.AddTrajectoryDataFile(&amps.ParticleTrajectory.IndividualTrajectories[n],cnt,"maxSourceRateTrajectories.dat");
        cnt++;

        break;
      }

    }


    //==================================================================
    //====== CORRECT VELOCITY OF THE DUST TRAJECTORIES: BEGING
/*    if (_MODE_==_DUST_MODE_) {
      int i,j;

      for (i=0;i<amps.ParticleTrajectory.nTotalTrajectories;i++) for (j=0;j<amps.ParticleTrajectory.IndividualTrajectories[i].nDataPoints;j++) {
        if (amps.ParticleTrajectory.IndividualTrajectories[i].Data[j][4]>1.0) amps.ParticleTrajectory.IndividualTrajectories[i].Data[j][4]=1.0;
      }
    }*/

/*    for (int i=0;i<amps.ParticleTrajectory.nTotalTrajectories;i++) {
      double a=2.0;
      int n0,n1,n=amps.ParticleTrajectory.IndividualTrajectories[i].nDataPoints-1;
      double l=0.0,t,s,v0,v1;

      for (int j=1;j<amps.ParticleTrajectory.IndividualTrajectories[i].nDataPoints;j++) {
        n0=j-1;
        n1=j;

        for (int ii=0;ii<3;ii++) l+=pow(amps.ParticleTrajectory.IndividualTrajectories[i].Data[n1][ii]-amps.ParticleTrajectory.IndividualTrajectories[i].Data[n0][ii],2);

        l=sqrt(l);
        v0=amps.ParticleTrajectory.IndividualTrajectories[i].Data[n0][4];
        v1=amps.ParticleTrajectory.IndividualTrajectories[i].Data[n1][4];

        t=(v1-v0)/a;
        s=v0*t+a*t*t/2.0;

        if (fabs((l-s)/l)>1.0E-5) {
          std::cout << "Error: in the trajectory integration" << endl;
        }
      }
    }*/



    //====== CORRECT VELOCITY OF THE DUST TRAJECTORIES: END
    //==================================================================


    //print out the limitab trajectories number 
    amps.PrintParticleTrajectory(300,_OUTPUT_MODE__UNIFORM_,NULL,"limited-trajectories.uniform.dat");
    amps.PrintParticleTrajectory(300,_OUTPUT_MODE__FLUX_,AcceptParticleTrajectory,"limited-trajectories.flux.dat"); 
//    return 1;


    //output location of the spacecraft
    if (amps.rank==0) {
      FILE *fSpacecraft=fopen("spacecraft-location.dat","w");
      fprintf(fSpacecraft,"VARIABLES=\"X\", \"Y\", \"Z\"\n");
      fprintf(fSpacecraft,"%e %e %e\n",xObservation[0],xObservation[1],xObservation[2]);
      fclose(fSpacecraft);
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
    IntegrationSet.PostProcessColumnIntegralVector=GAS::PostProcessColumnIntegralVector;
  }
  else {
    IntegrationSet.IntegrantVector=DUST::IntegrantVector;
    IntegrationSet.IntegrantVectorLength=DUST::IntegrantVectorLength;
    IntegrationSet.PrintVariableList=DUST::PrintVariableList;
    IntegrationSet.PostProcessColumnIntegralVector=DUST::PostProcessColumnIntegralVector;
  }


  //calculate the column integrals as seen by VIRTIS-M
  cVirtisM VirtisM;

  VirtisM.BlockNucleus.SetPixelLimits(200,240);
  VirtisM.BlockNucleus.SetBlock(et,CutCell::nBoundaryTriangleFaces,CutCell::BoundaryTriangleFaces);

  VirtisM.BlockNucleus.GetColumnIntegralMap("Virtis.map.dat",&IntegrationSet,&amps);



  amps.ColumnIntegral.Map.Circular(xObservation,xPrimary,xSecondary,2,200,200,"map.dat",&IntegrationSet);
 // amps.SaveDataFile("Res.dat", amps.data);


  amps.FinalizeMPI();
  return 1;
}


