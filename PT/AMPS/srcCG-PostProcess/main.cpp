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

#include "config.dfn"
#include "config.h"

//header of the functions for calculating of the dust scattering efficientcy




//post-processing object
cPostProcess3D amps;

//location of the Sun in the frame of reference related to the surface model
double xSun[3];

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

//================================================
//filter trajectories that are not consistent with the observations
//and determine the set of the surface faces the dust injection from which will be taken into account
namespace TRAJECTORY_FILTER {
  bool ProcessTrajectoriesInitFlag=false;
  bool *SelectedFaceTable=NULL;
  bool *SelectedTrajectoriesTable=NULL;

  double GetFaceCorrectionFactor(int nface) {
    return (SelectedFaceTable[nface]==true) ? 1.0 : 0.0;
  }

  void ProcessTrajectories(SpiceDouble et) {
    int nTrajectory,n,i,iPixel,jPixel;
    double x[3],c0,c1,l,lCorridor[2];
    cVirtisM Virtis;

    //the trajectory acceptance corridor
    struct cCorridor {
      double x0,y0;
      double x1,y1;
    };

/*    double x0PixelCorridor=80;
    double y0PixelCorridor=170.0;

    double x1PixelCorridor=90;
    double y1PixelCorridor=160.0;*/


    cCorridor IncludeCorridor={80,170.0,90,160.0};
    cCorridor ExcludeCorridor={180,94,190,83};


    double dyPixel=20.0;

    //set orientation of the axis of VIRTIS
    Virtis.SetFrameAxis(et);


    ProcessTrajectoriesInitFlag=true;

    SelectedTrajectoriesTable=new bool [amps.ParticleTrajectory.nTotalTrajectories];
    for (i=0;i<amps.ParticleTrajectory.nTotalTrajectories;i++) SelectedTrajectoriesTable[i]=false;

    SelectedFaceTable=new bool [amps.SurfaceData.nCells];
    for (i=0;i<amps.SurfaceData.nCells;i++) SelectedFaceTable[i]=false;

    //add surface faces into the allowed list
    lCorridor[0]=IncludeCorridor.x1-IncludeCorridor.x0;
    lCorridor[1]=IncludeCorridor.y1-IncludeCorridor.y0;

    for (nTrajectory=0;nTrajectory<=amps.ParticleTrajectory.nTotalTrajectories;nTrajectory++) {
      for (n=0;n<amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].nDataPoints;n++) {
        for (i=0;i<3;i++) x[i]=amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[n][i];

        //determine projection of 'x' as would be seen by VIRTIS-M
        for (c0=0.0,c1=0.0,l=0.0,i=0;i<3;i++) {
          c0+=(x[i]-Virtis.xRosetta[i])*Virtis.e0[i];
          c1+=(x[i]-Virtis.xRosetta[i])*Virtis.e1[i];
          l+=pow(x[i]-Virtis.xRosetta[i],2);
        }

        l=sqrt(l);

        iPixel=(Pi/2.0-acos(c0/l))/Virtis.dAnglePixel+Virtis.nFieldOfViewPixels/2;
        jPixel=(Pi/2.0-acos(c1/l))/Virtis.dAnglePixel+Virtis.nFieldOfViewPixels/2;

        //Trajectory is accepted is it has a point that falls into the corridor
        double r=(iPixel-IncludeCorridor.x0)/lCorridor[0]*lCorridor[1]+IncludeCorridor.y0;

        if ((fabs(jPixel-r)<dyPixel) && (iPixel>=IncludeCorridor.x0)  && (iPixel<=IncludeCorridor.x1)) {
          //trajectory falls within the corridor
          //add the face into the list of the faces from which the dust ejection is allowed
          SelectedTrajectoriesTable[nTrajectory]=true;
          SelectedFaceTable[(int)amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][7]]=true;
        }
      }
    }

    //exclude faces that produses jets seen on the reverse side
    lCorridor[0]=ExcludeCorridor.x1-ExcludeCorridor.x0;
    lCorridor[1]=ExcludeCorridor.y1-ExcludeCorridor.y0;

    for (nTrajectory=0;nTrajectory<=amps.ParticleTrajectory.nTotalTrajectories;nTrajectory++) {
      for (n=0;n<amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].nDataPoints;n++) {
        for (i=0;i<3;i++) x[i]=amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[n][i];

        //determine projection of 'x' as would be seen by VIRTIS-M
        for (c0=0.0,c1=0.0,l=0.0,i=0;i<3;i++) {
          c0+=(x[i]-Virtis.xRosetta[i])*Virtis.e0[i];
          c1+=(x[i]-Virtis.xRosetta[i])*Virtis.e1[i];
          l+=pow(x[i]-Virtis.xRosetta[i],2);
        }

        l=sqrt(l);

        iPixel=(Pi/2.0-acos(c0/l))/Virtis.dAnglePixel+Virtis.nFieldOfViewPixels/2;
        jPixel=(Pi/2.0-acos(c1/l))/Virtis.dAnglePixel+Virtis.nFieldOfViewPixels/2;

        //Trajectory is accepted is it has a point that falls into the corridor
        double r=(iPixel-ExcludeCorridor.x0)/lCorridor[0]*lCorridor[1]+ExcludeCorridor.y0;

        if ((fabs(jPixel-r)<dyPixel) && (iPixel>=ExcludeCorridor.x0)  && (iPixel<=ExcludeCorridor.x1)) {
          //trajectory falls within the corridor
          //add the face into the list of the faces from which the dust ejection is allowed
          SelectedTrajectoriesTable[nTrajectory]=false;
          SelectedFaceTable[(int)amps.ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][7]]=false;
        }
      }
    }

  }

  //output selected trajectories
  void PrintSelectedTrajectories(int nPrintedTrajectories,const char *fname) {
    int i,nt,cnt;
    vector<int> PrintedTrajectoriesTable;

    //output the header of the trajectory file
    amps.ParticleTrajectory.PrintDataFileHeader(fname);

    //output trajectories
    for (cnt=0;cnt<nPrintedTrajectories;cnt++) {
      bool foundflag=false;

      do {
        foundflag=false;

        do {
          nt=(int)(rnd()*amps.ParticleTrajectory.nTotalTrajectories);
        }
        while (SelectedTrajectoriesTable[nt]==false);

        for (i=0;i<PrintedTrajectoriesTable.size();i++) if (PrintedTrajectoriesTable[i]==nt) {
          foundflag=true;
          break;
        }
      }
      while (foundflag==true);

      PrintedTrajectoriesTable.push_back(nt);
      amps.ParticleTrajectory.AddTrajectoryDataFile(&amps.ParticleTrajectory.IndividualTrajectories[nt],nt,fname);
    }
  }

}


//=================================================
//print the surface additional data
namespace SURFACE {

   namespace ILLUMINATION {
     //calcualte the surface exposure time
     double etSunLocation=0.0,xSunLocation[3];
     double *ExposureTime=NULL;
     double *cosSolarZenithAngle=NULL;
     int *IlluminationMap=NULL;

     void InitDataBuffer() {
       ExposureTime=new double [CutCell::nBoundaryTriangleFaces];
       cosSolarZenithAngle=new double [CutCell::nBoundaryTriangleFaces];
       IlluminationMap=new int [CutCell::nBoundaryTriangleFaces];

       for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) ExposureTime[i]=0.0,cosSolarZenithAngle[i]=0.0,IlluminationMap[i]=false;
     }

     void GetIlliminationMapElement(int nface,double& cosSZA,int& IlluminationFlag,double* xLight) {
       double xCenter[3],*Normal,r2,l[3],t;
       int i;

       CutCell::BoundaryTriangleFaces[nface].GetCenterPosition(xCenter);
       Normal=CutCell::BoundaryTriangleFaces[nface].ExternalNormal;
       IlluminationFlag=true;

       //get the solar zenith angle
       for (i=0,r2=0.0,cosSZA=0.0;i<3;i++) {
         l[i]=xLight[i]-xCenter[i];

         r2+=pow(l[i],2),cosSZA+=l[i]*Normal[i];
       }

       cosSZA/=sqrt(r2);

       //determine whether the face is in a shadow
       if (cosSZA<=0.0) IlluminationFlag=false;
       else for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) if (i!=nface) {
         if (CutCell::BoundaryTriangleFaces[i].RayIntersection(xCenter,l,t,0.0)==true) if (t>0.0) {
           IlluminationFlag=false;
           break;
         }
       }
     }

     void CalculateIlliminationMap(SpiceDouble et) {
       int nface;

       etSunLocation=et;
       if (IlluminationMap==NULL) InitDataBuffer();

       //get the location of the Sun
       SpiceDouble StateSun[6],lt;
       spkezr_c("SUN",et,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateSun,&lt);
       for (int i=0;i<3;i++) xSunLocation[i]=1.0E3*StateSun[i];

       //init the solar zenith angle and illumination map tables
       for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) if (nface%amps.size==amps.rank) {
         double cosSZA;
         int IlluminationFlag;

         GetIlliminationMapElement(nface,cosSZA,IlluminationFlag,xSunLocation);
         cosSolarZenithAngle[nface]=cosSZA;
         IlluminationMap[nface]=IlluminationFlag;
       }

       //exchange the table between all processors
       if (amps.rank==0) {
         int thread;
         double tmpCosSolarZenithAngle[CutCell::nBoundaryTriangleFaces];
         int tmpIlluminationMap[CutCell::nBoundaryTriangleFaces];
         MPI_Status status;

         for (thread=1;thread<amps.size;thread++) {
           MPI_Recv(tmpCosSolarZenithAngle,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,thread,0,MPI_COMM_WORLD,&status);
           MPI_Recv(tmpIlluminationMap,CutCell::nBoundaryTriangleFaces,MPI_INT,thread,0,MPI_COMM_WORLD,&status);

           for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) if (nface%amps.size==thread) {
             IlluminationMap[nface]=tmpIlluminationMap[nface];
             cosSolarZenithAngle[nface]=tmpCosSolarZenithAngle[nface];
           }
         }
       }
       else {
         MPI_Send(cosSolarZenithAngle,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
         MPI_Send(IlluminationMap,CutCell::nBoundaryTriangleFaces,MPI_INT,0,0,MPI_COMM_WORLD);
       }

       MPI_Bcast(cosSolarZenithAngle,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,0,MPI_COMM_WORLD);
       MPI_Bcast(IlluminationMap,CutCell::nBoundaryTriangleFaces,MPI_INT,0,MPI_COMM_WORLD);
     }

     void SetIlliminationMap(SpiceDouble et) {
       FILE *fInMap=NULL;

       //init the data buffer if needed
       if (IlluminationMap==NULL) InitDataBuffer();

       //check whether the binary file with the map exists
       fInMap=fopen("illumination-map.bin","r");

       //if the file exists -> check the 'et' and the number of the cutsells in the surface model
       if (fInMap!=NULL) {
         SpiceDouble etSaved;
         int nBoundaryTriangleFacesSaved;

         fread(&nBoundaryTriangleFacesSaved,sizeof(int),1,fInMap);
         fread(&etSaved,sizeof(SpiceDouble),1,fInMap);

         if ((etSaved!=et)||(nBoundaryTriangleFacesSaved!=CutCell::nBoundaryTriangleFaces)) {
           //the saved file does not corresponds to the current calculations -> re-calculate the illumination map
           fclose(fInMap);
           fInMap=NULL;
         }
       }

       if (fInMap==NULL) {
         //there is no file that containes the useful data: calculate the map, and save it
         CalculateIlliminationMap(et);

         if (amps.rank==0) {
           FILE *fOutMap=fopen("illumination-map.bin","w");

           fwrite(&CutCell::nBoundaryTriangleFaces,sizeof(int),1,fOutMap);
           fwrite(&et,sizeof(SpiceDouble),1,fOutMap);

           fwrite(cosSolarZenithAngle,sizeof(double),CutCell::nBoundaryTriangleFaces,fOutMap);
           fwrite(IlluminationMap,sizeof(int),CutCell::nBoundaryTriangleFaces,fOutMap);

           fclose(fOutMap);
         }
       }
       else {
         //the file exists -> read it
         fread(cosSolarZenithAngle,sizeof(double),CutCell::nBoundaryTriangleFaces,fInMap);
         fread(IlluminationMap,sizeof(int),CutCell::nBoundaryTriangleFaces,fInMap);

         fclose(fInMap);
       }
     }

     //calculate the surface face exposure time
     void CalculateSurfaceExposureTime(SpiceDouble et,double SearchTimeInterval,double SearchTimeIncrement) {
       int nface;
       double AccumulatedSearchTime=0.0;

       //init the data buffer if needed
       if (IlluminationMap==NULL) InitDataBuffer();

       //init the data buffer and the current illumination state
       int localFaceIlluminationMap[CutCell::nBoundaryTriangleFaces],tmpBuffer[CutCell::nBoundaryTriangleFaces];
       double SunState[6],xLight[3],lt;

       SetIlliminationMap(et);
       memcpy(localFaceIlluminationMap,IlluminationMap,CutCell::nBoundaryTriangleFaces*sizeof(int));

       //setup the exposure time table
       for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) ExposureTime[nface]=(localFaceIlluminationMap[nface]==true) ? SearchTimeInterval : -1.0;

       //determine the exposure time
       int nFaceChangedState,FaceChangedStateList[CutCell::nBoundaryTriangleFaces];
       double ExposureTimeExchengeBuffer[CutCell::nBoundaryTriangleFaces];

       while (AccumulatedSearchTime<SearchTimeInterval) {
         AccumulatedSearchTime+=SearchTimeIncrement;
         nFaceChangedState=0;

         //get the location of the Sun
         SpiceDouble StateSun[6],lt;
         spkezr_c("SUN",et-AccumulatedSearchTime,"67P/C-G_CK","none","CHURYUMOV-GERASIMENKO",StateSun,&lt);
         for (int i=0;i<3;i++) xLight[i]=1.0E3*StateSun[i];

         //determine the new faces that was not in a shadow in time 'et-AccumulatedSearchTime'
         int cnt=0;

         for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) if (localFaceIlluminationMap[nface]==true) {
           if (cnt%amps.size==amps.rank) {
             double cosSZA;
             int IlluminationFlag;

             GetIlliminationMapElement(nface,cosSZA,IlluminationFlag,xLight);

             if (IlluminationFlag==false) {
               //the face was in shadow at time 'et-AccumulatedSearchTime'
               FaceChangedStateList[nFaceChangedState++]=nface;
               ExposureTime[nface]=AccumulatedSearchTime;
               localFaceIlluminationMap[nface]=false;
             }
           }

           cnt++;
         }

         //collect the time information from all processors
         int nFaceChangedStateGlobalList[amps.size];

         MPI_Gather(&nFaceChangedState,1,MPI_INT,nFaceChangedStateGlobalList,1, MPI_INT,0,MPI_COMM_WORLD);

         if (amps.rank==0) {
            int n,thread,buffer[CutCell::nBoundaryTriangleFaces];
            MPI_Status status;

            //get the list of new faces coming into shadow
            for (thread=1;thread<amps.size;thread++) {
              MPI_Recv(buffer,nFaceChangedStateGlobalList[thread],MPI_INT,thread,0,MPI_COMM_WORLD,&status);

              for (n=0;n<nFaceChangedStateGlobalList[thread];n++) {
                localFaceIlluminationMap[buffer[n]]=false;
                ExposureTime[buffer[n]]=AccumulatedSearchTime;
              }
            }
         }
         else {
           //the slave processor;
           MPI_Send(FaceChangedStateList,nFaceChangedState,MPI_INT,0,0,MPI_COMM_WORLD);
         }

         MPI_Bcast(localFaceIlluminationMap,CutCell::nBoundaryTriangleFaces,MPI_INT,0,MPI_COMM_WORLD);
       }

       //distribute the exposure time array
       MPI_Bcast(ExposureTime,CutCell::nBoundaryTriangleFaces,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }

     void SetSurfaceExposureTime(SpiceDouble et,double SearchTimeInterval,double SearchTimeIncrement) {
       FILE *fInMap=NULL;

       //init the data buffer if needed
       if (IlluminationMap==NULL) InitDataBuffer();

       //check whether the binary file with the map exists
       fInMap=fopen("exposure-time-map.bin","r");

       //if the file exists -> check the 'et' and the number of the cutsells in the surface model
       if (fInMap!=NULL) {
         SpiceDouble etSaved;
         int nBoundaryTriangleFacesSaved;
         double SearchTimeIntervalSaved,SearchTimeIncrementSaved;

         fread(&nBoundaryTriangleFacesSaved,sizeof(int),1,fInMap);
         fread(&etSaved,sizeof(SpiceDouble),1,fInMap);
         fread(&SearchTimeIntervalSaved,sizeof(double),1,fInMap);
         fread(&SearchTimeIncrementSaved,sizeof(double),1,fInMap);

         if ((etSaved!=et)||(nBoundaryTriangleFacesSaved!=CutCell::nBoundaryTriangleFaces)||(SearchTimeIntervalSaved!=SearchTimeInterval)||(SearchTimeIncrementSaved!=SearchTimeIncrement)) {
           //the saved file does not corresponds to the current calculations -> re-calculate the illumination map
           fclose(fInMap);
           fInMap=NULL;
         }
       }

       if (fInMap==NULL) {
         //there is no file that containes the useful data: calculate the map, and save it
         CalculateSurfaceExposureTime(et,SearchTimeInterval,SearchTimeIncrement);

         if (amps.rank==0) {
           FILE *fOutMap=fopen("exposure-time-map.bin","w");

           fwrite(&CutCell::nBoundaryTriangleFaces,sizeof(int),1,fOutMap);
           fwrite(&et,sizeof(SpiceDouble),1,fOutMap);
           fwrite(&SearchTimeInterval,sizeof(double),1,fOutMap);
           fwrite(&SearchTimeIncrement,sizeof(double),1,fOutMap);

           fwrite(ExposureTime,sizeof(double),CutCell::nBoundaryTriangleFaces,fOutMap);
           fclose(fOutMap);
         }
       }
       else {
         //the file exists -> read it
         fread(ExposureTime,sizeof(double),CutCell::nBoundaryTriangleFaces,fInMap);
         fclose(fInMap);
       }
     }
   }





  //print the surface data
  void PrintVariableList(FILE *fout) {fprintf(fout," \"Shadow Flag\", \"cos(Solar Zenith Angle)\","
      " \"Sun Exposure Time\", \"cosSZA/Exposure Time\", \"H2O Source Rate/Exposure Time\", \"H2O Source Rate\""); }
  int GetVariableNumber() {return 6;}

  void GetFaceDataVector(double *DataVector,CutCell::cTriangleFace *face,int nface) {
    if (ILLUMINATION::IlluminationMap==NULL) exit(__LINE__,__FILE__,"Error: the illumination map need to be initalized first");

    DataVector[0]=(ILLUMINATION::IlluminationMap[nface]==true) ? 1 : 0;
    DataVector[1]=ILLUMINATION::cosSolarZenithAngle[nface];
    DataVector[2]=ILLUMINATION::ExposureTime[nface];

    double FaceSourceRate_H2O=0.0;
    int i,nnode;

    for (i=0;i<3;i++) FaceSourceRate_H2O+=0.3*amps.SurfaceData.data[amps.SurfaceData.ConnectivityList[nface][i]][3];

    DataVector[3]=(ILLUMINATION::IlluminationMap[nface]==true) ? ILLUMINATION::cosSolarZenithAngle[nface]/ILLUMINATION::ExposureTime[nface] : -1;
    DataVector[4]=(ILLUMINATION::IlluminationMap[nface]==true) ? FaceSourceRate_H2O/ILLUMINATION::ExposureTime[nface] : -1;

    DataVector[5]=FaceSourceRate_H2O;
  }
}

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

      data[1]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][38]/2.99151E-26;
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
#elif _DUST_CASE_ == _DUST_CASE__4SPEC_1GROUP_
  int nRadii=1;
  int MeanDensityOffset=46;

  cVariablePair VariablePair[]={
      {51,0.5*(1.000000E-07+1.000000E-04)}
  };
#elif _DUST_CASE_ == _DUST_CASE__1SPEC_10GROUP_
  int nRadii=10;

  int MeanDensityOffset=43;

  cVariablePair VariablePair[]={
        {48,0.5*(1.000000E-07+1.995262E-07)},
        {53,0.5*(1.995262E-07+3.981072E-07)},
        {58,0.5*(3.981072E-07+7.943282E-07)},
        {63,0.5*(7.943282E-07+1.584893E-06)},
        {68,0.5*(1.584893E-06+3.162278E-06)},
        {73,0.5*(3.162278E-06+6.309573E-06)},
        {78,0.5*(6.309573E-06+1.258925E-05)},
        {83,0.5*(1.258925E-05+2.511886E-05)},
        {88,0.5*(2.511886E-05+5.011872E-05)},
        {93,0.5*(5.011872E-05+1.000000E-04)}
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



  int IntegrantVectorLength() {return 7+3*nRadii;}

  void PrintVariableList(FILE* fout) {
    fprintf(fout," \"Column Densty\", \"Total Brightness\", \"Number of Trajectory points\", \"Number of Independent Trajectories\", \"Mean Velocity\", \"Corrected Column Densty\", \"Cortrected Total Brightness\",");

    for (int i=0;i<nRadii;i++) fprintf(fout,", \"Column Density(%e)\", \"Brightness(%e)\", \"Velocity(%e)\"",VariablePair[i].GrainRadius,VariablePair[i].GrainRadius,VariablePair[i].GrainRadius);
  }

  void PostProcessColumnIntegralVector(double* data) {

    //total speed
    if (data[0]>0.0) data[4]/=data[0];

    //speed for individual dust groupls
    for (int iRadius=0;iRadius<nRadii;iRadius++) if (data[7+3*iRadius]>0.0) data[7+3*iRadius+2]/=data[7+3*iRadius];
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


      if ((data[0]<0.0)||(data[4]<0.0)) {
        double t=234;
        t+=30;

        std::cout << amps.data.data[Stencil.Node[i]][0] << "  " << amps.data.data[Stencil.Node[i]][1] << "  " << amps.data.data[Stencil.Node[i]][2] << endl; 
        std::cout << Stencil.Weight[i] << "  " << Stencil.Node[i] << "   " << MeanDensityOffset << "  " << amps.data.data[Stencil.Node[i]][MeanDensityOffset] << "  " << data[0] << endl;

        printf("Error: a value that must be positive is negative\n");
        exit(0);
      }
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
        data[7+3*iRadius]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
        data[7+3*iRadius+2]+=Stencil.Weight[i]*amps.data.data[Stencil.Node[i]][4+VariablePair[iRadius].nVar]*amps.data.data[Stencil.Node[i]][VariablePair[iRadius].nVar];
      }

      if (data[7+3*iRadius]<0.0) {
        double t=3435;
        t+=4053;

        printf("Error: a value that must be positive is negative\n");
        exit(0);
       } 

      data[7+3*iRadius+1]=ScatteringEfficentcy*data[7+3*iRadius];
      TotalBrightness+=data[7+3*iRadius+1];
    }

    data[1]=TotalBrightness;

    //calculate the total column density (data[5]) and brightness (data[6]) corrected with the surface exposure time
    int np,nt,nStatingFace;
    double CorrectionFactor;

    for (np=0;np<cl->TrajectoryPoints.size();np++) {
      nt=cl->TrajectoryPoints[np];
      nStatingFace=amps.ParticleTrajectory.IndividualTrajectories[nt].Data[0][7];

      double FaceSourceRate_H2O=0.0,c;
      int i,nnode;

      for (i=0;i<3;i++) FaceSourceRate_H2O+=0.3*amps.SurfaceData.data[amps.SurfaceData.ConnectivityList[nStatingFace][i]][3];

      CorrectionFactor=SURFACE::ILLUMINATION::cosSolarZenithAngle[nStatingFace]/SURFACE::ILLUMINATION::ExposureTime[nStatingFace];
      CorrectionFactor=FaceSourceRate_H2O/SURFACE::ILLUMINATION::ExposureTime[nStatingFace];

      CorrectionFactor*=TRAJECTORY_FILTER::GetFaceCorrectionFactor(nStatingFace);
//      if (TRAJECTORY_FILTER::SelectedTrajectoriesTable[nt]==false) CorrectionFactor=0.0;

      if (CorrectionFactor<0.0) CorrectionFactor=0.0;

      c=CorrectionFactor*amps.ParticleTrajectory.IndividualTrajectories[nt].Data[0][8]/
          (bl->xmax[0]-bl->xmin[0])/(bl->xmax[1]-bl->xmin[1])/(bl->xmax[2]-bl->xmin[2]);

      data[5]+=c;

      //particle size
      GrainRadius=amps.ParticleTrajectory.IndividualTrajectories[nt].Data[0][6];
      ScatteringEfficentcy=
              LK::GetScatteringEfficeintcy(GrainRadius,LK::Ice2Dust0_899999976__Porosity0_649122834::Data,LK::Ice2Dust0_899999976__Porosity0_649122834::nDataPoints);

      data[6]+=c*ScatteringEfficentcy;
    }
  }

}


int main(int argc,char **argv) {

  amps.InitMPI();
  amps.SetBlockSize(5,5,5);


/*
//    amps.ParticleTrajectory.LoadDataFile("amps.TrajectoryTracking.out=1.s=0.H2O.dat",".");

//amps.ParticleTrajectory.LoadDataFile("amps.TrajectoryTracking.out=1.s=2.DUST%0.dat",".");
amps.ParticleTrajectory.LoadDataFile("amps.TrajectoryTracking.out=5.s=5.DUST%3.dat",".");

    amps.PrintParticleTrajectory(1200,_OUTPUT_MODE__UNIFORM_,NULL,"limited-trajectories.gas.uniform.dat");
    amps.PrintParticleTrajectory(1200,_OUTPUT_MODE__FLUX_,AcceptParticleTrajectory,"limited-trajectories.gas.flux.dat"); 

MPI_Finalize();
return 1;
*/




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
    xSun[i]=1.0E3*StateSun[i];
  }


  std::cout << "Heliocentric Distance " << sqrt(StateSun[0]*StateSun[0]+StateSun[1]*StateSun[1]+StateSun[2]*StateSun[2])*1.0E3/_AU_ << std::endl;

  if (_MODE_==_GAS_MODE_) {
    std::string out,DataFileName;
    std::stringstream t;

    t << _OUTPUT_FILE_NUMBER_;
    out=t.str();

    DataFileName="pic.H2O.s=0.out="+out+".dat";

    amps.LoadDataFile(DataFileName.c_str(),".");
    amps.PrintVariableList();

    //load the surface data
    std::string SurfaceDataFileName="amps.cut-cell.surface-data.out="+out+".dat";
    amps.SurfaceData.LoadDataFile(SurfaceDataFileName.c_str(),".");
    amps.SurfaceData.PrintVariableList();


    //load particle trajectory file
     std::string TrajectoryFileName;
 
     t.str("");
     t << _OUTPUT_FILE_NUMBER_;
     out=t.str();
 
     TrajectoryFileName="amps.TrajectoryTracking.out="+out+".s=0.H2O.dat";
     amps.ParticleTrajectory.LoadDataFile(TrajectoryFileName.c_str(),".");

    amps.PrintParticleTrajectory(300,_OUTPUT_MODE__UNIFORM_,NULL,"limited-trajectories.gas.uniform.dat");
    amps.PrintParticleTrajectory(300,_OUTPUT_MODE__FLUX_,AcceptParticleTrajectory,"limited-trajectories.gas.flux.dat");
  }
  else if (_MODE_==_DUST_MODE_) {
    std::string TrajectoryFileName[_DUST_SPEC_NUMBER_];
    std::string s,out;
    std::stringstream t;

    t << _OUTPUT_FILE_NUMBER_;
    out=t.str();

    for (int i=0;i<_DUST_SPEC_NUMBER_;i++) {
      std::string s1;

      t.str("");
      t << _GAS_SPEC_NUMBER_+i;
      s=t.str();

      t.str("");
      t << i;
      s1=t.str();

      TrajectoryFileName[i]="amps.TrajectoryTracking.out="+out+".s="+s+".DUST%"+s1+".dat";
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

    for (int i=0;i<_DUST_SPEC_NUMBER_;i++) amps.ParticleTrajectory.LoadDataFile(TrajectoryFileName[i].c_str(),".");
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
  amps.ParticleTrajectory.PrintSurfaceData("surface.dat",NULL,NULL,NULL);


//  for (int n=0;n<amps.data.nNodes;n++) amps.data.data[n][38]=sqrt(pow(amps.data.data[n][0],2)+pow(amps.data.data[n][1],2)+pow(amps.data.data[n][2],2));


  //output surface data
  SURFACE::ILLUMINATION::SetSurfaceExposureTime(et,6.0*3600.0,10.0*60);
  SURFACE::ILLUMINATION::SetIlliminationMap(et);
  amps.ParticleTrajectory.PrintSurfaceData("surface-data.dat",SURFACE::GetVariableNumber,SURFACE::PrintVariableList,SURFACE::GetFaceDataVector);

//  MPI_Finalize();
//  return 1;

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


  //process all trajectories: select those that falls into the jet corridor
  TRAJECTORY_FILTER::ProcessTrajectories(et);
  TRAJECTORY_FILTER::PrintSelectedTrajectories(500,"selected-trajectories.dat");

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


