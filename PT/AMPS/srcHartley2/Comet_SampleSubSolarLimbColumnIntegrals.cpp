/*
 * Comet_SampleSubSolarLimbColumnIntegrals.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: vtenishe
 */

//functions that sample the column integrals at the subsolar point of the limb
//$Id$

#include "pic.h"

SpiceDouble Comet::Sampling::SubsolarLimbColumnIntegrals::etSampleBegin;
int Comet::Sampling::SubsolarLimbColumnIntegrals::SamplingPhase;
int Comet::Sampling::SubsolarLimbColumnIntegrals::firstPhaseRadialVelocityDirection;
int Comet::Sampling::SubsolarLimbColumnIntegrals::nOutputFile;

double Comet::Sampling::SubsolarLimbColumnIntegrals::rr=0.0;
long int Comet::Sampling::SubsolarLimbColumnIntegrals::nSampleAltitudeDistributionPoints=0;
//Comet::Sampling::SubsolarLimbColumnIntegrals::cSampleAltitudeDistrubutionBufferElement *Comet::Sampling::SubsolarLimbColumnIntegrals::SampleAltitudeDistrubutionBuffer=NULL;

Comet::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Comet::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_AntiSunwardMotion[nPhaseAngleIntervals];
Comet::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Comet::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_SunwardMotion[nPhaseAngleIntervals];
Comet::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Comet::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_AntiSunwardMotion__TotalModelRun[nPhaseAngleIntervals];
Comet::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Comet::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_SunwardMotion__TotalModelRun[nPhaseAngleIntervals];

void Comet::Sampling::SubsolarLimbColumnIntegrals::EmptyFunction() {}

//init the sampling block
void Comet::Sampling::SubsolarLimbColumnIntegrals::init() {
  int el;

  etSampleBegin=0.0;
  SamplingPhase=-1;
  firstPhaseRadialVelocityDirection=0;
  nOutputFile=0;

  //init the altitude distribution sampling buffer
  double R,dAlt;
  long int t;

  rr=(maxAltitude+dAltmax)/(maxAltitude+dAltmin);
  t=(long int)(log(dAltmax/dAltmin)/log(rr)-2.0);
  rr=pow(dAltmax/dAltmin,1.0/(t+2.0));

  R=dAltmin,dAlt=dAltmin;
  while (R<maxAltitude+1.0) {
    ++nSampleAltitudeDistributionPoints;

    R+=dAlt;
    dAlt*=rr;
  }


  //allocate the memory buffers
  int i,SamplingDataBufferOffset=0;
  int ColumnIntegralVectorLength=ColumnIntegral::GetVariableList(NULL);
  static double *SamplingDataBuffer=new double [4*nPhaseAngleIntervals*ColumnIntegralVectorLength];

  int AltitudeDistributionSamplingDataOffset=0;
  static double **AltitudeDistributionSamplingDataPointer=new double *[2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints];
  static double *AltitudeDistributionSamplingData=new double [2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints*ColumnIntegralVectorLength];


  for (i=0;i<4*nPhaseAngleIntervals*ColumnIntegralVectorLength;i++) SamplingDataBuffer[i]=0.0;
  for (i=0;i<2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints*ColumnIntegralVectorLength;i++) AltitudeDistributionSamplingData[i]=0.0;

  for (i=0;i<2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints;i++) AltitudeDistributionSamplingDataPointer[i]=AltitudeDistributionSamplingData+i*ColumnIntegralVectorLength;



  for (el=0;el<nPhaseAngleIntervals;el++) {
    SampleBuffer_AntiSunwardMotion[el].SampleColumnIntegrals=SamplingDataBuffer+SamplingDataBufferOffset;
    SamplingDataBufferOffset+=ColumnIntegralVectorLength;

    SampleBuffer_SunwardMotion[el].SampleColumnIntegrals=SamplingDataBuffer+SamplingDataBufferOffset;
    SamplingDataBufferOffset+=ColumnIntegralVectorLength;

    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].SampleColumnIntegrals=SamplingDataBuffer+SamplingDataBufferOffset;
    SamplingDataBufferOffset+=ColumnIntegralVectorLength;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals=AltitudeDistributionSamplingDataPointer+AltitudeDistributionSamplingDataOffset;
    AltitudeDistributionSamplingDataOffset+=nSampleAltitudeDistributionPoints;

    SampleBuffer_SunwardMotion__TotalModelRun[el].SampleColumnIntegrals=SamplingDataBuffer+SamplingDataBufferOffset;
    SamplingDataBufferOffset+=ColumnIntegralVectorLength;
    SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals=AltitudeDistributionSamplingDataPointer+AltitudeDistributionSamplingDataOffset;
    AltitudeDistributionSamplingDataOffset+=nSampleAltitudeDistributionPoints;
  }


/*  //OLDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

  if (SampleAltitudeDistrubutionBuffer==NULL) SampleAltitudeDistrubutionBuffer=new cSampleAltitudeDistrubutionBufferElement[2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints];

  for (t=0;t<2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints;t++) {
    SampleAltitudeDistrubutionBuffer[t].NA_ColumnDensity=0.0;
    SampleAltitudeDistrubutionBuffer[t].NA_EmissionIntensity__5891_58A=0.0;
    SampleAltitudeDistrubutionBuffer[t].NA_EmissionIntensity__5897_56A=0.0;
  }

  //flush the sampling biffer
  for (el=0,offset=0;el<nPhaseAngleIntervals;el++) {
    SampleBuffer_AntiSunwardMotion[el].NA_ColumnDensity=0.0;
    SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_AntiSunwardMotion[el].JulianDate=0.0;
    SampleBuffer_AntiSunwardMotion[el].nSamples=0;
    SampleBuffer_AntiSunwardMotion[el].AltitudeDistrubutionBuffer=NULL;

    SampleBuffer_SunwardMotion[el].NA_ColumnDensity=0.0;
    SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_SunwardMotion[el].JulianDate=0.0;
    SampleBuffer_SunwardMotion[el].nSamples=0;
    SampleBuffer_SunwardMotion[el].AltitudeDistrubutionBuffer=NULL;

    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_ColumnDensity=0.0;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].JulianDate=0.0;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples=0;
    SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer=SampleAltitudeDistrubutionBuffer+offset++;

    SampleBuffer_SunwardMotion__TotalModelRun[el].NA_ColumnDensity=0.0;
    SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_SunwardMotion__TotalModelRun[el].JulianDate=0.0;
    SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples=0;
    SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer=SampleAltitudeDistrubutionBuffer+offset++;
  }*/
}

//output a datafile
void Comet::Sampling::SubsolarLimbColumnIntegrals::PrintDataFile() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  char fname[200];
  FILE *fout;

  if (PIC::ThisThread==0) {
    //open the output data file
    sprintf(fname,"%s/pic.Comet.LimbIntegrals.out=%i.dat",PIC::OutputDataFileDirectory,nOutputFile);
    fout=fopen(fname,"w");

    //generate the title of the file
    const SpiceInt lenout=35;
    SpiceChar utcstr[lenout];

    et2utc_c (etSampleBegin,"C",6,lenout,utcstr);
    fprintf(fout,"TITLE=\"Column Integrals of the subsolar point of the limb. Sampling Start Time=%s",utcstr);

    et2utc_c (OrbitalMotion::et,"C",6,lenout,utcstr);
    fprintf(fout,", Sampling Finit Time=%s\"\n",utcstr);

    //the variable list
    char vlist[_MAX_STRING_LENGTH_PIC_];
    int StateVectorLength;

    StateVectorLength=ColumnIntegral::GetVariableList(vlist);
    fprintf(fout,"VARIABLES=\"Phase Angle [deg]\", \"Julian day\" %s \n",vlist);


    //output the sampled data
    int el,nSampledIntervals;

    //the Comet moves toward the Sun, current cycle sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Comet moves toward the Sun, current cycle sample\",I=%i\n",nSampledIntervals);

    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion[el].nSamples!=0) {
      fprintf(fout,"%e  %e",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_SunwardMotion[el].JulianDate);
      for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",SampleBuffer_SunwardMotion[el].SampleColumnIntegrals[i]/SampleBuffer_SunwardMotion[el].nSamples);
      fprintf(fout,"\n");
    }


    //the Comet moves outward from the Sun, current cycle sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Comet moves outward from the Sun, current cycle sample\",I=%i\n",nSampledIntervals);

    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion[el].nSamples!=0) {
      fprintf(fout,"%e  %e",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_AntiSunwardMotion[el].JulianDate);
      for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",SampleBuffer_AntiSunwardMotion[el].SampleColumnIntegrals[i]/SampleBuffer_AntiSunwardMotion[el].nSamples);
      fprintf(fout,"\n");
    }


    //the Comet moves toward the Sun, total run accumulated sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Comet moves toward the Sun, total run accumulated sample\",I=%i\n",nSampledIntervals);

    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) {
      fprintf(fout,"%e  %e",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_SunwardMotion__TotalModelRun[el].JulianDate);
      for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",SampleBuffer_SunwardMotion__TotalModelRun[el].SampleColumnIntegrals[i]/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples);
      fprintf(fout,"\n");
    }

    //the Comet moves outward from the Sun, total run accumulated sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Comet moves outward from the Sun, total run accumulated sample\",I=%i\n",nSampledIntervals);

    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) {
      fprintf(fout,"%e  %e\n",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_AntiSunwardMotion__TotalModelRun[el].JulianDate);
      for (int i=0;i<StateVectorLength;i++) fprintf(fout,"   %e",SampleBuffer_AntiSunwardMotion__TotalModelRun[el].SampleColumnIntegrals[i]/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples);
      fprintf(fout,"\n");
    }


    fclose(fout);

    if (_NA_SPEC_>=0) {
      //output the altitude distribution of the integrals
      FILE *fNA_ColumnDensity,*fNA_Intensity_5891_58A,*fNA_Intensity_5897_56A;

      fNA_ColumnDensity=fopen("pic.Comet.AltitudeDistribution.NA.ColumnDensity.dat","w");
      fNA_Intensity_5891_58A=fopen("pic.Comet.AltitudeDistribution.NA.Intensity.5891_58A.dat","w");
      fNA_Intensity_5897_56A=fopen("pic.Comet.AltitudeDistribution.NA.Intensity.5897_56A.dat","w");

      fprintf(fNA_ColumnDensity,"VARIABLES=\"Altitude\"");
      fprintf(fNA_Intensity_5891_58A,"VARIABLES=\"Altitude\"");
      fprintf(fNA_Intensity_5897_56A,"VARIABLES=\"Altitude\"");

      for (el=0;el<nPhaseAngleIntervals;el++) {
        fprintf(fNA_ColumnDensity,", \"Column Density (Pa=%e,sunward)\", \"Column Density (Pa=%e,anti-sunward)\"",(el+0.5)*dPhaseAngle/Pi*180.0,(el+0.5)*dPhaseAngle/Pi*180.0);
        fprintf(fNA_Intensity_5891_58A,", \"Intensity (5891_58A,Pa=%e,sunward)\", \"Intensity (5891_58A,Pa=%e,anti-sunward)\"",(el+0.5)*dPhaseAngle/Pi*180.0,(el+0.5)*dPhaseAngle/Pi*180.0);
        fprintf(fNA_Intensity_5897_56A,", \"Intensity (5897_56A,Pa=%e,sunward)\", \"Intensity (5897_56A,Pa=%e,anti-sunward)\"",(el+0.5)*dPhaseAngle/Pi*180.0,(el+0.5)*dPhaseAngle/Pi*180.0);
      }

      fprintf(fNA_ColumnDensity,"\n");
      fprintf(fNA_Intensity_5891_58A,"\n");
      fprintf(fNA_Intensity_5897_56A,"\n");

      double R=dAltmin,dAlt=dAltmin;
      long int nAltitudeSamplePoint=0;

      while (R<maxAltitude+1.0) {
        fprintf(fNA_ColumnDensity,"%e ",R*_RADIUS_(_TARGET_));
        fprintf(fNA_Intensity_5891_58A,"%e ",R*_RADIUS_(_TARGET_));
        fprintf(fNA_Intensity_5897_56A,"%e ",R*_RADIUS_(_TARGET_));

        for (el=0;el<nPhaseAngleIntervals;el++) {
          fprintf(fNA_ColumnDensity,"%e %e ",
              (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_COLUMN_DENSITY_OFFSET_]/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
              (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_COLUMN_DENSITY_OFFSET_]/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
              );

          fprintf(fNA_Intensity_5891_58A,"%e %e ",
              (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_EMISSION_5891_58A_SAMPLE_OFFSET_]/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
              (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_EMISSION_5891_58A_SAMPLE_OFFSET_]/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
              );

          fprintf(fNA_Intensity_5897_56A,"%e %e ",
              (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_EMISSION_5897_56A_SAMPLE_OFFSET_]/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
              (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][_NA_EMISSION_5897_56A_SAMPLE_OFFSET_]/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
              );
        }


        fprintf(fNA_ColumnDensity,"\n");
        fprintf(fNA_Intensity_5891_58A,"\n");
        fprintf(fNA_Intensity_5897_56A,"\n");

        R+=dAlt;
        dAlt*=rr;
        if (++nAltitudeSamplePoint>=nSampleAltitudeDistributionPoints) break;
      }

      fclose(fNA_ColumnDensity);
      fclose(fNA_Intensity_5891_58A);
      fclose(fNA_Intensity_5897_56A);
    }
  }

  //reset the sampling buffer
  int i,el,StateVectorLength=ColumnIntegral::GetVariableList(NULL);

  for (el=0;el<nPhaseAngleIntervals;el++) {
    SampleBuffer_AntiSunwardMotion[el].JulianDate=0.0;
    SampleBuffer_AntiSunwardMotion[el].nSamples=0;
    for (i=0;i<StateVectorLength;i++) SampleBuffer_AntiSunwardMotion[el].SampleColumnIntegrals[i]=0.0;

    SampleBuffer_SunwardMotion[el].JulianDate=0.0;
    SampleBuffer_SunwardMotion[el].nSamples=0;
    for (i=0;i<StateVectorLength;i++) SampleBuffer_SunwardMotion[el].SampleColumnIntegrals[i]=0.0;

    if (nOutputFile==0) {
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].JulianDate=0.0;
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples=0;
      for (i=0;i<StateVectorLength;i++) SampleBuffer_AntiSunwardMotion__TotalModelRun[el].SampleColumnIntegrals[i]=0.0;

      SampleBuffer_SunwardMotion__TotalModelRun[el].JulianDate=0.0;
      SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples=0;
      for (i=0;i<StateVectorLength;i++) SampleBuffer_SunwardMotion__TotalModelRun[el].SampleColumnIntegrals[i]=0.0;

      for (int t=0;t<nSampleAltitudeDistributionPoints;t++) {
        for (i=0;i<StateVectorLength;i++) {
          SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[t][i]=0.0;
          SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistributionColumnIntegrals[t][i]=0.0;
        }
      }

    }
  }


  etSampleBegin=OrbitalMotion::et;
  nOutputFile++;
#endif
}

void Comet::Sampling::SubsolarLimbColumnIntegrals::CollectSample(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int CometMotionDirection;

  //get the direction of the lunar motion in respect to the Sun
  SpiceDouble lt,SunState[6];
  spkezr_c("Sun",Exosphere::OrbitalMotion::et,SO_FRAME,"none","Comet",SunState,&lt);
  CometMotionDirection=(SunState[3]>0.0) ? -1 : 1;   //-1 -> the Comet moves toward the Sun, +1 -> the moon moves outward from the Sun



  //get the direction of motion from projection of the lunar position on the ecliptic
  SpiceDouble CometState[6];
  spkezr_c("Comet",Exosphere::OrbitalMotion::et,"GSE","none","Earth",CometState,&lt);
  CometMotionDirection=(CometState[1]>0.0) ? -1 : 1;



  //init the sampling the buffer at the first call of 'CollectSample'
  if (etSampleBegin==0.0) {
    etSampleBegin=OrbitalMotion::et;
    SamplingPhase=0;
    firstPhaseRadialVelocityDirection=CometMotionDirection;
  }


  //check if the moon has finished the full rotation around the Earth
  if ((firstPhaseRadialVelocityDirection==CometMotionDirection)&&(SamplingPhase==1)) {
    //the Comet has finished the full rotation around the Earth -> output the data file
    PrintDataFile();
    SamplingPhase=0;
  }
  else if ((firstPhaseRadialVelocityDirection!=CometMotionDirection)&&(SamplingPhase==0)) {
    //the Comet have changed the direction of motion in respect to the Sun
    SamplingPhase=1;
  }

  //sample the integrals
  int IntegralsVectorLength=ColumnIntegral::GetVariableList(NULL);
  double IntegralsVector[IntegralsVectorLength];
  double LimbDirection[3],EarthPosition[3];
  ColumnIntegral::GetSubsolarPointDirection(LimbDirection,EarthPosition);

  PIC::ColumnIntegration::GetCoulumnIntegral(IntegralsVector,IntegralsVectorLength,EarthPosition,LimbDirection,ColumnIntegral::CoulumnDensityIntegrant);
  ColumnIntegral::ProcessColumnIntegrationVector(IntegralsVector,IntegralsVectorLength);

  //determinethe phase angle and save the data
  double PhaseAngle;
  int i,iPhaseAngle;

  PhaseAngle=Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et);
  iPhaseAngle=int(PhaseAngle/dPhaseAngle);

  if (CometMotionDirection==1) {
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].nSamples++;
    for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_AntiSunwardMotion[iPhaseAngle].SampleColumnIntegrals[i]+=IntegralsVector[i];

    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].nSamples++;
    for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].SampleColumnIntegrals[i]+=IntegralsVector[i];
  }
  else {
    SampleBuffer_SunwardMotion[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_SunwardMotion[iPhaseAngle].nSamples++;
    for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_SunwardMotion[iPhaseAngle].SampleColumnIntegrals[i]+=IntegralsVector[i];

    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].nSamples++;
    for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].SampleColumnIntegrals[i]+=IntegralsVector[i];
  }



  SpiceInt n,nxpts;
  SpiceDouble rad[3],xpt0[3],xpt1[3],xLimb_IAU[3],xLimb_SO[3];
  SpiceDouble EarthState_IAU[6],EarthState_SO[6],SunState_IAU[6];
  SpiceEllipse limb;
  SpicePlane plane;

  //find the limb
  bodvrd_c(ObjectName, "RADII", 3, &n, rad );
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,IAU_FRAME,"none",ObjectName,EarthState_IAU,&lt);
  edlimb_c(rad[0],rad[1],rad[2],EarthState_IAU,&limb);

  //find intersection of the limb with the plane hat contains position of the Earth, center of the boly and the Sun
  spkezr_c("SUN",Exosphere::OrbitalMotion::et,IAU_FRAME,"none","EARTH",SunState_IAU,&lt);
  psv2pl_c (EarthState_IAU,EarthState_IAU,SunState_IAU,&plane);
  inelpl_c (&limb,&plane,&nxpts,xpt0,xpt1);

  //choose the limb
  double cosXpt0=0.0,cosXpt1=0.0;
  int idim=0;

  if (nxpts<2) {
    exit(__LINE__,__FILE__,"Error: cannot file point of intersection of the plane that contains centers of the Earth, the Object and the Sun with the elips of the limb");
  }

  for (idim=0;idim<3;idim++) {
    cosXpt0+=(xpt0[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
    cosXpt1+=(xpt1[idim]-EarthState_IAU[idim])*(SunState_IAU[idim]-EarthState_IAU[idim]);
  }

  if (cosXpt0>cosXpt1) memcpy(xLimb_IAU,xpt0,3*sizeof(SpiceDouble)); //'xpt0' is closer to the Sun
  else memcpy(xLimb_IAU,xpt1,3*sizeof(SpiceDouble)); //'xpt1' is closer to the Sun


  //recalcualte position of the Earth in the SO_FRAME
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,SO_FRAME,"none",ObjectName,EarthState_SO,&lt);

  //sample the altitude distribution of the column integrals
  double R=dAltmin,dAlt=dAltmin;
  long int nAltitudeSamplePoint=0;
  double l[3],rEarth_SO[3],c=0.0;

  while (R<maxAltitude+1.0) {
    for (idim=0,c=0.0;idim<3;idim++) {
      //recalculate position of the limb in the 'SO_FRAME'
      xLimb_SO[idim]=
          (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][0]*xLimb_IAU[0]*(R+1.0))+
          (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][1]*xLimb_IAU[1]*(R+1.0))+
          (OrbitalMotion::IAU_to_SO_TransformationMartix[idim][2]*xLimb_IAU[2]*(R+1.0));



      l[idim]=xLimb_SO[idim]-EarthState_SO[idim];
      c+=pow(l[idim],2);

      rEarth_SO[idim]=1.0E3*EarthState_SO[idim];
    }

    for (c=sqrt(c),idim=0;idim<3;idim++) l[idim]/=c;

    PIC::ColumnIntegration::GetCoulumnIntegral(IntegralsVector,IntegralsVectorLength,rEarth_SO,l,ColumnIntegral::CoulumnDensityIntegrant);
    ColumnIntegral::ProcessColumnIntegrationVector(IntegralsVector,IntegralsVectorLength);

    if (CometMotionDirection==1) {
 /*     SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity+=IntegralsVector[0];
      SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
      SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
 */
      for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][i]+=IntegralsVector[i];
    }
    else {
 /*     SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity+=IntegralsVector[0];
      SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
      SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
 */
      for (i=0;i<IntegralsVectorLength;i++) SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistributionColumnIntegrals[nAltitudeSamplePoint][i]+=IntegralsVector[i];
    }

    R+=dAlt;
    dAlt*=rr;
    if (++nAltitudeSamplePoint>=nSampleAltitudeDistributionPoints) break;
  }

#endif
}








//print into a file column integrals that correspond to observations of Kaguya's TVIS instrument
/* void Comet::Sampling::Kaguya::TVIS::OutputModelData(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int i;
  FILE *fout;

  int IntegralsVectorLength=ColumnIntegral::GetVariableList(NULL);
  double IntegralsVector[IntegralsVectorLength];

  if (_NA_SPEC_>=0) for (unsigned int nDataSet=0;nDataSet<TvisOrientationVector.size();nDataSet++) {
    int nPass=0;

    if (PIC::ThisThread==0) {
      char fname[100];

      sprintf(fname,"%s/pic.Comet.Kaguya.TVIS.set=%i.out=%i.dat",PIC::OutputDataFileDirectory,nDataSet,DataOutputFileNumber);
      fout=fopen(fname,"w");

      fprintf(fout,"VARIABLES= \"Equatorial Crossing Time\", \"Epoch\", \"NA ColumnDensity\", \"NA EmissionIntensity (5891_58A)\", \"NA EmissionIntensity (5897_56A)\"\n");
      fprintf(fout,"ZONE T=\"Pass=%i\"\n",nPass++);
    }

    //output the column integrals
    for (i=0;i<TvisOrientationVector[nDataSet].nTotalTvisOrientationElements;i++) {
      PIC::ColumnIntegration::GetCoulumnIntegral(IntegralsVector,IntegralsVectorLength,TvisOrientationVector[nDataSet].TvisOrientation[i].scPosition,TvisOrientationVector[nDataSet].TvisOrientation[i].fov,ColumnIntegral::CoulumnDensityIntegrant);

      if (PIC::ThisThread==0) {
        if (i!=0) if ((TvisOrientationVector[nDataSet].TvisOrientation[i].EquatorialCrossingTime<0.0)&&(TvisOrientationVector[nDataSet].TvisOrientation[i-1].EquatorialCrossingTime>0.0)) fprintf(fout,"ZONE T=\"Pass=%i\"\n",nPass++);

        fprintf(fout,"%e  %e  %e  %e  %e\n",TvisOrientationVector[nDataSet].TvisOrientation[i].EquatorialCrossingTime,TvisOrientationVector[nDataSet].TvisOrientation[i].Epoch,
            IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_],IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_],IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_]);
      }
    }

    if (PIC::ThisThread==0) fclose(fout);
  }

  //get the sampled information for the individual points
  if (nTotalIndividualPointsTvisOrientationElements!=0) {
    if (PIC::ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"\n\nUPI_TVIS Individual point sample:\nVARIABLES= \"UTC\", \"Equatorial Crossing Time\", \"Epoch\", \"NA ColumnDensity\", \"NA EmissionIntensity (5891_58A)\", \"NA EmissionIntensity (5897_56A)\"\n");

    for (i=0;i<nTotalIndividualPointsTvisOrientationElements;i++) {
      PIC::ColumnIntegration::GetCoulumnIntegral(IntegralsVector,IntegralsVectorLength,individualPointsTvisOrientation[i].scPosition,TVIS::individualPointsTvisOrientation[i].fov,ColumnIntegral::CoulumnDensityIntegrant);

      if (PIC::ThisThread==0) {
        fprintf(PIC::DiagnospticMessageStream,"%s %e  %e  %e  %e  %e\n",individualPointsTvisOrientation[i].UTC,individualPointsTvisOrientation[i].EquatorialCrossingTime,individualPointsTvisOrientation[i].Epoch,
            IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_COLUMN_DENSITY_OFFSET_],IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5891_58A_SAMPLE_OFFSET_],IntegralsVector[Comet::Sampling::SubsolarLimbColumnIntegrals::_NA_EMISSION_5897_56A_SAMPLE_OFFSET_]);
      }
    }
  }

#endif
}*/
