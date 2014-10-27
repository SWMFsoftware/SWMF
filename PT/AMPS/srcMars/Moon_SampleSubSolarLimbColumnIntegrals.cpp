/*
 * Moon_SampleSubSolarLimbColumnIntegrals.cpp
 *
 *  Created on: Nov 9, 2012
 *      Author: vtenishe
 */

//functions that sample the column integrals at the subsolar point of the limb
//$Id$

#include "pic.h"

SpiceDouble Moon::Sampling::SubsolarLimbColumnIntegrals::etSampleBegin;
int Moon::Sampling::SubsolarLimbColumnIntegrals::SamplingPhase;
int Moon::Sampling::SubsolarLimbColumnIntegrals::firstPhaseRadialVelocityDirection;
int Moon::Sampling::SubsolarLimbColumnIntegrals::nOutputFile;

double Moon::Sampling::SubsolarLimbColumnIntegrals::rr=0.0;
long int Moon::Sampling::SubsolarLimbColumnIntegrals::nSampleAltitudeDistributionPoints=0;
Moon::Sampling::SubsolarLimbColumnIntegrals::cSampleAltitudeDistrubutionBufferElement *Moon::Sampling::SubsolarLimbColumnIntegrals::SampleAltitudeDistrubutionBuffer=NULL;

Moon::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Moon::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_AntiSunwardMotion[nPhaseAngleIntervals];
Moon::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Moon::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_SunwardMotion[nPhaseAngleIntervals];
Moon::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Moon::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_AntiSunwardMotion__TotalModelRun[nPhaseAngleIntervals];
Moon::Sampling::SubsolarLimbColumnIntegrals::cSampleBufferElement Moon::Sampling::SubsolarLimbColumnIntegrals::SampleBuffer_SunwardMotion__TotalModelRun[nPhaseAngleIntervals];

void Moon::Sampling::SubsolarLimbColumnIntegrals::EmptyFunction() {}

//init the sampling block
void Moon::Sampling::SubsolarLimbColumnIntegrals::init() {
  int el;

  etSampleBegin=0.0;
  SamplingPhase=-1;
  firstPhaseRadialVelocityDirection=0;
  nOutputFile=0;

  //init the altitude distribution sampling buffer
  double R,dAlt;
  long int t,offset;

  rr=(maxAltitude+dAltmax)/(maxAltitude+dAltmin);
  t=(long int)(log(dAltmax/dAltmin)/log(rr)-2.0);
  rr=pow(dAltmax/dAltmin,1.0/(t+2.0));

  R=dAltmin,dAlt=dAltmin;
  while (R<maxAltitude+1.0) {
    ++nSampleAltitudeDistributionPoints;

    R+=dAlt;
    dAlt*=rr;
  }

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
  }
}

//output a datafile
void Moon::Sampling::SubsolarLimbColumnIntegrals::PrintDataFile() {
  char fname[200];
  FILE *fout;

  if (PIC::ThisThread==0) {
    //open the output data file
    sprintf(fname,"pic.Moon.LimbIntegrals.out=%i.dat",nOutputFile);
    fout=fopen(fname,"w");

    //generate the title of the file
    const SpiceInt lenout=35;
    SpiceChar utcstr[lenout];

    et2utc_c (etSampleBegin,"C",6,lenout,utcstr);
    fprintf(fout,"TITLE=\"Column Integrals of the subsolar point of the limb. Sampling Start Time=%s",utcstr);

    et2utc_c (OrbitalMotion::et,"C",6,lenout,utcstr);
    fprintf(fout,", Sampling Finit Time=%s\"\n",utcstr);

    //the variable list
    fprintf(fout,"VARIABLES=\"Phase Angle [deg]\", \"Julian day\",      \"Column Integral[m^-2]\", \"Intensity(5891_58A) [kR]\", \"Intensity(5897_56A) [kR]\" \n");


    //output the sampled data
    int el,nSampledIntervals;

    //the Moon moves toward the Sun, current cycle sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Moon moves toward the Sun, current cycle sample\",I=%i\n",nSampledIntervals);
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion[el].nSamples!=0) {
       fprintf(fout,"%e  %e  %e  %e  %e\n",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_SunwardMotion[el].JulianDate,
           SampleBuffer_SunwardMotion[el].NA_ColumnDensity/SampleBuffer_SunwardMotion[el].nSamples,
           SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5891_58A/SampleBuffer_SunwardMotion[el].nSamples,
           SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5897_56A/SampleBuffer_SunwardMotion[el].nSamples);
    }


    //the Moon moves outward from the Sun, current cycle sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Moon moves outward from the Sun, current cycle sample\",I=%i\n",nSampledIntervals);
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion[el].nSamples!=0) {
       fprintf(fout,"%e  %e  %e  %e  %e\n",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_AntiSunwardMotion[el].JulianDate,
           SampleBuffer_AntiSunwardMotion[el].NA_ColumnDensity/SampleBuffer_AntiSunwardMotion[el].nSamples,
           SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5891_58A/SampleBuffer_AntiSunwardMotion[el].nSamples,
           SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5897_56A/SampleBuffer_AntiSunwardMotion[el].nSamples);
    }


    //the Moon moves toward the Sun, total run accumulated sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Moon moves toward the Sun, total run accumulated sample\",I=%i\n",nSampledIntervals);
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) {
       fprintf(fout,"%e  %e  %e  %e  %e\n",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_SunwardMotion__TotalModelRun[el].JulianDate,
           SampleBuffer_SunwardMotion__TotalModelRun[el].NA_ColumnDensity/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples,
           SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples,
           SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples);
    }

    //the Moon moves outward from the Sun, total run accumulated sample
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) nSampledIntervals++;
    fprintf(fout,"ZONE T=\"the Moon moves outward from the Sun, total run accumulated sample\",I=%i\n",nSampledIntervals);
    for (el=0,nSampledIntervals=0;el<nPhaseAngleIntervals;el++) if (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) {
       fprintf(fout,"%e  %e  %e  %e  %e\n",(el+0.5)*dPhaseAngle/Pi*180.0,SampleBuffer_AntiSunwardMotion__TotalModelRun[el].JulianDate,
           SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_ColumnDensity/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples,
           SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples,
           SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples);
    }


    fclose(fout);

    //output the altitude distribution of the integrals
    FILE *fNA_ColumnDensity,*fNA_Intensity_5891_58A,*fNA_Intensity_5897_56A;

    fNA_ColumnDensity=fopen("pic.Moon.AltitudeDistribution.NA.ColumnDensity.dat","w");
    fNA_Intensity_5891_58A=fopen("pic.Moon.AltitudeDistribution.NA.Intensity.5891_58A.dat","w");
    fNA_Intensity_5897_56A=fopen("pic.Moon.AltitudeDistribution.NA.Intensity.5897_56A.dat","w");

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
            (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
            (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
            );

        fprintf(fNA_Intensity_5891_58A,"%e %e ",
            (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
            (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
            );

        fprintf(fNA_Intensity_5897_56A,"%e %e ",
            (SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_SunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A/SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples : 0.0,
            (SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples!=0) ? SampleBuffer_AntiSunwardMotion__TotalModelRun[el].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A/SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples : 0.0
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

  //reset the sampling buffer
  int el;

  for (el=0;el<nPhaseAngleIntervals;el++) {
    SampleBuffer_AntiSunwardMotion[el].NA_ColumnDensity=0.0;
    SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_AntiSunwardMotion[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_AntiSunwardMotion[el].JulianDate=0.0;
    SampleBuffer_AntiSunwardMotion[el].nSamples=0;

    SampleBuffer_SunwardMotion[el].NA_ColumnDensity=0.0;
    SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5891_58A=0.0;
    SampleBuffer_SunwardMotion[el].NA_EmissionIntensity__5897_56A=0.0;
    SampleBuffer_SunwardMotion[el].JulianDate=0.0;
    SampleBuffer_SunwardMotion[el].nSamples=0;

    if (nOutputFile==0) {
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_ColumnDensity=0.0;
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A=0.0;
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A=0.0;
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].JulianDate=0.0;
      SampleBuffer_AntiSunwardMotion__TotalModelRun[el].nSamples=0;

      SampleBuffer_SunwardMotion__TotalModelRun[el].NA_ColumnDensity=0.0;
      SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5891_58A=0.0;
      SampleBuffer_SunwardMotion__TotalModelRun[el].NA_EmissionIntensity__5897_56A=0.0;
      SampleBuffer_SunwardMotion__TotalModelRun[el].JulianDate=0.0;
      SampleBuffer_SunwardMotion__TotalModelRun[el].nSamples=0;
    }
  }

  if (nOutputFile==0) for (long int t=0;t<2*nPhaseAngleIntervals*nSampleAltitudeDistributionPoints;t++) {
    SampleAltitudeDistrubutionBuffer[t].NA_ColumnDensity=0.0;
    SampleAltitudeDistrubutionBuffer[t].NA_EmissionIntensity__5891_58A=0.0;
    SampleAltitudeDistrubutionBuffer[t].NA_EmissionIntensity__5897_56A=0.0;
  }

  etSampleBegin=OrbitalMotion::et;
  nOutputFile++;
}

void Moon::Sampling::SubsolarLimbColumnIntegrals::CollectSample(int DataOutputFileNumber) {
  int MoonMotionDirection;

  //get the direction of the lunar motion in respect to the Sun
  SpiceDouble lt,SunState[6];
  spkezr_c("Sun",Exosphere::OrbitalMotion::et,SO_FRAME,"none","Moon",SunState,&lt);
  MoonMotionDirection=(SunState[3]>0.0) ? -1 : 1;   //-1 -> the Moon moves toward the Sun, +1 -> the moon moves outward from the Sun



  //get the direction of motion from projection of the lunar position on the ecliptic
  SpiceDouble MoonState[6];
  spkezr_c("Moon",Exosphere::OrbitalMotion::et,"GSE","none","Earth",MoonState,&lt);
  MoonMotionDirection=(MoonState[1]>0.0) ? -1 : 1;



  //init the sampling the buffer at the first call of 'CollectSample'
  if (etSampleBegin==0.0) {
    etSampleBegin=OrbitalMotion::et;
    SamplingPhase=0;
    firstPhaseRadialVelocityDirection=MoonMotionDirection;
  }


  //check if the moon has finished the full rotation around the Earth
  if ((firstPhaseRadialVelocityDirection==MoonMotionDirection)&&(SamplingPhase==1)) {
    //the Moon has finished the full rotation around the Earth -> output the data file
    PrintDataFile();
    SamplingPhase=0;
  }
  else if ((firstPhaseRadialVelocityDirection!=MoonMotionDirection)&&(SamplingPhase==0)) {
    //the Moon have changed the direction of motion in respect to the Sun
    SamplingPhase=1;
  }

  //sample the integrals
  const int IntegralsVectorLength=3;
  double IntegralsVector[IntegralsVectorLength];

  Exosphere::GetColumnIntegralLimbSubsolarPoint(IntegralsVector,IntegralsVectorLength);

  //determinethe phase angle and save the data
  double PhaseAngle;
  int iPhaseAngle;

  PhaseAngle=Exosphere::OrbitalMotion::GetPhaseAngle(Exosphere::OrbitalMotion::et);
  iPhaseAngle=int(PhaseAngle/dPhaseAngle);

  if (MoonMotionDirection==1) {
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].NA_ColumnDensity+=IntegralsVector[0];
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_AntiSunwardMotion[iPhaseAngle].nSamples++;

    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].NA_ColumnDensity+=IntegralsVector[0];
    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].nSamples++;
  }
  else {
    SampleBuffer_SunwardMotion[iPhaseAngle].NA_ColumnDensity+=IntegralsVector[0];
    SampleBuffer_SunwardMotion[iPhaseAngle].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
    SampleBuffer_SunwardMotion[iPhaseAngle].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    SampleBuffer_SunwardMotion[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_SunwardMotion[iPhaseAngle].nSamples++;

    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].NA_ColumnDensity+=IntegralsVector[0];
    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].JulianDate=j2000_c()+Exosphere::OrbitalMotion::et/spd_c();
    SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].nSamples++;
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

    PIC::ColumnIntegration::GetCoulumnIntegral(IntegralsVector,IntegralsVectorLength,rEarth_SO,l,SodiumCoulumnDensityIntegrant);

    if (MoonMotionDirection==1) {
      SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity+=IntegralsVector[0];
      SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
      SampleBuffer_AntiSunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    }
    else {
      SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_ColumnDensity+=IntegralsVector[0];
      SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5891_58A+=IntegralsVector[1];
      SampleBuffer_SunwardMotion__TotalModelRun[iPhaseAngle].AltitudeDistrubutionBuffer[nAltitudeSamplePoint].NA_EmissionIntensity__5897_56A+=IntegralsVector[2];
    }

    R+=dAlt;
    dAlt*=rr;
    if (++nAltitudeSamplePoint>=nSampleAltitudeDistributionPoints) break;
  }


}
