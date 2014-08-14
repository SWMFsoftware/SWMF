//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
/*
 * mercury_ReferenceGroundBasedObservations.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: vtenishe
 */

//$Id$

#include "pic.h"
#include "Exosphere.h"

Exosphere::Sampling::ReferenceGroundBasedObservations::cObservationTag Exosphere::Sampling::ReferenceGroundBasedObservations::RemoteObservationList[nReferenceGroundBasedObservations];

void Exosphere::Sampling::ReferenceGroundBasedObservations::init() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int n;

  if (PIC::ThisThread==0) cout << "$PREFIX:Reference Ground Based Observations:" << endl << "$PREFIX:n\tObservation Time[degrees]\tTAA\tPhase Angle[degrees]" << endl;

  for (n=0;n<nReferenceGroundBasedObservations;n++) {
    RemoteObservationList[n].nOutputFile=0;

    sprintf(RemoteObservationList[n].TimeStamp,"%s",ReferenceGroundBasedObservationTime[n]);
    utc2et_c(RemoteObservationList[n].TimeStamp,&RemoteObservationList[n].et);
    RemoteObservationList[n].TAA=Exosphere::OrbitalMotion::GetTAA(RemoteObservationList[n].et);
    RemoteObservationList[n].PhaseAngle=Exosphere::OrbitalMotion::GetPhaseAngle(RemoteObservationList[n].et);

    if (PIC::ThisThread==0) cout << "$PREFIX:" << n << "\t" << ReferenceGroundBasedObservationTime[n] << "\t" << RemoteObservationList[n].TAA/Pi*180.0 << "\t" << RemoteObservationList[n].PhaseAngle/Pi*180.0 << endl;
  }
#endif
}

void Exosphere::Sampling::ReferenceGroundBasedObservations::OutputSampledData(SpiceDouble etStartInterval,SpiceDouble etFinishInterval,int nMercuryOutputFile) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  double dTAA,dTAAobservation,taaStartInterval,taaFinishInterval,taa;
  double paStartInterval,paFinishInterval;
  int n;

  double domainCharacteristicSize=0.0;
  for (int idim=0;idim<DIM;idim++) domainCharacteristicSize=max(max(fabs(PIC::Mesh::mesh.xGlobalMax[idim]),fabs(PIC::Mesh::mesh.xGlobalMin[idim])),domainCharacteristicSize);

  taaStartInterval=Exosphere::OrbitalMotion::GetTAA(etStartInterval);
  taaFinishInterval=Exosphere::OrbitalMotion::GetTAA(etFinishInterval);
  dTAA=(taaStartInterval<taaFinishInterval) ? taaFinishInterval-taaStartInterval : taaFinishInterval+2.0*Pi-taaStartInterval;

  paStartInterval=Exosphere::OrbitalMotion::GetPhaseAngle(etStartInterval);
  paFinishInterval=Exosphere::OrbitalMotion::GetPhaseAngle(etFinishInterval);

  for (n=0;n<nReferenceGroundBasedObservations;n++) {
    taa=RemoteObservationList[n].TAA;
    dTAAobservation=(taaStartInterval<taa) ? taa-taaStartInterval : taa+2.0*Pi-taaStartInterval;

    if ((dTAAobservation<dTAA) || ((min(paStartInterval,paFinishInterval)<=RemoteObservationList[n].PhaseAngle)&&(RemoteObservationList[n].PhaseAngle<=max(paStartInterval,paFinishInterval))) ) {
      //the taa of the remote observation falls in the range taaStartInterval - taaFinishInterval
      char fname[_MAX_STRING_LENGTH_PIC_];

      sprintf(fname,"%s/pic.GroundBasedObservation.n=%i.ColumnDensityMap.%s.nMercuryOutputFile=%i.out=%i.dat",PIC::OutputDataFileDirectory,n,RemoteObservationList[n].TimeStamp,nMercuryOutputFile,RemoteObservationList[n].nOutputFile);
      Exosphere::ColumnIntegral::CircularMap(fname,domainCharacteristicSize,0.05*_RADIUS_(_TARGET_),5*_RADIUS_(_TARGET_),80,RemoteObservationList[n].et);

      sprintf(fname,"%s/pic.GroundBasedObservation.n=%i.LimbColumnDensity.%s.nMercuryOutputFile=%i.out=%i.dat",PIC::OutputDataFileDirectory,n,RemoteObservationList[n].TimeStamp,nMercuryOutputFile,RemoteObservationList[n].nOutputFile);
      Exosphere::ColumnIntegral::Limb(fname);

      RemoteObservationList[n].nOutputFile++;
    }
  }
#endif
}
