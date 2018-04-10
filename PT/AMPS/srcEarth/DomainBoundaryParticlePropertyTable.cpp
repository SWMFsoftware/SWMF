//$Id$
//functionality of sampling of the distribution of the particle velocity direction and energy that car reach the near Earth' regions


/*
 * DomainBoundaryParticlePropertyTable.cpp
 *
 *  Created on: Sep 16, 2017
 *      Author: vtenishe
 */

#include "pic.h"
#include "specfunc.h"

double Earth::CutoffRigidity::DomainBoundaryParticleProperty::logEmax=log(1.0E1*GeV2J);
double Earth::CutoffRigidity::DomainBoundaryParticleProperty::logEmin=log(1.0E1*MeV2J);

int Earth::CutoffRigidity::DomainBoundaryParticleProperty::nAzimuthIntervals=10;
int Earth::CutoffRigidity::DomainBoundaryParticleProperty::nCosZenithIntervals=10;
int Earth::CutoffRigidity::DomainBoundaryParticleProperty::nLogEnergyLevels=10;

double Earth::CutoffRigidity::DomainBoundaryParticleProperty::dCosZenithAngle=1.0/Earth::CutoffRigidity::DomainBoundaryParticleProperty::nCosZenithIntervals;
double Earth::CutoffRigidity::DomainBoundaryParticleProperty::dAzimuthAngle=2.0*Pi/Earth::CutoffRigidity::DomainBoundaryParticleProperty::nAzimuthIntervals;
double Earth::CutoffRigidity::DomainBoundaryParticleProperty::dLogE=(logEmax-logEmin)/Earth::CutoffRigidity::DomainBoundaryParticleProperty::nLogEnergyLevels;

double Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX[6][2];
cBitwiseFlagTable *****Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable=NULL;
//cBitwiseFlagTable Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable[PIC::nTotalSpecies][6][Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection][Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection];

//parameters of the phase space sampling procedure
//ActiveFlag defines whether sampling of the phaswe space is turned on
//LastActiveOutputCycleNumber is the number of output cycles during which the phase space is sampled
//LastActiveOutputCycleResetParticleBuffer is the flag that defined whether the particle buffer need to be reset at the end of sampling of the phase space
bool Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::ActiveFlag=false;
int Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::LastActiveOutputCycleNumber=-1;
bool Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::LastActiveOutputCycleResetParticleBuffer=false;
double Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::OccupiedPhaseSpaceTableFraction=0.2;

bool Earth::CutoffRigidity::DomainBoundaryParticleProperty::ApplyInjectionPhaseSpaceLimiting=false;
bool Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;


//allocate the particle property sample table
void Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate() {
  int s,iface,i,j,iThreadOpenMP,iTestLocation;

  int offset,nlocs=std::max(1,Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength);

  SampleTable=new cBitwiseFlagTable**** [nlocs];
  SampleTable[0]=new cBitwiseFlagTable*** [PIC::nTotalSpecies*nlocs];

  for (offset=0,iTestLocation=0;iTestLocation<nlocs;iTestLocation++) {
    SampleTable[iTestLocation]=SampleTable[0]+offset;
    offset+=PIC::nTotalSpecies;
  }

  SampleTable[0][0]=new cBitwiseFlagTable** [6*PIC::nTotalSpecies*nlocs];

  for (offset=0,iTestLocation=0;iTestLocation<nlocs;iTestLocation++) for (s=0;s<PIC::nTotalSpecies;s++) {
    SampleTable[iTestLocation][s]=SampleTable[0][0]+offset;
    offset+=6;
  }

  SampleTable[0][0][0]=new cBitwiseFlagTable* [6*SampleMaskNumberPerSpatialDirection*PIC::nTotalSpecies*nlocs];

  for (offset=0,iTestLocation=0;iTestLocation<nlocs;iTestLocation++) for (s=0;s<PIC::nTotalSpecies;s++) for (iface=0;iface<6;iface++) {
    SampleTable[iTestLocation][s][iface]=SampleTable[0][0][0]+offset;
    offset+=SampleMaskNumberPerSpatialDirection;
  }

  SampleTable[0][0][0][0]=new cBitwiseFlagTable [6*SampleMaskNumberPerSpatialDirection*SampleMaskNumberPerSpatialDirection*PIC::nTotalSpecies*nlocs];

  for (offset=0,iTestLocation=0;iTestLocation<nlocs;iTestLocation++) for (s=0;s<PIC::nTotalSpecies;s++) for (iface=0;iface<6;iface++) for (i=0;i<SampleMaskNumberPerSpatialDirection;i++) {
    SampleTable[iTestLocation][s][iface][i]=SampleTable[0][0][0][0]+offset;
    offset+=SampleMaskNumberPerSpatialDirection;
  }

  for (iTestLocation=0;iTestLocation<nlocs;iTestLocation++) {
    for (s=0;s<PIC::nTotalSpecies;s++) for (iface=0;iface<6;iface++) for (i=0;i<SampleMaskNumberPerSpatialDirection;i++) for (j=0;j<SampleMaskNumberPerSpatialDirection;j++) {
      for (iThreadOpenMP=0;iThreadOpenMP<SampleTable[iTestLocation][s][iface][i][j].nThreadsOpenMP;iThreadOpenMP++) {
        SampleTable[iTestLocation][s][iface][i][j].AllocateTable(nAzimuthIntervals*nCosZenithIntervals*nLogEnergyLevels,iThreadOpenMP);
      }
    }
  }
}

//get the global index that corresponds to the velocity vector
int Earth::CutoffRigidity::DomainBoundaryParticleProperty::GetVelocityVectorIndex(int spec,double *v,int iface) {
  double t,l[3],Energy,Speed=0.0;
  int idim,iCosZenithInterval,iAzimuthInterval,iLogEnergyLevel;

  //getermine the speed
  for (idim=0;idim<3;idim++) Speed+=v[idim]*v[idim];

  Speed=sqrt(Speed);
  t=1.0/Speed;

  for (idim=0;idim<3;idim++) l[idim]=t*v[idim];

  //getermine the energy level
  Energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));
  iLogEnergyLevel=(int)(log(Energy)-logEmin)/nLogEnergyLevels;
  if (iLogEnergyLevel<0) iLogEnergyLevel=0;
  if (iLogEnergyLevel>=nLogEnergyLevels) iLogEnergyLevel=nLogEnergyLevels-1;

  //determine the azimuth angle level
  double c0=0.0,c1=0.0,c2=0.0;

  for (int i=0;i<3;i++) c0-=v[i]*e0FaceFrame[iface][i],c1-=v[i]*e1FaceFrame[iface][i],c2-=v[i]*InternalFaceNorm[iface][i];  //the reversed direction of the velocity vector is important

  //determine the zenith angle level
  if (c2>0.0) {
    iCosZenithInterval=c2/(Speed*dCosZenithAngle);
    if (iCosZenithInterval>=nCosZenithIntervals) iCosZenithInterval=nCosZenithIntervals-1;
  }
  else iCosZenithInterval=nCosZenithIntervals-1;

  //determine the azimuth angle level
  double phi=acos(c0/Speed);

  if (c1<0.0) phi=PiTimes2-phi;
  iAzimuthInterval=(int)(phi/dAzimuthAngle);
  if (iAzimuthInterval>=nAzimuthIntervals) iAzimuthInterval=nAzimuthIntervals-1;

  //return the index corresponding to the particle velocity
  return iAzimuthInterval+(iCosZenithInterval+iLogEnergyLevel*nCosZenithIntervals)*nAzimuthIntervals;
}

void Earth::CutoffRigidity::DomainBoundaryParticleProperty::ConvertVelocityVectorIndex2Velocity(int spec,double *v,int iface,int Index) {
  int idim,iCosZenithInterval,iAzimuthInterval,iLogEnergyLevel;

  iAzimuthInterval=Index%nAzimuthIntervals;
  Index/=nAzimuthIntervals;

  iCosZenithInterval=Index%nCosZenithIntervals;
  iLogEnergyLevel=Index/nCosZenithIntervals;

  double Speed,Energy,ZenithAngle,AzimuthAngle,cosZenithAngle,sinZenithAngle,cosAzimuthAngle,sinAzimuthAngle;

  Energy=exp((0.5+iLogEnergyLevel)*nLogEnergyLevels+logEmin);
  ZenithAngle=acos((iCosZenithInterval+0.5)*dCosZenithAngle);
  AzimuthAngle=(iAzimuthInterval+0.5)*dAzimuthAngle;

  cosZenithAngle=cos(ZenithAngle);
  sinZenithAngle=sin(ZenithAngle);

  cosAzimuthAngle=cos(AzimuthAngle);
  sinAzimuthAngle=sin(AzimuthAngle);

  Speed=Relativistic::E2Speed(Energy,PIC::MolecularData::GetMass(spec));

  for (idim=0;idim<3;idim++) v[idim]=Speed*(InternalFaceNorm[iface][idim]*sinZenithAngle+cosZenithAngle*(e0FaceFrame[iface][idim]*cosAzimuthAngle+e1FaceFrame[iface][idim]*sinAzimuthAngle));
}


//init the namespace
void Earth::CutoffRigidity::DomainBoundaryParticleProperty::Init() {
  int iface;

  for (iface=0;iface<6;iface++) {
    switch (iface) {
    case 0:case 1:
      dX[iface][0]=(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1])/SampleMaskNumberPerSpatialDirection;
      dX[iface][1]=(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2])/SampleMaskNumberPerSpatialDirection;
      break;
    case 2:case 3:
      dX[iface][0]=(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0])/SampleMaskNumberPerSpatialDirection;
      dX[iface][1]=(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2])/SampleMaskNumberPerSpatialDirection;
      break;
    case 4:case 5:
      dX[iface][0]=(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0])/SampleMaskNumberPerSpatialDirection;
      dX[iface][1]=(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1])/SampleMaskNumberPerSpatialDirection;
      break;
    }
  }
}

//register properties of the particles that cross the external boundary
void Earth::CutoffRigidity::DomainBoundaryParticleProperty::RegisterParticleProperties(int spec,double *x,double *v,int iOriginLocation,int iface) {
  int Index;

  if (EnableSampleParticleProperty==true) {
    //determine on which Samplking Table to refister the particle velocity
    int iTable,jTable;

    switch (iface) {
    case 0:case 1:
      iTable=(x[1]-PIC::Mesh::mesh.xGlobalMin[1])/dX[iface][1];
      jTable=(x[2]-PIC::Mesh::mesh.xGlobalMin[2])/dX[iface][2];
      break;
    case 2:case3:
      iTable=(x[0]-PIC::Mesh::mesh.xGlobalMin[0])/dX[iface][0];
      jTable=(x[2]-PIC::Mesh::mesh.xGlobalMin[2])/dX[iface][2];
      break;
    case 4:case5:
      iTable=(x[0]-PIC::Mesh::mesh.xGlobalMin[0])/dX[iface][0];
      jTable=(x[1]-PIC::Mesh::mesh.xGlobalMin[1])/dX[iface][1];
      break;
    }

    if (iTable>=SampleMaskNumberPerSpatialDirection) iTable=SampleMaskNumberPerSpatialDirection-1;
    if (jTable>=SampleMaskNumberPerSpatialDirection) jTable=SampleMaskNumberPerSpatialDirection-1;


    //register the particle velocity vector
    Index=GetVelocityVectorIndex(spec,v,iface);
    SampleTable[iOriginLocation][spec][iface][iTable][jTable].SetFlag(true,Index);
  }
}

//gather the flags from all processors
void Earth::CutoffRigidity::DomainBoundaryParticleProperty::Gather() {
  for (int iTestLocation=0;iTestLocation<std::max(1,Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength);iTestLocation++) for (int s=0;s<PIC::nTotalSpecies;s++) for (int iface=0;iface<6;iface++) {
    for (int i=0;i<SampleMaskNumberPerSpatialDirection;i++) for (int j=0;j<SampleMaskNumberPerSpatialDirection;j++) {
      SampleTable[iTestLocation][s][iface][i][j].Gather();
    }
  }
}

//test the particle properties
bool Earth::CutoffRigidity::DomainBoundaryParticleProperty::TestInjectedParticleProperties(int spec,double *x,double *v,int iTestLocation,int iface) {
  int Index,iTable,jTable;

  //determine on which Sampling Table to refister the particle velocity
  switch (iface) {
  case 0:case 1:
    iTable=(x[1]-PIC::Mesh::mesh.xGlobalMin[1])/dX[iface][1];
    jTable=(x[2]-PIC::Mesh::mesh.xGlobalMin[2])/dX[iface][2];
    break;
  case 2:case3:
    iTable=(x[0]-PIC::Mesh::mesh.xGlobalMin[0])/dX[iface][0];
    jTable=(x[2]-PIC::Mesh::mesh.xGlobalMin[2])/dX[iface][2];
    break;
  case 4:case5:
    iTable=(x[0]-PIC::Mesh::mesh.xGlobalMin[0])/dX[iface][0];
    jTable=(x[1]-PIC::Mesh::mesh.xGlobalMin[1])/dX[iface][1];
    break;
  }

  if (iTable>=SampleMaskNumberPerSpatialDirection) iTable=SampleMaskNumberPerSpatialDirection-1;
  if (jTable>=SampleMaskNumberPerSpatialDirection) jTable=SampleMaskNumberPerSpatialDirection-1;

  //Test Particle Properties
  Index=GetVelocityVectorIndex(spec,v,iface);
  return SampleTable[iTestLocation][spec][iface][iTable][jTable].Test(Index);
}

//Smooth the sampled distribution
void Earth::CutoffRigidity::DomainBoundaryParticleProperty::SmoothSampleTable() {
  int spec,i,j,k,iMin,iMax,jMin,jMax,Index,iTable,jTable,nSetPoints,IndexMax,iLogEnergyLevel,iface,iCosZenithInterval,iAzimuthInterval;

  IndexMax=nAzimuthIntervals*nCosZenithIntervals*nLogEnergyLevels;

  for (int iTestLocation=0;iTestLocation<std::max(1,Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength);iTestLocation++) for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iface=0;iface<6;iface++) {
    for (iTable=0;iTable<SampleMaskNumberPerSpatialDirection;iTable++) for (jTable=0;jTable<SampleMaskNumberPerSpatialDirection;jTable++) {
      //loop through all sampling tables

      for (Index=0,nSetPoints=0;Index<IndexMax;Index++) if (SampleTable[iTestLocation][spec][iface][iTable][jTable].Test(Index)==true) nSetPoints++;

      if (nSetPoints!=0) {
        //expand the set of the sampled points so a given fraction of the points is occupied
        bool ProcessFlagTable[nCosZenithIntervals][nAzimuthIntervals][nLogEnergyLevels];

        while (nSetPoints<SamplingParameters::OccupiedPhaseSpaceTableFraction*IndexMax) {
          nSetPoints=0;

          //reset the flag table
          for (iCosZenithInterval=0;iCosZenithInterval<nCosZenithIntervals;iCosZenithInterval++) for (iAzimuthInterval=0;iAzimuthInterval<nAzimuthIntervals;iAzimuthInterval++) {
            for (iLogEnergyLevel=0;iLogEnergyLevel<nLogEnergyLevels;iLogEnergyLevel++) {
              ProcessFlagTable[iCosZenithInterval][iAzimuthInterval][iLogEnergyLevel]=false;
            }
          }

          //scan through the sampled distribution
          for (iCosZenithInterval=0;iCosZenithInterval<nCosZenithIntervals;iCosZenithInterval++) for (iAzimuthInterval=0;iAzimuthInterval<nAzimuthIntervals;iAzimuthInterval++) {
            for (iLogEnergyLevel=0;iLogEnergyLevel<nLogEnergyLevels;iLogEnergyLevel++) {
              Index=iAzimuthInterval+(iCosZenithInterval+iLogEnergyLevel*nCosZenithIntervals)*nAzimuthIntervals;

              if ((SampleTable[iTestLocation][spec][iface][iTable][jTable].Test(Index)==true)&&(ProcessFlagTable[iCosZenithInterval][iAzimuthInterval][iLogEnergyLevel]==false)) {
                //check if the distribution souble be expanded
                ProcessFlagTable[iCosZenithInterval][iAzimuthInterval][iLogEnergyLevel]=true;
                nSetPoints++;

                //scan through the neibourhood of the found point
                for (i=-1;i<=1;i++) if ((iCosZenithInterval+i>=0)&&(iCosZenithInterval+i<nCosZenithIntervals)) {
                  for (j=-1;j<=1;j++)  if ((iAzimuthInterval+j>=0)&&(iAzimuthInterval+i<nAzimuthIntervals)) {
                    for (k=-1;k<=1;k++) if ((iLogEnergyLevel+k>=0)&&(iLogEnergyLevel+k<nLogEnergyLevels)) {
                      if (ProcessFlagTable[iCosZenithInterval+i][iAzimuthInterval+j][iLogEnergyLevel+k]==false) {
                        ProcessFlagTable[iCosZenithInterval+i][iAzimuthInterval+j][iLogEnergyLevel+j]=true;
                        Index=(iAzimuthInterval+j)+((iCosZenithInterval+i)+(iLogEnergyLevel+k)*nCosZenithIntervals)*nAzimuthIntervals;

                        if (SampleTable[iTestLocation][spec][iface][iTable][jTable].Test(Index)==false) {
                          SampleTable[iTestLocation][spec][iface][iTable][jTable].SetFlag(true,Index);
                        }

                        nSetPoints++;
                      }
                    }
                  }
                }
              }
            }
          }

        }
      }

    }
  }
}










