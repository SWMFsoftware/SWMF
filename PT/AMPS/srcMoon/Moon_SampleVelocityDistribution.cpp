/*
 * Moon_SampleVelocityDistribution.cpp
 *
 *  Created on: Mar 31, 2013
 *      Author: vtenishe
 */


//$Id$
//create a sample of the velocity distribution of the exospheric particles in anti-sunward direction


#include "pic.h"

Moon::Sampling::VelocityDistribution::cVelocitySampleBuffer *Moon::Sampling::VelocityDistribution::SampleBuffer=NULL;
int Moon::Sampling::VelocityDistribution::nTotalSampleDirections=0.0;
int Moon::Sampling::VelocityDistribution::nZenithPoints=0;


void Moon::Sampling::VelocityDistribution::Init() {
  int nSampleDirections=0;

  //calcualte the number of the sampling directions
  double dZ,rr,ZenithAngle,AzimuthAngle,dZenithAngle;

  dZ=dZenithAngleMin;
  rr=(maxZenithAngle+dZenithAngleMax)/(maxZenithAngle+dZenithAngleMin);
  nZenithPoints=(long int)(log(dZenithAngleMax/dZenithAngleMin)/log(rr)-2.0);
  rr=pow(dZenithAngleMax/dZenithAngleMin,1.0/(nZenithPoints+2.0));

  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {
    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;
  }


  nSampleDirections=nZenithPoints*(1+nAzimuthPoints);
  SampleBuffer=new cVelocitySampleBuffer[nSampleDirections];

  //calculate the pointing direction of the sampling
  dZ=dZenithAngleMin;
  rr=(maxZenithAngle+dZenithAngleMax)/(maxZenithAngle+dZenithAngleMin);
  nZenithPoints=(long int)(log(dZenithAngleMax/dZenithAngleMin)/log(rr)-2.0);
  rr=pow(dZenithAngleMax/dZenithAngleMin,1.0/(nZenithPoints+2.0));

  nZenithPoints=0,ZenithAngle=dZenithAngleMin,dZenithAngle=dZenithAngleMin;

  while (ZenithAngle<maxZenithAngle) {
    ZenithAngle+=dZenithAngle;
    dZenithAngle*=rr;
    nZenithPoints++;

    for (int iAzimuthPoint=0;iAzimuthPoint<1+nAzimuthPoints;iAzimuthPoint++) {
      AzimuthAngle=2.0*Pi*double(iAzimuthPoint)/double(nAzimuthPoints);

      SampleBuffer[nTotalSampleDirections].lGSE[0]=-cos(ZenithAngle);
      SampleBuffer[nTotalSampleDirections].lGSE[1]=sin(ZenithAngle)*sin(AzimuthAngle);
      SampleBuffer[nTotalSampleDirections].lGSE[2]=sin(ZenithAngle)*cos(AzimuthAngle);

      nTotalSampleDirections++;
    }
  }
}


void Moon::Sampling::VelocityDistribution::Sampling() {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  int nSampleDirection;
  int idim;
  SpiceDouble xform[6][6],EarthState[6],lt,SunState[6];
  double xLSO[3],dl,vEarthLSO[3],vSunLSO[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xNode;
  SpiceDouble lLSO[6]={0,0,0,0,0,0};  //only 3 first components of the vectors are used. all 6 components are needed in the definition in order SPICE routines work correctly

  //the ratio between the step of the integration procedure and the local cell size
  static const double IntegrationStep2CellSizeRatio=0.3;

  //calculate the currect position of the Earth and the Sun
  spkezr_c("Earth",Exosphere::OrbitalMotion::et,"LSO","none","Moon",EarthState,&lt);
  spkezr_c("Sun",Exosphere::OrbitalMotion::et,"LSO","none","Moon",SunState,&lt);

  for (idim=0;idim<3;idim++) vEarthLSO[idim]=1.0E3*EarthState[3+idim],vSunLSO[idim]=1.0E3*SunState[3+idim];

  //calculate the rotation matrix from 'GSE' to 'LSO'
  sxform_c("GSE","LSO",Exosphere::OrbitalMotion::et,xform);


  for (nSampleDirection=0;nSampleDirection<nTotalSampleDirections;nSampleDirection++) {
    //temporary sampling buffer and the buffer for the particle data
    cVelocitySampleBuffer tempSamplingBuffer;
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

    memcpy(&tempSamplingBuffer,SampleBuffer+nSampleDirection,sizeof(cVelocitySampleBuffer));

    //convert the pointing vector from 'GSE' to 'LSO'
    mxvg_c(xform,tempSamplingBuffer.lGSE,6,6,lLSO);
    for (idim=0;idim<3;idim++) xLSO[idim]=1.0E3*EarthState[idim];

    //determine the limits of integration
    double IntegrationPathLength,xStart[3],xFinish[3];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *xFinishNode;

    PIC::ColumnIntegration::FindIntegrationLimits(xLSO,lLSO,IntegrationPathLength,xStart,xNode,xFinish,xFinishNode);
    memcpy(xLSO,xStart,3*sizeof(double));

    while (xNode!=NULL) {
      dl=IntegrationStep2CellSizeRatio*xNode->GetCharacteristicCellSize();

      if (xNode->Thread==PIC::ThisThread) {
        //find the cell for sampling
        int i,j,k;
        long int ncell,ptr;
        PIC::ParticleBuffer::byte *ParticleData;
        int spec;
        double *v,*x,LocalParticleWeight,VelocityLineOfSight,HeliocentricRadialVelocity,c,cellMeasure,rHeliocentric,Speed;

        ncell=PIC::Mesh::mesh.fingCellIndex(xLSO,i,j,k,xNode);
        cellMeasure=xNode->block->GetCenterNode(ncell)->Measure;
//        dl=IntegrationStep2CellSizeRatio*xNode->GetCharacteristicCellSize();

        //sample particles
        ptr=xNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        while (ptr!=-1) {
          ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          memcpy((void*)tempParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);

          spec=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)tempParticleData);
          v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)tempParticleData);
          x=PIC::ParticleBuffer::GetX((PIC::ParticleBuffer::byte*)tempParticleData);

          LocalParticleWeight=xNode->block->GetLocalParticleWeight(spec)/cellMeasure;
          LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)tempParticleData);

          //calculate velocity component in the direction of the line of sight (in respent to the Earth)
          for (idim=0,VelocityLineOfSight=0.0;idim<DIM;idim++) VelocityLineOfSight+=(v[idim]-vEarthLSO[idim])*lLSO[idim];
          i=(int)((VelocityLineOfSight+maxVelocityLimit)/VelocityBinWidth);
          if ((i>=0.0)&&(i<nVelocitySamplePoints)) tempSamplingBuffer.VelocityLineOfSight[spec][i]+=LocalParticleWeight*dl;

          //calcualte the radial component of the heliocentric velocity
          for (idim=0,c=0.0;idim<DIM;idim++) c+=pow(x[idim]-Moon::xSun_SO[idim],2);

          for (Speed=0.0,rHeliocentric=sqrt(c),HeliocentricRadialVelocity=0.0,idim=0;idim<DIM;idim++) {
            HeliocentricRadialVelocity+=(v[idim]-vSunLSO[idim])*(x[idim]-Moon::xSun_SO[idim])/rHeliocentric;
            Speed+=pow(v[idim]-vSunLSO[idim],2);
          }

          Speed=sqrt(Speed);

          i=(int)((HeliocentricRadialVelocity+maxVelocityLimit)/VelocityBinWidth);
          if ((i>=0.0)&&(i<nVelocitySamplePoints)) tempSamplingBuffer.VelocityRadialHeliocentric[spec][i]+=LocalParticleWeight*dl;

          //sample particle speed (in respect to the Sun)
          i=(int)((Speed+maxVelocityLimit)/VelocityBinWidth);
          if ((i>=0.0)&&(i<nVelocitySamplePoints)) tempSamplingBuffer.Speed[spec][i]+=LocalParticleWeight*dl;

          //sample the mean value of the line-og-sight and heliocentric-radial velocities
          tempSamplingBuffer.meanVelocityLineOfSight[spec]+=VelocityLineOfSight*LocalParticleWeight*dl;
          tempSamplingBuffer.meanVelocityRadialHeliocentric[spec]+=HeliocentricRadialVelocity*LocalParticleWeight*dl;
          tempSamplingBuffer.meanSpeed[spec]+=Speed*LocalParticleWeight*dl;
          tempSamplingBuffer.ColumnDensityIntegral[spec]+=LocalParticleWeight*dl;

          //sample the exospehric brightness
          if (spec==_NA_SPEC_) if ( (Moon::EarthShadowCheck(x)==false) && ((x[0]>0.0)||(x[1]*x[1]+x[2]*x[2]>_RADIUS_(_MOON_)*_RADIUS_(_MOON_))) ) {
            tempSamplingBuffer.Brightness[spec]+=1.0E-10*
                (SodiumGfactor__5891_58A__Killen_2009_AJSS(HeliocentricRadialVelocity,rHeliocentric)+SodiumGfactor__5897_56A__Killen_2009_AJSS(HeliocentricRadialVelocity,rHeliocentric))*LocalParticleWeight*dl;
          }


          ptr=PIC::ParticleBuffer::GetNext((PIC::ParticleBuffer::byte*)tempParticleData);
        }
      }


      for (idim=0;idim<3;idim++) xLSO[idim]+=dl*lLSO[idim];

      //break the loop if the new point is within the Moon or outside of the computational domain
      if (xLSO[0]*xLSO[0]+xLSO[0]*xLSO[0]+xLSO[0]*xLSO[0]<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_)) break;

      xNode=PIC::Mesh::mesh.findTreeNode(xLSO,xNode);
    }


    memcpy(SampleBuffer+nSampleDirection,&tempSamplingBuffer,sizeof(cVelocitySampleBuffer));
  }
#endif
}


void Moon::Sampling::VelocityDistribution::OutputSampledData(int DataOutputFileNumber) {
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  char fname[_MAX_STRING_LENGTH_PIC_];
  int nZone,s,i,cnt;
  FILE *fout=NULL,*fMAP=NULL;

  const int mpiZoneExchangeBufferLength=PIC::nTotalSpecies*(5+3*nVelocitySamplePoints);
  double mpiZoneExchangeBuffer[mpiZoneExchangeBufferLength];


  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/pic.Moon.VelocityDistribution.AntisolarDirection.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"v\"");

    for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,", \"f Velocity along the Line of Sight (%s) (in respect to the Earth)\", \"f Heliocentric Radial Component Velocity (%s) (in respect to the Sun)\", \"f Speed (%s) (in respect to the Sun)\" ",
        PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s));
    fprintf(fout,"\n");


    //create a map of the zones
    sprintf(fname,"%s/pic.Moon.VelocityDistribution.AntisolarDirection.Map.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fMAP=fopen(fname,"w");

    fprintf(fMAP,"VARIABLES=\"l0GSE\", \"l1GSE\", \"RA (hourse)\", \"DEC (deg)\", \"nZone\"");
    for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fMAP,", \"Mean Velocity along the Line of Sight (%s) (in respect to the Earth)\", \"Mean Heliocentric Radial Component Velocity (%s) (in respect to the Sun)\", \"Mean Speed along the Line of Sight (%s) (in respect to the Sun)\", \"Column DensityIntegral (%s)\", \"Brightness along the line of sight (%s)\"",
        PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s),PIC::MolecularData::GetChemSymbol(s));
    fprintf(fMAP,"\n");

    fprintf(fMAP,"ZONE T=\"Velocity Distribution Zone Map\"\n");
    fprintf(fMAP,"I=%i, J=%i, K=1, ZONETYPE=Ordered\n",nAzimuthPoints+1,nZenithPoints);
    fprintf(fMAP,"DATAPACKING=POINT\n");
  }



  for (nZone=0;nZone<nTotalSampleDirections;nZone++) {
    SpiceDouble lJ2000[6],xform[6][6],r,ra,dec;

    sxform_c("GSE","J2000",Exosphere::OrbitalMotion::et,xform);
    mxvg_c(xform,SampleBuffer[nZone].lGSE,6,6,lJ2000);
    recrad_c (lJ2000,&r,&ra,&dec);
    if (ra>Pi) ra-=2.0*Pi;

    if (PIC::ThisThread==0) {
      fprintf(fout,"ZONE T=\"RA=%e (hours), DEC=%e (deg)\",I=%i\n",ra/Pi*12.0,dec/Pi*180.0,nVelocitySamplePoints);
    }

    //combine all sampled data on the root processor and print them into a file
    for (cnt=0,s=0;s<PIC::nTotalSpecies;s++) {
      mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].meanVelocityLineOfSight[s];
      mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].meanVelocityRadialHeliocentric[s];
      mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].meanSpeed[s];
      mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].ColumnDensityIntegral[s];
      mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].Brightness[s];

      for (i=0;i<nVelocitySamplePoints;i++) {
        mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].VelocityLineOfSight[s][i];
        mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].VelocityRadialHeliocentric[s][i];
        mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].Speed[s][i];
      }
    }


    if (PIC::ThisThread==0) {
      int thread;

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        MPI_Status status;

        MPI_Recv(mpiZoneExchangeBuffer,mpiZoneExchangeBufferLength,MPI_DOUBLE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (cnt=0,s=0;s<PIC::nTotalSpecies;s++) {
          SampleBuffer[nZone].meanVelocityLineOfSight[s]+=mpiZoneExchangeBuffer[cnt++];
          SampleBuffer[nZone].meanVelocityRadialHeliocentric[s]+=mpiZoneExchangeBuffer[cnt++];
          SampleBuffer[nZone].meanSpeed[s]+=mpiZoneExchangeBuffer[cnt++];
          SampleBuffer[nZone].ColumnDensityIntegral[s]+=mpiZoneExchangeBuffer[cnt++];
          SampleBuffer[nZone].Brightness[s]+=mpiZoneExchangeBuffer[cnt++];

          for (i=0;i<nVelocitySamplePoints;i++) {
            SampleBuffer[nZone].VelocityLineOfSight[s][i]+=mpiZoneExchangeBuffer[cnt++];
            SampleBuffer[nZone].VelocityRadialHeliocentric[s][i]+=mpiZoneExchangeBuffer[cnt++];
            SampleBuffer[nZone].Speed[s][i]+=mpiZoneExchangeBuffer[cnt++];
          }
        }
      }

      //normalize the distribution function
      for (s=0;s<PIC::nTotalSpecies;s++) {
        double w;

        for (w=0.0,i=0;i<nVelocitySamplePoints;i++) w+=SampleBuffer[nZone].VelocityLineOfSight[s][i];
        if (w>0.0) for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].VelocityLineOfSight[s][i]/=w*VelocityBinWidth;

        for (w=0.0,i=0;i<nVelocitySamplePoints;i++) w+=SampleBuffer[nZone].VelocityRadialHeliocentric[s][i];
        if (w>0.0) for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].VelocityRadialHeliocentric[s][i]/=w*VelocityBinWidth;

        for (w=0.0,i=0;i<nVelocitySamplePoints;i++) w+=SampleBuffer[nZone].Speed[s][i];
        if (w>0.0) for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].Speed[s][i]/=w*VelocityBinWidth;

        if (SampleBuffer[nZone].ColumnDensityIntegral[s]>0.0) {
          SampleBuffer[nZone].meanVelocityLineOfSight[s]/=SampleBuffer[nZone].ColumnDensityIntegral[s];
          SampleBuffer[nZone].meanVelocityRadialHeliocentric[s]/=SampleBuffer[nZone].ColumnDensityIntegral[s];
          SampleBuffer[nZone].meanSpeed[s]/=SampleBuffer[nZone].ColumnDensityIntegral[s];

          //devide the calcualted column density by the number of the iterations in the sample (MUST be the last operation within the loop)
          SampleBuffer[nZone].ColumnDensityIntegral[s]/=PIC::LastSampleLength;
          SampleBuffer[nZone].Brightness[s]/=PIC::LastSampleLength;
        }
      }


      //output sample data into a file
      for (i=0;i<nVelocitySamplePoints;i++) {
        fprintf(fout,"%e  ",-maxVelocityLimit+VelocityBinWidth*(i+0.5));

        for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,"  %e   %e  %e", SampleBuffer[nZone].VelocityLineOfSight[s][i],SampleBuffer[nZone].VelocityRadialHeliocentric[s][i],SampleBuffer[nZone].Speed[s][i]);
        fprintf(fout,"\n");
      }

      //print the map and the mean values
      fprintf(fMAP,"%e %e %e %e %i",SampleBuffer[nZone].lGSE[1],SampleBuffer[nZone].lGSE[2],ra/Pi*12.0,dec/Pi*180.0,nZone);
      for (s=0;s<PIC::nTotalSpecies;s++) {
        fprintf(fMAP,"  %e  %e  %e  %e  %e",SampleBuffer[nZone].meanVelocityLineOfSight[s],SampleBuffer[nZone].meanVelocityRadialHeliocentric[s],SampleBuffer[nZone].meanSpeed[s],SampleBuffer[nZone].ColumnDensityIntegral[s],SampleBuffer[nZone].Brightness[s]);
      }

      fprintf(fMAP,"\n");
    }
    else MPI_Send(mpiZoneExchangeBuffer,mpiZoneExchangeBufferLength,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);

    //clear the sampling buffer
    for (s=0;s<PIC::nTotalSpecies;s++) {
      SampleBuffer[nZone].meanVelocityLineOfSight[s]=0.0,SampleBuffer[nZone].meanVelocityRadialHeliocentric[s]=0.0,SampleBuffer[nZone].ColumnDensityIntegral[s]=0.0,SampleBuffer[nZone].meanSpeed[s]=0.0,SampleBuffer[nZone].Brightness[s]=0.0;

      for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].VelocityLineOfSight[s][i]=0.0,SampleBuffer[nZone].VelocityRadialHeliocentric[s][i]=0.0,SampleBuffer[nZone].Speed[s][i]=0.0;
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }


  //close the file
  if (PIC::Mesh::mesh.ThisThread==0) {
    fclose(fout);
    fclose(fMAP);
  }
#endif
}


