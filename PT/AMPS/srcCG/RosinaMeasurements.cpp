/*
 * RosinaMeasurements.cpp
 *
 *  Created on: Dec 3, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "RosinaMeasurements.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif


RosinaSample::cRosinaSamplingLocation RosinaSample::Rosina[RosinaSample::nPoints];

//the sampling time interval limit
SpiceDouble RosinaSample::SamplingTimeInterval::StartSamplingEt=0.0;
SpiceDouble RosinaSample::SamplingTimeInterval::EndSamplingEt=0.0;
char RosinaSample::SamplingTimeInterval::StartSamplingTime[_MAX_STRING_LENGTH_PIC_]="";
char RosinaSample::SamplingTimeInterval::EndSamplingTime[_MAX_STRING_LENGTH_PIC_]="";

char RosinaSample::SamplingTimeInterval::CurrentTimeWindowStamp[_MAX_STRING_LENGTH_PIC_]="";


//init the pointing directions and geomentry informationvoid



void RosinaSample::Init() {
  //init the start and end time of the particle sample
  //init the sampling point
  utc2et_c(SamplingTimeInterval::StartSamplingTime,&SamplingTimeInterval::StartSamplingEt);
  utc2et_c(SamplingTimeInterval::EndSamplingTime,&SamplingTimeInterval::EndSamplingEt);

  RosinaSample::Init(SamplingTimeInterval::StartSamplingEt,SamplingTimeInterval::EndSamplingEt);
}

void RosinaSample::Init(double etMin,double etMax) {
#ifndef _NO_SPICE_CALLS_
  int i,idim,iCell,jCell,kCell,nd;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  double dxSubCell,dySubCell,dzSubCell,xLocalCell,yLocalCell,zLocalCell;

  //the function call counter
  static int CallCounter=0;

  CallCounter++;

  //init line-of-sight vectors
  SpiceDouble lt,et,xRosetta[3],etStart;
  SpiceDouble       xform[6][6];

  FILE *fTrajectory=NULL;
  FILE *fTrajectoryZones=NULL;

  if (PIC::ThisThread==0) {
   char fname[1000];
   char TimeStampStart[100],TimeStampEnd[100];

   et2utc_c(etMin,"ISOC",4,100,TimeStampStart);
   et2utc_c(etMax,"ISOC",4,100,TimeStampEnd);
   printf("Setting sampling flags: data will be sampled for the time interval of %s - %s\n",TimeStampStart,TimeStampEnd);

   sprintf(fname,"%s/TrajectoryLastMeasurements.dat",PIC::OutputDataFileDirectory);
   fTrajectory=fopen(fname,"w");
   fprintf(fTrajectory,"VARIABLES=\"X\", \"Y\", \"Z\", \"iPoint\" \n");

   sprintf(fname,"%s/TrajectoryLastMeasurementsZones.dat",PIC::OutputDataFileDirectory);
   fTrajectoryZones=fopen(fname,"w");
   fprintf(fTrajectoryZones,"VARIABLES=\"X\", \"Y\", \"Z\" \n");
  }

  for (i=0;i<nPoints;i++) {
    utc2et_c(ObservationTime[i],&et);
    spkpos_c("ROSETTA",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xRosetta,&lt);
    for (idim=0;idim<3;idim++) xRosetta[idim]*=1.0E3;

    if (i==0) etStart=et;

    node=PIC::Mesh::mesh.findTreeNode(xRosetta,NULL);

    if (node!=NULL) {
      nd=PIC::Mesh::mesh.fingCellIndex(xRosetta,iCell,jCell,kCell,node);
    }
    else nd=-1,iCell=-1,jCell=-1,kCell-1;

    for (idim=0;idim<3;idim++) Rosina[i].x[idim]=xRosetta[idim];
    Rosina[i].nd=nd;
    Rosina[i].node=node;
    Rosina[i].iCell=iCell;
    Rosina[i].jCell=jCell;
    Rosina[i].kCell=kCell;

    Rosina[i].RadiusVectorLeangth=sqrt(xRosetta[0]*xRosetta[0]+xRosetta[1]*xRosetta[1]+xRosetta[2]*xRosetta[2]);
    Rosina[i].CometDistance=-1.0;
    Rosina[i].SecondsFromBegining=et-etStart;
    Rosina[i].CharacteristicCellSize=-1.0;

    if (Rosina[i].node==NULL) continue;

    dxSubCell=(Rosina[i].node->xmax[0]-Rosina[i].node->xmin[0])/_BLOCK_CELLS_X_/CellFractionationFactor;
    dySubCell=(Rosina[i].node->xmax[1]-Rosina[i].node->xmin[1])/_BLOCK_CELLS_Y_/CellFractionationFactor;
    dzSubCell=(Rosina[i].node->xmax[2]-Rosina[i].node->xmin[2])/_BLOCK_CELLS_Z_/CellFractionationFactor;

    Rosina[i].iLocalSubCell=(int)((xRosetta[0]-Rosina[i].node->xmin[0])/dxSubCell);
    Rosina[i].jLocalSubCell=(int)((xRosetta[1]-Rosina[i].node->xmin[1])/dySubCell);
    Rosina[i].kLocalSubCell=(int)((xRosetta[2]-Rosina[i].node->xmin[2])/dzSubCell);

    if (node==NULL) {
      Rosina[i].LocationCode=-1; //outside of the domain
    }
    else {
      if (CutCell::CheckPointInsideDomain(xRosetta,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0)==false) {
        Rosina[i].LocationCode=-2; //the point is inside the body
      }
      else {
        Rosina[i].LocationCode=0; //the point is within the domain
        Rosina[i].CharacteristicCellSize=sqrt(pow(node->xmax[0]-node->xmin[0],2)+pow(node->xmax[1]-node->xmin[1],2)+pow(node->xmax[2]-node->xmin[2],2));

        //determine the closest distance to the comet


        for (int ii=0;ii<CutCell::nBoundaryTriangleFaces;ii++) {
          double r,xCenter[3];
          int idim;

          CutCell::BoundaryTriangleFaces[ii].GetCenterPosition(xCenter);
          for (idim=0,r=0.0;idim<3;idim++) r+=pow(xCenter[idim]-xRosetta[idim],2);

          r=sqrt(r);
          if ((Rosina[i].CometDistance<0.0)||(r<Rosina[i].CometDistance)) Rosina[i].CometDistance=r;
        }
      }
    }

    if (PIC::ThisThread==0) {
      fprintf(fTrajectory,"%e %e %e %i\n",xRosetta[0],xRosetta[1],xRosetta[2],i);

      fprintf(fTrajectoryZones,"ZONE T=\"Point=%i\" F=POINT\n",i);
      fprintf(fTrajectoryZones,"%e  %e  %e\n",xRosetta[0],xRosetta[1],xRosetta[2]);
    }

    sxform_c("ROS_SPACECRAFT","67P/C-G_CK",et,xform);

    //pointing of COPS nude gauge: negative y
    SpiceDouble NudeGauge[6]={0.0,-1.0,0.0 ,0.0,0.0,0.0};
    SpiceDouble l[6];

    mxvg_c(xform,NudeGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[i].NudeGauge.LineOfSight[idim]=l[idim];


    //pointing of COPS ram gauge: positive z
    SpiceDouble RamGauge[6]={0.0,0.0,1.0 ,0.0,0.0,0.0};

    mxvg_c(xform,RamGauge,6,6,l);
    for (idim=0;idim<3;idim++) Rosina[i].RamGauge.LineOfSight[idim]=l[idim];


    if (CallCounter==1) {
      Rosina[i].NudeGauge.NucleusSolidAngle=0.0;
      Rosina[i].RamGauge.NucleusSolidAngle=0.0;

      //solid angle occupied by the nucleus by the nude gauge
      int t,nTotalTests=100000,iTest,iIntersectionFace,NucleusIntersectionCounter=0;
      double l[3],xIntersection[3];

      int iStartTest,iFinishTest,nTestThread;

      iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
      iFinishTest=(nTotalTests/PIC::nTotalThreads)*(PIC::ThisThread+1);

      nTestThread=nTotalTests/PIC::nTotalThreads;
      iStartTest=(nTotalTests/PIC::nTotalThreads)*PIC::ThisThread;
      iFinishTest=iStartTest+nTestThread;
      if (PIC::ThisThread==PIC::nTotalThreads-1) iFinishTest=nTotalTests;

      //get the solid angles, and estimate distance to the nucleus
      Liouville::GetSolidAngle(Rosina[i].NudeGauge.NucleusSolidAngle,Rosina[i].RamGauge.NucleusSolidAngle,i);
      Rosina[i].Altitude=PIC::Mesh::IrregularSurface::GetClosestDistance(Rosina[i].x,Rosina[i].xNucleusClosestPoint,Rosina[i].iNucleusClosestFace);
    }

    //init the sampling point
    SpiceDouble et;

    utc2et_c(ObservationTime[i],&et);
    Rosina[i].SamplingParticleDataFlag=((etMin<=et)&&(et<etMax)) ? true : false;
  }

  if (PIC::ThisThread==0) {
    fclose(fTrajectory);
    fclose(fTrajectoryZones);
    printf("$PREFIX: RosinaSample::Init is complete.....\n");
  }
#endif  
  
  Flush();
}

//flush the data buffers
void RosinaSample::Flush() {
  int spec,i;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) for (i=0;i<nPoints;i++) {
    Rosina[i].NudeGauge.Density[spec]=0.0,Rosina[i].RamGauge.Density[spec]=0.0;
    Rosina[i].NudeGauge.Flux[spec]=0.0,Rosina[i].RamGauge.Flux[spec]=0.0;
  }
}

//the sampling routine
void RosinaSample::SamplingProcessor() {
  int i,j,k,spec,idim;
  long int ptr;
  double ParticleWeight,*v,c,*x;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  int iCellFraction,jCellFraction,kCellFraction;
  double dxSubCell,dySubCell,dzSubCell;
  int iLocalSubCell,jLocalSubCell,kLocalSubCell;

  #ifdef _NO_SPICE_CALLS_
  return; //cannot sample the instrument related data when SPICE is not used
  #endif

  for (i=0;i<nPoints;i++) if ((Rosina[i].node!=NULL)&&(Rosina[i].SamplingParticleDataFlag==true)) if ((block=Rosina[i].node->block)!=NULL) {
    ptr=block->FirstCellParticleTable[Rosina[i].iCell+_BLOCK_CELLS_X_*(Rosina[i].jCell+_BLOCK_CELLS_Y_*Rosina[i].kCell)];

    dxSubCell=(Rosina[i].node->xmax[0]-Rosina[i].node->xmin[0])/_BLOCK_CELLS_X_/CellFractionationFactor;
    dySubCell=(Rosina[i].node->xmax[1]-Rosina[i].node->xmin[1])/_BLOCK_CELLS_Y_/CellFractionationFactor;
    dzSubCell=(Rosina[i].node->xmax[2]-Rosina[i].node->xmin[2])/_BLOCK_CELLS_Z_/CellFractionationFactor;

    while (ptr!=-1) {
      spec=PIC::ParticleBuffer::GetI(ptr);
      v=PIC::ParticleBuffer::GetV(ptr);
      x=PIC::ParticleBuffer::GetX(ptr);
      iLocalSubCell=(int)((x[0]-Rosina[i].node->xmin[0])/dxSubCell);
      jLocalSubCell=(int)((x[1]-Rosina[i].node->xmin[1])/dySubCell);
      kLocalSubCell=(int)((x[2]-Rosina[i].node->xmin[2])/dzSubCell);

      if ((Rosina[i].iLocalSubCell!=iLocalSubCell)||(Rosina[i].jLocalSubCell!=jLocalSubCell)||(Rosina[i].kLocalSubCell!=kLocalSubCell)) {
        ptr=PIC::ParticleBuffer::GetNext(ptr);
        continue;
      }

      ParticleWeight=block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
      if ((cell=block->GetCenterNode(Rosina[i].nd))!=NULL) ParticleWeight/=(cell->Measure/pow(CellFractionationFactor,3));

      //can the particle be registered by the ram gauge?
      if ((c=Rosina[i].RamGauge.LineOfSight[0]*v[0]+Rosina[i].RamGauge.LineOfSight[1]*v[1]+Rosina[i].RamGauge.LineOfSight[2]*v[2])<0.0) {
        //the particle can be registered by the instrument: sample the flux, and density dust to this particles
        Rosina[i].RamGauge.Density[spec]+=ParticleWeight;
        Rosina[i].RamGauge.Flux[spec]-=ParticleWeight*c;
      }

      //can the particle be registered by the nude gauge?
      if ((c=Rosina[i].NudeGauge.LineOfSight[0]*v[0]+Rosina[i].NudeGauge.LineOfSight[1]*v[1]+Rosina[i].NudeGauge.LineOfSight[2]*v[2])<0.0) {
        //the particle can be registered by the instrument: sample the flux, and density dust to this particles
        Rosina[i].NudeGauge.Density[spec]+=ParticleWeight;
        Rosina[i].NudeGauge.Flux[spec]-=ParticleWeight*c;
      }

      ptr=PIC::ParticleBuffer::GetNext(ptr);
    }
  }

  //remove all particles that would not be able to reach the spacecraft
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  int iPoint,ptrNext,iLocalNode;
  bool DirectAccessPathFound;

  for (iLocalNode=0;iLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;iLocalNode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[iLocalNode];

    if (node->block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          if (ptr!=-1) {
            v=PIC::ParticleBuffer::GetV(ptr);
            x=PIC::ParticleBuffer::GetX(ptr);
            DirectAccessPathFound=false;

            //check whether the particle can reach the nude gauge:
            //1 check whether the particle trajectory intersects a sphere with the radius of 2 cells, and centered in the observations point
            //2 the direction of the particle velocity is opposite to the line of sight ofthe instrument
            for (iPoint=0;iPoint<nPoints;iPoint++) {
              double a=0.0,b=0.0,c=0.0,d2,t;

              for (idim=0;idim<3;idim++) {
                a+=pow(v[idim],2);
                b+=2.0*(x[idim]-Rosina[iPoint].x[idim])*v[idim];
                c+=pow(x[idim]-Rosina[iPoint].x[idim],2);
              }

              c-=pow(3.0*node->GetCharacteristicCellSize(),2);
              d2=b*b-4.0*a*c;

              if (d2>0.0) {
                double d=sqrt(d2);

                if ( ((-b+d)/(2.0*a)>0.0) || ((-b-d)/(2.0*a)>0.0) ) {
                  if ( (Vector3D::DotProduct(Rosina[iPoint].RamGauge.LineOfSight,v)<0.0) || (Vector3D::DotProduct(Rosina[iPoint].NudeGauge.LineOfSight,v)<0.0) ) {
                    DirectAccessPathFound=true;
                    break;
                  }
                }
              }
            }

/*            //a particle can reach the gauge if
            //1) it is in the direction of the line of sight, and
            //2) the particle vbelocity is opposite to the line of sight
            for (iPoint=0;iPoint<nPoints;iPoint++) {
              for (c=0.0,idim=0;idim<3;idim++) c+=Rosina[iPoint].NudeGauge.LineOfSight[idim]*(x[idim]-Rosina[iPoint].x[idim]);

              if ((c>=0)&&(Vector3D::DotProduct(Rosina[iPoint].NudeGauge.LineOfSight,v)<=0.0)) {
                DirectAccessPathFound=true;
                break;
              }
            }

            //check whether the particle can reach the ram gauge
            if (DirectAccessPathFound==false) for (iPoint=0;iPoint<nPoints;iPoint++) {
              for (c=0.0,idim=0;idim<3;idim++) c+=Rosina[iPoint].RamGauge.LineOfSight[idim]*(x[idim]-Rosina[iPoint].x[idim]);

              if ((c>=0)&&(Vector3D::DotProduct(Rosina[iPoint].RamGauge.LineOfSight,v)<=0.0)) {
                DirectAccessPathFound=true;
                break;
              }
            }*/

            ptrNext=PIC::ParticleBuffer::GetNext(ptr);

            if (DirectAccessPathFound==false) {
              //remove the particle
              long int t=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              PIC::ParticleBuffer::DeleteParticle(ptr,t);
              node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=t;
            }

            ptr=ptrNext;
          }
        }
      }
    }
  }
}


//output sampled flux and density
void RosinaSample::PrintOutputFile(int nfile) {

  //collect all sampled data on the root processor
  double ExchangeBuffer[PIC::nTotalSpecies*nPoints];
  double DensityRamGauge[PIC::nTotalSpecies*nPoints],FluxRamGauge[PIC::nTotalSpecies*nPoints];
  double DensityNudeGauge[PIC::nTotalSpecies*nPoints],FluxNudeGauge[PIC::nTotalSpecies*nPoints];

  int i,spec;

  //collect all sampled data
  for (i=0;i<nPoints;i++) for (spec=0;spec<PIC::nTotalSpecies;spec++) ExchangeBuffer[spec+i*PIC::nTotalSpecies]=Rosina[i].RamGauge.Density[spec]/PIC::LastSampleLength;
  MPI_Reduce(ExchangeBuffer,DensityRamGauge,PIC::nTotalSpecies*nPoints,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  for (i=0;i<nPoints;i++) for (spec=0;spec<PIC::nTotalSpecies;spec++) ExchangeBuffer[spec+i*PIC::nTotalSpecies]=Rosina[i].RamGauge.Flux[spec]/PIC::LastSampleLength;
  MPI_Reduce(ExchangeBuffer,FluxRamGauge,PIC::nTotalSpecies*nPoints,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  for (i=0;i<nPoints;i++) for (spec=0;spec<PIC::nTotalSpecies;spec++) ExchangeBuffer[spec+i*PIC::nTotalSpecies]=Rosina[i].NudeGauge.Density[spec]/PIC::LastSampleLength;
  MPI_Reduce(ExchangeBuffer,DensityNudeGauge,PIC::nTotalSpecies*nPoints,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  for (i=0;i<nPoints;i++) for (spec=0;spec<PIC::nTotalSpecies;spec++) ExchangeBuffer[spec+i*PIC::nTotalSpecies]=Rosina[i].NudeGauge.Flux[spec]/PIC::LastSampleLength;
  MPI_Reduce(ExchangeBuffer,FluxNudeGauge,PIC::nTotalSpecies*nPoints,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

  //output the data file
  if (PIC::ThisThread==0) for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    char fname[1000];
    FILE *fout=NULL;

    //sampeled density
    sprintf(fname,"%s/RamNudeGaugeModelDensity.time=%s.spec=%s.out=%i.dat",PIC::OutputDataFileDirectory,RosinaSample::SamplingTimeInterval::CurrentTimeWindowStamp,PIC::MolecularData::GetChemSymbol(spec),nfile);
    fout=fopen(fname,"w");
    if (fout==0) exit(__LINE__,__FILE__,"Error: cannot open file");


    double SecondsFromBegining,RadiusVectorLeangth,CometDistance,beta;
    int LocationCode;


    fprintf(fout,"VARIABLES=\"i\", \"Ram Gauge Density\", \"Nude Gauge Density\", \"Ram Gauge Pressure\", \"Nude Gauge Pressure\", \"Total Ram Gauge Pressure\", \"Total Nude Gauge Pressure\", \"Seconds From The First Point\", \" Radius-Vector Length\", \"Distance to the Coment\", \"Location Code\", \"Characteristic Cell Size\", \"Ram Guage Nucleus Solid angle\", \"Nude Gauge Nucleus Solid Angle\", \"Altitude\" \n");

    for (int i=0;i<nPoints;i++) {
      double rgPressure=0.0;
      double ngPressure=0.0;

      double rgTotalPressure=0.0;
      double ngTotalPressure=0.0;

      const double rgTemperature=293.0;
      const double ngTemperature=293.0;

      beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*rgTemperature));
      rgPressure=Kbol*rgTemperature*Pi*beta/2.0*FluxRamGauge[spec+i*PIC::nTotalSpecies];


 //     rgPressure=4.0*Kbol*rgTemperature*FluxRamGauge[spec+i*PIC::nTotalSpecies]/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(spec)));


      ngPressure=DensityNudeGauge[spec+i*PIC::nTotalSpecies]*Kbol*ngTemperature;

      if (_H2O_SPEC_>=0) {
        beta=sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/(2.0*Kbol*rgTemperature));
        rgTotalPressure+=Kbol*rgTemperature*Pi*beta/2.0*FluxRamGauge[_H2O_SPEC_+i*PIC::nTotalSpecies];

//        rgTotalPressure+=4.0*Kbol*rgTemperature*FluxRamGauge[_H2O_SPEC_+i*PIC::nTotalSpecies]/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(_H2O_SPEC_)));
        ngTotalPressure+=DensityNudeGauge[_H2O_SPEC_+i*PIC::nTotalSpecies]*Kbol*ngTemperature;
      }

      if (_CO2_SPEC_>=0) {
        beta=sqrt(PIC::MolecularData::GetMass(_CO2_SPEC_)/(2.0*Kbol*rgTemperature));
        rgTotalPressure+=Kbol*rgTemperature*Pi*beta/2.0*FluxRamGauge[_CO2_SPEC_+i*PIC::nTotalSpecies];

//        rgTotalPressure+=4.0*Kbol*rgTemperature*FluxRamGauge[_CO2_SPEC_+i*PIC::nTotalSpecies]/sqrt(8.0*Kbol*rgTemperature/(Pi*PIC::MolecularData::GetMass(_CO2_SPEC_)));
        ngTotalPressure+=DensityNudeGauge[_CO2_SPEC_+i*PIC::nTotalSpecies]*Kbol*ngTemperature;
      }

      fprintf(fout,"%i %e %e %e %e %e %e %e %e %e %i %e %e %e %e\n",i, DensityRamGauge[spec+i*PIC::nTotalSpecies],DensityNudeGauge[spec+i*PIC::nTotalSpecies],
          rgPressure, ngPressure, rgTotalPressure, ngTotalPressure,
          Rosina[i].SecondsFromBegining,Rosina[i].RadiusVectorLeangth,Rosina[i].CometDistance,Rosina[i].LocationCode,Rosina[i].CharacteristicCellSize,
          Rosina[i].RamGauge.NucleusSolidAngle,Rosina[i].NudeGauge.NucleusSolidAngle,Rosina[i].Altitude);
    }

    fclose(fout);

    //sampeled flux
    sprintf(fname,"%s/RamNudeGaugeModelFlux.time=%s.spec=%s.out=%i.dat",PIC::OutputDataFileDirectory,RosinaSample::SamplingTimeInterval::CurrentTimeWindowStamp,PIC::MolecularData::GetChemSymbol(spec),nfile);
    fout=fopen(fname,"w");
    if (fout==0) exit(__LINE__,__FILE__,"Error: cannot open file");

    fprintf(fout,"VARIABLES=\"i\", \"Ram Gauge Flux\", \"Nude Gauge Flux\", \"Seconds From The First Point\", \" Radius-Vector Length\", \"Distance to the Coment\", \"Location Code\", \"Characteristic Cell Size\"\n");

    for (int i=0;i<nPoints;i++) {
      fprintf(fout,"%i %e %e %e %e %e %i %e\n",i,FluxRamGauge[spec+i*PIC::nTotalSpecies],FluxNudeGauge[spec+i*PIC::nTotalSpecies],
          Rosina[i].SecondsFromBegining,Rosina[i].RadiusVectorLeangth,Rosina[i].CometDistance,Rosina[i].LocationCode,Rosina[i].CharacteristicCellSize);
    }

    fclose(fout);
  }


  //flush sampled data buffers
 // Flush();
}



