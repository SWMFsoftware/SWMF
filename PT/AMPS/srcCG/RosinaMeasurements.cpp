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

//init the pointing directions and geomentry informationvoid
void RosinaSample::Init() {
#ifndef _NO_SPICE_CALLS_
  int i,idim,iCell,jCell,kCell,nd;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  //init line-of-sight vectors
  SpiceDouble lt,et,xRosetta[3];
  SpiceDouble       xform[6][6];

  FILE *fTrajectory=NULL;

  if (PIC::ThisThread==0) {
   char fname[1000];

   sprintf(fname,"%s/TrajectoryLastMeasuremetns.dat",PIC::OutputDataFileDirectory);
   fTrajectory=fopen(fname,"w");
   fprintf(fTrajectory,"VARIABLES=\"X\", \"Y\", \"Z\"\n");
  }

  for (i=0;i<nPoints;i++) {
    utc2et_c(ObservationTime[i],&et);
    spkpos_c("ROSETTA",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",xRosetta,&lt);
    for (idim=0;idim<3;idim++) xRosetta[idim]*=1.0E3;

    node=PIC::Mesh::mesh.findTreeNode(xRosetta,NULL);

    if (node!=NULL) {
      nd=PIC::Mesh::mesh.fingCellIndex(xRosetta,iCell,jCell,kCell,node);
    }
    else nd=-1,iCell=-1,jCell=-1,kCell-1;

    Rosina[i].nd=nd;
    Rosina[i].node=node;
    Rosina[i].iCell=iCell;
    Rosina[i].jCell=jCell;
    Rosina[i].kCell=kCell;

    if (PIC::ThisThread==0) fprintf(fTrajectory,"%e %e %e\n",xRosetta[0],xRosetta[1],xRosetta[2]);

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
  }

  if (PIC::ThisThread==0) {
    fclose(fTrajectory);
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
  int i,spec;
  long int ptr;
  double ParticleWeight,*v,c;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;

  #ifdef _NO_SPICE_CALLS_
  return; //cannot sample the instrument related data when SPICE is not used
  #endif

  for (i=0;i<nPoints;i++) if (Rosina[i].node!=NULL) if (Rosina[i].node->Thread==PIC::ThisThread) if ((block=Rosina[i].node->block)!=NULL) {
    ptr=block->FirstCellParticleTable[Rosina[i].iCell+_BLOCK_CELLS_X_*(Rosina[i].jCell+_BLOCK_CELLS_Y_*Rosina[i].kCell)];

    while (ptr!=-1) {
      spec=PIC::ParticleBuffer::GetI(ptr);
      v=PIC::ParticleBuffer::GetV(ptr);

      ParticleWeight=block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
      if ((cell=block->GetCenterNode(Rosina[i].nd))!=NULL) ParticleWeight/=cell->Measure;

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
}

//output sampled flux and density
void RosinaSample::PrintOutputFile(int nfile) {

  //collect all sampled data on the root processor
  double ExchangeBuffer[4*nPoints],SumBuffer[4*nPoints];
  int i,spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    for (i=0;i<nPoints;i++) {
      ExchangeBuffer[2*i]=Rosina[i].RamGauge.Density[spec]/PIC::LastSampleLength;
      ExchangeBuffer[2*i+1]=Rosina[i].RamGauge.Flux[spec]/PIC::LastSampleLength;

      ExchangeBuffer[2*i+2*nPoints]=Rosina[i].NudeGauge.Density[spec]/PIC::LastSampleLength;
      ExchangeBuffer[2*i+1+2*nPoints]=Rosina[i].NudeGauge.Flux[spec]/PIC::LastSampleLength;
    }

    MPI_Reduce(ExchangeBuffer,SumBuffer,4*nPoints,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);

    //output the data file
    if (PIC::ThisThread==0) {
      char fname[1000];
      FILE *fout=NULL;

      //sampeled density
      sprintf(fname,"%s/RamNudeGaugeModelDensity.spec=%s.out=%i.dat",PIC::OutputDataFileDirectory,PIC::MolecularData::GetChemSymbol(spec),nfile);
      fout=fopen(fname,"w");
      if (fout==0) exit(__LINE__,__FILE__,"Error: cannot open file");

      fprintf(fout,"VARIABLES=\"i\", \"Ram Gauge Density\", \"Nude Gauge Density\"\n");

      for (int i=0;i<nPoints;i++) {
        fprintf(fout,"%i %e %e\n",i,SumBuffer[2*i],SumBuffer[2*i+2*nPoints]);
      }

      fclose(fout);

      //sampeled flux
      sprintf(fname,"%s/RamNudeGaugeModelFlux.spec=%s.out=%i.dat",PIC::OutputDataFileDirectory,PIC::MolecularData::GetChemSymbol(spec),nfile);
      fout=fopen(fname,"w");
      if (fout==0) exit(__LINE__,__FILE__,"Error: cannot open file");

      fprintf(fout,"VARIABLES=\"i\", \"Ram Gauge Density\", \"Nude Gauge Density\"\n");

      for (int i=0;i<nPoints;i++) {
        fprintf(fout,"%i %e %e\n",i,SumBuffer[2*i+1],SumBuffer[2*i+1+2*nPoints]);
      }

      fclose(fout);
    }
  }

  //flush sampled data buffers
  Flush();
}



