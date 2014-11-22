
//$Id$

/*
 * pic_pt.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: vtenishe
 */

//contains routines that are used when AMPS is used in a 'particle tracker' mode

#include "pic.h"

long int PIC::ParticleTracker::ParticleDataRecordOffset=-1;

PIC::ParticleTracker::cTrajectoryRecord *PIC::ParticleTracker::TrajectoryDataBuffer::buffer=NULL;
unsigned long int PIC::ParticleTracker::TrajectoryDataBuffer::Size=500000;
unsigned long int PIC::ParticleTracker::TrajectoryDataBuffer::CurrentPosition=0;
unsigned long int PIC::ParticleTracker::TrajectoryDataBuffer::nfile=0;

unsigned long int PIC::ParticleTracker::TrajectoryList::Size=500000;
unsigned long int PIC::ParticleTracker::TrajectoryList::CurrentPosition=0;
PIC::ParticleTracker::cTrajectoryRecordReference *PIC::ParticleTracker::TrajectoryList::buffer=NULL;
unsigned long int PIC::ParticleTracker::TrajectoryList::nfile=0;


//init the particle tracker
void PIC::ParticleTracker::Init() {
  //reserve memory in the particle data
  PIC::ParticleBuffer::RequestDataStorage(ParticleDataRecordOffset,sizeof(cParticleData));

  //init the data buffers
  TrajectoryDataBuffer::buffer=new cTrajectoryRecord[TrajectoryDataBuffer::Size];
  TrajectoryList::buffer=new cTrajectoryRecordReference[TrajectoryList::Size];

  //remove old and create new directory for temporary files
  char cmd[_MAX_STRING_LENGTH_PIC_];
  sprintf(cmd,"rm -rf %s/ParticleTrackerTmp",PIC::OutputDataFileDirectory);
  system(cmd);

  sprintf(cmd,"mkdir -p %s/ParticleTrackerTmp",PIC::OutputDataFileDirectory);
  system(cmd);
}



//init the particle trajecotry record
void PIC::ParticleTracker::InitParticleID(void *ParticleData) {
  cParticleData *DataRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  DataRecord->TrajectoryTrackingFlag=false;

  DataRecord->lastRecordReference.file=0;
  DataRecord->lastRecordReference.offset=-1;
  DataRecord->lastRecordReference.thread=0;
}

void PIC::ParticleTracker::TrajectoryDataBuffer::flush() {
  FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];

  sprintf(fname,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryData.pt",PIC::OutputDataFileDirectory,PIC::ThisThread,nfile);
  fout=fopen(fname,"w");

  fwrite(buffer,sizeof(cTrajectoryRecord),CurrentPosition,fout);
  fclose(fout);

  CurrentPosition=0;
  ++nfile;
}

void PIC::ParticleTracker::TrajectoryList::flush() {
  FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];

  sprintf(fname,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryList.pt",PIC::OutputDataFileDirectory,PIC::ThisThread,nfile);
  fout=fopen(fname,"w");

  fwrite(&CurrentPosition,sizeof(unsigned long int),1,fout);
  fwrite(buffer,sizeof(cTrajectoryRecordReference),CurrentPosition,fout);
  fclose(fout);

  CurrentPosition=0;
  ++nfile;
}

void PIC::ParticleTracker::RecordTrajectoryPoint(double *x,double *v,int spec,void *ParticleData) {
  cParticleData *ParticleTrajectoryRecord;
  cTrajectoryRecord *TrajectoryRecord;

  //save the data buffer if full
  if (TrajectoryDataBuffer::CurrentPosition==TrajectoryDataBuffer::Size) TrajectoryDataBuffer::flush();

  //save trajectory point
  TrajectoryRecord=TrajectoryDataBuffer::buffer+TrajectoryDataBuffer::CurrentPosition;
  ParticleTrajectoryRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  //the trajectory is traced only if the particle trajecotry tracking flag is set 'true'
  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==false) return;

  //save physical data
  memcpy(TrajectoryRecord->x,x,3*sizeof(double));
  TrajectoryRecord->Speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  TrajectoryRecord->spec=spec;

  //save trajectory tracking data
  TrajectoryRecord->lastRecordReference=ParticleTrajectoryRecord->lastRecordReference;

  //update the trajecotry tracking data
  ParticleTrajectoryRecord->lastRecordReference.file=TrajectoryDataBuffer::nfile;
  ParticleTrajectoryRecord->lastRecordReference.thread=PIC::ThisThread;
  ParticleTrajectoryRecord->lastRecordReference.offset=TrajectoryDataBuffer::CurrentPosition;

  //update the trajectory data buffer pointer
  ++TrajectoryDataBuffer::CurrentPosition;
}

void PIC::ParticleTracker::FinilazeParticleRecord(void *ParticleData) {
  cParticleData *ParticleTrajectoryRecord;

  //the trajectory is traced only if the particle trajecotry tracking flag is set 'true'
  ParticleTrajectoryRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));
  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==false) return;

  //save the data buffer if full
  if (TrajectoryList::CurrentPosition==TrajectoryList::Size) TrajectoryList::flush();
  TrajectoryList::buffer[TrajectoryList::CurrentPosition]=ParticleTrajectoryRecord->lastRecordReference;

  //update the trajectory data buffer pointer
  ++TrajectoryList::CurrentPosition;
}

void PIC::ParticleTracker::CreateTrajectoryFile(const char *fname) {
  FILE *fout[PIC::nTotalSpecies];
  FILE *fTrajectoryList=NULL,*fTrajectoryData=NULL;
  char str[_MAX_STRING_LENGTH_PIC_];
  unsigned long int t;
  int i,spec,thread,nList;
  unsigned int TrajectoryCounter[PIC::nTotalSpecies];


  //the number of files, and size of the last file in the list
  unsigned long int TrajectoryList_nfile[PIC::nTotalThreads];

  //trajectory list
  t=TrajectoryList::nfile+1; //one more file will be written when the buffer is flushed
  MPI_Gather(&t,1,MPI_UNSIGNED_LONG,TrajectoryList_nfile,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);

  //trajectory buffer
  TrajectoryDataBuffer::flush();
  TrajectoryList::flush();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //the file will be created only by the root processor
  if (PIC::ThisThread==0) {
    //open nessesary files
    //output trajectory files
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      char ChemSymbol[_MAX_STRING_LENGTH_PIC_];

      TrajectoryCounter[spec]=0;

      PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);
      sprintf(str,"%s.s=%i.%s.dat",fname,spec,ChemSymbol);

      fout[spec]=fopen(str,"w");
      fprintf(fout[spec],"VARIABLES=\"x\",\"y\",\"z\",\"Speed\"\n");
    }

    //create particle trajectories
    unsigned long int lastRecordThread=-1,lastRecordFile=-1;

    for (thread=0;thread<PIC::nTotalThreads;thread++) for (nList=0;nList<TrajectoryList_nfile[thread];nList++) {
      unsigned long int tr,TotalListLegth;
      int StartTrajectorySpec; //during the particle motion the spec might be changed, which will change the trajectory output file
      unsigned long int RecordThread,RecordOffset,RecordFile;
      FILE *trOut;

      cTrajectoryRecordReference StartTrajectoryPoint;
      cTrajectoryRecord TrajectoryRecord;

      sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryList.pt",PIC::OutputDataFileDirectory,thread,nList);
      fTrajectoryList=fopen(str,"r");
      fread(&TotalListLegth,sizeof(unsigned long int),1,fTrajectoryList);

      for (tr=0;tr<TotalListLegth;tr++) {
        fread(&StartTrajectoryPoint,sizeof(cTrajectoryRecordReference),1,fTrajectoryList);
        StartTrajectorySpec=-1,trOut=NULL;

        do {
          if (StartTrajectorySpec==-1) {
            RecordThread=StartTrajectoryPoint.thread;
            RecordOffset=StartTrajectoryPoint.offset;
            RecordFile=StartTrajectoryPoint.file;
          }
          else {
            RecordThread=TrajectoryRecord.lastRecordReference.thread;
            RecordOffset=TrajectoryRecord.lastRecordReference.offset;
            RecordFile=TrajectoryRecord.lastRecordReference.file;
          }

          //open the trajectory data file if needed
          if ((RecordThread!=lastRecordThread)||(RecordFile!=lastRecordFile)) {
            //close previously opened file
            if (lastRecordFile!=-1) {
              fclose(fTrajectoryData);
            }

            //open a new trajecory data file
            lastRecordThread=RecordThread,lastRecordFile=RecordFile;

            sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryData.pt",PIC::OutputDataFileDirectory,RecordThread,RecordFile);
            fTrajectoryData=fopen(str,"r");
          }

          fseek (fTrajectoryData,RecordOffset*sizeof(cTrajectoryRecord),SEEK_SET);
          fread(&TrajectoryRecord,sizeof(cTrajectoryRecord),1,fTrajectoryData);

          #if _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE_ == _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE__ENTIRE_TRAJECTORY_
          if (StartTrajectorySpec==-1) {
            trOut=fout[TrajectoryRecord.spec];

            //print the header of the new trajectory
            fprintf(trOut,"ZONE T=\"Trajectory=%i\" F=POINT\n",TrajectoryCounter[TrajectoryRecord.spec]);
            ++TrajectoryCounter[TrajectoryRecord.spec];
            StartTrajectorySpec=TrajectoryRecord.spec;
          }
          #elif _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE_ == _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE__SPECIES_TYPE_SEGMENTS_
          if ((StartTrajectorySpec==-1)||(StartTrajectorySpec!=TrajectoryRecord.spec)) {
            trOut=fout[TrajectoryRecord.spec];

            //print the header of the new trajectory
            fprintf(trOut,"ZONE T=\"Trajectory=%i\" F=POINT\n",TrajectoryCounter[TrajectoryRecord.spec]);
            ++TrajectoryCounter[TrajectoryRecord.spec];
            StartTrajectorySpec=TrajectoryRecord.spec;
          }
          #else
          exit(__LINE__,__FILE__,"Error: unknown option");
          #endif

          fprintf(trOut,"%e  %e  %e  %e\n",TrajectoryRecord.x[0],TrajectoryRecord.x[1],TrajectoryRecord.x[2],TrajectoryRecord.Speed);
        }
        while (TrajectoryRecord.lastRecordReference.offset!=-1);
      }

      fclose(fTrajectoryList);
    }

    //close the trajectory files
    for (spec=0;spec<PIC::nTotalSpecies;spec++) fclose(fout[spec]);

    if (fTrajectoryData!=NULL) fclose(fTrajectoryData);
  }

  //clean up the trajectory data buffers and set the default value of the particle tracking flag
/*  TrajectoryDataBuffer::clean();
  TrajectoryList::clean();

  SetDefaultParticleTrackingFlag();*/
}


//set up the default tracking flag to all particles
void PIC::ParticleTracker::SetDefaultParticleTrackingFlag(void* StartNodeVoid) {
  int i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode,*StartNode;
  cParticleData *DataRecord;
  long int ptr;

  StartNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) StartNodeVoid;
  if (StartNode==NULL) StartNode=PIC::Mesh::mesh.rootTree;

  //serch the tree
  for (i=0;i<(1<<DIM);i++) if ((downNode=StartNode->downNode[i])!=NULL) SetDefaultParticleTrackingFlag(downNode);

  //reset the particles
  if (StartNode->block!=NULL) {
    long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

    memcpy(FirstCellParticleTable,StartNode->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++)  {
      ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      while (ptr!=-1) {
        DataRecord=(cParticleData*)(ParticleDataRecordOffset+PIC::ParticleBuffer::GetParticleDataPointer(ptr));
        DataRecord->TrajectoryTrackingFlag=false;

        ptr=PIC::ParticleBuffer::GetNext(ptr);
      }
    }
  }
}

//set the particle tracking flag "on"
void PIC::ParticleTracker::StartParticleTrajectoryTracking(void *ParticleData) {
  cParticleData *DataRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  DataRecord->TrajectoryTrackingFlag=true;
}

void PIC::ParticleTracker::StopParticleTrajectoryTracking(void *ParticleData) {
  cParticleData *DataRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  if (DataRecord->TrajectoryTrackingFlag==true) FinilazeParticleRecord(ParticleData);
  DataRecord->TrajectoryTrackingFlag=false;
}

//the default particle trajectory tracking condition
bool PIC::ParticleTracker::TrajectoryTrackingCondition_default(double *x,double *v,int spec,void *ParticleData) {
  return false;
}

//apply the particle tracking condition to a particle
void PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  bool flag;
  cParticleData *DataRecord;

  flag=_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_(x,v,spec,ParticleData);

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  DataRecord=(cParticleData*)(ParticleDataRecordOffset+(PIC::ParticleBuffer::byte*)ParticleData);
  DataRecord->TrajectoryTrackingFlag=flag;

  if (flag==true) RecordTrajectoryPoint(x,v,spec,ParticleData);
  #endif
}

//apply the particle tracking condition to all particles in the simulation
void PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(void* StartNodeVoid) {
  int i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode,*StartNode;
  long int ptr;

  double x[3],v[3];
  int spec;
  PIC::ParticleBuffer::byte* ParticleData;

  StartNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) StartNodeVoid;
  if (StartNode==NULL) StartNode=PIC::Mesh::mesh.rootTree;

  //serch the tree
  for (i=0;i<(1<<DIM);i++) if ((downNode=StartNode->downNode[i])!=NULL) SetDefaultParticleTrackingFlag(downNode);

  //reset the particles
  if (StartNode->block!=NULL) {
    long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

    memcpy(FirstCellParticleTable,StartNode->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++)  {
      ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        PIC::ParticleBuffer::GetX(x,ParticleData);
        PIC::ParticleBuffer::GetV(v,ParticleData);

        ApplyTrajectoryTrackingCondition(x,v,spec,ParticleData);

        ptr=PIC::ParticleBuffer::GetNext(ptr);
      }
    }
  }
}










