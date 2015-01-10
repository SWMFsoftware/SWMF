
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

PIC::ParticleTracker::cTrajectoryDataRecord *PIC::ParticleTracker::TrajectoryData::buffer=NULL;
unsigned long int PIC::ParticleTracker::TrajectoryData::Size=500000;
unsigned long int PIC::ParticleTracker::TrajectoryData::CurrentPosition=0;
unsigned long int PIC::ParticleTracker::TrajectoryData::nfile=0;

unsigned long int PIC::ParticleTracker::TrajectoryList::Size=500000;
unsigned long int PIC::ParticleTracker::TrajectoryList::CurrentPosition=0;
PIC::ParticleTracker::cTrajectoryListRecord *PIC::ParticleTracker::TrajectoryList::buffer=NULL;
unsigned long int PIC::ParticleTracker::TrajectoryList::nfile=0;

int PIC::ParticleTracker::maxSampledTrajectoryNumber=-1;
int *PIC::ParticleTracker::threadSampledTrajectoryNumber=NULL;
int *PIC::ParticleTracker::totalSampledTrajectoryNumber=NULL;
unsigned int PIC::ParticleTracker::SampledTrajectoryCounter=0;


//init the particle tracker
void PIC::ParticleTracker::Init() {
  //reserve memory in the particle data
  PIC::ParticleBuffer::RequestDataStorage(ParticleDataRecordOffset,sizeof(cParticleData));

  //init the data buffers
  TrajectoryData::buffer=new cTrajectoryDataRecord[TrajectoryData::Size];
  TrajectoryList::buffer=new cTrajectoryListRecord[TrajectoryList::Size];

  //init the trajectory counter
  threadSampledTrajectoryNumber=new int [PIC::nTotalSpecies];
  totalSampledTrajectoryNumber=new int [PIC::nTotalSpecies];
  for (int s=0;s<PIC::nTotalSpecies;s++) threadSampledTrajectoryNumber[s]=0,totalSampledTrajectoryNumber[s]=0;

  //remove old and create new directory for temporary files
  char cmd[_MAX_STRING_LENGTH_PIC_];
  sprintf(cmd,"rm -rf %s/ParticleTrackerTmp",PIC::OutputDataFileDirectory);
  system(cmd);

  sprintf(cmd,"mkdir -p %s/ParticleTrackerTmp",PIC::OutputDataFileDirectory);
  system(cmd);
}

//update the total number of samples trajectories
void PIC::ParticleTracker::UpdateTrajectoryCounter() {
  MPI_Allreduce(threadSampledTrajectoryNumber,totalSampledTrajectoryNumber,PIC::nTotalSpecies,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
}

//init the particle trajecotry record
void PIC::ParticleTracker::InitParticleID(void *ParticleData) {
  cParticleData *DataRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  DataRecord->TrajectoryTrackingFlag=false;
  DataRecord->nSampledTrajectoryPoints=0;

  DataRecord->Trajectory.StartingThread=0;
  DataRecord->Trajectory.id=0;
}

void PIC::ParticleTracker::TrajectoryData::flush() {
  FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];

  if (CurrentPosition!=0) {
    sprintf(fname,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryData.pt",PIC::OutputDataFileDirectory,PIC::ThisThread,nfile);
    fout=fopen(fname,"w");

    fwrite(&CurrentPosition,sizeof(unsigned long int),1,fout);
    fwrite(buffer,sizeof(cTrajectoryDataRecord),CurrentPosition,fout);
    fclose(fout);

    CurrentPosition=0;
    ++nfile;
  }
}

void PIC::ParticleTracker::TrajectoryList::flush() {
  FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];

  if (CurrentPosition!=0) {
    sprintf(fname,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryList.pt",PIC::OutputDataFileDirectory,PIC::ThisThread,nfile);
    fout=fopen(fname,"w");

    fwrite(&CurrentPosition,sizeof(unsigned long int),1,fout);
    fwrite(buffer,sizeof(cTrajectoryListRecord),CurrentPosition,fout);
    fclose(fout);

    CurrentPosition=0;
    ++nfile;
  }
}

void PIC::ParticleTracker::RecordTrajectoryPoint(double *x,double *v,int spec,void *ParticleData) {
  cParticleData *ParticleTrajectoryRecord;
  cTrajectoryDataRecord *TrajectoryRecord;

  //save the data buffer if full
  if (TrajectoryData::CurrentPosition==TrajectoryData::Size) TrajectoryData::flush();

  //pointers to the trajectory data buffer record and the particle trajectory data
  TrajectoryRecord=TrajectoryData::buffer+TrajectoryData::CurrentPosition;
  ParticleTrajectoryRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  //the trajectory is traced only if the particle trajecotry tracking flag is set 'true'
  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==false) return;

  //save physical data
  memcpy(TrajectoryRecord->data.x,x,3*sizeof(double));
  TrajectoryRecord->data.Speed=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  TrajectoryRecord->data.spec=spec;
  TrajectoryRecord->offset=ParticleTrajectoryRecord->nSampledTrajectoryPoints;

  //save the trajectory ID
  TrajectoryRecord->Trajectory=ParticleTrajectoryRecord->Trajectory;

  //increment the trajectory point counter
  ParticleTrajectoryRecord->nSampledTrajectoryPoints++;
  if (ParticleTrajectoryRecord->nSampledTrajectoryPoints==0) exit(__LINE__,__FILE__,"Error: ParticleTrajectoryRecord->nSampledTrajectoryPoints is zero");

  //update the trajectory data buffer pointer
  TrajectoryData::CurrentPosition++;
}

void PIC::ParticleTracker::FinilazeParticleRecord(void *ParticleData) {
  cParticleData *ParticleTrajectoryRecord;

  //the trajectory is traced only if the particle trajecotry tracking flag is set 'true'
  ParticleTrajectoryRecord=(cParticleData*)(ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));
  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==false) return;

  //save the data buffer if full
  if (TrajectoryList::CurrentPosition==TrajectoryList::Size) TrajectoryList::flush();
  TrajectoryList::buffer[TrajectoryList::CurrentPosition].Trajectory=ParticleTrajectoryRecord->Trajectory;
  TrajectoryList::buffer[TrajectoryList::CurrentPosition].nSampledTrajectoryPoints=ParticleTrajectoryRecord->nSampledTrajectoryPoints;

  //update the trajectory data buffer pointer
  ++TrajectoryList::CurrentPosition;

  //reset the tracking flag
  ParticleTrajectoryRecord->TrajectoryTrackingFlag=false;
}

//output sampled particle trajectories
//PIC::ParticleTracker::cTrajectoryRecordReference PIC::ParticleTracker::ParticleTrajectoryFile::StartTrajectoryPoint;
//int PIC::ParticleTracker::ParticleTrajectoryFile::lastRecordThread=-1,PIC::ParticleTracker::ParticleTrajectoryFile::lastRecordFile=-1;
/*FILE *PIC::ParticleTracker::ParticleTrajectoryFile::fout[PIC::nTotalSpecies];
FILE *PIC::ParticleTracker::ParticleTrajectoryFile::fTrajectoryData=NULL;
unsigned int PIC::ParticleTracker::ParticleTrajectoryFile::TrajectoryCounter[PIC::nTotalSpecies];
char PIC::ParticleTracker::ParticleTrajectoryFile::str[_MAX_STRING_LENGTH_PIC_];*/

/*
void PIC::ParticleTracker::ParticleTrajectoryFile::OutputParticleTrajectory() {
  unsigned long int RecordThread,RecordOffset,RecordFile;
  cTrajectoryDataRecord TrajectoryRecord;
  int StartTrajectorySpec; //during the particle motion the spec might be changed, which will change the trajectory output file
  FILE *trOut;

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

unsigned long int PIC::ParticleTracker::ParticleTrajectoryFile::nTrajectoryDataFiles[PIC::nTotalThreads];
unsigned long int PIC::ParticleTracker::ParticleTrajectoryFile::nTrajectoryListFiles[PIC::nTotalThreads];
unsigned long int PIC::ParticleTracker::ParticleTrajectoryFile::nTotalSampledTrajectories[PIC::nTotalThreads];

void PIC::ParticleTracker::ParticleTrajectoryFile::SaveCurrectParticlesTrajectoryList() {

  //output trajecotries of the particles from the root processor
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            memcpy(&ParticleTrajectoryRecord,ParticleDataRecordOffset+ParticleData,sizeof(cParticleData));

            if ((ParticleTrajectoryRecord.TrajectoryTrackingFlag==true)&&(ParticleTrajectoryRecord.lastRecordReference.offset!=-1)) {
              StartTrajectoryPoint=ParticleTrajectoryRecord.lastRecordReference;
              OutputParticleTrajectory(); //the starting point of the trajectory is in StartTrajectoryPoint
            }

            ptr=PIC::ParticleBuffer::GetNext(ParticleData);
          }

        }
      }
    }
  }
}
*/


void PIC::ParticleTracker::OutputTrajectory(const char *fname) {
  int thread;
  char str[_MAX_STRING_LENGTH_PIC_];

  //save a temporary file that contains the trajectory information of the particles that are currently in the system
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],ptr;
  long int i,j,k;
  PIC::ParticleBuffer::byte *ParticleData;
  cParticleData ParticleTrajectoryRecord;
  FILE *fTemporatyTrajectoryList;
  unsigned long int length=0;

  sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.TemporaryTrajectoryList.pt",PIC::OutputDataFileDirectory,PIC::ThisThread);
  fTemporatyTrajectoryList=fopen(str,"w");

  //calculate the number of the particles that will be placed into the list
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            memcpy(&ParticleTrajectoryRecord,ParticleDataRecordOffset+ParticleData,sizeof(cParticleData));

            if (ParticleTrajectoryRecord.TrajectoryTrackingFlag==true) {
              ++length;
            }

            ptr=PIC::ParticleBuffer::GetNext(ParticleData);
          }

        }
      }
    }
  }

  fwrite(&length,sizeof(unsigned long int),1,fTemporatyTrajectoryList);

  //save the list
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            memcpy(&ParticleTrajectoryRecord,ParticleDataRecordOffset+ParticleData,sizeof(cParticleData));

            if (ParticleTrajectoryRecord.TrajectoryTrackingFlag==true) {
              cTrajectoryListRecord record;

              record.Trajectory=ParticleTrajectoryRecord.Trajectory;
              record.nSampledTrajectoryPoints=ParticleTrajectoryRecord.nSampledTrajectoryPoints;
              if (record.nSampledTrajectoryPoints<=0) exit(__LINE__,__FILE__,"Error: the number of sampled points is out of range");

              fwrite(&record,sizeof(cTrajectoryListRecord),1,fTemporatyTrajectoryList);
            }

            ptr=PIC::ParticleBuffer::GetNext(ParticleData);
          }

        }
      }
    }
  }

  fclose(fTemporatyTrajectoryList);

  //flush trajectory buffer
  TrajectoryData::flush();
  TrajectoryList::flush();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //get the number of the trajectory list and data files, and the total number of the sampled trajectories
  unsigned long nTrajectoryListFiles[PIC::nTotalThreads],nTrajectoryDataFiles[PIC::nTotalThreads],nTotalSampledTrajectories[PIC::nTotalThreads];

  MPI_Gather(&TrajectoryList::nfile,1,MPI_UNSIGNED_LONG,nTrajectoryListFiles,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Gather(&TrajectoryData::nfile,1,MPI_UNSIGNED_LONG,nTrajectoryDataFiles,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Gather(&SampledTrajectoryCounter,1,MPI_UNSIGNED_LONG,nTotalSampledTrajectories,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);

  //output trajectories
  if (PIC::ThisThread==0) {
    //open sampled trajectory output files
    FILE *fout[PIC::nTotalSpecies];
    unsigned int spec, TrajectoryCounter[PIC::nTotalSpecies];

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      char ChemSymbol[_MAX_STRING_LENGTH_PIC_];

      TrajectoryCounter[spec]=0;

      PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);
      sprintf(str,"%s.s=%i.%s.dat",fname,spec,ChemSymbol);

      fout[spec]=fopen(str,"w");
      fprintf(fout[spec],"VARIABLES=\"x\",\"y\",\"z\",\"Speed\"\n");
    }

    //Create the array of sampled trajectory points for all trajectories (used also as a flag to determine trajectories that was already processed)
    long int nTotalTracedTrajectories=0;
    long int TrajectoryCounterOffset[PIC::nTotalThreads]; //the "global" trajectory number is TrajectoryCounterOffset[TrajectoryStartingThread]+ trajectory counting number at partuculat processor

    for (thread=0;thread<PIC::nTotalThreads;thread++) {
      TrajectoryCounterOffset[thread]=nTotalTracedTrajectories;
      nTotalTracedTrajectories+=nTotalSampledTrajectories[thread];
    }

    //2. Read the table of the sampled trajectory numbers
    int nfile;
    FILE *fTrajectoryList;
    PIC::ParticleTracker::cTrajectoryID Trajectory;
//    unsigned int nTrajectoryPoints;
    unsigned long int i,length,nReadTrajectoryNumber=0;

    int *nSampledTrajectoryPoints=new int [nTotalTracedTrajectories];
    int *SampledTrajectoryDataOffset=new int [nTotalTracedTrajectories];
    int *nReadSampledTrajectoryPoints=new int [nTotalTracedTrajectories];

    for (thread=0;thread<PIC::nTotalThreads;thread++) for (nfile=0;nfile<nTrajectoryListFiles[thread]+1;nfile++) {
      //scroll through all trajectory lists inclusing the temporary lists that contain the trajectory information regarding particles that are still in the simulation
      if (nfile<nTrajectoryListFiles[thread]) {
        sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryList.pt",PIC::OutputDataFileDirectory,thread,nfile);
      }
      else {
        sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.TemporaryTrajectoryList.pt",PIC::OutputDataFileDirectory,thread);
      }

      fTrajectoryList=fopen(str,"r");
      fread(&length,sizeof(unsigned long int),1,fTrajectoryList);

//      while (feof(fTrajectoryList)==false) {
      for (i=0;i<length;i++) {
        cTrajectoryListRecord Record;

        fread(&Record,sizeof(PIC::ParticleTracker::cTrajectoryListRecord),1,fTrajectoryList);

        nSampledTrajectoryPoints[Record.Trajectory.id+TrajectoryCounterOffset[Record.Trajectory.StartingThread]]=Record.nSampledTrajectoryPoints;
        nReadTrajectoryNumber++;
       }

      fclose(fTrajectoryList);
    }

    if (nReadTrajectoryNumber!=nTotalTracedTrajectories) {
      exit(__LINE__,__FILE__,"Error: the number of the read trajectories is different from the total number of the sampled trajectories");
    }

    //2. Read the particle trajectories
    nReadTrajectoryNumber=0;

    static const int TrajectoryPointBufferLength=10000000;
    cTrajectoryPhysicalData *TempTrajectoryBuffer=new cTrajectoryPhysicalData[TrajectoryPointBufferLength];
    int UsedTrajectoryPointBuffer,StartTrajectoryNumber;

    while (nReadTrajectoryNumber!=nTotalTracedTrajectories) {
      UsedTrajectoryPointBuffer=0;
      StartTrajectoryNumber=nReadTrajectoryNumber;

      //reset the offset array
      for (i=0;i<nTotalTracedTrajectories;i++) SampledTrajectoryDataOffset[i]=-1,nReadSampledTrajectoryPoints[i]=0;

      while (UsedTrajectoryPointBuffer+nSampledTrajectoryPoints[nReadTrajectoryNumber]<=TrajectoryPointBufferLength) {
        SampledTrajectoryDataOffset[nReadTrajectoryNumber]=UsedTrajectoryPointBuffer;
        UsedTrajectoryPointBuffer+=nSampledTrajectoryPoints[nReadTrajectoryNumber];
        nReadTrajectoryNumber++;

        if (nReadTrajectoryNumber==nTotalTracedTrajectories) break;
      }

      //scroll throught all trajectory files and locate all points that corresponds to this trajectories
      FILE *fTrajectoryData;
      cTrajectoryDataRecord TrajectoryRecord;
      int GlobalTrajectoryNumber,ReadTrajectoryPoints=0;

      for (thread=0;thread<PIC::nTotalThreads;thread++) {
        for (nfile=0;nfile<nTrajectoryDataFiles[thread];nfile++) {
          sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryData.pt",PIC::OutputDataFileDirectory,thread,nfile);

          fTrajectoryData=fopen(str,"r");
          fread(&length,sizeof(unsigned long int),1,fTrajectoryList);

//          while (feof(fTrajectoryData)==false) {
          for (i=0;i<length;i++) {
            fread(&TrajectoryRecord,sizeof(cTrajectoryDataRecord),1,fTrajectoryData);
            GlobalTrajectoryNumber=TrajectoryRecord.Trajectory.id+TrajectoryCounterOffset[TrajectoryRecord.Trajectory.StartingThread];

            if (SampledTrajectoryDataOffset[GlobalTrajectoryNumber]!=-1) {
              TempTrajectoryBuffer[SampledTrajectoryDataOffset[GlobalTrajectoryNumber]+TrajectoryRecord.offset]=TrajectoryRecord.data;
              ++nReadSampledTrajectoryPoints[GlobalTrajectoryNumber];
              ++ReadTrajectoryPoints;

              if (nReadSampledTrajectoryPoints[GlobalTrajectoryNumber]>nSampledTrajectoryPoints[GlobalTrajectoryNumber]) {
                exit(__LINE__,__FILE__,"Error: the number of the read trajectory number is out of range");
              }
            }
          }

          fclose(fTrajectoryData);

          //if all points are found stop reading the trajectory files
          if (ReadTrajectoryPoints==UsedTrajectoryPointBuffer) {
            break;
          }
        }

        //if all points are found stop reading the trajectory files
        if (ReadTrajectoryPoints==UsedTrajectoryPointBuffer) {
          break;
        }
      }

      //check whether all trajecotry points are found
      if (ReadTrajectoryPoints!=UsedTrajectoryPointBuffer) {
        exit(__LINE__,__FILE__,"Error: some trajectory points are not found");
      }

      //save the found trajectories
      long int tr,offset;
      int StartTrajectorySpec;
      cTrajectoryPhysicalData *TrajectoryData;
      FILE *trOut;

      for (tr=StartTrajectoryNumber;tr<nReadTrajectoryNumber;tr++) {
        offset=SampledTrajectoryDataOffset[tr];
        StartTrajectorySpec=-1;

        for (i=0;i<nSampledTrajectoryPoints[tr];i++) {
          TrajectoryData=TempTrajectoryBuffer+offset+i;

          #if _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE_ == _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE__ENTIRE_TRAJECTORY_
          if (StartTrajectorySpec==-1) {
            trOut=fout[TrajectoryData->spec];

            //print the header of the new trajectory
            fprintf(trOut,"ZONE T=\"Trajectory=%i\" F=POINT\n",TrajectoryCounter[TrajectoryData->spec]);
            ++TrajectoryCounter[TrajectoryData->spec];
            StartTrajectorySpec=TrajectoryData->spec;
          }
          #elif _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE_ == _PIC_PARTICLE_TRACKER__TRAJECTORY_OUTPUT_MODE__SPECIES_TYPE_SEGMENTS_
          if ((StartTrajectorySpec==-1)||(StartTrajectorySpec!=TrajectoryData->spec)) {
            trOut=fout[TrajectoryData->spec];

            //print the header of the new trajectory
            fprintf(trOut,"ZONE T=\"Trajectory=%i\" F=POINT\n",TrajectoryCounter[TrajectoryData->spec]);
            ++TrajectoryCounter[TrajectoryData->spec];
            StartTrajectorySpec=TrajectoryData->spec;
          }
          #else
          exit(__LINE__,__FILE__,"Error: unknown option");
          #endif

          fprintf(trOut,"%e  %e  %e  %e\n",TrajectoryData->x[0],TrajectoryData->x[1],TrajectoryData->x[2],TrajectoryData->Speed);
        }
      }

    }

    //close all open files and deallocate the temporatry data buffers
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      fclose(fout[spec]);
      fout[spec]=NULL;
    }

    delete [] TempTrajectoryBuffer;
    delete [] nSampledTrajectoryPoints;
    delete [] SampledTrajectoryDataOffset;
    delete [] nReadSampledTrajectoryPoints;
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

/*
void PIC::ParticleTracker::ParticleTrajectoryFile::Output(const char *fname) {
  FILE *fTrajectoryList=NULL;
//  unsigned long int t;
  int i,spec,thread,nList;

  //the number of files, and size of the last file in the list

  unsigned long int nTrajectoryDataFiles[PIC::nTotalThreads];
  unsigned long int nTrajectoryListFiles[PIC::nTotalThreads];
  unsigned long int nTotalSampledTrajectories[PIC::nTotalThreads];


  //flush trajectory buffer
  TrajectoryData::flush();
  TrajectoryList::flush();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //get the number of the trajectory list and data files, and the total number of the sampled trajectories
  MPI_Gather(&TrajectoryList::nfile,1,MPI_UNSIGNED_LONG,nTrajectoryListFiles,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Gather(&TrajectoryData::nfile,1,MPI_UNSIGNED_LONG,nTrajectoryDataFiles,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Gather(&TrajectoryData::nfile,1,MPI_UNSIGNED_LONG,nTotalSampledTrajectories,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);

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
    for (thread=0;thread<PIC::nTotalThreads;thread++) for (nList=0;nList<TrajectoryList_nfile[thread];nList++) {
      unsigned long int tr,TotalListLegth;
      cTrajectoryRecord TrajectoryRecord;

      sprintf(str,"%s/ParticleTrackerTmp/amps.ParticleTracker.thread=%i.out=%i.TrajectoryList.pt",PIC::OutputDataFileDirectory,thread,nList);
      fTrajectoryList=fopen(str,"r");
      fread(&TotalListLegth,sizeof(unsigned long int),1,fTrajectoryList);

      for (tr=0;tr<TotalListLegth;tr++) {
        fread(&StartTrajectoryPoint,sizeof(cTrajectoryRecordReference),1,fTrajectoryList);
        OutputParticleTrajectory(); //the starting point of the trajectory is in StartTrajectoryPoint
      }

      fclose(fTrajectoryList);
    }
  }

  //output trajectories of the particles that are still in the system
  if (PIC::ThisThread==0) {
    CMPI_channel pipe(1000000);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],ptr;
    int i,j,k;
    PIC::ParticleBuffer::byte *ParticleData;
    cParticleData ParticleTrajectoryRecord;


    //output trajecotries of the particles from the root processor
    for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
      memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              memcpy(&ParticleTrajectoryRecord,ParticleDataRecordOffset+ParticleData,sizeof(cParticleData));

              if ((ParticleTrajectoryRecord.TrajectoryTrackingFlag==true)&&(ParticleTrajectoryRecord.lastRecordReference.offset!=-1)) {
                StartTrajectoryPoint=ParticleTrajectoryRecord.lastRecordReference;
                OutputParticleTrajectory(); //the starting point of the trajectory is in StartTrajectoryPoint
              }

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }

          }
        }
      }

    }

    //output trajecotries of the particles stored on other processors
    pipe.openRecvAll();

    for (int thread=1;thread<PIC::nTotalThreads;thread++) {
      pipe.recv(ptr,thread);

      while (ptr!=-1) {
        pipe.recv(StartTrajectoryPoint,thread);
        OutputParticleTrajectory(); //the starting point of the trajectory is in StartTrajectoryPoint

        pipe.recv(ptr,thread);
      }
    }

    pipe.closeRecvAll();
  }
  else {
    CMPI_channel pipe(1000000);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],ptr;
    int i,j,k;
    PIC::ParticleBuffer::byte *ParticleData;
    cParticleData ParticleTrajectoryRecord;

    pipe.openSend(0);

    for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
      memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              memcpy(&ParticleTrajectoryRecord,ParticleDataRecordOffset+ParticleData,sizeof(cParticleData));

              if ((ParticleTrajectoryRecord.TrajectoryTrackingFlag==true)&&(ParticleTrajectoryRecord.lastRecordReference.offset!=-1)) {
                StartTrajectoryPoint=ParticleTrajectoryRecord.lastRecordReference;
                pipe.send(ptr);
                pipe.send(StartTrajectoryPoint);
              }

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }

          }
        }
      }
    }

    ptr=-1;
    pipe.send(ptr);

    pipe.closeSend();
  }

  //close open files
  if (PIC::ThisThread==0) {
    //close the trajectory files
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      fclose(fout[spec]);
      fout[spec]=NULL;
    }

    if (fTrajectoryData!=NULL) {
      fclose(fTrajectoryData);
      fTrajectoryData=NULL,lastRecordThread=-1,lastRecordFile=-1;
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}
*/


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
  if (maxSampledTrajectoryNumber>totalSampledTrajectoryNumber[spec]) {
    return true;
  }

  return false;
}

//apply the particle tracking condition to a particle
void PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  bool flag;
  cParticleData *DataRecord=(cParticleData*)(ParticleDataRecordOffset+(PIC::ParticleBuffer::byte*)ParticleData);

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  //check the tracking condition ONLY for partilces that was not tracked previously
  if (DataRecord->TrajectoryTrackingFlag==true) return;
  #endif

  flag=_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_(x,v,spec,ParticleData);

  #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  DataRecord->TrajectoryTrackingFlag=flag;

  if (flag==true) {
    DataRecord->nSampledTrajectoryPoints=0;
    DataRecord->Trajectory.StartingThread=PIC::ThisThread;
    DataRecord->Trajectory.id=SampledTrajectoryCounter;

    RecordTrajectoryPoint(x,v,spec,ParticleData);

    //increment the trajectory counter
    ++threadSampledTrajectoryNumber[spec];
    ++SampledTrajectoryCounter;
  }
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










