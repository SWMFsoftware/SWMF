//$Id$
//the restart procedures of AMPS

#include "pic.h"

int PIC::Restart::ParticleRestartAutosaveIterationInterval=1;
char PIC::Restart::SamplingDataRestartFileName[_MAX_STRING_LENGTH_PIC_]="SampleData.restart";

char PIC::Restart::saveParticleDataRestartFileName[_MAX_STRING_LENGTH_PIC_]="ParticleData.restart";
char PIC::Restart::recoverParticleDataRestartFileName[_MAX_STRING_LENGTH_PIC_]="ParticleData.restart";
bool PIC::Restart::ParticleDataRestartFileOverwriteMode=true;


//-------------------------------------- Save/Load Sampling Data Restart File ---------------------------------------------------------
//save sampling data
void PIC::Restart::SaveSamplingData(const char* fname) {
  FILE *fRestart=NULL;
  int BlockDataSize=0,mpiBufferSize=0;


  //determine the size of the data vector associated with a block
  BlockDataSize=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*PIC::Mesh::sampleSetDataLength;
  mpiBufferSize=max((int)10.0*BlockDataSize,(int)10000000);

  //init the MPI channel and open the restart file
  CMPI_channel pipe(mpiBufferSize);

  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.openRecvAll();

    fRestart=fopen(fname,"w");

    fwrite(&PIC::LastSampleLength,sizeof(PIC::LastSampleLength),1,fRestart);
    fwrite(&PIC::DataOutputFileNumber,sizeof(PIC::DataOutputFileNumber),1,fRestart);
  }
  else pipe.openSend(0);

  //save the restart information
  SaveSamplingDataBlock(PIC::Mesh::mesh.rootTree,&pipe,fRestart);

  //close the MPI channel and the restart file
  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.closeRecvAll();
    fclose(fRestart);
  }
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void PIC::Restart::SaveSamplingDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,CMPI_channel* pipe,FILE* fRestart) {
  int nAllocatedCells,i,j,k;

  //save the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((node->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      int LocalCellNumber;
      PIC::Mesh::cDataCenterNode *cell;
      char* SamplingData;

      //Calculate the number of the allocated cells in the block and send it to the root processor
      nAllocatedCells=0;

      if (node->Thread==PIC::ThisThread) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
          cell=node->block->GetCenterNode(LocalCellNumber);

          if (cell!=NULL) nAllocatedCells++;
        }

        if (node->Thread!=0) pipe->send(nAllocatedCells);
      }
      else if (PIC::ThisThread==0) {
        pipe->recv(&nAllocatedCells,1,node->Thread);
      }

      if (PIC::ThisThread==0) fwrite(&nAllocatedCells,sizeof(int),1,fRestart);

      //save the sampling data
      if (node->Thread==PIC::ThisThread) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
          cell=node->block->GetCenterNode(LocalCellNumber);

          if (cell==NULL) continue;

          SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

          if (PIC::ThisThread==0) {
            fwrite(SamplingData,sizeof(char),PIC::Mesh::sampleSetDataLength,fRestart);
          }
          else {
            //send the sampling information to the root processor
            pipe->send(SamplingData,PIC::Mesh::sampleSetDataLength);
          }
        }
      }
      else if (PIC::ThisThread==0) {
        //recieve the sampling information from other processor and save it into a file
        for (i=0;i<nAllocatedCells;i++) {
          SamplingData=pipe->recvPointer<char>(PIC::Mesh::sampleSetDataLength,node->Thread);
          fwrite(SamplingData,sizeof(char),PIC::Mesh::sampleSetDataLength,fRestart);
        }
      }

    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) SaveSamplingDataBlock(node->downNode[nDownNode],pipe,fRestart);
  }
}

void PIC::Restart::ReadSamplingData(const char* fname) {
  FILE *fRestart=NULL;

  fRestart=fopen(fname,"r");

  if (fRestart==NULL) {
    char msg[200];

    sprintf(msg,"Error: restart file %s is not found",fname);
    exit(__LINE__,__FILE__,msg);
  }

  fread(&PIC::LastSampleLength,sizeof(PIC::LastSampleLength),1,fRestart);
  fread(&PIC::DataOutputFileNumber,sizeof(PIC::DataOutputFileNumber),1,fRestart);

  ReadSamplingDataBlock(PIC::Mesh::mesh.rootTree,fRestart);
  fclose(fRestart);

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void PIC::Restart::ReadSamplingDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,FILE* fRestart) {
  int nAllocatedCells;

  //read the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    fread(&nAllocatedCells,sizeof(int),1,fRestart);

    if (node->block!=NULL) {
      //read the data for this block

      int i,j,k,LocalCellNumber;
      PIC::Mesh::cDataCenterNode *cell;
      char* SamplingData;

      //block is allocated -> march through the cells and save them into the restart file
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        cell=node->block->GetCenterNode(LocalCellNumber);

        if (cell!=NULL) {
          SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;
          fread(SamplingData,sizeof(char),PIC::Mesh::sampleSetDataLength,fRestart);
        }
      }
    }
    else {
      //skip the data for this block
      fseek(fRestart,PIC::Mesh::sampleSetDataLength*nAllocatedCells,SEEK_CUR);
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) ReadSamplingDataBlock(node->downNode[nDownNode],fRestart);
  }
}

//-------------------------------------- Save/Load Particle Data Restart File ---------------------------------------------------------
void PIC::Restart::SaveParticleData(const char* fname) {
  FILE *fRestart=NULL;

  //init the MPI channel and open the restart file
  CMPI_channel pipe(10000000);

  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.openRecvAll();
    fRestart=fopen(fname,"w");
  }
  else pipe.openSend(0);

  //save the restart information
  SaveParticleDataBlock(PIC::Mesh::mesh.rootTree,&pipe,fRestart);

  //close the MPI channel and the restart file
  if (PIC::Mesh::mesh.ThisThread==0) {
    pipe.closeRecvAll();
    fclose(fRestart);
  }
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  GetParticleDataCheckSum();
}

void PIC::Restart::SaveParticleDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,CMPI_channel* pipe,FILE* fRestart) {
  //save the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((node->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      int i,j,k,LocalCellNumber;
      PIC::Mesh::cDataCenterNode *cell;
      char* SamplingData;
      long int ptr;

      //determine the number of the model particle in the each cell of the block
      int nTotalParticleNumber=0;
      int ParticleNumberTable[_BLOCK_CELLS_X_][_BLOCK_CELLS_Y_][_BLOCK_CELLS_Z_];
      long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

      if (node->Thread==PIC::ThisThread) {
        memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          ParticleNumberTable[i][j][k]=0;
          ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleNumberTable[i][j][k]++;
            nTotalParticleNumber++;
            ptr=PIC::ParticleBuffer::GetNext(ptr);
          }
        }

        if (node->Thread!=0) {
          pipe->send(nTotalParticleNumber);
          if (nTotalParticleNumber!=0) pipe->send(&ParticleNumberTable[0][0][0],_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
        }
      }
      else if (PIC::ThisThread==0) {
        pipe->recv(&nTotalParticleNumber,1,node->Thread);
        if (nTotalParticleNumber!=0) pipe->recv(&ParticleNumberTable[0][0][0],_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,node->Thread);
      }

      //save the particle number into the restart file
      if (PIC::ThisThread==0) {
        fwrite(&nTotalParticleNumber,sizeof(int),1,fRestart);
        if (nTotalParticleNumber!=0) fwrite(&ParticleNumberTable[0][0][0],sizeof(int),_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,fRestart);
      }

      //save the particle data into the restart file
      //IMPORTANT: save the partilce data in the reverse order so, when thay are read back from the restart file they are in the seme order as
      //in the memory before saving --> the checksum of the particle data is the same
      if (nTotalParticleNumber!=0) for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
        int n;

        for (n=0;n<ParticleNumberTable[i][j][k];n++) {
          if (node->Thread==PIC::ThisThread) {
            if (n==0) {
              ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
              for (int t=0;t<ParticleNumberTable[i][j][k]-1;t++) ptr=PIC::ParticleBuffer::GetNext(ptr);
            }
            else ptr=PIC::ParticleBuffer::GetPrev(ptr);

           memcpy(tempParticleData,PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

           if (node->Thread!=0) pipe->send(tempParticleData,PIC::ParticleBuffer::ParticleDataLength);
          }
          else {
            pipe->recv(tempParticleData,PIC::ParticleBuffer::ParticleDataLength,node->Thread);
          }

          if (PIC::ThisThread==0) fwrite(tempParticleData,sizeof(char),PIC::ParticleBuffer::ParticleDataLength,fRestart);
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) SaveParticleDataBlock(node->downNode[nDownNode],pipe,fRestart);
  }
}


void PIC::Restart::ReadParticleDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,FILE* fRestart) {
  //read the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    int nTotalParticleNumber=0;
    int ParticleNumberTable[_BLOCK_CELLS_X_][_BLOCK_CELLS_Y_][_BLOCK_CELLS_Z_];
    long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];

    fread(&nTotalParticleNumber,sizeof(int),1,fRestart);

    if (nTotalParticleNumber!=0) {
      if (node->Thread==PIC::ThisThread) { ///(node->block!=NULL) {
        //read the data for this block
        fread(&ParticleNumberTable[0][0][0],sizeof(int),_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,fRestart);
        memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

        int i,j,k,np;
        char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];
        long int ptr;

        //block is allocated -> march through the cells and save them into the restart file
        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (np=0;np<ParticleNumberTable[i][j][k];np++) {
            fread(tempParticleData,sizeof(char),PIC::ParticleBuffer::ParticleDataLength,fRestart);
            ptr=PIC::ParticleBuffer::GetNewParticle(FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);

            PIC::ParticleBuffer::CloneParticle((PIC::ParticleBuffer::byte*) PIC::ParticleBuffer::GetParticleDataPointer(ptr),(PIC::ParticleBuffer::byte*) tempParticleData);

            //in case when particle tracking is on -> apply the particle tracking conditions if needed
            if (_PIC_PARTICLE_TRACKER_MODE_ ==_PIC_MODE_ON_) {
              PIC::ParticleBuffer::byte *ParticleData;
              PIC::ParticleTracker::cParticleData *DataRecord;

              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              DataRecord=(PIC::ParticleTracker::cParticleData*)(PIC::ParticleTracker::ParticleDataRecordOffset+(PIC::ParticleBuffer::byte*)ParticleData);
              DataRecord->TrajectoryTrackingFlag=false;


              if (_PIC_PARTICLE_TRACKER__RESTART_LOADED_PARTICLES__APPLY_TRACKING_CONDITION_MODE_==_PIC_MODE_ON_) {
                //apply the particle tracking condition if needed
                double x[3],v[3];
                int spec;

                PIC::ParticleBuffer::GetX(x,ParticleData);
                PIC::ParticleBuffer::GetV(v,ParticleData);
                spec=PIC::ParticleBuffer::GetI(ParticleData);

                PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,(void *)ParticleData);
              }

            }
          }
        }

        memcpy(node->block->FirstCellParticleTable,FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));
      }
      else {
        //skip the data for this block
        fseek(fRestart,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(int)+
            nTotalParticleNumber*PIC::ParticleBuffer::ParticleDataLength*sizeof(char),SEEK_CUR);
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) ReadParticleDataBlock(node->downNode[nDownNode],fRestart);
  }
}


void PIC::Restart::ReadParticleData(const char* fname) {
  FILE *fRestart=NULL;

  fRestart=fopen(fname,"r");

  if (fRestart==NULL) {
    char msg[200];

    sprintf(msg,"Error: restart file %s is not found",fname);
    exit(__LINE__,__FILE__,msg);
  }

  ReadParticleDataBlock(PIC::Mesh::mesh.rootTree,fRestart);
  fclose(fRestart);

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  GetParticleDataCheckSum();
}

//calculate the check sum of the particle data
//-------------------------------------- Calculate the checlsum of the particle buffer ------------------------------------------------
unsigned long PIC::Restart::GetParticleDataCheckSum() {
  CRC32 CheckSum;

  //the thread number that was processed last: the checkSum object is sent from PrevNodeThread directly to node->Thread
  //at the begining of the calculation CheckSum is located on the root thread
  int PrevNodeThread=0;

  GetParticleDataBlockCheckSum(PIC::Mesh::mesh.rootTree,&CheckSum,PrevNodeThread);
  MPI_Bcast(&CheckSum,sizeof(CRC32),MPI_CHAR,PrevNodeThread,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    printf("$PREFIX:The particle data CRC32 checksum=0x%lx (%i@%s):\n",CheckSum.checksum(),__LINE__,__FILE__);
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  return CheckSum.checksum();
}


void PIC::Restart::GetParticleDataBlockCheckSum(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,CRC32* CheckSum,int &PrevNodeThread) {
  //save the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    //recieve the CheckSum object

    if ((node->Thread==PIC::ThisThread)||(PIC::ThisThread==PrevNodeThread)) {
      if (node->Thread!=PrevNodeThread) {
        if (PIC::ThisThread==PrevNodeThread) {
          //send CheckSum to node->Thread
          MPI_Send(CheckSum,sizeof(CRC32),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR);
        }
        else {
          //recieve CheckSum from PrevNodeThread
          MPI_Status status;
          MPI_Recv(CheckSum,sizeof(CRC32),MPI_CHAR,PrevNodeThread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        }
      }
    }

    PrevNodeThread=node->Thread;

    //calculate the chekc sum
    if (node->Thread==PIC::ThisThread) {
      int i,j,k,LocalCellNumber;
      long int ptr;
      long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_];
      char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

      memcpy(FirstCellParticleTable,node->block->FirstCellParticleTable,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_*sizeof(long int));

      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        ptr=FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        while (ptr!=-1) {
          memcpy(tempParticleData,PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

          CheckSum->add(tempParticleData,sizeof(unsigned char));
          CheckSum->add(tempParticleData+sizeof(unsigned char)+2*sizeof(long int),
              PIC::ParticleBuffer::ParticleDataLength-sizeof(unsigned char)-2*sizeof(long int));

          ptr=PIC::ParticleBuffer::GetNext(ptr);
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) GetParticleDataBlockCheckSum(node->downNode[nDownNode],CheckSum,PrevNodeThread);
  }
}








