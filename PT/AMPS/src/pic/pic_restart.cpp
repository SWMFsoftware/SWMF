//$Id$
//the restart procedures of AMPS

#include "pic.h"

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
}

void PIC::Restart::SaveParticleDataBlock(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,CMPI_channel* pipe,FILE* fRestart) {

  //save the data
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->block!=NULL) {
      int i,j,k,LocalCellNumber;
      PIC::Mesh::cDataCenterNode *cell;
      char* SamplingData;

      //block is allocated -> march through the cells and save them into the restart file
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        cell=node->block->GetCenterNode(LocalCellNumber);

        SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset;

        if (PIC::ThisThread==0) {
          //save the data into the file
          if (node->Thread!=0) SamplingData=pipe->recvPointer<char>(PIC::Mesh::sampleSetDataLength,node->Thread);
          fwrite(SamplingData,sizeof(char),PIC::Mesh::sampleSetDataLength,fRestart);
        }
        else {
          pipe->send(SamplingData,PIC::Mesh::sampleSetDataLength);
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (node->downNode[nDownNode]!=NULL) SaveSamplingDataBlock(node->downNode[nDownNode],pipe,fRestart);
  }
}

















