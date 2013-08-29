//$Id$
//calculate integrated parameters that describes the "state" of the entire simulation

#include "pic.h"

//get the check sum of all particles
void PIC::RunTimeSystemState::GetParticleFieldCheckSum(char *msg) {
  CRC32 CheckSum;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  int i,j,k;
  long int ptr;

  //calcualte the check sum of the particles with in the domain
  while (node!=NULL) {
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++)  {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            CheckSum.add<char>((char*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);

            ptr=PIC::ParticleBuffer::GetNext(ParticleData);
          }

        }
      }
    }

    node=node->nextNodeThisThread;
  }

  //calcualte the check sum of particles in the 'ghost'cells
  for (int To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
    for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];node!=NULL;node=node->nextNodeThisThread) {
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
          for (i=0;i<_BLOCK_CELLS_X_;i++)  {
            ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              CheckSum.add<char>((char*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }

          }
        }
      }

    }
  }

  unsigned long CheckSumBuffer[PIC::nTotalThreads],t;

  t=CheckSum.checksum();
  MPI_Gather(&t,1,MPI_LONG,CheckSumBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
     printf("AMPS: Particle Field Check Sum:");
     if (msg!=NULL) printf("(message=\"%s\")",msg);
     for (int thread=0;thread<PIC::nTotalThreads;thread++) printf("  0x%lx",CheckSumBuffer[thread]);
     printf("\n");
  }
}

void PIC::RunTimeSystemState::GetParticleFieldCheckSum(long int nline,char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetParticleFieldCheckSum(msg);
}

//compare the domain decomposition: calcualte the chech sum of all block's TempID belong to the currect processor
void PIC::RunTimeSystemState::GetDomainDecompositionCheckSum(char *msg) {
  CRC32 CheckSum;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  while (node!=NULL) {
    CheckSum.add(node->Temp_ID);
    node=node->nextNodeThisThread;
  }

  unsigned long CheckSumBuffer[PIC::nTotalThreads],t;

  t=CheckSum.checksum();
  MPI_Gather(&t,1,MPI_LONG,CheckSumBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    printf("AMPS: Domain Decomposition Check Sum:");
    if (msg!=NULL) printf("(message=\"%s\") ",msg);

    for (int thread=0;thread<PIC::nTotalThreads;thread++) printf("  0x%lx",CheckSumBuffer[thread]);
    printf("\n");
  }
}

void PIC::RunTimeSystemState::GetDomainDecompositionCheckSum(long int nline,char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetDomainDecompositionCheckSum(msg);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(FILE* fout,char *msg) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  int i,j,k,idim;
  long int ptr;

  int s;
  double *v;
  double StatWeight,TotalStatWeight[PIC::nTotalSpecies],MeanSpeed[PIC::nTotalSpecies],MeanVelocity[3*PIC::nTotalSpecies];

  for (s=0;s<PIC::nTotalSpecies;s++) {
    MeanSpeed[s]=0.0;
    TotalStatWeight[s]=0.0;

    for (idim=0;idim<3;idim++) MeanVelocity[idim+3*s]=0.0;
  }

  //sample the mean values
  while (node!=NULL) {
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++)  {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            s=PIC::ParticleBuffer::GetI(ParticleData);

            StatWeight=node->block->GetLocalParticleWeight(s)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
            v=PIC::ParticleBuffer::GetV(ParticleData);

            TotalStatWeight[s]+=StatWeight;
            MeanSpeed[s]+=StatWeight*sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            for (idim=0;idim<3;idim++) MeanVelocity[idim+3*s]+=StatWeight*v[idim];

            ptr=PIC::ParticleBuffer::GetNext(ParticleData);
          }

        }
      }
    }

    node=node->nextNodeThisThread;
  }

  double TempBuffer[3*PIC::nTotalSpecies];

  MPI_Reduce(TotalStatWeight,TempBuffer,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(TotalStatWeight,TempBuffer,PIC::nTotalSpecies*sizeof(double));

  MPI_Reduce(MeanSpeed,TempBuffer,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(MeanSpeed,TempBuffer,PIC::nTotalSpecies*sizeof(double));

  MPI_Reduce(MeanVelocity,TempBuffer,3*PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(MeanVelocity,TempBuffer,3*PIC::nTotalSpecies*sizeof(double));

  if (PIC::ThisThread==0) {
    fprintf(fout,"AMPS: Averaged particles microscopic aprameters:");
    if (msg!=NULL) fprintf(fout,"(message=\"%s\") ",msg);
    fprintf(fout,"\nAMPS: spec\t<Speed>\t\t<v[0]>\t\t<v[1]>\t\t<v[2]>\n");

    for (s=0;s<PIC::nTotalSpecies;s++) {
      if (TotalStatWeight[s]>0.0) {
        MeanSpeed[s]/=TotalStatWeight[s];

        for (idim=0;idim<3;idim++) MeanVelocity[3*s+idim]/=TotalStatWeight[s];
      }

      fprintf(fout,"AMPS: %i\t%e\t%e\t%e\t%e\t\n",s,MeanSpeed[s],MeanVelocity[0+3*s],MeanVelocity[1+3*s],MeanVelocity[2+3*s]);
    }

    fprintf(fout,"\n");
  }

  fflush(fout);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(FILE* fout,long int nline,char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetMeanParticleMicroscopicParameters(fout,msg);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(const char *fname) {
  FILE *fout=NULL;

  if (PIC::ThisThread==0) fout=fopen(fname,"w");
  GetMeanParticleMicroscopicParameters(fout,NULL);
  if (PIC::ThisThread==0) fclose(fout);
}

