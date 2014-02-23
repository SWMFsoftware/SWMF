/*
 * pic_debugger.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: vtenishe
 */
//$Id$
//contains functions used for debugging AMPS

#include "pic.h"

//InfiniteLoop==false ==> no problem found; InfiniteLoop==true => the actual number of particles does not consider the that in teh particle buffer
bool PIC::Debugger::InfiniteLoop(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  bool res=false;
  long int ptr;
  static long int nAllCountedParticles=0;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;
  if (startNode==PIC::Mesh::mesh.rootTree) nAllCountedParticles=0;


  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InfiniteLoop(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block

    for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      ptr=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      while (ptr!=-1) {
        ptr=PIC::ParticleBuffer::GetNext(ptr);

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          FindDoubleReferencedParticle();

          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //collect the information from all processors;
    int Flag,FlagTable[PIC::nTotalThreads];

    Flag=(nAllCountedParticles==PIC::ParticleBuffer::NAllPart) ? true : false;
    MPI_Gather(FlagTable,1,MPI_INT,&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      for (int thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==false) {
        printf("$PREFIX: the actual number of particles does not coniside with that in the particle buffer (Thread=%i,line=%i,file=%s)\n",thread,__LINE__,__FILE__);
        res=true;
      }
    }

    MPI_Bcast(&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    nAllCountedParticles=0;
  }

  return res;
}



void PIC::Debugger::FindDoubleReferencedParticle(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  long int ptr,ptrNext;
  static long int nAllCountedParticles=0;

  static char* ParticleAllocationTable=NULL;
  static long int* ptrPrevTable=NULL;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;

  //create table of allocated/not-allocated flag table
  if (startNode==PIC::Mesh::mesh.rootTree) {
    long int cnt=0;

    nAllCountedParticles=0;
    ParticleAllocationTable=new char [PIC::ParticleBuffer::MaxNPart];
    ptrPrevTable=new long int [PIC::ParticleBuffer::MaxNPart];

    //the definition of the bits of ParticleAllocationTable[]
    //ParticleAllocationTable[] & 1: == 0 --> particle is not allocated; == 1 --> the particle is allocated
    //ParticleAllocationTable[] & 2: == 0 --> the particle is not lonked yet; == 2 --> the particle is already linked

    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) ptrPrevTable[ptr]=-1,ParticleAllocationTable[ptr]=1; //particle is allocated

    if (PIC::ParticleBuffer::FirstPBufferParticle!=-1) ParticleAllocationTable[PIC::ParticleBuffer::FirstPBufferParticle]=2; //particle is links and not allocated

    for (cnt=0,ptr=PIC::ParticleBuffer::FirstPBufferParticle;ptr!=-1;ptr=ptrNext) {
      ptrNext=PIC::ParticleBuffer::GetNext(ptr);

      if (ptrNext!=-1) {
        if ((ParticleAllocationTable[ptrNext]&2)==2) {
          printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);

          exit(__LINE__,__FILE__,"Error: an un-allocated particle is double-referenced");
        }
        else ParticleAllocationTable[ptrNext]=2,ptrPrevTable[ptrNext]=ptr;
      }

      if (++cnt>PIC::ParticleBuffer::MaxNPart) {
        exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
      }

      ptr=ptrNext;
    }
  }

  //check the particle lists
  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) FindDoubleReferencedParticle(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block

    for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (int npass=0;npass<2;npass++) {
      ptr=(npass==0) ? startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)] : startNode->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      if (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        if ((ParticleAllocationTable[ptr]&2)==2) {
          printf("Error: the first particle in the list is referenced: %ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptr],ptr);
          exit(__LINE__,__FILE__,"Error: the first particle in the list is referenced");
        }
        else ParticleAllocationTable[ptr]|=2;

        if (PIC::ParticleBuffer::GetPrev(ptr)>=0) {
          exit(__LINE__,__FILE__,"Error: the first particle in the list has non-negarive value of GetPrev");
        }
      }


      while (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        ptrNext=PIC::ParticleBuffer::GetNext(ptr);

        if (ptrNext!=-1) {
          if ((ParticleAllocationTable[ptrNext]&1)==0) {
            exit(__LINE__,__FILE__,"Error: reference to an un-allocated particles is found in the particle list");
          }

          if ((ParticleAllocationTable[ptrNext]&2)==2) {
            printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);
            exit(__LINE__,__FILE__,"Error: have found double-referenced particle in the list");
          }

          if (PIC::ParticleBuffer::GetPrev(ptrNext)!=ptr) {
            exit(__LINE__,__FILE__,"Error: PIC::ParticleBuffer::GetPrev(ptrNext) != ptr");
          }

          ParticleAllocationTable[ptrNext]|=2;
          ptrPrevTable[ptrNext]=ptr;
        }

        ptr=ptrNext;

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) {
      if ((ParticleAllocationTable[ptr]&2)==0) {
        exit(__LINE__,__FILE__,"Error: found particles that is not referenced at all");
      }
    }

    delete [] ParticleAllocationTable;
    delete [] ptrPrevTable;
  }
}



