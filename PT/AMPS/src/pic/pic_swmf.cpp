//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//the interface between AMPS and SWMF

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "pic.h"

using namespace std;

int PIC::CPLR::SWMF::MagneticFieldOffset=-1;
int PIC::CPLR::SWMF::TotalDataLength=0;


int PIC::CPLR::SWMF::RequestDataBuffer(int offset) {
  MagneticFieldOffset=offset;
  TotalDataLength+=3;

  return TotalDataLength*sizeof(double);
}

void PIC::CPLR::SWMF::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"gmBx\", \"gmBy\", \"gmBz\"");
}

void PIC::CPLR::SWMF::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double B[3]={0.0,0.0,0.0};
  int i,idim;
  char *SamplingBuffer,*CellNodeSamplingBuffer;

  CellNodeSamplingBuffer=CenterNode->GetAssociatedDataBufferPointer()+MagneticFieldOffset;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+MagneticFieldOffset;

    for (idim=0;idim<3;idim++) B[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];
  }

  memcpy(CellNodeSamplingBuffer,B,3*sizeof(double));
}

void PIC::CPLR::SWMF::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;

  for (idim=0;idim<3;idim++) {
    double t;

    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}


void PIC::CPLR::SWMF::init() {
  //request sampling buffer and particle fields
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);

  //print out of the otuput file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
}



void PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(signed int* iComm,signed int* iProc,signed int* nProc) {
  MPI_GLOBAL_COMMUNICATOR=MPI_Comm_f2c(*iComm);
  PIC::InitMPI();

  if (PIC::ThisThread==0) {
    printf("AMPS: MPI Communicatior is imported from SWMF, size=%i\n",PIC::nTotalThreads);
  }
}

void PIC::CPLR::SWMF::ResetCenterPointProcessingFlag() {
  int res=0,thread,i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
          if (cell!=NULL) cell->nodeDescriptor.nodeProcessedFlag=_OFF_AMR_MESH_;
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));
          if (cell!=NULL) cell->nodeDescriptor.nodeProcessedFlag=_OFF_AMR_MESH_;
        }
    }
  }
}


void PIC::CPLR::SWMF::GetCenterPointNumber(int *nCenterPoints) {
  int thread,i,j,k;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();
  *nCenterPoints=0;

  //count the number of the center points
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            (*nCenterPoints)++;
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            (*nCenterPoints)++;
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
          }
        }
    }
  }

}

void PIC::CPLR::SWMF::GetCenterPointCoordinates(double *x) {
  int thread,i,j,k,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            cell->GetX(x+3*cnt);
            cnt++;
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            cell->GetX(x+3*cnt);
            cnt++;
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;
          }
        }
    }
  }
}

void PIC::CPLR::SWMF::RecieveCenterPointData(double *data,int *index) {
  int thread,i,j,k,idim,offset,cnt=0;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataBlockAMR *block;
  PIC::Mesh::cDataCenterNode *cell;

  //init the cell processing flags
  ResetCenterPointProcessingFlag();

  //get coordinated of the center points
  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            offset=TotalDataLength*(index[cnt++]-1);
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;

            for (idim=0;idim<3;idim++) *((double*)(cell->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)))=data[idim+offset];
          }
        }
    }
  }

  for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)
        for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
          cell=block->GetCenterNode(PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k));

          if (cell!=NULL) if (cell->nodeDescriptor.nodeProcessedFlag==_OFF_AMR_MESH_) {
            offset=TotalDataLength*(index[cnt++]-1);
            cell->nodeDescriptor.nodeProcessedFlag=_ON_AMR_MESH_;

            for (idim=0;idim<3;idim++) *((double*)(cell->GetAssociatedDataBufferPointer()+MagneticFieldOffset+idim*sizeof(double)))=data[idim+offset];
          }
        }
    }
  }
}


