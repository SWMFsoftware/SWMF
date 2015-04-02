
/*
 * pic_batsrus.cpp
 *
 *  Created on: Apr 1, 2015
 *      Author: vtenishe
 */

//$Id$
//read of the BATSRUS output file
#include "pic.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

extern "C"
{

  void batsrus2amps_get_nvar_(int *nVar);
  void batsrus2amps_domain_limits_(double *xmin,double *xmax);
  void batsrus2amps_openfile_(char *FileName,int *FileNameLength);
  void batsrus2amps_read_file_header_(char *FileName,int *FileNameLength);
  void batsrus2amps_closefile_();
  void batsrus2amps_get_data_point_(double* x,double *res,int *FoundFlag);
  void batsrus2amps_get_namevardata_(char *varlist,int *varlistlength);
  void batsrus2amps_set_mpi_parameters_(int* ThisThread,int *nTotalThreads,int *Communicator);
}


void PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static double *StateLocal=NULL,*State=NULL;
  static int nVar;


  if (startNode==PIC::Mesh::mesh.rootTree) {
    //open data file
    char fullname[_MAX_STRING_LENGTH_PIC_],varlist[_MAX_STRING_LENGTH_PIC_];
    int length=_MAX_STRING_LENGTH_PIC_;

    //sprintf(fullname,"%s/%s",PIC::CPLR::DATAFILE::path,fname);



    sprintf(fullname,"%s",fname);

    for (length=0;length<_MAX_STRING_LENGTH_PIC_;length++) if (fullname[length]==0) {
//      length-=1;
      break;
    }

//    length++;

    int iComm=MPI_Comm_c2f(MPI_GLOBAL_COMMUNICATOR);

    batsrus2amps_set_mpi_parameters_(&PIC::ThisThread,&PIC::nTotalThreads,&iComm);

//    batsrus2amps_read_file_header_(fullname,&length);
//    batsrus2amps_get_namevardata_(varlist,&length);
//    batsrus2amps_get_nvar_(&nVar);

    batsrus2amps_openfile_(fullname,&length);

    batsrus2amps_get_namevardata_(varlist,&length);
    batsrus2amps_get_nvar_(&nVar);
    State=new double [nVar+1];
    StateLocal=new double [nVar+1];

    int IsFound;
    double x[3]={0.0,0.0,0.0};

    batsrus2amps_get_data_point_(x,StateLocal,&IsFound);
  }


  //perform the interpolation loop
  int i,j,k,nd,idim;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double x[3],T,n,p;

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
      x[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
      x[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);

      //recover the data from the BATSRUS data file
      int IsFound;

      batsrus2amps_get_data_point_(x,StateLocal,&IsFound);

      if (IsFound==true) {
        MPI_Allreduce(StateLocal,State,nVar+1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        //devide the interpolated values by the "total weight"
        for (int n=1;n<nVar+1;n++) State[n]/=State[0];
      }
      else for (int n=0;n<nVar+1;n++) State[n]=0.0;

      //locate the cell
      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
      offset=CenterNode->GetAssociatedDataBufferPointer();

      //save the interpolated values
      for (idim=0;idim<3;idim++) {
        *(idim+(double*)(offset+MagneticFieldOffset))=State[0];
        *(idim+(double*)(offset+PlasmaBulkVelocityOffset))=State[0];
      }

      p=State[0];
      n=State[0];
      T=(n>0.0) ? p/(n*Kbol) : 0.0;

      *((double*)(offset+PlasmaPressureOffset))=p;
      *((double*)(offset+PlasmaNumberDensityOffset))=n;
      *((double*)(offset+PlasmaTemperatureOffset))=T;
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(fname,startNode->downNode[nDownNode]);
  }


  if (startNode==PIC::Mesh::mesh.rootTree) {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    delete [] State;
    delete [] StateLocal;

    batsrus2amps_closefile_();
  }
}




