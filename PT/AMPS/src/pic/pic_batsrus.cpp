
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

//the definition of the fortran wrapper functinos from pic_batsrus_wrapper.f90
extern "C"
{
  //manipulate the file
  void batsrus2amps_read_file_header_(char *FileName,int *FileNameLength);
  void batsrus2amps_openfile_(char *FileName,int *FileNameLength);
  void batsrus2amps_closefile_();

  //get the side of the domain
  void batsrus2amps_domain_limits_(double *xmin,double *xmax);

  //set MPI data
  void batsrus2amps_set_mpi_parameters_(int* ThisThread,int *nTotalThreads,int *Communicator);

  //get number of the variable, the list of the variables and their uinits
  void batsrus2amps_get_nvar_(int *nVar);
  void batsrus2amps_get_namevardata_(char *varlist,int *varlistlength);
  void batsrus2amps_get_nameunitdata_(char *unitlist,int *unitlistlength);

  //get the value of the variables
  void batsrus2amps_get_data_point_(double* x,double *res,int *FoundFlag);
}

//definition of the variabled from BATSRUS namespace
double PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::PlasmaSpeciesAtomicMass=1.0;
double PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::UnitLength=1.0;
char PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::filename[_MAX_STRING_LENGTH_PIC_]="";
bool PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::InitFlag=false;

//init the BATSRUS data file AMR information
void PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::Init(const char *fname) {
  int length;
  
  InitFlag=true;

  //get the full file name
  sprintf(filename,"%s/%s",PIC::CPLR::DATAFILE::path,fname);

  //export MPI parameters into BATL
  int iComm=MPI_Comm_c2f(MPI_GLOBAL_COMMUNICATOR);
  batsrus2amps_set_mpi_parameters_(&PIC::ThisThread,&PIC::nTotalThreads,&iComm);
  
  //init BATSRUS' fiile reader
  for (length=0;length<_MAX_STRING_LENGTH_PIC_;length++) if (filename[length]==0) break;
  batsrus2amps_read_file_header_(filename,&length);
} 

void  PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::GetDomainLimits(double *xmin,double *xmax) {
  int idim;

  batsrus2amps_domain_limits_(xmin,xmax);
  for (idim=0;idim<3;idim++) xmin[idim]*=UnitLength,xmax[idim]*=UnitLength;
}

//read BATSRUS' .idl file
void PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::LoadDataFile(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static double *StateLocal=NULL,*State=NULL,*PhysicalVariableUnitConversionTable=NULL;
  static int nVar;

  //the offsets of the physical variables in the .idl file
  static int rhoBATSRUS2AMPS=-2;
  static int mxBATSRUS2AMPS=-2,myBATSRUS2AMPS=-2,mzBATSRUS2AMPS=-2;
  static int uxBATSRUS2AMPS=-2,uyBATSRUS2AMPS=-2,uzBATSRUS2AMPS=-2;
  static int bxBATSRUS2AMPS=-2,byBATSRUS2AMPS=-2,bzBATSRUS2AMPS=-2;
  static int pBATSRUS2AMPS=-2;

  if (InitFlag==false) exit(__LINE__,__FILE__,"Error: the reader needs to be initialized first! Call PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::InitAMR before PIC::CPLR::DATAFILE::BATSRUS::OUTPUT::LoadDataFile");

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //open data file
    char NameVar[_MAX_STRING_LENGTH_PIC_];
    int length=_MAX_STRING_LENGTH_PIC_;

    //open the data file
    for (length=0;length<_MAX_STRING_LENGTH_PIC_;length++) if (filename[length]==0) break;
    batsrus2amps_openfile_(filename,&length);

    //read the variable nabmer and string
    length=_MAX_STRING_LENGTH_PIC_;
    batsrus2amps_get_namevardata_(NameVar,&length);
    batsrus2amps_get_nvar_(&nVar);
    State=new double [nVar+1];
    StateLocal=new double [nVar+1];
    PhysicalVariableUnitConversionTable=new double [nVar+1];

    for (int i=0;i<nVar+1;i++) PhysicalVariableUnitConversionTable[i]=0.0;

    //parse the variable line
    char vname[200];
    int i0=0,i1=0;

    while ((NameVar[i0]!=0)&&(NameVar[i0]==' ')) i0++;

    for (int n=0;n<nVar;n++) {
      i1=i0;
      while (NameVar[i1]!=' ') {
        vname[i1-i0]=tolower(NameVar[i1]);
        i1++;
      }

      vname[i1-i0]=0;

      if (strcmp(vname,"rho")==0) rhoBATSRUS2AMPS=n;

      if (strcmp(vname,"mx")==0) mxBATSRUS2AMPS=n;
      if (strcmp(vname,"my")==0) myBATSRUS2AMPS=n;
      if (strcmp(vname,"mz")==0) mzBATSRUS2AMPS=n;

      if (strcmp(vname,"ux")==0) uxBATSRUS2AMPS=n;
      if (strcmp(vname,"uy")==0) uyBATSRUS2AMPS=n;
      if (strcmp(vname,"uz")==0) uzBATSRUS2AMPS=n;

      if (strcmp(vname,"bx")==0) bxBATSRUS2AMPS=n;
      if (strcmp(vname,"by")==0) byBATSRUS2AMPS=n;
      if (strcmp(vname,"bz")==0) bzBATSRUS2AMPS=n;

      if (strcmp(vname,"p")==0) pBATSRUS2AMPS=n;

      i0=i1;
      while ((NameVar[i0]!=0)&&(NameVar[i0]==' ')) i0++;
    }

    //check whether the state vector containes all nessesary physical quantaties
    if (rhoBATSRUS2AMPS==-1) exit(__LINE__,__FILE__,"Error: rho is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");

    if ((mxBATSRUS2AMPS==-1)&&(uxBATSRUS2AMPS==-1)) exit(__LINE__,__FILE__,"Error: Mx or ux is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");
    if ((myBATSRUS2AMPS==-1)&&(uxBATSRUS2AMPS==-1)) exit(__LINE__,__FILE__,"Error: My or uy is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");
    if ((mzBATSRUS2AMPS==-1)&&(uxBATSRUS2AMPS==-1)) exit(__LINE__,__FILE__,"Error: Mz or uz is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");

    if (bxBATSRUS2AMPS==-1) exit(__LINE__,__FILE__,"Error: Bx is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");
    if (byBATSRUS2AMPS==-1) exit(__LINE__,__FILE__,"Error: By is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");
    if (bzBATSRUS2AMPS==-1) exit(__LINE__,__FILE__,"Error: Bz is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");

    if (pBATSRUS2AMPS==-1) exit(__LINE__,__FILE__,"Error: p is not present in the BARSRUS .idl file. Please add this variable to the .idl file.");



    //parse the unit line
    char UnitVar[_MAX_STRING_LENGTH_PIC_],uname[200];

    length=_MAX_STRING_LENGTH_PIC_;
    batsrus2amps_get_nameunitdata_(UnitVar,&length);
    i0=0,i1=0;

    while ((UnitVar[i0]!=0)&&(UnitVar[i0]==' ')) i0++;

    for (int n=0;n<nVar;n++) {
      i1=i0;
      while (UnitVar[i1]!=' ') {
        uname[i1-i0]=tolower(UnitVar[i1]);
        i1++;
      }

      uname[i1-i0]=0;

      if (strcmp(uname,"r")==0) --n; //the special coordinate is not a part of the vector that is returned by the interpolation routine

      if (strcmp(uname,"mp/cc")==0) PhysicalVariableUnitConversionTable[n+1]=1.0E6/PlasmaSpeciesAtomicMass;
      if (strcmp(uname,"km/s")==0) PhysicalVariableUnitConversionTable[n+1]=1.0E3;
      if (strcmp(uname,"nt")==0) PhysicalVariableUnitConversionTable[n+1]=1.0E-9;
      if (strcmp(uname,"npa")==0) PhysicalVariableUnitConversionTable[n+1]=1.0E-9;

      i0=i1;
      while ((UnitVar[i0]!=0)&&(UnitVar[i0]==' ')) i0++;
    }


    //the first element in the state vector is the "weight" -> adjust the offsets
    rhoBATSRUS2AMPS++;
    mxBATSRUS2AMPS++,myBATSRUS2AMPS++,mzBATSRUS2AMPS++;
    uxBATSRUS2AMPS++,uyBATSRUS2AMPS++,uzBATSRUS2AMPS++;
    bxBATSRUS2AMPS++,byBATSRUS2AMPS++,bzBATSRUS2AMPS++;
    pBATSRUS2AMPS++;

    //check whether the unit conversion factors are defined for all physical variables
    if ((PhysicalVariableUnitConversionTable[rhoBATSRUS2AMPS]==0.0) || \
        (PhysicalVariableUnitConversionTable[bxBATSRUS2AMPS]==0.0) || (PhysicalVariableUnitConversionTable[byBATSRUS2AMPS]==0.0)  || (PhysicalVariableUnitConversionTable[bzBATSRUS2AMPS]==0.0) || \
        (PhysicalVariableUnitConversionTable[pBATSRUS2AMPS]==0.0)) {
      exit(__LINE__,__FILE__,"Error: the physical variable unit conversion factor is not defined");
    }

    if ( ((mxBATSRUS2AMPS>=0) && ((PhysicalVariableUnitConversionTable[mxBATSRUS2AMPS]==0.0) || (PhysicalVariableUnitConversionTable[myBATSRUS2AMPS]==0.0)  || (PhysicalVariableUnitConversionTable[mzBATSRUS2AMPS]==0.0))) || \
       ((uxBATSRUS2AMPS>=0) && ((PhysicalVariableUnitConversionTable[uxBATSRUS2AMPS]==0.0) || (PhysicalVariableUnitConversionTable[uyBATSRUS2AMPS]==0.0)  || (PhysicalVariableUnitConversionTable[uzBATSRUS2AMPS]==0.0))) ) {
      exit(__LINE__,__FILE__,"Error: the physical variable unit conversion factor is not defined");
    }


    cout << PlasmaNumberDensityOffset << "  " << PlasmaTemperatureOffset << "  " << PlasmaPressureOffset <<  "   "  <<  MagneticFieldOffset <<  "   " << PlasmaBulkVelocityOffset << endl;


  } //end of the initialization


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
    int IsFound;

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i))/UnitLength;
      x[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j))/UnitLength;
      x[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k))/UnitLength;

      //recover the data from the BATSRUS data file
      batsrus2amps_get_data_point_(x,StateLocal,&IsFound);

      //the found state vector could be located on onother processor
      if (IsFound==true) {
        MPI_Allreduce(StateLocal,State,nVar+1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
        //devide the interpolated values by the "total weight"
        for (int n=1;n<nVar+1;n++) State[n]/=State[0];
      }

      //locate the cell
      if (startNode->block==NULL) continue;

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
      offset=CenterNode->GetAssociatedDataBufferPointer();

      //save the interpolated values
      if (IsFound==true) {
        int idim;

        //if momentum is read that convert it to velocity
        if (mxBATSRUS2AMPS>=0) {
          exit(__LINE__,__FILE__,"Error: conversion of the momentun into the velocity is not implemented. The implementation should be simular to that in PIC::CPLR::SWMF::RecieveCenterPointData");
        }

        //the order of the state vector: number density, temperature
        *((double*)(offset+PlasmaNumberDensityOffset))=State[rhoBATSRUS2AMPS]*PhysicalVariableUnitConversionTable[rhoBATSRUS2AMPS];
        *((double*)(offset+PlasmaTemperatureOffset))=(State[rhoBATSRUS2AMPS]>0.0) ? PhysicalVariableUnitConversionTable[pBATSRUS2AMPS]*State[pBATSRUS2AMPS]/(Kbol*State[rhoBATSRUS2AMPS]*PhysicalVariableUnitConversionTable[rhoBATSRUS2AMPS]) : 0.0;

        //get pressure
        *((double*)(offset+PlasmaPressureOffset))=State[pBATSRUS2AMPS]*PhysicalVariableUnitConversionTable[pBATSRUS2AMPS];

        //bulk velocity and magnetic field
        for (idim=0;idim<3;idim++) {
          *((double*)(offset+MagneticFieldOffset+idim*sizeof(double)))=State[bxBATSRUS2AMPS+idim]*PhysicalVariableUnitConversionTable[bxBATSRUS2AMPS];

          if (uxBATSRUS2AMPS>=0) {
            *((double*)(offset+PlasmaBulkVelocityOffset+idim*sizeof(double)))=State[uxBATSRUS2AMPS+idim]*PhysicalVariableUnitConversionTable[uxBATSRUS2AMPS];
          }
          else {
            exit(__LINE__,__FILE__,"Error: saving of the velocity derived from the momentum is not implemented");
          }
        }

        //calculate the electric field
        double *E,*B,*v;

        v=(double*)(offset+PlasmaBulkVelocityOffset);
        B=(double*)(offset+MagneticFieldOffset);
        E=(double*)(offset+ElectricFieldOffset);

        E[0]=-(v[1]*B[2]-B[1]*v[2]);
        E[1]=+(v[0]*B[2]-B[0]*v[2]);
        E[2]=-(v[0]*B[1]-B[0]*v[1]);
      }
      else {
        *((double*)(offset+PlasmaNumberDensityOffset))=0.0;
        *((double*)(offset+PlasmaTemperatureOffset))=0.0;
        *((double*)(offset+PlasmaPressureOffset))=0.0;

        for (idim=0;idim<3;idim++) {
          *((double*)(offset+PlasmaBulkVelocityOffset+idim*sizeof(double)))=0.0;
          *((double*)(offset+MagneticFieldOffset+idim*sizeof(double)))=0.0;
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(startNode->downNode[nDownNode]);
  }


  if (startNode==PIC::Mesh::mesh.rootTree) {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    delete [] State;
    delete [] StateLocal;
    delete [] PhysicalVariableUnitConversionTable;

    batsrus2amps_closefile_();
  }
}



