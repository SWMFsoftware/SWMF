//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//The file containes function for using ICES


#include "pic.h"


char PIC::CPLR::DATAFILE::ICES::locationICES[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/ices/ICES/Models";
char PIC::CPLR::DATAFILE::ICES::ModeCaseSWMF[_MAX_STRING_LENGTH_PIC_]="";



int PIC::CPLR::DATAFILE::ICES::TotalAssociatedDataLength=0,PIC::CPLR::DATAFILE::ICES::AssociatedDataOffset=-1;

//offsets of the data loaded from the DSMC results
int PIC::CPLR::DATAFILE::ICES::NeutralBullVelocityOffset=-1,PIC::CPLR::DATAFILE::ICES::NeutralNumberDensityOffset=-1,PIC::CPLR::DATAFILE::ICES::NeutralTemperatureOffset=-1,PIC::CPLR::DATAFILE::ICES::DataStatusOffsetDSMC=-1;

//pre-processor of the data
PIC::CPLR::DATAFILE::ICES::fDSMCdataPreProcessor PIC::CPLR::DATAFILE::ICES::DSMCdataPreProcessor=NULL;
PIC::CPLR::DATAFILE::ICES::fSWMFdataPreProcessor PIC::CPLR::DATAFILE::ICES::SWMFdataPreProcessor=NULL;

//====================================================
void PIC::CPLR::DATAFILE::ICES::SetLocationICES(const char* loc) {
  sprintf(locationICES,"%s",loc);
}

//====================================================
void PIC::CPLR::DATAFILE::ICES::Init() {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  //reserve the data fields
  PIC::CPLR::DATAFILE::Offset::ElectricField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.allocate=true;
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  exit(__LINE__,__FILE__,"Error: the DSMC part is not updated for using ICES as PIC::CPLR::DATAFILE::ICES ");


  if (NeutralBullVelocityOffset!=-1) exit(__LINE__,__FILE__,"Error: the model is already initialied");

  if (NeutralBullVelocityOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    NeutralBullVelocityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
    TotalAssociatedDataLength+=3*sizeof(double);
  }

  if (NeutralNumberDensityOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    NeutralNumberDensityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  if (NeutralTemperatureOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    NeutralTemperatureOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  //!!The data center associated data can contain only 'double' type data -> memory reserved for DataStatusOffsetDSMC as for double
  DataStatusOffsetDSMC=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
  TotalAssociatedDataLength+=sizeof(double);
#endif

}

//====================================================
//retrive the data file from SWMF
void PIC::CPLR::DATAFILE::ICES::retriveSWMFdata(const char *DataFile) {
  char cCurrentPath[_MAX_STRING_LENGTH_PIC_],command[_MAX_STRING_LENGTH_PIC_],initDirectory[_MAX_STRING_LENGTH_PIC_];

  //check if the model is initialied
  if (PIC::CPLR::DATAFILE::Offset::ElectricField.offset==-1) exit(__LINE__,__FILE__,"Error: the model is not initialied");

  //create the directory for the trajectory file
  sprintf(command,"rm -f -r temp.ICES.thread=%i",PIC::Mesh::mesh.ThisThread);
  system(command);

  sprintf(command,"mkdir temp.ICES.thread=%i",PIC::Mesh::mesh.ThisThread);
  system(command);

  getcwd(initDirectory,_MAX_STRING_LENGTH_PIC_);
  sprintf(command,"%s/temp.ICES.thread=%i",initDirectory,PIC::Mesh::mesh.ThisThread);
  chdir(command);

  //copy the trajectory file in the currect directory
  sprintf(command,"cp ../icesCellCenterCoordinates.thread=%i .",PIC::Mesh::mesh.ThisThread);
  system(command);

  sprintf(command,"mv icesCellCenterCoordinates.thread=%i mhd_traj.dat",PIC::Mesh::mesh.ThisThread);
  system(command);

  //start the ices
  getcwd(cCurrentPath, sizeof(cCurrentPath));
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"Getting the MHD data from ICES.... \nThe current working directory is: %s\n", cCurrentPath);

  system("rm -f MHDRestart");

  sprintf(command,"ln -s %s/Data/%s/MHD/ MHDRestart",locationICES,DataFile);
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"%s\n",command);
  system(command);

  if (PIC::Mesh::mesh.ThisThread==0) {
    sprintf(command,"%s/MHD/comet_mhd.exe",locationICES);
    fprintf(PIC::DiagnospticMessageStream,"%s\n",command);
  }
  else sprintf(command,"%s/MHD/comet_mhd.exe > /dev/null",locationICES);

  system(command);

  system("rm MHDRestart");

  //copy the output file into the working directory
  sprintf(command,"mv mhd_values.dat ../icesCellCenterCoordinates.thread=%i.MHD.dat",PIC::Mesh::mesh.ThisThread);
  system(command);

  chdir(initDirectory);

  sprintf(command,"rm -f -r temp.ICES.thread=%i",PIC::Mesh::mesh.ThisThread);
  system(command);

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if (PIC::Mesh::mesh.ThisThread==0) {
    fprintf(PIC::DiagnospticMessageStream,"Combine individual trajectory files into a single file\n");

    FILE *fout;
    int thread;
    char fname[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_];
    CiFileOperations fin;

    fout=fopen("icesCellCenterCoordinates.MHD.dat","w");

    for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
      sprintf(fname,"icesCellCenterCoordinates.thread=%i.MHD.dat",thread);
      fin.openfile(fname);

      while (fin.eof()==false) {
        if (fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false)  {
          break;
        }

        fin.CutInputStr(str1,str);

        if ((strcmp("#VARIABLES",str1)==0)||(strcmp("#START",str1)==0)) {
          if (thread==0) fprintf(fout,"%s %s\n",str1,str);
        }
        else fprintf(fout,"%s %s\n",str1,str);
      }

      fin.closefile();


    }

    fclose(fout);
  }


  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"ICES (SWMF) done\n");
}


//====================================================
//retrive the data file from DSMC solver
void PIC::CPLR::DATAFILE::ICES::retriveDSMCdata(const char *Case,const char *DataFile,const char *MeshFile) {
  char command[_MAX_STRING_LENGTH_PIC_];

  exit(__LINE__,__FILE__,"Error: the DSMC part is not updated for using ICES as PIC::CPLR::DATAFILE::ICES ");

  //check if the model if initialied
  if (NeutralBullVelocityOffset==-1) exit(__LINE__,__FILE__,"Error: the model is already initialied");

  //start ICES
  sprintf(command,"%s/DSMC/ices-dsmc -extractdatapoints -grid %s/Data/%s/DSMC/%s -testpointslist icesCellCenterCoordinates.thread=%i -datafile %s/Data/%s/DSMC/%s -dim 2 -symmetry cylindrical",locationICES,   locationICES,Case,MeshFile, PIC::Mesh::mesh.ThisThread,   locationICES,Case,DataFile);

  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"ICES: retrive DSMC data \n Execute command: %s\n",command);
  system(command);

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if (PIC::Mesh::mesh.ThisThread==0) {
    fprintf(PIC::DiagnospticMessageStream,"Combine individual trajectory files into a single file\n");

    FILE *fout;
    int thread;
    char fname[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_];
    CiFileOperations fin;

    fout=fopen("icesCellCenterCoordinates.DSMC.dat","w");

    for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
      sprintf(fname,"icesCellCenterCoordinates.thread=%i.out",thread);
      fin.openfile(fname);

      while (fin.eof()==false) {
        if (fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false)  {
          break;
        }

        fin.CutInputStr(str1,str);

        if ((strcmp("#VARIABLES",str1)==0)||(strcmp("#START",str1)==0)) {
          if (thread==0) fprintf(fout,"%s %s\n",str1,str);
        }
        else fprintf(fout,"%s %s\n",str1,str);
      }

      fin.closefile();


    }

    fclose(fout);
  }


  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"ICES (DSMC) done\n");

}



//====================================================
//read and parse the data file from SWMF
void PIC::CPLR::DATAFILE::ICES::readSWMFdata(const double MeanIonMass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  cDataNodeSWMF dataSWMF;
  static CiFileOperations ices;
  long int nd;
  int status,idim,i,j,k;
  char *offset;

  PIC::Mesh::cDataCenterNode *CenterNode;
  char str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_],*endptr;

  //check if the model is initialied
  if (PIC::CPLR::DATAFILE::Offset::ElectricField.offset==-1) exit(__LINE__,__FILE__,"Error: the model is not initialied");


#if DIM == 3
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
#elif DIM == 2
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=0;
#elif DIM == 1
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=0;
  const int kMin=0,kMax=0;
#endif

  if (startNode==PIC::Mesh::mesh.rootTree) {
    ices.openfile("icesCellCenterCoordinates.MHD.dat");

    ices.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
    ices.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
  }

  //read the data file
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the read section
      ices.GetInputStr(str,sizeof(str));
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);

      status=strtol(str1,&endptr,10);
      dataSWMF.status=status;

      if (status==0) {
        //read the row data
        ices.CutInputStr(str1,str);
        dataSWMF.swNumberDensity=strtod(str1,&endptr);

        for (idim=0;idim<3;idim++) {
          ices.CutInputStr(str1,str);
          dataSWMF.swVel[idim]=strtod(str1,&endptr);
        }

        for (idim=0;idim<3;idim++) {
          ices.CutInputStr(str1,str);
          dataSWMF.B[idim]=strtod(str1,&endptr);
        }

        ices.CutInputStr(str1,str);
        dataSWMF.swPressure=strtod(str1,&endptr);
      }
      else dataSWMF.flush();


      if (startNode->block!=NULL) {
        //calculation of other plasma parameters
        /* ---------- calculation equations
            pV=NkT
            p=ne*k*Te+ni+k*Ti for single ion species + electrons
            p=2*n*k*T assuming quasi neutrality ni=ne=n and Ti=Te=T
            T=p/(2*n*T)
        ---------- end of calculation equations */



        dataSWMF.swNumberDensity/=MeanIonMass; // MHD output in [amu/m^3], mi = mean ion mass
        dataSWMF.swTemperature=dataSWMF.swPressure/(2.0*dataSWMF.swNumberDensity*Kbol);

        dataSWMF.E[0]=-(dataSWMF.swVel[1]*dataSWMF.B[2]-dataSWMF.swVel[2]*dataSWMF.B[1]);
        dataSWMF.E[1]=(dataSWMF.swVel[0]*dataSWMF.B[2]-dataSWMF.swVel[2]*dataSWMF.B[0]);
        dataSWMF.E[2]=-(dataSWMF.swVel[0]*dataSWMF.B[1]-dataSWMF.swVel[1]*dataSWMF.B[0]);

        //save the data on the mesh node
        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        CenterNode=startNode->block->GetCenterNode(nd);

        if (CenterNode==NULL) continue;

        //preprocess the data if needed
        if (SWMFdataPreProcessor!=NULL) SWMFdataPreProcessor(CenterNode->GetX(),dataSWMF);

        //save the data on the mesh
        offset=CenterNode->GetAssociatedDataBufferPointer();

/*
        if (status!=0) { //check the status of the reading
          ices.error("the point is not found");
          exit(__LINE__,__FILE__,"Error: the extracted point is not found");
        }
*/

        for (idim=0;idim<3;idim++) {
          *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.offset))=dataSWMF.E[idim];
          *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.offset))=dataSWMF.B[idim];
          *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.offset))=dataSWMF.swVel[idim];
        }

        *(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.offset)=dataSWMF.swPressure;
        *(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.offset)=dataSWMF.swNumberDensity;
        *(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.offset)=dataSWMF.swTemperature;

//        *(int*)(offset+DataStatusOffsetSWMF)=dataSWMF.status;
      }
    }
  }
  else {
    int nDownNode;

    for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) readSWMFdata(MeanIonMass,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) ices.closefile();
}

//====================================================
//read and parse the data file from SWMF
void PIC::CPLR::DATAFILE::ICES::readDSMCdata(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  cDataNodeDSMC dataDSMC;
  static CiFileOperations ices;
  long int nd;
  int status,idim,i,j,k;
  char *offset;

  exit(__LINE__,__FILE__,"Error: the DSMC part is not updated for using ICES as PIC::CPLR::DATAFILE::ICES ");

  PIC::Mesh::cDataCenterNode *CenterNode;
  char str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_],*endptr;

  //check if the model if initialied
  if (NeutralBullVelocityOffset==-1) exit(__LINE__,__FILE__,"Error: the model is already initialied");


#if DIM == 3
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
#elif DIM == 2
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=0;
#elif DIM == 1
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=0;
  const int kMin=0,kMax=0;
#endif

  if (startNode==PIC::Mesh::mesh.rootTree) {
    ices.openfile("icesCellCenterCoordinates.DSMC.dat");

    ices.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
    ices.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
  }

  //read the data file
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the read section
      ices.GetInputStr(str,sizeof(str));
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);
      ices.CutInputStr(str1,str);

      status=strtol(str1,&endptr,10);
      dataDSMC.status=status;

      if (status==0) {
        //read the row data
        for (idim=0;idim<3;idim++) {
          ices.CutInputStr(str1,str);
          dataDSMC.neutralVel[idim]=strtod(str1,&endptr);
        }

        ices.CutInputStr(str1,str);
        dataDSMC.neutralTemperature=strtod(str1,&endptr);

        ices.CutInputStr(str1,str);
        dataDSMC.neutralNumberDensity=strtod(str1,&endptr);
      }
      else dataDSMC.flush();


      if (startNode->block!=NULL) {
        //save the data on the mesh node
        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        CenterNode=startNode->block->GetCenterNode(nd);

        if (CenterNode==NULL) continue;

        //preprocess the data if needed
        if (DSMCdataPreProcessor!=NULL) DSMCdataPreProcessor(CenterNode->GetX(),dataDSMC);

        //save the data on the mesh
        offset=CenterNode->GetAssociatedDataBufferPointer();

/*
        if (status!=0) { //check the status of the reading
          ices.error("the point is not found");
          exit(__LINE__,__FILE__,"Error: the extracted point is not found");
        }
*/

        for (idim=0;idim<3;idim++) {
          *(idim+(double*)(offset+NeutralBullVelocityOffset))=dataDSMC.neutralVel[idim];
        }

        *(double*)(offset+NeutralNumberDensityOffset)=dataDSMC.neutralNumberDensity;
        *(double*)(offset+NeutralTemperatureOffset)=dataDSMC.neutralTemperature;

        *(int*)(offset+DataStatusOffsetDSMC)=dataDSMC.status;
      }
    }
  }
  else {
    int nDownNode;

    for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) readDSMCdata(startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) ices.closefile();
}

//====================================================
//the total number of cells on the mesh
long int PIC::CPLR::DATAFILE::ICES::getTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int res=0;
  int nDownNode;
  bool flag=false;

  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) {
    res+=getTotalCellNumber(startNode->downNode[nDownNode]);
    flag=true;
  }

  if (flag==false) res=1;

  return res;
}

//====================================================
void PIC::CPLR::DATAFILE::ICES::createCellCenterCoordinateList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];
  double x[3]={0.0,0.0,0.0},*xNodeMin,*xNodeMax;
  int idim,i,j,k;

  static long int nTotalCellNumber,startCellNumber,stopCellNumber,CellCounter;

#if DIM == 3
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
#elif DIM == 2
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=0;
#elif DIM == 1
  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=0;
  const int kMin=0,kMax=0;
#endif

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //generate the total file name and open the file
    sprintf(fname,"icesCellCenterCoordinates.thread=%i",PIC::Mesh::mesh.ThisThread);
    fout=fopen(fname,"w");
    PIC::Mesh::mesh.resetNodeProcessedFlag();

    fprintf(fout,"#START\n");

    //determine the limits of the cells that are printed
    long int nCellPerProcessor;

    nTotalCellNumber=getTotalCellNumber(PIC::Mesh::mesh.rootTree);
    nCellPerProcessor=nTotalCellNumber/PIC::Mesh::mesh.nTotalThreads;

    startCellNumber=PIC::Mesh::mesh.ThisThread*nCellPerProcessor;
    stopCellNumber=startCellNumber+nCellPerProcessor-1;
    CellCounter=0;

    if (PIC::Mesh::mesh.ThisThread==PIC::Mesh::mesh.nTotalThreads-1) stopCellNumber=nTotalCellNumber-1;
  }

  //output the cell's coordinates
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((startCellNumber<=CellCounter)&&(CellCounter<=stopCellNumber)) {
      xNodeMin=startNode->xmin;
      xNodeMax=startNode->xmax;

      for (k=kMin;k<=kMax;k++) {
        if (DIM==3) x[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);

        for (j=jMin;j<=jMax;j++) {
          if (DIM>=2) x[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);

          for (i=iMin;i<=iMax;i++) {
            x[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);

            for (idim=0;idim<DIM;idim++) fprintf(fout,"%15.10e ",x[idim]);
            fprintf(fout,"\n");
          }
        }
      }
    }

    CellCounter++;
  }
  else {
    int nDownNode;

    for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) createCellCenterCoordinateList(startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) fclose(fout);
}





























