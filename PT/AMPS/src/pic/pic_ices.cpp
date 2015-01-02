//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//The file containes function for using ICES


#include "pic.h"


char PIC::CPLR::ICES::locationICES[_MAX_STRING_LENGTH_PIC_]="/Users/vtenishe/ices/ICES/Models";
char PIC::CPLR::ICES::ModeCaseSWMF[_MAX_STRING_LENGTH_PIC_]="";

//background plasma parameter's offsets
int PIC::CPLR::ICES::ElectricFieldOffset=-1;
int PIC::CPLR::ICES::MagneticFieldOffset=-1;
int PIC::CPLR::ICES::PlasmaPressureOffset=-1;
int PIC::CPLR::ICES::PlasmaNumberDensityOffset=-1;
int PIC::CPLR::ICES::PlasmaTemperatureOffset=-1;
int PIC::CPLR::ICES::PlasmaBulkVelocityOffset=-1;
int PIC::CPLR::ICES::DataStatusOffsetSWMF=-1;

int PIC::CPLR::ICES::TotalAssociatedDataLength=0,PIC::CPLR::ICES::AssociatedDataOffset=-1;

//offsets of the data loaded from the DSMC results
int PIC::CPLR::ICES::NeutralBullVelocityOffset=-1,PIC::CPLR::ICES::NeutralNumberDensityOffset=-1,PIC::CPLR::ICES::NeutralTemperatureOffset=-1,PIC::CPLR::ICES::DataStatusOffsetDSMC=-1;

//pre-processor of the data
PIC::CPLR::ICES::fDSMCdataPreProcessor PIC::CPLR::ICES::DSMCdataPreProcessor=NULL;
PIC::CPLR::ICES::fSWMFdataPreProcessor PIC::CPLR::ICES::SWMFdataPreProcessor=NULL;

//====================================================
void PIC::CPLR::ICES::SetLocationICES(const char* loc) {
  sprintf(locationICES,"%s",loc);
}

//====================================================
void PIC::CPLR::ICES::Init() {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  if (ElectricFieldOffset!=-1) exit(__LINE__,__FILE__,"Error: the model is already initialied");

  //init the plasma parameters
  if (ElectricFieldOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    ElectricFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
    TotalAssociatedDataLength+=3*sizeof(double);
  }

  if (MagneticFieldOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    MagneticFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
    TotalAssociatedDataLength+=3*sizeof(double);
  }

  if (PlasmaPressureOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PlasmaPressureOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaNumberDensityOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PlasmaNumberDensityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaTemperatureOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PlasmaTemperatureOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaBulkVelocityOffset==-1) {
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PlasmaBulkVelocityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
    TotalAssociatedDataLength+=3*sizeof(double);
  }

  //!!The data center associated data can contain only 'double' type data -> memory reserved for DataStatusOffsetSWMF as for double
  DataStatusOffsetSWMF=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
  PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
  TotalAssociatedDataLength+=sizeof(double);
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
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


  //register the file output functions
  PIC::Mesh::fPrintVariableListCenterNode t0; //define the temporaty variables to satisfy the intell c++ compiler
  PIC::Mesh::fPrintDataCenterNode t1;
  PIC::Mesh::fInterpolateCenterNode t2;


  PIC::Mesh::PrintVariableListCenterNode.push_back(t0=PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(t1=PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(t2=Interpolate);


}

//====================================================
//retrive the data file from SWMF
void PIC::CPLR::ICES::retriveSWMFdata(const char *DataFile) {
  char cCurrentPath[_MAX_STRING_LENGTH_PIC_],command[_MAX_STRING_LENGTH_PIC_],initDirectory[_MAX_STRING_LENGTH_PIC_];

  //check if the model is initialied
  if (ElectricFieldOffset==-1) exit(__LINE__,__FILE__,"Error: the model is not initialied");

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
void PIC::CPLR::ICES::retriveDSMCdata(const char *Case,const char *DataFile,const char *MeshFile) {
  char command[_MAX_STRING_LENGTH_PIC_];

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
void PIC::CPLR::ICES::readSWMFdata(const double MeanIonMass,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  cDataNodeSWMF dataSWMF;
  static CiFileOperations ices;
  long int nd;
  int status,idim,i,j,k;
  char *offset;

  PIC::Mesh::cDataCenterNode *CenterNode;
  char str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_],*endptr;

  //check if the model is initialied
  if (ElectricFieldOffset==-1) exit(__LINE__,__FILE__,"Error: the model is not initialied");


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
          *(idim+(double*)(offset+ElectricFieldOffset))=dataSWMF.E[idim];
          *(idim+(double*)(offset+MagneticFieldOffset))=dataSWMF.B[idim];
          *(idim+(double*)(offset+PlasmaBulkVelocityOffset))=dataSWMF.swVel[idim];
        }

        *(double*)(offset+PlasmaPressureOffset)=dataSWMF.swPressure;
        *(double*)(offset+PlasmaNumberDensityOffset)=dataSWMF.swNumberDensity;
        *(double*)(offset+PlasmaTemperatureOffset)=dataSWMF.swTemperature;

        *(int*)(offset+DataStatusOffsetSWMF)=dataSWMF.status;
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
void PIC::CPLR::ICES::readDSMCdata(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  cDataNodeDSMC dataDSMC;
  static CiFileOperations ices;
  long int nd;
  int status,idim,i,j,k;
  char *offset;

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
//output background parameters loaded with ICES from SWMF output
void PIC::CPLR::ICES::PrintVariableList(FILE* fout,int DataSetNumber) {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  for (int idim=0;idim<3;idim++) fprintf(fout,", \"E%i\"",idim);
  for (int idim=0;idim<3;idim++) fprintf(fout,", \"B%i\"",idim);
  for (int idim=0;idim<3;idim++) fprintf(fout,", \"PlasmaVel%i\"",idim);

  fprintf(fout,", \"Plasma Pressure\", \"Plasma Temperature\", \"Plasma number Desnity\"");
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  for (int idim=0;idim<3;idim++) fprintf(fout,", \"NeutralVel%i\"",idim);
  fprintf(fout,", \"Neutral number Desnity\", \"Neutral Temperature\"");
#endif

}

void PIC::CPLR::ICES::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  cDataNodeSWMF dataSWMF;
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  cDataNodeDSMC dataDSMC;
#endif

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  char *offset=CenterNode->GetAssociatedDataBufferPointer();
#elif _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  char *offset=CenterNode->GetAssociatedDataBufferPointer();
#endif


  if (pipe->ThisThread==CenterNodeThread) {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
    dataSWMF.swNumberDensity=*((double*)(PlasmaNumberDensityOffset+offset));
    dataSWMF.swPressure=*((double*)(PlasmaPressureOffset+offset));
    dataSWMF.swTemperature=*((double*)(PlasmaTemperatureOffset+offset));

    for (int idim=0;idim<3;idim++) {
      dataSWMF.B[idim]=*(idim+(double*)(MagneticFieldOffset+offset));
      dataSWMF.E[idim]=*(idim+(double*)(ElectricFieldOffset+offset));
      dataSWMF.swVel[idim]=*(idim+(double*)(PlasmaBulkVelocityOffset+offset));
    }
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
    for (int idim=0;idim<3;idim++) dataDSMC.neutralVel[idim]=*(idim+(double*)(offset+NeutralBullVelocityOffset));

    char *offset=CenterNode->GetAssociatedDataBufferPointer();
    dataDSMC.neutralNumberDensity=*(double*)(offset+NeutralNumberDensityOffset);
    dataDSMC.neutralTemperature=*(double*)(offset+NeutralTemperatureOffset);
#endif

  }


  if (pipe->ThisThread==0) {
#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
    if (CenterNodeThread!=0) pipe->recv((char*)&dataSWMF,sizeof(dataSWMF),CenterNodeThread);

    fprintf(fout,"%e  %e  %e  ",dataSWMF.E[0],dataSWMF.E[1],dataSWMF.E[2]);
    fprintf(fout,"%e  %e  %e  ",dataSWMF.B[0],dataSWMF.B[1],dataSWMF.B[2]);
    fprintf(fout,"%e  %e  %e  ",dataSWMF.swVel[0],dataSWMF.swVel[1],dataSWMF.swVel[2]);

    fprintf(fout,"%e  %e  %e  ",dataSWMF.swPressure,dataSWMF.swTemperature,dataSWMF.swNumberDensity);
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
    if (CenterNodeThread!=0) pipe->recv((char*)&dataDSMC,sizeof(dataDSMC),CenterNodeThread);

    for (int idim=0;idim<3;idim++) fprintf(fout,"%e  ",dataDSMC.neutralVel[idim]);
    fprintf(fout,"%e  %e  ",dataDSMC.neutralNumberDensity,dataDSMC.neutralTemperature);
#endif

  }
  else {
#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
    pipe->send((char*)&dataSWMF,sizeof(dataSWMF));
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
    pipe->send((char*)&dataDSMC,sizeof(dataDSMC));
#endif

  }
}

void PIC::CPLR::ICES::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_ || _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  int i;
  double c;
  char *offset,*offsetCenterNode;

  #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if (CenterNode->GetAssociatedDataBufferPointer()==NULL) exit(__LINE__,__FILE__,"Error: The associated data buffer is not initialized");
  #endif

  //init the buffer of 'CenterNode'
  offsetCenterNode=CenterNode->GetAssociatedDataBufferPointer();

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  *((double*)(PlasmaNumberDensityOffset+offsetCenterNode))=0.0;
  *((double*)(PlasmaPressureOffset+offsetCenterNode))=0.0;
  *((double*)(PlasmaTemperatureOffset+offsetCenterNode))=0.0;

  for (int idim=0;idim<3;idim++) {
    *(idim+(double*)(MagneticFieldOffset+offsetCenterNode))=0.0;
    *(idim+(double*)(ElectricFieldOffset+offsetCenterNode))=0.0;
    *(idim+(double*)(PlasmaBulkVelocityOffset+offsetCenterNode))=0.0;
  }
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
  for (int idim=0;idim<3;idim++) {
    *(idim+(double*)(NeutralBullVelocityOffset+offsetCenterNode))=0.0;
  }

  *(double*)(NeutralNumberDensityOffset+offsetCenterNode)=0.0;
  *(double*)(NeutralTemperatureOffset+offsetCenterNode)=0.0;
#endif

  //interpolate the sampled data
  for (i=0;i<nInterpolationCoeficients;i++) {
    c=InterpolationCoeficients[i];
    offset=InterpolationList[i]->GetAssociatedDataBufferPointer();

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
    *((double*)(PlasmaNumberDensityOffset+offsetCenterNode))+=c*(*((double*)(PlasmaNumberDensityOffset+offset)));
    *((double*)(PlasmaPressureOffset+offsetCenterNode))+=c*(*((double*)(PlasmaPressureOffset+offset)));
    *((double*)(PlasmaTemperatureOffset+offsetCenterNode))+=c*(*((double*)(PlasmaTemperatureOffset+offset)));

    for (int idim=0;idim<3;idim++) {
      *(idim+(double*)(MagneticFieldOffset+offsetCenterNode))+=c*(*(idim+(double*)(MagneticFieldOffset+offset)));
      *(idim+(double*)(ElectricFieldOffset+offsetCenterNode))+=c*(*(idim+(double*)(ElectricFieldOffset+offset)));
      *(idim+(double*)(PlasmaBulkVelocityOffset+offsetCenterNode))+=c*(*(idim+(double*)(PlasmaBulkVelocityOffset+offset)));
    }
#endif

#if _PIC_ICES_DSMC_MODE_ == _PIC_ICES_MODE_ON_
    for (int idim=0;idim<3;idim++) {
      *(idim+(double*)(NeutralBullVelocityOffset+offsetCenterNode))+=c*(*(idim+(double*)(NeutralBullVelocityOffset+offset)));
    }

    *(double*)(NeutralNumberDensityOffset+offsetCenterNode)+=c*(*(double*)(NeutralNumberDensityOffset+offset));
    *(double*)(NeutralTemperatureOffset+offsetCenterNode)+=c*(*(double*)(NeutralTemperatureOffset+offset));
#endif


  }
#endif
}

//====================================================
//the total number of cells on the mesh
long int PIC::CPLR::ICES::getTotalCellNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
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
void PIC::CPLR::ICES::createCellCenterCoordinateList(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
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


//====================================================
//print the ion flux at a sphere
void PIC::CPLR::ICES::PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius) {
  FILE *fout=NULL;
  int nZenithAngle,nPolarAngle,LocalCellNumber,i,j,k;
  double ZenithAngle,PolarAngle,x[3],n,v[3],flux;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.rootTree;
  CMPI_channel pipe(1000000);

  const int nTotalZenithPoints=100;
  const double dZenithAngle=Pi/(nTotalZenithPoints-1);
  const int nTotalPolarPoints=2*nTotalZenithPoints;
  const double dPolarAngle=2.0*Pi/(nTotalPolarPoints-1);

  if (PIC::ThisThread==0) {
    pipe.openRecvAll();

    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Lon\", \"Lat\", \"Ion Flux\"\n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nTotalPolarPoints,nTotalZenithPoints);
  }
  else pipe.openSend(0);

  //calculate the flux distribution
  for (nZenithAngle=0,ZenithAngle=-Pi/2.0;nZenithAngle<nTotalZenithPoints;nZenithAngle++,ZenithAngle+=dZenithAngle) {
    for (nPolarAngle=0,PolarAngle=-Pi;nPolarAngle<nTotalPolarPoints;nPolarAngle++,PolarAngle+=dPolarAngle) {
      x[0]=SphereRadius*cos(PolarAngle)*cos(ZenithAngle);
      x[1]=SphereRadius*sin(PolarAngle)*cos(ZenithAngle);
      x[2]=SphereRadius*sin(ZenithAngle);

      node=PIC::Mesh::mesh.findTreeNode(x,node);
      if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find cell");

      if (node->Thread==PIC::ThisThread) {
        n=GetBackgroundPlasmaNumberDensity(x,LocalCellNumber,node);
        GetBackgroundPlasmaVelocity(v,x,LocalCellNumber,node);

        flux=-n*(v[0]*x[0]+v[1]*x[1]+v[2]*x[2]);

        if (PIC::ThisThread!=0) pipe.send(flux);
      }

      if (PIC::ThisThread==0) {
        if (node->Thread!=0) flux=pipe.recv<double>(node->Thread);

        fprintf(fout,"%e  %e  %e\n",PolarAngle/Pi*180.0,ZenithAngle/Pi*180.0,flux);
      }

    }
  }


  if (PIC::ThisThread==0) {
    pipe.closeRecvAll();
    fclose(fout);
  }
  else pipe.closeSend();

}


//====================================================
//evaluate the ion flux at the surface
void PIC::CPLR::ICES::EvaluateSurfaceIonFlux(double ShiftFactor) {
  int nTotalElements,el;
  list<cInternalBoundaryConditionsDescriptor>::iterator d;
  cInternalSphericalData* Sphere;
  double x[3],norm[3],l,n,v[3],flux;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.rootTree;
  int LocalCellNumber,i,j,k;


  for (d=PIC::Mesh::mesh.InternalBoundaryList.begin();d!=PIC::Mesh::mesh.InternalBoundaryList.end();d++) {
    if (d->BondaryType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
      Sphere=(cInternalSphericalData*)d->BoundaryElement;
      nTotalElements=Sphere->GetTotalSurfaceElementsNumber();

      //assemble the table of ion fluxes
      double FluxTable[nTotalElements],TotalSurfaceFlux;
      CMPI_channel pipe(1000000);

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else pipe.openSend(0);

      for (el=0;el<nTotalElements;el++) {
        Sphere->GetSurfaceElementNormal(norm,el);
        Sphere->GetSurfaceElementMiddlePoint(x,el);

        l=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]+=max(ShiftFactor-1.0,0.0)*l*norm[0],x[1]+=max(ShiftFactor-1.0,0.0)*l*norm[1],x[2]+=max(ShiftFactor-1.0,0.0)*l*norm[2];

        node=PIC::Mesh::mesh.findTreeNode(x,node);
        if ((LocalCellNumber=PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find cell");

        if (node->Thread==PIC::ThisThread) {
          n=GetBackgroundPlasmaNumberDensity(x,LocalCellNumber,node);
          GetBackgroundPlasmaVelocity(v,x,LocalCellNumber,node);

          flux=-n*(v[0]*x[0]+v[1]*x[1]+v[2]*x[2]);

          if (PIC::ThisThread!=0) pipe.send(flux);
        }

        if (PIC::ThisThread==0) {
          if (node->Thread!=0) flux=pipe.recv<double>(node->Thread);

          FluxTable[el]=flux;
        }
      }

      if (PIC::ThisThread==0) {
        pipe.closeRecvAll();
      }
      else pipe.closeSend();

      //broundcast the table across all processors
      MPI_Bcast(FluxTable,nTotalElements,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

      //save the table on the sphere
      for (el=0,TotalSurfaceFlux=0.0;el<nTotalElements;el++) {
        Sphere->SolarWindSurfaceFlux[el]=FluxTable[el];
        if (FluxTable[el]>0.0) TotalSurfaceFlux+=FluxTable[el]*Sphere->SurfaceElementArea[el];
      }

      Sphere->TotalSolarWindSurfaceFlux=TotalSurfaceFlux;

    }
    else {
      exit(__LINE__,__FILE__,"Error: unknown internal boundary type");
    }
  }
}




























