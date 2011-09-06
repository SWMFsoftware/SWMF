//====================================================
//$Id$
//====================================================
//The file containes function for using ICES


#include "pic.h"


char PIC::ICES::locationICES[_MAX_STRING_LENGTH_PIC_]; //location of the data and the dace cases

//background plasma parameter's offsets
int PIC::ICES::ElectricFieldOffset=-1;
int PIC::ICES::MagneticFieldOffset=-1;
int PIC::ICES::PlasmaPressureOffset=-1;
int PIC::ICES::PlasmaNumberDensityOffset=-1;
int PIC::ICES::PlasmaTemperatureOffset=-1;
int PIC::ICES::PlasmaBulkVelocityOffset=-1;

//====================================================
void PIC::ICES::SetLocationICES(const char* loc) {
  sprintf(locationICES,"%s",loc);
}

//====================================================
void PIC::ICES::Init() {

#if _PIC_ICES_SWMF_MODE_ == _PIC_ICES_MODE_ON_
  //init the plasma parameters
  if (ElectricFieldOffset==-1) {
    ElectricFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
  }

  if (MagneticFieldOffset==-1) {
    MagneticFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
  }

  if (PlasmaPressureOffset==-1) {
    PlasmaPressureOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaNumberDensityOffset==-1) {
    PlasmaNumberDensityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaTemperatureOffset==-1) {
    PlasmaTemperatureOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
  }

  if (PlasmaBulkVelocityOffset==-1) {
    PlasmaBulkVelocityOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
  }


  //register the file output functions
  PIC::Mesh::fPrintVariableListCenterNode t0; //define the temporaty variables to satisfy the intell c++ compiler
  PIC::Mesh::fPrintDataCenterNode t1;
  PIC::Mesh::fInterpolateCenterNode t2;


  PIC::Mesh::PrintVariableListCenterNode.push_back(t0=PrintVariableListSWMF);
  PIC::Mesh::PrintDataCenterNode.push_back(t1=PrintDataSWMF);
  PIC::Mesh::InterpolateCenterNode.push_back(t2=InterpolateSWMF);
#endif

}

//====================================================
//retrive the data file from SWMF
void PIC::ICES::retriveSWMFdata(const char *DataFile) {
  char cCurrentPath[_MAX_STRING_LENGTH_PIC_],command[_MAX_STRING_LENGTH_PIC_],initDirectory[_MAX_STRING_LENGTH_PIC_];

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
  if (PIC::Mesh::mesh.ThisThread==0) printf("Getting the MHD data from ICES.... \nThe current working directory is: %s\n", cCurrentPath);

  system("rm -f MHDRestart");

  sprintf(command,"ln -s %s/Data/%s/MHD/ MHDRestart",locationICES,DataFile);
  if (PIC::Mesh::mesh.ThisThread==0) printf("%s\n",command);
  system(command);

  if (PIC::Mesh::mesh.ThisThread==0) {
    sprintf(command,"%s/MHD/comet_mhd.exe",locationICES);
    printf("%s\n",command);
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

  MPI_Barrier(MPI_COMM_WORLD);
  if (PIC::Mesh::mesh.ThisThread==0) printf("ICES done\n");
}

//====================================================
//read and parse the data file from SWMF
void PIC::ICES::readSWMFdata(const double MeanIonMass) {
  cDataNodeSWMF dataSWMF;
  CiFileOperations ices;
  int status,idim;
  long int nd;
  char *offset;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_],*endptr;

  sprintf(str,"icesCellCenterCoordinates.thread=%i.MHD.dat",PIC::Mesh::mesh.ThisThread);
  ices.openfile(str);

  ices.GetInputStr(str,sizeof(str));
  ices.GetInputStr(str,sizeof(str));

  PIC::Mesh::mesh.resetNodeProcessedFlag();

  //read the data file
#if DIM == 3
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)*(2*_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_)*(2*_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_);
#elif DIM == 2
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)*(2*_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_);
#elif DIM == 1
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)
#endif

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    for (nd=0;nd<ndMax;nd++) {
      CenterNode=node->block->GetCenterNode(nd);

      if (CenterNode!=NULL) if (CenterNode->nodeDescriptor.nodeProcessedFlag==_AMR_FALSE_) {
        CenterNode->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;

        //the read section
        ices.GetInputStr(str,sizeof(str));

        ices.CutInputStr(str1,str);
        ices.CutInputStr(str1,str);
        ices.CutInputStr(str1,str);
        ices.CutInputStr(str1,str);

        status=strtol(str1,&endptr,10);

        if (status!=0) {
          ices.error("the point is not found");
          exit(__LINE__,__FILE__,"Error: the extracted point is not found");
        }

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
        offset=CenterNode->GetAssociatedDataBufferPointer();

        for (idim=0;idim<3;idim++) {
          *(idim+(double*)(offset+ElectricFieldOffset))=dataSWMF.E[idim];
          *(idim+(double*)(offset+MagneticFieldOffset))=dataSWMF.B[idim];
          *(idim+(double*)(offset+PlasmaBulkVelocityOffset))=dataSWMF.swVel[idim];
        }

        *(double*)(offset+PlasmaPressureOffset)=dataSWMF.swPressure;
        *(double*)(offset+PlasmaNumberDensityOffset)=dataSWMF.swNumberDensity;
        *(double*)(offset+PlasmaTemperatureOffset)=dataSWMF.swTemperature;


      }
    }
  }

  ices.closefile();
}

//====================================================
//output background parameters loaded with ICES from SWMF output
void PIC::ICES::PrintVariableListSWMF(FILE* fout,int DataSetNumber) {
  int idim;

  for (idim=0;idim<DIM;idim++) fprintf(fout,", \"E%i\"",idim);
  for (idim=0;idim<DIM;idim++) fprintf(fout,", \"B%i\"",idim);
  for (idim=0;idim<DIM;idim++) fprintf(fout,", \"PlasmaVel%i\"",idim);

  fprintf(fout,", \"Plasma Pressure\", \"Plasma Temperature\", \"Plasma number Desnity\"");
}

void PIC::ICES::PrintDataSWMF(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  cDataNodeSWMF dataSWMF;

  if (pipe->ThisThread==CenterNodeThread) {
    int idim;
    char *offset=CenterNode->GetAssociatedDataBufferPointer();

    dataSWMF.swNumberDensity=*((double*)(PlasmaNumberDensityOffset+offset));
    dataSWMF.swPressure=*((double*)(PlasmaPressureOffset+offset));
    dataSWMF.swTemperature=*((double*)(PlasmaTemperatureOffset+offset));

    for (idim=0;idim<3;idim++) {
      dataSWMF.B[idim]=*(idim+(double*)(MagneticFieldOffset+offset));
      dataSWMF.E[idim]=*(idim+(double*)(ElectricFieldOffset+offset));
      dataSWMF.swVel[idim]=*(idim+(double*)(PlasmaBulkVelocityOffset+offset));
    }
  }


  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv((char*)&dataSWMF,sizeof(dataSWMF),CenterNodeThread);

    fprintf(fout,"%e  %e  %e  ",dataSWMF.E[0],dataSWMF.E[1],dataSWMF.E[2]);
    fprintf(fout,"%e  %e  %e  ",dataSWMF.B[0],dataSWMF.B[1],dataSWMF.B[2]);
    fprintf(fout,"%e  %e  %e  ",dataSWMF.swVel[0],dataSWMF.swVel[1],dataSWMF.swVel[2]);

    fprintf(fout,"%e  %e  %e  ",dataSWMF.swPressure,dataSWMF.swTemperature,dataSWMF.swNumberDensity);
  }
  else pipe->send((char*)&dataSWMF,sizeof(dataSWMF));
}

void PIC::ICES::InterpolateSWMF(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim,i;
  double c;
  char *offset,*offsetCenterNode;

  #if _PIC_DEBUGGER_MODE_ ==  _PIC_DEBUGGER_MODE_ON_
  if (CenterNode->GetAssociatedDataBufferPointer()==NULL) exit(__LINE__,__FILE__,"Error: The associated data buffer is not initialized");
  #endif

  //init the buffer of 'CenterNode'
  offsetCenterNode=CenterNode->GetAssociatedDataBufferPointer();

  *((double*)(PlasmaNumberDensityOffset+offsetCenterNode))=0.0;
  *((double*)(PlasmaPressureOffset+offsetCenterNode))=0.0;
  *((double*)(PlasmaTemperatureOffset+offsetCenterNode))=0.0;

  for (idim=0;idim<3;idim++) {
    *(idim+(double*)(MagneticFieldOffset+offsetCenterNode))=0.0;
    *(idim+(double*)(ElectricFieldOffset+offsetCenterNode))=0.0;
    *(idim+(double*)(PlasmaBulkVelocityOffset+offsetCenterNode))=0.0;
  }

  //interpolate the sampled data
  for (i=0;i<nInterpolationCoeficients;i++) {
    c=InterpolationCoeficients[i];
    offset=InterpolationList[i]->GetAssociatedDataBufferPointer();

    *((double*)(PlasmaNumberDensityOffset+offsetCenterNode))+=c*(*((double*)(PlasmaNumberDensityOffset+offset)));
    *((double*)(PlasmaPressureOffset+offsetCenterNode))+=c*(*((double*)(PlasmaPressureOffset+offset)));
    *((double*)(PlasmaTemperatureOffset+offsetCenterNode))+=c*(*((double*)(PlasmaTemperatureOffset+offset)));

    for (idim=0;idim<3;idim++) {
      *(idim+(double*)(MagneticFieldOffset+offsetCenterNode))+=c*(*(idim+(double*)(MagneticFieldOffset+offset)));
      *(idim+(double*)(ElectricFieldOffset+offsetCenterNode))+=c*(*(idim+(double*)(ElectricFieldOffset+offset)));
      *(idim+(double*)(PlasmaBulkVelocityOffset+offsetCenterNode))+=c*(*(idim+(double*)(PlasmaBulkVelocityOffset+offset)));
    }

//===================  DEBUG =========================
/*
    cDataNodeSWMF dataSWMF;

    dataSWMF.swNumberDensity=*((double*)(PlasmaNumberDensityOffset+offset));
    dataSWMF.swPressure=*((double*)(PlasmaPressureOffset+offset));
    dataSWMF.swTemperature=*((double*)(PlasmaTemperatureOffset+offset));

    for (idim=0;idim<3;idim++) {
      dataSWMF.B[idim]=*(idim+(double*)(MagneticFieldOffset+offset));
      dataSWMF.E[idim]=*(idim+(double*)(ElectricFieldOffset+offset));
      dataSWMF.swVel[idim]=*(idim+(double*)(PlasmaBulkVelocityOffset+offset));
    }


    dataSWMF.swNumberDensity=*((double*)(PlasmaNumberDensityOffset+offsetCenterNode));
    dataSWMF.swPressure=*((double*)(PlasmaPressureOffset+offsetCenterNode));
    dataSWMF.swTemperature=*((double*)(PlasmaTemperatureOffset+offsetCenterNode));

    for (idim=0;idim<3;idim++) {
      dataSWMF.B[idim]=*(idim+(double*)(MagneticFieldOffset+offsetCenterNode));
      dataSWMF.E[idim]=*(idim+(double*)(ElectricFieldOffset+offsetCenterNode));
      dataSWMF.swVel[idim]=*(idim+(double*)(PlasmaBulkVelocityOffset+offsetCenterNode));
    }
*/
//=================  END DEBUG =======================


  }
}



//====================================================
void PIC::ICES::createCellCenterCoordinateList() {
  FILE *fout;
  char fname[_MAX_STRING_LENGTH_PIC_];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  double *x;
  int idim;

  //generate the total file name and open the file
  sprintf(fname,"icesCellCenterCoordinates.thread=%i",PIC::Mesh::mesh.ThisThread);
  fout=fopen(fname,"w");
  PIC::Mesh::mesh.resetNodeProcessedFlag();

  fprintf(fout,"#START\n");

  //output the cell's coordinates

#if DIM == 3
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)*(2*_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_)*(2*_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_);
#elif DIM == 2
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)*(2*_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_);
#elif DIM == 1
  const long int ndMax=(2*_GHOST_CELLS_X_+_BLOCK_CELLS_X_)
#endif

  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    for (nd=0;nd<ndMax;nd++) {
      CenterNode=node->block->GetCenterNode(nd);

      if (CenterNode!=NULL) if (CenterNode->nodeDescriptor.nodeProcessedFlag==_AMR_FALSE_) {
        CenterNode->nodeDescriptor.nodeProcessedFlag=_AMR_TRUE_;
        x=CenterNode->GetX();

        for (idim=0;idim<DIM;idim++) fprintf(fout,"%e ",x[idim]);
        fprintf(fout,"\n");
      }
    }
  }

  fclose(fout);
}
