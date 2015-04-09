
//$Id$

/*
 * pic_tecplot.cpp
 *
 *  Created on: Apr 8, 2015
 *      Author: vtenishe
 */

#include "pic.h"


//the procedures that created the script and extract data from tecplot data files
double PIC::CPLR::DATAFILE::TECPLOT::xDataMin[3]={0.0,0.0,0.0},PIC::CPLR::DATAFILE::TECPLOT::xDataMax[3]={0.0,0.0,0.0},PIC::CPLR::DATAFILE::TECPLOT::UnitLength=1.0;
int PIC::CPLR::DATAFILE::TECPLOT::maxScriptPointNumber=15000;
int PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=0;

PIC::CPLR::DATAFILE::TECPLOT::cLoadedVariableData PIC::CPLR::DATAFILE::TECPLOT::Velocity,PIC::CPLR::DATAFILE::TECPLOT::Pressure,PIC::CPLR::DATAFILE::TECPLOT::MagneticField,PIC::CPLR::DATAFILE::TECPLOT::Density;

//set the limits of the domain in the data file
void PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimits(double *xmin,double *xmax) {
  for (int idim=0;idim<3;idim++) xDataMin[idim]=xmin[idim],xDataMax[idim]=xmax[idim];
}

//calcualte the number of the points that need to be extracted and save them into a file
int PIC::CPLR::DATAFILE::TECPLOT::CountInterpolatedPointNumber(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static int nExtractedPoints=0;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    nExtractedPoints=0;
  }

  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double x[3];

    if (startNode->Thread==PIC::ThisThread) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i))/UnitLength;
      x[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j))/UnitLength;
      x[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k))/UnitLength;

      //locate the cell
      if (startNode->block==NULL) continue;

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if (startNode->block->GetCenterNode(nd)==NULL) continue;

      //count the point if it is within the limits of the domain
      if ( (xDataMin[0]<x[0])&&(x[0]<xDataMax[0]) && (xDataMin[1]<x[1])&&(x[1]<xDataMax[1]) && (xDataMin[2]<x[2])&&(x[2]<xDataMax[2]) ) {
        nExtractedPoints++;
      }

    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) CountInterpolatedPointNumber(startNode->downNode[nDownNode]);
  }

  return nExtractedPoints;
}

//create the TECPLOT Script, and extract and load the data
int PIC::CPLR::DATAFILE::TECPLOT::CreateScript(const char *ScriptBaseName,const char* DataFileTECPLOT,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static char ScriptName[_MAX_STRING_LENGTH_PIC_];
  static FILE *fScript=NULL;
  static int InterpolatedDataFileCounter=0;
  static int nDataPointLeft=0,nScriptPrintedPoints=0;


  //init the datapoint number
  if (startNode==PIC::Mesh::mesh.rootTree) {
    nDataPointLeft=CountInterpolatedPointNumber(startNode);
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    double *xNodeMin=startNode->xmin;
    double *xNodeMax=startNode->xmax;
    double x[3];

    if (startNode->Thread==PIC::ThisThread) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i))/UnitLength;
      x[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j))/UnitLength;
      x[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k))/UnitLength;

      //locate the cell
      if (startNode->block==NULL) continue;

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if (startNode->block->GetCenterNode(nd)==NULL) continue;

      //save the point if it is within the limits of the domain
      if ( (xDataMin[0]<x[0])&&(x[0]<xDataMax[0]) && (xDataMin[1]<x[1])&&(x[1]<xDataMax[1]) && (xDataMin[2]<x[2])&&(x[2]<xDataMax[2]) ) {
        //open new script file if needed
        if ((fScript==NULL) || (nScriptPrintedPoints==maxScriptPointNumber)) {
          nScriptPrintedPoints=0;

          if (fScript==NULL) {
            sprintf(ScriptName,"%s.thread=%i.mcr",ScriptBaseName,PIC::ThisThread);
            fScript=fopen(ScriptName,"w");
            fprintf(fScript,"#!MC 1400\n");

            //read the data set
            fprintf(fScript,"$!READDATASET  '\"%s\" '\n",DataFileTECPLOT);
            fprintf(fScript,"  READDATAOPTION = NEW\n");
            fprintf(fScript,"  RESETSTYLE = YES\n");
            fprintf(fScript,"  INCLUDETEXT = NO\n");
            fprintf(fScript,"  INCLUDEGEOM = NO\n");
            fprintf(fScript,"  INCLUDECUSTOMLABELS = NO\n");
            fprintf(fScript,"  VARLOADMODE = BYNAME\n");
            fprintf(fScript,"  ASSIGNSTRANDIDS = YES\n");
            fprintf(fScript,"  INITIALPLOTTYPE = CARTESIAN3D\n");
            fprintf(fScript,"  VARNAMELIST = '\"X [R]\" \"Y [R]\" \"Z [R]\" \"`r [amu/cm^3]\" \"U_x [km/s]\" \"U_y [km/s]\" \"U_z [km/s]\" \"B_x [nT]\" \"B_y [nT]\" \"B_z [nT]\" \"p [nPa]\" \"J_x [`mA/m^2]\" \"J_y [`mA/m^2]\" \"J_z [`mA/m^2]\"'\n");
          }

          fprintf(fScript,"$!EXTRACTFROMPOLYLINE\n");
          fprintf(fScript,"  EXTRACTTHROUGHVOLUME = YES\n");
          fprintf(fScript,"  EXTRACTLINEPOINTSONLY = YES\n");
          fprintf(fScript,"  INCLUDEDISTANCEVAR = NO\n");
          fprintf(fScript,"  NUMPTS = 20\n");
          fprintf(fScript,"  EXTRACTTOFILE = YES\n");
          fprintf(fScript,"  FNAME = '%s.thread=%i.%i.dat'\n",ScriptBaseName,PIC::ThisThread,InterpolatedDataFileCounter);
          fprintf(fScript,"  RAWDATA\n");

          InterpolatedDataFileCounter++;

          //the number of the pointed to be interpolated
          fprintf(fScript,"%i\n",std::min(maxScriptPointNumber,nDataPointLeft));
        }

        fprintf(fScript,"%e %e %e\n",x[0],x[1],x[2]);

        nScriptPrintedPoints++;
        nDataPointLeft--;
      }

    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) CreateScript(ScriptBaseName,DataFileTECPLOT,startNode->downNode[nDownNode]);
  }


  //close the script file
  if (startNode==PIC::Mesh::mesh.rootTree) {
    if (fScript!=NULL) fclose(fScript);
  }

  return InterpolatedDataFileCounter;
}


//read the data files
void PIC::CPLR::DATAFILE::TECPLOT::LoadDataFile(const char *fname,int nTotalOutputFiles,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static char DataFileName[_MAX_STRING_LENGTH_PIC_],str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_];
  static double *data=NULL;
  static int nLoadedDataFile=0;
  static int nLoadedDataPoints=0;
  static CiFileOperations fin;
  static bool FirstPassFlag=true;


  if (startNode==PIC::Mesh::mesh.rootTree) {
    data=new double [nTotalVarlablesTECPLOT];
    nLoadedDataFile=0,nLoadedDataPoints=0;
    FirstPassFlag=true;
  }

  //loop through all points
  //create the list of the points
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
    double x[3];

    if (startNode->Thread==PIC::ThisThread) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      //the interpolation location
      x[0]=(xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i))/UnitLength;
      x[1]=(xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j))/UnitLength;
      x[2]=(xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k))/UnitLength;

      //locate the cell
      if (startNode->block==NULL) continue;

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
      offset=CenterNode->GetAssociatedDataBufferPointer();

      //save the point if it is within the limits of the domain
      if ( (xDataMin[0]<x[0])&&(x[0]<xDataMax[0]) && (xDataMin[1]<x[1])&&(x[1]<xDataMax[1]) && (xDataMin[2]<x[2])&&(x[2]<xDataMax[2]) ) {
        //open new script file if needed
        if ((FirstPassFlag==true) || (nLoadedDataPoints==maxScriptPointNumber)) {
          nLoadedDataPoints=0;

          if (FirstPassFlag==false) {
            fin.closefile();
          }

          FirstPassFlag=false;
          sprintf(DataFileName,"%s.thread=%i.%i.dat",fname,PIC::ThisThread,nLoadedDataFile);
          fin.openfile(DataFileName);

          //skip the first 2+nTotalVarlablesTECPLOT lines
          for (int i=0;i<2+nTotalVarlablesTECPLOT;i++) fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);

          nLoadedDataFile++;
        }

        //read data points
        fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);

        for (int i=0;i<nTotalVarlablesTECPLOT;i++) {
          fin.CutInputStr(str1,str);
          data[i]=strtod(str1, NULL);
        }

        //save the data on the AMPS data buffers
        //the order of the state vector: number density, temperature
        *((double*)(offset+PlasmaNumberDensityOffset))=data[Density.offset]*Density.ScaleFactor;
        *((double*)(offset+PlasmaTemperatureOffset))=(data[Density.offset]>0.0) ? data[Pressure.offset]*Pressure.ScaleFactor/(Kbol*data[Density.offset]*Density.ScaleFactor) : 0.0;

        //get pressure
        *((double*)(offset+PlasmaPressureOffset))=data[Pressure.offset]*Pressure.ScaleFactor;

        //bulk velocity and magnetic field
        for (idim=0;idim<3;idim++) {
          *((double*)(offset+MagneticFieldOffset+idim*sizeof(double)))=data[idim+MagneticField.offset]*MagneticField.ScaleFactor;
          *((double*)(offset+PlasmaBulkVelocityOffset+idim*sizeof(double)))=data[idim+Velocity.offset]*Velocity.ScaleFactor;
        }

        //calculate the electric field
        double *E,*B,*v;

        v=(double*)(offset+PlasmaBulkVelocityOffset);
        B=(double*)(offset+MagneticFieldOffset);
        E=(double*)(offset+ElectricFieldOffset);

        E[0]=-(v[1]*B[2]-B[1]*v[2]);
        E[1]=+(v[0]*B[2]-B[0]*v[2]);
        E[2]=-(v[0]*B[1]-B[0]*v[1]);



        //increment the point counter
        nLoadedDataPoints++;
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
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(fname,nTotalOutputFiles,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    delete [] data;
    nLoadedDataFile=0,nLoadedDataPoints=0;
    fin.closefile();
  }

}


//save and load the binary data saved in the AMPS data structure
void PIC::CPLR::DATAFILE::SaveBinaryFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static CMPI_channel pipe;
  static FILE *fout=NULL;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    pipe.init(1000000);

    if (PIC::Mesh::mesh.ThisThread==0) {
      pipe.openRecvAll();
      fout=fopen(fname,"w");
    }
    else pipe.openSend(0);
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char SendCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((PIC::ThisThread==0)||(startNode->Thread==PIC::ThisThread)) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      SendCellFlag=true;

      //determine whether the cell data neede to be saved
      if (startNode->Thread==PIC::ThisThread) {
        //locate the cell
        if (startNode->block==NULL) SendCellFlag=false,offset=NULL;

        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
        if (SendCellFlag==true) {
          if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) offset=CenterNode->GetAssociatedDataBufferPointer();
          else SendCellFlag=false,offset=NULL;
        }

        if (startNode->Thread!=0) pipe.send(SendCellFlag);
      }
      else {
        pipe.recv(SendCellFlag,startNode->Thread);
      }

      //save the cell data saving flag
      if (PIC::ThisThread==0) fwrite(&SendCellFlag,sizeof(char),1,fout);

      //save the cell data
      if (SendCellFlag==true) {
        if (startNode->Thread==PIC::ThisThread) {
          if (startNode->Thread==0) {
            fwrite(offset+PlasmaNumberDensityOffset,sizeof(double),1,fout);
            fwrite(offset+PlasmaTemperatureOffset,sizeof(double),1,fout);
            fwrite(offset+PlasmaPressureOffset,sizeof(double),1,fout);
            fwrite(offset+MagneticFieldOffset,3*sizeof(double),1,fout);
            fwrite(offset+ElectricFieldOffset,3*sizeof(double),1,fout);
            fwrite(offset+PlasmaBulkVelocityOffset,3*sizeof(double),1,fout);
          }
          else {
            pipe.send((double*)(offset+PlasmaNumberDensityOffset),1);
            pipe.send((double*)(offset+PlasmaTemperatureOffset),1);
            pipe.send((double*)(offset+PlasmaPressureOffset),1);
            pipe.send((double*)(offset+MagneticFieldOffset),3);
            pipe.send((double*)(offset+ElectricFieldOffset),3);
            pipe.send((double*)(offset+PlasmaBulkVelocityOffset),3);
          }
        }
        else {
          double data[3];

          pipe.recv(data,1,startNode->Thread);
          fwrite(data,sizeof(double),1,fout);

          pipe.recv(data,1,startNode->Thread);
          fwrite(data,sizeof(double),1,fout);

          pipe.recv(data,1,startNode->Thread);
          fwrite(data,sizeof(double),1,fout);

          pipe.recv(data,3,startNode->Thread);
          fwrite(data,sizeof(double),3,fout);

          pipe.recv(data,3,startNode->Thread);
          fwrite(data,sizeof(double),3,fout);

          pipe.recv(data,3,startNode->Thread);
          fwrite(data,sizeof(double),3,fout);
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SaveBinaryFile(fname,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    if (PIC::Mesh::mesh.ThisThread==0) {
      pipe.closeRecvAll();
      fclose(fout);
    }
    else pipe.closeSend();

    pipe.remove();
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }

}


void PIC::CPLR::DATAFILE::LoadBinaryFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static FILE *fData=NULL;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    fData=fopen(fname,"r");
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char LoadCellFlag,savedLoadCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->Thread!=PIC::ThisThread) {
       //the block belongs to a other processor -> skip all data
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
        fread(&savedLoadCellFlag,sizeof(char),1,fData);

        if (savedLoadCellFlag==true) {
          //the cell data is saved -> skip it
          fseek(fData,(1+1+1+3+3+3)*sizeof(double),SEEK_CUR);
        }
      }

    }
    else for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      LoadCellFlag=true;

      //determine whether the cell data needed to be read
      //locate the cell
      if (startNode->block==NULL) LoadCellFlag=false,offset=NULL;

      nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      if (LoadCellFlag==true) {
        if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) offset=CenterNode->GetAssociatedDataBufferPointer();
        else LoadCellFlag=false,offset=NULL;
      }

      //read and compare the saved 'cell read' flag
      fread(&savedLoadCellFlag,sizeof(char),1,fData);
      if (savedLoadCellFlag!=LoadCellFlag) exit(__LINE__,__FILE__,"LoadCellFlag is not consistent with the saved one. Someing is wrong with saving/reading of the binary data");

      if (LoadCellFlag==true) {
        //read the data
        fread(offset+PlasmaNumberDensityOffset,sizeof(double),1,fData);
        fread(offset+PlasmaTemperatureOffset,sizeof(double),1,fData);
        fread(offset+PlasmaPressureOffset,sizeof(double),1,fData);
        fread(offset+MagneticFieldOffset,3*sizeof(double),1,fData);
        fread(offset+ElectricFieldOffset,3*sizeof(double),1,fData);
        fread(offset+PlasmaBulkVelocityOffset,3*sizeof(double),1,fData);
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadBinaryFile(fname,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    fclose(fData);
  }

}



