
//$Id$
//functions that are used for interpolation of the dataset obtainedwith other models into AMPS' grid

/*
 * pic_datafile.cpp
 *
 *  Created on: May 5, 2015
 *      Author: vtenishe
 */


#include "pic.h"


//path to the location of the datafiles
char PIC::CPLR::DATAFILE::path[_MAX_STRING_LENGTH_PIC_]=".";

//the offset from the cell->AssociatedData()
int PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin=-1;
int PIC::CPLR::DATAFILE::nTotalBackgroundVariables=0;
bool PIC::CPLR::DATAFILE::Offset::InitFlag=false;

//Physical quantaties offsets that could be read and srored
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity={false,false,1,"\"Plasma number density\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity={false,false,3,"\"vPlasmaX\", \"vPlasmaY\", \"vPlasmaZ\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaTemperature={false,false,1,"\"Plasma temperature\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure={false,false,1,"\"Plasma pressure\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaElectronPressure={false,false,1,"\"Plasma electron pressure\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::MagneticField={false,false,3,"\"Bx\", \"By\", \"Bz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::ElectricField={false,false,3,"\"Ex\", \"Ey\", \"Ez\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::GradientMagneticField={false,false,3,"\"d|B|/dx\", \"d|B|/dy\", \"d|B|/dz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::AbsoluteValueMagneticField={false,false,1,"\"|B|\"",-1};

//load new data file
//IMPORTANT! The list of the data that are loaded has to be indicated before PIC::Init_BeforeParser
void PIC::CPLR::DATAFILE::ImportData(const char *fname) {

  //the particular  reader
  if (_PIC_COUPLER_DATAFILE_READER_MODE_==_PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_) {
    TECPLOT::ImportData(fname);
  }
  else exit(__LINE__,__FILE__,"Error: the option is unknown");

}

//initialize the data reading namespace
void PIC::CPLR::DATAFILE::Init() {

  //set the initialization flag and call the init procedures of the particular file reader
  Offset::InitFlag=true;

  if (_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__DATAFILE_) {
    switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
    case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
      TECPLOT::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__ARMS_:
      ARMS::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__ICES_:
      ICES::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__KAMELEON_:
      KAMELEON::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
      BATSRUS::Init();
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  }
  else exit(__LINE__,__FILE__,"Error: wrong option");

  //init the offset table and request memory
  CenterNodeAssociatedDataOffsetBegin=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

  //init the data offsets
  if (Offset::PlasmaNumberDensity.allocate==true) {
    Offset::PlasmaNumberDensity.active=true;
    Offset::PlasmaNumberDensity.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::PlasmaNumberDensity.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaNumberDensity.nVars;
  }


  if (Offset::PlasmaBulkVelocity.allocate==true) {
    Offset::PlasmaBulkVelocity.active=true;
    Offset::PlasmaBulkVelocity.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::PlasmaBulkVelocity.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaBulkVelocity.nVars;
  }


  if (Offset::PlasmaTemperature.allocate==true) {
    Offset::PlasmaTemperature.active=true;
    Offset::PlasmaTemperature.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::PlasmaTemperature.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaTemperature.nVars;
  }


  if (Offset::PlasmaIonPressure.allocate==true) {
    Offset::PlasmaIonPressure.active=true;
    Offset::PlasmaIonPressure.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::PlasmaIonPressure.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaIonPressure.nVars;
  }


  if (Offset::PlasmaElectronPressure.allocate==true) {
    Offset::PlasmaElectronPressure.active=true;
    Offset::PlasmaElectronPressure.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::PlasmaElectronPressure.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaElectronPressure.nVars;
  }

  if (Offset::MagneticField.allocate==true) {
    Offset::MagneticField.active=true;
    Offset::MagneticField.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::MagneticField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::MagneticField.nVars;
  }

  if (Offset::ElectricField.allocate==true) {
    Offset::ElectricField.active=true;
    Offset::ElectricField.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::ElectricField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::ElectricField.nVars;
  }

  if (Offset::GradientMagneticField.allocate==true) {
    Offset::GradientMagneticField.active=true;
    Offset::GradientMagneticField.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::GradientMagneticField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::GradientMagneticField.nVars;
  }

  if (Offset::AbsoluteValueMagneticField.allocate==true) {
    Offset::AbsoluteValueMagneticField.active=true;
    Offset::AbsoluteValueMagneticField.offset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;

    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=Offset::AbsoluteValueMagneticField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::AbsoluteValueMagneticField.nVars;
  }




  if (nTotalBackgroundVariables==0) {
    exit(__LINE__,__FILE__,"Error: no background variables will be loaded. The background variables to load has to be indicated before call of the PIC::Init_BeforeParser()");
  }

  //print out of the output file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);

}


//save and load the binary data saved in the AMPS data structure
bool PIC::CPLR::DATAFILE::BinaryFileExists(const char *fNameBase) {
  FILE *fData=NULL;
  char fname[400];

  sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",PIC::Mesh::mesh.getMeshSignature(),fNameBase);
  fData=fopen(fname,"r");

  if (fData!=NULL) {
    fclose(fData);
    return true;
  }

  return false;

}


void PIC::CPLR::DATAFILE::SaveBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static CMPI_channel pipe;
  static FILE *fout=NULL;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    char fname[400];
    sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",PIC::Mesh::mesh.getMeshSignature(),fNameBase);

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
            fwrite(offset+CenterNodeAssociatedDataOffsetBegin,nTotalBackgroundVariables*sizeof(double),1,fout);

            #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
            PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin),nTotalBackgroundVariables,__LINE__,__FILE__);
            #endif
          }
          else {
            pipe.send((double*)(offset+CenterNodeAssociatedDataOffsetBegin),nTotalBackgroundVariables);
          }
        }
        else {
          double data[nTotalBackgroundVariables];

          pipe.recv(data,nTotalBackgroundVariables,startNode->Thread);
          fwrite(data,nTotalBackgroundVariables*sizeof(double),1,fout);

          #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
          PIC::Debugger::CatchOutLimitValue(data,nTotalBackgroundVariables,__LINE__,__FILE__);
          #endif
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SaveBinaryFile(NULL,startNode->downNode[nDownNode]);
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


void PIC::CPLR::DATAFILE::LoadBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static FILE *fData=NULL;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    char fname[400];
    sprintf(fname,"amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",PIC::Mesh::mesh.getMeshSignature(),fNameBase);

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

  char savedLoadCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block==NULL) {
       //the block belongs to a other processor -> skip all data
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
        fread(&savedLoadCellFlag,sizeof(char),1,fData);

        if (savedLoadCellFlag==true) {
          //the cell data is saved -> skip it
          fseek(fData,nTotalBackgroundVariables*sizeof(double),SEEK_CUR);
        }
      }

    }
    else for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      fread(&savedLoadCellFlag,sizeof(char),1,fData);

      if (savedLoadCellFlag==true) {
        //determine whether the cell data needed to be read
        //locate the cell
        nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);

        if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) {
          offset=CenterNode->GetAssociatedDataBufferPointer();

          //read the data
          fread(offset+CenterNodeAssociatedDataOffsetBegin,nTotalBackgroundVariables*sizeof(double),1,fData);

          #if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
          PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin),nTotalBackgroundVariables,__LINE__,__FILE__);
          #endif
        }
        else {
          //the cell data is saved -> skip it
          fseek(fData,nTotalBackgroundVariables*sizeof(double),SEEK_CUR);
        }
      }

    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadBinaryFile(NULL,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    fclose(fData);
  }
}

//====================================================
//output background parameters loaded with ICES from SWMF output
void PIC::CPLR::DATAFILE::PrintVariableList(FILE* fout,int DataSetNumber) {
  int idim;


  if (Offset::PlasmaNumberDensity.active) fprintf(fout,", %s",Offset::PlasmaNumberDensity.VarList);
  if (Offset::PlasmaBulkVelocity.active) fprintf(fout,", %s",Offset::PlasmaBulkVelocity.VarList);
  if (Offset::PlasmaTemperature.active) fprintf(fout,", %s",Offset::PlasmaTemperature.VarList);
  if (Offset::PlasmaIonPressure.active) fprintf(fout,", %s",Offset::PlasmaIonPressure.VarList);
  if (Offset::PlasmaElectronPressure.active) fprintf(fout,", %s",Offset::PlasmaElectronPressure.VarList);
  if (Offset::MagneticField.active) fprintf(fout,", %s",Offset::MagneticField.VarList);
  if (Offset::ElectricField.active) fprintf(fout,", %s",Offset::ElectricField.VarList);
  if (Offset::GradientMagneticField.active) fprintf(fout,", %s",Offset::GradientMagneticField.VarList);
  if (Offset::AbsoluteValueMagneticField.active) fprintf(fout,", %s",Offset::AbsoluteValueMagneticField.VarList);
}

void PIC::CPLR::DATAFILE::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int iData,iCoefficient;
  double c,data[nTotalBackgroundVariables];
  double *offset;

  //init the interpolation buffer
  for (iData=0;iData<nTotalBackgroundVariables;iData++) data[iData]=0.0;

  //interpolation loop
  for (iCoefficient=0;iCoefficient<nInterpolationCoeficients;iCoefficient++) {
    c=InterpolationCoeficients[iCoefficient];
    offset=(double*)(CenterNodeAssociatedDataOffsetBegin+InterpolationList[iCoefficient]->GetAssociatedDataBufferPointer());

    for (iData=0;iData<nTotalBackgroundVariables;iData++) data[iData]+=c*offset[iData];
  }

  //copy the interpolated data
  memcpy(CenterNodeAssociatedDataOffsetBegin+CenterNode->GetAssociatedDataBufferPointer(),data,nTotalBackgroundVariables*sizeof(double));
}

void PIC::CPLR::DATAFILE::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double data[nTotalBackgroundVariables];
  int i;

  //get the interpolated data
  if (pipe->ThisThread==CenterNodeThread) {
    memcpy(data,CenterNodeAssociatedDataOffsetBegin+CenterNode->GetAssociatedDataBufferPointer(),nTotalBackgroundVariables*sizeof(double));
  }

  //send the data to the root processor if needed and print them into a file
  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(data,nTotalBackgroundVariables,CenterNodeThread);

    for (i=0;i<nTotalBackgroundVariables;i++) fprintf(fout,"%e ",data[i]);
  }
  else pipe->send(data,nTotalBackgroundVariables);
}

//====================================================
//print the ion flux at a sphere
void PIC::CPLR::DATAFILE::PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius) {
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
void PIC::CPLR::DATAFILE::EvaluateSurfaceIonFlux(double ShiftFactor) {
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


