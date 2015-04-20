/*
 * pic_arms.cpp
 *
 *  Created on: Mar 18, 2015
 *      Author: vtenishe
 */

//$Id$
//reading of ARMS' output file

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "pic.h"

//path to the location of the data files
char   PIC::CPLR::DATAFILE::path[_MAX_STRING_LENGTH_PIC_]=".";
double PIC::CPLR::DATAFILE::ARMS::OUTPUT::TimeCurrent   =-1.0;
double PIC::CPLR::DATAFILE::ARMS::OUTPUT::TimeCoupleNext=-1.0;
int    PIC::CPLR::DATAFILE::ARMS::OUTPUT::GradientMagneticFieldOffset=-1;
int    PIC::CPLR::DATAFILE::ARMS::OUTPUT::AbsoluteValueMagneticFieldOffset=-1;



void PIC::CPLR::DATAFILE::ARMS::OUTPUT::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,",\"|B|\",\"d|B|/dx\",\"d|B|/dy\",\"d|B|/dz\"");
}

void PIC::CPLR::DATAFILE::ARMS::OUTPUT::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  double dB[3]={0.0,0.0,0.0},B=0.0;
  int i,idim;
  char *SamplingBuffer;

  for (i=0;i<nInterpolationCoeficients;i++) {

    for (idim=0,SamplingBuffer=InterpolationList[i]->GetAssociatedDataBufferPointer()+GradientMagneticFieldOffset;idim<3;idim++) dB[idim]+=(*((double*)(SamplingBuffer+idim*sizeof(double))))*InterpolationCoeficients[i];

    B+=(*((double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+AbsoluteValueMagneticFieldOffset)))*InterpolationCoeficients[i];
  }

  memcpy(CenterNode->GetAssociatedDataBufferPointer()+GradientMagneticFieldOffset,dB,3*sizeof(double));
  memcpy(CenterNode->GetAssociatedDataBufferPointer()+AbsoluteValueMagneticFieldOffset,&B,sizeof(double));
}

void PIC::CPLR::DATAFILE::ARMS::OUTPUT::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  int idim;
  double t;

  // |B|
  if (pipe->ThisThread==CenterNodeThread) {
    t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+AbsoluteValueMagneticFieldOffset));
  }

  if (pipe->ThisThread==0) {
    if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

    fprintf(fout,"%e ",t);
  }
  else pipe->send(t);

  //Magnetic Field
  for (idim=0;idim<3;idim++) {
    if (pipe->ThisThread==CenterNodeThread) {
      t= *((double*)(CenterNode->GetAssociatedDataBufferPointer()+GradientMagneticFieldOffset+idim*sizeof(double)));
    }

    if (pipe->ThisThread==0) {
      if (CenterNodeThread!=0) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}


// initialize reading ARMS data: set additional offsets
void PIC::CPLR::DATAFILE::ARMS::OUTPUT::Init(){

  if(AbsoluteValueMagneticFieldOffset==-1){
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    AbsoluteValueMagneticFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=sizeof(double);
    TotalAssociatedDataLength+=sizeof(double);
  }

  if(GradientMagneticFieldOffset==-1){
    if (AssociatedDataOffset==-1) AssociatedDataOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    GradientMagneticFieldOffset=PIC::Mesh::cDataCenterNode::totalAssociatedDataLength;
    PIC::Mesh::cDataCenterNode::totalAssociatedDataLength+=3*sizeof(double);
    TotalAssociatedDataLength+=3*sizeof(double);
  }

  //print out of the output file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);


}

//read ARM's output file
void PIC::CPLR::DATAFILE::ARMS::OUTPUT::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  // size of the uniform grid
  static int nX = -1, nZ = -1;
  // container for the data and grid
  const int nvar = 9;
  const int b_ = 0, v_ = 3, n_ = 6, t_ = 7, p_ = 8; 
  static double *Xpos, *Zpos;
  static double ***Data;

  if (startNode==PIC::Mesh::mesh.rootTree) {
    // read data file using class CiFileOperations (see src/general/ifileopr.h)
    CiFileOperations fin;
    char fullname[_MAX_STRING_LENGTH_PIC_];
    sprintf(fullname,"%s/%s",PIC::CPLR::DATAFILE::path,fname);
    fin.openfile(fullname);
    
    // start reading: read essential constants time & size of the grid
    char str[_MAX_STRING_LENGTH_PIC_],str1[_MAX_STRING_LENGTH_PIC_];
    // reset parameters
    nX = -1; nZ=-1; 
    TimeCurrent   =-1.0;
    TimeCoupleNext=-1.0;
    while(fin.eof()==false && (nX < 0 || nZ < 0|| TimeCurrent < 0 || TimeCoupleNext < 0) ){
      if(fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false){
	exit(__LINE__,__FILE__, "ERROR: the size of the grid couldn't be read.");
	break;
      }
      fin.CutInputStr(str1,str);
      
      if((strcmp("#TIME",str1)==0)){
	fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	TimeCurrent = strtod(str1, NULL);
      }
      else
	if((strcmp("#TIME_NEXT",str1)==0)){
	  fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	  TimeCoupleNext=strtod(str1, NULL);
	  if(TimeCoupleNext < 0)
	    exit(__LINE__,__FILE__,"Reached the last ARMS data file; exit");
	}
	else
	  if((strcmp("#NX",str1)==0)){
	    fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	    nX = strtol(str1, NULL, 10);
	  }
	  else
	    if((strcmp("#NZ",str1)==0)){
	      fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	      nZ = strtol(str1, NULL, 10);
	    }
    }
    // allocate container for the data and the grid
    Data = new double** [nvar];
    for(int ivar = 0; ivar < nvar; ivar++){
      Data[ivar] = new double* [nX];
      for(int iX = 0; iX < nX; iX++){
	Data[ivar][iX] = new double [nZ];
      }
    }
    Xpos = new double [nX];
    Zpos = new double [nZ];
    // keep reading the actual data
    while(fin.eof()==false){
      if(fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false){
	break;
      }
      fin.CutInputStr(str1,str);
      
      if((strcmp("#XPOS",str1)==0)){
	for(int iX = 0; iX < nX; iX++){
	  fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	  Xpos[iX] = strtod(str1, NULL)*cm2m;
	}
      }
      else
	if((strcmp("#ZPOS",str1)==0)){
	  for(int iZ = 0; iZ < nZ; iZ++){
	    fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	    Zpos[iZ] = strtod(str1, NULL)*cm2m;
	  }
	}
	else{
	  int offset = -1; 
	  double convert = -1.0;
	  //           N_E is electron mass density [cm^-3]
	  if((strcmp("#N_E",str1)==0)){offset=n_;convert=pow(cm2m,-3);}
	  else//         PRE is thermal pressure [dyn/cm^-2]
	    if((strcmp("#PRE",str1)==0)){offset=p_;convert=0.1;}
	    else//         TEM is temperature [K]
	      if((strcmp("#TEM",str1)==0)){offset=t_;convert=1.0;}
	      else//         VEL_RAD [cm/s]
		if((strcmp("#VEL_RAD",str1)==0)){offset=v_;convert=cm2m;}
		else//         VEL_THETA [cm/s]
		  if((strcmp("#VEL_THETA",str1)==0)){offset=v_+1;convert=cm2m;}
		  else//         VEL_PHI [cm/s]
		    if((strcmp("#VEL_PHI",str1)==0)){offset=v_+2;convert=cm2m;}
		    else//         B_RAD [Gauss]
		      if((strcmp("#B_RAD",str1)==0)){offset=b_;convert=1E-4;}
		      else//         B_THETA [Gauss]
			if((strcmp("#B_THETA",str1)==0)){offset=b_+1;convert=1E-4;}
			else//         B_PHI [Gauss]
			  if((strcmp("#B_PHI",str1)==0)){offset=b_+2;convert=1E-4;}
	  if(offset >= 0)
	    for(int iZ = 0; iZ < nZ; iZ++){
	      for(int iX = 0; iX < nX; iX++){
		fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
		Data[offset][iX][iZ] = strtod(str1, NULL) * convert;
	      }
	    }
	}
    }
  }
  // now the data has been read

  // compute gradient of the magnetic field
  // allocate container for the magnetic field data
  double **absB, **dBdR, **dBdZ;
  absB = new double* [nX];
  dBdR = new double* [nX-1];
  dBdZ = new double* [nX];
  for(int iX = 0; iX < nX; iX++){
    absB[iX] = new double [nZ];
    if(iX < nX-1)  dBdR[iX] = new double [nZ];
    dBdZ[iX] = new double [nZ-1];
  }

  {    
    double dX = Xpos[1] - Xpos[0];
    double dZ = Zpos[1] - Zpos[0];
    // absolute value of B
    for(int iX = 0; iX < nX; iX++)
      for(int iZ = 0; iZ < nZ; iZ++)
	absB[iX][iZ] = pow(Data[b_  ][iX][iZ]*Data[b_  ][iX][iZ]+
			   Data[b_+1][iX][iZ]*Data[b_+1][iX][iZ]+
			   Data[b_+2][iX][iZ]*Data[b_+2][iX][iZ],0.5);
    // components of gradient
    for(int iX = 0; iX < nX; iX++)
      for(int iZ = 0; iZ < nZ; iZ++){
	if(iX < nX-1) dBdR[iX][iZ] = (absB[iX+1][iZ]-absB[iX][iZ]) / dX;
	if(iZ < nZ-1) dBdZ[iX][iZ] = (absB[iX][iZ+1]-absB[iX][iZ]) / dZ;
	
      }
  }
  
  // perform the interpolation
  {
    int nd;
    PIC::Mesh::cDataCenterNode *CenterNode;
    char *offset;
    double dX = Xpos[1] - Xpos[0], dX2 = dX*dX, twodX = 2.0 * dX;
    double dZ = Zpos[1] - Zpos[0];
    
    const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
    const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
    const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;
    
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      double *xNodeMin=startNode->xmin;
      double *xNodeMax=startNode->xmax;
      double x[3],T,n,p;
      
      for (int k=kMin;k<=kMax;k++)for(int j=jMin;j<=jMax;j++)for (int i=iMin;i<=iMax;i++) {
	    //the interpolation location
	    x[0]=xNodeMin[0]+(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_*(0.5+i);
	    x[1]=xNodeMin[1]+(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_*(0.5+j);
	    x[2]=xNodeMin[2]+(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_*(0.5+k);
	    
	    // locate corresponding cell in the Data container
	    // note: initial Data is 2.5D, i.e. has cylindrical symmetry
	    double radpol2 = x[0]*x[0] + x[1]*x[1]+1E-15, radpol = pow(radpol2, 0.5);
	    double radius  = pow(radpol2 + x[2]*x[2],0.5);
	    int xCell = floor( (radpol - Xpos[0]) / dX);
	    int zCell = floor( (x[2]   - Zpos[0]) / dZ);
	    if(xCell <= 0||zCell <= 0||xCell >= nX-2||zCell >= nZ-2) continue;
	    
	    // interpolation weights
	    double wX = (pow(radpol,2)-pow(Xpos[xCell],2)) / (twodX*Xpos[xCell]+dX2);
	    // threshold weight for interpolating magnetic field gradient
	    double wXthresh = 0.25 + .5 / (2.0 + dX/Xpos[xCell]);
	    double wZ = (x[2] - Zpos[zCell]) / dZ;
	    
	    //interpolate values
	    double DataInterp[nvar] = {0.0}, GrdBInterp[3]= {0.0}, absBInterp=0.0;
	    for(int ivar=0;ivar<nvar;ivar++)for(int ii=0;ii<2;ii++) for(int jj=0;jj<2;jj++){
		  DataInterp[ivar] += 
		    Data[ivar][xCell+ii][zCell+jj] * ((1-wX)*(1-ii) + wX*ii) * ((1-wZ)*(1-jj) + wZ*jj);
		  absBInterp    += 
		    absB[xCell+ii][zCell+jj]       * ((1-wX)*(1-ii) + wX*ii) * ((1-wZ)*(1-jj) + wZ*jj);
		  GrdBInterp[0] += 
		    (wX > wXthresh) ? 
		    dBdR[xCell+ii][zCell+jj] * ((1-wX+wXthresh)*(1-ii) + (wX-wXthresh)*ii)* ((1-wZ)*(1-jj) + wZ*jj):
		    dBdR[xCell-ii][zCell+jj] * ((0.5-wX+wXthresh)*(1-ii) + (0.5+wX-wXthresh)*ii)* ((1-wZ)*(1-jj) + wZ*jj);
		  GrdBInterp[2] += 
		    (wZ > 0.5) ? 
		    dBdZ[xCell+ii][zCell+jj] * ((1-wX)*(1-ii) + wX*ii)* ((1.5-wZ)*(1-jj) + (wZ-0.5)*jj):
		    dBdZ[xCell+ii][zCell-jj] * ((1-wX)*(1-ii) + wX*ii)* ((0.5+wZ)*(1-jj) + (0.5-wZ)*jj);
		}
	    {
	      // transform data to cartesian: 
	      // initial data has vectors' components in spherical coordinates
	      double cosPhi    = x[0] / radpol, sinPhi   = x[1]   / radpol;
	      double cosTheta  = x[2] / radius, sinTheta = radpol / radius;
	      double tmp;
	      // transform velocity vector
	      tmp = DataInterp[v_]*cosTheta-DataInterp[v_+1]*sinTheta;
	      DataInterp[v_]   = DataInterp[v_]*sinTheta+DataInterp[v_+1]*cosTheta;
	      DataInterp[v_+1] = DataInterp[v_]*sinPhi+DataInterp[v_+2]*cosPhi;
	      DataInterp[v_  ] = DataInterp[v_]*cosPhi-DataInterp[v_+2]*sinPhi;
	      DataInterp[v_+2] = tmp;
	      // transfrom magnetic field vector
	      tmp = DataInterp[b_]*cosTheta-DataInterp[b_+1]*sinTheta;
	      DataInterp[b_]   = DataInterp[b_]*sinTheta+DataInterp[b_+1]*cosTheta;
	      DataInterp[b_+1] = DataInterp[b_]*sinPhi+DataInterp[b_+2]*cosPhi;
	      DataInterp[b_  ] = DataInterp[b_]*cosPhi-DataInterp[b_+2]*sinPhi;
	      DataInterp[b_+2] = tmp;

	      // gradient is computed in cylindrical coordinates
	      // transform gradient of magnetic field
	      GrdBInterp[1] = GrdBInterp[b_]*sinPhi;
	      GrdBInterp[0] = GrdBInterp[b_]*cosPhi;
	    }

	    
	    
	    //locate the cell
	    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
	    if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
	    offset=CenterNode->GetAssociatedDataBufferPointer();
	    
	    //save the interpolated values
	    for (int idim=0;idim<3;idim++) {
	      *(idim+(double*)(offset+MagneticFieldOffset))        =DataInterp[b_+idim];
	      *(idim+(double*)(offset+GradientMagneticFieldOffset))=GrdBInterp[idim];
	      *(idim+(double*)(offset+PlasmaBulkVelocityOffset))   =DataInterp[v_+idim];
	    }
	    // E = -VxB
	    *(0+(double*)(offset+ElectricFieldOffset)) = 
	      DataInterp[b_+1]*DataInterp[v_+2] - DataInterp[b_+2]*DataInterp[v_+1];
	    *(1+(double*)(offset+ElectricFieldOffset)) = 
	      DataInterp[b_+2]*DataInterp[v_+0] - DataInterp[b_+0]*DataInterp[v_+2];
	    *(2+(double*)(offset+ElectricFieldOffset)) = 
	      DataInterp[b_+0]*DataInterp[v_+1] - DataInterp[b_+1]*DataInterp[v_+0];
    
	    *((double*)(offset+PlasmaPressureOffset))            =DataInterp[p_];
	    *((double*)(offset+PlasmaNumberDensityOffset))       =DataInterp[n_];
	    *((double*)(offset+PlasmaTemperatureOffset))         =DataInterp[t_];
	    *((double*)(offset+AbsoluteValueMagneticFieldOffset))=absBInterp;
	  }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadDataFile(fname,startNode->downNode[nDownNode]);
    }
  }
  

  if (startNode==PIC::Mesh::mesh.rootTree) {
    // deallocate data containers
    delete [] Xpos;
    delete [] Zpos;  
    for(int iX = 0; iX < nX; iX++){
      delete [] absB[iX];
      if(iX < nX-1) delete [] dBdR[iX];
      delete [] dBdZ[iX];
    }
    delete [] absB;
    delete [] dBdR;
    delete [] dBdZ;
    for(int ivar = 0; ivar < nvar; ivar++){
      for(int iX = 0; iX < nX; iX++){
	delete [] Data[ivar][iX];
      }
      delete [] Data[ivar];
    }
    delete [] Data;
  }
}
