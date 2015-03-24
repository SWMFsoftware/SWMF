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
char PIC::CPLR::DATAFILE::path[_MAX_STRING_LENGTH_PIC_]=".";
double PIC::CPLR::DATAFILE::ARMS::OUTPUT::TimeCoupleNext=0.0;

//read ARM's output file
void PIC::CPLR::DATAFILE::ARMS::OUTPUT::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  // current time
  static double time = -1.0;
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
    while(fin.eof()==false && (nX < 0 || nZ < 0) ){
      if(fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false){
	exit(__LINE__,__FILE__, "ERROR: the size of the grid couldn't be read.");
	break;
      }
      fin.CutInputStr(str1,str);
      
      if((strcmp("#TIME",str1)==0)){
	fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	time = strtod(str1, NULL);
      }
      else
	if((strcmp("#TIME_NEXT",str1)==0)){
	  fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	  PIC::CPLR::DATAFILE::ARMS::OUTPUT::TimeCoupleNext=strtod(str1, NULL);
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
	  if((strcmp("#N_E",str1)==0)){offset=n_;convert=pow(cm2m,-3)/ ElectronMass;}
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
		      if((strcmp("#B_RAD",str1)==0)){offset=b_;convert=1E+5;}
		      else//         B_THETA [Gauss]
			if((strcmp("#B_THETA",str1)==0)){offset=b_+1;convert=1E+5;}
			else//         B_PHI [Gauss]
			  if((strcmp("#B_PHI",str1)==0)){offset=b_+2;convert=1E+5;}
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
	    if(xCell < 0||zCell < 0||xCell >= nX-1||zCell >= nZ-1) continue;
	    
	    // interpolation weights
	    double wX = (pow(radpol,2)-pow(Xpos[xCell],2)) / (twodX*Xpos[xCell]+dX2);
	    double wZ = (x[2] - Zpos[zCell]) / dZ;
	    
	    //interpolate values
	    double DataInterp[nvar] = {0.0}, E[nvar] = {0.0};
	    for(int ivar=0;ivar<nvar;ivar++)for(int ii=0;ii<2;ii++) for(int jj=0;jj<2;jj++)
		  DataInterp[ivar] += 
		    Data[ivar][xCell+ii][zCell+jj] * ((1-wX)*(1-ii) + wX*ii) * ((1-wZ)*(1-jj) + wZ*jj);
	    {
	      // transform data to cartesian: 
	      // initial data has vectors' components in spherical coordinates
	      double cosPhi    = x[0] / radpol, sinPhi   = x[1]   / radpol;
	      double cosTheta  = x[2] / radius, sinTheta = radpol / radius;
	      double tmp;
	      // transform velocity vector
	      tmp = DataInterp[v_]*cosTheta+DataInterp[v_+1]*sinTheta;
	      DataInterp[v_]   = DataInterp[v_]*sinTheta-DataInterp[v_+1]*cosTheta;
	      DataInterp[v_+1] = DataInterp[v_]*sinPhi+DataInterp[v_+2]*cosPhi;
	      DataInterp[v_  ] = DataInterp[v_]*cosPhi-DataInterp[v_+2]*sinPhi;
	      DataInterp[v_+2] = tmp;
	      // transfrom magnetic field vector
	      tmp = DataInterp[b_]*cosTheta+DataInterp[b_+1]*sinTheta;
	      DataInterp[b_]   = DataInterp[b_]*sinTheta-DataInterp[b_+1]*cosTheta;
	      DataInterp[b_+1] = DataInterp[b_]*sinPhi+DataInterp[b_+2]*cosPhi;
	      DataInterp[b_  ] = DataInterp[b_]*cosPhi-DataInterp[b_+2]*sinPhi;
	      DataInterp[b_+2] = tmp;
	    }
	    
	    
	    //locate the cell
	    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
	    if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
	    offset=CenterNode->GetAssociatedDataBufferPointer();
	    
	    //save the interpolated values
	    for (int idim=0;idim<3;idim++) {
	      *(idim+(double*)(offset+MagneticFieldOffset))     =DataInterp[b_+idim];
	      *(idim+(double*)(offset+PlasmaBulkVelocityOffset))=DataInterp[v_+idim];
	    }
	    
	    *((double*)(offset+PlasmaPressureOffset))     =DataInterp[p_];
	    *((double*)(offset+PlasmaNumberDensityOffset))=DataInterp[n_];
	    *((double*)(offset+PlasmaTemperatureOffset))  =DataInterp[t_];
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
    for(int ivar = 0; ivar < nvar; ivar++){
      for(int iX = 0; iX < nX; iX++){
	delete [] Data[ivar][iX];
      }
      delete [] Data[ivar];
    }
    delete [] Data;
  }
}
