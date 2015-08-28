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

// initialize reading ARMS data: set additional offsets
void PIC::CPLR::DATAFILE::ARMS::Init() {

  //reserve place for the interpolated data
  PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::ElectricField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.allocate=true;

  PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate=true;
}

//read ARM's output file
void PIC::CPLR::DATAFILE::ARMS::LoadDataFile(const char *fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  // size of the uniform grid
  static int nX = -1, nZ = -1;
  // container for the data and grid
  const int nvar = 22;
  //original data
  const int b_ = 0, v_ = 3, n_ = 6, t_ = 7, p_ = 8; 
  //additional data (gradients come in the END)
  const int Db_ = 10;
  // vectors in the data
  const int nvec       = 6;
  const int vecs[nvec] = {b_, v_, Db_, Db_+3, Db_+6};
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
    MULTIFILE::TimeCurrent   =-1.0;
    MULTIFILE::TimeCoupleNext=-1.0;
    while(fin.eof()==false && (nX < 0 || nZ < 0|| 
			       MULTIFILE::TimeCurrent < 0 || MULTIFILE::TimeCoupleNext < 0) ){
      if(fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_)==false){
	exit(__LINE__,__FILE__, "ERROR: the size of the grid couldn't be read.");
	break;
      }
      fin.CutInputStr(str1,str);
      
      if((strcmp("#TIME",str1)==0)){
	fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	MULTIFILE::TimeCurrent = strtod(str1, NULL);
      }
      else
	if((strcmp("#TIME_NEXT",str1)==0)){
	  fin.GetInputStr(str1,_MAX_STRING_LENGTH_PIC_);
	  MULTIFILE::TimeCoupleNext=strtod(str1, NULL);
	  if(MULTIFILE::TimeCoupleNext < 0)
	    return;
	    //exit(__LINE__,__FILE__,"Reached the last ARMS data file; exit");
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
    // keep reading the actual data -------------------------------------------
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
  // now the data has been read -----------------------------------------------

  //convert velocity and magnetic field vectors to cartesian coordinates ------
  for(int iZ = 0; iZ < nZ; iZ++){
    for(int iX = 0; iX < nX; iX++){
      double tmp;
      double radius   = pow(Xpos[iX]*Xpos[iX] + Zpos[iZ]*Zpos[iZ],0.5);
      double cosTheta = Zpos[iZ] / radius, sinTheta = Xpos[iX] / radius;
      // transform velocity and magnetic field vectors:
      // data is a slice of y=0 plane => need only Theta angle
      for(int i = b_; i==b_ || i==v_; i += (v_-b_) ){
	tmp              =Data[i][iX][iZ]*cosTheta-Data[i+1][iX][iZ]*sinTheta;
	Data[i  ][iX][iZ]=Data[i][iX][iZ]*sinTheta+Data[i+1][iX][iZ]*cosTheta;
	Data[i+1][iX][iZ]=Data[i+2][iX][iZ]*1;
	Data[i+2][iX][iZ]=tmp;
      }
    }
  }
  //conversion is finished ----------------------------------------------------

  //compute additional data: 
  //gradient of magnetic field, abs value of magnetic field and its gradient
  {    
    double dX = Xpos[1] - Xpos[0];
    double dZ = Zpos[1] - Zpos[0];

    //components of gradient in the y=0 slice, i.e. d(*)/dy = 0
    for(int iX = 0; iX < nX; iX++)
      for(int iZ = 0; iZ < nZ; iZ++){	
	if(iX<nX-1){
	  for(int i=0; i<3; i++){
	    Data[Db_+  3*i][iX][iZ]=(Data[b_+i][iX+1][iZ]-Data[b_+i][iX][iZ])/dX;
	  }
	}
	if(iZ<nZ-1){
	  for(int i=0; i<3; i++){
	    Data[Db_+2+3*i][iX][iZ]=(Data[b_+i][iX][iZ+1]-Data[b_+i][iX][iZ])/dZ;
	  }
	}
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
    
    if ((startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) && (startNode->block!=NULL)) {
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
	    double radpol2 = x[0]*x[0] + x[1]*x[1]+1E-15;
	    double radpol  = pow(radpol2,            0.5);
	    double radius  = pow(radpol2 + x[2]*x[2],0.5);
	    int xCell = floor( (radpol - Xpos[0]) / dX);
	    int zCell = floor( (x[2]   - Zpos[0]) / dZ);
	    if(xCell <= 0||zCell <= 0||xCell >= nX-2||zCell >= nZ-2) continue;
	    
	    // interpolation weights
	    double wX = (radpol-Xpos[xCell]) / dX;
	    double wZ = (x[2] - Zpos[zCell]) / dZ;

	    // threshold weight to switch interpolation stencils for gradients
	    double wXth = 0.5;
	    double wZth = 0.5;
	    
	    //interpolate values
	    double DataInterp[nvar] = {0};
	    //double GrdabsBInterp[3]= {0}, absBInterp=0, GrdBInterp[3][3]={0};
	    for(int ii=0;ii<2;ii++) 
	      for(int jj=0;jj<2;jj++){
		//interpolate variables stored at nodes
		for(int ivar=0; ivar<nvar && ivar<Db_; ivar++)
		  DataInterp[ivar] += 
		    Data[ivar][xCell+ii][zCell+jj] * 
		    ((1-wX)*(1-ii)+wX*ii) * ((1-wZ)*(1-jj)+wZ*jj);
		//interpolate gradients (stored at edges)
		// d(*)/dx
		for(int ivar=Db_; ivar<nvar; ivar+=3)
		  DataInterp[ivar] += (wX > wXth) 
		    ? 
		    Data[ivar][xCell+ii][zCell+jj] * 
		    ((1-wX+wXth)*(1-ii)+(wX-wXth)*ii) * ((1-wZ)*(1-jj)+wZ*jj)
		    :
		    Data[ivar][xCell-ii][zCell+jj] * 
		    ((1+wX-wXth)*(1-ii)+(wXth-wX)*ii) * ((1-wZ)*(1-jj)+wZ*jj);
		// d(*)/dz
		for(int ivar=Db_+2; ivar<nvar; ivar+=3)
		  DataInterp[ivar] += (wZ > wZth) 
		    ?
		    Data[ivar][xCell+ii][zCell+jj] * 
		    ((1-wX)*(1-ii)+wX*ii) * ((1-wZ+wZth)*(1-jj)+(wZ-wZth)*jj)
		    :
		    Data[ivar][xCell+ii][zCell-jj] * 
		    ((1-wX)*(1-ii)+wX*ii) * ((1+wZ-wZth)*(1-jj)+(wZth-wZ)*jj);
	      }
	    {
	      // rotate vetors to the current plane
	      // initial data is in slice y=0
	      double cosPhi    = x[0] / radpol, sinPhi   = x[1]   / radpol;
	      double tmp;
	      for(int ivec=0, i=vecs[ivec]; ivec<nvec; i=vecs[++ivec] ){
		tmp             = DataInterp[i]*cosPhi-DataInterp[i+1]*sinPhi;
		DataInterp[i+1] = DataInterp[i]*sinPhi+DataInterp[i+1]*cosPhi;
		DataInterp[i  ] = tmp;
	      }
	    }
  
	    //locate the cell
	    nd=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
	    if ((CenterNode=startNode->block->GetCenterNode(nd))==NULL) continue;
	    offset=CenterNode->GetAssociatedDataBufferPointer();
	    
	    //save the interpolated values
	    for (int idim=0;idim<3;idim++) {
	      *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticField.offset))        =DataInterp[b_+idim];
	      *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.offset))=DataInterp[Db_+idim];
	      *(3+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.offset))=DataInterp[Db_+3+idim];
	      *(6+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.offset))=DataInterp[Db_+6+idim];
	      *(idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.offset))   =DataInterp[v_+idim];
	    }
	    // E = -VxB
	    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.offset)) =
	      DataInterp[b_+1]*DataInterp[v_+2] - DataInterp[b_+2]*DataInterp[v_+1];
	    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.offset)) =
	      DataInterp[b_+2]*DataInterp[v_+0] - DataInterp[b_+0]*DataInterp[v_+2];
	    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::ElectricField.offset)) =
	      DataInterp[b_+0]*DataInterp[v_+1] - DataInterp[b_+1]*DataInterp[v_+0];
    
	    *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.offset))     =DataInterp[p_];
	    *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.offset))   =DataInterp[n_];
	    *((double*)(offset+PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.offset))     =DataInterp[t_];
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
