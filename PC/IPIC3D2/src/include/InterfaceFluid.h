/*******************************************************
InterfaceFluid_H handels all interaction between
a fluid model (BATSRUS/SWMF) and the particle model
in iPic3D

Writen by Lars Daldorff (daldorff@umich.edu) 15 Jan 2013

********************************************************/


#ifndef InterfaceFluid_H
#define InterfaceFluid_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "ConfigFile.h"
#include "input_array.h"
#include "Alloc.h"
#include "VCtopology3D.h"
#include "asserts.h"
#include "Grid3DCU.h"

#ifdef BATSRUS
#include "my_mpi.h"
#endif

using namespace std;

class InterfaceFluid
{ 
 private:
  
  //Error codes
  static const int readEOF=11;

  // BATSRUS output header file name 
  string FilenameFluid,FilenameRoot;
  // variables to read data from header file
  string TempString;  
  int INnProc;      // number of processors in the fluid code
  int INnStep;      // iteration number
  int INnPlotvar;   // number of variables in the fluid file
  int INnEqpar;     // number of equation parameters
  int INnRoot_D[3]; // number of root blocks
  int INdim;        // number of dimentions
  int INstep;       // from which step we read the file
  int INnIJK_D[3];  // number of cells per block
  double INtime,INxRange_I[6],INplotDx_D[3],INgridDx_D[3],INdt;

  double SItime; // time in SI units

  int nIJK_D[3];    // number of cells in the passed uniform grid
  double XYZ_D[3];  // physical size of the passed domain
  
  double ****Bc_GD;    // cell centered B 
  double ****State_GV; // node centered state variables

  int nVarFluid, nIonFluid, nSpecies, nIon, nVarCoupling;
  bool useMultiSpecies, useMultiFluid;

  int StartIdx_D[3], EndIdx_D[3];  // storage for starting/ending grid indexes of this processor
  
  static const int NG = 1; // number of ghost cell

  bool useAnisoP;  // Use anisotripic pressure

  bool useMhdPe; 
  
  double dt;            // scaled time step from fluid model
  double ParamDt;       // dt read in from input file

  double tUnitPic;      // conversenfactor to time used in IPIC3D 
  double invtUnitPic;   // 1/tUnitPic 
  double *Si2No_V; // array storing unit conversion factors
  double *No2Si_V; // array storing inverse unit conversion factors
  double Si2NoM, Si2NoV, Si2NoRho, Si2NoB, Si2NoP, Si2NoJ, Si2NoL;
  double No2SiV, No2SiL;
  double Lnorm, Unorm, Mnorm, Qnorm; // normalization units for length, velocity, mass and charge
                                     // Normalized q/m ==1 for proton in CGS units  
  long    nS; // number of particle species
  double *MoMi_S; // masses for the particles species
  double *QoQi_S; // charge for each particle species
  double PeRatio; // temperature ratio for electrons: PeRatio = Pe/Ptotal
  double SumMass; // Sum of masses of each particle species
  int nOverlap; // Number of grid point from boundary where we have MHD+PIC solution
  int nOverlapP; // Nomber of overlap cells for redistrebute particles
  int nIsotropic, nCharge; // Intrpolation region for curents/pressure (nIsotropoc) and charge 
                           // nCharge. <0 : do nothing; ==0 only ghost region; >0 interpolate inside domain

  unsigned long iSyncStep; // Iterator for sync with fluid
  long nSync; //

  int myrank; // this process mpi rank
  int thisrank; // my mpi rank ony used here! this objects

  VCtopology3D *_vct; // internal pointer to processor topolegy object

  bool isFirstTime;

  int iRegion;

  // Do not include ghost cells.
  int nxcLocal, nycLocal, nzcLocal;

  // Nodes, include ghost cells.
  int nxnLG, nynLG, nznLG;

  // Unless it is the first step to initilize PC, otherwise, only boundary
  // information is needed from GM. 
  bool doNeedBCOnly;

  int iCycle;

 protected:
  // Variables for IDL format output.
  static const int nDimMax=3;
  int nPlotFile;
  int *dnOutput_I;
  double *dtOutput_I, *plotDx_I;
  //The second dimension: xmin, xmax, ymin, ymax, zmin, zmax. 
  double **plotRange_ID; 
  string *plotString_I;
  string *plotVar_I;
  bool doSaveBinary;
  
 public:
  // These variables are also used in PSKOutput.h
  int *iRho_I, *iRhoUx_I, *iRhoUy_I, *iRhoUz_I, iBx,iBy,iBz,
    iPe,*iPpar_I,*iP_I,iJx,iJy,iJz, *iUx_I, *iUy_I, *iUz_I, iRhoTotal;

  int nBCLayer;
  bool useRandomPerCell;
  bool doUseOldRestart;
  string testFuncs;
  int iTest, jTest,kTest;
  
 private:
  
  /** Get nomal and pendicular vector to magnetic field */
  inline double **MagneticBaseVectors(const double Bx, const double By, const double Bz)const{
    double **norm_DD, inv;
    int Norm_,Perp1_, Perp2_, X_, Y_, Z_;
    Norm_  = 0;
    Perp1_ = 1;
    Perp2_ = 2;
    X_ = 0;
    Y_ = 1;
    Z_ = 2;
   
    norm_DD = newArr2(double,3,3);
  
    inv = 1.0/sqrt(Bx*Bx+By*By+Bz*Bz);
    norm_DD[Norm_][X_] = Bx*inv;
    norm_DD[Norm_][Y_] = By*inv;
    norm_DD[Norm_][Z_] = Bz*inv;
  
    if(norm_DD[Norm_][Z_] < 0.5) {
      norm_DD[Perp1_][X_] =  norm_DD[Norm_][Y_];
      norm_DD[Perp1_][Y_] = -norm_DD[Norm_][X_];
      norm_DD[Perp1_][Z_] =  0.0;
      norm_DD[Perp2_][X_] =  norm_DD[Norm_][Z_]*norm_DD[Norm_][X_];
      norm_DD[Perp2_][Y_] =  norm_DD[Norm_][Z_]*norm_DD[Norm_][Y_]; 
      norm_DD[Perp2_][Z_] = -pow(norm_DD[Norm_][X_],2) - pow(norm_DD[Norm_][Y_],2);
    }else{
      norm_DD[Perp1_][X_] =  0.0;
      norm_DD[Perp1_][Y_] =  norm_DD[Norm_][Z_];
      norm_DD[Perp1_][Z_] = -norm_DD[Norm_][Y_];
      norm_DD[Perp2_][X_] = -pow(norm_DD[Norm_][Y_],2) - pow(norm_DD[Norm_][Z_],2);
      norm_DD[Perp2_][Y_] =  norm_DD[Norm_][Y_]*norm_DD[Norm_][X_];
      norm_DD[Perp2_][Z_] =  norm_DD[Norm_][Z_]*norm_DD[Norm_][X_]; 
    }
  	
    inv = 1.0/sqrt(norm_DD[Perp1_][X_]*norm_DD[Perp1_][X_] +
  	           norm_DD[Perp1_][Y_]*norm_DD[Perp1_][Y_] +  
  	           norm_DD[Perp1_][Z_]*norm_DD[Perp1_][Z_]);  
    norm_DD[Perp1_][X_] *= inv;
    norm_DD[Perp1_][Y_] *= inv;
    norm_DD[Perp1_][Z_] *= inv;
  
    inv = 1.0/sqrt(norm_DD[Perp2_][X_]*norm_DD[Perp2_][X_] +
  	           norm_DD[Perp2_][Y_]*norm_DD[Perp2_][Y_] +  
  	           norm_DD[Perp2_][Z_]*norm_DD[Perp2_][Z_]);  
    norm_DD[Perp2_][X_] *= inv;
    norm_DD[Perp2_][Y_] *= inv;
    norm_DD[Perp2_][Z_] *= inv;
  
    return(norm_DD);
  }
  
  
  /** Convert from local processor index to full domain index */
 public:  inline void getGlobalIndex(const int il, const int jl, const int kl,
                             int *ig, int *jg, int *kg)const
  {
    if(StartIdx_D[0] == -1){
      *ig=il;
      *jg=jl;
      *kg=kl;
    } 
    else{
      // convert from local cpu index il to global domain index ig
      *ig =StartIdx_D[0] + il;
      
      // solution is constant across ignored dimensions indicated by nIJK_D == 1
      if(nIJK_D[1] == 1) *jg = 0; else *jg =StartIdx_D[1] + jl;
      if(nIJK_D[2] == 1) *kg = 0; else *kg =StartIdx_D[2] + kl;
    }
  }

 public:  inline void getLocalIndex(int *i, int *j, int *k)const
  {
    *i = *i - StartIdx_D[0];
      
    // solution is constant across ignored dimensions indicated by nIJK_D == 1
    if(nIJK_D[1] == 1) *j = 0; else *j = *j - StartIdx_D[1];
    if(nIJK_D[2] == 1) *k = 0; else *k = *k - StartIdx_D[2];
  }

  
  /** Convert from position to full domain index */
 public: inline void getGlobalIndex(const double x, const double y, const double  z,
			     int *ig, int *jg, int *kg)const
  {
    *ig = max(floor(x/INgridDx_D[0]),0.0);
    if(nIJK_D[1] == 1) *jg =0 ; else *jg = max(floor(y/INgridDx_D[1]),0.0);
    if(nIJK_D[2] == 1) *kg =0 ; else *kg = max(floor(z/INgridDx_D[2]),0.0);
  }
  
  
  /** Second order interpolation  for a given position */

  inline void getInterpolatedValue(const int i, const int j, const int k,
				   double *Var, int iVar)const
  {
    int i1,j1,k1;
    i1=i;j1=j;k1=k;
    if(nIJK_D[1] == 1) j1=0;
    if(nIJK_D[2] == 1) k1=0;
    *Var = State_GV[i1][j1][k1][iVar];
  }
  
  	
  /** Second order interpolation  for a given position */
  inline void getInterpolatedValue(const double x, const double y, const double  z, double *Var, int iVar)const
  {
    int i1,j1,k1,i2,j2,k2;
    double dx1,dx2,dy1,dy2,dz1,dz2;
  
    // Get the index of the for the nodes surounding the cell
    i1 = max(floor(x/INgridDx_D[0]),0.0)+1; i2=i1+1;
    if(nIJK_D[1] == 1){j1=0; j2=0;} else{j1 = max(floor(y/INgridDx_D[1]),0.0)+1; j2=j1+1;}
    if(nIJK_D[2] == 1){k1=0; k2=0;} else{k1 = max(floor(z/INgridDx_D[2]),0.0)+1; k2=k1+1;}
  
    // Get distenc form the cell walls
    dx1 = x/INgridDx_D[0] - i1+1; dx2 = 1.0-dx1;	
    if(nIJK_D[1] == 1){dy1=0.5; dy2=0.5;} else{dy1 = y/INgridDx_D[1] - j1+1; dy2 = 1.0-dy1;}
    if(nIJK_D[2] == 1){dz1=0.5; dz2=0.5;} else{dz1 = z/INgridDx_D[2] - k1+1; dz2 = 1.0-dz1;}

    getLocalIndex(&i1,&j1,&k1);
    getLocalIndex(&i2,&j2,&k2);
    
    *Var = dz2*(dy2*(dx2*State_GV[i1][j1][k1][iVar]
		     +           dx1*State_GV[i2][j1][k1][iVar])
		+      dy1*(dx2*State_GV[i1][j2][k1][iVar]
			    +           dx1*State_GV[i2][j2][k1][iVar]))
      + dz1*(dy2*(dx2*State_GV[i1][j1][k2][iVar]
		  +           dx1*State_GV[i2][j1][k2][iVar])
	     +      dy1*(dx2*State_GV[i1][j2][k2][iVar]
			 +           dx1*State_GV[i2][j2][k2][iVar]));    
  }
  
  
  void InitData()
  {
  
    // normalization variables
    double RHOnorm,Bnorm,Jnorm,Pnorm;
    
    Qnorm   = sqrt(Mnorm*Lnorm)*Unorm;
    RHOnorm = Mnorm/(Lnorm*Lnorm*Lnorm);
    Bnorm   = sqrt(RHOnorm)*Unorm;
    Pnorm   = RHOnorm*Unorm*Unorm;
  
    // Jnorm = Unorm*Qnorm/(Lnorm*Lnorm*Lnorm) is the correct formula
    // but the Unorm cancels out when we convert the 
    // curlB_SI/mu0 into curlB_CGS*c_CGS/4pi because c_CGS=Unorm.
    // Since Unorm is not necessarily defined below, we remove it here:
    Jnorm = Qnorm/(Lnorm*Lnorm*Lnorm);
  
    if(thisrank == 0) {
      cout<<"============= Normalization factors ============="<<endl;
      cout<<"Mnorm   = "<<Mnorm<<endl;
      cout<<"Unorm   = "<<Unorm<<endl;
      cout<<"Qnorm   = "<<Qnorm<<endl;
      cout<<"Lnorm   = "<<Lnorm<<endl;
      cout<<"RhoNorm = "<<RHOnorm<<endl;
      cout<<"Bnorm   = "<<Bnorm<<endl;
      cout<<"Pnorm   = "<<Pnorm<<endl;
      cout<<"Jnorm   = "<<Jnorm<<endl;
      //cout<<"========================================="<<endl;
    }


    // SI -> CGS conversion 
    Si2NoRho = 0.001;  // [kg/m^3] -> [g/cm^3] and get numnber desity
    Si2NoV   = 100.0;  // [m/s] -> [cm/s]
    Si2NoB   = 1.0e4;  // [Tesla] -> [gauss] 
    Si2NoP   = 10.0;   // [Pa] -> [Ba]
    Si2NoJ   = 1.0e-5; // [T/m]/mu0 -> [G/cm] c_CGS/4pi and
    Si2NoL   = 100;    // [m] ->[cm] 1.0/100.0*sqrt(1e9); 
    
    // Normalization: CGS -> non dimensional cgs    
    Si2NoRho /= RHOnorm;
    Si2NoV   /= Unorm;
    Si2NoB   /= Bnorm;
    Si2NoP   /= Pnorm;
    Si2NoJ   /= Jnorm;
    Si2NoL   /= Lnorm;

    Si2NoM = Si2NoRho*Si2NoV;

    No2SiV = 1./Si2NoV;
    No2SiL = 1./Si2NoL;
    
    Si2No_V[iBx] = Si2NoB;
    Si2No_V[iBy] = Si2NoB;
    Si2No_V[iBz] = Si2NoB;
    Si2No_V[iJx] = Si2NoJ;
    Si2No_V[iJy] = Si2NoJ;
    Si2No_V[iJz] = Si2NoJ;
    if(useMhdPe)Si2No_V[iPe] = Si2NoP;
    if(useMultiSpecies) Si2No_V[iRhoTotal] = Si2NoRho;
    for(int iIon=0; iIon<nIon; ++iIon){
      Si2No_V[iRho_I[iIon]]     = Si2NoRho;
      Si2No_V[iRhoUx_I[iIon]]   = Si2NoM;
      Si2No_V[iRhoUy_I[iIon]]   = Si2NoM;
      Si2No_V[iRhoUz_I[iIon]]   = Si2NoM;
      Si2No_V[iP_I[iIon]]       = Si2NoP;
      if(useAnisoP)Si2No_V[iPpar_I[iIon]] = Si2NoP;
    }

    // Get back to SI units
    for(int iVar=0;iVar<nVarCoupling;iVar++)
      No2Si_V[iVar] = 1.0/Si2No_V[iVar];
    }
 
  void ReNormLength()
  {
    // Normalization
    for(int i =0; i<3; i++){
      INgridDx_D[i]   *= Si2NoL;
      XYZ_D[i]        *= Si2NoL;
    }
  
    for(int i =0; i<6; i++){ INxRange_I[i]*=Si2NoL;}
  }	
  
  
  /** Convert to internal units in IPIC3D */
  void ReNormVariables()
  {
    // Convert all units to internal IPIC3D units
    for(int i=0;i<nxnLG;i++)
      for(int j=0;j<nynLG;j++)
	for(int k=0;k<nznLG;k++)	  
	  if(doGetFromGM(i,j,k)){
	    for(int iVar=0;iVar<=iJz;iVar++)
	      State_GV[i][j][k][iVar] *= Si2No_V[iVar];  
	  }
  }

  void Moment2Velocity(){
    int i,j,k;
    double Rhot; 

    if(INdim == 1) {
      k = 0;
      j  =0;
      for(i = 0; i < nxnLG; i++){
	if(doGetFromGM(i,0,0)){
	  if(useMultiSpecies){
	    Rhot=0;
	    for(int iIon=0; iIon<nIon; ++iIon){
	      // Rho = sum(Rhoi) + Rhoe;
	      Rhot +=
		State_GV[i][j][k][iRho_I[iIon]]*(1+MoMi_S[0]/MoMi_S[iIon+1]);
	    }// iIon

	    State_GV[i][j][k][iUx_I[0]] /= Rhot;
	    State_GV[i][j][k][iUy_I[0]] /= Rhot;
	    State_GV[i][j][k][iUz_I[0]] /= Rhot;	  
	  }else{
	    for(int iIon=0; iIon<nIon; ++iIon){
	      State_GV[i][j][k][iUx_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
	      State_GV[i][j][k][iUy_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
	      State_GV[i][j][k][iUz_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
	    }// iIon
	  } // else
	}
      } // i
    }
    else if (INdim == 2){
      k = 0;
      for(i = 0; i < nxnLG; i++){
	for(j = 0; j < nynLG; j++){
	  if(doGetFromGM(i,j,0)){
	    if(useMultiSpecies){
	      Rhot=0;
	      for(int iIon=0; iIon<nIon; ++iIon){
		// Rho = sum(Rhoi) + Rhoe;
		Rhot +=
		  State_GV[i][j][k][iRho_I[iIon]]*(1+MoMi_S[0]/MoMi_S[iIon+1]);
	      }// iIon
	      State_GV[i][j][k][iUx_I[0]] /= Rhot;
	      State_GV[i][j][k][iUy_I[0]] /= Rhot;
	      State_GV[i][j][k][iUz_I[0]] /= Rhot;	  
	    }else{
	      for(int iIon=0; iIon<nIon; ++iIon){
		State_GV[i][j][k][iUx_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
		State_GV[i][j][k][iUy_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
		State_GV[i][j][k][iUz_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
	      }// iIon
	    } // else
	  }
	}
      }
    }
    else{ 
      for(i = 0; i < nxnLG; i++){
	for(j = 0; j < nynLG; j++){
	  for(k = 0; k < nznLG; k++){
	    if(doGetFromGM(i,j,k)){
	      if(useMultiSpecies){
		Rhot=0;
		for(int iIon=0; iIon<nIon; ++iIon){
		  // Rho = sum(Rhoi) + Rhoe;
		  Rhot +=
		    State_GV[i][j][k][iRho_I[iIon]]*(1+MoMi_S[0]/MoMi_S[iIon+1]);
		}// iIon

		State_GV[i][j][k][iUx_I[0]] /= Rhot;
		State_GV[i][j][k][iUy_I[0]] /= Rhot;
		State_GV[i][j][k][iUz_I[0]] /= Rhot;	  
	      }else{
		for(int iIon=0; iIon<nIon; ++iIon){
		  State_GV[i][j][k][iUx_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
		  State_GV[i][j][k][iUy_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
		  State_GV[i][j][k][iUz_I[iIon]] /= State_GV[i][j][k][iRho_I[iIon]];
		}// iIon
	      } // else
	    }
	  }
	}
      }
    }    
  }
   
   
  /** for debugging, printing stat var matrix by variable */
  void PrintStateVar()
  {
    ofstream myfile;
    cout<<" PrintStateVar thisrank = "<<thisrank<<endl;
    if(thisrank == 0) myfile.open("state_vg_0.txt");
    if(thisrank == 1) myfile.open("state_vg_1.txt");
    if(thisrank == 2) myfile.open("state_vg_2.txt");
    if(thisrank == 3) myfile.open("state_vg_3.txt");
    if(thisrank == 4) myfile.open("state_vg_4.txt");
    if(thisrank == 5) myfile.open("state_vg_5.txt");
    if(thisrank == 6) myfile.open("state_vg_6.txt");
    if(thisrank == 7) myfile.open("state_vg_7.txt");
    if(thisrank == 8) myfile.open("state_vg_8.txt");
    if(thisrank == 9) myfile.open("state_vg_9.txt");
    if(thisrank == 10) myfile.open("state_vg_10.txt");
    if(thisrank == 11) myfile.open("state_vg_11.txt");
    if(thisrank == 12) myfile.open("state_vg_12.txt");
    if(thisrank == 13) myfile.open("state_vg_13.txt");
    if(thisrank == 14) myfile.open("state_vg_14.txt");
    if(thisrank == 15) myfile.open("state_vg_15.txt");

    myfile.setf(ios::fixed,ios::floatfield);
    myfile.precision(10);
    myfile.width(10);
    for(int k=0;k<nznLG;k++)
      {
	for(int j=0;j<nynLG;j++)
	  { 
	    for(int i=0;i<nxnLG; i++)
	      { 
		myfile<<setw(16)<<scientific<< i*INgridDx_D[0] + INxRange_I[0]<<" "<<j*INgridDx_D[1] + INxRange_I[2]<<"  ";
		for(int iVar=0;iVar<INnPlotvar +3;iVar++)
		  {
  		    myfile<<setw(16)<<scientific<<State_GV[i][j][k][iVar]<<"  ";
  		  }
		myfile<<"\n";
  	      }	
  	  }
  
      }
    myfile.close();
  }
  
  
  /** Help function for reading input variables */
  void InputToArray(ConfigFile *config,int ns, string Name, double *out)
  {
    array_double tmp = config->read<array_double>( Name);
    for(int i=0;i<ns;i++)
      switch(i){
      case 0:
        out[0] = tmp.a;
        break;
      case 1:
        out[1] = tmp.b;
        break;
      case 2:
        out[2] = tmp.c;
        break;
      case 3:
        out[3] = tmp.d;
        break;
      case 4:
        out[4] = tmp.e;
        break;
      case 5:
        out[5] = tmp.f;
        break;
      default:
  	{ 
  	  cout<<"Only 5 Species allowed!"<<flush;
  	  abort();
  	}
      }
  
  }
    
   
 public:

  /** constructor */
  InterfaceFluid(){
    
    SItime = 0.0; 
    MPI_Comm_rank(MPI_COMM_MYSIM, &thisrank);    
    for(int iDim=0;iDim<3;iDim++){
       StartIdx_D[iDim] = -1;
       EndIdx_D[iDim] = -1;
    }
    dt      = 0.0;
    ParamDt = 0.0; 
  
    // form normalizationof length C/Wp 100;  //[m] ->[cm]
    nSync = 1;
    iSyncStep = 1; // first sync time, 0 is at init

    // For now
    Lnorm = 1.0;
    Mnorm = 1.0;
    Unorm = 1.0;
   
    isFirstTime = true;
    doNeedBCOnly = false;
  }
  
  /** destructor */
  ~InterfaceFluid(){
    delete MoMi_S;
    delete QoQi_S;
    delArr4(State_GV,nxnLG,nynLG,nznLG);
    delArr4(Bc_GD,nxnLG,nynLG,nznLG);

    delete iRho_I;
    delete iRhoUx_I;
    delete iRhoUy_I;
    delete iRhoUz_I;
    delete iUx_I;
    delete iUy_I;
    delete iUz_I;
    delete iPpar_I;
    delete iP_I;

    if(nPlotFile>0){
      delete [] dnOutput_I;
      delete [] dtOutput_I;
      delete [] plotDx_I;
      delete [] plotString_I;
      delete [] plotVar_I;
      delArr2(plotRange_ID,2*nDimMax);
    }
  }

  /** return min X for fluid grid without ghostcell */
  inline double getFluidStartX(){return(INxRange_I[0]+INgridDx_D[0]);}
  
  
  /** return min Y for fluid grid without ghostcell */
  inline double getFluidStartY(){return(INxRange_I[2]+INgridDx_D[1]);}
  
  
  /** return min Z for fluid grid without ghostcell */
  inline double getFluidStartZ(){return(INxRange_I[4]+INgridDx_D[2]);}
  
  
  /** Average over nodes in ignored dimentions. Used by the original 2D file coupling. */
  inline double getValueNode(double ****vec, int i, int j, int k, int si)
  {
    double Var,N;
    
    N = 1;	
    Var = vec[si][i][j][k];
    if(nIJK_D[1] == 1){
      Var += vec[si][i][j+1][k];
      N++;
    }
    if(nIJK_D[2] == 1){
      Var += vec[si][i][j][k+1];
      N++;
      if(nIJK_D[1] == 1){
  	Var += vec[si][i][j+1][k+1];
  	N++;
      }
    }
    return(Var/N);
  }
  
  
  /** Average over nodes in ignored dimentions. Used by the original 2D file coupling. */
  inline double getValueNode(double ***vec, int i, int j, int k)
  {
    double Var,N;

    N = 1;	
    Var = vec[i][j][k];
    if(nIJK_D[1] == 1){
      Var += vec[i][j+1][k];
      N++;
    }
    if(nIJK_D[2] == 1){
      Var += vec[i][j][k+1];
      N++;
      if(nIJK_D[1] == 1){
  	Var += vec[i][j+1][k+1];
  	N++;
      }
    }
    return(Var/N);
  }
  
  
  /** return the cell value form the node */
  inline double getN2C(double ***vec, int i, int j, int k)
  {
    return(.125*(vec[i][j][k] + vec[i+1][j][k] + vec[i][j+1][k] + vec[i][j][k+1] 
  		 + vec[i+1][j+1][k] + vec[i+1][j][k+1] + vec[i][j+1][k+1] + vec[i+1][j+1][k+1]));
  }
  
  
  // return the cell value form the node
  inline double getN2C(double ****vec, int i, int j, int k,int is)
  {
    return(.125*(vec[is][i][j][k] + vec[is][i+1][j][k] + vec[is][i][j+1][k] + vec[is][i][j][k+1] 
  		 + vec[is][i+1][j+1][k] + vec[is][i+1][j][k+1] + vec[is][i][j+1][k+1] + vec[is][i+1][j+1][k+1]));
  }
  
  
  /** Find out if we will sync with fluid solution at this time/cycle */
  inline unsigned long doSyncWithFluid(unsigned long cycle)
  {
    if(cycle==iSyncStep*nSync){
      iSyncStep++;
      return((iSyncStep-1));
    }
    else 
      return(0);
  }

  /** Finding the start index of each processors subdomain */
  void setGlobalStartIndex(VCtopology3D *vct)
  {
    bool isCrowdedX, isCrowdedY, isCrowdedZ;
        
    // If vct == NULL it has bin set earlyer
    if( vct != NULL) _vct = vct;

    nxcLocal=1;
    nycLocal=1;
    nzcLocal=1;
    
    nxcLocal = getFluidNxc()/(double)_vct->getXLEN();
    if(INdim > 1) nycLocal = getFluidNyc()/(double)_vct->getYLEN();
    if(INdim > 2) nzcLocal = getFluidNzc()/(double)_vct->getZLEN();

    myrank = _vct-> getCartesian_rank();

    StartIdx_D[0] = _vct->getCoordinates(0)*nxcLocal;    
    StartIdx_D[1] = _vct->getCoordinates(1)*nycLocal;    
    StartIdx_D[2] = _vct->getCoordinates(2)*nzcLocal;    

    isCrowdedX =							\
      _vct->getCoordinates(0) < (getFluidNxc() - nxcLocal*_vct->getXLEN());
    isCrowdedY =							\
      _vct->getCoordinates(1) < (getFluidNyc() - nycLocal*_vct->getYLEN());
    isCrowdedZ =							\
      _vct->getCoordinates(2) < (getFluidNzc() - nzcLocal*_vct->getZLEN());
    
    if(isCrowdedX) {
      ++nxcLocal;
      StartIdx_D[0] += _vct->getCoordinates(0);
    }else{
	StartIdx_D[0] += getFluidNxc() - nxcLocal*_vct->getXLEN();	
    }    
    if(isCrowdedY) {
      ++nycLocal;
      StartIdx_D[1] += _vct->getCoordinates(1);
    }else{
      StartIdx_D[1] += (getFluidNyc() - nycLocal*_vct->getYLEN());
    }    
    if(isCrowdedZ) {
      ++nzcLocal;
      StartIdx_D[2] += _vct->getCoordinates(2);
    }else{
      StartIdx_D[2] += (getFluidNzc() - nzcLocal*_vct->getZLEN());
    }

    EndIdx_D[0] = StartIdx_D[0];
    EndIdx_D[1] = StartIdx_D[1];
    EndIdx_D[2] = StartIdx_D[2];

    EndIdx_D[0] += nxcLocal;    
    if(INdim > 1) EndIdx_D[1] += nycLocal;    
    if(INdim > 2) EndIdx_D[2] += nzcLocal;

    nxnLG = nxcLocal + 3;
    nynLG = nycLocal + 3;
    nznLG = nzcLocal + 3; 
  }	

  /** find which processor have with dordinate */
  void findProcForPoint(VCtopology3D *vct, int nPoint, double *Xyz_I, int *iProc_I){

    double invLx, invLy, invLz;
    double rank;
    int iPx, iPy, iPz;
    const double Epsilon = 0.0;//1.0e-13;
    myrank = vct-> getCartesian_rank();

    for(int i =0; i<nPoint*INdim; i++)
      Xyz_I[i]*=Si2NoL; 


    // Length of each processors sub domain for each dimentions.
    invLx = (double)vct->getXLEN()/(INgridDx_D[0]*getFluidNxn());
    invLy = (double)vct->getYLEN()/(INgridDx_D[1]*getFluidNyn());
    invLz = (double)vct->getZLEN()/(INgridDx_D[2]*getFluidNzn());

    if(INdim == 1){
      for(int i = 0; i < nPoint; i++){
	if(isThisRun(&Xyz_I[i*INdim])){
	  iProc_I[i] = (int)((Xyz_I[i*INdim    ]-INxRange_I[0]-Epsilon)*invLx);
	}
      }    
    }
    else if(INdim == 2){
      for(int i = 0; i < nPoint; i++){
	// Find out whether the point is in this simulation run region or not. 
	if(isThisRun(&Xyz_I[i*INdim])){
	  iPx = (int)((Xyz_I[i*INdim    ]-INxRange_I[0]-Epsilon)*invLx);
	  iPy = (int)((Xyz_I[i*INdim + 1]-INxRange_I[2]-Epsilon)*invLy);
	  
	  rank = iPy + vct->getYLEN()*iPx; 
	  iProc_I[i] = rank;
	}
      }  
    }
    else{
      for(int i = 0; i < nPoint; i++){
	if(isThisRun(&Xyz_I[i*INdim])){
	  iPx = (int)((Xyz_I[i*3    ]-INxRange_I[0]-Epsilon)*invLx);
	  iPy = (int)((Xyz_I[i*3 + 1]-INxRange_I[2]-Epsilon)*invLy);
	  iPz = (int)((Xyz_I[i*3 + 2]-INxRange_I[4]-Epsilon)*invLz);
	  
	  rank = iPz + vct->getZLEN()*(iPy + iPx*vct->getYLEN()); 
	  
	  iProc_I[i] = rank;
	}
      }
    }
    // if its used later ( may be removed)
    for(int i =0; i<nPoint*INdim; i++)
      Xyz_I[i]*=No2SiL; 

  }


  /* Is Pos_D inside this PIC domain? */
  bool isThisRun(double *Pos_D){
    bool isThisRunReturn;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double csmall; 
    csmall = 1e-7;
      
    xmin = getPhyXMin() - csmall;
    ymin = getPhyYMin() - csmall;
    zmin = getPhyZMin() - csmall;
    xmax = getPhyXMax() + csmall;
    ymax = getPhyYMax() + csmall;
    zmax = getPhyZMax() + csmall;

    isThisRunReturn = true;
    if(Pos_D[0] < xmin || Pos_D[0] > xmax) isThisRunReturn = false;
    if(INdim > 1){
      if(Pos_D[1] < ymin || Pos_D[1] > ymax) isThisRunReturn = false;
    }
    if(INdim > 2){
      if(Pos_D[2] < zmin || Pos_D[2] > zmax) isThisRunReturn = false;
    }
    
    return isThisRunReturn;
  }

  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNxc() const{ return(nIJK_D[0]-3*NG); }
  
  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNyc() const{ if(nIJK_D[1] > 3*NG) return(nIJK_D[1]-3*NG); else return(1); }
  
  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNzc() const{ if(nIJK_D[2] > 3*NG) return(nIJK_D[2]-3*NG); else return(1); }
  
  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNxn() { return(nIJK_D[0]-2*NG); }
  
  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNyn() { if(nIJK_D[1] > 3*NG) return(nIJK_D[1]-2*NG); else return(1); }
  
  /** nIJK_D includes 1 guard/ghost cell layer... */
  inline double getFluidNzn() { if(nIJK_D[2] > 3*NG) return(nIJK_D[2]-2*NG); else return(1); }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLx(){ return(XYZ_D[0]-NG*2*INgridDx_D[0]); }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLy()
  { 
    if(nIJK_D[1] > 3*NG) return(XYZ_D[1]-NG*2*INgridDx_D[1]); 
    else return(XYZ_D[1]); 
  }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLz()
  { 
    if(nIJK_D[2] > 3*NG) return(XYZ_D[2]-NG*2*INgridDx_D[2]); 
    else return(XYZ_D[2]); 
  }
  
  /** place holder for all variables like charge density that we know is 0.0 */
  template <typename Type>	
    inline double getFluidZero(const Type x, const Type y, const Type  z, const int is)const{
    return(0.0);
  }
  
  /** get Pxx from fluid */
  template <typename Type>	
    inline double getFluidPxx(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*(getFluidP(x,y,z,is)/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUx(x,y,z,is),2)));
  }
  
  /** get Pyy from fluid */
  template <typename Type>	
    inline double getFluidPyy(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*(getFluidP(x,y,z,is)/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUy(x,y,z,is),2)));
  }
  
  /** get Pzz from fluid */
  template <typename Type>	
    inline double getFluidPzz(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*(getFluidP(x,y,z,is)/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUz(x,y,z,is),2)));
  }
  
  /** get Pxy from fluid */
  template <typename Type>	
    inline double getFluidPxy(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidRhoNum(x,y,z,is)*getFluidUx(x,y,z,is)*getFluidUy(x,y,z,is));
  }
  
  /** get Pxz from fluid */
  template <typename Type>	
    inline double getFluidPxz(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidRhoNum(x,y,z,is)*getFluidUx(x,y,z,is)*getFluidUz(x,y,z,is));
  }
  
  /** get Pyz from fluid */
  template <typename Type>	
    inline double getFluidPyz(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidRhoNum(x,y,z,is)*getFluidUy(x,y,z,is)*getFluidUz(x,y,z,is));
  }
  
  /** get Jx from fluid */
  template <typename Type>	
    inline double getFluidJx(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidU(x,y,z,is,iUx_I,iJx)*getFluidRhoNum(x,y,z,is));
  }
  
  /** get Jy from fluid */
  template <typename Type>	
    inline double getFluidJy(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidU(x,y,z,is,iUy_I,iJy)*getFluidRhoNum(x,y,z,is));
  }
  
  /** get Jz from fluid */
  template <typename Type>	
    inline double getFluidJz(const Type x, const Type y,const Type z, const int is)const{
    return(QoQi_S[is]*getFluidU(x,y,z,is,iUz_I,iJz)*getFluidRhoNum(x,y,z,is));	
  }
  
  /** get bulk velocity Ux from fluid */
  template <typename Type>	
    inline double getFluidUx(const Type x, const Type y,const Type z, const int is)const{
    return(getFluidU(x,y,z,is,iUx_I,iJx)); 
  }
  
  /** get bulk velocity Uy from fluid */
  template <typename Type>	
    inline double getFluidUy(const Type x, const Type y,const Type z, const int is)const{
    return(getFluidU(x,y,z,is,iUy_I,iJy));
  }
  
  /** get bulk velocity Uz from fluid */
  template <typename Type>	
    inline double getFluidUz(const Type x, const Type y,const Type z, const int is)const{
    return(getFluidU(x,y,z,is,iUz_I,iJz)); 
  }
  
  /** get bulk velocity U in x, y or z from fluid */
  template <typename Type>	
    inline double getFluidU(const Type x, const Type y,const Type z, 
			    const int is, const int * iU_I, const int iJ)const
    {
      // Assume qe = -qi;
      double U,Rho,J, Rhoit, Qit, Rhot;
      
      if(useMultiFluid){
	if(is==0){
	  // Electron
	  /** Ue = (J - sum(ni*qi*Ui))/(ne*qe) 
	         = J/(ne*qe) + sum(ni*Ui)/ne */
	  double Ui, ni, ne; 
	  getInterpolatedValue(x,y,z,&J,iJ);
	  ne = getFluidRhoNum(x,y,z,0);
	  U  = J/(QoQi_S[0]*ne);

	  for(int iIon=0; iIon<nIon; ++iIon){
	    Ui = getFluidU(x,y,z,iIon+1,iU_I,iJ);
	    ni = getFluidRhoNum(x,y,z,iIon+1);
	    U += ni*Ui/ne;
	  }	  
	}else{
	  // Ion
	  getInterpolatedValue(x,y,z,&U, iU_I[is-1]);
	}
      }else{
	// Single fluid or multi-species.

	// Ui = U_{MHD} - me/qe*J/Rho;
	// where Rho is total density include electrons.
	// Ue = U_{MHD} - Rhoit/Qit*J/Rho;
	// where Rhoit = sum(ni*Mi), Qit = sum(ni*Qi).	

	double moq,Numi;
	if(is==0){
	  Rhoit = 0; Qit = 0;
	  for(int iIon=0; iIon<nIon; ++iIon){
	    Numi   = getFluidRhoNum(x,y,z,iIon+1);
	    Rhoit += Numi*MoMi_S[iIon+1];
	    Qit   += Numi*QoQi_S[iIon+1];
	    }
	  moq = Rhoit/Qit;
	}else moq = MoMi_S[0]/QoQi_S[0];

	Rhot = 0;
	for(int is0=0; is0<nS; ++is0)
	  Rhot += MoMi_S[is0]*getFluidRhoNum(x,y,z,is0);
	
	getInterpolatedValue(x,y,z,&U,iU_I[0]);
	getInterpolatedValue(x,y,z,&J,iJ);	
	U -= moq*J/Rhot;
      }
      return(U);
    }
  
  /** get thermal velocity at location x, y, z for fluid is */
  template <typename Type>
    inline double getFluidUth(const Type x, const Type y ,const Type z, const int is)const
    {
      double Uth, p, ni;

      p  = getFluidP(x,y,z,is);
      ni = getFluidRhoNum(x,y,z,is);
      Uth = sqrt(p/(ni*MoMi_S[is]));
      return(Uth);
    }
  
  /** set thermal velocity in magnetic cordinates for a given position */
  template <typename Type>
    inline void setFluidansioUth(const Type x, const Type y ,const Type z,
				 double *u, double *v, double *w, 
				 const double rand1, const double rand2, 
				 const double rand3, const double rand4, 
				 const int is)const
    {
      double Rho, Bx,By,Bz,B,P,Ppar,Pperp,Uthperp,Uthpar,Uthperp1, Uthperp2;
      double harvest, prob, theta;
      double **norm_DD;
      // indexes for the norm_DD matix
      int Norm_, Perp1_, Perp2_, X_, Y_, Z_;

      if(useMultiFluid || useMultiSpecies){
	cout<<" setFluidanisoUth has not implemented for multifluid/multispecies!!!"<<endl;
	abort();	
      }


      Norm_  = 0;
      Perp1_ = 1;
      Perp2_ = 2;
      X_ = 0;
      Y_ = 1;
      Z_ = 2;

      // Get rho and B at the particle position
      getInterpolatedValue(x,y,z,&Rho,iRho_I[0]);	
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);	
  
      // Get Parallel and perpendicular presure
      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);
     
      // Get 3 vertors spaning the vector space
      norm_DD = MagneticBaseVectors(Bx,By,Bz);
  
      //Get the termail verlocity's
      prob  = sqrt(-2.0*log(1.0-.999999999*rand1));
      theta = 2.0*M_PI*rand2;
      Uthpar = sqrt(Ppar*SumMass/(MoMi_S[is]*Rho))*prob*cos(theta);
  
      prob  = sqrt(-2.0*log(1.0-.999999999*rand3));
      theta = 2.0*M_PI*rand4;
      Uthperp  = sqrt(Pperp*SumMass/(MoMi_S[is]*Rho))*prob;
      Uthperp1 = Uthperp*cos(theta);
      Uthperp2 = Uthperp*sin(theta);
  
      // Set particle thermal velocity
      (*u) = Uthpar*norm_DD[Norm_][X_] + Uthperp1*norm_DD[Perp1_][X_] + Uthperp2*norm_DD[Perp2_][X_];
      (*v) = Uthpar*norm_DD[Norm_][Y_] + Uthperp1*norm_DD[Perp1_][Y_] + Uthperp2*norm_DD[Perp2_][Y_];
      (*w) = Uthpar*norm_DD[Norm_][Z_] + Uthperp1*norm_DD[Perp1_][Z_] + Uthperp2*norm_DD[Perp2_][Z_];
  
      //Cleanig
      delArr2(norm_DD,3);
    }
  
  /** get total pressure for species from fluid */
  template <typename Type>	
    inline double getFluidP(const Type x, const Type y,const Type z, const int is)const{
    double P, Pe, Pi;

    
    if(useMultiFluid){
      // Multi-fluid.
      if(is==0)
	getInterpolatedValue(x,y,z,&P,iPe); // Electron
      else	
	getInterpolatedValue(x,y,z,&P,iP_I[is-1]); // Ion
    }else{
      // Single-fluid and multi-species.
      if(! useMhdPe){
	getInterpolatedValue(x,y,z,&P,iP_I[0]);
	if(is == 0) P *= PeRatio;
	else if(is > 0) P *= (1 - PeRatio);	
      }else{
	if(is ==0 ) getInterpolatedValue(x,y,z,&P,iPe); // Electron
	else if(is >0) getInterpolatedValue(x,y,z,&P,iP_I[0]); // Ion
      }

      // Split pressure among ions.
      if(useMultiSpecies && is>0){
	double Numit;
	Numit = 0; // Number of all ions.
	for(int iIon=0; iIon<nIon; ++iIon)
	  Numit += getFluidRhoNum(x,y,z,iIon+1);

	P *= getFluidRhoNum(x,y,z,is)/Numit;
      }      
    }

    return(P);
  }
  
  /** get parallel thermal presure for species from a fluid*/ 
  template <typename Type>	
    inline double getFluidPpar(const Type x, const Type y,const Type z, const int is)const{
    // Need to check whether this function works correctly! -- Yuxi
    double P;
    if(useMultiSpecies || useMultiFluid){
      cout<<" getFluidPpar has not implemented for multifluid/multispecies!!"<<endl;
      abort();      
    }

    getInterpolatedValue(x,y,z,&P,iPpar_I[0]);
    if(is == 0) P *= PeRatio;
    else if(is > 0) P *= (1 - PeRatio);		
    return(P);
  }
  
  /** Get max thermal velocity in the domain. Result is on root processor only */
  double getMaxFluidUth(const int is, const int dir)const
  {
    double maxVth, allmaxVth;
    // Helping variable for looping
    // one step over ignored dimention
    int fix_D[3] = {0};
    if(INdim < 2) fix_D[1] = 1;
    if(INdim < 3) fix_D[2] = 1;

    maxVth=0.0;
    for(int i = 0; i<EndIdx_D[0] - StartIdx_D[0] + fix_D[0]; i++)
      for(int j = 0; j< EndIdx_D[1] - StartIdx_D[1] + fix_D[1]; j++)
  	for(int k = 0; k< EndIdx_D[2] - StartIdx_D[2] + fix_D[2]; k++)
	  maxVth = max(maxVth, getFluidUth(i,j,k,is));
    allmaxVth = maxVth; 
    MPI_Reduce(&maxVth, &allmaxVth, 1, MPI_DOUBLE, MPI_MAX, 0 , MPI_COMM_MYSIM);
    return(allmaxVth);
  }


  /** this should return NUMBER DENSITY from fluid */
  template <typename Type>
    inline double getFluidRhoNum(const Type x,const Type  y, const Type z, const int is)const{
    double Rho, NumDens;
    
    if(useMultiFluid || useMultiSpecies){
      if(is == 0){
	// Electron
	NumDens = 0; 
	for(int iIon=0; iIon<nIon; ++iIon){
	  getInterpolatedValue(x,y,z,&Rho,iRho_I[iIon]);
	  NumDens += Rho/MoMi_S[iIon+1];
	}
      }else{
	// Ion
	getInterpolatedValue(x,y,z,&Rho,iRho_I[is-1]);
	NumDens = Rho/MoMi_S[is];
      }      
    }else{
      // Electrons and iones have same density, ignoring is
      getInterpolatedValue(x,y,z,&Rho,iRho_I[0]);
      NumDens = Rho/SumMass;
    }

    return(NumDens);
  }
  
  /** Get the Electic field as from the fluid description */
  inline void setFluidFieldsNode(
  				 double *Ex, double *Ey, double *Ez, 
  				 double *Bx, double *By, double *Bz, 
  				 const int ii,const int jj, const int kk)
  {
    if(doGetFromGM(ii,jj,kk)){
      int i,j,k;
      i=ii;j=jj;k=kk;
      if(nIJK_D[1] == 1) j=0;
      if(nIJK_D[2] == 1) k=0;

      //cout<<"mhd: E = - Ue x B"<<endl;
      // Ue = Ui - J/ne
      (*Ex) = (getFluidUz(ii,jj,kk,0)*State_GV[i][j][k][iBy] - getFluidUy(ii,jj,kk,0)*State_GV[i][j][k][iBz]);
      (*Ey) = (getFluidUx(ii,jj,kk,0)*State_GV[i][j][k][iBz] - getFluidUz(ii,jj,kk,0)*State_GV[i][j][k][iBx]);
      (*Ez) = (getFluidUy(ii,jj,kk,0)*State_GV[i][j][k][iBx] - getFluidUx(ii,jj,kk,0)*State_GV[i][j][k][iBy]);
      (*Bx) = State_GV[i][j][k][iBx];
      (*By) = State_GV[i][j][k][iBy];
      (*Bz) = State_GV[i][j][k][iBz];
    }
  }
  
  
  // Data recived from SWMF coupler
  void ReadFromGMinit(int *paramint, double *griddim, double *paramreal, stringstream *ss){

    INdim       = paramint[0];
    nVarFluid   = paramint[2]; 
    nIonFluid   = paramint[3];
    nSpecies    = paramint[4];      
    nS          = nIonFluid + nSpecies; // nFluid+nSpecies-1+1 (electrons)
    nIon        = nS - 1;
    nVarCoupling= nVarFluid + 3; // nVarFluid + (Jx, Jy, Jz)
    useMultiFluid   = nIonFluid > 1;
    useMultiSpecies = nSpecies > 1;
    
    iRho_I   = new int[nIon];
    iRhoUx_I = new int[nIon];
    iRhoUy_I = new int[nIon];
    iRhoUz_I = new int[nIon];
    iUx_I    = new int[nIon];
    iUy_I    = new int[nIon];
    iUz_I    = new int[nIon];
    iPpar_I  = new int[nIon];
    iP_I     = new int[nIon];      

    // c++ index starts from 0. So, minus 1. 
    iPe      = paramint[5] - 1;
    iBx      = paramint[6] - 1;
    iBy      = iBx + 1;
    iBz      = iBy + 1;

    int n = 7;
    if(useMultiSpecies){
      // MultiSpecies. Densities of each species are known. Total velocity 
      // and total pressure are known. 
      iRhoTotal   = paramint[n++] - 1;
      iRho_I[0]   = iRhoTotal + 1; 
      iRhoUx_I[0] = paramint[n++] - 1; iUx_I[0] = iRhoUx_I[0];
      iRhoUy_I[0] = iRhoUx_I[0] + 1;   iUy_I[0]  = iRhoUy_I[0];
      iRhoUz_I[0] = iRhoUx_I[0] + 2;   iUz_I[0]  = iRhoUz_I[0];      
      iPpar_I[0]  = paramint[n++] - 1;
      iP_I[0]     = paramint[n++] - 1;

      for(int iIon=1; iIon<nIon; ++iIon){
	iRho_I[iIon]      = iRho_I[0] + iIon;
	iRhoUx_I[iIon]    = iRhoUx_I[0]; iUx_I[iIon] = iUx_I[0];
	iRhoUy_I[iIon]    = iRhoUy_I[0]; iUy_I[iIon] = iUy_I[0];
	iRhoUz_I[iIon]    = iRhoUz_I[0]; iUz_I[iIon] = iUz_I[0];
	iPpar_I[iIon]     = iPpar_I[0];
	iP_I[iIon]        = iP_I[0];
      }      
    }else{
      // Single fluid or multi-fluid. 
      for (int iIon = 0; iIon < nIon; ++iIon)
	iRho_I[iIon]     = paramint[n++] - 1;
      for (int iIon = 0; iIon < nIon; ++iIon){
	iRhoUx_I[iIon]   = paramint[n++] - 1;  iUx_I[iIon]  = iRhoUx_I[iIon];
	iRhoUy_I[iIon]   = iRhoUx_I[iIon] + 1; iUy_I[iIon]  = iRhoUy_I[iIon];
	iRhoUz_I[iIon]   = iRhoUx_I[iIon] + 2; iUz_I[iIon]  = iRhoUz_I[iIon];
      }
      for (int iIon = 0; iIon < nIon; ++iIon)
	iPpar_I[iIon]    = paramint[n++] - 1;
      for (int iIon = 0; iIon < nIon; ++iIon)
	iP_I[iIon]       = paramint[n++] - 1;
    }

    // See GM/BATSRUS/src/ModExtraVariables.f90. 
    useAnisoP = iPpar_I[0] != 0;
    useMhdPe  = iPe != 0;

    iJx = nVarFluid;
    iJy = iJx + 1;
    iJz = iJx + 2;

    /* cout<<" iPe = "<<iPe<<" iBx = "<<iBx<<" iBy = "<<iBy<<" iBz = "<<iBz<<endl; */
    /* for (int iIon = 0; iIon<nIon; ++iIon){ */
    /*   cout<<" iIon = "<<iIon<<endl; */
    /*   cout<<" iRho = "<<iRho_I[iIon]<<" iRhoUx = "<<iRhoUx_I[iIon] */
    /* 	  <<" iRhoUy = "<<iRhoUy_I[iIon]<<" iRhoUz = "<<iRhoUz_I[iIon] */
    /* 	  <<" iPpar = "<<iPpar_I[iIon]<<" iP = "<<iP_I[iIon]<<endl; */
    /* } */
    /* cout<<" iJx = "<<iJx<<" iJy = "<<iJy<<" iJz = "<<iJz<<endl; */

    Si2No_V = new double[nVarCoupling];
    No2Si_V = new double[nVarCoupling];
    for(int i =0; i<nVarCoupling; i++) Si2No_V[i] = 1;
    
    n = 0; 
    for(int i =0; i<3; i++){
      INxRange_I[2*i]     = griddim[n++];      // Lmin
      INxRange_I[2*i + 1] = griddim[n++]; // Lmax
      INgridDx_D[i]       = griddim[n++]; //dx
      XYZ_D[i] = (INxRange_I[2*i+1] - INxRange_I[2*i]);      
    }

    for(int i =0; i<3; i++)
      nIJK_D[i] = (int)(XYZ_D[i]/INgridDx_D[i] +0.5);

    // Add 1 as we count the nodes not cells
    for(int i =0; i<INdim; i++)
      nIJK_D[i] += 1;
    
    QoQi_S = new double[nS];
    MoMi_S = new double[nS];

    // electron charge is always -1
    QoQi_S[0] = -1.0;

    /** Do not change the order of the following lines. */
    n = 0; 
    for(int i = 1; i < nS; ++i){
      QoQi_S[i] = paramreal[n++];
      MoMi_S[i] = paramreal[n++];
    }
    // Electron pressure ratio: Pe/Ptotal
    PeRatio = paramreal[n++];
    Lnorm   = paramreal[n++]; 	
    Unorm   = paramreal[n++];	
    Mnorm   = paramreal[n++];
    /** Do not change the order of above lines. */
 
    // Normalization units converted [SI] -> [cgs]
    if(Lnorm > 0){ Lnorm *=100.0;}  else {Lnorm =1.0;}
    if(Unorm > 0){ Unorm *=100.0;}  else {Unorm =1.0;}
    if(Mnorm > 0){ Mnorm *=1000.0;} else {Mnorm =1.0;}
    
    InitData();
    ReNormLength();
    checkParam();
  }

  /** Check the parameters passed or calculated from BATSRUS*/
  void checkParam(){
    assert_ge(XYZ_D[0],0.0);
    assert_ge(XYZ_D[1],0.0);
    assert_ge(XYZ_D[2],0.0);

    assert_ge(INgridDx_D[0],0.0);
    assert_ge(INgridDx_D[1],0.0);
    assert_ge(INgridDx_D[2],0.0);

    assert_eq(QoQi_S[0],-1.0*QoQi_S[1]);
    for(int iIon = 1; iIon<nIon; ++iIon){
      if(QoQi_S[iIon+1] != QoQi_S[1]){
	cout<<" So far, -qe = qi1 = qi2 = qi3... is assumed. QoQi_S[1]= "
	    <<QoQi_S[1]<<" and QoQi_S["<<iIon+1<<"]= "<<QoQi_S[iIon+1]
	    <<" are not valided!!"<<endl;
	abort();
      }
    }

    if(useMultiFluid && !useMhdPe) {
      cout<<" Use multi-fluid but do not use electron pressure. This case is not supported so far!!!"<<endl;
      abort();
    }    
  }

  bool doGetFromGM(int i,int j,int k){
    // 1) The first step initialize from GM: all nodes need information from GM;
    // 2) During updates, only boundary nodes needs to be overwritten by GM.
    // Input: i, j, k are the local indexes. 
    bool doGetGM;
    int minDn;
    if(doNeedBCOnly){
      doGetGM = false;
      int ig = getGlobalStartIndex(0) + i; int nxg = getFluidNxc()+3;
      int jg = getGlobalStartIndex(1) + j; int nyg = getFluidNyc()+3;
      int kg = getGlobalStartIndex(2) + k; int nzg = getFluidNzc()+3;
      minDn = min(ig, nxg -ig -1);
      if(INdim>1) minDn = min(minDn, min(jg, nyg - jg - 1));
      if(INdim>2) minDn = min(minDn, min(kg, nzg - kg - 1));
      if(minDn <nBCLayer) doGetGM = true;
    }else{
      doGetGM = true;
    }
    return doGetGM;
  }
  
  /** get the number of cordinate points used by the simulation. recived from GM */
  void GetNgridPnt(int *nPoint){
    int i,j,k;
    *nPoint = 0;
    if(INdim == 1) {
      for(i = 0; i <= EndIdx_D[0] - StartIdx_D[0] + 2*NG; i++)
	if(doGetFromGM(i,0,0)) *nPoint +=1;
    }
    else if (INdim == 2){
      for(j = 0; j <=  EndIdx_D[1] - StartIdx_D[1] + 2*NG; j++)
	for(i = 0; i <= EndIdx_D[0] - StartIdx_D[0] + 2*NG; i++){
	  if(doGetFromGM(i,j,0)) *nPoint +=1;
	}
    }
    else{
      for(k = 0; k <=  EndIdx_D[2] - StartIdx_D[2] + 2*NG; k++)
	for(j = 0; j <=  EndIdx_D[1] - StartIdx_D[1] + 2*NG; j++)
	  for(i = 0; i <= EndIdx_D[0] - StartIdx_D[0] + 2*NG; i++){
	    if(doGetFromGM(i,j,k)) *nPoint +=1;
	  }
    }
    *nPoint *= INdim;
  }

  /** Get index of node points recived from GM */ 
  void GetGridPnt(double *Pos_I){
    int n, i, j, k;

    n=0;
    if(INdim == 1) {
      for(i = StartIdx_D[0]; i <= EndIdx_D[0] + 2*NG; i++)
	if(doGetFromGM(i-StartIdx_D[0],0,0))
	  Pos_I[i] = i*INgridDx_D[0]  + INxRange_I[0];
    }
    else if (INdim == 2){
      for(j = StartIdx_D[1]; j <=  EndIdx_D[1] + 2*NG; j++)
	for(i = StartIdx_D[0]; i <= EndIdx_D[0] + 2*NG; i++){
	  if(doGetFromGM(i-StartIdx_D[0],j-StartIdx_D[1],0)){
	    Pos_I[n++] = i*INgridDx_D[0] + INxRange_I[0];
	    Pos_I[n++] = j*INgridDx_D[1] + INxRange_I[2];
	  }
	}
    }
    else{ 
      for(k = StartIdx_D[2]; k <=  EndIdx_D[2] + 2*NG; k++)
        for(j = StartIdx_D[1]; j <=  EndIdx_D[1] + 2*NG; j++)
	  for(i = StartIdx_D[0]; i <= EndIdx_D[0] + 2*NG; i++){
	    if(doGetFromGM(i-StartIdx_D[0],j-StartIdx_D[1],k-StartIdx_D[2])){
	      Pos_I[n++] = i*INgridDx_D[0] + INxRange_I[0];
	      Pos_I[n++] = j*INgridDx_D[1] + INxRange_I[2];
	      Pos_I[n++] = k*INgridDx_D[2] + INxRange_I[4];
	    }
	  }
    }
    //to SI
    for(i = 0; i<n; i++){
      Pos_I[i] *= No2SiL;
    }
  }

  void setStateVar(double *State_I, int *iPoint_I){
    const int x_=0, y_=1, z_=2;
    int i, j, k, iVar;
    int ii, jj, kk;
    int idx, Nx, Ny, Nz; 

    if(isFirstTime){
      InitData();
      // Node values. 
      State_GV = newArr4(double,nxnLG,nynLG,nznLG,nVarCoupling);
      // Cell center values.
      Bc_GD = newArr4(double,nxnLG-1,nynLG-1,nznLG-1,3);
    }

    Nx = EndIdx_D[0]-StartIdx_D[0]+1 + 2*NG;
    Ny = EndIdx_D[1]-StartIdx_D[1]+1 + 2*NG;
    Nz = EndIdx_D[2]-StartIdx_D[2]+1 + 2*NG;

    // The order to loop through State_GV should be the same as the
    // order to find out the posion, which is in GetGridPnt().
    if(INdim == 1) {
      k = 0;
      j  =0;
      ii = 0;
      for(i = 0; i < nxnLG; i++){
	if(doGetFromGM(i,0,0)){
	  for(iVar=0; iVar < nVarFluid; iVar++){
 	    // convert 3D index to 1D
 	    idx = iVar + nVarFluid*(iPoint_I[ii] - 1);
	    State_GV[i][j][k][iVar] = State_I[idx];
	  }
          ii++;
	}
      }
    }
    else if (INdim == 2){
      k = 0;
      ii = 0;
      for(j = 0; j < nynLG; j++){
	for(i = 0; i < nxnLG; i++){
	  if(doGetFromGM(i,j,0)){
	    for(iVar=0; iVar < nVarFluid; iVar++){
	      idx = iVar + nVarFluid*(iPoint_I[ii] - 1);
	      State_GV[i][j][k][iVar] = State_I[idx];
	    }
	    ii++;
	  }
	}
      }
    }
    else{
      ii = 0;
      for(k = 0; k < nznLG; k++){
	for(j = 0; j < nynLG; j++){           
	  for(i = 0; i < nxnLG; i++){
	    if(doGetFromGM(i,j,k)){
	      for(iVar=0; iVar < nVarFluid; iVar++){
		idx = iVar + nVarFluid*(iPoint_I[ii] - 1);
		State_GV[i][j][k][iVar] = State_I[idx];
	      }
	      ii++;
	    }
	  }
	}
      }
    }


    // Get magnetic field B in the cell center form node values   

    if(INdim == 1){
      for(i = 0; i < nxnLG-1; i++){
	if(doGetFromGM(i,0,0)){
	Bc_GD[i][0][0][x_] = 0.5 *(State_GV[i][0][0][iBx] + State_GV[i+1][0][0][iBx]);
           
	Bc_GD[i][0][0][y_] = 0.5 *(State_GV[i][0][0][iBy] + State_GV[i+1][0][0][iBy]);
           
	Bc_GD[i][0][0][z_] = 0.5 *(State_GV[i][0][0][iBz] + State_GV[i+1][0][0][iBz]);
      } 
    }
    }
    else if(INdim == 2){
      for(i = 0; i < nxnLG-1; i++)
	for(j = 0; j < nynLG-1; j++){
	  if(doGetFromGM(i,j,0)){
	  Bc_GD[i][j][0][x_] = 0.25 *(
				      State_GV[i][j][0][iBx] + State_GV[i+1][j][0][iBx] +
				      State_GV[i][j+1][0][iBx] + State_GV[i+1][j+1][0][iBx]);
           
	  Bc_GD[i][j][0][y_] = 0.25 *(
				      State_GV[i][j][0][iBy] + State_GV[i+1][j][0][iBy] +
				      State_GV[i][j+1][0][iBy] + State_GV[i+1][j+1][0][iBy]);
           
	  Bc_GD[i][j][0][z_] = 0.25 *( 
				      State_GV[i][j][0][iBz] + State_GV[i+1][j][0][iBz] +
				      State_GV[i][j+1][0][iBz] + State_GV[i+1][j+1][0][iBz]);
	}          
    }
    }
    else if(INdim == 3){
      for(i = 0; i < nxnLG-1; i++)
        for(j = 0; j < nynLG-1; j++)
	  for(k = 0; k < nznLG-1; k++) {
	    if(doGetFromGM(i,j,k)){
	    Bc_GD[i][j][k][x_] = 0.125 *(
					 State_GV[i][j][k][iBx] + State_GV[i+1][j][k][iBx] +
					 State_GV[i][j+1][k][iBx] + State_GV[i+1][j+1][k][iBx] +
					 State_GV[i][j][k+1][iBx] + State_GV[i+1][j][k+1][iBx] +
					 State_GV[i][j+1][k+1][iBx] + State_GV[i+1][j+1][k+1][iBx]);
           
	    Bc_GD[i][j][k][y_] = 0.125 *(
					 State_GV[i][j][k][iBy] + State_GV[i+1][j][k][iBy] +
					 State_GV[i][j+1][k][iBy] + State_GV[i+1][j+1][k][iBy] +
					 State_GV[i][j][k+1][iBy] + State_GV[i+1][j][k+1][iBy] +
					 State_GV[i][j+1][k+1][iBy] + State_GV[i+1][j+1][k+1][iBy]);
           
	    Bc_GD[i][j][k][z_] = 0.125 *(
					 State_GV[i][j][k][iBz] + State_GV[i+1][j][k][iBz] +
					 State_GV[i][j+1][k][iBz] + State_GV[i+1][j+1][k][iBz] +
					 State_GV[i][j][k+1][iBz] + State_GV[i+1][j][k+1][iBz] +
					 State_GV[i][j+1][k+1][iBz] + State_GV[i+1][j+1][k+1][iBz]);
	  }          
    }
    }

    // J = curl B/ mu0
    double invdx, invdy, invdz;
    double compZdx, compXdy, compYdz, compXdz, compZdy, compYdx;


    // 1/(dx*mu0)  
    invdx = 1.0/(INgridDx_D[x_]*No2SiL*4.0*3.14159265359*1.0e-7);
    invdy = 1.0/(INgridDx_D[y_]*No2SiL*4.0*3.14159265359*1.0e-7);
    invdz = 1.0/(INgridDx_D[z_]*No2SiL*4.0*3.14159265359*1.0e-7);

    if(INdim == 1){
      j = 0;
      k = 0;
      for(i = 0; i < nxnLG-2; i++){
	if(doGetFromGM(i,0,0)){
	compZdx = invdx*( 
			 Bc_GD[i+1][j  ][k  ][z_] - Bc_GD[i  ][j  ][k  ][z_]);

	compYdx = invdx*( 
			 Bc_GD[i+1][j  ][k  ][y_] - Bc_GD[i  ][j  ][k  ][y_]);

	State_GV[i+1][j][k][iJx] =  0.0;
	State_GV[i+1][j][k][iJy] =  - compZdx;
	State_GV[i+1][j][k][iJz] =  compYdx;
      }
      }
      // fill in the ghost cell, constant value
      // Only the "real" ghost cells, which are at the boundaries are useful.
      State_GV[0][j][k][iJx] = State_GV[1][j][k][iJx];
      State_GV[0][j][k][iJy] = State_GV[1][j][k][iJy];
      State_GV[0][j][k][iJz] = State_GV[1][j][k][iJz];

      State_GV[nxnLG-1][j][k][iJx] = State_GV[nynLG-2][j][k][iJx];
      State_GV[nxnLG-1][j][k][iJy] = State_GV[nxnLG-2][j][k][iJy];
      State_GV[nxnLG-1][j][k][iJz] = State_GV[nxnLG-2][j][k][iJz];
    }
    else if(INdim == 2){
      k = 0;      
      for(i = 0; i < nxnLG-2; i++)
	for(j = 0; j < nynLG-2; j++){
	  if(doGetFromGM(i,j,0)){
	  compZdy = 0.5*invdy*( 
			       Bc_GD[i  ][j+1][k  ][z_] - Bc_GD[i  ][j][k  ][z_] +
			       Bc_GD[i+1][j+1][k  ][z_] - Bc_GD[i+1][j][k  ][z_]);

	  compXdy = 0.5*invdy*( 
			       Bc_GD[i  ][j+1][k  ][x_] - Bc_GD[i  ][j][k  ][x_] +
			       Bc_GD[i+1][j+1][k  ][x_] - Bc_GD[i+1][j][k  ][x_]);

	  compZdx = 0.5*invdx*( 
			       Bc_GD[i+1][j  ][k  ][z_] - Bc_GD[i  ][j  ][k  ][z_] +
			       Bc_GD[i+1][j+1][k  ][z_] - Bc_GD[i  ][j+1][k  ][z_]);

	  compYdx = 0.5*invdx*( 
			       Bc_GD[i+1][j  ][k  ][y_] - Bc_GD[i  ][j  ][k  ][y_] +
			       Bc_GD[i+1][j+1][k  ][y_] - Bc_GD[i  ][j+1][k  ][y_]);

	  State_GV[i+1][j+1][k][iJx] =   compZdy;
	  State_GV[i+1][j+1][k][iJy] = - compZdx;
	  State_GV[i+1][j+1][k][iJz] =   compYdx - compXdy;
	 
	  }
	}

      // fill in the ghost cell, constant value
      for(j=1;j<nynLG-1;j++){
	State_GV[0][j][k][iJx] = State_GV[1][j][k][iJx];
	State_GV[0][j][k][iJy] = State_GV[1][j][k][iJy];
	State_GV[0][j][k][iJz] = State_GV[1][j][k][iJz];

	State_GV[nxnLG-1][j][k][iJx] = State_GV[nxnLG-2][j][k][iJx];
	State_GV[nxnLG-1][j][k][iJy] = State_GV[nxnLG-2][j][k][iJy];
	State_GV[nxnLG-1][j][k][iJz] = State_GV[nxnLG-2][j][k][iJz];
      }

      for(i=0;i<nxnLG;i++){
	State_GV[i][0][k][iJx] = State_GV[i][1][k][iJx];
	State_GV[i][0][k][iJy] = State_GV[i][1][k][iJy];
	State_GV[i][0][k][iJz] = State_GV[i][1][k][iJz];

	State_GV[i][nynLG-1][k][iJx] = State_GV[i][nynLG-2][k][iJx];
	State_GV[i][nynLG-1][k][iJy] = State_GV[i][nynLG-2][k][iJy];
	State_GV[i][nynLG-1][k][iJz] = State_GV[i][nynLG-2][k][iJz];
      }
  
    }
    else if(INdim == 3){
      for(i = 0; i < nxnLG-2; i++)
        for(j = 0; j < nynLG-2; j++)
	  for(k = 0; k < nznLG-2; k++){
	    if(doGetFromGM(i,j,k)){
	    compZdy = 0.25*invdy*( 
				  Bc_GD[i  ][j+1][k  ][z_] - Bc_GD[i  ][j][k  ][z_] +
				  Bc_GD[i+1][j+1][k  ][z_] - Bc_GD[i+1][j][k  ][z_] +
				  Bc_GD[i  ][j+1][k+1][z_] - Bc_GD[i  ][j][k+1][z_] +
				  Bc_GD[i+1][j+1][k+1][z_] - Bc_GD[i+1][j][k+1][z_]);

	    compXdy = 0.25*invdy*( 
				  Bc_GD[i  ][j+1][k  ][x_] - Bc_GD[i  ][j][k  ][x_] +
				  Bc_GD[i+1][j+1][k  ][x_] - Bc_GD[i+1][j][k  ][x_] +
				  Bc_GD[i  ][j+1][k+1][x_] - Bc_GD[i  ][j][k+1][x_] +
				  Bc_GD[i+1][j+1][k+1][x_] - Bc_GD[i+1][j][k+1][x_]);

	    compYdz = 0.25*invdz*( 
				  Bc_GD[i  ][j  ][k+1][y_] - Bc_GD[i  ][j  ][k  ][y_] +
				  Bc_GD[i+1][j  ][k+1][y_] - Bc_GD[i+1][j  ][k  ][y_] +
				  Bc_GD[i  ][j+1][k+1][y_] - Bc_GD[i  ][j+1][k  ][y_] +
				  Bc_GD[i+1][j+1][k+1][y_] - Bc_GD[i+1][j+1][k  ][y_]);

	    compXdz = 0.25*invdz*( 
				  Bc_GD[i  ][j  ][k+1][x_] - Bc_GD[i  ][j  ][k  ][x_] +
				  Bc_GD[i+1][j  ][k+1][x_] - Bc_GD[i+1][j  ][k  ][x_] +
				  Bc_GD[i  ][j+1][k+1][x_] - Bc_GD[i  ][j+1][k  ][x_] +
				  Bc_GD[i+1][j+1][k+1][x_] - Bc_GD[i+1][j+1][k  ][x_]);

	    compZdx = 0.25*invdx*( 
				  Bc_GD[i+1][j  ][k  ][z_] - Bc_GD[i  ][j  ][k  ][z_] +
				  Bc_GD[i+1][j+1][k  ][z_] - Bc_GD[i  ][j+1][k  ][z_] +
				  Bc_GD[i+1][j  ][k+1][z_] - Bc_GD[i  ][j  ][k+1][z_] +
				  Bc_GD[i+1][j+1][k+1][z_] - Bc_GD[i  ][j+1][k+1][z_]);

	    compYdx = 0.25*invdx*( 
				  Bc_GD[i+1][j  ][k  ][y_] - Bc_GD[i  ][j  ][k  ][y_] +
				  Bc_GD[i+1][j+1][k  ][y_] - Bc_GD[i  ][j+1][k  ][y_] +
				  Bc_GD[i+1][j  ][k+1][y_] - Bc_GD[i  ][j  ][k+1][y_] +
				  Bc_GD[i+1][j+1][k+1][y_] - Bc_GD[i  ][j+1][k+1][y_]);

	    State_GV[i+1][j+1][k+1][iJx] =  compZdy - compYdz;
	    State_GV[i+1][j+1][k+1][iJy] =  compXdz - compZdx;
	    State_GV[i+1][j+1][k+1][iJz] =  compYdx - compXdy;
	    }   
	  }
      // fill in the ghost cell, constant value

      for(i=1;i<nxnLG-1;i++)
	for(j=1;j<nynLG-1;j++){
	  State_GV[i][j][0][iJx] = State_GV[i][j][1][iJx];
	  State_GV[i][j][0][iJy] = State_GV[i][j][1][iJy];
	  State_GV[i][j][0][iJz] = State_GV[i][j][1][iJz];

	  State_GV[i][j][nznLG-1][iJx] = State_GV[i][j][nznLG-2][iJx];
	  State_GV[i][j][nznLG-1][iJy] = State_GV[i][j][nznLG-2][iJy];
	  State_GV[i][j][nznLG-1][iJz] = State_GV[i][j][nznLG-2][iJz];
        }

      for(i=1;i<nxnLG-1;i++)
	for(k=0;k<nznLG;k++){
	  State_GV[i][0][k][iJx] = State_GV[i][1][k][iJx];
	  State_GV[i][0][k][iJy] = State_GV[i][1][k][iJy];
	  State_GV[i][0][k][iJz] = State_GV[i][1][k][iJz];

	  State_GV[i][nynLG-1][k][iJx] = State_GV[i][nynLG-2][k][iJx];
	  State_GV[i][nynLG-1][k][iJy] = State_GV[i][nynLG-2][k][iJy];
	  State_GV[i][nynLG-1][k][iJz] = State_GV[i][nynLG-2][k][iJz];
	}
      
      for(j=0;j<nynLG;j++)
	for(k=0;k<nznLG;k++){
	  State_GV[0][j][k][iJx] = State_GV[1][j][k][iJx];
	  State_GV[0][j][k][iJy] = State_GV[1][j][k][iJy];
	  State_GV[0][j][k][iJz] = State_GV[1][j][k][iJz];

	  State_GV[nxnLG-1][j][k][iJx] = State_GV[nxnLG-2][j][k][iJx];
	  State_GV[nxnLG-1][j][k][iJy] = State_GV[nxnLG-2][j][k][iJy];
	  State_GV[nxnLG-1][j][k][iJz] = State_GV[nxnLG-2][j][k][iJz];
	}      
    }

    ReNormVariables();    
    Moment2Velocity();
    isFirstTime = false;
  }
 

  /** print info for coupling */
  void PrintInterfaceFluid(){

    if(thisrank == 0){
      cout << "nS = " << nS <<endl;
      cout << "Sum all particle masses = " << SumMass <<endl;
      cout << "useMultiFluid   = "<<useMultiFluid
	   <<" nIonFluid = "<<nIonFluid<<endl;
      cout << "useMultiSpecies = "<<useMultiSpecies
	   <<" nSpecies ="<<nSpecies<<endl;
      for(int is=0; is < nS; is++){
          cout << "Q/Qi[" << is << "] = "<< QoQi_S[is] << endl;
          cout << "M/Mi[" << is << "] = "<< MoMi_S[is] << endl;
      }
      if(!useMhdPe)cout << "Pe/Ptotal = " << PeRatio<<endl;
      cout << "========== Unit conversion factors ==============" << endl;
      cout << " Si2NoRho = " << Si2NoRho << endl;
      cout << " Si2NoV   = " << Si2NoV   << endl;
      cout << " Si2NoB   = " << Si2NoB   << endl;
      cout << " Si2NoP   = " << Si2NoP   << endl;
      cout << " Si2NoJ   = " << Si2NoJ   << endl;
      cout << " Si2NoL   = " << Si2NoL << endl;
      cout << "===================================================" << endl;

    }

  }


  /** electron mass given by IPIC3D params while ion mass comes form BATSRUS */
  void fixPARAM(double*& qom, int*& npcelx, int*& npcely, int*& npcelz, int *ns){

    double qom_el;
    int npx, npy, npz;

    SumMass = 0.0;

    //tmp store old values 
    qom_el = qom[0];
    npx = npcelx[0];
    npy = npcely[0];
    npz = npcelz[0];

    delete[] qom;
    delete[] npcelx;
    delete[] npcely;
    delete[] npcelz;

    *ns =nS;
    qom    = new double[nS];
    npcelx = new int[nS];
    npcely = new int[nS];
    npcelz = new int[nS];

    // electron
    MoMi_S[0] = QoQi_S[0]/qom_el;
 
    //iones
    for(int is = 0; is < nS; is++){
      qom[is]    = QoQi_S[is]/MoMi_S[is];
      npcelx[is] = npx;
      npcely[is] = npy;
      npcelz[is] = npz;
    }

    for(int is = 0; is < nS; is++)
      SumMass   += MoMi_S[is];


    if(thisrank == 0 ){
      cout << "============= fixPARAM =============" << endl;
      for(int is = 0; is < nS; is++){
        cout << "qom[" << is<< "]  =  " << qom[is] << endl;
        cout << "npcelx[" << is << "] = " << npcelx[is] << endl;
        cout << "npcely[" << is << "] = " << npcely[is] << endl;
        cout << "npcelz[" << is << "] = " << npcelz[is] << endl;
        cout << "Q/Qi[" << is << "] = "<< QoQi_S[is] << endl;
        cout << "M/Mi[" << is << "] = "<< MoMi_S[is] << endl;
      }
      if(!useMhdPe)cout << "Pe/Ptotal = " << PeRatio<<endl;      
      //cout << "===================================" << endl;
    }	
  

  }

  /** get start index in for the total domain */
  int getGlobalStartIndex(int dir){return(StartIdx_D[dir]);}
  
  
  /** get number of cell in overlap region for field solver, used by fields */
  int getnOverlap() {return(nOverlap);}
  void setnOverlap(double param) {nOverlap = param;}
  
  
  /** get number of cell in overlap region for quintities relativ for paricle initialization  */
  int getnOverlapP() {return(nOverlapP);}
  void setnOverlapP(double param) {nOverlapP = param;}
  
  
  /** get number of cell from the boundary where we push the solution closer to a fluid solution, 
      used by the Jxyz and P tensor */
  int getnIsotropic() const{return(nIsotropic);}
  void setnIsotropic(double param) {nIsotropic = param;}
  
  
  /** get number of cell from the boundary where we push the solution closer to a fluid solution, 
      used by charge density */
  int getnCharge() const{return(nCharge);}
  void setnCharge(double param) {nCharge = param;}
  
  
  /** Number of steps between syncrinosation of the fluid solution */
  int getnSync() {return(nSync);}
  void setnSync(int SyncDn) {nSync = SyncDn;}
  
  /** Use anisotropisc pressure when seting upt the particle distribution */
  bool getUseAnisoP()const{return(useAnisoP);}
  //void setAnisoP(bool useAnisoPIn){useAnisoP = useAnisoPIn;}

  /** Whether electron pressure is passed between PIC and BATSRUS */
  bool getUseMhdPe()const{return(useMhdPe);}

  bool getUseMultiFluid()const{return(useMultiFluid);}
  bool getUseMultiSpecies()const{return(useMultiSpecies);}
  
  /** Get convertion factor to from IPIC3D internal units */
  inline double getNo2Si_V(int idx){return(No2Si_V[idx]);}
  
  /** Get convertion factor to from IPIC3D internal units */
  inline double getSi2No_V(int idx){return(Si2No_V[idx]);}

  /** Get time convertion units to internal IPIC3D units */
  double gettUnitPic() {return(tUnitPic);}
  
  
  /** Get time convertion units from internal IPIC3D units */
  double getinvtUnitPic() {return(invtUnitPic);}; 
  
  /** Store pointer to processor topolegy object internaly */
  void setVCTpointer(VCtopology3D *vct){ _vct = vct;};


  /** Get time convertion units from internal IPIC3D units */
  void setSItime(double time) {
    SItime = time;
  }; 
  
  /** Get time convertion units from internal IPIC3D units */
  double getSItime() {return(SItime);}; 
  
  /** set normalized dt */
  void setNormDt(double normDt){
    dt   = normDt;
    INdt = normDt *(No2SiL/No2SiV);
    SItime += INdt;
  }
  
  /** get fluid time step in IPIC3D units */
  /** set SI dt */
  void setSIDt(double SIDt){
    INdt = SIDt;
    SItime += INdt;
    dt      = INdt*(Si2NoL/Si2NoV);
  }
  
  /** get fluid time step in IPIC3D units */
  double getFluidDt() const{
    if(dt == 0.0){
      //cout<<"getFluidDt : "<<dt<<", "<<ParamDt<<endl;
      return(ParamDt);
    }
    else{
      //cout<<"getFluidDt : "<<dt<<endl;
      return(dt);
    }
  }

  inline void divide_processors(int &npx, int &npy, int &npz, int nprocs){
    int divisor1, divisor2, divisor3, nprocs0, npmax, nDim;
    nDim = INdim;

    if(nDim<3){
      divisor1 = 1;
    }else{
      npmax = (int) pow(nprocs+1,1./3);
      for(int i = 1; i<=npmax; ++i){
	if(nprocs % i == 0) divisor1 = i;

      }
    }

    if(nDim == 1){
      divisor2 = 1;
      nprocs0 = nprocs;
    }else{
      nprocs0 = nprocs/divisor1;
      npmax = (int) pow(nprocs0+1,0.5);
      for(int i = 1; i<=npmax; ++i){
	if(nprocs0 % i == 0) divisor2 = i;
      }

    }
    divisor3 = nprocs0/divisor2;

    // divisor3 is alway the largest one;
    npx = divisor3;
    npy = divisor1 > divisor2 ? divisor1 : divisor2;
    npz = divisor1 > divisor2 ? divisor2 : divisor1;
  }
  
  // The begining 'physical' point of this IPIC region. Assume there is one 
  // layer PIC ghost cell.
  double getPhyXMin(){return(INxRange_I[0] + INgridDx_D[0]);}
  double getPhyYMin(){return(INxRange_I[2] + INgridDx_D[1]);}
  double getPhyZMin(){return(INxRange_I[4] + INgridDx_D[2]);}
  double getPhyXMax(){return(INxRange_I[1] - INgridDx_D[0]);}
  double getPhyYMax(){return(INxRange_I[3] - INgridDx_D[1]);}
  double getPhyZMax(){return(INxRange_I[5] - INgridDx_D[2]);}


  void setiRegion(int i){iRegion = i;}
  int getiRegion()const{return(iRegion);}
  
  int getnDim()const{return(INdim);}

  int getnVarFluid(){return(nVarFluid);}

  int getnIon(){return(nIon);}

  double getSi2NoL(){return(Si2NoL);}
  double getNo2SiL(){return(No2SiL);}
  double getNo2SiRho(){return(1./Si2NoRho);}
  double getNo2SiV(){return(1./Si2NoV);}
  double getNo2SiB(){return(1./Si2NoB);}
  double getNo2SiP(){return(1./Si2NoP);}
  double getNo2SiJ(){return(1./Si2NoJ);}

  int getNxcLocal(){return nxcLocal;}
  int getNycLocal(){return nycLocal;}
  int getNzcLocal(){return nzcLocal;}

  bool getdoNeedBCOnly(){return(doNeedBCOnly);};
  void setdoNeedBCOnly(bool doBCOnly){doNeedBCOnly=doBCOnly;};
  void setCycle(int iCycleIn){iCycle=iCycleIn;};
  int getCycle()const{return iCycle;};
  bool getUseRandomPerCell()const{return useRandomPerCell;};
  string getTestFunc()const{return testFuncs;};
  int getiTest()const{return iTest;};
  int getjTest()const{return jTest;};
  int getkTest()const{return kTest;};

  // IDL format output.
  int getnPlotFile()const{return nPlotFile;};
  int getdnOutput(int i)const{return dnOutput_I[i];};
  double getdtOutput(int i)const{return dtOutput_I[i];};
  string getplotString(int i)const{return plotString_I[i];};
  double getplotDx(int i)const{return plotDx_I[i];};
  string getplotVar(int i)const{return plotVar_I[i];};
  double getplotRange(int iPlot,int i)const{
    // plotRange_ID is only set from PARAM.in for 'cut'!!!!
    return plotRange_ID[iPlot][i];
  };
  bool getdoSaveBinary()const{return doSaveBinary;};
};


#endif 
