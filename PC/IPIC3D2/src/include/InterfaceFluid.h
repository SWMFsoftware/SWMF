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
#include <vector>
#include <array>
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

  static const int iErr = 11; 
  
  int nDim;        // number of dimentions

  // The meaning of INdt is not clear -Yuxi
  double INdt;

  // Min and Max of the physical domain in normalized PIC units. 
  double phyMin_D[3], phyMax_D[3];

  // In normalized PIC units, but in MHD coordinates.
  double gstMin_D[3], gstMax_D[3];

  // The length of the computational domain in normalized PIC units, including
  // the ghost cell layers. 
  double lenGst_D[3];
  
  // Cell Size
  double dx_D[3];

  // Rotation matrix. 
  double R_DD[3][3];

  bool doRotate;     

  // Number of cells/nodes in each direction, including the ghost cell layers. 
  int nCellGst_D[3], nNodeGst_D[3];
    
  double SItime; // time in SI units
  
  double ****Bc_GD;    // cell centered B 
  double ****State_GV; // node centered state variables

  // Number of variables passing between MHD and PIC. 
  int nVarFluid;

  // Number of fluid at the MHD side. One 'fluid' has its own density,
  // velocity and pressure. Electron can be one fluid. 
  int nFluid;

  // Number of ion fluid at the MHD side. 
  int nIonFluid;

  // Number of species at the MHD side. One 'species' only has its own density.
  int nSpecies;
  
  // Total number of ion/electron species exit in the fluid code. 
  int nIon;

  int nVarCoupling;
  
  bool useMultiSpecies, useMultiFluid, useElectronFluid;

  // storage for starting/ending physical (not include ghost cell)
  // cell indexes of this processor
  int StartIdx_D[3], EndIdx_D[3];  
  
  static const int NG = 1; // number of ghost cell

  bool useAnisoP;  // Use anisotripic pressure

  bool useMhdPe;

  double rPlanetSi; 
  
  double dt;            // scaled time step from fluid model
  double ParamDt;       // dt read in from input file

  double tUnitPic;      // conversenfactor to time used in IPIC3D 
  double invtUnitPic;   // 1/tUnitPic 
  double *Si2No_V; // array storing unit conversion factors
  double *No2Si_V; // array storing inverse unit conversion factors
  double Si2NoM, Si2NoV, Si2NoRho, Si2NoB, Si2NoP, Si2NoJ, Si2NoL, Si2NoE;
  double No2SiV, No2SiL;
  double MhdNo2SiL; // Length in BATSRUS normalized unit -> Si
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
  string sRegion;

  // Do not include ghost cells.
  int nxcLocal, nycLocal, nzcLocal;

  // Nodes, include ghost cells.
  int nxnLG, nynLG, nznLG;

  // The range of the computtational domain on this processor.
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  double *xStart_I, *xEnd_I, *yStart_I, *yEnd_I, *zStart_I, *zEnd_I;
  
  // Unless it is the first step to initilize PC, otherwise, only boundary
  // information is needed from GM. 
  bool doNeedBCOnly;

  int iCycle;

 protected:
  static const int x_=0, y_=1, z_=2;

  bool doSubCycling;
  
  // Variables for IDL format output.
  static const int nDimMax=3;
  int nPlotFile;
  int *dnOutput_I;
  double *dtOutput_I, *plotDx_I;
  //The second dimension: xmin, xmax, ymin, ymax, zmin, zmax. 
  double **plotRangeMin_ID, **plotRangeMax_ID; 
  string *plotString_I;
  string *plotVar_I;
  bool doSaveBinary;
  double drSat; // A particle within drSat*dx from a satellite point will be wrote out. 
  
  vector< vector< array<double,4> > > satInfo_III;

  // Simulation start time.
  int iYear, iMonth, iDay, iHour, iMinute, iSecond;
  
  // 'CFL' condition: uth*dt/dx < 1 needs to be satisfied for all species.
  double cflLimit;
  double maxDt; // maxDt = min(dxi/uth, dyi/uth, dzi/uth), i=0...nspecies-1

  // If the maximum thermal velocity of one node exceeds maxUth, which is in
  // normalized PIC unit, then save the output and stop runing. 
  double maxUth; //

  // 1) If useSWMFDt is true, use the dt given by coupling frequency.
  // 2) If useSWMFDt is false and useFixedDt is true, use fixedDt, which is set with
  //    command #TIMESTEP.
  // 3) If both useSWMFDt and useFixedDt are false, calculate dt using the
  // 'CFL' condition.
  bool useSWMFDt, useFixedDt;
  double fixedDt; // In SI unit

  bool isPeriodicX, isPeriodicY, isPeriodicZ; // Use periodic BC in one direction?

  // Variables for test setup.
  bool doTestEMWave;
  double waveVec_D[3], phase0, amplE_D[3]; 
  
 public:
  // These variables are also used in PSKOutput.h
  int *iRho_I, *iRhoUx_I, *iRhoUy_I, *iRhoUz_I, iBx,iBy,iBz, iEx,iEy,iEz,
    iPe,*iPpar_I,*iP_I,iJx,iJy,iJz, *iUx_I, *iUy_I, *iUz_I, iRhoTotal;

  int nBCLayer;
  bool useRandomPerCell;
  bool doUseOldRestart;
  string testFuncs;
  int iTest, jTest,kTest;

  int nPartGhost;

  // Change smooth coefficient near the boundary.
  // Parameters 'SmoothNiter' and 'Smooth' are declared in Colective.h
  bool doSmoothAll; // Smooth jh and rhoh?. 
  double innerSmoothFactor, boundarySmoothFactor;
  double nBoundarySmooth;

  // At most 10 vectors are supported during the coupling.
  static const int nVecMax = 10; 
  int vecIdxStart_I[nVecMax], nVec;
  
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
  
  
  /** Convert from local index to global index. Only correct for State_GV because of the ignored
   dimension. */
 public:  inline void getGlobalIndex(const int il, const int jl, const int kl,
                             int *ig, int *jg, int *kg)const{
    *ig =StartIdx_D[0] + il - 1 ;      
    if(nNodeGst_D[1] == 1) *jg = 0; else *jg =StartIdx_D[1] + jl - 1;
    if(nNodeGst_D[2] == 1) *kg = 0; else *kg =StartIdx_D[2] + kl - 1;
  }

  /** Only correct for State_GV because of the ignored dimension. */
 public:  inline void getLocalIndex(int *i, int *j, int *k)const
  {
    *i = *i - StartIdx_D[0] + 1;      
    if(nNodeGst_D[1] == 1) *j = 0; else *j = *j - StartIdx_D[1] + 1;
    if(nNodeGst_D[2] == 1) *k = 0; else *k = *k - StartIdx_D[2] + 1;
  }
  
  inline void getInterpolatedValue(int i, int j,  int k, double *Var, int iVar)const{
    if(nNodeGst_D[1] == 1) j=0;
    if(nNodeGst_D[2] == 1) k=0;
    *Var = State_GV[i][j][k][iVar];
  }
  
  	
  /** Second order interpolation  for a given position */
  inline void getInterpolatedValue(const double x, const double y, const double  z, double *Var, int iVar)const
  {
    int i1,j1,k1,i2,j2,k2;
    double dx1,dx2,dy1,dy2,dz1,dz2;
  
    // Get the index of the for the nodes surounding the cell
    i1 = max(floor(x/dx_D[0]),0.0)+1; i2=i1+1;
    if(nNodeGst_D[1] == 1){j1=0; j2=0;} else{j1 = max(floor(y/dx_D[1]),0.0)+1; j2=j1+1;}
    if(nNodeGst_D[2] == 1){k1=0; k2=0;} else{k1 = max(floor(z/dx_D[2]),0.0)+1; k2=k1+1;}
  
    // Get distenc form the cell walls
    dx1 = x/dx_D[0] - i1+1; dx2 = 1.0-dx1;	
    if(nNodeGst_D[1] == 1){dy1=0.5; dy2=0.5;} else{dy1 = y/dx_D[1] - j1+1; dy2 = 1.0-dy1;}
    if(nNodeGst_D[2] == 1){dz1=0.5; dz2=0.5;} else{dz1 = z/dx_D[2] - k1+1; dz2 = 1.0-dz1;}

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
      cout.precision(15);
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
    Si2NoE   = 1e6; // 1 V/m = 1e6  statV/cm/c. 'c' is the speed of light in cgs unit.  
    
    // Normalization: CGS -> non dimensional cgs    
    Si2NoRho /= RHOnorm;
    Si2NoV   /= Unorm;
    Si2NoB   /= Bnorm;
    Si2NoP   /= Pnorm;
    Si2NoJ   /= Jnorm;
    Si2NoL   /= Lnorm;
    Si2NoE   /= (Bnorm*Unorm);

    Si2NoM = Si2NoRho*Si2NoV;

    No2SiV = 1./Si2NoV;
    No2SiL = 1./Si2NoL;
    
    Si2No_V[iBx] = Si2NoB;
    Si2No_V[iBy] = Si2NoB;
    Si2No_V[iBz] = Si2NoB;


    Si2No_V[iJx] = Si2NoJ;
    Si2No_V[iJy] = Si2NoJ;
    Si2No_V[iJz] = Si2NoJ;

    Si2No_V[iEx] = Si2NoE;
    Si2No_V[iEy] = Si2NoE;
    Si2No_V[iEz] = Si2NoE; 

    if(useMhdPe)Si2No_V[iPe] = Si2NoP;
    if(useMultiSpecies) Si2No_V[iRhoTotal] = Si2NoRho;
    
    for(int iFluid=0; iFluid<nFluid; ++iFluid){
      Si2No_V[iRho_I[iFluid]]     = Si2NoRho;
      Si2No_V[iRhoUx_I[iFluid]]   = Si2NoM;
      Si2No_V[iRhoUy_I[iFluid]]   = Si2NoM;
      Si2No_V[iRhoUz_I[iFluid]]   = Si2NoM;
      Si2No_V[iP_I[iFluid]]       = Si2NoP;
      if(useAnisoP)Si2No_V[iPpar_I[iFluid]] = Si2NoP;
    }

    // Get back to SI units
    for(int iVar=0;iVar<nVarCoupling;iVar++)
      No2Si_V[iVar] = 1.0/Si2No_V[iVar];
    }
 
  void ReNormLength()
  {
    // Normalization
    for(int i =0; i<3; i++){
      dx_D[i]         *= Si2NoL;
      lenGst_D[i]     *= Si2NoL;
      gstMin_D[i]     *= Si2NoL;
      gstMax_D[i]     *= Si2NoL;
      phyMin_D[i] = gstMin_D[i] + dx_D[i];
      phyMax_D[i] = gstMax_D[i] - dx_D[i];
    }
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

    for(i = 0; i < nxnLG; i++)
      for(j = 0; j < nynLG; j++)
	for(k = 0; k < nznLG; k++){
	  if(doGetFromGM(i,j,k)){
	    if(useMultiSpecies){
	      Rhot=0;
	      for(int iFluid=0; iFluid<nFluid; ++iFluid){
		// Rho = sum(Rhoi) + Rhoe;
		Rhot +=
		  State_GV[i][j][k][iRho_I[iFluid]]*(1+MoMi_S[0]/MoMi_S[iFluid+1]);
	      }// iFluid

	      State_GV[i][j][k][iUx_I[0]] /= Rhot;
	      State_GV[i][j][k][iUy_I[0]] /= Rhot;
	      State_GV[i][j][k][iUz_I[0]] /= Rhot;	  
	    }else{
	      for(int iFluid=0; iFluid<nFluid; ++iFluid){
		State_GV[i][j][k][iUx_I[iFluid]] /= State_GV[i][j][k][iRho_I[iFluid]];
		State_GV[i][j][k][iUy_I[iFluid]] /= State_GV[i][j][k][iRho_I[iFluid]];
		State_GV[i][j][k][iUz_I[iFluid]] /= State_GV[i][j][k][iRho_I[iFluid]];
	      }// iFluid
	    } // else
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

    delete xStart_I;
    delete xEnd_I;
    delete yStart_I;
    delete yEnd_I;
    delete zStart_I;
    delete zEnd_I;
    
    if(nPlotFile>0){
      delete [] dnOutput_I;
      delete [] dtOutput_I;
      delete [] plotDx_I;
      delete [] plotString_I;
      delete [] plotVar_I;
      delArr2(plotRangeMin_ID,nDimMax);
      delArr2(plotRangeMax_ID,nDimMax);
    }    
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
  void setGlobalStartIndex()
  {
    /**
       StartIdx_D could be either cell index or node index. 

       Example: A 2d case. 
       
       7 physical cells in x direction, and 2 processors are used. 
       1 cell and 1 processor in z direction. Do not consider y 
       direction in this case. 
       
      y|
       |
       --->x
       
       proc 1 
      ________________________
      |   |   |   |   |   |  |
      |___|___|___|___|___|__|
      |   |   |   |   |   |  |
      |___|___|___|___|___|__|
      |   |   |   |   |   |  |
      |___|___|___|___|___|__|                                
      0   1   2   3   4   5      global node index
        0   1   2   3   4        global cell index
      0   1   2   3   4   5      local node index
        0   1   2   3   4   5    local cell index

	x: StartIdx_D[0] = 1, nxcLocal = 4, EndIdx_D[0] = 4
	z: StartIdx_D[2] = 1, nxcLocal = 1, EndIdx_D[0] = 1
	


	proc 2
      ____________________
      |   |   |   |   |  |
      |___|___|___|___|__|
      |   |   |   |   |  |
      |___|___|___|___|__|
      |   |   |   |   |  |
      |___|___|___|___|__|
          5   6   7   8     global node index
            5   6   7       global cell index
      0   1   2   3   4  5  local node index
        0   1   2   3  4    local cell index


	x: StartIdx_D[0] = 5, nxcLocal = 3, EndIdx_D[0] = 7
	z: StartIdx_D[2] = 1, nxcLocal = 1, EndIdx_D[0] = 1


     **/
    
    bool isCrowdedX, isCrowdedY, isCrowdedZ;
        
    nxcLocal=1;
    nycLocal=1;
    nzcLocal=1;
    
    nxcLocal = getFluidNxc()/(double)_vct->getXLEN();
    if(nDim > 1) nycLocal = getFluidNyc()/(double)_vct->getYLEN();
    if(nDim > 2) nzcLocal = getFluidNzc()/(double)_vct->getZLEN();

    myrank = _vct-> getCartesian_rank();

    StartIdx_D[0] = _vct->getCoordinates(0)*nxcLocal + 1;    
    StartIdx_D[1] = _vct->getCoordinates(1)*nycLocal + 1;    
    StartIdx_D[2] = _vct->getCoordinates(2)*nzcLocal + 1;    

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

    EndIdx_D[0] += nxcLocal - 1;    
    if(nDim > 1) EndIdx_D[1] += nycLocal - 1;    
    if(nDim > 2) EndIdx_D[2] += nzcLocal - 1;

    nxnLG = nxcLocal + 3;
    nynLG = nycLocal + 3;
    nznLG = nzcLocal + 3; 
  }	


  void setSubDomainRange(){
    int xlen, ylen, zlen;
    xlen = _vct->getXLEN();
    ylen = _vct->getYLEN();
    zlen = _vct->getZLEN();

    xStart_I = new double[xlen];
    xEnd_I   = new double[xlen];
    yStart_I = new double[ylen];
    yEnd_I   = new double[ylen];
    zStart_I = new double[zlen];
    zEnd_I   = new double[zlen];

    double *xStart0_I, *xEnd0_I, *yStart0_I, *yEnd0_I, *zStart0_I, *zEnd0_I;
    xStart0_I = new double[xlen];
    xEnd0_I   = new double[xlen];
    yStart0_I = new double[ylen];
    yEnd0_I   = new double[ylen];
    zStart0_I = new double[zlen];
    zEnd0_I   = new double[zlen];

    for(int ix=0; ix<xlen; ix++){
      xStart_I[ix] = -1e99;
      xEnd_I[ix]   = -1e99;
      xStart0_I[ix] = -1e99;
      xEnd0_I[ix]   = -1e99; 

    }

    for(int iy=0; iy<ylen; iy++){
      yStart_I[iy] = -1e99;
      yEnd_I[iy]   = -1e99;
      yStart0_I[iy] = -1e99;
      yEnd0_I[iy]   = -1e99; 

    }

    for(int iz=0; iz<zlen; iz++){
      zStart_I[iz] = -1e99;
      zEnd_I[iz]   = -1e99;      
      zStart0_I[iz] = -1e99;
      zEnd0_I[iz]   = -1e99;
    }


    xStart0_I[_vct->getCoordinates(0)] = xStart;
    xEnd0_I[_vct->getCoordinates(0)] = xEnd;
    yStart0_I[_vct->getCoordinates(1)] = yStart;
    yEnd0_I[_vct->getCoordinates(1)] = yEnd;
    zStart0_I[_vct->getCoordinates(2)] = zStart;
    zEnd0_I[_vct->getCoordinates(2)] = zEnd;

    MPI_Allreduce(xStart0_I,xStart_I,xlen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);
    MPI_Allreduce(xEnd0_I,xEnd_I,xlen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);
    MPI_Allreduce(yStart0_I,yStart_I,ylen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);
    MPI_Allreduce(yEnd0_I,yEnd_I,ylen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);
    MPI_Allreduce(zStart0_I,zStart_I,zlen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);
    MPI_Allreduce(zEnd0_I,zEnd_I,zlen,MPI_DOUBLE,MPI_MAX,MPI_COMM_MYSIM);

    delete xStart0_I;
    delete xEnd0_I;
    delete yStart0_I;
    delete yEnd0_I;
    delete zStart0_I;
    delete zEnd0_I;

    bool doTest;
    doTest=false;
    if(doTest && myrank==0){
    for(int ix=0; ix<xlen; ix++){
      cout<<"ix= "<<ix<<" xstart= "<<xStart_I[ix]<<endl;
    }

    for(int iy=0; iy<ylen; iy++){
      cout<<"iy= "<<iy<<" ystart= "<<yStart_I[iy]<<endl;
    }

    for(int iz=0; iz<zlen; iz++){
      cout<<"iz= "<<iz<<" zstart= "<<zStart_I[iz]<<endl;
    }
    }
  }


  
  /** find which processor have with dordinate */
  void findProcForPoint(VCtopology3D *vct, int nPoint, double *Xyz_I, int *iProc_I){

    double invLx, invLy, invLz;
    int iPx, iPy, iPz;
    myrank = vct-> getCartesian_rank();

    for(int i =0; i<nPoint*nDim; i++)
      Xyz_I[i]*=Si2NoL; 
       
    // Estimated inversed length of a sub domain for each dimension. 
    invLx = (double)vct->getXLEN()/(dx_D[0]*getFluidNxc());
    invLy = (double)vct->getYLEN()/(dx_D[1]*getFluidNyc());
    invLz = (double)vct->getZLEN()/(dx_D[2]*getFluidNzc());


    for(int i = 0; i < nPoint; i++){
      double mhd_D[3],pic_D[3];
      mhd_D[0] = Xyz_I[i*nDim ];
      mhd_D[1] = Xyz_I[i*nDim+1];
      mhd_D[2] = Xyz_I[i*nDim+2];
      mhd_to_Pic_Vec(mhd_D,pic_D);

      if(isThisRun(pic_D)){
	double xp = pic_D[0];
	// Estimated process index in x direction. It can be shifted because the
	// domain size on different processor may different, so the
	// while loop is used to do the correction.
	iPx = (int)(xp*invLx); 
	while( xp > xEnd_I[iPx] || xp < xStart_I[iPx] ){
	    if(xp > xEnd_I[iPx]) iPx++;
	    if(xp < xStart_I[iPx])iPx--;
	  }

	double yp = pic_D[1];
	iPy = (int)(yp*invLy);
	while( yp > yEnd_I[iPy] || yp < yStart_I[iPy] ){
	    if(yp > yEnd_I[iPy]) iPy++;
	    if(yp < yStart_I[iPy])iPy--;
	  }

	  if(nDim<3)
	    iPz = 0;
	  else{
	    double zp = pic_D[2];
	    iPz = (int)(zp*invLz);
	    while( zp > zEnd_I[iPz] || zp < zStart_I[iPz] ){
	      if(zp > zEnd_I[iPz]) iPz++;
	      if(zp < zStart_I[iPz])iPz--;
	  }	  

	  }
	  
    	  iProc_I[i] = iPz + vct->getZLEN()*(iPy + iPx*vct->getYLEN());	  
    	}
      }

    bool doTest;
    doTest=false;
    if(doTest && myrank==0){
      for(int iPoint=0; iPoint<nPoint; iPoint++){
	cout<<"ipoint= "<<iPoint
	    <<" x = "<<Xyz_I[iPoint*nDim  ] - phyMin_D[x_]
	    <<" y = "<<Xyz_I[iPoint*nDim+1] - phyMin_D[y_]
	    <<" z = "<<Xyz_I[iPoint*nDim+2] - phyMin_D[z_]
	    <<" proc= "<<iProc_I[iPoint]
	    <<endl;
	  }
    }

    for(int i =0; i<nPoint*nDim; i++)
      Xyz_I[i]*=No2SiL;     

  }

  bool isThisProc(double *Pos_D, bool doDebug=false){
    bool isThisProcReturn;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double csmall; 
    csmall = 0;
      
    xmin = xStart - csmall;
    ymin = yStart - csmall;
    zmin = zStart - csmall;
    xmax = xEnd + csmall;
    ymax = yEnd + csmall;
    zmax = zEnd + csmall;

    isThisProcReturn = true;
    if(Pos_D[0] < xmin || Pos_D[0] > xmax) isThisProcReturn = false;
    if(nDim > 1){
      if(Pos_D[1] < ymin || Pos_D[1] > ymax) isThisProcReturn = false;
    }
    if(nDim > 2){
      if(Pos_D[2] < zmin || Pos_D[2] > zmax) isThisProcReturn = false;
    }

    if(doDebug && !isThisProcReturn){
      cout<<" xmin= "<<xmin<<" ymin= "<<ymin<<" zmin= "<<zmin
	  <<" xmax= "<<xmax<<" ymax= "<<ymax<<" zmax= "<<zmax
	  <<" x= "<<Pos_D[0]<<" y= "<<Pos_D[1]<<" z= "<<Pos_D[2]<<endl;
    }
    
    return isThisProcReturn;
  }
  

  /* Is Pos_D inside this PIC domain? */
  bool isThisRun(double *Pos_D){
    bool isThisRunReturn;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double csmall; 
    csmall = 1e-7;

    xmin = - csmall;
    ymin = - csmall;
    zmin = - csmall;
    xmax = lenGst_D[x_] - 2*dx_D[x_] + csmall;
    ymax = lenGst_D[y_] - 2*dx_D[y_] + csmall;
    zmax = lenGst_D[z_] - 2*dx_D[z_] + csmall;

    isThisRunReturn = true;
    if(Pos_D[0] < xmin || Pos_D[0] > xmax) isThisRunReturn = false;
    if(nDim > 1){
      if(Pos_D[1] < ymin || Pos_D[1] > ymax) isThisRunReturn = false;
    }
    if(nDim > 2){
      if(Pos_D[2] < zmin || Pos_D[2] > zmax) isThisRunReturn = false;
    }
    
    return isThisRunReturn;
  }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNxc() const{ return(nNodeGst_D[0]-3*NG); }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNyc() const{ if(nNodeGst_D[1] > 3*NG) return(nNodeGst_D[1]-3*NG); else return(1); }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNzc() const{ if(nNodeGst_D[2] > 3*NG) return(nNodeGst_D[2]-3*NG); else return(1); }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNxn() { return(nNodeGst_D[0]-2*NG); }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNyn() { if(nNodeGst_D[1] > 3*NG) return(nNodeGst_D[1]-2*NG); else return(1); }
  
  /** nNodeGst_D includes 1 guard/ghost cell layer... */
  inline double getFluidNzn() { if(nNodeGst_D[2] > 3*NG) return(nNodeGst_D[2]-2*NG); else return(1); }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLx(){ return(lenGst_D[0]-NG*2*dx_D[0]); }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLy()
  { 
    if(nNodeGst_D[1] > 3*NG) return(lenGst_D[1]-NG*2*dx_D[1]); 
    else return(lenGst_D[1]); 
  }
  
  /** Get physical dimetntions to the simulation domain, without ghost cells */
  inline double getFluidLz()
  { 
    if(nNodeGst_D[2] > 3*NG) return(lenGst_D[2]-NG*2*dx_D[2]); 
    else return(lenGst_D[2]); 
  }
  
  /** place holder for all variables like charge density that we know is 0.0 */
  template <typename Type>	
    inline double getFluidZero(const Type x, const Type y, const Type  z, const int is)const{
    return(0.0);
  }
  
  /** get Pxx from fluid */
  template <typename Type>	
    inline double getFluidPxx(const Type x, const Type y,const Type z, const int is)const{
    double Pxx;
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pxx = Pperp + (Ppar - Pperp)*Bx*Bx/Bt2;
      
    }else{
      Pxx = getFluidP(x,y,z,is);
    }
    return(QoQi_S[is]*(Pxx/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUx(x,y,z,is),2)));
  }
  
  /** get Pyy from fluid */
  template <typename Type>	
    inline double getFluidPyy(const Type x, const Type y,const Type z, const int is)const{
    double Pyy; 
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pyy = Pperp + (Ppar - Pperp)*By*By/Bt2;      
    }else{
      Pyy = getFluidP(x,y,z,is);
    }
    return(QoQi_S[is]*(Pyy/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUy(x,y,z,is),2)));	
  }
  
  /** get Pzz from fluid */
  template <typename Type>	
    inline double getFluidPzz(const Type x, const Type y,const Type z, const int is)const{
    double Pzz;
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pzz = Pperp + (Ppar - Pperp)*Bz*Bz/Bt2;
      
    }else{
      Pzz = getFluidP(x,y,z,is);
    }

    return(QoQi_S[is]*(Pzz/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*pow(getFluidUz(x,y,z,is),2)));
    
  }
  
  /** get Pxy from fluid */
  template <typename Type>	
    inline double getFluidPxy(const Type x, const Type y,const Type z, const int is)const{
    double Pxy;
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pxy = (Ppar - Pperp)*Bx*By/Bt2;
      
    }else{
      Pxy = 0; 
    }

    return(QoQi_S[is]*(Pxy/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*getFluidUx(x,y,z,is)*getFluidUy(x,y,z,is)));
  }
  
  /** get Pxz from fluid */
  template <typename Type>	
    inline double getFluidPxz(const Type x, const Type y,const Type z, const int is)const{
    double Pxz;
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pxz = (Ppar - Pperp)*Bx*Bz/Bt2;
      
    }else{
      Pxz = 0; 
    }

    return(QoQi_S[is]*(Pxz/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*getFluidUx(x,y,z,is)*getFluidUz(x,y,z,is)));
  }
  
  /** get Pyz from fluid */
  template <typename Type>	
    inline double getFluidPyz(const Type x, const Type y,const Type z, const int is)const{
    double Pyz;
    if(useAnisoP){
      double Bx, By, Bz, Bt2, Ppar, P, Pperp; 
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);
      Bt2 = Bx*Bx + By*By + Bz*Bz;

      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);

      Pyz = (Ppar - Pperp)*By*Bz/Bt2;
      
    }else{
      Pyz = 0; 
    }

    return(QoQi_S[is]*(Pyz/MoMi_S[is] + getFluidRhoNum(x,y,z,is)*getFluidUy(x,y,z,is)*getFluidUz(x,y,z,is)));
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

      if(useElectronFluid){
	getInterpolatedValue(x,y,z,&U, iU_I[is]);      
      }else if(useMultiFluid){
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
      double Bx,By,Bz,B,P,Ppar,Pperp,Uthperp,Uthpar,Uthperp1, Uthperp2;
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

      // Get number density and B at the particle position
      double ni = getFluidRhoNum(x,y,z,is);
      getInterpolatedValue(x,y,z,&Bx,iBx);	
      getInterpolatedValue(x,y,z,&By,iBy);	
      getInterpolatedValue(x,y,z,&Bz,iBz);	
  
      // Get Parallel and perpendicular presure
      Ppar = getFluidPpar(x,y,z,is);
      P    = getFluidP(x,y,z,is);
      Pperp = 0.5*(3.0*P - Ppar);
     
      // Get 3 vertors spaning the vector space
      norm_DD = MagneticBaseVectors(Bx,By,Bz);
  
      //Get the thermal verlocities
      prob  = sqrt(-2.0*log(1.0-.999999999*rand1));
      theta = 2.0*M_PI*rand2;
      Uthpar = sqrt(Ppar/(MoMi_S[is]*ni))*prob*cos(theta);
  
      prob  = sqrt(-2.0*log(1.0-.999999999*rand3));
      theta = 2.0*M_PI*rand4;
      Uthperp  = sqrt(Pperp/(MoMi_S[is]*ni))*prob;
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

    if(useElectronFluid){
      getInterpolatedValue(x,y,z,&P,iP_I[is]); 
    }else if(useMultiFluid){
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

    if(useElectronFluid){
      getInterpolatedValue(x,y,z,&P,iPpar_I[is]);          
    }else if(useMhdPe){
      if(is==0) getInterpolatedValue(x,y,z,&P,iPe); // Electron
      if(is==1) getInterpolatedValue(x,y,z,&P,iPpar_I[0]); //Ion
    }else{
      getInterpolatedValue(x,y,z,&P,iPpar_I[0]);
      if(is == 0) P *= PeRatio;
      else if(is == 1) P *= (1 - PeRatio);	
    }

    return(P);
  }
  
  /** Get max thermal velocity in the domain. Result is on root processor only */
  double getMaxFluidUth(const int is, const int dir)const
  {
    double maxVth, allmaxVth;

    maxVth=0.0;
    // The loop range seems not right. --Yuxi
    for(int i = 1; i<=nxcLocal; i++)
      for(int j = 1; j<=nycLocal; j++)
  	for(int k = 1; k<=nzcLocal; k++)
	  maxVth = max(maxVth, getFluidUth(i,j,k,is));
    allmaxVth = maxVth; 
    MPI_Reduce(&maxVth, &allmaxVth, 1, MPI_DOUBLE, MPI_MAX, 0 , MPI_COMM_MYSIM);
    return(allmaxVth);
  }


  /** this should return NUMBER DENSITY from fluid */
  template <typename Type>
    inline double getFluidRhoNum(const Type x,const Type  y, const Type z, const int is)const{
    double Rho, NumDens;
    
    if(useElectronFluid){
      getInterpolatedValue(x,y,z,&Rho,iRho_I[is]);
      NumDens = Rho/MoMi_S[is];	      
    }else if(useMultiFluid || useMultiSpecies){
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
      if(nNodeGst_D[1] == 1) j=0;
      if(nNodeGst_D[2] == 1) k=0;

      //cout<<"mhd: E = - Ue x B"<<endl;
      // Ue = Ui - J/ne
      if(useElectronFluid){
      (*Ex) = State_GV[i][j][k][iEx];
      (*Ey) = State_GV[i][j][k][iEy];
      (*Ez) = State_GV[i][j][k][iEz];      
      }else{
	(*Ex) = (getFluidUz(ii,jj,kk,0)*State_GV[i][j][k][iBy] -
		 getFluidUy(ii,jj,kk,0)*State_GV[i][j][k][iBz]);
	(*Ey) = (getFluidUx(ii,jj,kk,0)*State_GV[i][j][k][iBz] -
		 getFluidUz(ii,jj,kk,0)*State_GV[i][j][k][iBx]);
	(*Ez) = (getFluidUy(ii,jj,kk,0)*State_GV[i][j][k][iBx] -
		 getFluidUx(ii,jj,kk,0)*State_GV[i][j][k][iBy]);
      }
      (*Bx) = State_GV[i][j][k][iBx];
      (*By) = State_GV[i][j][k][iBy];
      (*Bz) = State_GV[i][j][k][iBz];
    }
  }
  
  
  // Data recived from SWMF coupler
  void ReadFromGMinit(int *paramint, double *griddim, double *paramreal, stringstream *ss){

    nDim       = paramint[0];
    nVarFluid   = paramint[2]; 
    nFluid   = paramint[3];
    nSpecies    = paramint[4];

    // c++ index starts from 0. So, minus 1. 
    iPe      = paramint[5] - 1;
    iBx      = paramint[6] - 1;
    iBy      = iBx + 1;
    iBz      = iBy + 1;

    iEx      = paramint[7] - 1;
    iEy      = iEx + 1;
    iEz      = iEy + 1;

    useElectronFluid = iEx > 1;

    if(useElectronFluid){
      nIonFluid = -1; // Do not distinguish between electrons and ions.
      nIon = -1; 
      nS = nFluid; 
    }else{
      nIonFluid = nFluid;
      nIon = nFluid + nSpecies - 1; // Assuming one electron species. 
      nS = nIon + 1; // + electron
    }
    
    useMultiFluid   = nIonFluid > 1;
    useMultiSpecies = nSpecies > 1;
    
    nVarCoupling= nVarFluid + 3; // nVarFluid + (Jx, Jy, Jz)
    
    iRho_I   = new int[nFluid];
    iRhoUx_I = new int[nFluid];
    iRhoUy_I = new int[nFluid];
    iRhoUz_I = new int[nFluid];
    iUx_I    = new int[nFluid];
    iUy_I    = new int[nFluid];
    iUz_I    = new int[nFluid];
    iPpar_I  = new int[nFluid];
    iP_I     = new int[nFluid];      

    int n = 8;
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
      // Not multi-species
      for (int iFluid = 0; iFluid < nFluid; ++iFluid)
	iRho_I[iFluid]     = paramint[n++] - 1;
      for (int iFluid = 0; iFluid < nFluid; ++iFluid){
	iRhoUx_I[iFluid]   = paramint[n++] - 1;  iUx_I[iFluid]  = iRhoUx_I[iFluid];
	iRhoUy_I[iFluid]   = iRhoUx_I[iFluid] + 1; iUy_I[iFluid]  = iRhoUy_I[iFluid];
	iRhoUz_I[iFluid]   = iRhoUx_I[iFluid] + 2; iUz_I[iFluid]  = iRhoUz_I[iFluid];
      }

      for (int iFluid = 0; iFluid < nFluid; ++iFluid)
	iPpar_I[iFluid]    = paramint[n++] - 1;
      for (int iFluid = 0; iFluid < nFluid; ++iFluid)
	iP_I[iFluid]       = paramint[n++] - 1;
    }

    nVec = nFluid + 1;    
    if(useElectronFluid) nVec ++; // + E field.
    if(nVec>nVecMax){
      if(myrank==0) cout<<"Error: nVec > nVecMax!!!!"<<endl;
      MPI_Abort(MPI_COMM_MYSIM,iErr);
    }
    for (int iVec = 0; iVec<nFluid; iVec++) vecIdxStart_I[iVec] = iRhoUx_I[iVec];
    vecIdxStart_I[nFluid] = iBx;
    if(useElectronFluid) vecIdxStart_I[nFluid+1] = iEx;

    // See GM/BATSRUS/src/ModExtraVariables.f90. 
    useAnisoP = iPpar_I[0] != 0;
    useMhdPe  = iPe != 0;

    iJx = nVarFluid;
    iJy = iJx + 1;
    iJz = iJx + 2;

    Si2No_V = new double[nVarCoupling];
    No2Si_V = new double[nVarCoupling];
    for(int i =0; i<nVarCoupling; i++) Si2No_V[i] = 1;
    
    n = 0; 
    for(int i =0; i<3; i++){
      gstMin_D[i] = griddim[n++];      // Lmin
      lenGst_D[i] = griddim[n++]; 
      dx_D[i]       = griddim[n++]; //dx
      gstMax_D[i] =  gstMin_D[i] + lenGst_D[i];
    }

    for (int i=0; i<3; i++){
      for(int j=0; j<3; j++){
    	R_DD[i][j] = griddim[n++];
      }      
    }

    doRotate = false;
    double csmall = 1e-7;
    for (int i=0; i<nDim; i++)
      if(fabs(R_DD[i][i]-1) > csmall){
	doRotate = true;
      }

    for(int i =0; i<3; i++)
      nNodeGst_D[i] = (int)(lenGst_D[i]/dx_D[i] +0.5);

    // Add 1 as we count the nodes not cells
    for(int i =0; i<nDim; i++)
      nNodeGst_D[i] += 1;
    
    QoQi_S = new double[nS];
    MoMi_S = new double[nS];


    /** Do not change the order of the following lines. */
    n = 0;
    if(useElectronFluid){
      for(int i = 0; i < nS; ++i){
	QoQi_S[i] = paramreal[n++];
	MoMi_S[i] = paramreal[n++];
      }      
    }else{
      QoQi_S[0] = -1.0;
      for(int i = 1; i < nS; ++i){
	QoQi_S[i] = paramreal[n++];
	MoMi_S[i] = paramreal[n++];
      }      
    }
    
    // Electron pressure ratio: Pe/Ptotal
    PeRatio   = paramreal[n++];
    Lnorm     = paramreal[n++]; 	
    Unorm     = paramreal[n++];
    Mnorm     = paramreal[n++];
    rPlanetSi = paramreal[n++];
    MhdNo2SiL = paramreal[n++];
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
    assert_ge(lenGst_D[0],0.0);
    assert_ge(lenGst_D[1],0.0);
    assert_ge(lenGst_D[2],0.0);

    assert_ge(dx_D[0],0.0);
    assert_ge(dx_D[1],0.0);
    assert_ge(dx_D[2],0.0);

    if(!useElectronFluid){
      assert_eq(QoQi_S[0],-1.0*QoQi_S[1]);
      for(int iIon = 1; iIon<nIon; ++iIon){
	if(QoQi_S[iIon+1] != QoQi_S[1]){
	  cout<<" So far, -qe = qi1 = qi2 = qi3... is assumed. QoQi_S[1]= "
	      <<QoQi_S[1]<<" and QoQi_S["<<iIon+1<<"]= "<<QoQi_S[iIon+1]
	      <<" are not valided!!"<<endl;
	  abort();
	}
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
      int ig = StartIdx_D[0] + i -1; int nxg = getFluidNxc()+3;
      int jg = StartIdx_D[1] + j -1; int nyg = getFluidNyc()+3;
      int kg = StartIdx_D[2] + k -1; int nzg = getFluidNzc()+3;
      minDn = min(ig, nxg -ig -1);
      if(nDim>1) minDn = min(minDn, min(jg, nyg - jg - 1));
      if(nDim>2) minDn = min(minDn, min(kg, nzg - kg - 1));
      if(minDn <nBCLayer) doGetGM = true;
    }else{
      doGetGM = true;
    }
    return doGetGM;
  }
  
  /** get the number of cordinate points used by the simulation. recived from GM */
  void GetNgridPnt(int &nPoint){
    double *Pos_I=nullptr, *state_I=nullptr;
    int *iPoint_I=nullptr;
    bool doCount, doGetPos, doSetData;
    
    doCount = true; doGetPos = false; doSetData = false;    
    get_region_points(doCount, doGetPos, doSetData, nPoint, Pos_I,
		      state_I, iPoint_I);
  }

  void get_region_points(bool doCount, bool doGetPos, bool doSetData,
			 int &nPoint, double *pos_I, double *state_I, int *iPoint_I){
    int i, j, k, iMin, iMax, jMin, jMax, kMin, kMax;
    double pic_D[3], mhd_D[3];

    iMin = 0; iMax = 1 ; jMin = 0; jMax = 1; kMin = 0; kMax = 1; 
    iMax = nxnLG;
    if(nDim >= 2) jMax = nynLG;
    if(nDim == 3) kMax = nznLG;

    int n = 0, ii = 0; 
    nPoint = 0;
    for(k = kMin; k < kMax; k++)
      for(j = jMin; j < jMax; j++)
	for(i = iMin; i < iMax; i++){
	  if(doGetFromGM(i,j,k)){

	    if(doCount) nPoint +=nDim;

	    if(doGetPos){
	      pic_D[x_] = (i + StartIdx_D[x_] - 2)*dx_D[x_];
	      pic_D[y_] = (j + StartIdx_D[y_] - 2)*dx_D[y_];
	      pic_D[z_] = (k + StartIdx_D[z_] - 2)*dx_D[z_]; 
	      pic_to_Mhd_Vec(pic_D,mhd_D);
	      for (int iDim=0; iDim<nDim; iDim++){
		pos_I[n++] = mhd_D[iDim]*No2SiL;
	      }

	    }// doGetPos

	    if(doSetData){
	      for(int iVar=0; iVar < nVarFluid; iVar++){
		int idx; 
		idx = iVar + nVarFluid*(iPoint_I[ii] - 1);
		State_GV[i][j][k][iVar] = state_I[idx];		
	      }	      
	      ii++;

	      
	      // Convert vectors from MHD coordinates to PIC coordinates----------
	      for (int iVec = 0; iVec<nVec; iVec++){
		int idx0, idx1;
		idx0 = vecIdxStart_I[iVec];
		idx1 = idx0 + nDim;
		for (int iVar = idx0; iVar< idx1; iVar++)
		  mhd_D[iVar - idx0] = State_GV[i][j][k][iVar];
		mhd_to_Pic_Vec(mhd_D,pic_D,true);
		for (int iVar = idx0; iVar< idx1; iVar++){
		  State_GV[i][j][k][iVar] = pic_D[iVar-idx0];
		}
	      }
	      //------------------------------------------------------------------
	      
	    } // doSetData
	    
	  }
	}// i j k 

    if(doSetData && nDim<3){
      // Fill in the data in the ignored dimension.
      if(nDim==1){
	int j0 = 0;      
	for(i = iMin; i < iMax; i++)
	  for(j = 1; j < nynLG; j++)
	    for(k = kMin; k < kMax; k++)
	      for( int iVar=0; iVar < nVarFluid; iVar++){
		State_GV[i][j][k][iVar] = State_GV[i][j0][k][iVar];
	      }      
      }
      
      int k0 = 0;      
      for(i = iMin; i < iMax; i++)
	for(j = iMin; j < nynLG; j++)
	  for(k = 1; k < nznLG; k++)
	    for( int iVar=0; iVar < nVarFluid; iVar++){
	      State_GV[i][j][k][iVar] = State_GV[i][j][k0][iVar];
	    }      
    }
    
  }
  
  /** Get index of node points recived from GM */ 
  void GetGridPnt(double *Pos_I){
    double *state_I=nullptr;
    int *iPoint_I=nullptr, nPoint;
    bool doCount, doGetPos, doSetData;

    doCount = false; doGetPos = true; doSetData = false;    
    get_region_points(doCount, doGetPos, doSetData, nPoint, Pos_I,
		      state_I, iPoint_I);
  }


  void setStateVar(double *state_I, int *iPoint_I){
    const int x_=0, y_=1, z_=2;
    int i, j, k, iVar;
    int ii, jj, kk;
    int idx;

    if(isFirstTime){
      InitData();
      // Node values. 
      State_GV = newArr4(double,nxnLG,nynLG,nznLG,nVarCoupling);
      // Cell center values.
      Bc_GD = newArr4(double,nxnLG-1,nynLG-1,nznLG-1,3);
    }

    //-----------------Put data to State_GV--------------------
    double *Pos_I=nullptr;
    int nPoint;
    bool doCount, doGetPos, doSetData;
    doCount = false; doGetPos = false; doSetData = true;    
    get_region_points(doCount, doGetPos, doSetData, nPoint, Pos_I,
		      state_I, iPoint_I);
    //-----------------------------------------------------------

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

    // J = curl B/ mu0
    double invdx, invdy, invdz;
    double compZdx, compXdy, compYdz, compXdz, compZdy, compYdx;

    // 1/(dx*mu0)  
    invdx = 1.0/(dx_D[x_]*No2SiL*4.0*3.14159265359*1.0e-7);
    invdy = 1.0/(dx_D[y_]*No2SiL*4.0*3.14159265359*1.0e-7);
    invdz = 1.0/(dx_D[z_]*No2SiL*4.0*3.14159265359*1.0e-7);

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

	    if(nDim == 2 ){
	      for(int kk = 0; kk<=nznLG-1; kk += nznLG-1)
		for (int iVar = iJx; iVar <= iJz; iVar++){
		  State_GV[i+1][j+1][kk][iVar] = State_GV[i+1][j+1][k+1][iVar];
		}	      
	    }
	    
	  }   
	}

    // Fill in ghost cells. This part can be further simplified. - Yuxi
    if(nDim == 2){
      k = 0;
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
    else if(nDim == 3){
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
      cout <<" nS = " << nS <<endl;
      cout <<" Sum all particle masses = " << SumMass <<endl;
      cout << "useMultiFluid   = "<<useMultiFluid
	   <<" nFluid = "<<nFluid<<endl;
      cout <<" useMultiSpecies = "<<useMultiSpecies
	   <<" nSpecies ="<<nSpecies<<endl;
      cout <<" useElectronFluid = "<<useElectronFluid<<endl;
      for(int is=0; is < nS; is++){
          cout << "Q/Qi[" << is << "] = "<< QoQi_S[is] << endl;
          cout << "M/Mi[" << is << "] = "<< MoMi_S[is] << endl;
      }
      if(!useMhdPe)cout << "Pe/Ptotal = " << PeRatio<<endl;
      cout << "========== Unit conversion factors ==============" << endl;
      cout << " Si2NoRho = " << Si2NoRho << endl;
      cout << " Si2NoV   = " << Si2NoV   << endl;
      cout << " Si2NoB   = " << Si2NoB   << endl;
      cout << " Si2NoE   = " << Si2NoE   << endl;
      cout << " Si2NoP   = " << Si2NoP   << endl;
      cout << " Si2NoJ   = " << Si2NoJ   << endl;
      cout << " Si2NoL   = " << Si2NoL << endl;
      cout << "===================================================" << endl;

    }

  }


  // Read Satellite files. 
  void read_satellite_file(string filename){
    ifstream file;
    string line;
    bool doStartRead;
    vector< array<double,4> > satInfo_II;
    file.open(filename.c_str());
    int len;
    while(getline(file, line)){
      len = line.length();
      if(line.substr(0,6)=="#START") break;
    }

    int yr, mo, dy, hr, mn, sc;
    double msc; 
    array<double, 4> data;
    // A better way??
    while(file >> yr &&     //yr
	  file >> mo &&     //Mo
	  file >> dy &&     //Dy
	  file >> hr &&     //Hr
	  file >> mn &&     //Mn
	  file >> sc &&     //Sc
	  file >> msc &&     //Msc
	  file >> data[1] &&     //X
	  file >> data[2] &&     //Y
	  file >> data[3]        //Z
	  ){
      data[0] = convert_time(yr,mo,dy,hr,mn,sc,msc);

      // Assume the location of satellite is in normalized BATSRUS unit.It
      // is usually planet radius. 
      data[1] = data[1]*getMhdNo2NoL() - getFluidStartX();
      data[2] = data[2]*getMhdNo2NoL() - getFluidStartY();
      data[3] = data[3]*getMhdNo2NoL() - getFluidStartZ();

      satInfo_II.push_back(data);
    }
    
    bool doTestMe;
    doTestMe = true;
    int nline;
    nline = satInfo_II.size();
    for(int i=0; i<nline; i++){
      cout.precision(10);
    }            
    file.close();
    satInfo_III.push_back(satInfo_II);
  }

  void find_sat_points(double **pointList_ID, long &nPoint, int nPointMax,
		       double plotRange_I[6],
		       double xStart,double xEnd,double yStart,double yEnd,
		       double zStart,double zEnd){
    // Need to know which satellite. Here always use the last satellite
    // information read by read_satellite_file(), since this method is called
    // after read_sat_points() in iPic3dlib.cpp. Not a good approach!!!!
    
    nPoint = 0;
    int const iSat = satInfo_III.size()-1;
    int const nLine = satInfo_III[iSat].size();
    double x, y, z, xm, ym, zm;

    double dl; // The distance between two virual satellite output points.

    bool isFirstPoint;
    
    // dl = 0.25*min(dx,dy,dz); The parameter 0.25 can be changed. 0.5 may be
    // good enough. 
    dl = dx_D[0]<dx_D[1]? dx_D[0]:dx_D[1];
    if(nDim>2) dl = dl<dx_D[2]? dl:dx_D[2];
    dl *=0.25;

    isFirstPoint=true;
    for(int iLine = 0; iLine<nLine; iLine++){
      x = satInfo_III[iSat][iLine][1];
      y = satInfo_III[iSat][iLine][2];
      z = satInfo_III[iSat][iLine][3];
      if(nDim<3) z = 0;


      
      if(x>=0 && x<=getFluidLx() &&
	 y>=0 && y<=getFluidLy() &&
	 z>=0 && z<=getFluidLz()){ // inside the computational domain?
	
	if(isFirstPoint){
	  isFirstPoint=false;
	  xm = x;
	  ym = y;
	  zm = z;

	  plotRange_I[0] = xm; plotRange_I[1] = xm;
	  plotRange_I[2] = ym; plotRange_I[3] = ym;
	  plotRange_I[4] = zm; plotRange_I[5] = zm;	  
	  
	  if(xm>=xStart && xm<xEnd && ym>=yStart && ym<yEnd &&
	     zm>=zStart && zm<zEnd){	    
	    // This position is on this processor.
	    pointList_ID[nPoint][0] = xm;
	    pointList_ID[nPoint][1] = ym;
	    pointList_ID[nPoint][2] = zm;
	    nPoint++;
	  }

	}else{
	  double dl0;
	  dl0 = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym)+(z-zm)*(z-zm));
	  while(dl0 >=dl ){	    
	    xm = dl/dl0*x + (1-dl/dl0)*xm;
	    ym = dl/dl0*y + (1-dl/dl0)*ym;
	    zm = dl/dl0*z + (1-dl/dl0)*zm;

	    if(xm<plotRange_I[0]) plotRange_I[0] = xm;
	    if(xm>plotRange_I[1]) plotRange_I[1] = xm;
	    if(ym<plotRange_I[2]) plotRange_I[2] = ym;
	    if(ym>plotRange_I[3]) plotRange_I[3] = ym; 
	    if(zm<plotRange_I[4]) plotRange_I[4] = zm;
	    if(zm>plotRange_I[5]) plotRange_I[5] = zm; 
	    
	    if(xm>=xStart && xm<xEnd && ym>=yStart && ym<yEnd &&
	       zm>=zStart && zm<zEnd){
	      // This position is on this processor.
	      pointList_ID[nPoint][0] = xm;
	      pointList_ID[nPoint][1] = ym;
	      pointList_ID[nPoint][2] = zm;

	      nPoint++;
	      if(nPoint>nPointMax){
		cout<<"Error: nPoint = "<<nPoint<<" nPointMax= "<<nPointMax<<endl;
		abort();
	      }
	    }
	    
	    dl0 = sqrt((x-xm)*(x-xm)+(y-ym)*(y-ym)+(z-zm)*(z-zm));
	    
	  } // while
	}// if(nPoint==0)
      }// if
    }//for

    
    /* for(int iPoint=0; iPoint<nPoint;iPoint++){ */
    /*   cout<<"iPoint = "<<iPoint */
    /* 	  <<" x = "<<pointList_ID[iPoint][0] */
    /* 	  <<" y = "<<pointList_ID[iPoint][1] */
    /* 	  <<" z = "<<pointList_ID[iPoint][2] */
    /* 	  <<endl; */
    /* } */
    
  }

  int second_to_clock_time(int second){
    int iHr, iMn, iSc, time;
    iHr = floor(second/3600);
    second -= iHr*3600;
    iMn = floor(second/60);
    iSc = second - iMn*60;

    time = iSc + 100*iMn + 10000*iHr; 
    return time;
  }
    
  /** Convert real time into simulation time in second. **/
  double convert_time(int yr, int mo, int dy, int hr, int mn,int sc,double msc){
    double time; 
    double doyStart, doyNow; 
    doyStart = doy(iYear, iMonth, iDay);
    doyNow   = doy(yr, mo, dy);

    time = (((doyNow-doyStart)*24 + (hr - iHour))*60 + (mn - iMinute))*60 + sc - iSecond + msc/1000.0;

    bool doTestMe;
    doTestMe = false;
    if(doTestMe){
      cout<<" yr= "<<yr<<" mo= "<<mo<<" dy= "<<dy<<" hr= "<<hr<<" mn= "<<mn
	  <<" sc= "<<sc<<" msc="<<msc<<" iyr= "<<iYear<<" imo= "<<iMonth
	  <<" idy= "<<iDay<<" ihr= "<<iHour<<" imn= "<<iMinute
	  <<" isc= "<<iSecond<<" time= "<<time<<endl;
    }    
    return time; 
  }

  /** day of year **/
  int doy(int yr, int mo, int dy){
    int nDayInMonth_I [] ={31,28,31,30,31,30,31,31,30,31,30,31};
    int doy; 

    if(yr % 4 == 0) nDayInMonth_I[1]++;

    doy = 0;
    for(int i=0; i< mo-1; i++) doy += nDayInMonth_I[i];
    doy += dy;
    
    return doy;
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

    if(!useElectronFluid){
      MoMi_S[0] = QoQi_S[0]/qom_el;
    }
 
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
      if(!useMhdPe && !useElectronFluid)cout << "Pe/Ptotal = " << PeRatio<<endl;      
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
  bool get_useElectronFluid()const{return useElectronFluid;}
  
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
  
  /** set SI dt */
  void setSIDt(double SIDt, bool isSWMFDt){
    if(isSWMFDt && !useSWMFDt) return;
    
    INdt = SIDt;
    dt   = INdt*(Si2NoL/Si2NoV);
    return;
  }

  double calSIDt(){
    double dt0; 
    if(useSWMFDt){
      dt0 = INdt;
    }else{
      if(useFixedDt){
	dt0 = fixedDt;
      }else{
	dt0 = maxDt*cflLimit;
      }
    }
    return dt0;
  }
  
  void updateSItime(){
    SItime += INdt;
    if(myrank==0) {
      cout<<"SItime = "<<SItime<<" dt (s) = "<<INdt
	  <<" , normalized dt = "<<INdt*(Si2NoL/Si2NoV)<<endl;
    }	
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
    int divisor1, divisor2, divisor3, nprocs0, npmax;

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

  double getSmoothFactor(int i, int j, int k)const{
    double smooth; 
    int iMin,jMin,kMin,idxMin;

    int ig, jg, kg;
    getGlobalIndex(i,j,k,&ig,&jg,&kg);

    // The edge for PIC domain is node with index 1. It is different from
    // standalone PIC
    iMin = min(ig-1,nNodeGst_D[0]-ig-2);
    jMin = min(jg-1,nNodeGst_D[1]-jg-2);
    kMin = min(kg-1,nNodeGst_D[2]-kg-2);
    idxMin = min(iMin,jMin);

    if(nDim>2) idxMin = min(idxMin,kMin);
    
    double ix;
    if(idxMin<nBoundarySmooth){
      ix = (nBoundarySmooth-idxMin)/nBoundarySmooth; 
      smooth = ix*boundarySmoothFactor + (1-ix)*innerSmoothFactor;
    }else{
      smooth = innerSmoothFactor;
    }

    if(smooth > 1) smooth = 1;
    return smooth;    
  }

  string expandVariable(string inVars)const{
    // Expand the plot variables inside { };
    // Only support {fluid} so far.
    string::size_type pos1, pos2;
    string var0;
    
    pos1 = inVars.find_first_of("{");
    while(pos1 !=string::npos){
      pos2 = inVars.find_first_of("}");
      if(pos2 == string::npos){
  	cout<<"Variables should be inside { }: "<<inVars<<endl;
  	abort();
      }

      var0 = inVars.substr(pos1+1,pos2-pos1-1);
      inVars.erase(pos1,pos2-pos1+1);
      if(var0=="fluid"){
  	inVars += " rhoS0 rhoS1 Bx By Bz Ex Ey Ez uxS0 uyS0 uzS0 uxS1 uyS1 uzS1 pS0 pS1 pXXS0 pYYS0 pZZS0 pXYS0 pXZS0 pYZS0 pXXS1 pYYS1 pZZS1 pXYS1 pXZS1 pYZS1";
      }
      pos1 = inVars.find_first_of("{");
    }
    return inVars;
  }

  void pic_to_Mhd_Vec(double const *vecIn_D, double *vecOut_D, bool isZeroOrigin=false){
    /** 1) Change a vector in coupling PIC coordinates to MHD coordinates. 
	If not isZeroOrigin, then shifting the origin of the coupling 
	PIC coordinates to INxRange_I.
     **/

    for(int iDim=0; iDim < nDim; iDim++){
      vecOut_D[iDim] = 0;
      for(int jDim=0; jDim<nDim; jDim++){
	vecOut_D[iDim] +=R_DD[iDim][jDim]*vecIn_D[jDim];
      }
      if(!isZeroOrigin) vecOut_D[iDim] += phyMin_D[iDim];
    }    
  }
  
  void mhd_to_Pic_Vec(double const *vecIn_D, double *vecOut_D, bool isZeroOrigin=false){    
    /** 1) Change a vector in MHD coordinates to coupling PIC coordinates. 
	If not isZeroOrigin, then shifting the origin of the coupling 
	PIC coordinates to phyMin_D.
    **/
    
    double vec_D[3];
    if(!isZeroOrigin){
      for(int iDim=0; iDim < nDim; iDim++)
	vec_D[iDim] = vecIn_D[iDim] - phyMin_D[iDim];
    }else{
      for(int iDim=0; iDim < nDim; iDim++)
	vec_D[iDim] = vecIn_D[iDim];
    }
    
    for(int iDim=0; iDim < nDim; iDim++){
      vecOut_D[iDim] = 0;
      for(int jDim=0; jDim<nDim; jDim++){
	vecOut_D[iDim] +=R_DD[jDim][iDim]*vec_D[jDim];
      }
    } // iDim   
  }
  
  
  // The begining 'physical' point of this IPIC region. Assume there is one 
  // layer PIC ghost cell.
  /* double getPhyXMin()const{return phyMin_D[0];} */
  /* double getPhyYMin()const{return phyMin_D[1];} */
  /* double getPhyZMin()const{return phyMin_D[2];} */
  /* double getPhyXMax()const{return phyMax_D[0];} */
  /* double getPhyYMax()const{return phyMax_D[1];} */
  /* double getPhyZMax()const{return phyMax_D[2];} */
  double getphyMin(int i)const{return phyMin_D[i];}
  bool getdoRotate()const{return doRotate;}
  double getRDD(int i, int j)const{return R_DD[i][j];}

  void setiRegion(int i){
    iRegion = i;
    stringstream ss;
    sRegion=ss.str();
  }
  int getiRegion()const{return(iRegion);}
  string getsRegion()const{return sRegion;}
  
  int getnDim()const{return(nDim);}

  int getnVarFluid()const{return(nVarFluid);}
  int getnIon()const{return(nIon);}
  int get_nFluid()const{return(nFluid);}

  double getSi2NoL()const{return(Si2NoL);}
  double getNo2SiL()const{return(No2SiL);}
  double getNo2SiRho()const{return(1./Si2NoRho);}
  double getNo2SiV()const{return(1./Si2NoV);}
  double getNo2SiB()const{return(1./Si2NoB);}
  double getNo2SiP()const{return(1./Si2NoP);}
  double getNo2SiJ()const{return(1./Si2NoJ);}

  int getNxcLocal()const{return nxcLocal;}
  int getNycLocal()const{return nycLocal;}
  int getNzcLocal()const{return nzcLocal;}

  bool getdoNeedBCOnly()const{return(doNeedBCOnly);};
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
  double getplotDx(int i)const{return plotDx_I[i]>0 ? plotDx_I[i]:1;};
  string getplotVar(int i)const{return plotVar_I[i];};
  double getplotRangeMin(int iPlot,int i)const{
    // plotRangeMin_ID is only set from PARAM.in for 'cut'!!!!
    return plotRangeMin_ID[iPlot][i];
  };
  double getplotRangeMax(int iPlot,int i)const{
    // plotRangeMin_ID is only set from PARAM.in for 'cut'!!!!
    return plotRangeMax_ID[iPlot][i];
  };

  bool getdoSaveBinary()const{return doSaveBinary;};
  double getSatRadius()const{return drSat*dx_D[0];};
  void setmaxDt(double dt){maxDt = dt/(Si2NoL/Si2NoV);};
  void setxStart(double v){xStart=v;};
  void setxEnd(double v){xEnd=v;};
  void setyStart(double v){yStart=v;};
  void setyEnd(double v){yEnd=v;};
  void setzStart(double v){zStart=v;};
  void setzEnd(double v){zEnd=v;};
  int getnPartGhost()const{return nPartGhost;};
  bool getdoSmoothAll()const{return doSmoothAll;};
  double getMiSpecies(int i)const{return MoMi_S[i];};
  double getcLightSI()const{return Unorm/100;/*Unorm is in cgs unit*/};
  
  bool getdoTestEMWave()const{return doTestEMWave;};
  double getwaveVec_D(int i)const{return waveVec_D[i]/Si2NoL;};
  double getphase0deg()const{return phase0;};
  double getamplE_D(int i)const{return amplE_D[i];};

  // Replace these functions. --Yuxi
    /** return min X for fluid grid without ghostcell */  
  inline double getFluidStartX()const{return phyMin_D[x_];}
  /** return min Y for fluid grid without ghostcell */
  inline double getFluidStartY()const{return phyMin_D[y_];}
  /** return min Z for fluid grid without ghostcell */
  inline double getFluidStartZ()const{return phyMin_D[z_];}

  
  // return planet radius in SI unit.
  inline double getrPlanet()const{return(rPlanetSi);}
  // return MhdNo2SiL
  inline double getMhdNo2SiL()const{return(MhdNo2SiL);}
  // BATSRUS normalized unit -> PIC normalized unit;
  inline double getMhdNo2NoL()const{return(MhdNo2SiL*Si2NoL);}

  bool getdoSubCycling()const{return doSubCycling;}

  double get_maxUth()const{return maxUth;}
};


#endif 
