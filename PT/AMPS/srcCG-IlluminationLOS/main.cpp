//$Id$

#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>


#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "Comet.h"
#include "SpiceUsr.h"


double BulletLocalResolution(double *x) {
  int idim,i;
  double res,r,l[3];
  double SubsolarAngle;



  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  if (r<1600) return 300*12.0;

  SubsolarAngle=acos(l[0])*180.0/Pi;
  if (SubsolarAngle>=89.5) return 2000.0*pow(r/1600.0,2.0);

  return 300*pow(r/1600.0,2.0)*12.0;
 }

int SurfaceBoundaryCondition(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace) {
  int spec=PIC::ParticleBuffer::GetI(ptr);

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)  return _PARTICLE_DELETED_ON_THE_FACE_;
  else {
  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;
  //  return _PARTICLE_DELETED_ON_THE_FACE_;
  }
#else

  double c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];

  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];

  return _PARTICLE_REJECTED_ON_THE_FACE_;

#endif

}


double SurfaceResolution(CutCell::cTriangleFace* t) {
  return max(1.0,t->CharacteristicSize()*36.0); //18.0
}

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
    if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      ElectricallyChargedDust::EvaluateLocalTimeStep(spec,dt,startNode); //CharacteristicSpeed=3.0;
      return dt*3.0; //3.0
    }else CharacteristicSpeed=3.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
#else
    CharacteristicSpeed=3.0e2*sqrt(PIC::MolecularData::GetMass(_H2O_SPEC_)/PIC::MolecularData::GetMass(spec));
#endif

    CellSize=startNode->GetCharacteristicCellSize();
    return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double res=0.0;

  return res;
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  return false;
}

//injection of model particles through the faces of the bounding box
long int BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  long int nInjectedParticles=0;

  return nInjectedParticles;
}

long int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    long int nInjectedParticles=0;

  return nInjectedParticles;
}


double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

  return res;
}


double Comet::GetTotalProduction(int spec,void *BoundaryElement) {
#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__OFF_
  return Comet::GetTotalProductionRateBjornNASTRAN(spec)+Comet::GetTotalProductionRateUniformNASTRAN(spec)+Comet::GetTotalProductionRateJetNASTRAN(spec);
#else
  if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
   static bool InitFlag=false;
   static double MeanDustGrainMass=0.0;

   if (InitFlag==false) {
     InitFlag=true;

     //get mean mass of a dust grain
     const int nTotalTest=50;
     double r,c,cTotal=0.0,m=0.0;


     for (int ntest=0;ntest<nTotalTest;ntest++) {
       ElectricallyChargedDust::SizeDistribution::GenerateGrainRandomRadius(r,c);
       m+=pow(r,3.0)*c;
       cTotal+=c;
     }

     MeanDustGrainMass=m*4.0/3.0*Pi*ElectricallyChargedDust::MeanDustDensity/cTotal;
   }
   return ElectricallyChargedDust::TotalMassDustProductionRate/MeanDustGrainMass*500.0; //500.0
  }
  else return Comet::GetTotalProductionRateBjornNASTRAN(spec)+Comet::GetTotalProductionRateUniformNASTRAN(spec)+Comet::GetTotalProductionRateJetNASTRAN(spec);   
#endif

}

double Comet::GetTotalProduction(int spec,int BoundaryElementType,void *BoundaryElement) {
  return GetTotalProduction(spec,BoundaryElement);
}


int readDATAlength(long int& Data_length, long int File_Header, FILE *fH){
  char str[10000];
  Data_length=0;
  rewind(fH);
  while (!feof(fH)){
    fgets(str,10000,fH);
    Data_length++;
  }
  Data_length -= (1+File_Header);
  return 1;
}

int main(int argc,char **argv) {


  FILE *fTime;
  fTime = fopen("time_table.txt","r");
 
  string TimeTable[40000];

  long int timeFileLength=0;
  long int timeFileHeader=0;

  readDATAlength(timeFileLength, timeFileHeader, fTime);
  if (PIC::Mesh::mesh.ThisThread==0) printf("timeFileLength=%li  \n",timeFileLength);

  ifstream file("time_table.txt");
  int ct=0;
  int ichar=0;
  if(file.is_open()) {
    for (int i=0;i<timeFileLength;i++) {
      file >> TimeTable[i];
    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) cout << "TimeTable[0]" << TimeTable[0] << ", " << "TimeTable[1]" << TimeTable[1] << "\n";

  fclose(fTime);

  //init the particle solver
  PIC::InitMPI();
  PIC::Init_BeforeParser();
  Comet::Init_BeforeParser();

  double xmin[3]={0.0,-1.0,1.0};
  double xmax[3]={1.0,1.0,2.0};

  //load the NASTRAN mesh
  PIC::Mesh::IrregularSurface::ReadNastranSurfaceMeshLongFormat_km("cg-spc-shap5-v1.1-cheops_mod.bdf");

  PIC::Mesh::IrregularSurface::GetSurfaceSizeLimits(xmin,xmax);
  PIC::Mesh::IrregularSurface::PrintSurfaceTriangulationMesh("SurfaceTriangulation.dat",PIC::OutputDataFileDirectory);


#if _PIC_MODEL__3DGRAVITY__MODE_ == _PIC_MODEL__3DGRAVITY__MODE__ON_
  //Computation of the gravity field for an irregular nucleus shape
  nucleusGravity::readMesh_longformat("cg.RMOC-volume.nas");
  nucleusGravity::setDensity(430);
#endif   

  //set up the nastran object
  //init the nucleus
  cInternalBoundaryConditionsDescriptor CGSurfaceDescriptor;
  cInternalNastranSurfaceData *CG;

  CGSurfaceDescriptor=PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  CG=(cInternalNastranSurfaceData*) CGSurfaceDescriptor.BoundaryElement;

  CG->InjectionRate=Comet::GetTotalProduction;
  CG->faceat=0;

  Comet::GetNucleusNastranInfo(CG);

  for (int i=0;i<3;i++) xmin[i]=-0.5e5,xmax[i]=0.5e5;

  PIC::Mesh::mesh.CutCellSurfaceLocalResolution=SurfaceResolution;
  PIC::Mesh::mesh.AllowBlockAllocation=false;
  PIC::Mesh::mesh.init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh.memoryAllocationReport();

  //generate mesh or read from file                                                                                                                  
  char mesh[_MAX_STRING_LENGTH_PIC_]="amr.sig=0x73712d6492b3417e.mesh.bin";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";                                                              

  //generate mesh or read from file                                                                                                                  
  if (PIC::Mesh::mesh.ThisThread==0) printf("Generate Mesh \n");
  bool NewMeshGeneratedFlag=false;

  char fullname[STRING_LENGTH];
  sprintf(fullname,"%s",mesh);

  FILE *fmesh=NULL;

  fmesh=fopen(fullname,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh.readMeshFile(fullname);
  }
  else {
    NewMeshGeneratedFlag=true;
    if (PIC::Mesh::mesh.ThisThread==0) {
      PIC::Mesh::mesh.buildMesh();    
      PIC::Mesh::mesh.saveMeshFile("mesh.msh");
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
      PIC::Mesh::mesh.readMeshFile("mesh.msh");
    }
  }

  PIC::Mesh::mesh.SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh.CreateNewParallelDistributionLists();



  //initialize the blocks
  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh.AllowBlockAllocation=true;
  PIC::Mesh::mesh.AllocateTreeBlocks();

  PIC::Mesh::mesh.memoryAllocationReport();
  PIC::Mesh::mesh.GetMeshTreeStatistics();


  //init the volume of the cells'
  PIC::Mesh::IrregularSurface::CheckPointInsideDomain=PIC::Mesh::IrregularSurface::CheckPointInsideDomain_default;

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin                                                             
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh.getMeshSignature();

    if (PIC::Mesh::mesh.ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  PIC::Mesh::IrregularSurface::InitExternalNormalVector();

  //Define regions
  int nDev=25;
  int i,j,idim,harm;
  double x[3],xCenter[3],colatitude[200000],longitude[200000];

  for (i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
    CutCell::BoundaryTriangleFaces[i].GetCenterPosition(xCenter);
    colatitude[i]=acos(xCenter[2]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]+xCenter[2]*xCenter[2]));
    if(xCenter[1]>0.0) {
      longitude[i]=acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
    }else{
      longitude[i]=-acos(xCenter[0]/sqrt(xCenter[0]*xCenter[0]+xCenter[1]*xCenter[1]));
    }
  }

  double phiX=3.6669/2.0; 
  double phiY=3.6669/2.0;
  int nPixelsX=256;
  int nPixelsY=256;
  
  double lx=0.0,ly=0.0;
  double rPointingInstrument[3][nPixelsX*nPixelsY];
  double rPointing[3][nPixelsX*nPixelsY];
  double normPointing=0.0;
  double rotMat[3][3];

  
  double *sphericalHarmonic;
  sphericalHarmonic = new double [nDev];

  //Use illumination model
  if (PIC::Mesh::mesh.ThisThread==0) printf("Initialize Arrays \n");
  static double **C;

  C = new double* [nPixelsX*nPixelsY];
  C[0] = new double [nPixelsX*nPixelsY*nDev];
  for (i=0;i<nPixelsX*nPixelsY;i++) {
    C[i]=C[0]+i*nDev;
    for (j=0;j<nDev;j++) C[i][j]=0.0;
  }


  static double *Cvec;
  Cvec = new double [nPixelsX*nPixelsY*nDev];

  static double *buffer;
  buffer = new double [nPixelsX*nPixelsY*nDev];

  int ctGeom=0;

  for (i=0;i<nPixelsX*nPixelsY;i++) for (j=0;j<nDev;j++) C[i][j]=0.0,Cvec[i*nDev+j]=0.0,buffer[i*nDev+j]=0.0; //initialize C
  
  double a=0.02; 
  double g;

  MPI_Status status;
  int iStart,iFinish,nFaceThread,thread;
  nFaceThread=nPixelsX*nPixelsY/PIC::nTotalThreads;
  iStart=nFaceThread*PIC::ThisThread;
  iFinish=iStart+nFaceThread;
  if (PIC::ThisThread==PIC::nTotalThreads-1) iFinish=nPixelsX*nPixelsY;


  printf("PIC::ThisThread=%i iStart=%i iFinish=%i \n",PIC::ThisThread,iStart,iFinish); 

  FILE *fout;
  fout = fopen("rAU.dat","w"); 

  //init SPICE                                                                                                              
  SpiceDouble StateRosetta[6],StateSun[6],et,lt;
  double xObservation[3]={1.0E6,0,0},xPrimary[3]={0,0,0},xSecondary[3]={0,1.0E6,0};
  double xSun[3]={0,0,0};
  double xLightSource[3];
  double sunDistance=0.0, scDistance=0.0, rAU=0.0;
  

  furnsh_c("operational_1.tm");
  furnsh_c("operational_2.tm");
  furnsh_c("operational_3.tm"); 

  for (i=0;i<timeFileLength;i++) {
    
    utc2et_c(TimeTable[i].c_str(),&et);  
    spkpos_c("ROSETTA",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",StateRosetta,&lt);
    spkpos_c("SUN",et,"67P/C-G_CK","NONE","CHURYUMOV-GERASIMENKO",StateSun,&lt);
    
    for (idim=0;idim<3;idim++) {
      xObservation[idim]=1.0E3*StateRosetta[idim];
      xPrimary[idim]=0.0;
      xSecondary[idim]=1.0E3*StateSun[idim];
      xSun[idim]=1.0E3*StateSun[idim];
    }
    
    scDistance=sqrt(xObservation[0]*xObservation[0]+xObservation[1]*xObservation[1]+xObservation[2]*xObservation[2]);
    sunDistance=sqrt(xSun[0]*xSun[0]+xSun[1]*xSun[1]+xSun[2]*xSun[2]);
    rAU=sunDistance/_AU_;
    if (PIC::Mesh::mesh.ThisThread==0) {
      printf("xSun[0]=%e xSun[1]=%e xSun[2]=%e rAU=%e scDistance=%e\n",xSun[0],xSun[1],xSun[2],rAU,scDistance);
      printf("xObservation[0]=%e xObservation[1]=%e xObservation[2]=%e\n",xObservation[0],xObservation[1],xObservation[2]);
    }
    
    //shadow procedure    
    for (idim=0;idim<3;idim++) xLightSource[idim]=xSun[idim];
    PIC::RayTracing::SetCutCellShadowAttribute(xLightSource,false);    
    
    fprintf(fout,"% e\n",rAU);
    
    pxform_c("ROS_VIRTIS-M", "67P/C-G_CK", et, rotMat);
    
    //create pointing vectors
    lx=2*sin(phiX/180.0*Pi);
    ly=2*sin(phiY/180.0*Pi);
    
    int ct=0,nx=0,ny=0;
    double dx,dy;
    for (ny=0;ny<nPixelsY;ny++) {
      dy=-ly/2.0 + (double)ny/nPixelsY*ly;
      for(nx=0;nx<nPixelsX;nx++) {
	dx=-lx/2.0 + (double)nx/nPixelsX*lx;
	normPointing=sqrt(dx*dx+dy*dy+1.0);
	rPointingInstrument[0][ct]=dx/normPointing;
	rPointingInstrument[1][ct]=dy/normPointing;
	rPointingInstrument[2][ct]=1.0/normPointing;
	ct++;
      }
    }
    
    //rotate pointing vector to the nucleus frame
    for (ct=0;ct<nPixelsX*nPixelsY;ct++) {
      for(idim=0;idim<3;idim++) rPointing[idim][ct]=rotMat[idim][0]*rPointingInstrument[0][ct]+rotMat[idim][1]*rPointingInstrument[1][ct]+rotMat[idim][2]*rPointingInstrument[2][ct];
    }
  } //end of the time loop will have to be at the very end if we want to allow several times, remember to change iteration variable
  
  long int iTriangle;
  
  double origin[3]={0.0,0.0,0.0};  
  double lMax=1e6,rNucleus=0.0;
  double dl=0.0;
  bool intersection=false;
  for (i=iStart;i<iFinish;i++) {
    double scX[3]={xObservation[0],xObservation[1],xObservation[2]};
    rNucleus=sqrt(xObservation[0]*xObservation[0]+xObservation[1]*xObservation[1]+xObservation[2]*xObservation[2]);
    ct=0;
    intersection=true;
    //define the points of interest for the line of sight calculations
    while (rNucleus<lMax && intersection==true){
      intersection=false;
      dl=max(10.0,rNucleus/10.0); // /30
      for (idim=0;idim<3;idim++) scX[idim]+=rPointing[idim][i]*dl;
      rNucleus=sqrt(scX[0]*scX[0]+scX[1]*scX[1]+scX[2]*scX[2]);
      PIC::RayTracing::SetCutCellSCAttribute(scX,false);    
      
      //Check if the point is inside the nucleus (line of sight intersecting the nucleus)
      for (iTriangle=0;iTriangle<CutCell::nBoundaryTriangleFaces;iTriangle++){
	if(CutCell::BoundaryTriangleFaces[iTriangle].IntervalIntersection(scX,origin,0.0)==true){
	  intersection=true;  
	}
      }
      
      if(intersection==true) {
	for (j=0;j<CutCell::nBoundaryTriangleFaces;j++) {
	  CutCell::BoundaryTriangleFaces[j].GetCenterPosition(xCenter);
	  //Sun angle and SC angle
	  double sunAngle=0.0,scAngle=0.0;
	  double cSun=0.0,cSC=0.0,rCenter=0.0;
	  
	  for (idim=0;idim<3;idim++) {
	    cSun+=xLightSource[idim]*xCenter[idim];
	    cSC+=scX[idim]*xCenter[idim];
	    rCenter+=xCenter[idim]*xCenter[idim];	
	  }
	  
	  rCenter=sqrt(rCenter);
	  
	  scDistance=sqrt(pow(scX[0]-xCenter[0],2.0)+pow(scX[1]-xCenter[1],2.0)+pow(scX[2]-xCenter[2],2.0));
	  
	  sunAngle=acos(cSun/(sunDistance*rCenter));
	  scAngle=acos(cSC/(scDistance*rCenter));
	  
	  //compute C
	  if(scAngle<90*Pi/180.0) {
	    
	    sphericalHarmonic[0]=1; //Y00
	    sphericalHarmonic[1]=sin(colatitude[j])*sin(longitude[j]);//Y1-1
	    sphericalHarmonic[2]=cos(colatitude[j]); //Y10
	    sphericalHarmonic[3]=sin(colatitude[j])*cos(longitude[j]); //Y11
	    sphericalHarmonic[4]=pow(sin(colatitude[j]),2.0)*sin(2*longitude[j]);//Y2-2
	    sphericalHarmonic[5]=cos(colatitude[j])*sin(colatitude[j])*sin(longitude[j]);//Y2-1
	    sphericalHarmonic[6]=3*pow(cos(colatitude[j]),2.0)-1.0;//Y20
	    sphericalHarmonic[7]=cos(colatitude[j])*sin(colatitude[j])*cos(longitude[j]);//Y21
	    sphericalHarmonic[8]=pow(sin(colatitude[j]),2.0)*cos(2*longitude[j]);//Y22
	    sphericalHarmonic[9]=pow(sin(colatitude[j]),3.0)*sin(3*longitude[j]);//Y3-3                                        
	    sphericalHarmonic[10]=pow(sin(colatitude[j]),2.0)*cos(colatitude[j])*sin(2*longitude[j]);//Y3-2                    
	    sphericalHarmonic[11]=(5*pow(cos(colatitude[j]),2.0)-1)*sin(colatitude[j])*sin(longitude[j]);//Y3-1                
	    sphericalHarmonic[12]=5*pow(cos(colatitude[j]),3.0)-3*cos(colatitude[j]);//Y30                                     
	    sphericalHarmonic[13]=(5*pow(cos(colatitude[j]),2.0)-1)*sin(colatitude[j])*cos(longitude[j]);//Y31                 
	    sphericalHarmonic[14]=pow(sin(colatitude[j]),2.0)*cos(colatitude[j])*cos(2*longitude[j]);//Y32                     
	    sphericalHarmonic[15]=pow(sin(colatitude[j]),3.0)*cos(3*longitude[j]);//Y33              
	    sphericalHarmonic[16]=pow(sin(colatitude[j]),4.0)*sin(4*longitude[j]);//Y4-4              
	    sphericalHarmonic[17]=pow(sin(colatitude[j]),3.0)*cos(colatitude[j])*sin(3*longitude[j]);//Y4-3              
	    sphericalHarmonic[18]=pow(sin(colatitude[j]),2.0)*(7*pow(cos(colatitude[j]),2.0)-1)*sin(2*longitude[j]);//Y4-2              
	    sphericalHarmonic[19]=sin(colatitude[j])*(7*pow(cos(colatitude[j]),3.0)-3*cos(colatitude[j]))*sin(longitude[j]);//Y4-1              
	    sphericalHarmonic[20]=35*pow(cos(colatitude[j]),4.0)-30*pow(cos(colatitude[j]),2.0)+3;//Y40              
	    sphericalHarmonic[21]=sin(colatitude[j])*(7*pow(cos(colatitude[j]),3.0)-3*cos(colatitude[j]))*cos(longitude[j]);//Y41              
	    sphericalHarmonic[22]=pow(sin(colatitude[j]),2.0)*(7*pow(cos(colatitude[j]),2.0)-1)*cos(2*longitude[j]);//Y42              
	    sphericalHarmonic[23]=pow(sin(colatitude[j]),3.0)*cos(colatitude[j])*cos(3*longitude[j]);//Y43              
	    sphericalHarmonic[24]=pow(sin(colatitude[j]),4.0)*cos(4*longitude[j]);//Y44           
	    /*      sphericalHarmonic[25]=pow(sin(colatitude[j]),5.0)*sin(5*longitude[j]);//Y5-5              
		    sphericalHarmonic[26]=pow(sin(colatitude[j]),4.0)*cos(colatitude[j])*sin(4*longitude[j]);//Y5-4              
		    sphericalHarmonic[27]=pow(sin(colatitude[j]),3.0)*(9.0*pow(cos(colatitude[j]),2.0)-1.0)*sin(3*longitude[j]);//Y5-3
		    sphericalHarmonic[28]=pow(sin(colatitude[j]),2.0)*(3.0*pow(cos(colatitude[j]),3.0)-cos(colatitude[j]))*sin(2*longitude[j]);//Y5-2   
		    sphericalHarmonic[29]=sin(colatitude[j])*(21.0*pow(cos(colatitude[j]),4.0)-14.0*pow(cos(colatitude[j]),2.0)+1.0)*sin(longitude[j]);//Y5-1       
		    sphericalHarmonic[30]=65*pow(cos(colatitude[j]),5.0)-70*pow(cos(colatitude[j]),3.0)+15*cos(colatitude[j]);//Y50  
		    sphericalHarmonic[31]=sin(colatitude[j])*(21.0*pow(cos(colatitude[j]),4.0)-14.0*pow(cos(colatitude[j]),2.0)+1.0)*cos(longitude[j]);//Y51    
		    sphericalHarmonic[32]=pow(sin(colatitude[j]),2.0)*(3.0*pow(cos(colatitude[j]),3.0)-cos(colatitude[j]))*cos(2*longitude[j]);//Y52   
		    sphericalHarmonic[33]=pow(sin(colatitude[j]),3.0)*(9.0*pow(cos(colatitude[j]),2.0)-1.0)*cos(3*longitude[j]);//Y53
		    sphericalHarmonic[34]=pow(sin(colatitude[j]),4.0)*cos(colatitude[j])*cos(4*longitude[j]);//Y54              
		    sphericalHarmonic[35]=pow(sin(colatitude[j]),5.0)*cos(5*longitude[j]);//Y55              
		    sphericalHarmonic[36]=pow(sin(colatitude[j]),6.0)*sin(6*longitude[j]);//Y6-6              
		    sphericalHarmonic[37]=pow(sin(colatitude[j]),5.0)*cos(colatitude[j])*sin(5*longitude[j]);//Y6-5              
		    sphericalHarmonic[38]=pow(sin(colatitude[j]),4.0)*(11.0*pow(cos(colatitude[j]),2.0)-1.0)*sin(4*longitude[j]);//Y6-4
		    sphericalHarmonic[39]=pow(sin(colatitude[j]),3.0)*(11.0*pow(cos(colatitude[j]),3.0)-3.0*cos(colatitude[j]))*sin(3*longitude[j]);//Y6-3   
		    sphericalHarmonic[40]=pow(sin(colatitude[j]),2.0)*(33.0*pow(cos(colatitude[j]),4.0)-18.0*pow(cos(colatitude[j]),2.0)+1.0)*sin(2*longitude[j]);//6-2       
		    sphericalHarmonic[41]=sin(colatitude[j])*(33*pow(cos(colatitude[j]),5.0)-30*pow(cos(colatitude[j]),3.0)+5*cos(colatitude[j]))*sin(longitude[j]);//Y6-1  
		    sphericalHarmonic[42]=231*pow(cos(colatitude[j]),6.0)-315*pow(cos(colatitude[j]),4.0)+105*pow(cos(colatitude[j]),2.0)-5;//Y60 
		    sphericalHarmonic[43]=sin(colatitude[j])*(33*pow(cos(colatitude[j]),5.0)-30*pow(cos(colatitude[j]),3.0)+5*cos(colatitude[j]))*cos(longitude[j]);//Y61  
		    sphericalHarmonic[44]=pow(sin(colatitude[j]),2.0)*(33.0*pow(cos(colatitude[j]),4.0)-18.0*pow(cos(colatitude[j]),2.0)+1.0)*cos(2*longitude[j]);//62       
		    sphericalHarmonic[45]=pow(sin(colatitude[j]),3.0)*(11.0*pow(cos(colatitude[j]),3.0)-3.0*cos(colatitude[j]))*cos(3*longitude[j]);//Y63   
		    sphericalHarmonic[46]=pow(sin(colatitude[j]),4.0)*(11.0*pow(cos(colatitude[j]),2.0)-1.0)*cos(4*longitude[j]);//Y64
		    sphericalHarmonic[47]=pow(sin(colatitude[j]),5.0)*cos(colatitude[j])*cos(5*longitude[j]);//Y65              
		    sphericalHarmonic[48]=pow(sin(colatitude[j]),6.0)*cos(6*longitude[j]);//Y66    */          
	    
	    
	    //test if the nucleus blocks the view from the SC
	    if(CutCell::BoundaryTriangleFaces[j].pic__schidden_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__FALSE_) {
	      if(sunAngle>90*Pi/180.0) {
		g=a;
	      }else {
		if(CutCell::BoundaryTriangleFaces[j].pic__shadow_attribute==_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_) {
		  g=a;
		}else{
		  g=max(a,cos(sunAngle));
		}
	      }
	      
	      
	      for(harm=0;harm<nDev;harm++) Cvec[i*nDev+harm]+=g*CutCell::BoundaryTriangleFaces[j].SurfaceArea*cos(scAngle)/pow(scDistance,2.0)*sphericalHarmonic[harm]/pow(rAU,1.0)*dl; 
	    }
	  }
	  
	}
      }
    }
    if (PIC::ThisThread==0) printf("%i \n",i);
  }
  if (PIC::ThisThread==0) {
    for (i=0;i<nPixelsX*nPixelsY;i++) for (j=0;j<nDev;j++) C[i][j]+=Cvec[i*nDev+j];
    for (thread=1;thread<PIC::nTotalThreads;thread++) {
      MPI_Recv(buffer, nPixelsX*nPixelsY*nDev, MPI_DOUBLE, thread, 0, MPI_GLOBAL_COMMUNICATOR,&status);
      for (i=0;i<nPixelsX*nPixelsY;i++) for (j=0;j<nDev;j++) C[i][j]+=buffer[i*nDev+j];
    }
  }else if (PIC::ThisThread !=0) {
    for (i=0;i<nPixelsX*nPixelsY;i++) for (j=0;j<nDev;j++) buffer[i*nDev+j]+=Cvec[i*nDev+j];
    MPI_Send(buffer, nPixelsX*nPixelsY*nDev, MPI_DOUBLE, 0, 0, MPI_GLOBAL_COMMUNICATOR);
  } 
  
  if (PIC::ThisThread==0) {
    FILE *fC;
    fC = fopen("Cmatrix.dat","w"); 
    for (i=0;i<nPixelsX*nPixelsY;i++) {
      for (j=0;j<nDev;j++) {
	fprintf(fC,"%e ",C[i][j]);
      } 
      fprintf(fC,"\n");
    }
    fclose(fC);  
  }
  
  fclose(fout);  
  
  MPI_Finalize();
  
  return 1;
}

