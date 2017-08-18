//the function for sampling of the particle's velocity/evergy distribution function according to origin/population within the heliosphere

#include "pic.h"
#include "OH.h"

const int OH::Sampling::DistributionFunctionSample::_LINEAR_SAMPLING_SCALE_=0,OH::Sampling::DistributionFunctionSample::_LOGARITHMIC_SAMPLING_SCALE_=1;
const bool OH::Sampling::DistributionFunctionSample::Use = false;
int OH::Sampling::DistributionFunctionSample::v2SamplingMode=_LINEAR_SAMPLING_SCALE_,OH::Sampling::DistributionFunctionSample::speedSamplingMode=_LINEAR_SAMPLING_SCALE_;
double OH::Sampling::DistributionFunctionSample::vMin=-1000.0;
double OH::Sampling::DistributionFunctionSample::vMax=1000.0;
long int OH::Sampling::DistributionFunctionSample::nSampledFunctionPoints=100;
double** OH::Sampling::DistributionFunctionSample::SamplingBuffer=NULL;
double OH::Sampling::DistributionFunctionSample::SamplingLocations[][3]={{0.0,0.0,0.0}};
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>** OH::Sampling::DistributionFunctionSample::SampleNodes=NULL;
double OH::Sampling::DistributionFunctionSample::dV=0.0,OH::Sampling::DistributionFunctionSample::dV2=0.0,OH::Sampling::DistributionFunctionSample::dSpeed=0.0;
long int *OH::Sampling::DistributionFunctionSample::SampleLocalCellNumber=NULL;
int OH::Sampling::DistributionFunctionSample::nSampleLocations=0;
bool OH::Sampling::DistributionFunctionSample::SamplingInitializedFlag=false;

int OH::Sampling::DistributionFunctionSample::Sample_Velocity_Offset=0,OH::Sampling::DistributionFunctionSample::Sample_Speed_Offset=0;
int OH::Sampling::DistributionFunctionSample::Sample_V2_Offset=0,OH::Sampling::DistributionFunctionSample::SampleDataLength=0;

//2d phase space distribution function variables
OH::Sampling::DistributionFunctionSample::cSampled2DFunction **OH::Sampling::DistributionFunctionSample::Sampled2DFunction=NULL;
int OH::Sampling::DistributionFunctionSample::nSampledPlanes=DIM;

//====================================================
//init the sampling buffers
void OH::Sampling::DistributionFunctionSample::Init() {
  int nProbe,i,j,k;

  if (SamplingInitializedFlag==true) exit(__LINE__,__FILE__,"Error: DistributionFunctionSample is already initialized");

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_OFF_
  if (PIC::Mesh::mesh.ThisThread==0) fprintf(PIC::DiagnospticMessageStream,"WARNING: Sampling of the distribution function is prohibited in the settings of the model");
  return;
#endif


//  nSampleLocations=nProbeLocations;
  SamplingInitializedFlag=true;

  //get the lengths of the sampling intervals
  double t0,t1; //temporary variables to satisfy intel c++ compiler

  dV=1.05*(vMax-vMin)/(nSampledFunctionPoints-1);
  dSpeed=((speedSamplingMode==_LINEAR_SAMPLING_SCALE_) ? max((t0=fabs(vMax)),(t1=fabs(vMin)))*sqrt(2.0) : log(max((t0=fabs(vMax)),(t1=fabs(vMin)))*sqrt(2.0)))/(nSampledFunctionPoints-1);
  dV2=((v2SamplingMode==_LINEAR_SAMPLING_SCALE_) ? max(t0=vMax*vMax,t1=vMin*vMin)*sqrt(2.0) : log(max(t0=vMax*vMax,t1=vMin*vMin)*sqrt(2.0)))/(nSampledFunctionPoints-1);


  Sample_Velocity_Offset=0;
  SampleDataLength=DIM;

  Sample_Speed_Offset=SampleDataLength;
  SampleDataLength++;

  Sample_V2_Offset=SampleDataLength;
  SampleDataLength++;

  //allocate the sampling buffers
  SampleLocalCellNumber=new long int [nSampleLocations];
  SampleNodes=new cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* [nSampleLocations];
//  SamplingLocations=new double* [nProbeLocations];
//  SamplingLocations[0]=new double [DIM*nProbeLocations];

  // OriginTag+1, the 1 is the place in the buffer where we calculate the total solution across all origin tags
  SamplingBuffer=new double* [nSampleLocations];
  SamplingBuffer[0]=new double [nSampleLocations*(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::nPhysSpec*SampleDataLength*(nSampledFunctionPoints-1)];

  for (nProbe=1;nProbe<nSampleLocations;nProbe++) {
//    SamplingLocations[nProbe]=SamplingLocations[nProbe-1]+DIM;
    SamplingBuffer[nProbe]=SamplingBuffer[nProbe-1]+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::nPhysSpec*SampleDataLength*(nSampledFunctionPoints-1);
  }

  // initialization of the array of objects used to caclulate the 2d distribution functions- OriginTag+1, the 1 is the place in the buffer where we calculate the total solution across all origin tags
  Sampled2DFunction=new OH::Sampling::DistributionFunctionSample::cSampled2DFunction* [nSampleLocations];
  Sampled2DFunction[0]=new OH::Sampling::DistributionFunctionSample::cSampled2DFunction [nSampleLocations*(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::nPhysSpec*nSampledPlanes];

  for (nProbe=1;nProbe<nSampleLocations;nProbe++) {
    Sampled2DFunction[nProbe]=Sampled2DFunction[nProbe-1]+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::nPhysSpec*nSampledPlanes;
  }

  //initializing each object with correct direction of vectors describing plane
  int Sample_2D_Offset;
  char fname[_MAX_STRING_LENGTH_PIC_];
  char zoneName[_MAX_STRING_LENGTH_PIC_];
  char vx[_MAX_STRING_LENGTH_PIC_];
  char vy[_MAX_STRING_LENGTH_PIC_];
  char vz[_MAX_STRING_LENGTH_PIC_];
  char vr[_MAX_STRING_LENGTH_PIC_];
  char vtheta[_MAX_STRING_LENGTH_PIC_];
  char vphi[_MAX_STRING_LENGTH_PIC_];
  
  // cartesian coordinates
  double DirX[3]={1.0,0.0,0.0},DirY[3]={0.0,1.0,0.0},DirZ[3]={0.0,0.0,1.0};
  sprintf(vx,"Vx"); sprintf(vy,"Vy"); sprintf(vz,"Vz");
  
  //spherical coordinate vectors
  double DirR[3]={1.0,0.0,0.0},DirTheta[3]={0.0,1.0,0.0},DirPhi[3]={0.0,0.0,1.0};
  sprintf(vr,"Vr"); sprintf(vtheta,"Vtheta"); sprintf(vphi,"Vphi");
  double x=0.0,y=0.0,z=0.0,r=0.0,xy=0.0;
  
  for (nProbe=0;nProbe<nSampleLocations;nProbe++) {
    x=SamplingLocations[nProbe][0];
    if(DIM==3){ 
      y=SamplingLocations[nProbe][1]; 
      z=SamplingLocations[nProbe][2];
    }
    r=sqrt(x*x+y*y+z*z);xy=sqrt(x*x+y*y);
    //setting directions to calculate the velocities in spherical coords
    DirR[0]=x/r; DirR[1]=y/r; DirR[2]=x/xy;
    DirTheta[0]=-y/xy; DirTheta[1]=x/xy; DirTheta[2]=0;
    DirPhi[0]=z*x/(r*xy); DirPhi[0]=z*y/(r*xy); DirPhi[0]=-xy/r;
    
    for (int physpec=0;physpec<OH::nPhysSpec; physpec++){
      for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	Sample_2D_Offset=iSource*nSampledPlanes+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*physpec*nSampledPlanes;
	sprintf(fname,"PhysSpec=%i.2DVelocityDistFunc.nSampledPoint=%i",physpec,nProbe);
	if (iSource==0){
	  //rotating cartesian coordinates into spherical coordinates more useful for solar wind
	  sprintf(zoneName,"R-Theta plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][0+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirR,DirTheta,vr,vtheta,fname,zoneName);
	  sprintf(zoneName,"R-Phi plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][1+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirR,DirPhi,vr,vphi,fname,zoneName);
	  sprintf(zoneName,"Theta-Phi plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][2+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirTheta,DirPhi,vtheta,vphi,fname,zoneName);
	}
	else if(iSource == OH::Sampling::OriginLocation::nSampledOriginLocations){
	  // use default cartesian coordinates
	  sprintf(zoneName,"XY plane (Total Solution)");
	  Sampled2DFunction[nProbe][0+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirX,DirY,vx,vy,fname,zoneName);
	  sprintf(zoneName,"XZ plane (Total Solution)");
	  Sampled2DFunction[nProbe][1+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirX,DirZ,vx,vz,fname,zoneName);
	  sprintf(zoneName,"YZ plane (Total Solution)");
	  Sampled2DFunction[nProbe][2+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirY,DirZ,vy,vz,fname,zoneName);
	}
	else{
	  // use default cartesian coordinates
	  sprintf(zoneName,"XY plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][0+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirX,DirY,vx,vy,fname,zoneName);
	  sprintf(zoneName,"XZ plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][1+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirX,DirZ,vx,vz,fname,zoneName);
	  sprintf(zoneName,"YZ plane (OriginID=%i)",iSource);
	  Sampled2DFunction[nProbe][2+Sample_2D_Offset].Init(nSampledFunctionPoints,nSampledFunctionPoints,dV,dV,vMin,vMin,DirY,DirZ,vy,vz,fname,zoneName);
	}					   						    
      }
    }
  }

  //init the sampling informations
  for (nProbe=0;nProbe<nSampleLocations;nProbe++) {
//    for (idim=0;idim<DIM;idim++) SamplingLocations[nProbe][idim]=ProbeLocations[nProbe][idim];

    SampleNodes[nProbe]=PIC::Mesh::mesh.findTreeNode(SamplingLocations[nProbe]);
    if (SampleNodes[nProbe]==NULL) exit(__LINE__,__FILE__,"Error: the point is outside of the domain");

    SampleLocalCellNumber[nProbe]=PIC::Mesh::mesh.fingCellIndex(SamplingLocations[nProbe],i,j,k,SampleNodes[nProbe],false);
    if (SampleLocalCellNumber[nProbe]==-1) exit(__LINE__,__FILE__,"Error: cannot find the cell");
  }

  flushSamplingBuffers();
}

//====================================================
//flush the sampling buffer
void OH::Sampling::DistributionFunctionSample::flushSamplingBuffers() {
  long int i,TotalDataLength=nSampleLocations*(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::nPhysSpec*SampleDataLength*(nSampledFunctionPoints-1);
  double *ptr=SamplingBuffer[0];

  for (i=0;i<TotalDataLength;i++,ptr++) *ptr=0.0;
}
//====================================================
//return the offset where the sample data for the particular species, sampling interval and the sampling point are located
long int OH::Sampling::DistributionFunctionSample::GetSampleDataOffset(int spec,int OriginID,int SampleVariableOffset) {
  long int offset;

  offset=(OriginID+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::PhysSpec[spec])*SampleDataLength*(nSampledFunctionPoints-1);
  offset+=SampleVariableOffset*(nSampledFunctionPoints-1);
  
  return offset;
}
//====================================================
//Sample the distribution function
void OH::Sampling::DistributionFunctionSample::SampleDistributionFunction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,idim,offset,offsetTotal;
  double LocalParticleWeight,v2;

/*
    for (node=NULL,nProbe=0;nProbe<nSampleLocations;nProbe++) if (sampleNode==SampleNodes[nProbe]) {
      node=sampleNode;
      break;
    }

    if (node!=NULL) if (node->block!=NULL) {
    */

  for (node=SampleNodes[0],nProbe=0;nProbe<nSampleLocations;node=SampleNodes[++nProbe]) if (node->Thread==PIC::ThisThread) {
      double *v;
      PIC::ParticleBuffer::byte *ParticleData;
      int i,j,k, OriginID;

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
      ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

//      if (PIC::Mesh::mesh.fingCellIndex(SampleLocalCellNumber[nProbe],i,j,k,node,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");
//      ptr=node->block->GetCenterNode(SampleLocalCellNumber[nProbe])->FirstCellParticle;

      // loop through all particles in the cell
      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        v=PIC::ParticleBuffer::GetV(ParticleData);
	OriginID=OH::GetOriginTag(ParticleData);

        LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
        LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        for (v2=0.0,idim=0;idim<DIM;idim++) {
          i=(int)((v[idim]-vMin)/dV);
          v2+=pow(v[idim],2);
          offset=GetSampleDataOffset(spec,OriginID,Sample_Velocity_Offset+idim);
	  
	  // offset to calculate the total solution
	  offsetTotal=GetSampleDataOffset(spec,OH::Sampling::OriginLocation::nSampledOriginLocations,Sample_Velocity_Offset+idim);

          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight, SamplingBuffer[nProbe][offsetTotal+i]+=LocalParticleWeight; 
        }

        if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
          i=(int)(sqrt(v2)/dSpeed);
          offset=GetSampleDataOffset(spec,OriginID,Sample_Speed_Offset);
	  // offset to calculate the total solution
	  offsetTotal=GetSampleDataOffset(spec,OH::Sampling::OriginLocation::nSampledOriginLocations,Sample_Speed_Offset);
          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight, SamplingBuffer[nProbe][offsetTotal+i]+=LocalParticleWeight;

          i=(int)(v2/dV2);
	  // offset to calculate the total solution
	  offsetTotal=GetSampleDataOffset(spec,OH::Sampling::OriginLocation::nSampledOriginLocations,Sample_V2_Offset);
          offset=GetSampleDataOffset(spec,OriginID,Sample_V2_Offset);
          if ((i>=0)&&(i<nSampledFunctionPoints-1)) SamplingBuffer[nProbe][offset+i]+=LocalParticleWeight, SamplingBuffer[nProbe][offsetTotal+i]+=LocalParticleWeight;
        }
        else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
          exit(__LINE__,__FILE__,"not implemented");
        }
        else exit(__LINE__,__FILE__,"Error: the option is unknown");

        ptr=PIC::ParticleBuffer::GetNext(ParticleData);
      }
    }

}


//====================================================
//print the distribution function into a file
void OH::Sampling::DistributionFunctionSample::printDistributionFunction(int DataOutputFileNumber) {
  long int idim,nProbe,i,nVariable,thread,offset;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0,dInterval=0.0;
  char str[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);

  // loop over all  species
  for (int spec=0;spec<PIC::nTotalSpecies;spec++){
    // print out 1 file per physical species
    if (OH::DoPrintSpec[spec]){
      PIC::MolecularData::GetChemSymbol(ChemSymbol,spec);
      // loop over number of points (nSampledLocation) want distribution function at
      for (nProbe=0;nProbe<nSampleLocations;nProbe++) {
	if (PIC::Mesh::mesh.ThisThread==0) {
	  sprintf(str,"%s/pic.%s.s=%i.VelocityDistributionFunction.nSamplePoint=%ld.out=%i.dat",PIC::OutputDataFileDirectory,ChemSymbol,spec,nProbe,DataOutputFileNumber);
	  fout=fopen(str,"w");
	  
	  fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);
	  
	  fprintf(fout,"\"TITLE=Distribution function at x=%e",SamplingLocations[nProbe][0]);
	  for (idim=1;idim<DIM;idim++) fprintf(fout,", %E",SamplingLocations[nProbe][idim]);
	  fprintf(fout,", sampled over velocity range %E to %E [m/s]",vMin,vMax);	  

	  fprintf(fout,"\"\nVARIABLES=\"V\",\"|V|\",\"V2\"");
	  for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations;iSource++) {
	    fprintf(fout,", \"f(Vx) (Source Region ID=%i)\"",iSource);
	    if (DIM>1) fprintf(fout,",\"f(Vy) (Source Region ID=%i)\"",iSource);
	    if (DIM>2) fprintf(fout,",\"f(Vz) (Source Region ID=%i)\"",iSource);
	    
	    fprintf(fout,",\"f(|V|) (Source Region ID=%i)\"",iSource);
	    fprintf(fout,", \"f(V2) (Source Region ID=%i)\"",iSource);
	  }
    
	  fprintf(fout,", \"f(Vx) (Total Solution)\"");
	  if (DIM>1) fprintf(fout,",\"f(Vy) (Total Solution)\"");
	  if (DIM>2) fprintf(fout,",\"f(Vz) (Total Solution)\"");
	  fprintf(fout,",\"f(|V|) (Total Solution)\",\"f(V2) (Total Solution)\"\n");
	
	  //collect the sampled information from other processors
	  for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++){
	    for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	      for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
		offset=GetSampleDataOffset(spec,iSource,nVariable);
	      
		for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]+=pipe.recv<double>(thread);
	      }
	    }
	  }

	  //normalize the distribution functions
	  for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	    for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
	      norm=0.0;
	      offset=GetSampleDataOffset(spec,iSource,nVariable);
	      
	      for (i=0;i<nSampledFunctionPoints-1;i++) {
		if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
		  if (nVariable<DIM) dInterval=dV;
		  else if (nVariable==Sample_Speed_Offset) dInterval=dSpeed;
		  else if (nVariable==Sample_V2_Offset) dInterval=dV2;
		  else exit(__LINE__,__FILE__,"Error: unknown option");
		}
		else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
		  exit(__LINE__,__FILE__,"not implemented");
		}
		else exit(__LINE__,__FILE__,"Error: the option is unknown");
		
		norm+=SamplingBuffer[nProbe][i+offset]*dInterval;
	      }
	      
	      if (fabs(norm)>0.0) for (i=0;i<nSampledFunctionPoints-1;i++) SamplingBuffer[nProbe][i+offset]/=norm;
	    }
	  }
	  
	  //print the output file
	  for (i=0;i<nSampledFunctionPoints-1;i++) {
	    double v=0.0,v2=0.0,Speed=0.0;
	      
	    if (speedSamplingMode==_LINEAR_SAMPLING_SCALE_) {
	      v=vMin+i*dV,v2=i*dV2,Speed=i*dSpeed;
	    }
	    else if (speedSamplingMode==_LOGARITHMIC_SAMPLING_SCALE_) {
	      exit(__LINE__,__FILE__,"not implemented");
	    }
	    else exit(__LINE__,__FILE__,"Error: the option is unknown");
	    
	    fprintf(fout,"%e  %e  %e",v, Speed, v2);
	    
	    for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	      for (idim=0;idim<DIM;idim++) {
		offset=GetSampleDataOffset(spec,iSource,idim);
		fprintf(fout,"  %e",SamplingBuffer[nProbe][i+offset]);
	      }
	      
	      offset=GetSampleDataOffset(spec,iSource,Sample_Speed_Offset);
	      fprintf(fout,"  %e",SamplingBuffer[nProbe][i+offset]);
	      
	      offset=GetSampleDataOffset(spec,iSource,Sample_V2_Offset);
	      fprintf(fout,"  %e",SamplingBuffer[nProbe][i+offset]);	      
	    }
	    fprintf(fout,"\n"); // returns line
	  }
	  
	  //close the output file
	  fclose(fout);
	  fprintf(PIC::DiagnospticMessageStream,"done.\n");
	}
	else {
	  for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	    for (nVariable=0;nVariable<SampleDataLength;nVariable++) {
	      offset=GetSampleDataOffset(spec,iSource,nVariable);
	      
	      for (i=0;i<nSampledFunctionPoints-1;i++) {
		pipe.send(SamplingBuffer[nProbe][i+offset]);
		SamplingBuffer[nProbe][i+offset]=0.0; //this sampled information is stored by the root processor
	      }
	    }
	  }
	}
      }
    }
  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

//====================================================
// the print function for the calss used to sample 2D distribution functions
void OH::Sampling::DistributionFunctionSample::cSampled2DFunction::Print(int DataOutputFileNumber, char* printCharIN){
  long int idim,thread;
  FILE *fout=NULL;
  CMPI_channel pipe(1000000);
  double norm=0.0;
  char str[_MAX_STRING_LENGTH_PIC_];

  sprintf(printChar, "%s", printCharIN);

  if (PIC::Mesh::mesh.ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);


  if (PIC::Mesh::mesh.ThisThread==0) {
    sprintf(str,"%s/pic.%s.out=%i.dat",PIC::OutputDataFileDirectory,
	    this->fname,DataOutputFileNumber);
    fout=fopen(str,printChar);
	  
    fprintf(PIC::DiagnospticMessageStream,"printing output file: %s.........         ",str);
	  
    //fprintf(fout,"\"TITLE=2D distribution function at x=%e",SamplingLocations[nProbe][0]);
    //for (idim=1;idim<DIM;idim++) fprintf(fout,", %E",SamplingLocations[nProbe][idim]);
    //fprintf(fout,", sampled over velocity range %E to %E [m/s]",vMin,vMax);	  

    fprintf(fout,"\nVARIABLES=\"%s\",\"%s\", \"f(%s,%s)\"\n", this->Dir1name, this->Dir2name,this->Dir1name, this->Dir2name);
    fprintf(fout,"ZONE T=\"%s\" I=%i, J=%i, DATAPACKING=POINT\n",this->zoneName,this->N1-1,this->N2-1);
    
    //collect the sampled information from other processors
    for (thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++)
      for (int i1=0;i1<this->N1-1;i1++)      
	for (int i2=0;i2<this->N2-1;i2++) 
	  this->Buffer[i1][i2]+=pipe.recv<double>(thread);
    
    //normalize the distribution functions
    norm=0.0;
    for (int i1=0;i1<this->N1-1;i1++)      
      for (int i2=0;i2<this->N2-1;i2++) 
	norm += this->Buffer[i1][i2] * this->dV1 * this->dV2;
    
    if (fabs(norm)>0.0) 
      for (int i1=0;i1<this->N1-1;i1++)      
	for (int i2=0;i2<this->N2-1;i2++) 
	  this->Buffer[i1][i2]/=norm;
	  
    //print the output file
    for (int i1=0;i1<this->N1-1;i1++)      
      for (int i2=0;i2<this->N2-1;i2++){
	double v1=0.0,v2=0.0;
	v1 = this->VMin1 + i1*dV1;
	v2 = this->VMin2 + i2*dV2;
	
	fprintf(fout,"%e  %e  %e", v1, v2, this->Buffer[i1][i2]);
	fprintf(fout,"\n"); // returns line
      }
    //close the output file
    fclose(fout);
    fprintf(PIC::DiagnospticMessageStream,"done.\n");
  }
  else {
	
    for (int i1=0;i1<this->N1-1;i1++)      
      for (int i2=0;i2<this->N2-1;i2++){
	pipe.send(this->Buffer[i1][i2]);
	this->Buffer[i1][i2]=0.0; //this sampled information is stored by the root processor
    }

  }

  if (PIC::Mesh::mesh.ThisThread==0) pipe.closeRecvAll();
  else pipe.closeSend();

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

}

//====================================================
//Sample the 2D distribution function
void OH::Sampling::DistributionFunctionSample::Sample2dDistributionFunction() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node;
  long int ptr,nProbe,spec,idim,offset,offsetTotal;
  double LocalParticleWeight;

  // only do this for 3D
  if (DIM != 3) return;

  for (node=SampleNodes[0],nProbe=0;nProbe<nSampleLocations;node=SampleNodes[++nProbe]) if (node->Thread==PIC::ThisThread) {
      double *v;
      PIC::ParticleBuffer::byte *ParticleData;
      int i,j,k, OriginID;

      PIC::Mesh::mesh.convertCenterNodeLocalNumber2LocalCoordinates(SampleLocalCellNumber[nProbe],i,j,k);
      ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      // loop through all particles in the cell
      while (ptr!=-1) {
        ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        spec=PIC::ParticleBuffer::GetI(ParticleData);
        v=PIC::ParticleBuffer::GetV(ParticleData);
	OriginID=OH::GetOriginTag(ParticleData);

        LocalParticleWeight=node->block->GetLocalParticleWeight(spec);
        LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
	
	// add particle to each plane of the corresponding origin tag
	offset=OriginID*nSampledPlanes+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::PhysSpec[spec]*nSampledPlanes;
	Sampled2DFunction[nProbe][0+offset].Sample(v,LocalParticleWeight);
	Sampled2DFunction[nProbe][1+offset].Sample(v,LocalParticleWeight);
	Sampled2DFunction[nProbe][2+offset].Sample(v,LocalParticleWeight);

	// add particle to each plane of the total solution for that physical species
	offsetTotal=OH::Sampling::OriginLocation::nSampledOriginLocations*nSampledPlanes+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*OH::PhysSpec[spec]*nSampledPlanes;
	Sampled2DFunction[nProbe][0+offsetTotal].Sample(v,LocalParticleWeight);
	Sampled2DFunction[nProbe][1+offsetTotal].Sample(v,LocalParticleWeight);
	Sampled2DFunction[nProbe][2+offsetTotal].Sample(v,LocalParticleWeight);

	ptr=PIC::ParticleBuffer::GetNext(ParticleData);
      }
    }
}

//====================================================
//print the distribution function into a file
void OH::Sampling::DistributionFunctionSample::print2dDistributionFunction(int DataOutputFileNumber) {
  long int offset;

  if (DIM != 3) return;

  for (int nProbe=0;nProbe<nSampleLocations;nProbe++) {
    for (int physpec=0;physpec<OH::nPhysSpec; physpec++){
      for (int iSource=0;iSource<OH::Sampling::OriginLocation::nSampledOriginLocations+1;iSource++) {
	offset=iSource*nSampledPlanes+(OH::Sampling::OriginLocation::nSampledOriginLocations+1)*physpec*nSampledPlanes;
	
	// printing out the three planes per tag per physical species
	if (iSource ==0){
	  Sampled2DFunction[nProbe][0+offset].Print(DataOutputFileNumber,"w");
	}
	else Sampled2DFunction[nProbe][0+offset].Print(DataOutputFileNumber,"a");
	Sampled2DFunction[nProbe][1+offset].Print(DataOutputFileNumber,"a");
	Sampled2DFunction[nProbe][2+offset].Print(DataOutputFileNumber,"a");
      }
    }
  }

}
