//$Id$
//routines for operating with the particle trajectory file

/*
 * PostProcess3D_ParticleTrajectory.cpp
 *
 *  Created on: Jan 12, 2016
 *      Author: vtenishe
 */

#include "PostProcess3D.h"
#include "ifileopr.h"

int cPostProcess3D::cParticleTrajectory::nTrajectoryVariables=0;

//allocate individual trajectory data buffer
void cPostProcess3D::cParticleTrajectory::cIndividualTrajectoryData::AllocateDataBuffer(int n) {
  nDataPoints=n;

  if (nDataPoints!=0) {
    Data=new double* [nDataPoints];
    Data[0]=new double [nDataPoints*cPostProcess3D::cParticleTrajectory::nTrajectoryVariables];
    for (int i=1;i<nDataPoints;i++) Data[i]=Data[i-1]+cPostProcess3D::cParticleTrajectory::nTrajectoryVariables;
  }
}

//=============================================================
//add trajectory data to the list of the trajectories
void cPostProcess3D::cParticleTrajectory::AddIndividualTrajectoryData(int& nDataPoints,std::vector<double>& data) {
  cIndividualTrajectoryData t;
  int i;

  t.AllocateDataBuffer(nDataPoints);
  for (i=0;i<nDataPoints*nTrajectoryVariables;i++) t.Data[0][i]=data[i];

  //add the trajectory data
  IndividualTrajectories.push_back(t);
  nTotalTrajectories++;

  //clean the tamporaty data
  nDataPoints=0;
  data.clear();
}

//=============================================================
//load the trajectory data file
void cPostProcess3D::cParticleTrajectory::LoadDataFile(const char *fname,const char* path) {
  std::vector<double> TrajectoryData;
  char FullName[50000],str[50000],str1[50000];
  CiFileOperations ifile;
  FILE *fBinaryIn=NULL,*fBinaryOut=NULL;
  int nDataPoints=0;
  int nLoadedTrajectories=0,nOriginalTrajectories=nTotalTrajectories;


  //the slave processers open the files first
  if (PostProcess3D->rank==0) MPI_Barrier(MPI_COMM_WORLD);

  //get the full name of the data file and open the file
  sprintf(FullName,"%s/%s",path,fname);

  if (access(FullName,R_OK)!=0) {
     printf("Cannot find the file:%s\n",FullName);
     exit(__LINE__,__FILE__);
  }

  ifile.openfile(FullName);

  //check whether the binary file that corresponds to the data file 'fname' exists. If the file exists -> compare the time
  char BinaryFullName[5000];

  sprintf(BinaryFullName,"%s.trajectory-post-processing.tmp.bin",fname);

  if (access(BinaryFullName,R_OK)==0) {
    //the binary file exists: comapre the creation time
    struct stat t_stat;
    struct tm *timeinfo;
    time_t BinaryFileCreationTime,DataFileCreationFile;

    stat(BinaryFullName,&t_stat);
    timeinfo = localtime(&t_stat.st_ctime);
    BinaryFileCreationTime=mktime(timeinfo);

    stat(FullName,&t_stat);
    timeinfo = localtime(&t_stat.st_ctime);
    DataFileCreationFile=mktime(timeinfo);

    if (difftime(BinaryFileCreationTime,DataFileCreationFile)>0.0) {
      //the binary files created later than the data file -> use the binary file
      fBinaryIn=fopen(BinaryFullName,"r");
    }
  }

  if (fBinaryIn==NULL) {
    //the binary file either do not exists or created after after the data file -> create a new binary file
    if (PostProcess3D->rank==0) {
      fBinaryOut=fopen(BinaryFullName,"w");
      system("rm -f post-processing.trajectory-cell-distribution.tmp.bin"); //remove the file that assignes particle trajectories to the cells
    }
  }

  //read the variable line
  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);

  if (VariableList.size()==0) {
    while (strcmp(str1,"")!=0) {
      std::string var(str1);
      VariableList.push_back(var);
      nTrajectoryVariables++;

      ifile.CutInputStr(str1,str);
    }
  }

  //all processors are sincronize at this point
  if (PostProcess3D->rank!=0) MPI_Barrier(MPI_COMM_WORLD);

  //read the trajectory information
  if (fBinaryIn==NULL) {
    bool StartFileReading=false;

    while (ifile.eof()==false) {
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(str1,str);

      if (strcmp(str1,"ZONE")==0) {
        //beginig of the new trajectory: close the previous trajectory and prepare the buffer for the new trajectory
        if (StartFileReading==true) {
          AddIndividualTrajectoryData(nDataPoints,TrajectoryData);
          nLoadedTrajectories++;
        }
        else StartFileReading=true;
      }
      else {
        if (nLoadedTrajectories%PostProcess3D->size==PostProcess3D->rank) {
          //the line is a trajectory point
          for (int i=0;i<nTrajectoryVariables;i++) {
            TrajectoryData.push_back(atof(str1));
            ifile.CutInputStr(str1,str);
          }

          nDataPoints++;
        }
      }
    }

    //save the last trajectory
    AddIndividualTrajectoryData(nDataPoints,TrajectoryData);
    nLoadedTrajectories++;

    //collect the trajectory data stored on different processors
    for (int n=0;n<nLoadedTrajectories;n++) {
      if (n%PostProcess3D->size==PostProcess3D->rank) {
        //send the trajectory data
        MPI_Bcast(&IndividualTrajectories[n+nOriginalTrajectories].nDataPoints,1,MPI_INT,PostProcess3D->rank,MPI_COMM_WORLD);
        MPI_Bcast(IndividualTrajectories[n+nOriginalTrajectories].Data[0],IndividualTrajectories[n+nOriginalTrajectories].nDataPoints*nTrajectoryVariables,MPI_DOUBLE,PostProcess3D->rank,MPI_COMM_WORLD);
      }
      else {
        //recieve the trajectory data
        MPI_Bcast(&nDataPoints,1,MPI_INT,n%PostProcess3D->size,MPI_COMM_WORLD);
        IndividualTrajectories[n+nOriginalTrajectories].AllocateDataBuffer(nDataPoints);
        MPI_Bcast(IndividualTrajectories[n+nOriginalTrajectories].Data[0],IndividualTrajectories[n+nOriginalTrajectories].nDataPoints*nTrajectoryVariables,MPI_DOUBLE,n%PostProcess3D->size,MPI_COMM_WORLD);
      }
    }

    //save the binary file with the trajectory data
    if (PostProcess3D->rank==0) {
      fwrite(&nLoadedTrajectories,sizeof(int),1,fBinaryOut);

      for (int nTrajectory=0;nTrajectory<nLoadedTrajectories;nTrajectory++) {
        int n=IndividualTrajectories[nTrajectory+nOriginalTrajectories].nDataPoints;
        double *d=IndividualTrajectories[nTrajectory+nOriginalTrajectories].Data[0];

        fwrite(&n,sizeof(int),1,fBinaryOut);
        fwrite(d,sizeof(double),n*nTrajectoryVariables,fBinaryOut);
      }

      fclose(fBinaryOut);
    }
  }
  else {
    //read the trajectory data from a binary file
    int nReadTrajectories;

    fread(&nReadTrajectories,sizeof(int),1,fBinaryIn);
    nTotalTrajectories+=nReadTrajectories;

    for (int nTrajectory=0;nTrajectory<nReadTrajectories;nTrajectory++) {
      cIndividualTrajectoryData t;
      int nDataPoints;

      fread(&nDataPoints,sizeof(int),1,fBinaryIn);
      t.AllocateDataBuffer(nDataPoints);
      fread(t.Data[0],sizeof(double),t.nDataPoints*nTrajectoryVariables,fBinaryIn);

      IndividualTrajectories.push_back(t);
    }
  }

  //close the binary files
  if (fBinaryIn!=NULL) fclose(fBinaryIn);
  if (fBinaryOut!=NULL) fclose(fBinaryOut);

  //close the trajectory data file and save the binary trajectory data
  ifile.closefile();
}


//=====================================================================================
//save trajectory data file
void cPostProcess3D::cParticleTrajectory::PrintDataFileHeader(const char* fname) {
  if (PostProcess3D->rank!=0) return;
  FILE *fout=fopen(fname,"w");

  fprintf(fout,"VARIABLES=\"%s\"",VariableList[0].c_str());
  for (int nvar=1;nvar<VariableList.size();nvar++) fprintf(fout,", \"%s\"",VariableList[nvar].c_str());

  fprintf(fout,"\n");
  fclose(fout);
}

void cPostProcess3D::cParticleTrajectory::AddTrajectoryDataFile(cIndividualTrajectoryData* Trajectory,int TrajectoryNumber,const char* fname) {
  if (PostProcess3D->rank!=0) return;
  FILE *fout=fopen(fname,"a");

  fprintf(fout,"\nZONE T=\"Trajectory=%i\" F=POINT\n",TrajectoryNumber);

  for (int nline=0;nline<Trajectory->nDataPoints;nline++) {
    for (int i=0;i<nTrajectoryVariables;i++) fprintf(fout," %e",Trajectory->Data[nline][i]);
    fprintf(fout,"\n");
  }

  fclose(fout);
}


//======================================================================================
//print the surface data
void cPostProcess3D::cParticleTrajectory::PrintSurfaceData(const char *fname) {

  if (PostProcess3D->rank!=0) return;

  FILE *fout=fopen(fname,"w");
  int i,nface;

  struct cNodelInterpolationBase {
    double TotalInterpolationWeight;
    int StencilLength;
    int StencilOffset;
  };

  struct cFaceInterpolationStencil {
    double Weight;
    double *SurfaceData;
  };

  //prepare the list of the faces that contains the same node
  cNodelInterpolationBase *InterpolationBase=new cNodelInterpolationBase[CutCell::nBoundaryTriangleNodes];
  cFaceInterpolationStencil *Stencil=new cFaceInterpolationStencil[3*CutCell::nBoundaryTriangleFaces];
  int maxStencilLength=0,offset=0;

  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) {
    CutCell::BoundaryTriangleNodes[i].id=i;
    InterpolationBase[i].TotalInterpolationWeight=0.0;
    InterpolationBase[i].StencilLength=0;
    InterpolationBase[i].StencilOffset=0;
  }

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) for (i=0;i<3;i++) {
    InterpolationBase[CutCell::BoundaryTriangleFaces[nface].node[i]->id].StencilLength++;
  }


  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) {
    InterpolationBase[i].StencilOffset=offset;

    offset+=InterpolationBase[i].StencilLength;
    if (maxStencilLength<InterpolationBase[i].StencilLength) maxStencilLength=InterpolationBase[i].StencilLength;
  }

  //populate the interpolation base list
  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) InterpolationBase[i].StencilLength=0;

  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) for (i=0;i<3;i++) {
    double weight=CutCell::BoundaryTriangleFaces[nface].SurfaceArea;
    int nnode=CutCell::BoundaryTriangleFaces[nface].node[i]->id;

//    Stencil[InterpolationBase[nnode].StencilOffset+InterpolationBase[nnode].StencilLength].UserData=&(BoundaryTriangleFaces[nface].UserData);
    Stencil[InterpolationBase[nnode].StencilOffset+InterpolationBase[nnode].StencilLength].Weight=weight;
    InterpolationBase[nnode].TotalInterpolationWeight+=weight;

    InterpolationBase[nnode].StencilLength++;
  }

  //normalize the interpolation weights
  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) for (nface=0;nface<InterpolationBase[i].StencilLength;nface++) {
    Stencil[InterpolationBase[i].StencilOffset+nface].Weight/=InterpolationBase[i].TotalInterpolationWeight;
  }

  //print the variable list
  fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"");


  fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",CutCell::nBoundaryTriangleNodes,CutCell::nBoundaryTriangleFaces);

  //print the list of the data points
  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) {
    fprintf(fout,"%e %e %e ",CutCell::BoundaryTriangleNodes[i].x[0],CutCell::BoundaryTriangleNodes[i].x[1],CutCell::BoundaryTriangleNodes[i].x[2]);


/*    //prepare and output averaged interpolated face data

    //Prepare the interpolation stencil
    for (nface=0;nface<InterpolationBase[i].StencilLength;nface++) {
      InterpolationWeightList[nface]=Stencil[InterpolationBase[i].StencilOffset+nface].Weight;
      InterpolationFaceList[nface]=Stencil[InterpolationBase[i].StencilOffset+nface].UserData;
    }

    //output the user-defined surface data
    BoundaryTriangleFaces[0].UserData.Print(fout,InterpolationWeightList,InterpolationFaceList,InterpolationBase[i].StencilLength);*/


    fprintf(fout,"\n");
  }

  //print the connectovoty list
  for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++) {
    fprintf(fout,"%i %i %i\n",1+CutCell::BoundaryTriangleFaces[nface].node[0]->id,
        1+CutCell::BoundaryTriangleFaces[nface].node[1]->id,
        1+CutCell::BoundaryTriangleFaces[nface].node[2]->id);
  }


  fclose(fout);
}

//=============================================================================
//print the variable list
void cPostProcess3D::cParticleTrajectory::PrintVariableList() {
  if (PostProcess3D->rank!=0) return;
  printf("Particle Trajectory: Variable List: BEGIN\n");

  for (int nvar=0;nvar<nTrajectoryVariables;nvar++) printf("%i:\t%s\n",nvar,VariableList[nvar].c_str());
  printf("Particle Trajectory: Variable List: END\n\n");
}


//=============================================================================
//print particle trajectories
void cPostProcess3D::PrintParticleTrajectory(int nTrajectories,int OutputMode,double (*TrajectoryAcceptableProbability)(int),const char* fname) {
  int n;
  double MaxProbability=-1.0;

  if (rank!=0) return;

  //determine the probability table of the chosing a particulat trajectory to be printed
  double ProbabilityTable[ParticleTrajectory.nTotalTrajectories];

  if (OutputMode==_OUTPUT_MODE__UNIFORM_) {
    for (n=0;n<ParticleTrajectory.nTotalTrajectories;n++) ProbabilityTable[n]=1.0;
    MaxProbability=1.0;
  }
  else if (OutputMode==_OUTPUT_MODE__FLUX_) {
    if (TrajectoryAcceptableProbability==NULL) exit(__LINE__,__FILE__,"Error: TrajectoryAcceptableProbability must be defined");

    //the weight is ~ ((w/td)/vel)
    for (n=0;n<ParticleTrajectory.nTotalTrajectories;n++) {
      ProbabilityTable[n]=TrajectoryAcceptableProbability(n);
      if (MaxProbability<ProbabilityTable[n]) MaxProbability=ProbabilityTable[n];
    }

  }
  else exit(__LINE__,__FILE__,"Error: the option is unknown");

  //output the header of the trajectory file
  ParticleTrajectory.PrintDataFileHeader(fname);

  for (int i=0;i<nTrajectories;i++) {
    //determine the trajectory to output
    vector <int> PrintedTrajectories;
    bool found=false;
    int ii;

    do {
      found=false;

      do {
        n=rnd()*ParticleTrajectory.nTotalTrajectories;
      }
      while (ProbabilityTable[n]/MaxProbability<rnd());

      //check whether the trajecory is not chosen already
      for (ii=0;ii<PrintedTrajectories.size();ii++) if (PrintedTrajectories[ii]==n) {
        found=true;
        break;
      }
    }
    while (found==true);

    //print out new trajectory
    PrintedTrajectories.push_back(n);
    ParticleTrajectory.AddTrajectoryDataFile(&ParticleTrajectory.IndividualTrajectories[n],n,fname);
  }
}

//=============================================================================
//assign individual particle trajectories to the cells
void cPostProcess3D::AssignParticleTrajectoriesToCells() {
  double dt=0.0,vmax=-1.0,dxCellMin=-1.0,dtIntegration=0.0;
  int iInterval,nTrajectory,iBlock;

  //check whether the binary file that corresponds to the data file 'fname' exists. If the file exists -> compare the time
  if (access("post-processing.trajectory-cell-distribution.tmp.bin",R_OK)==0) {
    //the binary file exists: read the trajectory assignement from that file
    FILE *fBinaryIn=fopen("post-processing.trajectory-cell-distribution.tmp.bin","r"); //the file is removed when a new trajectory or a data file is read
    int iBlock,i,j,k;
    int nCellTrajectories,nTrajectory,n;

    for (iBlock=0;iBlock<nBlocks;iBlock++) for (i=0;i<nBlockCellX;i++) for (j=0;j<nBlockCellY;j++) for (k=0;k<nBlockCellZ;k++) {

      //read trajectory points in the cell
      fread(&nCellTrajectories,sizeof(int),1,fBinaryIn);

      for (nTrajectory=0;nTrajectory<nCellTrajectories;nTrajectory++) {
        fread(&n,sizeof(int),1,fBinaryIn);
        Block[iBlock].cell[i][j][k].TrajectoryPoints.push_back(n);
      }

      //read individual trajectories
      fread(&nCellTrajectories,sizeof(int),1,fBinaryIn);

      for (nTrajectory=0;nTrajectory<nCellTrajectories;nTrajectory++) {
        fread(&n,sizeof(int),1,fBinaryIn);
        Block[iBlock].cell[i][j][k].IndividualTrajectories.push_back(n);
      }
    }

    fclose(fBinaryIn);
    return;
  }

  //sinchronize
  MPI_Barrier(MPI_COMM_WORLD);

  //determine the trajectory path integration time step -> determine the minimum cell size, and the maximum particle velocity
  //the minimum cell size:
  for (iBlock=0;iBlock<nBlocks;iBlock++) {
    cBlock* bl;
    double l;

    bl=Block+iBlock;
    l=sqrt(pow(bl->xmax[0]-bl->xmin[0],2)+pow(bl->xmax[1]-bl->xmin[1],2)+pow(bl->xmax[2]-bl->xmin[2],2));
    if ((dxCellMin<0.0)||(dxCellMin>l)) dxCellMin=l;
  }

  dxCellMin/=std::min(std::min(nBlockCellX,nBlockCellY),nBlockCellZ);

  //the maximum particle speed
  for (nTrajectory=0;nTrajectory<ParticleTrajectory.nTotalTrajectories;nTrajectory++) {
    for (iInterval=0;iInterval<ParticleTrajectory.IndividualTrajectories[nTrajectory].nDataPoints;iInterval++) {
      if (vmax<ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[iInterval][4]) {
        vmax=ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[iInterval][4];
      }
    }
  }

  dtIntegration=dxCellMin/vmax;

  //distribute the trajectories
  double x0[3],x1[3],xParticle[3],ParticleSpeed,TrajectoryIntervalLength,dtLeftIntegrationStep,IntervalTravelTime,l[3];
  int i;
  int nTrajectoriesPerThread,nTrajectoryStart,nTrajectoryFinish;

  nTrajectoriesPerThread=ParticleTrajectory.nTotalTrajectories/size;

  nTrajectoryStart=rank*nTrajectoriesPerThread;
  nTrajectoryFinish=nTrajectoryStart+nTrajectoriesPerThread-1;

  if (rank==size-1) nTrajectoryFinish=ParticleTrajectory.nTotalTrajectories-1;

  for (nTrajectory=nTrajectoryStart;nTrajectory<=nTrajectoryFinish;nTrajectory++) if (ParticleTrajectory.IndividualTrajectories[nTrajectory].nDataPoints!=1) {
    for (i=0;i<3;i++) {
      x0[i]=ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][i];
      xParticle[i]=x0[i];
    }

    dtLeftIntegrationStep=dtIntegration;
    iInterval=1;

    //register the first point
    GetCell(x0)->TrajectoryPoints.push_back(nTrajectory);

    //get the information for the first segment
    for (i=0;i<3;i++) x1[i]=ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[1][i];

    TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
    for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;

    ParticleSpeed=0.5*(ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[0][4]+
        ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[1][4]);

    IntervalTravelTime=TrajectoryIntervalLength/ParticleSpeed;

    //start the integration loop
    do {
      if (IntervalTravelTime>dtLeftIntegrationStep) {
        IntervalTravelTime-=dtLeftIntegrationStep;

        //update the location of the particle and the remaining size of the interval
        for (i=0;i<3;i++) xParticle[i]+=dtLeftIntegrationStep*ParticleSpeed*l[i];

        //the time integratio step is finished -> register the point
        cCell* cl=GetCell(xParticle);
        cl->TrajectoryPoints.push_back(nTrajectory);
        dtLeftIntegrationStep=dtIntegration;

        //determine whether the trajectory is already saved with in the cell
        if (cl->LastTrajectoryProcessed!=nTrajectory) {
          //the trajectory is not saved in the cell yet -> save it
          cl->IndividualTrajectories.push_back(nTrajectory);
          cl->LastTrajectoryProcessed=nTrajectory;
        }
      }
      else {
        //the partice reached the end of the trajectory segment before the end of the time integraion step
        dtLeftIntegrationStep-=IntervalTravelTime;

        for (i=0;i<3;i++) xParticle[i]=x1[i];

        //move to the next interval
        iInterval++;
        if (iInterval==ParticleTrajectory.IndividualTrajectories[nTrajectory].nDataPoints) {
          //the integration has reached the end of the particle trajectory line -> stop the integration
          break;
        }

        for (i=0;i<3;i++) {
          x0[i]=x1[i];
          x1[i]=ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[iInterval][i];
        }

        TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
        for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;

        ParticleSpeed=0.5*(ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[iInterval-1][4]+
            ParticleTrajectory.IndividualTrajectories[nTrajectory].Data[iInterval][4]);

        IntervalTravelTime=TrajectoryIntervalLength/ParticleSpeed;
      }
    }
    while (true);
  }

  //output the trtajectory data into a file
  //the binary file exists: read the trajectory assignement from that file
  FILE *fBinaryOut=NULL;
  int j,k;
  int nCellTrajectories,n,nTotalTrajectories,thread;

  if (rank==0) fBinaryOut=fopen("post-processing.trajectory-cell-distribution.tmp.bin","w");

  for (iBlock=0;iBlock<nBlocks;iBlock++) for (i=0;i<nBlockCellX;i++) for (j=0;j<nBlockCellY;j++) for (k=0;k<nBlockCellZ;k++) {
    int nCellTrajectoriesVector[size];

    //save the trajectory points
    nCellTrajectories=Block[iBlock].cell[i][j][k].TrajectoryPoints.size();
    MPI_Allgather(&nCellTrajectories,1,MPI_INT,nCellTrajectoriesVector,1,MPI_INT,MPI_COMM_WORLD);

    for (nTotalTrajectories=0,thread=0;thread<size;thread++) nTotalTrajectories+=nCellTrajectoriesVector[thread];
    if (rank==0) fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (thread=0;thread<size;thread++) {
      int TrajectoryList[nCellTrajectoriesVector[thread]];

      if (thread==rank) {
        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          TrajectoryList[nTrajectory]=Block[iBlock].cell[i][j][k].TrajectoryPoints[nTrajectory];

        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);
      }
      else {
        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);

        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          Block[iBlock].cell[i][j][k].TrajectoryPoints.push_back(TrajectoryList[nTrajectory]);
      }

      if (rank==0) fwrite(TrajectoryList,sizeof(int),nCellTrajectoriesVector[thread],fBinaryOut);
    }

    //save the individual trajectories
    nCellTrajectories=Block[iBlock].cell[i][j][k].IndividualTrajectories.size();
    MPI_Allgather(&nCellTrajectories,1,MPI_INT,nCellTrajectoriesVector,1,MPI_INT,MPI_COMM_WORLD);

    for (nTotalTrajectories=0,thread=0;thread<size;thread++) nTotalTrajectories+=nCellTrajectoriesVector[thread];
    if (rank==0) fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (thread=0;thread<size;thread++) {
      int TrajectoryList[nCellTrajectoriesVector[thread]];

      if (thread==rank) {
        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          TrajectoryList[nTrajectory]=Block[iBlock].cell[i][j][k].IndividualTrajectories[nTrajectory];

        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);
      }
      else {
        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);

        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          Block[iBlock].cell[i][j][k].IndividualTrajectories.push_back(TrajectoryList[nTrajectory]);
      }

      if (rank==0) fwrite(TrajectoryList,sizeof(int),nCellTrajectoriesVector[thread],fBinaryOut);
    }

  }

  if (rank==0) fclose(fBinaryOut);
  MPI_Barrier(MPI_COMM_WORLD);
}












































