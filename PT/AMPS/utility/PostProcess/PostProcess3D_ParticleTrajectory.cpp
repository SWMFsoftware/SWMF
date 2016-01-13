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


//=============================================================
//add trajectory data to the list of the trajectories
void cPostProcess3D::cParticleTrajectory::AddIndividualTrajectoryData(int& nDataPoints,std::vector<double>& data) {
  cIndividualTrajectoryData t;

  if (nDataPoints!=0) {
    t.nDataPoints=nDataPoints;
    t.Data=new double [nDataPoints*nTrajectoryVariables];

    for (int i=0;i<nDataPoints*nTrajectoryVariables;i++) t.Data[i]=data[i];

    //add the trajectory data
    IndividualTrajectories.push_back(t);
    nTotalTrajectories++;
  }

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

  //all slave processers will read the mesh only after the root processor finish reading
  if (PostProcess3D->rank!=0) MPI_Barrier(MPI_COMM_WORLD);

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
    if (PostProcess3D->rank==0) fBinaryOut=fopen(BinaryFullName,"w");
  }

  //read the variable line
  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);

  while (strcmp(str1,"")!=0) {
    std::string var(str1);
    VariableList.push_back(var);
    nTrajectoryVariables++;

    ifile.CutInputStr(str1,str);
  }

  //read the trajectory information
  if (fBinaryOut!=NULL) {
    while (ifile.eof()==false) {
      ifile.GetInputStr(str,sizeof(str));
      ifile.CutInputStr(str1,str);

      if (strcmp(str1,"ZONE")==0) {
        //beginig of the new trajectory: close the previous trajectory and prepare the buffer for the new trajectory
        AddIndividualTrajectoryData(nDataPoints,TrajectoryData);
      }
      else {
        //the line is a trajectory point
        for (int i=0;i<nTrajectoryVariables;i++) {
          TrajectoryData.push_back(atof(str1));
          ifile.CutInputStr(str1,str);
        }

        nDataPoints++;
      }
    }

    //save the last trajectory
    AddIndividualTrajectoryData(nDataPoints,TrajectoryData);

    //save the binary file with the trajectory data
    fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (int nTrajectory=0;nTrajectory<nTotalTrajectories;nTrajectory++) {
      cIndividualTrajectoryData t;

      t=IndividualTrajectories[nTrajectory];
      fwrite(&t.nDataPoints,sizeof(int),1,fBinaryOut);
      fwrite(t.Data,sizeof(double),t.nDataPoints*nTrajectoryVariables,fBinaryOut);
    }

    fclose(fBinaryOut);
  }
  else {
    //read the trajectory data from a binary file
    fread(&nTotalTrajectories,sizeof(int),1,fBinaryIn);

    for (int nTrajectory=0;nTrajectory<nTotalTrajectories;nTrajectory++) {
      cIndividualTrajectoryData t;

      fread(&t.nDataPoints,sizeof(int),1,fBinaryIn);
      t.Data=new double [t.nDataPoints*nTrajectoryVariables];

      fread(t.Data,sizeof(double),t.nDataPoints*nTrajectoryVariables,fBinaryIn);

      IndividualTrajectories.push_back(t);
    }
  }

  //close the binary files
  if (fBinaryIn!=NULL) fclose(fBinaryIn);
  if (fBinaryOut!=NULL) fclose(fBinaryOut);

  //the root processor will read the mesh before all other processors
  if (PostProcess3D->rank==0) MPI_Barrier(MPI_COMM_WORLD);

  //close the trajectory data file and save the binary trajectory data
  ifile.closefile();
}


//=====================================================================================
//save trajectory data file
void cPostProcess3D::cParticleTrajectory::PrintDataFileHeader(const char* fname) {
  FILE *fout=fopen(fname,"w");

  fprintf(fout,"VARIABLES=\"%s\"",VariableList[0].c_str());
  for (int nvar=1;nvar<VariableList.size();nvar++) fprintf(fout,", %s",VariableList[nvar].c_str());

  fprintf(fout,"\n");
  fclose(fout);
}

void cPostProcess3D::cParticleTrajectory::AddTrajectoryDataFile(cIndividualTrajectoryData* Trajectory,int TrajectoryNumber,const char* fname) {
  FILE *fout=fopen(fname,"a");

  fprintf(fout,"ZONE T=\"Trajectory=%i\" F=POINT",TrajectoryNumber);

  for (int nline=0;nline<Trajectory->nDataPoints;nline++) {
    for (int i=0;i<nTrajectoryVariables;i++) fprintf(fout," %e",Trajectory->Data[i+nline*nTrajectoryVariables]);
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

























