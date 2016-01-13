//$Id$

#include "PostProcess3D.h"
#include "ifileopr.h"


//load the data file
void cPostProcess3D::cSurfaceData::LoadDataFile(const char *fname,const char* path) {
  char FullName[50000],str[50000],str1[50000];
  CiFileOperations ifile;
  int nline;
  FILE *fBinaryIn=NULL,*fBinaryOut=NULL;

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

  sprintf(BinaryFullName,"%s.post-processing.tmp.bin",fname);

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
    nVariables++;

    ifile.CutInputStr(str1,str);
  }

  //read the number of the nodes and cells
  char *endptr;
  ifile.GetInputStr(str,sizeof(str));

  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);
  nNodes=strtol(str1,&endptr,10);


  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);
  nCells=strtol(str1,&endptr,10);

  //init the data buffer
  data=new double* [nNodes];
  data[0]=new double [nNodes*nVariables];

  for (nline=1;nline<nNodes;nline++) data[nline]=data[nline-1]+nVariables;

  //reas the data
  if (fBinaryIn==NULL) {
    for (nline=0;nline<nNodes;nline++) {
      ifile.GetInputStr(str,sizeof(str));

      for (int nvar=0;nvar<nVariables;nvar++) {
        ifile.CutInputStr(str1,str);
        data[nline][nvar]=atof(str1);
      }
    }

    fwrite(data[0],sizeof(double),nNodes*nVariables,fBinaryOut);
  }
  else {
    fread(data[0],sizeof(double),nNodes*nVariables,fBinaryIn);
  }

  //init the connectivity list
  int ncell,nd;

  ConnectivityList=new int* [nCells];
  ConnectivityList[0]=new int [3*nCells];

  for (ncell=1;ncell<nCells;ncell++) ConnectivityList[ncell]=ConnectivityList[ncell-1]+3;

  //read the connectivity list
  if (fBinaryIn==NULL) {
    for (ncell=0;ncell<nCells;ncell++) {
      ifile.GetInputStr(str,sizeof(str));

      for (nd=0;nd<3;nd++) {
        ifile.CutInputStr(str1,str);
        ConnectivityList[ncell][nd]=strtol(str1,&endptr,10)-1;
      }
    }

    fwrite(ConnectivityList[0],sizeof(int),3*nCells,fBinaryOut);
  }
  else {
    fread(ConnectivityList[0],sizeof(int),3*nCells,fBinaryIn);
  }

  //close the binary files
  if (fBinaryIn!=NULL) fclose(fBinaryIn);
  if (fBinaryOut!=NULL) fclose(fBinaryOut);

  //the root processor will read the mesh before all other processors
  if (PostProcess3D->rank==0) MPI_Barrier(MPI_COMM_WORLD);
}

//=============================================================================
//print the variable list
void cPostProcess3D::cSurfaceData::PrintVariableList() {
  if (PostProcess3D->rank!=0) return;
  printf("Surface Data: Variable List: BEGIN\n");

  for (int nvar=0;nvar<nVariables;nvar++) printf("%i:\t%s\n",nvar,VariableList[nvar].c_str());
  printf("Surface Data: Variable List: END\n\n");
}
   
