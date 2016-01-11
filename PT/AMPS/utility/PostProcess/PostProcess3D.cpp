//$Id$
//the library for reading of AMPS' output files

/*
 * PostProcess.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: vtenishe
 */

#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "PostProcess3D.h"
#include "ifileopr.h"

//init the data buffer
void cPostProcess3D::sDataStructure::InitDataBuffer(int nnodes,int nvars) {
  nNodes=nnodes,nVariables=nvars;

  //alllocate the memory
  data=new double* [nNodes];
  data[0]=new double [nNodes*nVariables];

  for (int i=1;i<nNodes;i++) data[i]=data[i-1]+nVariables;
}

void cPostProcess3D::PrintVariableList() {
  int i;
  if (rank!=0) return;

  for (i=0;i<data.nVariables;i++) printf("%i:\t%s\n",i,data.VariableList[i].c_str());
}

void cPostProcess3D::sDataStructure::CopyVariable(int iTarget,sDataStructure* Source,int iSource) {
  int i;

  for (i=0;i<nNodes;i++) data[i][iTarget]=Source->data[i][iSource];
}

void cPostProcess3D::cBlock::Init(int nx,int ny,int nz) {
  nCellX=nx,nCellY=ny,nCellZ=nz;

  cell=new cPostProcess3D::cCell**[nCellX];
  cell[0]=new cPostProcess3D::cCell* [nCellX*nCellY];
  cell[0][0]=new cPostProcess3D::cCell[nCellX*nCellY*nCellZ];

  int i,j,k,cnt;

  for (i=1;i<nCellX;i++) cell[i]=cell[i-1]+nCellY;

  for (cnt=0,i=0;i<nCellX;i++) for (j=0;j<nCellY;j++) {
    cell[i][j]=cell[0][0]+cnt;
    cnt+=nCellZ;
  }
}

//load the data file
void cPostProcess3D::LoadDataFile(const char *fname,const char* path) {
  char FullName[50000],str[50000],str1[50000];
  CiFileOperations ifile;
  int nvars=0,nline;
  FILE *fBinaryIn=NULL,*fBinaryOut=NULL;

  //all slave processers will read the mesh only after the root processor finish reading
  if (rank!=0) MPI_Barrier(MPI_COMM_WORLD);

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
    if (rank==0) fBinaryOut=fopen(BinaryFullName,"w");
  }


  //read the variable line
  ifile.GetInputStr(str,sizeof(str));
  ifile.CutInputStr(str1,str);
  ifile.CutInputStr(str1,str);

  while (strcmp(str1,"")!=0) {
    std::string var(str1);
    data.VariableList.push_back(var);
    nvars++;

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

  //read the data
  data.InitDataBuffer(nNodes,nvars);

  if (fBinaryIn==NULL) {
    //read the data file
    for (nline=0;nline<nNodes;nline++) {
      ifile.GetInputStr(str,sizeof(str));

      for (int i=0;i<nvars;i++) {
        ifile.CutInputStr(str1,str);
        data.data[nline][i]=atof(str1);
      }
    }

    //save the binary file
    if (rank==0) fwrite(data.data[0],sizeof(double),nNodes*nvars,fBinaryOut);
  }
  else {
    fread(data.data[0],sizeof(double),nNodes*nvars,fBinaryIn);
  }


  //save the connectivity list and initialize blocks
  int nblock,iBlock,jBlock,kBlock;
  unsigned long int cnt;
  ConnectivityList=new cCell[nCells];

  if (nCells%(nBlockCellX*nBlockCellY*nBlockCellZ)) {
    printf("Error: somthing is wrong...\n");
    exit(__LINE__,__FILE__);
  }

  nBlocks=nCells/(nBlockCellX*nBlockCellY*nBlockCellZ);
  Block=new cBlock[nBlocks];

  for (nline=0;nline<nCells;nline++) {
    if (nline%(nBlockCellX*nBlockCellY*nBlockCellZ)==0) {
      //new block;
      nblock=nline/(nBlockCellX*nBlockCellY*nBlockCellZ);
      cnt=0;
      Block[nblock].Init(nBlockCellX,nBlockCellY,nBlockCellZ);
    }

    //determine the coordinate of the cell in the block
    int r;

    kBlock=cnt/(nBlockCellX*nBlockCellY);

    r=cnt%(nBlockCellX*nBlockCellY);
    jBlock=r/nBlockCellX;
    iBlock=r%nBlockCellX;

    cnt++;

    if (fBinaryIn==NULL) {
      ifile.GetInputStr(str,sizeof(str));

      for (int i=0;i<8;i++) {
        ifile.CutInputStr(str1,str);
        ConnectivityList[nline].n[i]=strtol(str1,&endptr,10);
        Block[nblock].cell[iBlock][jBlock][kBlock].n[i]=strtol(str1,&endptr,10)-1;
      }

      Block[nblock].id=nblock;
      if (rank==0) fwrite(ConnectivityList[nline].n,sizeof(int),8,fBinaryOut);
    }
    else {
      Block[nblock].id=nblock;
      fread(ConnectivityList[nline].n,sizeof(int),8,fBinaryIn);
      for (int i=0;i<8;i++) Block[nblock].cell[iBlock][jBlock][kBlock].n[i]=ConnectivityList[nline].n[i]-1;
    }
  }

  //close the binary files
  if (fBinaryIn!=NULL) fclose(fBinaryIn);
  if (fBinaryOut!=NULL) fclose(fBinaryOut);

  //the root processor will read the mesh before all other processors
  if (rank==0) MPI_Barrier(MPI_COMM_WORLD);

  //find xmin amd xmax and dx for each block
  int i,j,k;

  for (nblock=0;nblock<nBlocks;nblock++) {
    for (i=0;i<3;i++) {
      Block[nblock].xmin[i]=data.data[Block[nblock].cell[0][0][0].n[0]][i];
      Block[nblock].xmax[i]=data.data[Block[nblock].cell[nBlockCellX-1][nBlockCellY-1][nBlockCellZ-1].n[6]][i];

      if ((xmin[i]>Block[nblock].xmin[i]) || (nblock==0)) xmin[i]=Block[nblock].xmin[i];
      if ((xmax[i]<Block[nblock].xmax[i]) || (nblock==0)) xmax[i]=Block[nblock].xmax[i];
    }

    //the size of the cells in the block
    Block[nblock].dx[0]=(Block[nblock].xmax[0]-Block[nblock].xmin[0])/nBlockCellX;
    Block[nblock].dx[1]=(Block[nblock].xmax[1]-Block[nblock].xmin[1])/nBlockCellY;
    Block[nblock].dx[2]=(Block[nblock].xmax[2]-Block[nblock].xmin[2])/nBlockCellZ;

    //the minimum size of the block in each dimension
    for (i=0;i<3;i++) {
      double l=Block[nblock].xmax[i]-Block[nblock].xmin[i];

      if ((dxBlockMin[i]>l) || (nblock==0)) dxBlockMin[i]=l;
    }
  }

  //Init Mesh


  dxBlockMin[0]=(xmax[0]-xmin[0])/_POST_PROCESS_MAP_SIZE_;
  dxBlockMin[1]=(xmax[1]-xmin[1])/_POST_PROCESS_MAP_SIZE_;
  dxBlockMin[2]=(xmax[2]-xmin[2])/_POST_PROCESS_MAP_SIZE_;

  //assign blocks to the Mesh elements
  for (nblock=0;nblock<nBlocks;nblock++) {
    int imin[3],imax[3],i,j,k;

    for (i=0;i<3;i++) {
      imin[i]=(Block[nblock].xmin[i]-xmin[i])/dxBlockMin[i];
      imax[i]=(Block[nblock].xmax[i]-xmin[i])/dxBlockMin[i];
    }

    for (i=imin[0];i<=imax[0];i++) for (j=imin[1];j<=imax[1];j++) for (k=imin[2];k<=imax[2];k++) {
      if ((i<_POST_PROCESS_MAP_SIZE_)&&(j<_POST_PROCESS_MAP_SIZE_)&&(k<_POST_PROCESS_MAP_SIZE_)) {
        Mesh[i][j][k].push_back(Block+nblock);
      }
    }
  }

  //create vector that containes only coordinated of the nodes
  X.InitDataBuffer(nNodes,3);
  for (i=0;i<3;i++) X.CopyVariable(i,&data,i);

  //close the data file
  ifile.closefile();

  //init the procedure of the column integration
  ColumnIntegral.Init(xmin,xmax);
}

//save the data file
void cPostProcess3D::SaveDataFile(const char *fname, sDataStructure &dat) {
  if (rank!=0) return;
  FILE *fout=fopen(fname,"w");

  //save the variable list, and the number of the nodes and cells
  int i,ncell,nnode=0;

  fprintf(fout,"VARIABLES=");
  for (std::vector<std::string>::iterator v=dat.VariableList.begin();v!=dat.VariableList.end();v++) {
    if (nnode!=0) fprintf(fout,", ");
    nnode++;

    fprintf(fout,"\"%s\"",v->c_str());
  }

  fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",nNodes,nCells);

  //output the data
  for (nnode=0;nnode<nNodes;nnode++) {
    for (i=0;i<dat.nVariables;i++) fprintf(fout,"%e ",dat.data[nnode][i]);
    fprintf(fout,"\n");
  }

  //output the connectovity list
  for (ncell=0;ncell<nCells;ncell++) {
    for (i=0;i<8;i++) fprintf(fout,"%i ",ConnectivityList[ncell].n[i]);
    fprintf(fout,"\n");
  }

  //close the data file
  fclose(fout);
}


//getermine whether a point is within the domain
bool cPostProcess3D::IsInsideDomain(double *x) {
  for (int i=0;i<3;i++) if ((x[i]<xmin[i])||(xmax[i]<x[i])) return false;

  return true;
}

//determine the block
cPostProcess3D::cBlock* cPostProcess3D::GetBlock(double *x) {
  int i,j,k,n,idim,nblocks;
  cBlock*  res=NULL;

  i=(x[0]-xmin[0])/dxBlockMin[0];
  j=(x[1]-xmin[1])/dxBlockMin[1];
  k=(x[2]-xmin[2])/dxBlockMin[2];

  if (i==_POST_PROCESS_MAP_SIZE_) i--;
  if (j==_POST_PROCESS_MAP_SIZE_) j--;
  if (k==_POST_PROCESS_MAP_SIZE_) k--;

  nblocks=Mesh[i][j][k].size();

  for (n=0;n<nblocks;n++) {
    bool found=true;

    for (idim=0;idim<3;idim++) if ((x[idim]<Mesh[i][j][k][n]->xmin[idim]) || (Mesh[i][j][k][n]->xmax[idim]<x[idim])) {
      found=false;
      goto ExitSearchLoop;
    }

ExitSearchLoop:
    if (found==true) {
      res=Mesh[i][j][k][n];
      break;
    }
  }

  if (res==NULL) exit(__LINE__,__FILE__,"Error: the block is not found");

  return res;
}

//characteristic size of the cell
double cPostProcess3D::CharacteristicCellSize(double *x) {
  cBlock* block;

  block=GetBlock(x);
  return sqrt(pow(block->dx[0],2)+pow(block->dx[1],2)+pow(block->dx[2],2));
}

//determine the interpolation stencil
void cPostProcess3D::GetInterpolationStencil(double *x,cStencil* stencil) {
  int i,j,k;
  cBlock *bl;

  //determine the block
  bl=GetBlock(x);

  //determine the cell in the block
  int iCell[3];
  double w[3],InterpolationWeight;

  for (i=0;i<3;i++) {
    iCell[i]=(x[i]-bl->xmin[i])/bl->dx[i];
    w[i]=(x[i]-bl->xmin[i])/bl->dx[i]-iCell[i];
  }

  for (i=0;i<2;i++) for (j=0;j<2;j++) for (k=0;k<2;k++) {
    switch (i+2*j+4*k) {
    case 0+0*2+0*4:
      InterpolationWeight=(1.0-w[0])*(1.0-w[1])*(1.0-w[2]);

      stencil->Node[0]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[0];
      stencil->Weight[0]=InterpolationWeight;

      break;
    case 1+0*2+0*4:
      InterpolationWeight=w[0]*(1.0-w[1])*(1.0-w[2]);

      stencil->Node[1]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[1];
      stencil->Weight[1]=InterpolationWeight;

      break;
    case 0+1*2+0*4:
      InterpolationWeight=(1.0-w[0])*w[1]*(1.0-w[2]);

      stencil->Node[3]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[3];
      stencil->Weight[3]=InterpolationWeight;

      break;
    case 1+1*2+0*4:
      InterpolationWeight=w[0]*w[1]*(1.0-w[2]);

      stencil->Node[2]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[2];
      stencil->Weight[2]=InterpolationWeight;

      break;

    case 0+0*2+1*4:
      InterpolationWeight=(1.0-w[0])*(1.0-w[1])*w[2];

      stencil->Node[4]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[4];
      stencil->Weight[4]=InterpolationWeight;

      break;
    case 1+0*2+1*4:
      InterpolationWeight=w[0]*(1.0-w[1])*w[2];

      stencil->Node[5]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[5];
      stencil->Weight[5]=InterpolationWeight;

      break;
    case 0+1*2+1*4:
      InterpolationWeight=(1.0-w[0])*w[1]*w[2];

      stencil->Node[7]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[7];
      stencil->Weight[7]=InterpolationWeight;

      break;
    case 1+1*2+1*4:
      InterpolationWeight=w[0]*w[1]*w[2];

      stencil->Node[6]=bl->cell[iCell[0]][iCell[1]][iCell[2]].n[6];
      stencil->Weight[6]=InterpolationWeight;

      break;
    }

  }

}
























