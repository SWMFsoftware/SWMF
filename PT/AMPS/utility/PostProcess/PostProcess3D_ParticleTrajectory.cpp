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
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <omp.h>
#include <boost/iostreams/device/mapped_file.hpp> // for mmap                           
#include <algorithm>  // for std::find                                                      
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>



int cPostProcess3D::cParticleTrajectory::nTrajectoryVariables=0;
int cPostProcess3D::cParticleTrajectory::ParticleWeightOverTimeStepOffset=-1;

//allocate individual trajectory data buffer
void cPostProcess3D::cParticleTrajectory::cIndividualTrajectoryData::AllocateDataBuffer(int n) {
  
  if(Data!=NULL){
    delete [] Data[0];
    delete [] Data;
  }
  
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

void cPostProcess3D::cParticleTrajectory::InitTrajectoryVarList(const char* filename){
  std::string line;
  std::ifstream infile(filename);
  std::string var_name,skip;
 
  
  std::getline(infile,line);
  std::istringstream ss(line);
  while (std::getline(ss,skip,'\"')) {
      std::getline(ss,var_name,'\"');
    VariableList.push_back(var_name); 
  }
  nTrajectoryVariables = VariableList.size();

  printf("Initialize Trajecotry Variable of Size %d\n",VariableList.size());
  
  int i;
  for ( i=0; i<VariableList.size()-1; i++) std::cout<<VariableList[i]<<',';
  std::cout<<VariableList[i]<<'\n';
}

void cPostProcess3D::cParticleTrajectory::GetTrajectoryLineNumber(int & traj_num, std::vector<int> & ZoneLineNum, const char* filename){
  std::string line;
  std::ifstream infile(filename);

  int cnt=0;
  while (std::getline(infile,line) ){
    cnt++;
    if (line[0]=='Z') {
      ZoneLineNum.push_back(cnt);
    }
  }
 
  cnt++; // to get ( total number of lines +1)
  ZoneLineNum.push_back(cnt);  

  traj_num=ZoneLineNum.size()-1;
  printf("the total number of trajecotries:%d\n",traj_num);

}


void cPostProcess3D::cParticleTrajectory::GetTrajectoryLocation_boost(int & traj_num, std::vector<const char *> & ZoneCharLoc,  boost::iostreams::mapped_file & mmap){

  const char * f = mmap.const_data();
  const char * l = f + mmap.size(); // point to end of the file
   
  while (f && f!=l){
    if ((f = static_cast<const char*>(memchr(f, 'Z', l-f)))){
      ZoneCharLoc.push_back(f);
      f++;
    }

  }

  traj_num=ZoneCharLoc.size();
  ZoneCharLoc.push_back(NULL);
  printf("the total number of trajecotries:%d\n",traj_num);

}

//head_traj_num: the index number of the first buffer traj
//buffer_traj_size: the total number of traj in the buffer 
void cPostProcess3D::cParticleTrajectory::WriteBinaryBuffer_boost(const char * bufferName, int head_traj_num, int buffer_traj_size, const char ** ZoneCharLocArr, boost::iostreams::mapped_file & mmap){
  
  FILE* fBinaryOut=fopen(bufferName,"w");
  char BufferIndexName[500];
  sprintf(BufferIndexName,"%s.index",bufferName);
  FILE* fBinaryOutOffset=fopen(BufferIndexName,"w");
  long int offset;
  
  fwrite(&buffer_traj_size,sizeof(int),1,fBinaryOut);
  for(int i=0; i<buffer_traj_size; i++){
    offset= ftell(fBinaryOut);
    fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
    WriteIndividualTrajectory_boost(ZoneCharLocArr[head_traj_num+i],ZoneCharLocArr[head_traj_num+i+1], mmap, fBinaryOut);
  }
    
  fclose(fBinaryOutOffset);
  fclose(fBinaryOut);

}



void cPostProcess3D::cParticleTrajectory::AppendBinaryBuffer_boost(const char * bufferName, int head_traj_num, int add_buffer_traj_size, const char ** ZoneCharLocArr, boost::iostreams::mapped_file & mmap){
  
  FILE* fBinaryOut=fopen(bufferName,"ab");
  char BufferIndexName[500];
  sprintf(BufferIndexName,"%s.index",bufferName);
  FILE* fBinaryOutOffset=fopen(BufferIndexName,"ab");
  long int offset;
  
  // fwrite(&buffer_traj_size,sizeof(int),1,fBinaryOut);
  for(int i=0; i<add_buffer_traj_size; i++){
    offset= ftell(fBinaryOut);
    fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
    WriteIndividualTrajectory_boost(ZoneCharLocArr[head_traj_num+i], ZoneCharLocArr[head_traj_num+i+1], mmap, fBinaryOut);
  }
    
  fclose(fBinaryOutOffset);
  fclose(fBinaryOut);

  int nPrevTraj, nCurrTraj;
  FILE* fBinaryUpdate=fopen(bufferName,"r+b");
  fread(&nPrevTraj,sizeof(int),1, fBinaryUpdate);
  fseek(fBinaryUpdate, 0, 0);
  nCurrTraj = nPrevTraj + add_buffer_traj_size;
  fwrite(&nCurrTraj, sizeof(int),1, fBinaryUpdate);
  fclose(fBinaryUpdate);

  if (nCurrTraj==BufferFileSize){
    LastBufferFull=true;
  }else{
    LastBufferFull=false;
    AvailBufferSize = BufferFileSize - nCurrTraj;
  }
  
}

void cPostProcess3D::cParticleTrajectory::WriteIndividualTrajectory_boost(const char * ZoneBegin, const char * ZoneEnd,  boost::iostreams::mapped_file & mmap, FILE* fBinaryOut){
 
  std::string line;
  int nDataPoint=0;
  int nVar = nTrajectoryVariables;
  if (ZoneEnd==NULL) ZoneEnd=mmap.const_data()+mmap.size();
  
  
  const char * charPtr=static_cast<const char*>(memchr(ZoneBegin, '\n', ZoneEnd-ZoneBegin));
  ZoneBegin = charPtr++; //first char after "Zone....\n"
 

  while(charPtr && charPtr!=ZoneEnd)
    if ((charPtr = static_cast<const char*>(memchr(charPtr, '\n', ZoneEnd-charPtr))))
      nDataPoint++, charPtr++;
    
  fwrite(&nDataPoint,sizeof(int),1,fBinaryOut);
   
  double d[nDataPoint*nVar];
  double *ptr=d;
  char * endPtr;
  charPtr = ZoneBegin;
 
  for (int i=0; i<nDataPoint*nVar; i++){
    *ptr = strtod(charPtr, &endPtr);
    charPtr = static_cast<const char*>(endPtr);
    ptr++;
  }
  
  fwrite(d,sizeof(double),nDataPoint*nVar,fBinaryOut);
}

//initialize variables and store trajectory info in several binary buffer files
void cPostProcess3D::cParticleTrajectory::InitLoadBufferFile_boost(const char *fname,const char* path, int buffer_traj_num) {
  char FullName[50000]; 
  int nTrajectoryInFile=0;// number of traj in this individual datafile
  //get the full name of the data file and open the file
  sprintf(FullName,"%s/%s",path,fname);
  printf("%s\n",FullName);
  printf("InitLoadBufferFile\n");
  
  if (!initBufferFileSizeFlag){
    BufferFileSize = buffer_traj_num;
    initBufferFileSizeFlag = true;
  }
  
  if (access(FullName,R_OK)!=0) {
     printf("Cannot find the file:%s\n",FullName);
     exit(__LINE__,__FILE__);
  }
  
  if (VariableList.size()==0) InitTrajectoryVarList(FullName);
  printf("InitLoadBufferFile inittrajecotryvarlist\n");
  //To store line number of "ZONE" for each trajectory in the ascii file 
  std::vector <const char *> ZoneCharLoc;
  // get total trajectory number and trajectory's line number in file

  boost::iostreams::mapped_file mmap(FullName, boost::iostreams::mapped_file::readonly);
  GetTrajectoryLocation_boost(nTrajectoryInFile, ZoneCharLoc, mmap);

  nTotalTrajectories += nTrajectoryInFile;
 
  printf("InitLoadBufferFile gettrajectorylinenumber done\n");
  if (nTrajectoryInFile==0) return;

  const char ** ZoneCharLocArr=&ZoneCharLoc[0]; //pointer of the head of ZoneLineNumVec
  int TrajNumToLoad=nTrajectoryInFile;

  //for test below
  // nTotalTrajectories = 3*PostProcess3D->BufferFileSize;
  //TrajNumToLoad = nTotalTrajectories;
  

  int head_traj_num_buffer=0;

  int BufferFileNum = BufferNameArray.size();
 
  printf("InitLoadBufferFile begin loadbuffer \n");
  //traj num in the buffer
  int bufferSize;
  //fill available buffer first
  if (TrajNumToLoad>0 && !LastBufferFull){
    bufferSize = AvailBufferSize<TrajNumToLoad?AvailBufferSize:TrajNumToLoad;
    AppendBinaryBuffer_boost(BufferNameArray.back().c_str(), head_traj_num_buffer, bufferSize, ZoneCharLocArr, mmap);
    head_traj_num_buffer += bufferSize;
    TrajNumToLoad -= bufferSize;
  }
  
  if(TrajNumToLoad>0){
    int nBuff = (int)ceil((double)TrajNumToLoad/buffer_traj_num);
    
#pragma omp parallel for ordered default(none)  shared(nBuff, TrajNumToLoad,buffer_traj_num, BufferFileNum,path,head_traj_num_buffer,ZoneCharLocArr, mmap) 
    for (int i=0; i<nBuff; i++){
      char BufferName[500];
      sprintf(BufferName,"%s/BufferFile%d.bin",".",BufferFileNum+i);
#pragma omp ordered
      {
#pragma omp critical(buffernameArray)
      BufferNameArray.push_back(BufferName);
      }
      int bufferSize = buffer_traj_num;
      if (i== nBuff-1 && TrajNumToLoad%buffer_traj_num!=0) bufferSize = TrajNumToLoad%buffer_traj_num;
      WriteBinaryBuffer_boost(BufferName, head_traj_num_buffer+i*buffer_traj_num, bufferSize, ZoneCharLocArr, mmap);

#pragma omp critical(OffsetInBuffVec)
      TrajOffsetInBuffer.push_back(NULL);
    }
    if (TrajNumToLoad%buffer_traj_num!=0) {
      LastBufferFull = false;
    }else{
      LastBufferFull = true;
    }
    
    BufferFileNum += nBuff;
  }
  printf("InitLoadBufferFile  loadbuffer done\n");
  
}

//initialize variables and store trajectory info in several binary buffer files 
// with inputnTraj for test cases
void cPostProcess3D::cParticleTrajectory::InitLoadBufferFile_boost(const char *fname,const char* path, int buffer_traj_num, int inputnTraj) {
  char FullName[50000]; 
  int nTrajectoryInFile=0;// number of traj in this individual datafile
  //get the full name of the data file and open the file
  sprintf(FullName,"%s/%s",path,fname);
  printf("InitLoadBufferFile\n");
  
  if (!initBufferFileSizeFlag){
    BufferFileSize = buffer_traj_num;
    initBufferFileSizeFlag = true;
  }
  
  if (access(FullName,R_OK)!=0) {
     printf("Cannot find the file:%s\n",FullName);
     exit(__LINE__,__FILE__);
  }
  
  if (VariableList.size()==0) InitTrajectoryVarList(FullName);
  printf("InitLoadBufferFile inittrajecotryvarlist\n");
  //To store line number of "ZONE" for each trajectory in the ascii file 
  std::vector <const char *> ZoneCharLoc;
  // get total trajectory number and trajectory's line number in file

  // boost::iostreams::mapped_file mmap(FullName, boost::iostreams::mapped_file::readonly);
  //  GetTrajectoryLocation_boost(nTrajectoryInFile, ZoneCharLoc, mmap);
  
  nTrajectoryInFile = inputnTraj;
  nTotalTrajectories += nTrajectoryInFile;
  
 
  printf("InitLoadBufferFile gettrajectorylinenumber done\n");
  if (nTrajectoryInFile==0) return;

  const char ** ZoneCharLocArr=&ZoneCharLoc[0]; //pointer of the head of ZoneLineNumVec
  int TrajNumToLoad=nTrajectoryInFile;

  //for test below
  // nTotalTrajectories = 3*PostProcess3D->BufferFileSize;
  //TrajNumToLoad = nTotalTrajectories;
  

  int head_traj_num_buffer=0;

  int BufferFileNum = BufferNameArray.size();
 
  printf("InitLoadBufferFile begin loadbuffer \n");
  //traj num in the buffer
  int bufferSize;
  //fill available buffer first
  if (TrajNumToLoad>0 && !LastBufferFull){
    bufferSize = AvailBufferSize<TrajNumToLoad?AvailBufferSize:TrajNumToLoad;
    //  AppendBinaryBuffer_boost(BufferNameArray.back().c_str(), head_traj_num_buffer, bufferSize, ZoneCharLocArr, mmap);
    head_traj_num_buffer += bufferSize;
    TrajNumToLoad -= bufferSize;
  }
  
  if(TrajNumToLoad>0){
    int nBuff = (int)ceil((double)TrajNumToLoad/buffer_traj_num);
    
    //#pragma omp parallel for default(none)  shared(nBuff, TrajNumToLoad,buffer_traj_num, BufferFileNum,path,head_traj_num_buffer,ZoneCharLocArr, mmap)
#pragma omp parallel for ordered default(none)  shared(nBuff, TrajNumToLoad,buffer_traj_num, BufferFileNum,path,head_traj_num_buffer,ZoneCharLocArr)
    for (int i=0; i<nBuff; i++){
      char BufferName[500];
      sprintf(BufferName,"%s/BufferFile%d.bin",".",BufferFileNum+i);
#pragma omp ordered
      {
#pragma omp critical(buffernameArray)
      BufferNameArray.push_back(BufferName);
      }
      int bufferSize = buffer_traj_num;
      if (i== nBuff-1 && TrajNumToLoad%buffer_traj_num!=0) bufferSize = TrajNumToLoad%buffer_traj_num;
      //  WriteBinaryBuffer_boost(BufferName, head_traj_num_buffer+i*buffer_traj_num, bufferSize, ZoneCharLocArr, mmap);
#pragma omp critical(OffsetInBuffVec)
      TrajOffsetInBuffer.push_back(NULL);
    }
    if (TrajNumToLoad%buffer_traj_num!=0) {
      LastBufferFull = false;
    }else{
      LastBufferFull = true;
    }
    
    BufferFileNum += nBuff;
  }
  printf("InitLoadBufferFile  loadbuffer done\n");
  
}

//head_traj_num: the index number of the first buffer traj
//buffer_traj_size: the total number of traj in the buffer 
void cPostProcess3D::cParticleTrajectory::WriteBinaryBuffer(const char * bufferName, int head_traj_num, int buffer_traj_size, int* ZoneLineNum, std::ifstream &infile, int & currentLineNum){
  
  FILE* fBinaryOut=fopen(bufferName,"w");
  char BufferIndexName[500];
  sprintf(BufferIndexName,"%s.index",bufferName);
  FILE* fBinaryOutOffset=fopen(BufferIndexName,"w");
  long int offset;
  
  fwrite(&buffer_traj_size,sizeof(int),1,fBinaryOut);
  for(int i=0; i<buffer_traj_size; i++){
    offset= ftell(fBinaryOut);
    fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
    WriteIndividualTrajectory(ZoneLineNum[head_traj_num+i], ZoneLineNum[head_traj_num+i+1], infile, currentLineNum, fBinaryOut);
  }
    
  fclose(fBinaryOutOffset);
  fclose(fBinaryOut);
  if (buffer_traj_size==BufferFileSize){
    LastBufferFull=true;
  }else{
    LastBufferFull=false;
    AvailBufferSize = BufferFileSize - buffer_traj_size;
  }
}

void cPostProcess3D::cParticleTrajectory::AppendBinaryBuffer(const char * bufferName, int head_traj_num, int add_buffer_traj_size, int* ZoneLineNum, std::ifstream &infile, int & currentLineNum){
  
  FILE* fBinaryOut=fopen(bufferName,"ab");
  char BufferIndexName[500];
  sprintf(BufferIndexName,"%s.index",bufferName);
  FILE* fBinaryOutOffset=fopen(BufferIndexName,"ab");
  long int offset;
  
  // fwrite(&buffer_traj_size,sizeof(int),1,fBinaryOut);
  for(int i=0; i<add_buffer_traj_size; i++){
    offset= ftell(fBinaryOut);
    fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
    WriteIndividualTrajectory(ZoneLineNum[head_traj_num+i], ZoneLineNum[head_traj_num+i+1], infile, currentLineNum, fBinaryOut);
  }
    
  fclose(fBinaryOutOffset);
  fclose(fBinaryOut);

  int nPrevTraj, nCurrTraj;
  FILE* fBinaryUpdate=fopen(bufferName,"r+b");
  fread(&nPrevTraj,sizeof(int),1, fBinaryUpdate);
  fseek(fBinaryUpdate, 0, 0);
  nCurrTraj = nPrevTraj + add_buffer_traj_size;
  fwrite(&nCurrTraj, sizeof(int),1, fBinaryUpdate);
  fclose(fBinaryUpdate);

  if (nCurrTraj==BufferFileSize){
    LastBufferFull=true;
  }else{
    LastBufferFull=false;
    AvailBufferSize = BufferFileSize - nCurrTraj;
  }

}





void cPostProcess3D::cParticleTrajectory::WriteIndividualTrajectory(int LineBegin, int LineEnd, std::ifstream &infile, int& cnt, FILE* fBinaryOut){
 
  std::string line;
  int nDataPoint= LineEnd-LineBegin-1;
  int nVar = nTrajectoryVariables;
 
  fwrite(&nDataPoint,sizeof(int),1,fBinaryOut);
   
  double d[nDataPoint*nVar];
  double *ptr=d;
 
  while (std::getline(infile,line)){
    cnt++;
    if (cnt>LineBegin && cnt<LineEnd) {
      std::stringstream ss(line);
      while (ss >> *ptr) ptr++;
    }
    if (cnt>=LineEnd) break;
  }
  
  // if !infile is true, the end of file is reached 
  if(cnt==LineEnd || !infile) fwrite(d,sizeof(double),nDataPoint*nVar,fBinaryOut);
}

void cPostProcess3D::cParticleTrajectory::LoadIndividualTrajData(cIndividualTrajectoryData &traj, int iTraj){
  
  long int offset;
  int iBuff = iTraj/BufferFileSize;

  GetIndividualTrajOffset(iTraj, offset);

  ReadIndividualTrajFromBuffer(traj, BufferNameArray[iBuff].c_str(), offset);
  
}

void cPostProcess3D::cParticleTrajectory::LoadIndividualTrajDataOneLine(cIndividualTrajectoryData &traj, int iTraj){
  
  long int offset;
  int iBuff = iTraj/BufferFileSize;

  GetIndividualTrajOffset(iTraj, offset);

  ReadIndividualTrajOneLineFromBuffer(traj, BufferNameArray[iBuff].c_str(), offset);
  
}



void cPostProcess3D::cParticleTrajectory::GetIndividualTrajOffset(int iTraj, long int & offset){
  int iBuff = iTraj/BufferFileSize;
  int jTraj = iTraj%BufferFileSize;
  
  if (iTraj>=nTotalTrajectories) {
    printf("%d is larger than total trajectory number %d.\n",iTraj,nTotalTrajectories);
    exit(__LINE__,__FILE__);
  }
  
  if (TrajOffsetInBuffer[iBuff]==NULL) {
    int nTrajInBuffer = BufferFileSize;
    if (iBuff==BufferNameArray.size()-1 && nTotalTrajectories%BufferFileSize!=0)
      nTrajInBuffer = nTotalTrajectories%BufferFileSize;
    TrajOffsetInBuffer[iBuff] = new long int [nTrajInBuffer];
    char BufferIndexName[1000];
    sprintf(BufferIndexName,"%s.index",BufferNameArray[iBuff].c_str());
    FILE* fBinaryInOffset=fopen(BufferIndexName,"r");
    fread(TrajOffsetInBuffer[iBuff],sizeof(long int), nTrajInBuffer,fBinaryInOffset);
    fclose(fBinaryInOffset);  
  }
  
  offset = TrajOffsetInBuffer[iBuff][jTraj];
}


void cPostProcess3D::cParticleTrajectory::ReadIndividualTrajFromBuffer(cIndividualTrajectoryData &traj, const char * bufferName, long int offset){
 
  int nDataPoint;
  FILE* fBinaryIn=fopen(bufferName,"rb");
  int nVar = nTrajectoryVariables;
  
 
  //FILE *fout=fopen("./test.dat","w");
  
  
  fseek(fBinaryIn,offset,0);
  fread(&nDataPoint,sizeof(int),1,fBinaryIn);
  traj.nDataPoints = nDataPoint;
  // fprintf(fout,"%d\n",nDataPoint);
  traj.AllocateDataBuffer(nDataPoint);
  // for (int jDataPoint=0; jDataPoint<nDataPoint; jDataPoint++)
    fread(traj.Data[0], sizeof(double), nVar*nDataPoint,fBinaryIn);
  
  /*
  for(int j=0; j<nDataPoint; j++){
    for(int k=0; k<nVar; k++){
      fprintf(fout," %e",data[j][k]);
    }
    fprintf(fout,"\n");
  }
  */
  
  fclose(fBinaryIn);
  // fclose(fout);

}


void cPostProcess3D::cParticleTrajectory::ReadIndividualTrajOneLineFromBuffer(cIndividualTrajectoryData &traj, const char * bufferName, long int offset){
 
  int nDataPoint;
  FILE* fBinaryIn=fopen(bufferName,"rb");
  int nVar = nTrajectoryVariables;
  
 
  //FILE *fout=fopen("./test.dat","w");
  
  
  fseek(fBinaryIn,offset,0);
  fread(&nDataPoint,sizeof(int),1,fBinaryIn);
  traj.nDataPoints = nDataPoint;
  // fprintf(fout,"%d\n",nDataPoint);
  traj.AllocateDataBuffer(1);
  fread(traj.Data[0], sizeof(double), nVar,fBinaryIn);
  
  /*
  for(int j=0; j<nDataPoint; j++){
    for(int k=0; k<nVar; k++){
      fprintf(fout," %e",data[j][k]);
    }
    fprintf(fout,"\n");
  }
  */
  
  fclose(fBinaryIn);
  // fclose(fout);

}


void cPostProcess3D::cParticleTrajectory::LoadAllBufferOneLine(){

  int nBuffer =  BufferNameArray.size();
  IndividualTrajectories.clear();
  IndividualTrajectories.reserve(nTotalTrajectories);
  IndividualTrajectories.resize(nTotalTrajectories);
  for (int iBuffer=0; iBuffer < nBuffer; iBuffer++){
    LoadBufferDataOneLine(iBuffer,iBuffer*BufferFileSize);
  }

  printf("LoadAllBufferOneLine done\n");
}


void cPostProcess3D::cParticleTrajectory::LoadBufferDataOneLine(int iBuffer, int head_traj_num){
  const char * BufferName = BufferNameArray[iBuffer].c_str(); 
  //head_traj_num is the first buffer trajectory's id in the individualtrajectories  
  
  //IndividualTrajectories.clear();
  
  int buffer_traj_size;
  int nDataPoint;
  int nVar =nTrajectoryVariables;
  printf("nVar:%d\n",nVar);
  FILE* fBinaryIn=fopen(BufferName,"rb");
  printf("load buffer bufferName:%s\n",BufferName);
  fread(&buffer_traj_size,sizeof(int),1,fBinaryIn);
  printf("buffer_size:%d\n",buffer_traj_size);
  // IndividualTrajectories.resize(buffer_traj_size);
  //IndividualTrajectories.clear();
  // IndividualTrajectories.swap(IndividualTrajectories);
  // IndividualTrajectories.reserve(BufferFileSize);
  //IndividualTrajectories.resize(BufferFileSize); 


  for(int i=0; i<buffer_traj_size; i++){
    
    fread(&nDataPoint,sizeof(int),1,fBinaryIn);
   
    cIndividualTrajectoryData * TrajPtr = &IndividualTrajectories[i+head_traj_num];
    //  printf("i+head_traj_num:%d\n",i+head_traj_num);
    if (TrajPtr==NULL) {
      printf("null traj in list:%d\n ",i);
    }
    TrajPtr->nDataPoints = nDataPoint;
    // printf("nDataPoint:%d\n",nDataPoint);
    if (nDataPoint>0){
      TrajPtr->AllocateDataBuffer(1);
      fread(TrajPtr->Data[0], sizeof(double), nVar,fBinaryIn);
      fseek(fBinaryIn, (nDataPoint-1)*nVar*sizeof(double),SEEK_CUR);
    }
   
  
  }
  
  fclose(fBinaryIn);    
  printf("Load Buffer Done\n");
 
}

void cPostProcess3D::cParticleTrajectory::LoadBufferHeader(int iBuffer, FILE * & fBinaryIn){
  const char * BufferName = BufferNameArray[iBuffer].c_str(); 

  //IndividualTrajectories.clear();
  
  int nDataPoint, buffer_traj_size;
  int nVar =nTrajectoryVariables;
  printf("nVar:%d\n",nVar);
  fBinaryIn=fopen(BufferName,"rb");
  printf("load buffer bufferName:%s\n",BufferName);
  fread(&buffer_traj_size,sizeof(int),1,fBinaryIn);
  printf("buffer_size:%d\n",buffer_traj_size);
 

  IndividualTrajectories.reserve(BufferFileSize);
  IndividualTrajectories.resize(BufferFileSize); 

}

void cPostProcess3D::cParticleTrajectory::ReadSequentialTrajFromBuffer(int nTraj, FILE * & fBinaryIn){
  
  if (nTraj%BufferFileSize==0) LoadBufferHeader(nTraj/BufferFileSize,fBinaryIn);
  int nVar =nTrajectoryVariables;
  int nDataPoint;
  fread(&nDataPoint,sizeof(int),1,fBinaryIn);
  cIndividualTrajectoryData* TrajPtr = &IndividualTrajectories[nTraj%BufferFileSize];
  
  TrajPtr->nDataPoints = nDataPoint;
  TrajPtr->AllocateDataBuffer(nDataPoint);
  
  /*    for (int j=1;j<nDataPoint;j++){
	if((TrajPtr->Data[j]-TrajPtr->Data[j-1])!=9) printf("error Data[j] address\n");
	}
  */
  
  //    for (int j=0; j<nDataPoint; j++){
  int cnt=fread(TrajPtr->Data[0], sizeof(double), nVar*nDataPoint,fBinaryIn);
  
  
  // }
  if (cnt!=nVar*nDataPoint){
    printf("fread error nTraj:%d\n",nTraj);
    // exit(__LINE__,__FILE__);
  }
  
  if (nTraj%BufferFileSize==BufferFileSize-1 || nTraj==nTotalTrajectories-1 ) fclose(fBinaryIn);    
  
}


void cPostProcess3D::cParticleTrajectory::LoadBufferData(int iBuffer){
  const char * BufferName = BufferNameArray[iBuffer].c_str(); 

  //IndividualTrajectories.clear();
  
  int buffer_traj_size;
  int nDataPoint;
  int nVar =nTrajectoryVariables;
  printf("nVar:%d\n",nVar);
  FILE* fBinaryIn=fopen(BufferName,"rb");
  printf("load buffer bufferName:%s\n",BufferName);
  fread(&buffer_traj_size,sizeof(int),1,fBinaryIn);
  printf("buffer_size:%d\n",buffer_traj_size);
  // IndividualTrajectories.resize(buffer_traj_size);
  IndividualTrajectories.clear();
  // IndividualTrajectories.swap(IndividualTrajectories);
  IndividualTrajectories.reserve(BufferFileSize);
  IndividualTrajectories.resize(BufferFileSize); 


  for(int i=0; i<buffer_traj_size; i++){
    
    fread(&nDataPoint,sizeof(int),1,fBinaryIn);
    cIndividualTrajectoryData* TrajPtr = &IndividualTrajectories[i];
    if (TrajPtr==NULL) {
      printf("null traj in list:%d\n ",i);
    }
    TrajPtr->nDataPoints = nDataPoint;
    TrajPtr->AllocateDataBuffer(nDataPoint);
    
    /*    for (int j=1;j<nDataPoint;j++){
      if((TrajPtr->Data[j]-TrajPtr->Data[j-1])!=9) printf("error Data[j] address\n");
    }
    */

    //    for (int j=0; j<nDataPoint; j++){
   int cnt=fread(TrajPtr->Data[0], sizeof(double), nVar*nDataPoint,fBinaryIn);
   
 
      // }
   if (cnt!=nVar*nDataPoint){
     printf("fread error iBuffer:%d, jTraj:%d\n",iBuffer,i);
     // exit(__LINE__,__FILE__);
   }
  }
  
  fclose(fBinaryIn);    
  printf("Load Buffer Done\n");
 
}


void cPostProcess3D::cParticleTrajectory::ReadBinaryBuffer(const char * bufferName){

  int buffer_traj_size;
  int nDataPoint;
  int nVar =nTrajectoryVariables;
  FILE* fBinaryIn=fopen(bufferName,"rb");
  fread(&buffer_traj_size,sizeof(int),1,fBinaryIn);

  FILE *fout=fopen("./test.dat","w");
  fprintf(fout,"%d\n",buffer_traj_size);
  
  for(int i=0; i<buffer_traj_size; i++){
    fread(&nDataPoint,sizeof(int),1,fBinaryIn);
    fprintf(fout,"%d\n",nDataPoint);
   
    double ** data;
    data = new double* [nDataPoint];
    data[0] = new double [nDataPoint*nVar];
    
    for (int jDataPoint=1; jDataPoint<nDataPoint; jDataPoint++)
      data[jDataPoint]=data[jDataPoint-1]+nVar;
    for (int jDataPoint=0; jDataPoint<nDataPoint; jDataPoint++)
      fread(data[jDataPoint], sizeof(double), nVar,fBinaryIn);
    
    for(int j=0; j<nDataPoint; j++){
      for(int k=0; k<nVar; k++){
	fprintf(fout," %e",data[j][k]);
      }
      fprintf(fout,"\n");
    }
    
  }
  
  fclose(fBinaryIn);
  fclose(fout);
  
}

//initialize variables and store trajectory info in several binary buffer files
void cPostProcess3D::cParticleTrajectory::InitLoadBufferFile(const char *fname,const char* path, int buffer_traj_num, int inputnTraj) {
  char FullName[50000]; 
  int nTrajectoryInFile=0;// number of traj in this individual datafile
  //get the full name of the data file and open the file
  sprintf(FullName,"%s/%s",path,fname);
  printf("InitLoadBufferFile\n");
  
  if (!initBufferFileSizeFlag){
    BufferFileSize = buffer_traj_num;
    initBufferFileSizeFlag = true;
  }
  
  if (access(FullName,R_OK)!=0) {
     printf("Cannot find the file:%s\n",FullName);
     exit(__LINE__,__FILE__);
  }
  
  if (VariableList.size()==0) InitTrajectoryVarList(FullName);
  printf("InitLoadBufferFile inittrajecotryvarlist\n");
  //To store line number of "ZONE" for each trajectory in the ascii file 
  std::vector <int> ZoneLineNumVec;
  // get total trajectory number and trajectory's line number in file
  //GetTrajectoryLineNumber(nTrajectoryInFile, ZoneLineNumVec, FullName);
  
  nTrajectoryInFile = inputnTraj;
  nTotalTrajectories += nTrajectoryInFile;
  
 
  printf("InitLoadBufferFile gettrajectorylinenumber done\n");
  if (nTrajectoryInFile==0) return;

  int* ZoneLineNum=&ZoneLineNumVec[0]; //pointer of the head of ZoneLineNumVec
  int TrajNumToLoad=nTrajectoryInFile;

  //for test below
  // nTotalTrajectories = 3*PostProcess3D->BufferFileSize;
  //TrajNumToLoad = nTotalTrajectories;
  

  int head_traj_num_buffer=0;

  std::ifstream infile(FullName); // ifstream object of the ascii file
  int currentLineNum=0;// currentLineNum of the pointer in the file
  int BufferFileNum = BufferNameArray.size();
 
  printf("InitLoadBufferFile begin loadbuffer \n");
  
  while(TrajNumToLoad>0){
    //trajectory number in this buffer
    int bufferSize;
    
    if(!LastBufferFull){
      bufferSize = AvailBufferSize<TrajNumToLoad?AvailBufferSize:TrajNumToLoad;
      //  AppendBinaryBuffer(BufferNameArray.back().c_str(), head_traj_num_buffer, bufferSize, ZoneLineNum, infile, currentLineNum);    
    }else{
      char BufferName[500];
      sprintf(BufferName,"%s/BufferFile%d.bin",path,BufferFileNum);   
      BufferNameArray.push_back(BufferName);
      BufferFileNum++;
      
      bufferSize = buffer_traj_num<TrajNumToLoad?buffer_traj_num:TrajNumToLoad;
      //  WriteBinaryBuffer(BufferName, head_traj_num_buffer, bufferSize, ZoneLineNum, infile, currentLineNum);
      TrajOffsetInBuffer.push_back(NULL);
    }
    head_traj_num_buffer += bufferSize;
    TrajNumToLoad -= bufferSize;
  
  }
  printf("InitLoadBufferFile  loadbuffer done\n");
  
}



//initialize variables and store trajectory info in several binary buffer files
void cPostProcess3D::cParticleTrajectory::InitLoadBufferFile(const char *fname,const char* path, int buffer_traj_num) {
  char FullName[50000]; 
  int nTrajectoryInFile=0;// number of traj in this individual datafile
  //get the full name of the data file and open the file
  sprintf(FullName,"%s/%s",path,fname);
  printf("InitLoadBufferFile\n");
  
  if (!initBufferFileSizeFlag){
    BufferFileSize = buffer_traj_num;
    initBufferFileSizeFlag = true;
  }
  
  if (access(FullName,R_OK)!=0) {
     printf("Cannot find the file:%s\n",FullName);
     exit(__LINE__,__FILE__);
  }
  
  if (VariableList.size()==0) InitTrajectoryVarList(FullName);
  printf("InitLoadBufferFile inittrajecotryvarlist\n");
  //To store line number of "ZONE" for each trajectory in the ascii file 
  std::vector <int> ZoneLineNumVec;
  // get total trajectory number and trajectory's line number in file
  GetTrajectoryLineNumber(nTrajectoryInFile, ZoneLineNumVec, FullName);
  nTotalTrajectories += nTrajectoryInFile;
  
  printf("InitLoadBufferFile gettrajectorylinenumber done\n");
  if (nTrajectoryInFile==0) return;

  int* ZoneLineNum=&ZoneLineNumVec[0]; //pointer of the head of ZoneLineNumVec
  int TrajNumToLoad=nTrajectoryInFile;

  //for test below
  // nTotalTrajectories = 3*PostProcess3D->BufferFileSize;
  //TrajNumToLoad = nTotalTrajectories;
  

  int head_traj_num_buffer=0;

  std::ifstream infile(FullName); // ifstream object of the ascii file
  int currentLineNum=0;// currentLineNum of the pointer in the file
  int BufferFileNum = BufferNameArray.size();
 
  printf("InitLoadBufferFile begin loadbuffer \n");
  
  while(TrajNumToLoad>0){
    //trajectory number in this buffer
    int bufferSize;
    
    if(!LastBufferFull){
      bufferSize = AvailBufferSize<TrajNumToLoad?AvailBufferSize:TrajNumToLoad;
      AppendBinaryBuffer(BufferNameArray.back().c_str(), head_traj_num_buffer, bufferSize, ZoneLineNum, infile, currentLineNum);    
    }else{
      char BufferName[500];
      sprintf(BufferName,"%s/BufferFile%d.bin",path,BufferFileNum);   
      BufferNameArray.push_back(BufferName);
      BufferFileNum++;
      
      bufferSize = buffer_traj_num<TrajNumToLoad?buffer_traj_num:TrajNumToLoad;
      WriteBinaryBuffer(BufferName, head_traj_num_buffer, bufferSize, ZoneLineNum, infile, currentLineNum);
      TrajOffsetInBuffer.push_back(NULL);
    }
    head_traj_num_buffer += bufferSize;
    TrajNumToLoad -= bufferSize;
  
  }
  printf("InitLoadBufferFile  loadbuffer done\n");
  
}

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
void cPostProcess3D::cParticleTrajectory::PrintSurfaceData(const char *fname,int (*GetVariableNumber)(),void (*PrintVariableList)(FILE*),void (*GetFaceDataVector)(double*,CutCell::cTriangleFace*,int)) {
  FILE *fout=NULL;
  int i,nface;

  if (PostProcess3D->rank==0) fout=fopen(fname,"w");

  struct cNodelInterpolationBase {
    double TotalInterpolationWeight;
    int StencilLength;
    int StencilOffset;
  };

  struct cFaceInterpolationStencil {
    double Weight;
    int nface;
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

    Stencil[InterpolationBase[nnode].StencilOffset+InterpolationBase[nnode].StencilLength].Weight=weight;
    Stencil[InterpolationBase[nnode].StencilOffset+InterpolationBase[nnode].StencilLength].nface=nface;
    InterpolationBase[nnode].TotalInterpolationWeight+=weight;

    InterpolationBase[nnode].StencilLength++;
  }

  //normalize the interpolation weights
  for (i=0;i<CutCell::nBoundaryTriangleNodes;i++) for (nface=0;nface<InterpolationBase[i].StencilLength;nface++) {
    Stencil[InterpolationBase[i].StencilOffset+nface].Weight/=InterpolationBase[i].TotalInterpolationWeight;
  }

  //print the variable list
  if (PostProcess3D->rank==0) {
    fprintf(fout,"VARIABLES=\"X\",\"Y\",\"Z\"");
    if (PrintVariableList!=NULL) PrintVariableList(fout);

    fprintf(fout,"\nZONE N=%i, E=%i, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n",CutCell::nBoundaryTriangleNodes,CutCell::nBoundaryTriangleFaces);
  }

  //print the list of the data points
  int iMeshNode,iInterplationElement;


  int ivar,nvars=0;
  double *data=NULL,*InterpolatedDataLocal=NULL,*InterpolatedDataGlobal=NULL;

  if (GetVariableNumber!=NULL) {
    nvars=GetVariableNumber();

    data=new double[nvars];
    InterpolatedDataLocal=new double[nvars];
    InterpolatedDataGlobal=new double[nvars];
  }

  for (iMeshNode=0;iMeshNode<CutCell::nBoundaryTriangleNodes;iMeshNode++) {
    if (PostProcess3D->rank==0)
      fprintf(fout,"%e %e %e ",CutCell::BoundaryTriangleNodes[iMeshNode].x[0],CutCell::BoundaryTriangleNodes[iMeshNode].x[1],CutCell::BoundaryTriangleNodes[iMeshNode].x[2]);

    if (GetVariableNumber!=NULL) {
    //prepare and output averaged interpolated face data
  //    int ivar,nvars=GetVariableNumber();
  //    double data[nvars],InterpolatedDataLocal[nvars],InterpolatedDataGlobal[nvars];

      for (ivar=0;ivar<nvars;ivar++) InterpolatedDataLocal[ivar]=0.0;

      //interpolate the data vector
      for (iInterplationElement=0;iInterplationElement<InterpolationBase[iMeshNode].StencilLength;iInterplationElement++) if (iInterplationElement%PostProcess3D->size==PostProcess3D->rank) {
        int nface=Stencil[InterpolationBase[iMeshNode].StencilOffset+iInterplationElement].nface;
        double w=Stencil[InterpolationBase[iMeshNode].StencilOffset+iInterplationElement].Weight;

        GetFaceDataVector(data,CutCell::BoundaryTriangleFaces+nface,nface);
        for (ivar=0;ivar<nvars;ivar++) InterpolatedDataLocal[ivar]+=w*data[ivar];
      }

      //collect the interpolted data from all processors
      MPI_Reduce(InterpolatedDataLocal,InterpolatedDataGlobal,nvars,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      //print the data into a file
      if (PostProcess3D->rank==0) for (ivar=0;ivar<nvars;ivar++) fprintf(fout," %e",InterpolatedDataGlobal[ivar]); 
    }

    if (PostProcess3D->rank==0) fprintf(fout,"\n");
  }

  //print the connectivity list
  if (PostProcess3D->rank==0) {
     for (nface=0;nface<CutCell::nBoundaryTriangleFaces;nface++)
       fprintf(fout,"%i %i %i\n",1+CutCell::BoundaryTriangleFaces[nface].node[0]->id,
         1+CutCell::BoundaryTriangleFaces[nface].node[1]->id,
         1+CutCell::BoundaryTriangleFaces[nface].node[2]->id);


     fclose(fout);
   }

  //deallocate the fata buffers
  delete [] InterpolationBase;
  delete [] Stencil;

  if (GetVariableNumber!=NULL) {
    delete [] data;
    delete [] InterpolatedDataLocal;
    delete [] InterpolatedDataGlobal;
  }
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
      
      if(ParticleTrajectory.BufferFileSize!=0){
	/*	int iBuffer = n/ParticleTrajectory.BufferFileSize;
	int jTraj = n%ParticleTrajectory.BufferFileSize;
	if (jTraj==0)  ParticleTrajectory.LoadBufferData(iBuffer);
	ProbabilityTable[n]=TrajectoryAcceptableProbability(jTraj);
	*/
	if (n==0) ParticleTrajectory.LoadAllBufferOneLine();
	ProbabilityTable[n]=TrajectoryAcceptableProbability(n);

	
	if (MaxProbability<ProbabilityTable[n]) MaxProbability=ProbabilityTable[n];
	
      }else{
      ProbabilityTable[n]=TrajectoryAcceptableProbability(n);
      if (MaxProbability<ProbabilityTable[n]) MaxProbability=ProbabilityTable[n];
      }
    }
    
  }
  else exit(__LINE__,__FILE__,"Error: the option is unknown");

  //output the header of the trajectory file
  ParticleTrajectory.PrintDataFileHeader(fname);
  printf("In the print trajectory\n");
  
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
//    PrintedTrajectories.push_back(n);
    
    if(ParticleTrajectory.BufferFileSize!=0){
      cParticleTrajectory::cIndividualTrajectoryData t;
      ParticleTrajectory.LoadIndividualTrajData(t,n);
      ParticleTrajectory.AddTrajectoryDataFile(&t,n,fname);
    }else{
    ParticleTrajectory.AddTrajectoryDataFile(&ParticleTrajectory.IndividualTrajectories[n],n,fname);
    }
    
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

  static int *TrajectoryList=NULL;
  static int TrajectoryListLength=0;

  if (rank==0) fBinaryOut=fopen("post-processing.trajectory-cell-distribution.tmp.bin","w");

  for (iBlock=0;iBlock<nBlocks;iBlock++) for (i=0;i<nBlockCellX;i++) for (j=0;j<nBlockCellY;j++) for (k=0;k<nBlockCellZ;k++) {
    int nCellTrajectoriesVector[size];

    //save the trajectory points
    nCellTrajectories=Block[iBlock].cell[i][j][k].TrajectoryPoints.size();
    MPI_Allgather(&nCellTrajectories,1,MPI_INT,nCellTrajectoriesVector,1,MPI_INT,MPI_COMM_WORLD);

    for (nTotalTrajectories=0,thread=0;thread<size;thread++) nTotalTrajectories+=nCellTrajectoriesVector[thread];
    if (rank==0) fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (thread=0;thread<size;thread++) {
//      int TrajectoryList[nCellTrajectoriesVector[thread]];

      if (nCellTrajectoriesVector[thread]>TrajectoryListLength) {
        TrajectoryListLength=nCellTrajectoriesVector[thread];
        if (TrajectoryList!=NULL) delete [] TrajectoryList;
        TrajectoryList=new int [TrajectoryListLength];
      }

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

//=============================================================================
//assign individual particle trajectories to the cells
void cPostProcess3D::AssignParticleTrajectoriesToCellsOptimize() {
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
	    
	    cCell* CellPtr = &Block[iBlock].cell[i][j][k];
	    CellPtr->TrajectoryPoints.resize(nCellTrajectories);
	    int tempData[nCellTrajectories];
	    fread(tempData,sizeof(int),nCellTrajectories,fBinaryIn);
	    for (nTrajectory=0;nTrajectory<nCellTrajectories;nTrajectory++)
	      CellPtr->TrajectoryPoints[nTrajectory] = tempData[nTrajectory];
	    
	    //read individual trajectories
	    fread(&nCellTrajectories,sizeof(int),1,fBinaryIn);
	    CellPtr->IndividualTrajectories.resize(nCellTrajectories);
	    fread(tempData,sizeof(int),nCellTrajectories,fBinaryIn);
	    for (nTrajectory=0;nTrajectory<nCellTrajectories;nTrajectory++)
	      CellPtr->IndividualTrajectories[nTrajectory] = tempData[nTrajectory];
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
    int iBuffer = nTrajectory/ParticleTrajectory.BufferFileSize;// i-th Buffer
    int jTraj = nTrajectory%ParticleTrajectory.BufferFileSize;// j-th trajectory
    
    if (jTraj==0)  ParticleTrajectory.LoadBufferData(iBuffer);
    
    for (iInterval=0;iInterval<ParticleTrajectory.IndividualTrajectories[jTraj].nDataPoints;iInterval++) {
      if (vmax<ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][4]) {
        vmax=ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][4];
      }
    }
  }

  dtIntegration=dxCellMin/vmax;

  //distribute the trajectories
  double x0[3],x1[3],xParticle[3],ParticleSpeed,TrajectoryIntervalLength,dtLeftIntegrationStep,IntervalTravelTime,l[3];
  int i;
  int nTrajectoriesPerThread,nTrajectoryStart,nTrajectoryFinish;

  // nTrajectoriesPerThread=ParticleTrajectory.nTotalTrajectories/size;
  nTrajectoriesPerThread=1/size;
  nTrajectoryStart=rank*nTrajectoriesPerThread;
  nTrajectoryFinish=nTrajectoryStart+nTrajectoriesPerThread-1;

  if (rank==size-1) nTrajectoryFinish=ParticleTrajectory.nTotalTrajectories-1;

  for (nTrajectory=nTrajectoryStart;nTrajectory<=nTrajectoryFinish;nTrajectory++){
    int iBuffer = nTrajectory/ParticleTrajectory.BufferFileSize;// i-th Buffer
    int jTraj = nTrajectory%ParticleTrajectory.BufferFileSize;// j-th trajectory
    
    if (jTraj==0)  ParticleTrajectory.LoadBufferData(iBuffer);
  
    if (ParticleTrajectory.IndividualTrajectories[jTraj].nDataPoints!=1) {
      for (i=0;i<3;i++) {
	x0[i]=ParticleTrajectory.IndividualTrajectories[jTraj].Data[0][i];
	xParticle[i]=x0[i];
      }
      
      dtLeftIntegrationStep=dtIntegration;
      iInterval=1;
      
      //register the first point
      GetCell(x0)->TrajectoryPoints.push_back(nTrajectory);
      
      //get the information for the first segment
      for (i=0;i<3;i++) x1[i]=ParticleTrajectory.IndividualTrajectories[jTraj].Data[1][i];
      
      TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
      for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;
      
      ParticleSpeed=0.5*(ParticleTrajectory.IndividualTrajectories[jTraj].Data[0][4]+
			 ParticleTrajectory.IndividualTrajectories[jTraj].Data[1][4]);

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
	  if (iInterval==ParticleTrajectory.IndividualTrajectories[jTraj].nDataPoints) {
	    //the integration has reached the end of the particle trajectory line -> stop the integration
	    break;
	  }

	  for (i=0;i<3;i++) {
	    x0[i]=x1[i];
	    x1[i]=ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][i];
	  }

	  TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
	  for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;

	  ParticleSpeed=0.5*(ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval-1][4]+
			     ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][4]);

	  IntervalTravelTime=TrajectoryIntervalLength/ParticleSpeed;
	}
      }
      while (true);
    } // if (ParticleTrajectory.IndividualTrajectories[nTtraj]!=1
  } //for (nTrajectory=nTrajectoryStart;nTrajectory<=nTrajectoryFinish;nTrajectory++) 
  
  //output the trtajectory data into a file
  //the binary file exists: read the trajectory assignement from that file
  FILE *fBinaryOut=NULL;
  int j,k;
  int nCellTrajectories,n,nTotalTrajectories,thread;

  static int *TrajectoryList=NULL;
  static int TrajectoryListLength=0;

  if (rank==0) fBinaryOut=fopen("post-processing.trajectory-cell-distribution.tmp.bin","w");

  for (iBlock=0;iBlock<nBlocks;iBlock++) for (i=0;i<nBlockCellX;i++) for (j=0;j<nBlockCellY;j++) for (k=0;k<nBlockCellZ;k++) {
    int nCellTrajectoriesVector[size];

    cCell* CellPtr = &Block[iBlock].cell[i][j][k];
    //save the trajectory points
    nCellTrajectories=CellPtr->TrajectoryPoints.size();
    MPI_Allgather(&nCellTrajectories,1,MPI_INT,nCellTrajectoriesVector,1,MPI_INT,MPI_COMM_WORLD);

    for (nTotalTrajectories=0,thread=0;thread<size;thread++) nTotalTrajectories+=nCellTrajectoriesVector[thread];
    if (rank==0) fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (thread=0;thread<size;thread++) {
//      int TrajectoryList[nCellTrajectoriesVector[thread]];

      if (nCellTrajectoriesVector[thread]>TrajectoryListLength) {
        TrajectoryListLength=nCellTrajectoriesVector[thread];
        if (TrajectoryList!=NULL) delete [] TrajectoryList;
        TrajectoryList=new int [TrajectoryListLength];
      }

      if (thread==rank) {
        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          TrajectoryList[nTrajectory]=CellPtr->TrajectoryPoints[nTrajectory];

        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);
      }
      else {
        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);

        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          CellPtr->TrajectoryPoints.push_back(TrajectoryList[nTrajectory]);
      }

      if (rank==0) fwrite(TrajectoryList,sizeof(int),nCellTrajectoriesVector[thread],fBinaryOut);
    }

    //save the individual trajectories
    nCellTrajectories=CellPtr->IndividualTrajectories.size();
    MPI_Allgather(&nCellTrajectories,1,MPI_INT,nCellTrajectoriesVector,1,MPI_INT,MPI_COMM_WORLD);

    for (nTotalTrajectories=0,thread=0;thread<size;thread++) nTotalTrajectories+=nCellTrajectoriesVector[thread];
    if (rank==0) fwrite(&nTotalTrajectories,sizeof(int),1,fBinaryOut);

    for (thread=0;thread<size;thread++) {
      int TrajectoryList[nCellTrajectoriesVector[thread]];

      if (thread==rank) {
        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          TrajectoryList[nTrajectory]=CellPtr->IndividualTrajectories[nTrajectory];

        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);
      }
      else {
        MPI_Bcast(TrajectoryList,nCellTrajectoriesVector[thread],MPI_INT,thread,MPI_COMM_WORLD);

        for (nTrajectory=0;nTrajectory<nCellTrajectoriesVector[thread];nTrajectory++)
          CellPtr->IndividualTrajectories.push_back(TrajectoryList[nTrajectory]);
      }

      if (rank==0) fwrite(TrajectoryList,sizeof(int),nCellTrajectoriesVector[thread],fBinaryOut);
    }

  }

  if (rank==0) fclose(fBinaryOut);
  MPI_Barrier(MPI_COMM_WORLD);
}

void cPostProcess3D::GetCellTrajInfo(){

 
  CellTrajMap = new std::map<int,int> [nCells];
  printf("nCells:%d\n",nCells);
  for (int iBuffer=0; iBuffer<ParticleTrajectory.BufferNameArray.size();iBuffer++){
    ReadTrajCellData(iBuffer);
   }
  
}

int  cPostProcess3D::sumTrajPoints(int iCell){
  std::map<int,int> * mapPtr = &CellTrajMap[iCell];
  if (mapPtr->size()==0) return 0;
  int sum=0;

  for (std::map<int,int>::iterator it=mapPtr->begin(); it!=mapPtr->end(); it++){
    sum += it->second;
  }

  return sum;
}

void cPostProcess3D::ReadTrajCellData(int iBuffer){
  char BufferTrajCellInfo[500];
  sprintf(BufferTrajCellInfo,"%s.OneTrajNCell",ParticleTrajectory.BufferNameArray[iBuffer].c_str());
  FILE* fBinaryIn=fopen(BufferTrajCellInfo,"rb");

  fseek(fBinaryIn,0, SEEK_END);
  long int endpos = ftell(fBinaryIn);
  
  rewind(fBinaryIn);
  printf("filename:%s\n",BufferTrajCellInfo);

  int buffer_traj_size = ParticleTrajectory.BufferFileSize;
  if (iBuffer == ParticleTrajectory.BufferNameArray.size()-1 && ParticleTrajectory.nTotalTrajectories%ParticleTrajectory.BufferFileSize!=0)
    buffer_traj_size= ParticleTrajectory.nTotalTrajectories%ParticleTrajectory.BufferFileSize;

  printf("buffer_traj_size:%d\n",buffer_traj_size);

  std::map<int, int> * CellTrajMapTmp = CellTrajMap;
  for(int i=0; i<buffer_traj_size; i++){
    int iTraj = iBuffer*ParticleTrajectory.BufferFileSize+i;
    int n;
    fread(&n, sizeof(int),1,fBinaryIn);

    if (n>0){
      int mapArr[2*n];
      int cnt = fread(mapArr, sizeof(int),2*n,fBinaryIn);
      if (cnt!=2*n) {
	printf("reading error\n");
	printf("error:%d\n",ferror(fBinaryIn));
	printf("eof:%d\n",feof(fBinaryIn));
	printf("diff:%d\n", endpos- ftell(fBinaryIn) );
      }

    
#pragma omp parallel for default(none) shared(n, mapArr, CellTrajMapTmp,iTraj) 
      for (int j=0; j<n; j++){
	int iCell = mapArr[2*j];
	int nDataPoints = mapArr[2*j+1];
	
	CellTrajMapTmp[iCell][iTraj]=nDataPoints;
      }
    }
    
  }
  
  fclose(fBinaryIn);

}


void cPostProcess3D::ComputeDtIntegration(double & dtIntegration){
  
  double dxCellMin=-1.0, vmax=-1.0;
  //determine the trajectory path integration time step -> determine the minimum cell size, and the maximum particle velocity
  //the minimum cell size:
  for (int iBlock=0;iBlock<nBlocks;iBlock++) {
    cBlock* bl;
    double l;

    bl=Block+iBlock;
    l=sqrt(pow(bl->xmax[0]-bl->xmin[0],2)+pow(bl->xmax[1]-bl->xmin[1],2)+pow(bl->xmax[2]-bl->xmin[2],2));
    if ((dxCellMin<0.0)||(dxCellMin>l)) dxCellMin=l;
  }

  dxCellMin/=std::min(std::min(nBlockCellX,nBlockCellY),nBlockCellZ);
  printf("dxCellMin:%f\n",dxCellMin);
  //the maximum particle speed
  /* 
     for (int nTrajectory=0;nTrajectory<ParticleTrajectory.nTotalTrajectories;nTrajectory++) {
    int iBuffer = nTrajectory/ParticleTrajectory.BufferFileSize;// i-th Buffer
    int jTraj = nTrajectory%ParticleTrajectory.BufferFileSize;// j-th trajectory
    
    if (jTraj==0)  ParticleTrajectory.LoadBufferData(iBuffer);
    
    for (int iInterval=0;iInterval<ParticleTrajectory.IndividualTrajectories[jTraj].nDataPoints;iInterval++) {
      if (vmax<ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][4]) {
        vmax=ParticleTrajectory.IndividualTrajectories[jTraj].Data[iInterval][4];
      }
    }
  }
  */

  int nBuff = ceil((double)ParticleTrajectory.nTotalTrajectories/ParticleTrajectory.BufferFileSize);
 
  for (int iBuffer = 0; iBuffer< nBuff; iBuffer++){
    ParticleTrajectory.LoadBufferData(iBuffer);
    int bufferSize = ParticleTrajectory.BufferFileSize;
    if (iBuffer==nBuff-1 && ParticleTrajectory.nTotalTrajectories%ParticleTrajectory.BufferFileSize!=0)   bufferSize = ParticleTrajectory.nTotalTrajectories%ParticleTrajectory.BufferFileSize;

    std::vector<cParticleTrajectory::cIndividualTrajectoryData> & TrajPtr= ParticleTrajectory.IndividualTrajectories;
    
    double vmax_traj=-1.0;
#pragma omp parallel for reduction(max:vmax_traj) default(none) shared(bufferSize, TrajPtr)
    for (int jTraj=0; jTraj<bufferSize; jTraj++){ 
      //   vmax_traj = -1.0; // re-init for openmp
      for (int iInterval=0;iInterval<TrajPtr[jTraj].nDataPoints;iInterval++) {
	if (vmax_traj<TrajPtr[jTraj].Data[iInterval][4]) {
	  vmax_traj=TrajPtr[jTraj].Data[iInterval][4];
	}
      }
    }
    
    if (vmax_traj>vmax) vmax = vmax_traj;
    
  }
  
  printf("vmax:%f\n",vmax);
  dtIntegration=dxCellMin/vmax;
  printf("dtIntegration:%f\n",dtIntegration);

}


/*
void cPostProcess3D::AddOneTrajToCell(int iCell, int iTraj, int nTrajPoint){


  
  CellTrajMap[iCell][iTraj]=nTrajPoint;

}
*/


void cPostProcess3D::AddOneCellTrajPointToTraj(std::map<int, int>  &TrajCellMap, double* x){

  int iCell = CellId(x);
 
  if (iCell<0 || iCell>=nCells)  exit(__LINE__,__FILE__); 
  std::map<int, int>::iterator  it=TrajCellMap.find(iCell);
  if (it!=TrajCellMap.end())  {
    (it->second)++;
  }else{
    TrajCellMap[iCell]=1;
  }

}

void cPostProcess3D::cParticleTrajectory::SaveAllTrajsCellInfo(){

  double dtIntegration;
  
  PostProcess3D->ComputeDtIntegration(dtIntegration);
 
  PostProcess3D->TrajCellMap = new std::map<int,int> [BufferFileSize];

  int nBuff = ceil((double)nTotalTrajectories/BufferFileSize);
 
  for (int iBuffer = 0; iBuffer< nBuff; iBuffer++){
    printf("SaveAllTrajsCellInfo,iBuffer:%d\n",iBuffer);
    LoadBufferData(iBuffer);
    int bufferSize = BufferFileSize;
    if (iBuffer==nBuff-1 && nTotalTrajectories%BufferFileSize!=0) 
      bufferSize = nTotalTrajectories%BufferFileSize;
 
    char BufferTrajCellInfo[500], BufferTrajCellIndexName[500];
    sprintf(BufferTrajCellInfo,"%s.OneTrajNCell",BufferNameArray[iBuffer].c_str());
    FILE* fBinaryOut=fopen(BufferTrajCellInfo,"wb");
    sprintf(BufferTrajCellIndexName,"%s.index",BufferTrajCellInfo);
    FILE* fBinaryOutOffset=fopen(BufferTrajCellIndexName,"wb");
 

    bool IsTrajDone[bufferSize];
    for (int jTraj=0; jTraj < bufferSize; jTraj++) IsTrajDone[jTraj]=false;


    std::map<int,int> * TrajCellMapTmp = PostProcess3D->TrajCellMap; 
    // cPostProcess3D * PostProcess3DTmp=PostProcess3D;
#pragma omp parallel default(none)  \
  shared(  IsTrajDone, fBinaryOut, fBinaryOutOffset,TrajCellMapTmp,dtIntegration, bufferSize, iBuffer)
    {
    
      if (omp_get_thread_num()!=0){  // slave threads processing data
	for (int jTraj=0; jTraj < bufferSize; jTraj++){	
	  if (jTraj%(omp_get_num_threads()-1)==omp_get_thread_num()-1){
	    ComputeOneTrajsCellInfo(jTraj, dtIntegration, TrajCellMapTmp[jTraj]);

#pragma omp critical (TrajDone)
	    IsTrajDone[jTraj] = true;
     	    
	  }
	}
      }else{  // master thread write data
	
	for (int iTraj=0; iTraj < bufferSize; iTraj++) {
	  //#pragma omp critical(TrajDone)

	  while (!IsTrajDone[iTraj]){
	    bool tmpFlag;
#pragma omp critical (TrajDone)
	    tmpFlag = IsTrajDone[iTraj];
	    if (tmpFlag) break; 
	  }

	  WriteTrajsCellInfoTraj(iTraj,fBinaryOut,fBinaryOutOffset, TrajCellMapTmp[iTraj]);  	 
	}
	fclose(fBinaryOutOffset);
	fclose(fBinaryOut);
	printf("TrajCell Buffer:%d done\n",iBuffer);
      }
      #pragma omp barrier
    }
     
    
  }

  delete [] PostProcess3D->TrajCellMap;
  
}

void cPostProcess3D::cParticleTrajectory::ReadIndividualTrajCellInfo(int nTraj){
  int iBuffer = nTraj/BufferFileSize;
  int jTraj = nTraj%BufferFileSize;

  char BufferTrajCellInfo[500], BufferTrajCellIndexName[500];
  sprintf(BufferTrajCellInfo,"%s.OneTrajNCell",BufferNameArray[iBuffer].c_str());
  FILE* fBinaryIn=fopen(BufferTrajCellInfo,"rb");
  sprintf(BufferTrajCellIndexName,"%s.index",BufferTrajCellInfo);
  FILE* fBinaryInOffset=fopen(BufferTrajCellIndexName,"rb");
  
  long int offset;
  fseek(fBinaryInOffset,jTraj*sizeof(long int),0);
  fread(&offset,sizeof(long int),1,fBinaryInOffset);
  fseek(fBinaryIn,offset,0);


  int n;
  fread(&n, sizeof(int), 1,fBinaryIn);
  
  int mapArr[2*n];
  fread(mapArr, sizeof(int),2*n,fBinaryIn);

  printf("Traj %d, jTraj %d Cell Info:\n", nTraj, jTraj);
  printf("Cell number %d\n", n);
  for (int j =0; j<n; j++){
    printf("%d, %d\n", mapArr[2*j],mapArr[2*j+1]);
  }

  fclose(fBinaryIn);
  fclose(fBinaryInOffset);
}

void cPostProcess3D::cParticleTrajectory::WriteTrajsCellInfoTraj(int jTraj, FILE * fBinaryOut, FILE * fBinaryOutOffset, std::map<int,int> & TrajCellMap){

  long int offset;
  
  offset= ftell(fBinaryOut);
  // printf("offset:%d\n",offset);
  fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
  std::map<int,int> * mapPtr;
  mapPtr = & TrajCellMap;
  int n = mapPtr->size();
  // printf("the cell number of traj:%d",n);
  fwrite(&n, sizeof(int), 1,fBinaryOut);
  if (n>0){
    int mapArr[2*n];
    int j=0;
    for (std::map<int,int>::iterator it=mapPtr->begin(); it != mapPtr->end();it++){
      mapArr[j++]=it->first;
      mapArr[j++]=it->second;
    }
    if (j!=2*n) {
      printf("jTraj:MapSize!=DataPoints\n");
      exit(__LINE__,__FILE__);
    }

    fwrite(mapArr, sizeof(int),2*n,fBinaryOut);
  }
  mapPtr->clear();
 
}




void cPostProcess3D::cParticleTrajectory::WriteTrajsCellInfoBuffer(int iBuffer){
  char BufferTrajCellInfo[500], BufferTrajCellIndexName[500];
  sprintf(BufferTrajCellInfo,"%s.OneTrajNCell",BufferNameArray[iBuffer].c_str());
  FILE* fBinaryOut=fopen(BufferTrajCellInfo,"w");
  sprintf(BufferTrajCellIndexName,"%s.index",BufferTrajCellInfo);
  FILE* fBinaryOutOffset=fopen(BufferTrajCellIndexName,"w");
  long int offset;
  
  int buffer_traj_size = BufferFileSize;
  if (iBuffer == BufferNameArray.size()-1 && nTotalTrajectories%BufferFileSize!=0 )
    buffer_traj_size = nTotalTrajectories%BufferFileSize;
  
  for(int i=0; i<buffer_traj_size; i++){
    offset= ftell(fBinaryOut);
    fwrite(&offset,sizeof(long int),1,fBinaryOutOffset); // record the offset of trajectory data
    std::map<int,int> * mapPtr;
    mapPtr = & PostProcess3D->TrajCellMap[i];
    int n = mapPtr->size();
    fwrite(&n, sizeof(long int), 1,fBinaryOut);
    if (n>0){
      int mapArr[2*n];
      int j=0;
      for (std::map<int,int>::iterator it=mapPtr->begin(); it != mapPtr->end();it++){
	mapArr[j++]=it->first;
	mapArr[j++]=it->second;
      }
      fwrite(mapArr, sizeof(int),2*n,fBinaryOut);
    }
    mapPtr->clear();
  }
    
  fclose(fBinaryOutOffset);
  fclose(fBinaryOut);
 
}



void cPostProcess3D::cParticleTrajectory::ComputeOneTrajsCellInfo(int nTraj, double dtIntegration, std::map<int, int> & TrajCellMap){
  
  double x0[3],x1[3], xParticle[3];
 

  //clear previous existing cell info
  //TrajCellMap stores the map as  (CellId, TrajPnts)
  //TrajPnts tells how many TrajPnts of Traj nTraj are in the cell with the id CellID. 
  int i;


  TrajCellMap.clear();
  
  // printf("nTraj:%d\n",nTraj);
  cIndividualTrajectoryData * t = &IndividualTrajectories[nTraj];
  // printf("t->nDataPoints:%d\n",t->nDataPoints);
  if (t->nDataPoints >1) {
    for (i=0;i<3;i++) {
      x0[i]=t->Data[0][i];
      xParticle[i]=x0[i];
    }
    
    int  iInterval=1;
    
    double dtLeftIntegration = dtIntegration;
    //register the first point
    
    //add 1 to TrajPnts in the tuple (CellId(x0), TrajPnts).
    PostProcess3D->AddOneCellTrajPointToTraj(TrajCellMap, x0);
    
    //get the information for the first segment
    for (i=0;i<3;i++) x1[i]=t->Data[1][i];
    
    double TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
    
    double l[3];
    for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;
    
    double ParticleSpeed=0.5*(t->Data[0][4]+t->Data[1][4]);
    
    double IntervalTravelTime=TrajectoryIntervalLength/ParticleSpeed;
    
    //start the integration loop
    do {
      if (IntervalTravelTime>dtLeftIntegration) {
	IntervalTravelTime-=dtLeftIntegration;

	//update the location of the particle and the remaining size of the interval
	for (i=0;i<3;i++) xParticle[i]+=dtLeftIntegration*ParticleSpeed*l[i];
	
	//the time integratio step is finished -> register the point
	//	cPostProcess3D::cCell* cl=PostProcess3D->GetCell(xParticle);
	  
	//add 1 to TrajPnts in the tuple (CellId(xParticle), TrajPnts).
	PostProcess3D->AddOneCellTrajPointToTraj(TrajCellMap, xParticle);
	dtLeftIntegration=dtIntegration;
      }
      else {
	  //the partice reached the end of the trajectory segment before the end of the time integraion step
	  dtLeftIntegration-=IntervalTravelTime;
	  
	  for (i=0;i<3;i++) xParticle[i]=x1[i];
	  
	  //move to the next interval
	  iInterval++;
	  //  printf("iInterval:%d\n",iInterval);
	
	  if (iInterval==t->nDataPoints) {
	    //the integration has reached the end of the particle trajectory line -> stop the integration
	    break;
	  }
	  
	  for (i=0;i<3;i++) {
	    x0[i]=x1[i];
	    x1[i]=t->Data[iInterval][i];
	  }
	  
	  TrajectoryIntervalLength=sqrt(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2));
	  for (i=0;i<3;i++) l[i]=(x1[i]-x0[i])/TrajectoryIntervalLength;
	  
	  ParticleSpeed=0.5*(t->Data[iInterval-1][4]+t->Data[iInterval][4]);
	  
	  IntervalTravelTime=TrajectoryIntervalLength/ParticleSpeed;
      }
    }
    while (true);
  }
  


}









































