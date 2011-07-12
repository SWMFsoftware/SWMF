//==========================================================
//$Id$
//==========================================================
//the particle data buffer
#include <stdio.h>
#include <stdlib.h>

#include "pic.h"

long int PIC::ParticleBuffer::ParticleDataLength=_PIC_PARTICLE_DATA_POSITION_OFFSET_+sizeof(double)*DIM;
PIC::ParticleBuffer::byte *PIC::ParticleBuffer::ParticleDataBuffer=NULL;
long int PIC::ParticleBuffer::MaxNPart=0;
long int PIC::ParticleBuffer::NAllPart=0;
long int PIC::ParticleBuffer::FirstPBufferParticle=-1;

//==========================================================
//init the buffer
void PIC::ParticleBuffer::Init(long int BufrerLength) {

  if ((ParticleDataBuffer!=NULL)||(MaxNPart!=0)) exit(__LINE__,__FILE__,"Reallocation of the particle data buffer");
  if (sizeof(byte)!=1) exit(__LINE__,__FILE__,"The size of 'byte' is diferent from 1");
  if (BufrerLength<=0) exit(__LINE__,__FILE__,"BufrerLength is less that zero");

  MaxNPart=BufrerLength;
  ParticleDataBuffer=(PIC::ParticleBuffer::byte*) malloc(ParticleDataLength*MaxNPart);

  //init the list of particles in the buffer
  for (long int ptr=0;ptr<MaxNPart-1;ptr++) SetNext(ptr+1,ptr);
  FirstPBufferParticle=0;

}

//==========================================================
//Request additional data for a particle
void PIC::ParticleBuffer::RequestDataStorage(long int &offset,int TotalDataLength) {
  offset=ParticleDataLength;
  ParticleDataLength+=TotalDataLength;
}

//==========================================================
//the basic data access functions for a particle
PIC::ParticleBuffer::byte *PIC::ParticleBuffer::GetParticleDataPointer(long int ptr) {
  return ParticleDataBuffer+ptr*ParticleDataLength;
}

//==========================================================
//get the particle position
double *PIC::ParticleBuffer::GetX(long int ptr) {
  return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
}

double *PIC::ParticleBuffer::GetX(byte *ParticleDataStart) {
  return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
}

void PIC::ParticleBuffer::GetX(double* x,long int ptr) {
  register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
}

void PIC::ParticleBuffer::GetX(double* x,byte *ParticleDataStart) {
  register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) x[idim]=xptr[idim];
}

void PIC::ParticleBuffer::SetX(double* x,long int ptr) {
  register double *xptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];
}

void PIC::ParticleBuffer::SetX(double* x,byte *ParticleDataStart) {
  register double *xptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_POSITION_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) xptr[idim]=x[idim];
}

//==========================================================
//get the particle velocity
double *PIC::ParticleBuffer::GetV(long int ptr) {
  return (double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
}

double *PIC::ParticleBuffer::GetV(byte *ParticleDataStart) {
  return (double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
}

void PIC::ParticleBuffer::GetV(double* v,long int ptr) {
  register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
}

void PIC::ParticleBuffer::GetV(double* v,byte *ParticleDataStart) {
  register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) v[idim]=vptr[idim];
}

void PIC::ParticleBuffer::SetV(double* v,long int ptr) {
  register double *vptr=(double*) (ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
}

void PIC::ParticleBuffer::SetV(double* v,byte *ParticleDataStart) {
  register double *vptr=(double*) (ParticleDataStart+_PIC_PARTICLE_DATA_VELOCITY_OFFSET_);
  register int idim;

  for (idim=0;idim<DIM;idim++) vptr[idim]=v[idim];
}


//==========================================================
//get the particle's species ID

unsigned int PIC::ParticleBuffer::GetI(byte* ParticleDataStart) {
  return *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
}

unsigned int PIC::ParticleBuffer::GetI(long int ptr) {
  return *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_));
}

void PIC::ParticleBuffer::SetI(unsigned int spec,byte* ParticleDataStart) {
  *((unsigned char*)(ParticleDataStart+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=spec;
}

void PIC::ParticleBuffer::SetI(unsigned int spec,long int ptr) {
  *((unsigned char*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_SPECIEDID_OFFSET_))=spec;
}

//==========================================================
//get/set prev
long int PIC::ParticleBuffer::GetPrev(long int ptr) {
  return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_PREV_OFFSET_));
}

long int PIC::ParticleBuffer::GetPrev(byte* ParticleDataStart) {
  return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_PREV_OFFSET_));
}

void PIC::ParticleBuffer::SetPrev(long int prev,long int ptr) {
  *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_PREV_OFFSET_))=prev;
}

void PIC::ParticleBuffer::SetPrev(long int prev,byte* ParticleDataStart) {
  *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_PREV_OFFSET_))=prev;
}

//get/set next
long int PIC::ParticleBuffer::GetNext(long int ptr) {
  return *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_NEXT_OFFSET_));
}

long int PIC::ParticleBuffer::GetNext(byte* ParticleDataStart) {
  return *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_NEXT_OFFSET_));
}

void PIC::ParticleBuffer::SetNext(long int next,long int ptr) {
  *((long int*)(ParticleDataBuffer+ptr*ParticleDataLength+_PIC_PARTICLE_DATA_NEXT_OFFSET_))=next;
}

void PIC::ParticleBuffer::SetNext(long int next,byte* ParticleDataStart) {
  *((long int*)(ParticleDataStart+_PIC_PARTICLE_DATA_NEXT_OFFSET_))=next;
}


//==========================================================
//the functions that controls the particle buffer
long int PIC::ParticleBuffer::GetMaxNPart() {return MaxNPart;}
long int PIC::ParticleBuffer::GetAllPartNum() {return NAllPart;}
long int PIC::ParticleBuffer::GetParticleDataLength() {return ParticleDataLength;}

long int PIC::ParticleBuffer::GetNewParticle() {
  long int newptr;
  byte *pdataptr;

  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);

  FirstPBufferParticle=GetNext(pdataptr);
  SetPrev(-1,pdataptr);
  SetNext(-1,pdataptr);

  return newptr;
}

long int PIC::ParticleBuffer::GetNewParticle(long int &ListFirstParticle) {
  long int newptr;
  byte *pdataptr;

  if (MaxNPart==NAllPart) exit(__LINE__,__FILE__,"The particle buffer is full");

  NAllPart++;
  newptr=FirstPBufferParticle;
  pdataptr=GetParticleDataPointer(newptr);

  FirstPBufferParticle=GetNext(pdataptr);

  if (ListFirstParticle>=0) {
    byte *listFirstPData=GetParticleDataPointer(ListFirstParticle);

    SetPrev(GetPrev(listFirstPData),pdataptr);
    SetNext(ListFirstParticle,pdataptr);
    SetPrev(newptr,listFirstPData);
  }
  else {
    SetPrev(ListFirstParticle,pdataptr);
    SetNext(ListFirstParticle,pdataptr);
  }

  ListFirstParticle=newptr;
  return newptr;
}

void PIC::ParticleBuffer::ExcludeParticleFromList(long int ptr,long int& ListFirstParticle) {
  byte *pdataptr=GetParticleDataPointer(ptr);
  long int prev,next;

  //exclude the particle from the list
  prev=GetPrev(pdataptr);
  next=GetNext(pdataptr);

  if (ptr==ListFirstParticle) {
    SetPrev(prev,next);
    ListFirstParticle=next;
  }
  else {
    SetNext(next,prev);
    SetPrev(prev,next);
  }
}


void PIC::ParticleBuffer::DeleteParticle(long int ptr) {
  NAllPart--;
  SetNext(FirstPBufferParticle,ptr);
  FirstPBufferParticle=ptr;
}


void PIC::ParticleBuffer::DeleteParticle(long int ptr,long int& ListFirstParticle) {
  ExcludeParticleFromList(ptr,ListFirstParticle);
  DeleteParticle(ptr);
}



void PIC::ParticleBuffer::CloneParticle(long int copy,long int source) {
  byte *SourceData,*CopyData;
  long int i,next,prev;

  SourceData=GetParticleDataPointer(source);
  CopyData=GetParticleDataPointer(copy);

  prev=GetPrev(CopyData);
  next=GetNext(CopyData);

  for (i=0;i<ParticleDataLength;i++) CopyData[i]=SourceData[i];

  SetPrev(prev,CopyData);
  SetNext(next,CopyData);
}

//==========================================================
//save the particle buffer in a restart file
void PIC::ParticleBuffer::SaveImageFile(int fd) {
  exit(__LINE__,__FILE__,"Not implemented");
}

void PIC::ParticleBuffer::LoadImageFile(int fd) {
  exit(__LINE__,__FILE__,"not implemented");
}


//==========================================================
//pack the particle data

void PIC::ParticleBuffer::PackParticleData(char* buffer,long int ptr) {
  byte *SourceData=GetParticleDataPointer(ptr);
  long int i;

  for (i=0;i<ParticleDataLength;i++) buffer[i]=SourceData[i];
}


void PIC::ParticleBuffer::UnPackParticleData(char* buffer,long int ptr) {
  byte *pdata;
  long int i,next,prev;

  pdata=GetParticleDataPointer(ptr);
  prev=GetPrev(pdata);
  next=GetNext(pdata);

  for (i=0;i<ParticleDataLength;i++) pdata[i]=buffer[i];

  SetPrev(prev,pdata);
  SetNext(next,pdata);
}

//==========================================================
//get the checksum of the particle buffer
unsigned long PIC::ParticleBuffer::GetChecksum() {
  CRC32 sum;

  //save the particle's buffer internal data
  sum.add(&ParticleDataLength,1);
  sum.add(&MaxNPart,1);
  sum.add(&NAllPart,1);
  sum.add(&FirstPBufferParticle,1);

  //save the particle's data
  sum.add(ParticleDataBuffer,MaxNPart*ParticleDataLength);

  unsigned long int *buffer=new unsigned long int[TotalThreadsNumber];
  char str[10*_MAX_STRING_LENGTH_PIC_];

  buffer[0]=sum.checksum();
  MPI_Gather(buffer,1,MPI_UNSIGNED_LONG,buffer,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);

  if (ThisThread==0) {
    sprintf(str,"Cdsmc::pbuffer CRC32 checksum: ");

    for (long int thread=0;thread<TotalThreadsNumber;thread++) sprintf(str,"%s 0x%lx ",str,buffer[thread]);

    printf("%s\n",str);
    PrintErrorLog(str);
  }

  delete [] buffer;
  return sum.checksum();
}

