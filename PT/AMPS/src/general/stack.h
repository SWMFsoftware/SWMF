//$Id$

//declaration of general purpose stack class

#ifndef STACK
#define STACK

#include "specfunc.h"

#define _GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_ 10000
#define _GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_ 10000

//the stack class
template<class T>
class cStack {
public: 
  //data element's stack
  long int nMaxElements; 
  T*** elementStackList;
  long int elementStackPointer;

  //array pointers on the allocates data blocks
  T** dataBufferList;
  long int dataBufferListSize,dataBufferListPointer;  


  //control allocated memory
  long int MemoryAllocation;   

  long int getAllocatedMemory() {
    return MemoryAllocation;
  }
  
  void initMemoryBlock() {
    long int i,j;

    if (sizeof(T)==0) return;

    //check available space in the dataBufferList list: if needed increment the size of 'elementStackList' and 'dataBufferList' 
    if (dataBufferListPointer==dataBufferListSize) {
      T** tmpDataList=new T*[dataBufferListSize+_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_];
      MemoryAllocation+=sizeof(T*)*(dataBufferListSize+_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_);

      T*** tmpStackList=new T**[dataBufferListSize+_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_];
      MemoryAllocation+=sizeof(T**)*(dataBufferListSize+_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_); 

      i=0;
     
      //copy the content of the old lists to the new ones 
      if (dataBufferList!=NULL) for (;i<dataBufferListSize;i++) tmpDataList[i]=dataBufferList[i],tmpStackList[i]=elementStackList[i];
      for (j=0;j<_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_;j++,i++) tmpDataList[i]=NULL,tmpStackList[i]=NULL;

      if (dataBufferList!=NULL) {
        delete [] dataBufferList;
        MemoryAllocation-=sizeof(T*)*dataBufferListSize;

        delete [] elementStackList;
        MemoryAllocation-=sizeof(T**)*dataBufferListSize;
      }

      dataBufferList=tmpDataList;
      elementStackList=tmpStackList;
      dataBufferListSize+=_GENERAL_STACK_DEFAULT_BUFFER_LIST_SIZE_;
    }

    //allocate a new memory chunk for the element's data and update the stack list
    dataBufferList[dataBufferListPointer]=new T[_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
    elementStackList[dataBufferListPointer]=new T*[_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_];

    MemoryAllocation+=sizeof(T)*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    MemoryAllocation+=sizeof(T*)*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    for (i=0;i<_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;i++) 
      elementStackList[dataBufferListPointer][i]=dataBufferList[dataBufferListPointer]+i; 
 
    nMaxElements+=_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    dataBufferListPointer++;
  }

  void resetStack() {
    long int databank,offset;

    for (databank=0;databank<dataBufferListPointer;databank++) 
      for (offset=0;offset<_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;offset++) 
	elementStackList[databank][offset]=dataBufferList[databank]+offset; 

    elementStackPointer=0;
  }
 


  //get the entry pointer and counting number
  long int GetEntryCountingNumber(T* ptr) {
    long int nMemoryBank,res=-1;

    if (ptr!=NULL) {
      for (nMemoryBank=0;nMemoryBank<dataBufferListPointer;nMemoryBank++) if ((ptr>=dataBufferList[nMemoryBank])&&(ptr<dataBufferList[nMemoryBank]+_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_)) if ((res=ptr-dataBufferList[nMemoryBank])<_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_) {
        return res+nMemoryBank*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
      }
 
      exit(__LINE__,__FILE__,"Error: canot find the entry number");  
    }
 
    return -1;
  }

  T* GetEntryPointer(long int countingNumber) {
    long int nMemoryBank,offset;

    if ((countingNumber<0.0)||(countingNumber>=nMaxElements)) return NULL;

    nMemoryBank=countingNumber/_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=countingNumber-nMemoryBank*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    return dataBufferList[nMemoryBank]+offset;
  }

  void clear() {
    for (int i=0;i<dataBufferListPointer;i++) {
      delete [] dataBufferList[i];
      delete [] elementStackList[i];

      MemoryAllocation-=(sizeof(T)+sizeof(T*))*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    }

    if (dataBufferList!=NULL) {
      delete [] dataBufferList;
      delete [] elementStackList; 

      MemoryAllocation-=(sizeof(T*)+sizeof(T**))*dataBufferListSize;
    }

    nMaxElements=0,elementStackList=NULL,elementStackPointer=0;
    dataBufferList=0,dataBufferList=NULL,dataBufferListSize=0,dataBufferListPointer=0;
  } 


  //save and load the allocation of the stack
  void saveAllocationParameters(FILE *fout) {
    fwrite(&nMaxElements,sizeof(long int),1,fout);
    fwrite(&elementStackPointer,sizeof(long int),1,fout);
    fwrite(&dataBufferListSize,sizeof(long int),1,fout);
    fwrite(&dataBufferListPointer,sizeof(long int),1,fout);

    //save the elementStack 
    long int i,j,elementCountingNumber;

    for (i=0;i<dataBufferListPointer;i++) for (j=0;j<_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;j++) {    
      elementCountingNumber=GetEntryCountingNumber(elementStackList[i][j]);
      fwrite(&elementCountingNumber,sizeof(long int),1,fout); 
    }  
  }


  void readAllocationParameters(FILE *fout) {
    clear();

    fread(&nMaxElements,sizeof(long int),1,fout);
    fread(&elementStackPointer,sizeof(long int),1,fout);
    fread(&dataBufferListSize,sizeof(long int),1,fout);
    fread(&dataBufferListPointer,sizeof(long int),1,fout);

    //allocate the stack's buffers
    long int i,j,elementCountingNumber;

    elementStackList=new T** [dataBufferListSize];
    MemoryAllocation+=sizeof(T**)*dataBufferListSize;
    for (i=0;i<dataBufferListSize;i++) elementStackList[i]=NULL;

    dataBufferList=new T*[dataBufferListSize];
    MemoryAllocation+=sizeof(T*)*dataBufferListSize;
    for (i=0;i<dataBufferListSize;i++) dataBufferList[i]=NULL;

    for (i=0;i<dataBufferListPointer;i++) {
      dataBufferList[i]=new T[_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
      MemoryAllocation+=sizeof(T)*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

      elementStackList[i]=new T*[_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_]; 
      MemoryAllocation+=sizeof(T*)*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    }

    //read the elementStack
    for (i=0;i<dataBufferListPointer;i++) for (j=0;j<_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;j++) {
      fread(&elementCountingNumber,sizeof(long int),1,fout);
      elementStackList[i][j]=GetEntryPointer(elementCountingNumber);
    } 
  }
   

  
  void init() {
    clear();
    initMemoryBlock();
  }
    

  void explicitConstructor() {
    MemoryAllocation=0;

    nMaxElements=0,elementStackList=NULL,elementStackPointer=0;
    dataBufferList=0,dataBufferList=NULL,dataBufferListSize=0,dataBufferListPointer=0;

  }

  cStack() {
    explicitConstructor();
  }

  long int capasity() {return nMaxElements-elementStackPointer;}
  long int totalSize() {return nMaxElements;}
  long int usedElements() {return elementStackPointer;}


  T* newElement() {
    T* res;

    if (sizeof(T)==0) return NULL;
    if (elementStackPointer==nMaxElements) initMemoryBlock(); 

    long int elementStackBank,offset;

    elementStackBank=elementStackPointer/_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=elementStackPointer-elementStackBank*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    res=elementStackList[elementStackBank][offset];
    elementStackPointer++;

    return res;
  }

  void deleteElement(T* delElement) {
    if (sizeof(T)==0) return;

    if (elementStackPointer==0) {
      cout << "$PREFIX:ERROR: stack pointer is 0 (line=" << __LINE__ << ", file=" << __FILE__ << ")" << endl;
    } 

    long int elementStackBank,offset;
    --elementStackPointer;

    elementStackBank=elementStackPointer/_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=elementStackPointer-elementStackBank*_GENERAL_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    elementStackList[elementStackBank][offset]=delElement;
  } 
};


#endif
