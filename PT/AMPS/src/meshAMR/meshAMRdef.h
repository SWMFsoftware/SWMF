
//$Id$
//global macroscopic variables for the AMR mesh (no cut cells) 


#include "mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>

#include "global.h"
#include "constants.h"

#ifndef _AMR_MESH_DEFINITION_
#define _AMR_MESH_DEFINITION_

/*
#define Pi 3.14159265358979323846264338327950288419716939937510582
#define sqrtPi 1.7724538509055160272981674833411
#define PiTimes2 6.28318530717958647692528676655900576839433879875021
*/

#define _ON_AMR_MESH_    1
#define _OFF_AMR_MESH_   0


#define _BLOCK_CELLS_X_    5
#define _BLOCK_CELLS_Y_    5
#define _BLOCK_CELLS_Z_    5

#define _GHOST_CELLS_X_       2
#define _GHOST_CELLS_Y_       2
#define _GHOST_CELLS_Z_       2

#define _TOTAL_BLOCK_CELLS_X_ (_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)
#define _TOTAL_BLOCK_CELLS_Y_ (_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)
#define _TOTAL_BLOCK_CELLS_Z_ (_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)

//use cut-cells in the mesh
#define _AMR__CUT_CELL__MODE__ON_      0
#define _AMR__CUT_CELL__MODE__OFF_     1

#define _AMR__CUT_CELL__MODE_  _AMR__CUT_CELL__MODE__ON_


//the existance of the internal boundaries (bodies inside the compurational domain)
#define _INTERNAL_BOUNDARY_MODE_ON_  0
#define _INTERNAL_BOUNDARY_MODE_OFF_ 1

#define _INTERNAL_BOUNDARY_MODE_ _INTERNAL_BOUNDARY_MODE_ON_


//allow user to add data to the definitions of the internal boundaries
#define _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_     0
#define _USED_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_OFF_    1

#define _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_

//load the user difinitions
#define _AMR__LOAD_USER_DEFINITION__MODE__ON_      0
#define _AMR__LOAD_USER_DEFINITION__MODE__OFF_     1

#define _AMR__LOAD_USER_DEFINITION__MODE_  _AMR__LOAD_USER_DEFINITION__MODE__OFF_


//the types of the internal boundaries
#define _INTERNAL_BOUNDARY_TYPE_UNDEFINED_         0
#define _INTERNAL_BOUNDARY_TYPE_SPHERE_            1
#define _INTERNAL_BOUNDARY_TYPE_CIRCLE_            2
#define _INTERNAL_BOUNDARY_TYPE_1D_SPHERE_         3
#define _INTERNAL_BOUNDARY_TYPE_BODY_OF_ROTATION_  4

//the types of symmetries used with the mesh
#define _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_    0
#define _AMR_SYMMETRY_MODE_AXIAL_SYMMETRY_     1
#define _AMR_SYMMETRY_MODE_SPHERICAL_SYMMETRY_ 2

#define _AMR_SYMMETRY_MODE_ _AMR_SYMMETRY_MODE_PLANAR_SYMMETRY_

//Zenith angle distribution on the internal boundary sphere
#define _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_  0
#define _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_   1

#define _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_ _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_
//#define _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_


//the maximum number and othe length of the variable that containes the number of a mesh node's connection (the number of how many blocks containes the corner node)
#define _MAX_CORNER_NODE_CONNECTION_       15
#define _MAX_CORNER_NODE_CONNECTION_BITS_   6 

//the maximum value of the size of the variable that contains the block's refinment level
#define _MAX_REFINMENT_LEVEL_BITS_ 5
#define _MAX_REFINMENT_LEVEL_  15


//save additional information when the mesh generator is used in the debugger mode
#define _AMR_DEBUGGER_MODE_ON_    1
#define _AMR_DEBUGGER_MODE_OFF_   0

#define _AMR_DEBUGGER_MODE_ _AMR_DEBUGGER_MODE_ON_ 


//check the consistancy of the mesh on each level of refinments during creation of the mesh
#define _CHECK_MESH_CONSISTANCY_ _ON_AMR_MESH_


//the number of bits reserved to store the counting number of the mesh elements
#define _MESH_ELEMENTS_NUMBERING_BITS_   24
#define _MAX_MESH_ELEMENT_NUMBER_        (-1+(1<<_MESH_ELEMENTS_NUMBERING_BITS_)) 


//define the internal structure of the block 
#define _AMR_CENTER_NODE_  _ON_AMR_MESH_
#define _AMR_FACE_NODE_    _OFF_AMR_MESH_
#define _AMR_EDGE_NODE_    _OFF_AMR_MESH_
 

//definition for the default face attributes
#define _INTERNAL_FACE_NO_NEIBOUR_FACEAT_      -2
#define _ROOT_BLOCK_EXTERNAL_BOUNDARY_FACEAT_  -1

//the location of of the block regarding to the boundary of the mesh
#define _AMR_BLOCK_INSIDE_DOMAIN_              0
#define _AMR_BLOCK_OUTSIDE_DOMAIN_             1
#define _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_  2

//keep the blocks that are entirely covered by the internal surfaces installed into the mesh
#define _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_REMOVE_ 0
#define _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_KEEP_   1

#define _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_ _AMR_BLOCK_OUTSIDE_DOMAIN_MODE_KEEP_

//the parallel mode used on the mesh
#define _AMR_PARALLEL_MODE_ON_    1
#define _AMR_PARALLEL_MODE_OFF_   0

#define _AMR_PARALLEL_MODE_ _AMR_PARALLEL_MODE_ON_


//the type of the data exchange in the mesh's parallel mode : GhostCells -> only that data, which is stored in ghost cells is transfered; BlockBoundaryLayer -> a layer of blocks is created on the domain's boundary, and the whole blocks is transfered
#define _AMR_PARALLEL_DATA_EXCHANGE_MODE__GHOST_CELLS_           0
#define _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_ 1

#define _AMR_PARALLEL_DATA_EXCHANGE_MODE_ _AMR_PARALLEL_DATA_EXCHANGE_MODE__DOMAIN_BOUNDARY_LAYER_


//the length of symbolic strings used by the mesh generator
#define STRING_LENGTH 200


//a stack class for the mesh elements
#define _STACK_DEFAULT_BUFFER_BUNK_SIZE_ 10000
#define _STACK_DEFAULT_BUFFER_LIST_SIZE_ 10000

//macro definition for the real and ghost blocks 
#define _GHOST_BLOCK_ true
#define _REAL_BLOCK_  false 

//the mode of output of the connectivity list:
//_AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_: each block is printed into a file independently from other blocks -> nodes with the same coordinate can appear in the file
//_AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_: global node numbering is used -> no nodes with the same locations in the output file
#define _AMR_CONNECTIVITY_LIST_MODE_INTERNAL_BLOCK_    0
#define _AMR_CONNECTIVITY_LIST_MODE_GLOBAL_NODE_NUMBER_ 1

using namespace std;

//the exist function with printing of the error message
class cAMRexit {
public:
void exit(const long int nline, const char* fname,const char* msg=NULL) {
  char str[1000];
  int mpiInitFlag,ThisThread;


  if (msg==NULL) sprintf(str," exit: line=%ld, file=%s\n",nline,fname);
  else sprintf(str," exit: line=%ld, file=%s, message=%s\n",nline,fname,msg);

  FILE* errorlog=fopen("$ERROR","a+");

  time_t TimeValue=time(0);
  tm *ct=localtime(&TimeValue);

  MPI_Initialized(&mpiInitFlag);

  if (mpiInitFlag==true)  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  else ThisThread=0;

  fprintf(errorlog,"Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  fprintf(errorlog,"file=%s, line=%ld\n",fname,nline);
  fprintf(errorlog,"%s\n\n",msg);

  printf("$PREFIX:Thread=%i: (%i/%i %i:%i:%i)\n",ThisThread,ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);
  printf("$PREFIX:file=%s, line=%ld\n",fname,nline);
  printf("$PREFIX:%s\n\n",msg);

  fclose(errorlog);
  ::exit(1);
}

virtual ~cAMRexit() { }
};

//the stack class to store the data structure of the mesh
template<class T>
class cAMRstack : public cAMRexit {
public: 
  //data element's stack
  long int nMaxElements; 
  T*** elementStackList;
  long int elementStackPointer;

  //array pointers on the allocates data blocks
  T** dataBufferList;
  long int dataBufferListSize,dataBufferListPointer;  


  //the element temporary ID (used only in the debugger mode)
  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID_counter;
  #endif

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
      T** tmpDataList=new T*[dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_];
      MemoryAllocation+=sizeof(T*)*(dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_);

      T*** tmpStackList=new T**[dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_];
      MemoryAllocation+=sizeof(T**)*(dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_); 

      i=0;
     
      //copy the content of the old lists to the new ones 
      if (dataBufferList!=NULL) for (;i<dataBufferListSize;i++) tmpDataList[i]=dataBufferList[i],tmpStackList[i]=elementStackList[i];
      for (j=0;j<_STACK_DEFAULT_BUFFER_LIST_SIZE_;j++,i++) tmpDataList[i]=NULL,tmpStackList[i]=NULL;

      if (dataBufferList!=NULL) {
        delete [] dataBufferList;
        MemoryAllocation-=sizeof(T*)*dataBufferListSize;

        delete [] elementStackList;
        MemoryAllocation-=sizeof(T**)*dataBufferListSize;
      }

      dataBufferList=tmpDataList;
      elementStackList=tmpStackList;
      dataBufferListSize+=_STACK_DEFAULT_BUFFER_LIST_SIZE_;
    }

    //allocate a new memory chunk for the element's data and update the stack list
    dataBufferList[dataBufferListPointer]=new T[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
    elementStackList[dataBufferListPointer]=new T*[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];

    MemoryAllocation+=sizeof(T)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    MemoryAllocation+=sizeof(T*)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    for (i=0;i<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;i++) elementStackList[dataBufferListPointer][i]=dataBufferList[dataBufferListPointer]+i; 
 
    nMaxElements+=_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    dataBufferListPointer++;
  }

  void resetStack() {
    long int databank,offset;

    for (databank=0;databank<dataBufferListPointer;databank++) for (offset=0;offset<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;offset++) elementStackList[databank][offset]=dataBufferList[databank]+offset; 

    elementStackPointer=0;
  }
 


  //get the entry pointer and counting number
  long int GetEntryCountingNumber(T* ptr) {
    long int nMemoryBank,res=-1;

    if (ptr!=NULL) {
      for (nMemoryBank=0;nMemoryBank<dataBufferListPointer;nMemoryBank++) if ((ptr>=dataBufferList[nMemoryBank])&&(ptr<dataBufferList[nMemoryBank]+_STACK_DEFAULT_BUFFER_BUNK_SIZE_)) if ((res=ptr-dataBufferList[nMemoryBank])<_STACK_DEFAULT_BUFFER_BUNK_SIZE_) {
        return res+nMemoryBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
      }
 
      exit(__LINE__,__FILE__,"Error: canot find the entry number");  
    }
 
    return -1;
  }

  T* GetEntryPointer(long int countingNumber) {
    long int nMemoryBank,offset;

    if ((countingNumber<0.0)||(countingNumber>=nMaxElements)) return NULL;

    nMemoryBank=countingNumber/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=countingNumber-nMemoryBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    return dataBufferList[nMemoryBank]+offset;
  }

  void clear() {
    for (int i=0;i<dataBufferListPointer;i++) {
      delete [] dataBufferList[i];
      delete [] elementStackList[i];

      MemoryAllocation-=(sizeof(T)+sizeof(T*))*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
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

    for (i=0;i<dataBufferListPointer;i++) for (j=0;j<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;j++) {    
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
      dataBufferList[i]=new T[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
      MemoryAllocation+=sizeof(T)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

      elementStackList[i]=new T*[_STACK_DEFAULT_BUFFER_BUNK_SIZE_]; 
      MemoryAllocation+=sizeof(T*)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    }

    //read the elementStack
    for (i=0;i<dataBufferListPointer;i++) for (j=0;j<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;j++) {
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

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID_counter=0;
    #endif

  }

  cAMRstack() {
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

    elementStackBank=elementStackPointer/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=elementStackPointer-elementStackBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    res=elementStackList[elementStackBank][offset];
    elementStackPointer++;

    if (usedElements()>_MAX_MESH_ELEMENT_NUMBER_) exit(__LINE__,__FILE__,"The number of the requster mesh elements exeeds the limit -> increase _MESH_ELEMENTS_NUMBERING_BITS_ "); 


    res->cleanDataBuffer();

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    res->Temp_ID=Temp_ID_counter++;
    #endif
    
    return res;
  }

  void deleteElement(T* delElement) {
    if (sizeof(T)==0) return;

    delElement->cleanDataBuffer();

    if (elementStackPointer==0) {
      cout << "$PREFIX:ERROR: stack pointer is 0 (line=" << __LINE__ << ", file=" << __FILE__ << ")" << endl;
    } 

    long int elementStackBank,offset;
    --elementStackPointer;

    elementStackBank=elementStackPointer/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=elementStackPointer-elementStackBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    elementStackList[elementStackBank][offset]=delElement;
  } 
};


//=============================================================
//the stack class that is capable to store additional data associated with the objects stored in the stack
template<class T>
class cAssociatedDataAMRstack : public cAMRstack <T> {
public:
  //the stack's structure to store the associated data
  char*** associatedDataStackList;
  char** associatedDataBufferList;

  long int getAllocatedMemory() {
    T t;
	return cAMRstack <T>::dataBufferListSize*(sizeof(T)+sizeof(T*)+sizeof(char)*t.AssociatedDataLength());
  }

  void initMemoryBlock() {
	T t;

    if (t.AssociatedDataLength()!=0) {
      long int i=0,j=0;

	  //check available space in the dataBufferList list: if needed increment the size of 'elementStackList' and 'dataBufferList'
	  if (cAMRstack <T>::dataBufferListPointer==cAMRstack <T>::dataBufferListSize) {
	    char** tmpDataList=new char*[cAMRstack <T>::dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_];
	    char*** tmpStackList=new char**[cAMRstack <T>::dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_];

	    cAMRstack <T>::MemoryAllocation+=sizeof(char*)*(cAMRstack <T>::dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_);
	    cAMRstack <T>::MemoryAllocation+=sizeof(char**)*(cAMRstack <T>::dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_);

	    //copy the content of the old lists to the new ones
	    if (associatedDataBufferList!=NULL) for (;i<cAMRstack <T>::dataBufferListSize;i++) tmpDataList[i]=associatedDataBufferList[i],tmpStackList[i]=associatedDataStackList[i];
	    for (j=0;j<_STACK_DEFAULT_BUFFER_LIST_SIZE_;j++,i++) tmpDataList[i]=NULL,tmpStackList[i]=NULL;

	    if (associatedDataBufferList!=NULL) {
	      delete [] associatedDataBufferList;
	      delete [] associatedDataStackList;

	      cAMRstack <T>::MemoryAllocation-=(sizeof(char*)+sizeof(char**))*cAMRstack <T>::dataBufferListSize;
	    }

	    associatedDataBufferList=tmpDataList;
	    associatedDataStackList=tmpStackList;
	  }

	  //allocate a new memory chunk for the element's data and update the stack list
	  long int offset=t.AssociatedDataLength();

	  associatedDataBufferList[cAMRstack <T>::dataBufferListPointer]=new char[_STACK_DEFAULT_BUFFER_BUNK_SIZE_*offset];
	  associatedDataStackList[cAMRstack <T>::dataBufferListPointer]=new char*[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];

	  cAMRstack <T>::MemoryAllocation+=(offset*sizeof(char)+sizeof(char*))*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

	  for (i=0;i<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;i++) associatedDataStackList[cAMRstack <T>::dataBufferListPointer][i]=associatedDataBufferList[cAMRstack <T>::dataBufferListPointer]+i*offset;
	}

   //init the buffer for the stack object itself
   cAMRstack <T>::initMemoryBlock() ;
  }



  T* newElement() {
    T* res;

    if (sizeof(T)==0) return NULL;
    if (cAMRstack <T>::elementStackPointer==cAMRstack <T>::nMaxElements) initMemoryBlock();

    long int elementStackBank,offset;

    elementStackBank=cAMRstack <T>::elementStackPointer/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=cAMRstack <T>::elementStackPointer-elementStackBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    res=cAMRstack <T>::elementStackList[elementStackBank][offset];
    if (associatedDataStackList!=NULL) res->SetAssociatedDataBufferPointer(associatedDataStackList[elementStackBank][offset]);
    cAMRstack <T>::elementStackPointer++;

    if (cAMRstack <T>::usedElements()>_MAX_MESH_ELEMENT_NUMBER_) cAMRexit::exit(__LINE__,__FILE__,"The number of the requested mesh elements exeeds the limit -> increase _MESH_ELEMENTS_NUMBERING_BITS_ ");

    res->cleanDataBuffer();

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    res->Temp_ID=cAMRstack <T>::Temp_ID_counter++;
    #endif




//================   DEBUG =========================
/*

    int t;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&t);

    if (t==3) if ((res->Temp_ID==2708279)||(res->Temp_ID==536172)) {
      cout << __FILE__ << __LINE__ << endl;
    }
*/
//================ END DEBUG =======================


    return res;
  }


  void clear() {
    T t;
    long int offset=t.AssociatedDataLength();

    if (associatedDataBufferList!=NULL) {
      for (int i=0;i<cAMRstack <T>::dataBufferListPointer;i++) {
        delete [] associatedDataBufferList[i];
        delete [] associatedDataStackList[i];

        cAMRstack <T>::MemoryAllocation-=(sizeof(char)*offset+sizeof(char*))*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
      }

      delete [] associatedDataBufferList;
      delete [] associatedDataStackList;

      cAMRstack <T>::MemoryAllocation-=(sizeof(char*)+sizeof(char**))*cAMRstack <T>::dataBufferListSize;
    }

    associatedDataBufferList=NULL,associatedDataStackList=NULL;

    cAMRstack <T>::clear();
  }

  void init() {
    clear();
    initMemoryBlock();

    cAMRstack <T>::init();
  }

  void explicitConstructor() {
	associatedDataStackList=NULL;
	associatedDataBufferList=NULL;

	cAMRstack <T>::explicitConstructor();
  }


   cAssociatedDataAMRstack() : cAMRstack<T>() {
	 explicitConstructor();
   }

   void deleteElement(T* delElement) {


	 if (delElement->AssociatedDataLength()!=0) {
       long int elementStackBank,offset;
       long int localElementStackPointer=cAMRstack<T>::elementStackPointer-1;

       elementStackBank=localElementStackPointer/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
       offset=localElementStackPointer-elementStackBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

       associatedDataStackList[elementStackBank][offset]=delElement->GetAssociatedDataBufferPointer();
	 }




	 cAMRstack<T>::deleteElement(delElement);
	 delElement->SetAssociatedDataBufferPointer(NULL);
   }

   //save and load the allocation of the stack
   void saveAllocationParameters(FILE *fout) {
     cAMRstack<T>::saveAllocationParameters(fout);
   }

   void readAllocationParameters(FILE *fout) {
     cAMRstack<T>::readAllocationParameters(fout);

     if (cAMRstack<T>::dataBufferListPointer!=0) cAMRstack<T>::exit(__LINE__,__FILE__,"not implemented");
   }

};

//the heap class to store the data structure of the mesh: the data are stored linearly, no data can be deleter from the heap
template<class T>
class cAMRheap : public cAMRexit {
public:
  //data element's stack
  long int nMaxElements;
  long int elementHeapPointer;

  //array pointers on the allocates data blocks
  T** dataBufferList;
  long int dataBufferListSize,dataBufferListPointer;


  //the element temporary ID (used only in the debugger mode)
  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID_counter;
  #endif

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
      T** tmpDataList=new T*[dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_];
      MemoryAllocation+=sizeof(T*)*(dataBufferListSize+_STACK_DEFAULT_BUFFER_LIST_SIZE_);

      i=0;

      //copy the content of the old lists to the new ones
      if (dataBufferList!=NULL) for (;i<dataBufferListSize;i++) tmpDataList[i]=dataBufferList[i];
      for (j=0;j<_STACK_DEFAULT_BUFFER_LIST_SIZE_;j++,i++) tmpDataList[i]=NULL;

      if (dataBufferList!=NULL) {
        delete [] dataBufferList;
        MemoryAllocation-=sizeof(T*)*dataBufferListSize;
      }

      dataBufferList=tmpDataList;
      dataBufferListSize+=_STACK_DEFAULT_BUFFER_LIST_SIZE_;
    }

    //allocate a new memory chunk for the element's data and update the stack list
    dataBufferList[dataBufferListPointer]=new T[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
    MemoryAllocation+=sizeof(T)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;


    nMaxElements+=_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    dataBufferListPointer++;
  }





  //get the entry pointer and counting number
  long int GetEntryCountingNumber(T* ptr) {
    long int nMemoryBank,res=-1;

    if (ptr!=NULL) {
      for (nMemoryBank=0;nMemoryBank<dataBufferListPointer;nMemoryBank++) if ((ptr>=dataBufferList[nMemoryBank])&&(ptr<dataBufferList[nMemoryBank]+_STACK_DEFAULT_BUFFER_BUNK_SIZE_)) if ((res=ptr-dataBufferList[nMemoryBank])<_STACK_DEFAULT_BUFFER_BUNK_SIZE_) {
        return res+nMemoryBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
      }

      exit(__LINE__,__FILE__,"Error: canot find the entry number");
    }

    return -1;
  }

  T* GetEntryPointer(long int countingNumber) {
    long int nMemoryBank,offset;

    if ((countingNumber<0.0)||(countingNumber>=nMaxElements)) return NULL;

    nMemoryBank=countingNumber/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=countingNumber-nMemoryBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    return dataBufferList[nMemoryBank]+offset;
  }

  void clear() {
    for (int i=0;i<dataBufferListPointer;i++) {
      delete [] dataBufferList[i];

      MemoryAllocation-=sizeof(T)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    }

    if (dataBufferList!=NULL) {
      delete [] dataBufferList;

      MemoryAllocation-=sizeof(T*)*dataBufferListSize;
    }

    nMaxElements=0,elementHeapPointer=0;
    dataBufferList=0,dataBufferList=NULL,dataBufferListSize=0,dataBufferListPointer=0;
  }


  //save and load the allocation of the stack
  void saveAllocationParameters(FILE *fout) {
    fwrite(&nMaxElements,sizeof(long int),1,fout);
    fwrite(&elementHeapPointer,sizeof(long int),1,fout);
    fwrite(&dataBufferListSize,sizeof(long int),1,fout);
    fwrite(&dataBufferListPointer,sizeof(long int),1,fout);


    exit(__LINE__,__FILE__,"Check the implementetion!");

  }


  void readAllocationParameters(FILE *fout) {
    clear();

    fread(&nMaxElements,sizeof(long int),1,fout);
    fread(&elementHeapPointer,sizeof(long int),1,fout);
    fread(&dataBufferListSize,sizeof(long int),1,fout);
    fread(&dataBufferListPointer,sizeof(long int),1,fout);

    //allocate the stack's buffers
    long int i;

    dataBufferList=new T*[dataBufferListSize];
    MemoryAllocation+=sizeof(T*)*dataBufferListSize;
    for (i=0;i<dataBufferListSize;i++) dataBufferList[i]=NULL;

    for (i=0;i<dataBufferListPointer;i++) {
      dataBufferList[i]=new T[_STACK_DEFAULT_BUFFER_BUNK_SIZE_];
      MemoryAllocation+=sizeof(T)*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    }

    exit(__LINE__,__FILE__,"Check the implementation");
  }





  void init() {
    clear();
    initMemoryBlock();
  }


  void explicitConstructor() {
    MemoryAllocation=0;

    nMaxElements=0,elementHeapPointer=0;
    dataBufferList=0,dataBufferList=NULL,dataBufferListSize=0,dataBufferListPointer=0;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID_counter=0;
    #endif

  }

  cAMRheap() {
    explicitConstructor();
  }

  long int capasity() {return nMaxElements-elementHeapPointer;}
  long int totalSize() {return nMaxElements;}
  long int usedElements() {return elementHeapPointer;}


  T* newElement() {
    T* res;

    if (sizeof(T)==0) return NULL;
    if (elementHeapPointer==nMaxElements) initMemoryBlock();

    long int elementHeapBank,offset;

    elementHeapBank=elementHeapPointer/_STACK_DEFAULT_BUFFER_BUNK_SIZE_;
    offset=elementHeapPointer-elementHeapBank*_STACK_DEFAULT_BUFFER_BUNK_SIZE_;

    res=dataBufferList[elementHeapBank]+offset;
    elementHeapPointer++;

    if (usedElements()>_MAX_MESH_ELEMENT_NUMBER_) exit(__LINE__,__FILE__,"The number of the requested mesh elements exceeds the limit -> increase _MESH_ELEMENTS_NUMBERING_BITS_ ");


    res->cleanDataBuffer();

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    res->Temp_ID=Temp_ID_counter++;
    #endif

    return res;
  }

};



#endif
