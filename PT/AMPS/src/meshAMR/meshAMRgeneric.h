
//$Id$
//the AMR mesh in 1,2 and 3D without cut cells

#ifndef _AMR_MESH_GENERIC_ 
#define _AMR_MESH_GENERIC_ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>

#include <sys/time.h>
#include <sys/resource.h>


#include <string>
#include <iterator>
#include <utility>
#include <algorithm>
#include <list>
#include <vector>

#include "mpi.h"


#include "meshAMRdef.h"
#include "meshAMRinternalSurface.h"

#include "specfunc.h"


//define the boolian macro variables used in the mesh 
#define _FALSE_ 0
#define _TRUE_ 1



//the limits of the comlutational domain
extern double _MESH_AMR_XMAX_[3],_MESH_AMR_XMIN_[3];

class cBasicNode : public cAMRexit {
public:
  //the place holder for the structure that contained the associated data
  int AssociatedDataLength() {return 0;}
  void SetAssociatedDataBufferPointer(char* ptr) {}
  char* GetAssociatedDataBufferPointer() {return NULL;}

  //placeholder for the printing procedured
  void PrintData(FILE* fout,int DataSetNumber) {}
  void PrintVariableList(FILE* fout,int DataSetNumber) {}
  void PrintFileDescriptior(FILE* fout,int DataSetNumber) {}

protected:
  #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
  double x[_MESH_DIMENSION_];
  #endif

public:
  #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  struct cNodeDescriptor {
    unsigned long int nodeno : _MESH_ELEMENTS_NUMBERING_BITS_;
    unsigned nodeProcessedFlag : 1;
    unsigned maxRefinmentLevel : _MAX_REFINMENT_LEVEL_BITS_;
    unsigned internalMeshNode : 1;
    unsigned int nNodeConnections : _MAX_CORNER_NODE_CONNECTION_BITS_;
  } nodeDescriptor;


  //get and set the nodes' positions
  inline double *GetX() {
    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    return x;
  }

  inline void GetX(double *l) {
    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    for (int idim=0;idim<_MESH_DIMENSION_;idim++) l[idim]=x[idim];
  }

  inline void SetX(double *l) {


//########  DEBUG ##########
	  if (Temp_ID==88861) {
		  cout << __LINE__ << endl;
	  }

//######## END DEBUG #######


    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    exit(__LINE__,__FILE__,"The operation is allowed only in the debugger mode");
    #endif

    for (int idim=0;idim<_MESH_DIMENSION_;idim++) x[idim]=l[idim];
  }


  //clean the data buffers
  void cleanDataBuffer() {
    nodeDescriptor.nodeno=0;
    nodeDescriptor.nodeProcessedFlag=_FALSE_;
    nodeDescriptor.maxRefinmentLevel=0;
    nodeDescriptor.internalMeshNode=_TRUE_;
    nodeDescriptor.nNodeConnections=0;

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  //increment and decrement the node conenction's counter
  void incrementConnectionCounter() {
    if (nodeDescriptor.nNodeConnections==_MAX_CORNER_NODE_CONNECTION_) exit(__LINE__,__FILE__,"the node's connections exeeds _MAX_CORNER_NODE_CONNECTION_");

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    nodeDescriptor.nNodeConnections++;
    #endif
  }

  int getNodeConnectionNumber() {

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    return nodeDescriptor.nNodeConnections;
    #else
    return 1;
    #endif
  }

  int decrementConnectionCounter() {

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    if (nodeDescriptor.nNodeConnections==0) exit(__LINE__,__FILE__,"the node's connections becomes negative");
    nodeDescriptor.nNodeConnections--;
    #endif

    return getNodeConnectionNumber();
  }

}; 

class cBasicCornerNode : public cBasicNode {
public: 
  //double x[_MESH_DIMENSION_];
};


class cBasicCenterNode : public cBasicNode {
public: 


};



template <class cBlockAMR> class cTreeNodeAMR;
template <class cCornerNode,class cCenterNode> class cBasicBlockAMR; 

//===================================================================

#define _BOTTOM_BRANCH_TREE_ false

#define _GLOBAL_POSITION_REAL_NODE_UNKNOWN_ 0
#define _GLOBAL_POSITION_REAL_NODE_UP_      1
#define _GLOBAL_POSITION_REAL_NODE_CURRENT_ 2
#define _GLOBAL_POSITION_REAL_NODE_DOWN_    3

template <class cBlockAMR>
class cTreeNodeAMR : public cAMRexit {
public:
  cTreeNodeAMR *upNode,*downNode[1<<_MESH_DIMENSION_],*neibNode[2*_MESH_DIMENSION_];
  cBlockAMR *block;

  #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  //the pointers to members of the list of the nodes that "belongs" to this processor (or some other lists)
  cTreeNodeAMR *nextNode,*prevNode;

  double xmin[_MESH_DIMENSION_],xmax[_MESH_DIMENSION_];
  int RefinmentLevel;

  //the list of the descriptors of the internal boundaries installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  cInternalBoundaryConditionsDescriptor *InternalBoundaryDescriptorList;
  #endif

  //get the external (directed outward of the computational domain) normal
  void GetExternalNormal(double *norm,int nface) {
    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");

    #elif _MESH_DIMENSION_ == 3
    switch(nface) {
    case 0: case 1:
      norm[0]=(nface==0) ? -1.0 : 1.0;
      norm[1]=0.0,norm[2]=0.0;
      break;
    case 2: case 3:
      norm[2]=(nface==2) ? -1.0 : 1.0;
      norm[0]=0.0,norm[2]=0.0;
      break;
    case 4: case 5:
      norm[3]=(nface==4) ? -1.0 : 1.0;
      norm[0]=0.0,norm[1]=0.0;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: 'nface' is out of the range");
    }
    #endif

  }

  double GetCharacteristicCellSize() {
    double CellSize;

    CellSize=pow((xmax[0]-xmin[0])/_BLOCK_CELLS_X_,2);
    if (_MESH_DIMENSION_>1) CellSize+=pow((xmax[1]-xmin[1])/_BLOCK_CELLS_Y_,2);
    if (_MESH_DIMENSION_>2) CellSize+=pow((xmax[2]-xmin[2])/_BLOCK_CELLS_Z_,2);

    return sqrt(CellSize);
  }

  double GetBlockFaceSurfaceArea(int nface) {
    double res;

    if (nface>=2*_MESH_DIMENSION_) exit(__LINE__,__FILE__,"Error: 'nface' is out of range");

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3

    const static int faceNodeCoordinateFrame[6][3][3]={ {{0,0,0},{0,1,0},{0,0,1}},{{1,0,0},{1,1,0},{1,0,1}},
                                                        {{0,0,0},{1,0,0},{0,0,1}},{{0,1,0},{1,1,0},{0,1,1}},
                                                        {{0,0,0},{1,0,0},{0,1,0}},{{0,0,1},{1,0,1},{0,1,1}} };



    double x0[3],x1[3],x2[3],l0,l1;
    int idim;

    x0[0]=(faceNodeCoordinateFrame[nface][0][0]==0) ? xmin[0] : xmax[0];
    x0[1]=(faceNodeCoordinateFrame[nface][0][1]==0) ? xmin[1] : xmax[1];
    x0[2]=(faceNodeCoordinateFrame[nface][0][2]==0) ? xmin[2] : xmax[2];

    x1[0]=(faceNodeCoordinateFrame[nface][1][0]==0) ? xmin[0] : xmax[0];
    x1[1]=(faceNodeCoordinateFrame[nface][1][1]==0) ? xmin[1] : xmax[1];
    x1[2]=(faceNodeCoordinateFrame[nface][1][2]==0) ? xmin[2] : xmax[2];

    x2[0]=(faceNodeCoordinateFrame[nface][2][0]==0) ? xmin[0] : xmax[0];
    x2[1]=(faceNodeCoordinateFrame[nface][2][1]==0) ? xmin[1] : xmax[1];
    x2[2]=(faceNodeCoordinateFrame[nface][2][2]==0) ? xmin[2] : xmax[2];

    for (idim=0,l0=0.0,l1=0.0;idim<_MESH_DIMENSION_;idim++) l0+=pow(x1[idim]-x0[idim],2),l1+=pow(x2[idim]-x0[idim],2);
    res=sqrt(l0*l1);

    #endif

    return res;
  }

  double GetCellFaceSurfaceArea(int nface) {
    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3

    const static double CellSurfaceAreaMultiplyer[6]={_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_,
                                                      _BLOCK_CELLS_X_*_BLOCK_CELLS_Z_,_BLOCK_CELLS_X_*_BLOCK_CELLS_Z_,
                                                      _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_,_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_};
    #endif

    return GetBlockFaceSurfaceArea(nface)/CellSurfaceAreaMultiplyer[nface];
  }

  struct cTreeNodeDescriptor {
    unsigned int GlobalPositionRealTreeNode : 2;
  } treeNodeDescriptor; 


  void cleanDataBuffer() {
    int i;

    block=NULL,upNode=NULL;

    for (i=0;i<(1<<_MESH_DIMENSION_);i++) downNode[i]=NULL;
    for (i=0;i<2*_MESH_DIMENSION_;i++) neibNode[i]=NULL;
  
    treeNodeDescriptor.GlobalPositionRealTreeNode=_GLOBAL_POSITION_REAL_NODE_UNKNOWN_; 

    RefinmentLevel=-1;
    nextNode=NULL,prevNode=NULL;

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif

    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    InternalBoundaryDescriptorList=NULL;
    #endif
  }

  //find the neibour of the tree node: this version of the function is used only in 1D case
  cTreeNodeAMR *GetNeib(int i) {
    cTreeNodeAMR *res=NULL;

    #if _MESH_DIMENSION_ == 1 
    res=neibNode[i];
    #else
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif
  }

  cTreeNodeAMR *GetNeib(int i,int j,int k) {
    #if _MESH_DIMENSION_ == 1 
    return GetNeib(i);
    #else 
    exit(__LINE__,__FILE__,"not implemented yet");
    #endif
  }

  bool allocated() {return (block==NULL) ? false : true;};

  int GetGlobalPositinoRealNode() {return treeNodeDescriptor.GlobalPositionRealTreeNode;}
  void SettGlobalPositinoRealNode(int RealNodePosition) {treeNodeDescriptor.GlobalPositionRealTreeNode=RealNodePosition;} 


  cTreeNodeAMR *GetDownNode(int i,int j,int k) {return downNode[i+2*(j+k*2)];} 
  void SetDownNode(cTreeNodeAMR* node,int i,int j,int k) {downNode[i+2*(j+k*2)]=node;}  

  bool lastBranchFlag() {return (downNode[0]==NULL) ? _BOTTOM_BRANCH_TREE_ : !_BOTTOM_BRANCH_TREE_;} //if _BOTTOM_BRANCH_TREE_ -> the node is on the bottom of the tree   
};





//=======================================================================
template <class cCornerNode,class cCenterNode> 
class cBasicBlockAMR : public cAMRexit {
public:
  //the place holder for the structure that contained the associated data
  int AssociatedDataLength() {return 0;}
  void SetAssociatedDataBufferPointer(char* ptr) {}
  char* GetAssociatedDataBufferPointer() {return NULL;}

  //place holder for print function that outputs the general tree node/blocks values into a file
  void PrintData(FILE* fout,int DataSetNumber) {}
  void PrintVariableList(FILE* fout) {}

protected:
  //nodes of the block

  #if _MESH_DIMENSION_ == 1
  cCornerNode *cornerNodes[1+_TOTAL_BLOCK_CELLS_X_];

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  exit(__LINE__,__FILE__,"not implemented yet");
  #endif

  #elif _MESH_DIMENSION_ == 2
  cCornerNode *cornerNodes[(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)];

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  exit(__LINE__,__FILE__,"not implemented yet");
  #endif

  #elif _MESH_DIMENSION_ == 3
  cCornerNode *cornerNodes[(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_)]; 

  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *centerNodes[_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_];
  #endif


  #endif

public:

  //set and get the pointers to the center nodes of the block
  cCenterNode **GetCenterNodeBuffer() {
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    return centerNodes;
    #else
    return NULL;
    #endif
  }

  inline cCenterNode *GetCenterNode(long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;
    #endif

    if ((nd<0)||(nd>=nMaxCenterNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    return centerNodes[nd];
    #else
    return NULL;
    #endif
  }


  inline void SetCenterNode(cCenterNode* nodeptr,long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_;
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCenterNodes=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_;
    #endif

    if ((nd<0)||(nd>=nMaxCenterNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    centerNodes[nd]=nodeptr;
    #endif
  }



  //set and get the pointers to the corner nodes of the block
  cCornerNode **GetCornerNodeBuffer() {
	 return cornerNodes;
  }

  inline cCornerNode *GetCornerNode(long int nd) {

    #if _MESH_DIMENSION_ == 1
    static const int nMaxCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    if ((nd<0)||(nd>=nMaxCornerNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    return cornerNodes[nd];
  }


  inline void SetCornerNode(cCornerNode* nodeptr,long int nd) {
    #if _MESH_DIMENSION_ == 1
    static const int nMaxCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    static const int nMaxCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    if ((nd<0)||(nd>=nMaxCornerNodes)) exit(__LINE__,__FILE__,"The value is outside of the limit");

    cornerNodes[nd]=nodeptr;
  }


  #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif

  void *nextBlock; //the next block in the list of allocated blocks of the mesh


  struct cBlockDescriptor {
    unsigned RefinmentLevel : _MAX_REFINMENT_LEVEL_BITS_;
    unsigned ghostBlock : 1;
  } blockDescriptor; 

  double CellMeasure(int i,int j, int k) {
    register double vol;

    vol=(_MESH_AMR_XMAX_[0]-_MESH_AMR_XMIN_[0])*(_MESH_AMR_XMAX_[1]-_MESH_AMR_XMIN_[1])*(_MESH_AMR_XMAX_[2]-_MESH_AMR_XMIN_[2]);

    #if _MESH_DIMENSION_ == 1
    vol/=_TOTAL_BLOCK_CELLS_X_*(1<<blockDescriptor.RefinmentLevel);
    #elif _MESH_DIMENSION_ == 2
    vol/=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*(1<<(1+blockDescriptor.RefinmentLevel));
    #elif _MESH_DIMENSION_ == 3
    vol/=_TOTAL_BLOCK_CELLS_X_*_TOTAL_BLOCK_CELLS_Y_*_TOTAL_BLOCK_CELLS_Z_*(1<<(2+blockDescriptor.RefinmentLevel));
    #endif

    return vol; 
  } 

  double CellCharacteristicSize() {
    register double res;
    register int idim;

    res=pow((_MESH_AMR_XMAX_[0]-_MESH_AMR_XMIN_[0])/_TOTAL_BLOCK_CELLS_X_,2);
    if (_MESH_DIMENSION_>1) res+=pow((_MESH_AMR_XMAX_[1]-_MESH_AMR_XMIN_[1])/_TOTAL_BLOCK_CELLS_Y_,2);
    if (_MESH_DIMENSION_>2) res+=pow((_MESH_AMR_XMAX_[2]-_MESH_AMR_XMIN_[2])/_TOTAL_BLOCK_CELLS_Z_,2); 

    res=sqrt(res)/(1<<blockDescriptor.RefinmentLevel);

    return res;
  }


  void cleanDataBuffer() {
    blockDescriptor.RefinmentLevel=0;
    blockDescriptor.ghostBlock=0;

    nextBlock=NULL;

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  } 

  cBasicBlockAMR() {
    cleanDataBuffer();
  }


  bool GetGhostFlag() {return (blockDescriptor.ghostBlock==0) ? _REAL_BLOCK_ : _GHOST_BLOCK_;}
  void SetGhostFlag(bool flag) {blockDescriptor.ghostBlock=(flag==_REAL_BLOCK_) ? 0 : 1;} //if false -> the block is real, if true -> the block is a ghost block 

  int GetRefinmentLevel() {return blockDescriptor.RefinmentLevel;}
  void SetRefinmentLevel(int level) {blockDescriptor.RefinmentLevel=level;}

};



  

//======================================================================================
template <class cCornerNode,class cCenterNode,class cBlockAMR>
class cMeshAMRgeneric : public cAMRexit {
public:
  cTreeNodeAMR<cBlockAMR>  *rootTree;

  //the function that calculates the interpolation coefficients to get an interpolated value for the block's nodes
  //return the number of the interpolation coefficients that was used in the stencil. if the return value <=0 -> the operation is not succesful
  typedef int (*cGetCornerNodesInterpolationCoefficients)(double *x,double *CoefficientsList,cCornerNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients);
  cGetCornerNodesInterpolationCoefficients GetCornerNodesInterpolationCoefficients;

  typedef int (*cGetCenterNodesInterpolationCoefficients)(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients);
  cGetCenterNodesInterpolationCoefficients GetCenterNodesInterpolationCoefficients;

  //accept tree node function
  typedef bool (*cAcceptBlockFunc)(double*,double*);
  cAcceptBlockFunc accepltTreeNodeFunction;
  void SetAcceptBlockFunction(cAcceptBlockFunc t) {accepltTreeNodeFunction=t;}

  //limit of the accuracy of calculation of the mesh parameters
  double EPS;

  //the name and fignature of the mesh: the name contains the data when the mesh has been created; all files associated with the mesh will starts with the mesh name
  char MeshName[STRING_LENGTH];
  unsigned long int MeshSignature;


  //mesh modyfied flag -> is set to true each time the mesh is modified, the flag is set to false when the number of elements and the connectivity list are prepared
  bool meshModifiedFlag;

  //the total number of mesh nodes and blocks 
  long int meshNodesNumber,meshBlocksNumber;

  //the upper (maximum) refinment level used on the mesh 
  int meshMaximumRefinmentLevel;

  //the limits of the mesh
  double xGlobalMin[_MESH_DIMENSION_],xGlobalMax[_MESH_DIMENSION_];

  //the user defined function for the local mesh resolution 
  double (*localResolution)(double*);

  //the stacks for the tree, blocks and nodes
  cAssociatedDataAMRstack<cCenterNode> CenterNodes;
  cAssociatedDataAMRstack<cCornerNode> CornerNodes;
  cAssociatedDataAMRstack<cBlockAMR> blocks;
  cAMRstack<cTreeNodeAMR<cBlockAMR> > treeNodes;

  //the stack containing the descriptors for all internal surfaces installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  cAMRstack<cInternalBoundaryConditionsDescriptor> InternalBoundaryDescriptors;
  #endif

  //the root block;
  cBlockAMR *rootBlock;
  double dxRootBlock[3];

  //the list of blocks of the mesh
  void *blockList;

  //the list of the nodes that "belongs" to this processor
  cTreeNodeAMR<cBlockAMR> *LocalNodes;

  //the list of the descriptors of all internal boundaries installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  list<cInternalBoundaryConditionsDescriptor> InternalBoundaryList;
  #endif

  //generate the mesh signeture: the signature contained the time of the mesh creeation, the user name and the computer name where the lesh is created
  void generateMeshSignature() {
    CRC32 meshSign;
    char str[STRING_LENGTH];
    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);
    char *username;

    //add the time stamp
    meshSign.add(TimeValue);

    //add the hostname
    gethostname(str,STRING_LENGTH);
    meshSign.add(str,STRING_LENGTH);

    //add the mesh name
    meshSign.add(MeshName,STRING_LENGTH);

    //add the user name
    username=getenv("USER");
    if (username!=NULL) meshSign.add(username,20);

    MeshSignature=meshSign.checksum();

    //distribute the signature among all processors
    int mpiInitFlag;

    MPI_Initialized(&mpiInitFlag);

    if (mpiInitFlag==true) {
      MPI_Bcast(&MeshSignature,1,MPI_UNSIGNED_LONG,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  unsigned long getMeshSignature() {
    if (meshModifiedFlag==true) generateMeshSignature();

    return MeshSignature;
  }

  //set the mesh name
  void setMeshName(const char *str) {
    meshModifiedFlag=true;

    sprintf(MeshName,"%s",str);
    generateMeshSignature();
  }

  void getMeshName(char *str) {
    sprintf(str,"%s",MeshName);
  }

  void generateMeshName() {
    char mname[STRING_LENGTH];
    time_t TimeValue=time(0);
    tm *ct=localtime(&TimeValue);

    meshModifiedFlag=true;

    sprintf(mname,"noname.%i-%i-%i.AMR.mesh",ct->tm_mon+1,ct->tm_mday,ct->tm_year+1900);
    setMeshName(mname);
  }

  inline long int getCornerNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+(1+_TOTAL_BLOCK_CELLS_X_)*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*(1+_TOTAL_BLOCK_CELLS_Y_));
    #endif

    return nd;
  }

  inline long int getCenterNodeLocalNumber(int i,int j,int k) {
    long int nd;

    #if _MESH_DIMENSION_ == 1
    nd=i+_GHOST_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_);
    #elif _MESH_DIMENSION_ == 3
    nd=i+_GHOST_CELLS_X_+_TOTAL_BLOCK_CELLS_X_*(j+_GHOST_CELLS_Y_+(k+_GHOST_CELLS_Z_)*_TOTAL_BLOCK_CELLS_Y_);
    #endif

    return nd;
  }

  long int findCornerNodeIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode) { 
    double dx;
    
    //if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) return -1;
    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    i=(x[0]-startNode->xmin[0])/dx;

    if (fabs(startNode->xmin[0]+dx*i-x[0])>EPS) {
      i++; 
      if (fabs(startNode->xmin[0]+dx*i-x[0])>EPS) return -1;
    }

    if (_MESH_DIMENSION_>=2) {
      //if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) return -1;
      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      j=(x[1]-startNode->xmin[1])/dx;

      if (fabs(startNode->xmin[1]+dx*j-x[1])>EPS) {
        j++;
        if (fabs(startNode->xmin[1]+dx*j-x[1])>EPS) return -1;
      }
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      //if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) return -1;
      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      k=(x[2]-startNode->xmin[2])/dx;

      if (fabs(startNode->xmin[2]+dx*k-x[2])>EPS) {
        k++;
        if (fabs(startNode->xmin[2]+dx*k-x[2])>EPS) return -1;
      }
    }
    else k=0;

    return getCornerNodeLocalNumber(i,j,k);
  }  


  long int findCenterNodeIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode) {
    double dx,dx2;
  
    //if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) return -1;
    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    dx2=dx/2.0;
    i=(x[0]-dx2-startNode->xmin[0])/dx;

    if (fabs(startNode->xmin[0]+dx2+dx*i-x[0])>EPS) {
      i++;
      if (fabs(startNode->xmin[0]+dx2+dx*i-x[0])>EPS) return -1;
    }

    if (_MESH_DIMENSION_>=2) {
      //if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) return -1;
      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      dx2=dx/2.0;
      j=(x[1]-dx2-startNode->xmin[1])/dx;

      if (fabs(startNode->xmin[1]+dx2+dx*j-x[1])>EPS) {
        j++;
        if (fabs(startNode->xmin[1]+dx2+dx*j-x[1])>EPS) return -1;
      }
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      //if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) return -1;
      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      dx2=dx/2.0;
      k=(x[2]-dx2-startNode->xmin[2])/dx;

      if (fabs(startNode->xmin[2]+dx2+dx*k-x[2])>EPS) {
        k++;
        if (fabs(startNode->xmin[2]+dx2+dx*k-x[2])>EPS) return -1;
      }
    }
    else k=0;

    return getCenterNodeLocalNumber(i,j,k);
  }

  //find the index of the cell where the point 'x' is located
  long int fingCellIndex(double *x,int &i,int &j,int &k,cTreeNodeAMR<cBlockAMR>* startNode) {
    double dx;

    if ((x[0]<startNode->xmin[0])||(startNode->xmax[0]<x[0])) exit(__LINE__,__FILE__,"x is outside of the block");
    dx=dxRootBlock[0]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
    i=(int)((x[0]-startNode->xmin[0])/dx);

    if (_MESH_DIMENSION_>=2) {
      if ((x[1]<startNode->xmin[1])||(startNode->xmax[1]<x[1])) exit(__LINE__,__FILE__,"x is outside of the block");
      dx=dxRootBlock[1]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_);
      j=(int)((x[1]-startNode->xmin[1])/dx);
    }
    else j=0;

    if (_MESH_DIMENSION_==3) {
      if ((x[2]<startNode->xmin[2])||(startNode->xmax[2]<x[2])) exit(__LINE__,__FILE__,"x is outside of the block");
      dx=dxRootBlock[2]/(1<<startNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_);
      k=(int)((x[2]-startNode->xmin[2])/dx);
    }
    else k=0;

    return getCenterNodeLocalNumber(i,j,k);
  }
 

  //the constructor
  void init (double *xMin,double *xMax,double (*localResolutionFunction)(double*)) {
    int i,j,k,nd,idim;

    //check the number of bits reserved to store the number of connection's of a corner node
    if ((1<<_MAX_CORNER_NODE_CONNECTION_BITS_)<_MAX_CORNER_NODE_CONNECTION_+1) exit(__LINE__,__FILE__,"The value of _MAX_CORNER_NODE_CONNECTION_BITS_ in not sufficient to store _MAX_CORNER_NODE_CONNECTION_");

    //check if the number of bits reserved to store the blocks refinment level is sufficient to store all levels' numbers
    if ((1<<_MAX_REFINMENT_LEVEL_BITS_)<=_MAX_REFINMENT_LEVEL_) exit(__LINE__,__FILE__,"The value of _MAX_REFINMENT_LEVEL_BITS_ is not sufficient to store all refinment levels");
    if ((_MESH_DIMENSION_<1)||(_MESH_DIMENSION_>3)) exit(__LINE__,__FILE__,"The mesh dimension is wrong");

    xGlobalMin[0]=xMin[0],xGlobalMax[0]=xMax[0],dxRootBlock[0]=(xMax[0]-xMin[0]);
    if (2*_GHOST_CELLS_X_>=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
    EPS=0.0001*dxRootBlock[0]/double(_BLOCK_CELLS_X_)/(1<<_MAX_REFINMENT_LEVEL_); 

    if (_MESH_DIMENSION_>1) {
      xGlobalMin[1]=xMin[1],xGlobalMax[1]=xMax[1],dxRootBlock[1]=(xMax[1]-xMin[1]); 
      if (2*_GHOST_CELLS_Y_>=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
      if (EPS>0.0001*dxRootBlock[1]/double(_BLOCK_CELLS_Y_)/(1<<_MAX_REFINMENT_LEVEL_)) EPS=0.0001*dxRootBlock[1]/double(_BLOCK_CELLS_Y_)/(1<<_MAX_REFINMENT_LEVEL_); 
    }

    if (_MESH_DIMENSION_>2) {
      xGlobalMin[2]=xMin[2],xGlobalMax[2]=xMax[2],dxRootBlock[2]=(xMax[2]-xMin[2]);
      if (2*_GHOST_CELLS_Z_>=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__,"The mesh dimension is wrong");
      if (EPS>0.0001*dxRootBlock[2]/double(_BLOCK_CELLS_Z_)/(1<<_MAX_REFINMENT_LEVEL_)) EPS=0.0001*dxRootBlock[2]/double(_BLOCK_CELLS_Z_)/(1<<_MAX_REFINMENT_LEVEL_);
    }

    for (idim=0;idim<_MESH_DIMENSION_;idim++) _MESH_AMR_XMAX_[idim]=xMax[idim],_MESH_AMR_XMIN_[idim]=xMin[idim]; 
 
    //set the default value for the 'interpolation functions'
    GetCenterNodesInterpolationCoefficients=NULL;
    GetCornerNodesInterpolationCoefficients=NULL;

    localResolution=localResolutionFunction;

    //set up the tree and the root block
    rootBlock=blocks.newElement();
    rootTree=treeNodes.newElement();


    meshNodesNumber=0,meshBlocksNumber=0,meshModifiedFlag=true,meshMaximumRefinmentLevel=0;
    rootTree->upNode=NULL;
    rootTree->block=rootBlock;
    for (i=0;i<(1<<_MESH_DIMENSION_);i++) rootTree->downNode[i]=NULL; 

    rootTree->RefinmentLevel=0;
    for (i=0;i<_MESH_DIMENSION_;i++) rootTree->xmin[i]=xMin[i],rootTree->xmax[i]=xMax[i];

    accepltTreeNodeFunction=NULL;
    MeshName[0]='\0',MeshSignature=0;
    blockList=NULL;

    //set the corner nodes of the block
    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"Adjust as in 3D");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"Adjust as in 3D");
    #elif _MESH_DIMENSION_ == 3
    int iMax=_BLOCK_CELLS_X_,jMax=_BLOCK_CELLS_Y_,kMax=_BLOCK_CELLS_Z_;
    #else
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif

    //set up the corner nodes
    cCornerNode *newCornerNode;
    double xNode[3];

    for (k=0;k<=kMax;k++) for (j=0;j<=jMax;j++) for (i=0;i<=iMax;i++) {
      nd=getCornerNodeLocalNumber(i,j,k);

      newCornerNode=CornerNodes.newElement();
      rootBlock->SetCornerNode(newCornerNode,nd);
      newCornerNode->incrementConnectionCounter();

      xNode[0]=xGlobalMin[0]+i*dxRootBlock[0]/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>=2) xNode[1]=xGlobalMin[1]+j*dxRootBlock[1]/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_==3) xNode[2]=xGlobalMin[2]+k*dxRootBlock[2]/_BLOCK_CELLS_Z_;

      newCornerNode->SetX(xNode);
    }


    //set the center node 
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    cCenterNode *newCenterNode;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      nd=getCenterNodeLocalNumber(i,j,k);

      newCenterNode=CenterNodes.newElement();
      rootBlock->SetCenterNode(newCenterNode,nd);
      newCenterNode->incrementConnectionCounter();

      xNode[0]=xGlobalMin[0]+(i+0.5)*dxRootBlock[0]/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>=2) xNode[1]=xGlobalMin[1]+(j+0.5)*dxRootBlock[1]/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_==3) xNode[2]=xGlobalMin[2]+(k+0.5)*dxRootBlock[2]/_BLOCK_CELLS_Z_;

      newCenterNode->SetX(xNode);

    }
    #endif
 
    //install the internal surfaces into the rooBlock
    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    list<cInternalBoundaryConditionsDescriptor>::iterator ptr;
    cInternalBoundaryConditionsDescriptor *newDescriptor;

    for (ptr=InternalBoundaryList.begin();ptr!=InternalBoundaryList.end();ptr++) {
      newDescriptor=InternalBoundaryDescriptors.newElement();

      //copy the content of the descriptor
      newDescriptor->BondaryType=ptr->BondaryType;
      newDescriptor->BoundaryElement=ptr->BoundaryElement;

      //add the descriptor to the rootBlock
      newDescriptor->nextInternalBCelement=rootTree->InternalBoundaryDescriptorList;
      rootTree->InternalBoundaryDescriptorList=newDescriptor;
    }

    #endif
  }  

  cMeshAMRgeneric (double *xMin,double *xMax,double (*localResolutionFunction)(double*)) {
    init(xMin,xMax,localResolutionFunction);
  }


  cMeshAMRgeneric (){
    for (int idim=0;idim<_MESH_DIMENSION_;idim++) _MESH_AMR_XMAX_[idim]=0.0,_MESH_AMR_XMIN_[idim]=0.0;

     //set the default value for the 'interpolation functions'
     GetCenterNodesInterpolationCoefficients=NULL;
     GetCornerNodesInterpolationCoefficients=NULL;
     localResolution=NULL;

     //set up the tree and the root block
     rootBlock=NULL;
     rootTree=NULL;
     meshNodesNumber=0,meshBlocksNumber=0,meshModifiedFlag=true,meshMaximumRefinmentLevel=0;

     accepltTreeNodeFunction=NULL;
     MeshName[0]='\0',MeshSignature=0;
     blockList=NULL;
     LocalNodes=NULL;
  }

  //register the 'internal boundary' (the surface determining cut cells)
  void RegisterInternalBoundary(cInternalBoundaryConditionsDescriptor Descriptor) {
    #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
    if (rootTree!=NULL) exit(__LINE__,__FILE__,"Error: all internal surface must be registered before initialization of the mesh");

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"Error: internal surfaces are not allowed in _MESH_DIMENSION_ == 1");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 3
    if (Descriptor.BondaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_) exit(__LINE__,__FILE__,"Error: Attempted boundary type is not allowed for _MESH_DIMENSION_ == 3");
    #else
    exit(__LINE__,__FILE__,"Error: unknown value of _MESH_DIMENSION_");
    #endif


    InternalBoundaryList.push_back(Descriptor);
    #else
    exit(__LINE__,__FILE__,"Error: internal boundary is allowed only when _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_");
    #endif
  }


  cTreeNodeAMR<cBlockAMR>* findTreeNode(double *x,cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    cTreeNodeAMR<cBlockAMR> *res=NULL,*t=NULL;
    double blockBasedCoordinates[3]={0.0,0.0,0.0};
    int i,j=0,k=0;

    if (startNode==NULL) startNode=rootTree;

    blockBasedCoordinates[0]=(x[0]-startNode->xmin[0])/dxRootBlock[0]*(1<<startNode->RefinmentLevel);
    if (_MESH_DIMENSION_>1) blockBasedCoordinates[1]=(x[1]-startNode->xmin[1])/dxRootBlock[1]*(1<<startNode->RefinmentLevel);
    if (_MESH_DIMENSION_>2) blockBasedCoordinates[2]=(x[2]-startNode->xmin[2])/dxRootBlock[2]*(1<<startNode->RefinmentLevel);

    if ((0.0<=blockBasedCoordinates[0])&&(blockBasedCoordinates[0]<=1.0)&&(0.0<=blockBasedCoordinates[1])&&(blockBasedCoordinates[1]<=1.0)&&(0.0<=blockBasedCoordinates[2])&&(blockBasedCoordinates[2]<=1.0)) {
      i=(blockBasedCoordinates[0]<0.5) ? 0 : 1;
      if (_MESH_DIMENSION_>1) j=(blockBasedCoordinates[1]<0.5) ? 0 : 1;
      if (_MESH_DIMENSION_>2) k=(blockBasedCoordinates[2]<0.5) ? 0 : 1;

      t=startNode->downNode[i+2*(j+2*k)]; 

      res=(t!=NULL) ? findTreeNode(x,t) : startNode; 
    }
    else res=(startNode->upNode!=NULL) ? findTreeNode(x,startNode->upNode) : NULL;

    return res;
  }



  void reconnectDownTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
    cTreeNodeAMR<cBlockAMR> *neibNode,*neibDownNode; 

    #if _MESH_DIMENSION_ == 1
    const static int neibNodeNumber[2][2]={ {0,-1}, {-1,1}};
    const static int neibDownNodeNumber[2][2]={ {1,-1}, {-1,0}};
    const static int reverseFaceMap[2]={1,0};

    #else

    #if _MESH_DIMENSION_ == 2
 
    //the face number of the parent block that the downBlock is connecter through with the surrounding blocks
    const static int neibNodeNumber[4][4]={ {0,1,-1,-1}, {-1,1,2,-1}, {0,-1,-1,3}, {-1,-1,2,3}};  

    //the number of the downNodes in the neibNode that the startNode->downNode[?] is connected with
    const static int neibDownNodeNumber[4][4]={ {1,2,-1,-1}, {-1,3,0,-1}, {3,-1,-1,0}, {-1,-1,2,1}}; 

    //the corresponding face number that neibNode is connected to startNode
    const static int reverseFaceMap[4]={2,3,0,1};

    #else

    //the face number of the parent block that the downBlock is connecter through with the surrounding blocks
    const static int neibNodeNumber[8][6]={ {0,-1,2,-1,4,-1}, {-1,1,2,-1,4,-1}, {0,-1,-1,3,4,-1}, {-1,1,-1,3,4,-1},
                                            {0,-1,2,-1,-1,5}, {-1,1,2,-1,-1,5}, {0,-1,-1,3,-1,5}, {-1,1,-1,3,-1,5},   };

    const static int neibDownNodeNumber[8][6]={ {1,-1,2,-1,4,-1,}, {-1,0,3,-1,5,-1}, {3,-1,-1,0,6,-1}, {-1,2,-1,1,7,-1},
                                                {5,-1,6,-1,-1,0}, {-1,4,7,-1,-1,1,}, {7,-1,-1,4,-1,2}, {-1,6,-1,5,-1,3}  }; 

    const static int reverseFaceMap[6]={1,0,3,2,5,4};



    #endif

    int nDownNode,nface;



    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if (neibNodeNumber[nDownNode][nface]!=-1) {
      neibNode=startNode->neibNode[neibNodeNumber[nDownNode][nface]];  

      if (neibNode!=NULL) {
        neibDownNode=neibNode->downNode[neibDownNodeNumber[nDownNode][nface]];

        startNode->downNode[nDownNode]->neibNode[nface]=neibDownNode;
        if (neibDownNode!=NULL) neibDownNode->neibNode[ reverseFaceMap[nface] ]=startNode->downNode[nDownNode]; 
      }
    }


    #endif
  }




  cTreeNodeAMR<cBlockAMR> *getNeibNode_2D(int i,int j,cTreeNodeAMR<cBlockAMR>* startNode) {
    double x[3]={0.0,0.0,0.0};

    x[0]=startNode->xmin[0]+0.5*dxRootBlock[0]/(1<<startNode->RefinmentLevel)*(1.0+1.01*i);  
    x[1]=startNode->xmin[1]+0.5*dxRootBlock[1]/(1<<startNode->RefinmentLevel)*(1.0+1.01*j);

    return findTreeNode(x,startNode); 
  }

  cTreeNodeAMR<cBlockAMR> *getNeibNode(int i,int j,int k,cTreeNodeAMR<cBlockAMR>* startNode) {
    double x[3]={0.0,0.0,0.0};
    cTreeNodeAMR<cBlockAMR> *res;

    x[0]=startNode->xmin[0]+0.5*dxRootBlock[0]/(1<<startNode->RefinmentLevel)*(1.0+1.01*i);
    if (_MESH_DIMENSION_>1) x[1]=startNode->xmin[1]+0.5*dxRootBlock[1]/(1<<startNode->RefinmentLevel)*(1.0+1.01*j);
    if (_MESH_DIMENSION_>2) x[2]=startNode->xmin[2]+0.5*dxRootBlock[2]/(1<<startNode->RefinmentLevel)*(1.0+1.01*k);

    res=findTreeNode(x,startNode);
    
    return res;
  }
 

 void collectNeibNodeMap_1D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
   int i,j,nface,iDownBlock,iNeib,iOffset;
   cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;


exit(__LINE__,__FILE__,"update as in 3D collect"); 


   //if the neibNode->RefinmentLevel == startNodes->RefinementLevel -> the 'neibNode' does not contribute to the 'newCornerNodeMap'

   //check if the face neibours can contribute to the lists
   for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=startNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
     //transfer the nodes

     static const int _X_FACE_NEIB_OFFSET_[2]={-2*_BLOCK_CELLS_X_,2*_BLOCK_CELLS_X_};

     for (iDownBlock=0;iDownBlock<2;iDownBlock++) { 
       downNode=neibNode->downNode[iDownBlock];
       if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"the branch is too deap");
       iOffset=_X_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_X_*iDownBlock;

       for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
         newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][0][0]=downNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,0,0)];
       }

     }
   }
 }


  void collectNeibNodeMap_2D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int i,j,nface,iDownBlock,jDownBlock,iNeib,jNeib,iOffset,jOffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;

exit(__LINE__,__FILE__,"update as in 3D collect");

    //if the neibNode->RefinmentLevel == startNodes->RefinementLevel -> the 'neibNode' does not contribute to the 'newCornerNodeMap' 

    //check if the face neibours can contribute to the lists
    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=startNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    //transfer the nodes

    static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
    static const int jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;

    
    static const int _X_FACE_NEIB_OFFSET_[4]={-2*_BLOCK_CELLS_X_,0,2*_BLOCK_CELLS_X_,0};
    static const int _Y_FACE_NEIB_OFFSET_[4]={0,-2*_BLOCK_CELLS_Y_,0,2*_BLOCK_CELLS_Y_}; 
 
    for (jDownBlock=0;jDownBlock<2;jDownBlock++) for (iDownBlock=0;iDownBlock<2;iDownBlock++) {
      downNode=neibNode->downNode[iDownBlock+2*jDownBlock]; 

      iOffset=_X_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_X_*iDownBlock;
      jOffset=_Y_FACE_NEIB_OFFSET_[nface]+_BLOCK_CELLS_Y_*jDownBlock; 

      //chech if 'downNode' can overlap 'startNode'
      int nIntersection=0;

      if (iOffset>0) {
        if (iOffset<=2*_BLOCK_CELLS_X_) nIntersection++;
      }
      else if (iOffset<0) {
        if (iOffset==-_BLOCK_CELLS_X_) nIntersection++;
      }
      else nIntersection++;

      if (jOffset>0) {
        if (jOffset<=2*_BLOCK_CELLS_Y_) nIntersection++;
      }
      else if (jOffset<0) {
        if (jOffset==-_BLOCK_CELLS_Y_) nIntersection++;
      }
      else nIntersection++;

      if (nIntersection!=2) continue;

      if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) exit(__LINE__,__FILE__,"the branch is too deap");

      for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {   
        for (jNeib=-_GHOST_CELLS_Y_;jNeib<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jNeib++) if ((-_GHOST_CELLS_Y_<=jNeib+jOffset)&&(jNeib+jOffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {  
          newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][jNeib+jOffset+_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,jNeib,0)]; 
        }
      }
    }
  }

  //check the corner blocks
  int iNeibBlock,jNeibBlock;

  for (iNeibBlock=-1;iNeibBlock<2;iNeibBlock+=2) for (jNeibBlock=-1;jNeibBlock<2;jNeibBlock+=2) if ((neibNode=getNeibNode_2D(iNeibBlock,jNeibBlock,startNode))!=NULL) if (neibNode->RefinmentLevel>startNode->RefinmentLevel) { 
    //transfer the nodes
    if ((neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_)||(neibNode->RefinmentLevel!=startNode->RefinmentLevel+1)) exit(__LINE__,__FILE__,"the branch is too deap");
     
    iOffset=_BLOCK_CELLS_X_*((iNeibBlock==-1) ? -1 : 2);
    jOffset=_BLOCK_CELLS_Y_*((jNeibBlock==-1) ? -1 : 2); 

    for (iNeib=-_GHOST_CELLS_X_;iNeib<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;iNeib++) if ((-_GHOST_CELLS_X_<=iNeib+iOffset)&&(iNeib+iOffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
      for (jNeib=-_GHOST_CELLS_Y_;jNeib<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jNeib++) if ((-_GHOST_CELLS_Y_<=jNeib+jOffset)&&(jNeib+jOffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
        newCornerNodeMap[iNeib+iOffset+_GHOST_CELLS_X_][jNeib+jOffset+_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[getCornerNodeLocalNumber(iNeib,jNeib,0)];
      }
    }

  }

}
 





  void collectNeibCornerNodeMap_3D(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int nIntersection,ii,jj,kk,iNeibNode,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode,*upNode;

    //init the node's map
    for (ii=0;ii<1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_);ii++) for (jj=0;jj<1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_);jj++) for (kk=0;kk<1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_);kk++) newCornerNodeMap[ii][jj][kk]=NULL;

    //transfer nodes to the 'newCornerNodeMap' from the 'cCornerNode'
    for (ii=-(_GHOST_CELLS_X_/2);ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_/2;ii++) for (jj=-(_GHOST_CELLS_Y_/2);jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_/2;jj++) for (kk=-(_GHOST_CELLS_Z_/2);kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_/2;kk++) {
      newCornerNodeMap[2*ii+_GHOST_CELLS_X_][2*jj+_GHOST_CELLS_Y_][2*kk+_GHOST_CELLS_Z_]=startNode->block->GetCornerNode(getCornerNodeLocalNumber(ii,jj,kk));
    }

    //transfer nodes from the neibouring blocks
    for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode; 

        for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
          downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

          ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
          joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
          koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

          nIntersection=0;

          //chech if 'downNode' can overlap 'startNode'
          if ((ioffset==0)||(ioffset==2*_BLOCK_CELLS_X_)||(ioffset==_BLOCK_CELLS_X_)||(ioffset==-_BLOCK_CELLS_X_)) nIntersection++; 
          if ((joffset==0)||(joffset==2*_BLOCK_CELLS_Y_)||(joffset==_BLOCK_CELLS_Y_)||(joffset==-_BLOCK_CELLS_Y_)) nIntersection++;
          if ((koffset==0)||(koffset==2*_BLOCK_CELLS_Z_)||(koffset==_BLOCK_CELLS_Z_)||(koffset==-_BLOCK_CELLS_Z_)) nIntersection++;

          if (nIntersection!=3) continue;

          for (ii=-_GHOST_CELLS_X_;ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
            for (jj=-_GHOST_CELLS_Y_;jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
              for (kk=-_GHOST_CELLS_Z_;kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
                newCornerNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(getCornerNodeLocalNumber(ii,jj,kk));
              }
            }
          }

        }
      }
    }
 
  }




	
  void collectNeibCenterNodeMap_3D(cTreeNodeAMR<cBlockAMR> *startNode,cCenterNode *newCenterNodeMap[2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)]) {
    int nIntersection,ii,jj,kk,iNeibNode,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
    cTreeNodeAMR<cBlockAMR> *neibNode,*downNode,*upNode;

    //init the node's map
    for (ii=0;ii<2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_);ii++) for (jj=0;jj<2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_);jj++) for (kk=0;kk<2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_);kk++) newCenterNodeMap[ii][jj][kk]=NULL;

    //transfer nodes from the neibouring blocks
    for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode;

       //transfer nodes from the neibouring blocks
       for (iNeibNode=-1;iNeibNode<=1;iNeibNode++) for (jNeibNode=-1;jNeibNode<=1;jNeibNode++) for (kNeibNode=-1;kNeibNode<=1;kNeibNode++) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
        if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
          upNode=neibNode->upNode;

          for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
            downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

            ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
            joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
            koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

            nIntersection=0;

            //chech if 'downNode' can overlap 'startNode'
            if ((ioffset==0)||(ioffset==2*_BLOCK_CELLS_X_)||(ioffset==_BLOCK_CELLS_X_)||(ioffset==-_BLOCK_CELLS_X_)) nIntersection++;
            if ((joffset==0)||(joffset==2*_BLOCK_CELLS_Y_)||(joffset==_BLOCK_CELLS_Y_)||(joffset==-_BLOCK_CELLS_Y_)) nIntersection++;
            if ((koffset==0)||(koffset==2*_BLOCK_CELLS_Z_)||(koffset==_BLOCK_CELLS_Z_)||(koffset==-_BLOCK_CELLS_Z_)) nIntersection++;

            if (nIntersection!=3) continue;

            for (ii=-_GHOST_CELLS_X_;ii<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_-1)) {
              for (jj=-_GHOST_CELLS_Y_;jj<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_-1)) {
                for (kk=-_GHOST_CELLS_Z_;kk<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_-1)) {
                  newCenterNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=downNode->block->GetCenterNode(getCenterNodeLocalNumber(ii,jj,kk));
                }
              }
            }

          }
        }
      }

    }
  }
}



	

void checkMeshConsistency(cTreeNodeAMR<cBlockAMR> *startNode) {
  int i,j,k,iDownBlock,jDownBlock,kDownBlock,iStartNode,iNeibNode;
  cTreeNodeAMR<cBlockAMR> *neibNode,*downNode;
  bool flag;

  #if _MESH_DIMENSION_ == 1
  static const int iDownNodeMax=1,jDownNodeMax=0,kDownNodeMax=0;
  static const int iMin=-1,iMax=2,jMin=0,jMax=1,kMin=0,kMax=1;
  static const int nCornerNodes=1+_TOTAL_BLOCK_CELLS_X_;
  #elif _MESH_DIMENSION_ == 2
  static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=0;
  static const int iMin=-1,iMax=2,jMin=-1,jMax=2,kMin=0,kMax=1;
  static const int nCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);  
  #else 
  static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=1;
  static const int iMin=-1,iMax=2,jMin=-1,jMax=2,kMin=-1,kMax=2;
  static const int nCornerNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
  #endif



  if (startNode==rootTree) {
    //reset the 'processed flag' for the nodes
    resetNodeProcessedFlag();
  }

  //check the 'xmax' and 'xmin' value of the node
  if (startNode!=rootTree) {
    double *xminUpNode,*xmaxUpNode,*xmin,*xmax;
    int kDownBlock,jDownBlock,iDownBlock;
    bool flag;

    flag=false;
    xminUpNode=startNode->upNode->xmin,xmaxUpNode=startNode->upNode->xmax;
    xmin=startNode->xmin,xmax=startNode->xmax;

    for (kDownBlock=0;kDownBlock<=kDownNodeMax;kDownBlock++) for (jDownBlock=0;jDownBlock<=jDownNodeMax;jDownBlock++) for (iDownBlock=0;iDownBlock<=iDownNodeMax;iDownBlock++) {
      if (startNode->upNode->downNode[iDownBlock+2*(jDownBlock+2*kDownBlock)]==startNode) {
        flag=true; //there is a connection between the upNode and the startNode;

        if (fabs(xminUpNode[0]+(xmaxUpNode[0]-xminUpNode[0])*iDownBlock/2.0-xmin[0])>EPS) exit(__LINE__,__FILE__,"Block's xmin value if not correct");
        if (_MESH_DIMENSION_>=2) if (fabs(xminUpNode[1]+(xmaxUpNode[1]-xminUpNode[1])*jDownBlock/2.0-xmin[1])>EPS) exit(__LINE__,__FILE__,"Block's xmin value if not correct");
        if (_MESH_DIMENSION_==3) if (fabs(xminUpNode[2]+(xmaxUpNode[2]-xminUpNode[2])*kDownBlock/2.0-xmin[2])>EPS) exit(__LINE__,__FILE__,"Block's xmin value if not correct");

        if (fabs(xminUpNode[0]+(xmaxUpNode[0]-xminUpNode[0])*(iDownBlock+1)/2.0-xmax[0])>EPS) exit(__LINE__,__FILE__,"Block's xmax value if not correct");
        if (_MESH_DIMENSION_>=2) if (fabs(xminUpNode[1]+(xmaxUpNode[1]-xminUpNode[1])*(jDownBlock+1)/2.0-xmax[1])>EPS) exit(__LINE__,__FILE__,"Block's xmax value if not correct");
        if (_MESH_DIMENSION_==3) if (fabs(xminUpNode[2]+(xmaxUpNode[2]-xminUpNode[2])*(kDownBlock+1)/2.0-xmax[2])>EPS) exit(__LINE__,__FILE__,"Block's xmax value if not correct");
      }
    }

    if (flag==false) exit(__LINE__,__FILE__,"There is a disconnection between the upNode and the startNode");
  }


  //check the tree node 
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {

     //check the neibours 
     for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++)  if ((neibNode=getNeibNode(i,j,k,startNode))!=NULL) {
       if (fabs(neibNode->RefinmentLevel-startNode->RefinmentLevel)>1) {
         exit(__LINE__,__FILE__,"The resolution change exseeds the limit"); 
       }

       //check if the startNode and neibNode share a cornerNode 
       for (flag=false,iStartNode=0;(iStartNode<nCornerNodes)&&(flag==false);iStartNode++) for (iNeibNode=0;iNeibNode<nCornerNodes;iNeibNode++) {
         if (startNode->block->GetCornerNode(iStartNode)==neibNode->block->GetCornerNode(iNeibNode)) {
            flag=true;
            break;
         }
       }  

       if (flag==false) {
         long int nd,ii,jj,kk,idim;

         cout << "ERROR: Two neiubouring blocks doesn't share nodes" << endl;  
         cout << "Block 1: ID=" << startNode->Temp_ID << endl << "Corner Nodes:" << endl;

         for (ii=0;ii<iMax;ii++) for (jj=0;jj<jMax;jj++) for (kk=0;kk<kMax;kk++) {
           nd=getCornerNodeLocalNumber(ii*_BLOCK_CELLS_X_,jj*_BLOCK_CELLS_Y_,kk*_BLOCK_CELLS_Z_);
  
           cout << "(i,j,k)=" << ii <<"," << jj << "," << kk << " (nd=" << nd << "), Temp_ID=" << startNode->block->GetCornerNode(nd)->Temp_ID << " x=";
           for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << startNode->block->GetCornerNode(nd)->GetX()[idim] << " ";
           cout << endl; 
         }

         cout << "Block 2: ID=" << neibNode->Temp_ID << endl << "Corner Nodes:";
         cout << "neib coordinates=" << i << ", " << j << ", " << k << endl; 

         for (ii=0;ii<iMax;ii++) for (jj=0;jj<jMax;jj++) for (kk=0;kk<kMax;kk++) {
           nd=getCornerNodeLocalNumber(ii*_BLOCK_CELLS_X_,jj*_BLOCK_CELLS_Y_,kk*_BLOCK_CELLS_Z_);
  
           cout << "(i,j,k)=" << ii <<"," << jj << "," << kk << " (nd=" << nd << "), Temp_ID=" << neibNode->block->GetCornerNode(nd)->Temp_ID << " x=";
           for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << neibNode->block->GetCornerNode(nd)->GetX()[idim] << " ";
           cout << endl;
         }

         exit(__LINE__,__FILE__,"The blocks do not share a node"); 
       }

     }

     //check the value of the corner nodes in the block
     int ii,jj,kk;
     long int nd;
     double *xmin,*xmax,dx[3]={0.0,0.0,0.0},*xnode;

     #if _MESH_DIMENSION_ == 1
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=0,jjCornerMax=0,kkCornerMin=0,kkCornerMax=0; 
     #elif _MESH_DIMENSION_ == 2
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=-_GHOST_CELLS_Y_,jjCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkCornerMin=0,kkCornerMax=0;
     #else
     static const int iiCornerMin=-_GHOST_CELLS_X_,iiCornerMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjCornerMin=-_GHOST_CELLS_Y_,jjCornerMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkCornerMin=-_GHOST_CELLS_Z_,kkCornerMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
     #endif

     xmin=startNode->xmin,xmax=startNode->xmax;
     dx[0]=(xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
     if (_MESH_DIMENSION_>=2) dx[1]=(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_; 
     if (_MESH_DIMENSION_==3) dx[2]=(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

     for (kk=kkCornerMin;kk<=kkCornerMax;kk++) for (jj=jjCornerMin;jj<=jjCornerMax;jj++)  for (ii=iiCornerMin;ii<=iiCornerMax;ii++) {
       nd=getCornerNodeLocalNumber(ii,jj,kk);

       if (startNode->block->GetCornerNode(nd)==NULL) {
         if ((ii>=0)&&(ii<=_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<=_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<=_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
         else {
           int iNeib=0,jNeib=0,kNeib=0;
 
           if (ii<0) iNeib=-1;
           if (ii>_BLOCK_CELLS_X_) iNeib=1;
 
           if (jj<0) jNeib=-1;
           if (jj>_BLOCK_CELLS_Y_) jNeib=1;
 
           if (kk<0) kNeib=-1;
           if (kk>_BLOCK_CELLS_Z_) kNeib=1; 
 
           if (getNeibNode(iNeib,jNeib,kNeib,startNode)!=NULL) {
             exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
           }
           else continue; 
         }
       } 

       xnode=startNode->block->GetCornerNode(nd)->GetX();

       if (fabs(xmin[0]+dx[0]*ii-xnode[0])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct"); 
       if (_MESH_DIMENSION_>=2) if (fabs(xmin[1]+dx[1]*jj-xnode[1])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct"); 
       if (_MESH_DIMENSION_==3) if (fabs(xmin[2]+dx[2]*kk-xnode[2])>EPS) exit(__LINE__,__FILE__,"The position of the block's corner nodes is not correct");
     }  

     //check if the nodes that the block shares with its neibours the have the same coordinates points to the same location in physical memory
     int ioffset,joffset,koffset,iNeibNode,jNeibNode,kNeibNode,iDownNode,jDownNode,kDownNode;
     cTreeNodeAMR<cBlockAMR> *upNode,*downNode;
     double *xStartNode,*xNeibNode; 
     long int ndStartNode,ndNeibNode;
     int iiDownCornerNode=0,jjDownCornerNode=0,kkDownCornerNode=0;

     #if _MESH_DIMENSION_ == 1
     exit(__LINE__,__FILE__,"not implemented"); 
     #elif _MESH_DIMENSION_ == 2 
     exit(__LINE__,__FILE__,"not implemented");
     #else 
     static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;
     static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
     static const int kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
     static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
     static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=-1,kNeibNodeMax=1;
     #endif 

     for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        for (ii=iiMin;ii<=iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndNeibNode=getCornerNodeLocalNumber(ii,jj,kk);
              ndStartNode=getCornerNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);

              xStartNode=startNode->block->GetCornerNode(ndStartNode)->GetX();
              xNeibNode=neibNode->block->GetCornerNode(ndNeibNode)->GetX();

              if ((xNeibNode==NULL)&&(xStartNode==NULL)) continue;
              else if ((xNeibNode==NULL)||(xStartNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
              else {
                if (xStartNode!=xNeibNode) exit(__LINE__,__FILE__,"Error: two different corner nodes have found that has the same coordinate 'x'");
              }
            }
          }
        }
      }
      else if (neibNode->RefinmentLevel>startNode->RefinmentLevel) {
        upNode=neibNode->upNode;

        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        for (ii=iiMin;ii<=iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndStartNode=getCornerNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);
              xStartNode=startNode->block->GetCornerNode(ndStartNode)->GetX();

              for (iDownNode=0;iDownNode<iDownNodeMax;iDownNode++) for (jDownNode=0;jDownNode<jDownNodeMax;jDownNode++)  for (kDownNode=0;kDownNode<kDownNodeMax;kDownNode++) {
                downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

                //find the index of the node neibNode->block->cornerNodes(ii,jj,kk) in the downNode 
                iiDownCornerNode=2*ii-_BLOCK_CELLS_X_*iDownNode; 
                if (_MESH_DIMENSION_>=2) jjDownCornerNode=2*jj-_BLOCK_CELLS_Y_*jDownNode; 
                if (_MESH_DIMENSION_==3) kkDownCornerNode=2*kk-_BLOCK_CELLS_Z_*kDownNode;

                if ((-_GHOST_CELLS_X_<=iiDownCornerNode)&&(iiDownCornerNode<=_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
                  if ((-_GHOST_CELLS_Y_<=jjDownCornerNode)&&(jjDownCornerNode<=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
                    if ((-_GHOST_CELLS_Z_<=kkDownCornerNode)&&(kkDownCornerNode<=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
                      ndNeibNode=getCornerNodeLocalNumber(iiDownCornerNode,jjDownCornerNode,kkDownCornerNode);
                      xNeibNode=downNode->block->GetCornerNode(ndNeibNode)->GetX();

                      if ((xNeibNode==NULL)&&(xStartNode==NULL)) continue;
                      else if ((xNeibNode==NULL)||(xStartNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
                      else {
                        if (xStartNode!=xNeibNode) exit(__LINE__,__FILE__,"Error: two different corner nodes have found that has the same coordinate 'x'");
                      }
                    }
                  }
                }
           
              }

            }
          }
        }

      }
    }

    #if _AMR_CENTER_NODE_  == _ON_AMR_MESH_
    cCenterNode *neibCenterNode,*startCenterNode;


    //check if all centeral nodes (including in 'ghost' cells are defined)
    for (k=((_MESH_DIMENSION_==3) ? -_GHOST_CELLS_Z_ : 0);k<((_MESH_DIMENSION_==3) ? _BLOCK_CELLS_Z_+_GHOST_CELLS_Z_ : 1);k++) for (j=((_MESH_DIMENSION_>=2) ? -_GHOST_CELLS_Y_ : 0);j<((_MESH_DIMENSION_>=2) ? _BLOCK_CELLS_Y_+_GHOST_CELLS_Y_ : 1);j++) for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) { 
      nd=getCenterNodeLocalNumber(i,j,k);

      if (startNode->block->GetCenterNode(nd)==NULL) {
        if ((i>=0)&&(i<_BLOCK_CELLS_X_)&&(j>=0)&&(j<_BLOCK_CELLS_Y_)&&(k>=0)&&(k<_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        else {
          int iNeib=0,jNeib=0,kNeib=0;

          if (i<0) iNeib=-1;
          if (i>=_BLOCK_CELLS_X_) iNeib=1;

          if (j<0) jNeib=-1;
          if (j>=_BLOCK_CELLS_Y_) jNeib=1;

          if (k<0) kNeib=-1;
          if (k>=_BLOCK_CELLS_Z_) kNeib=1;

          if (getNeibNode(iNeib,jNeib,kNeib,startNode)!=NULL) {
            exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
          }
          else continue;
        }
      }
    }


    //check if the internal nodes of the block are used only ones on the mesh
    for (k=0;k<((_MESH_DIMENSION_==3) ? _BLOCK_CELLS_Z_ : 1);k++) for (j=0;j<((_MESH_DIMENSION_>=2) ? _BLOCK_CELLS_Y_ : 1);j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      nd=getCenterNodeLocalNumber(i,j,k);
      startCenterNode=startNode->block->GetCenterNode(nd);

      if (startCenterNode==NULL) exit(__LINE__,__FILE__,"Error: a center node is not defined"); 

      if (startCenterNode->nodeDescriptor.nodeProcessedFlag==_TRUE_) exit(__LINE__,__FILE__,"Error: a center node is used twice on the mesh");
      else startCenterNode->nodeDescriptor.nodeProcessedFlag=_TRUE_; 
    }

    //check the location of the center nodes: only blocks that have the same resolution level can share 'center nodes'
     for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
      if (neibNode->RefinmentLevel==startNode->RefinmentLevel) {
        ioffset=_BLOCK_CELLS_X_*iNeibNode;
        joffset=_BLOCK_CELLS_Y_*jNeibNode;
        koffset=_BLOCK_CELLS_Z_*kNeibNode;

        for (ii=iiMin;ii<iiMax;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<jjMax;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<kkMax;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
              ndNeibNode=getCenterNodeLocalNumber(ii,jj,kk);
              ndStartNode=getCenterNodeLocalNumber(ii+ioffset,jj+joffset,kk+koffset);

              startCenterNode=startNode->block->GetCenterNode(ndStartNode);
              neibCenterNode=neibNode->block->GetCenterNode(ndNeibNode);


              if ((neibCenterNode==NULL)&&(startCenterNode==NULL)) continue;
              else if ((neibCenterNode==NULL)||(startCenterNode==NULL)) exit(__LINE__,__FILE__,"Error: one of the nodes is not defined");
              else {
                if (startCenterNode!=neibCenterNode) {
                	long int nd;
                	int iii[3]={0,0,0},idim;

                	//recalculate the position of the neibBlock
                	neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode);

                	//recalculate the indexes of the 'center nodes'
                	nd=findCenterNodeIndex(startCenterNode->GetX(),iii[0],iii[1],iii[2],startNode);
                	cout << "The center node is not found:x=";
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << startCenterNode->GetX()[idim] << "  ";
                	cout << "; (i,j,k,nd)=";
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << iii[idim] << "  ";
                	cout << nd << endl;

                	//recalculate the indexes of the 'center nodes'
                	nd=findCenterNodeIndex(neibCenterNode->GetX(),iii[0],iii[1],iii[2],neibNode);
                	cout << "The center node is not found:x=";
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << neibCenterNode->GetX()[idim] << "  ";
                	cout << "; (i,j,k,nd)=";
                	for (idim=0;idim<_MESH_DIMENSION_;idim++) cout << iii[idim] << "  ";
                	cout << nd << endl;


                	exit(__LINE__,__FILE__,"Error: two different center nodes have found that has the same coordinate 'x'");
                }



              }
            }
          }
        }

      }
    }

    #endif

  }
  else {
    //check the downBlocks
    for (kDownBlock=0;kDownBlock<=kDownNodeMax;kDownBlock++) for (jDownBlock=0;jDownBlock<=jDownNodeMax;jDownBlock++) for (iDownBlock=0;iDownBlock<=iDownNodeMax;iDownBlock++) {
      if ((downNode=startNode->downNode[iDownBlock+2*(jDownBlock+2*kDownBlock)])!=NULL) checkMeshConsistency(downNode);  
    }
  }
}

  

//=================================================================================



/*

void collectneibCornerNodes_1D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int nface,ii,nDownNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode; 
  long int nd;

  //collect the nodes from the all upNode->downNodes[:]  
  static const int iiOffset[2]={0,_TOTAL_BLOCK_CELLS_X_}; 
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_; 
  
  upNode=startNode->upNode;

  for (nDownNode=0,nDownNode<2;nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    for (ii=iiMin;ii<=iiMax;ii++) {
      nd=getCornerNodeLocalNumber(ii,0,0);
      newCornerNodeMap[ii+iiOffset[nDownNode]-iiMin][0][0]=downNode->cornerNodes[nd];    
    }
  }

  //search the neibours
  static int iiNeibOffset[2]=[-2*_BLOCK_CELLS_X_,2*_BLOCK_CELLS_X_];

  for (nface=0;nface=2*_MESH_DIMENSION_;nface++) if ((neibNode=upNode->neibNode[nface])!=NULL) {
    if (neibNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      int ioffset=iiNeibOffset[nface]; 

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) { 
        nd=getCornerNodeLocalNumber(ii,0,0);
        newCornerNodeMap[2*ii+ioffset-iiMin][0][0]=downNode->cornerNodes[nd];
      }
    }
    else {
      int ioffset=iiNeibOffset[nface]+nDownNode*_BLOCK_CELLS_X_;

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {                        
        nd=getCornerNodeLocalNumber(ii,0,0);
        newCornerNodeMap[ii+ioffset-iiMin][0][0]=downNode->cornerNodes[nd];
      }
    }
  }
}

*/

void collectNeibCornerNodes_2D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int nface,ii,jj,iDownNode,jDownNode,ioffset,joffset;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode;
  long int nd;

  //collect the nodes from the all upNode->downNodes[:]
  static const int iiOffset[2]={0,_BLOCK_CELLS_X_};
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;

  static const int jjOffset[2]={0,_BLOCK_CELLS_Y_};
  static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;
  
  upNode=startNode->upNode;

  for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) {
    downNode=upNode->downNode[iDownNode+2*jDownNode];

    for (ii=iiMin;ii<=iiMax;ii++) for (jj=jjMin;jj<=jjMax;jj++) {
      nd=getCornerNodeLocalNumber(ii,jj,0);
      newCornerNodeMap[ii+iiOffset[iDownNode]+2*_GHOST_CELLS_X_][jj+jjOffset[jDownNode]+2*_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[nd];
    }
  }

  //search the neibours
  static const int _X_FACE_NEIB_OFFSET_[4]={-2*_BLOCK_CELLS_X_,0,2*_BLOCK_CELLS_X_,0};
  static const int _Y_FACE_NEIB_OFFSET_[4]={0,-2*_BLOCK_CELLS_Y_,0,2*_BLOCK_CELLS_Y_};

  for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=upNode->neibNode[nface])!=NULL) { 
    if (neibNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      ioffset=_X_FACE_NEIB_OFFSET_[nface],joffset=_Y_FACE_NEIB_OFFSET_[nface]; 

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
        for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=2*jj+joffset)&&(2*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
          nd=getCornerNodeLocalNumber(ii,jj,0);
          newCornerNodeMap[2*ii+ioffset+2*_GHOST_CELLS_X_][2*jj+joffset+2*_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[nd];
        }
      }
    }
    else {
      for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) {
        downNode=neibNode->downNode[iDownNode+2*jDownNode];
        ioffset=_X_FACE_NEIB_OFFSET_[nface]+iDownNode*_BLOCK_CELLS_X_;
        joffset=_Y_FACE_NEIB_OFFSET_[nface]+jDownNode*_BLOCK_CELLS_Y_; 

        for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
            nd=getCornerNodeLocalNumber(ii,jj,0);
            newCornerNodeMap[ii+ioffset+2*_GHOST_CELLS_X_][jj+joffset+2*_GHOST_CELLS_Y_][0]=downNode->block->cornerNodes[nd];
          }
        }
      }
    }
  } 

  
  //search the corner neibours 
  int stepFactor,iNeibNode,jNeibNode; 

  for (iNeibNode=-1;iNeibNode<2;iNeibNode+=2) for (jNeibNode=-1;jNeibNode<2;jNeibNode+=2) if (abs(iNeibNode)+abs(jNeibNode)>1) if ((neibNode=getNeibNode_2D(iNeibNode,jNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      stepFactor=2;
      ioffset=_BLOCK_CELLS_X_*((iNeibNode==-1) ? -2 : 2);
      joffset=_BLOCK_CELLS_Y_*((jNeibNode==-1) ? -2 : 2);
    }
    else {
      stepFactor=1;
      ioffset=_BLOCK_CELLS_X_*((iNeibNode==-1) ? -1 : 2);
      joffset=_BLOCK_CELLS_Y_*((jNeibNode==-1) ? -1 : 2);
    } 

    for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=stepFactor*ii+ioffset)&&(stepFactor*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
      for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=stepFactor*jj+joffset)&&(stepFactor*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
        nd=getCornerNodeLocalNumber(ii,jj,0);
        newCornerNodeMap[stepFactor*ii+ioffset+2*_GHOST_CELLS_X_][stepFactor*jj+joffset+2*_GHOST_CELLS_Y_][0]=neibNode->block->cornerNodes[nd];
      }
    }
  }

  
}


void collectNeibCornerNodes_3D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)]) {
  int ii,jj,kk,iDownNode,jDownNode,kDownNode,ioffset,joffset,koffset;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode;
  long int nd;

  for (ii=0;ii<1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_);ii++) for (jj=0;jj<1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_);jj++) for (kk=0;kk<1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_);kk++) newCornerNodeMap[ii][jj][kk]=NULL;

  //collect the nodes from the all upNode->downNodes[:]
  static const int iiOffset[2]={0,_BLOCK_CELLS_X_};
  static const int iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_;

  static const int jjOffset[2]={0,_BLOCK_CELLS_Y_};
  static const int jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;

  static const int kkOffset[2]={0,_BLOCK_CELLS_Z_};
  static const int kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;

  upNode=startNode->upNode;

  for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++) for (kDownNode=0;kDownNode<2;kDownNode++) {
    downNode=upNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

    for (ii=iiMin;ii<=iiMax;ii++) for (jj=jjMin;jj<=jjMax;jj++) for (kk=kkMin;kk<=kkMax;kk++) {
      nd=getCornerNodeLocalNumber(ii,jj,kk);
      newCornerNodeMap[ii+iiOffset[iDownNode]+2*_GHOST_CELLS_X_][jj+jjOffset[jDownNode]+2*_GHOST_CELLS_Y_][kk+kkOffset[kDownNode]+2*_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(nd);
    }
  }

  //search the neibours
  int iNeibNode,jNeibNode,kNeibNode;

  for (iNeibNode=-1;iNeibNode<2;iNeibNode++) for (jNeibNode=-1;jNeibNode<2;jNeibNode++) for (kNeibNode=-1;kNeibNode<2;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      ioffset=_BLOCK_CELLS_X_*2*iNeibNode;
      joffset=_BLOCK_CELLS_Y_*2*jNeibNode;
      koffset=_BLOCK_CELLS_Z_*2*kNeibNode;

      for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=2*ii+ioffset)&&(2*ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
        for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=2*jj+joffset)&&(2*jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
          for (kk=kkMin;kk<=kkMax;kk++) if ((-2*_GHOST_CELLS_Z_<=2*kk+koffset)&&(2*kk+koffset<=2*_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)) {
            nd=getCornerNodeLocalNumber(ii,jj,kk);
            newCornerNodeMap[2*ii+ioffset+2*_GHOST_CELLS_X_][2*jj+joffset+2*_GHOST_CELLS_Y_][2*kk+koffset+2*_GHOST_CELLS_Z_]=neibNode->block->GetCornerNode(nd);
          }
        }
      }
    }
    else {
      neibNode=neibNode->upNode;

      for (iDownNode=0;iDownNode<2;iDownNode++) for (jDownNode=0;jDownNode<2;jDownNode++)  for (kDownNode=0;kDownNode<2;kDownNode++) {
        downNode=neibNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)];

        ioffset=_BLOCK_CELLS_X_*(iDownNode+2*iNeibNode);
        joffset=_BLOCK_CELLS_Y_*(jDownNode+2*jNeibNode);
        koffset=_BLOCK_CELLS_Z_*(kDownNode+2*kNeibNode);

        for (ii=iiMin;ii<=iiMax;ii++) if ((-2*_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<=2*_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)) {
          for (jj=jjMin;jj<=jjMax;jj++) if ((-2*_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<=2*_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)) {
            for (kk=kkMin;kk<=kkMax;kk++) if ((-2*_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<=2*_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)) {
              nd=getCornerNodeLocalNumber(ii,jj,kk);
              newCornerNodeMap[ii+ioffset+2*_GHOST_CELLS_X_][jj+joffset+2*_GHOST_CELLS_Y_][kk+koffset+2*_GHOST_CELLS_Z_]=downNode->block->GetCornerNode(nd);
            }
          }
        }
      }
    }
  }


}


void collectNeibCenterNodes_3D_deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode,cCenterNode *newCenterNodeMap[_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_][_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_][_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_]) {
  int ii,jj,kk,ioffset,joffset,koffset,iNeibNode,jNeibNode,kNeibNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*neibNode;
  long int nd;

  for (ii=0;ii<_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_;ii++) for (jj=0;jj<_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_;jj++) for (kk=0;kk<_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_;kk++) newCenterNodeMap[ii][jj][kk]=NULL; 

  upNode=startNode->upNode;

  for (iNeibNode=-1;iNeibNode<2;iNeibNode++) for (jNeibNode=-1;jNeibNode<2;jNeibNode++) for (kNeibNode=-1;kNeibNode<2;kNeibNode++)  if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) if ((neibNode=getNeibNode(iNeibNode,jNeibNode,kNeibNode,upNode))!=NULL) {
    if (neibNode->RefinmentLevel==upNode->RefinmentLevel) {
      ioffset=_BLOCK_CELLS_X_*iNeibNode;
      joffset=_BLOCK_CELLS_Y_*jNeibNode;
      koffset=_BLOCK_CELLS_Z_*kNeibNode;

      for (ii=-_GHOST_CELLS_X_;ii<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;ii++) if ((-_GHOST_CELLS_X_<=ii+ioffset)&&(ii+ioffset<_BLOCK_CELLS_X_+_GHOST_CELLS_X_)) {
        for (jj=-_GHOST_CELLS_Y_;jj<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;jj++) if ((-_GHOST_CELLS_Y_<=jj+joffset)&&(jj+joffset<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)) {
          for (kk=-_GHOST_CELLS_Z_;kk<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;kk++) if ((-_GHOST_CELLS_Z_<=kk+koffset)&&(kk+koffset<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)) {
            nd=getCenterNodeLocalNumber(ii,jj,kk);
            newCenterNodeMap[ii+ioffset+_GHOST_CELLS_X_][jj+joffset+_GHOST_CELLS_Y_][kk+koffset+_GHOST_CELLS_Z_]=neibNode->block->GetCenterNode(nd);
          }
        }
      }
    }
  }
}
  


bool deleteTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
  int nface,i,j,k,iMax,jMax,kMax,iMin,jMin,kMin,nDownNode;
  cTreeNodeAMR<cBlockAMR> *upNode,*downNode,*neibNode; 


//######## DEBUG #############

if (startNode->Temp_ID==1169) {
cout << __LINE__ << endl;
}

//######## END DEBUG #########


  upNode=startNode->upNode;
  if (upNode==NULL) return false;

  //reset the mesh modified flag 
  meshModifiedFlag=true;

  //check if the neib nodes need to be removed first
  //1. check the face neibbours
  //2. check the neibours connected throught edges and corners

  //check the face neibours of all downNodes of hte 'upNode' 
  for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=downNode->neibNode[nface])!=NULL) if (neibNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
      deleteTreeNode(neibNode->downNode[0]); 
    }

    if (downNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) deleteTreeNode(downNode->downNode[0]);
  

    //check the neibours connected throught edges and corners
    #if _MESH_DIMENSION_ == 1
    iMin=-1,iMax=1,jMin=0,jMax=0,kMin=0,kMax=0; 
    #elif _MESH_DIMENSION_ == 2
    iMin=-1,iMax=1,jMin=-1,jMax=1,kMin=0,kMax=0; 
    #elif _MESH_DIMENSION_ == 3 
    iMin=-1,iMax=1,jMin=-1,jMax=1,kMin=-1,kMax=1;
    #else 
    exit(__LINE__,__FILE__,"The wrong dimension");
    #endif

    for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) if (abs(i)+abs(j)+abs(k)>1) {     
      neibNode=getNeibNode(i,j,k,downNode);

      if (neibNode!=NULL) if (neibNode->RefinmentLevel>startNode->RefinmentLevel) deleteTreeNode(neibNode);
    }
  }


  //create a new block for the upNode
  cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_)];
  long int ii,iiMin,iiMax,jj,jjMin,jjMax,kk,kkMin,kkMax,nd;

  #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_ 
  upNode->block=blocks.newElement();
  #endif


  //transfer the corner nodes of the block to the upNode->block->cornerNodes[]
  #if _MESH_DIMENSION_ == 1 
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 2 
  iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkMin=0,kkMax=0;
  collectNeibCornerNodes_2D_deleteTreeNode(startNode,newCornerNodeMap); 
  #elif _MESH_DIMENSION_ == 3 
  iiMin=-_GHOST_CELLS_X_,iiMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jjMin=-_GHOST_CELLS_Y_,jjMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kkMin=-_GHOST_CELLS_Z_,kkMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
  collectNeibCornerNodes_3D_deleteTreeNode(startNode,newCornerNodeMap);
  #else 
  exit(__LINE__,__FILE__,"wrong dimension");
  #endif 

  for (kk=kkMin;kk<=kkMax;kk++) for (jj=jjMin;jj<=jjMax;jj++) for (ii=iiMin;ii<=iiMax;ii++) {  
    nd=getCornerNodeLocalNumber(ii,jj,kk);
    upNode->block->SetCornerNode(newCornerNodeMap[2*(ii-iiMin)][2*(jj-jjMin)][2*(kk-kkMin)],nd);
    if (upNode->block->GetCornerNode(nd)!=NULL) upNode->block->GetCornerNode(nd)->incrementConnectionCounter();

    //check if the node is defined
    if (upNode->block->GetCornerNode(nd)==NULL) {

      if ((ii>=0)&&(ii<=_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<=_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<=_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
      else {
        int iNeib=0,jNeib=0,kNeib=0;

        if (ii<0) iNeib=-1;
        if (ii>_BLOCK_CELLS_X_) iNeib=1;

        if (jj<0) jNeib=-1;
        if (jj>_BLOCK_CELLS_Y_) jNeib=1;

        if (kk<0) kNeib=-1;
        if (kk>_BLOCK_CELLS_Z_) kNeib=1; 

        if ((neibNode=getNeibNode(iNeib,jNeib,kNeib,upNode))!=NULL) {
          double xMissing[3]={0.0,0.0,0.0};
          cTreeNodeAMR<cBlockAMR>* xMissingNode;

          cout << "neibNode->Temp_ID=" << neibNode->Temp_ID << ", neibNode->RefinmentLevel=" << neibNode->RefinmentLevel << endl;

          //coordinates of the missing node
          cout << "Coordinates of the missing node: x=" << (xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*ii); 
          if (_MESH_DIMENSION_>=2) cout << ", " << (xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*jj);
          if (_MESH_DIMENSION_==3) cout << ", " << (xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*kk); 
          cout << endl;

          //find the node that contains xMissing and compare it with neibNode
          xMissingNode=findTreeNode(xMissing);

          if ((xMissing[0]<xMissingNode->xmin[0])||(xMissingNode->xmax[0]<xMissing[0])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");  
          if (_MESH_DIMENSION_>=2) if ((xMissing[1]<xMissingNode->xmin[1])||(xMissingNode->xmax[1]<xMissing[1])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_==3) if ((xMissing[2]<xMissingNode->xmin[2])||(xMissingNode->xmax[2]<xMissing[2])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");

          if (neibNode!=xMissingNode) { //recalculate the value of 'neibBlock' for debugging perposes
            cTreeNodeAMR<cBlockAMR>* recalculatedNeibNode;

            recalculatedNeibNode=getNeibNode(iNeib,jNeib,kNeib,upNode); 
            //exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          }
   
          //the position of the missed node in the 'xMissingNode' 
          long int ndMissingNode;
          int iMissingNode=0,jMissingNode=0,kMissingNode=0;

          ndMissingNode=findCornerNodeIndex(xMissing,iMissingNode,jMissingNode,kMissingNode,xMissingNode); 
          cout << "The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=" << iMissingNode << ", " << jMissingNode << ", " << kMissingNode << ", ndMissingNode=" << ndMissingNode << endl;


          exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        }
      }
    }
  }  

 
  //set up central nodes of the upNode->block
  #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
  cCenterNode *newCenterNode,*newCenterNodeMap[_BLOCK_CELLS_X_+2*_GHOST_CELLS_X_][_BLOCK_CELLS_Y_+2*_GHOST_CELLS_Y_][_BLOCK_CELLS_Z_+2*_GHOST_CELLS_Z_];
  double xNode[3];
  int idim;


  //###########  DEBUG  #######
  if (startNode->Temp_ID==729) {
  	cout << __LINE__ << endl;
  }

  if (upNode->Temp_ID==261) {
  	cout << __LINE__ << endl;
  }

  //########### END DEBUG ######

 
  #if _MESH_DIMENSION_ == 1
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 2
  exit(__LINE__,__FILE__,"not implemented");
  #elif _MESH_DIMENSION_ == 3 
  collectNeibCenterNodes_3D_deleteTreeNode(startNode,newCenterNodeMap); 
  #else 
  exit(__LINE__,__FILE__,"The value is not defined");
  #endif


  //create the missing nodes
  for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
    nd=getCenterNodeLocalNumber(ii,jj,kk); 

    if (newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin]==NULL) {
      newCenterNode=CenterNodes.newElement();
      newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin]=newCenterNode;

//###########  DEBUG  #######
if (newCenterNode->Temp_ID==88861) {
	cout << __LINE__ << endl;
}

//########### END DEBUG ######

      xNode[0]=upNode->xmin[0]+(ii+0.5)*dxRootBlock[0]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_X_);
      if (_MESH_DIMENSION_>=2) xNode[1]=upNode->xmin[1]+(jj+0.5)*dxRootBlock[1]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_Y_); 
      if (_MESH_DIMENSION_==3) xNode[2]=upNode->xmin[2]+(kk+0.5)*dxRootBlock[2]/(1<<upNode->RefinmentLevel)/double(_BLOCK_CELLS_Z_); 

      #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
      newCenterNode->SetX(xNode);
      #endif

      for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-EPS)||(xNode[idim]>xGlobalMax[idim]+EPS)) {
        if (newCenterNode->getNodeConnectionNumber()==0) CenterNodes.deleteElement(newCenterNode);
        newCenterNodeMap[i-iiMin][j-jjMin][k-kkMin]=NULL;
        break;
      }
    }
  }


  //transfer the node table 
  for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
    nd=getCenterNodeLocalNumber(ii,jj,kk);
    upNode->block->SetCenterNode(newCenterNodeMap[ii-iiMin][jj-jjMin][kk-kkMin],nd);
    if (upNode->block->GetCenterNode(nd)!=NULL) upNode->block->GetCenterNode(nd)->incrementConnectionCounter();

    //check if the node is defined
    if (upNode->block->GetCenterNode(nd)==NULL) {
      if ((ii>=0)&&(ii<_BLOCK_CELLS_X_)&&(jj>=0)&&(jj<_BLOCK_CELLS_Y_)&&(kk>=0)&&(kk<_BLOCK_CELLS_Z_)) exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
      else {
        int iNeib=0,jNeib=0,kNeib=0;

        if (ii<0) iNeib=-1;
        if (ii>=_BLOCK_CELLS_X_) iNeib=1;

        if (jj<0) jNeib=-1;
        if (jj>=_BLOCK_CELLS_Y_) jNeib=1;

        if (kk<0) kNeib=-1;
        if (kk>=_BLOCK_CELLS_Z_) kNeib=1;

        if ((neibNode=getNeibNode(iNeib,jNeib,kNeib,upNode))!=NULL) {
          double xMissing[3]={0.0,0.0,0.0};
          cTreeNodeAMR<cBlockAMR>* xMissingNode;

          cout << "neibNode->Temp_ID=" << neibNode->Temp_ID << ", neibNode->RefinmentLevel=" << neibNode->RefinmentLevel << endl;

          //coordinates of the missing node
          cout << "Coordinates of the missing node: x=" << (xMissing[0]=upNode->xmin[0]+(upNode->xmax[0]-upNode->xmin[0])/double(_BLOCK_CELLS_X_)*(ii+0.5));
          if (_MESH_DIMENSION_>=2) cout << ", " << (xMissing[1]=upNode->xmin[1]+(upNode->xmax[1]-upNode->xmin[1])/double(_BLOCK_CELLS_Y_)*(jj+0.5));
          if (_MESH_DIMENSION_==3) cout << ", " << (xMissing[2]=upNode->xmin[2]+(upNode->xmax[2]-upNode->xmin[2])/double(_BLOCK_CELLS_Z_)*(kk+0.5));
          cout << endl;

          //find the node that contains xMissing and compare it with neibNode
          xMissingNode=findTreeNode(xMissing);

          if ((xMissing[0]<xMissingNode->xmin[0])||(xMissingNode->xmax[0]<xMissing[0])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_>=2) if ((xMissing[1]<xMissingNode->xmin[1])||(xMissingNode->xmax[1]<xMissing[1])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          if (_MESH_DIMENSION_==3) if ((xMissing[2]<xMissingNode->xmin[2])||(xMissingNode->xmax[2]<xMissing[2])) exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");

          if (neibNode!=xMissingNode) { //recalculate the value of 'neibBlock' for debugging perposes
            cTreeNodeAMR<cBlockAMR>* recalculatedNeibNode;

            recalculatedNeibNode=getNeibNode(iNeib,jNeib,kNeib,upNode);
            //exit(__LINE__,__FILE__,"Cannot find 'xMissingNode'");
          }
  
          //the position of the missed node in the 'xMissingNode'
          long int ndMissingNode;
          int iMissingNode=0,jMissingNode=0,kMissingNode=0;

          ndMissingNode=findCenterNodeIndex(xMissing,iMissingNode,jMissingNode,kMissingNode,xMissingNode);
          cout << "The (i,j,k) index of the missing node in 'xMissingNode': (i,j,k)=" << iMissingNode << ", " << jMissingNode << ", " << kMissingNode << ", ndMissingNode=" << ndMissingNode << endl;


          exit(__LINE__,__FILE__,"Blocks's corner node is not defined");
        }
      }
    }
  }
  #endif 


  //delete the tree node
  for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
    downNode=upNode->downNode[nDownNode];

    //remove the corner nodes
    for (kk=kkMin;kk<=kkMax;kk++) for (jj=jjMin;jj<=jjMax;jj++) for (ii=iiMin;ii<=iiMax;ii++) {
      nd=getCornerNodeLocalNumber(ii,jj,kk);
      if (downNode->block->GetCornerNode(nd)->decrementConnectionCounter()==0) CornerNodes.deleteElement(downNode->block->GetCornerNode(nd));
    }

    //remove the center nodes
    #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
    for (kk=kkMin;kk<kkMax;kk++) for (jj=jjMin;jj<jjMax;jj++) for (ii=iiMin;ii<iiMax;ii++) {
      nd=getCenterNodeLocalNumber(ii,jj,kk);
      if (downNode->block->GetCenterNode(nd)->decrementConnectionCounter()==0) CenterNodes.deleteElement(downNode->block->GetCenterNode(nd));
    }
    #endif


    //remove the connections with the neibours
    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) if ((neibNode=downNode->neibNode[nface])!=NULL) {
      for (int nNeibFace=0;nNeibFace<2*_MESH_DIMENSION_;nNeibFace++) if (neibNode->neibNode[nNeibFace]==downNode) {
        neibNode->neibNode[nNeibFace]=NULL;
        break;
      }

      downNode->neibNode[nface]=NULL;  
    }

    //remove the connection with the parent tree node
    upNode->downNode[nDownNode]=NULL;

    //empty the memoty 
    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_OFF_
    blocks.deleteElement(downNode->block);
    downNode->block=NULL;
    #endif
  }

  return true; 
}


    
//=================================================================================================
 bool splitTreeNode(cTreeNodeAMR<cBlockAMR> *startNode) {
   int idim,i,j,k;
   long int nd,rLevel,nDownNode;
   double dxBlock[3];


if (startNode->Temp_ID==3) {
cout << __LINE__ << endl;
}



   if (startNode->RefinmentLevel>=_MAX_REFINMENT_LEVEL_) return false;
   meshModifiedFlag=true;

   //create the daugher blocks
   if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //there is no daugher blocks
     cTreeNodeAMR<cBlockAMR> *newTreeNode,*neibNode;
     cCornerNode *newCornerNodeMap[1+2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][1+2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][1+2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)];
     cCornerNode *newCornerNode;

     //check if the neibNodes need to be refined  
     for (i=-1;i<2;i++) for (j=-1;j<2;j++) for (k=-1;k<2;k++) if ((neibNode=getNeibNode(i,j,k,startNode))!=NULL) if (neibNode->RefinmentLevel<startNode->RefinmentLevel) { //refine the neibNode   
       splitTreeNode(neibNode);
     }

     #if _MESH_DIMENSION_ == 1
     collectNeibNodeMap_1D(startNode,newCornerNodeMap);
     #elif _MESH_DIMENSION_ == 2
     collectNeibNodeMap_2D(startNode,newCornerNodeMap); 
     #else 
     collectNeibCornerNodeMap_3D(startNode,newCornerNodeMap);
     #endif 

     //transfer the nodes from 'startNode' to 'newCornerNodeMap' and create the mission nodes 
     #if _MESH_DIMENSION_ == 1
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=0,kCellMax=0;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=0,kCellMin=0;
     #elif _MESH_DIMENSION_ == 2
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kCellMax=0;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=-_GHOST_CELLS_Y_,kCellMin=0;
     #else
     static const int iCellMax=2*_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jCellMax=2*_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kCellMax=2*_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
     static const int iCellMin=-_GHOST_CELLS_X_,jCellMin=-_GHOST_CELLS_Y_,kCellMin=-_GHOST_CELLS_Z_; 
     #endif

      rLevel=startNode->RefinmentLevel+1;
      dxBlock[0]=(xGlobalMax[0]-xGlobalMin[0])/(1<<rLevel)/_BLOCK_CELLS_X_;
      if (_MESH_DIMENSION_>1) dxBlock[1]=(xGlobalMax[1]-xGlobalMin[1])/(1<<rLevel)/_BLOCK_CELLS_Y_;
      if (_MESH_DIMENSION_>2) dxBlock[2]=(xGlobalMax[2]-xGlobalMin[2])/(1<<rLevel)/_BLOCK_CELLS_Z_;

      for (k=kCellMin;k<=kCellMax;k++) for (j=jCellMin;j<=jCellMax;j++) for (i=iCellMin;i<=iCellMax;i++) if (newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) {
        //check if the new node has to be created or transfered from the startBlock
        if ((2*(i/2)!=i)||(2*(j/2)!=j)||(2*(k/2)!=k)) { // create new nodes
          if (newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) {
            newCornerNode=CornerNodes.newElement();
            newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=newCornerNode;  

            double xNode[3];

            xNode[0]=startNode->xmin[0]+i*dxBlock[0];
            if (_MESH_DIMENSION_>=2) xNode[1]=startNode->xmin[1]+j*dxBlock[1];
            if (_MESH_DIMENSION_==3) xNode[2]=startNode->xmin[2]+k*dxBlock[2];

            newCornerNode->SetX(xNode);

            for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-0.001*dxBlock[idim])||(xNode[idim]>xGlobalMax[idim]+0.001*dxBlock[idim])) {
              if (newCornerNode->getNodeConnectionNumber()==0) CornerNodes.deleteElement(newCornerNode);
              newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=NULL; 
              break;
            }
          }
        }
        else {
          nd=getCornerNodeLocalNumber(i/2,j/2,k/2);
          newCornerNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=startNode->block->GetCornerNode(nd);
        }
      }

      //collect and create the block's center nodes
      #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_ 
      cCenterNode *newCenterNodeMap[2*(_BLOCK_CELLS_X_+_GHOST_CELLS_X_)][2*(_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_)][2*(_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_)];
      cCenterNode *newCenterNode;

      #if _MESH_DIMENSION_ == 1
      exit(__LINE__,__FILE__,"not implemented"); 
      #elif _MESH_DIMENSION_ == 2
      exit(__LINE__,__FILE__,"not implemented");
      #else
      collectNeibCenterNodeMap_3D(startNode,newCenterNodeMap);
      #endif 

      //create new 'center nodes' when its needed
      double xNode[3];

      for (k=kCellMin;k<kCellMax;k++) for (j=jCellMin;j<jCellMax;j++) for (i=iCellMin;i<iCellMax;i++) if (newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]==NULL) { //create bew new node 
        newCenterNode=CenterNodes.newElement();
        newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=newCenterNode;


//################  DEBUG ##############

if (newCenterNode->Temp_ID==88861) {
	cout << __LINE__ << endl;
}

//################ END DEBUG ###########


        xNode[0]=startNode->xmin[0]+(i+0.5)*dxBlock[0];
        if (_MESH_DIMENSION_>=2) xNode[1]=startNode->xmin[1]+(j+0.5)*dxBlock[1];
        if (_MESH_DIMENSION_==3) xNode[2]=startNode->xmin[2]+(k+0.5)*dxBlock[2];

        #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
        newCenterNode->SetX(xNode);
        #endif

        for (idim=0;idim<_MESH_DIMENSION_;idim++) if ((xNode[idim]<xGlobalMin[idim]-0.001*dxBlock[idim])||(xNode[idim]>xGlobalMax[idim]+0.001*dxBlock[idim])) {
          if (newCenterNode->getNodeConnectionNumber()==0) CenterNodes.deleteElement(newCenterNode);
          newCenterNodeMap[i-iCellMin][j-jCellMin][k-kCellMin]=NULL;
          break;
        }
      } 
      #endif


      //create the downNodes
      int iOffset,jOffset,kOffset,kDownNode,jDownNode,iDownNode;

      #if _MESH_DIMENSION_ == 1
      static const int iDownNodeMax=2,jDownNodeMax=1,kDownNodeMax=1;
      #elif _MESH_DIMENSION_ == 2
      static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=1;
      #else
      static const int iDownNodeMax=2,jDownNodeMax=2,kDownNodeMax=2;
      static const int iMin=-_GHOST_CELLS_X_,iMax=_BLOCK_CELLS_X_+_GHOST_CELLS_X_,jMin=-_GHOST_CELLS_Y_,jMax=_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_,kMin=-_GHOST_CELLS_Z_,kMax=_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;
      #endif


      for (kDownNode=0;kDownNode<kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<iDownNodeMax;iDownNode++) if (startNode->GetDownNode(iDownNode,jDownNode,kDownNode)==NULL) {
        newTreeNode=treeNodes.newElement();
        startNode->SetDownNode(newTreeNode,iDownNode,jDownNode,kDownNode);

        //set block into the tree
        newTreeNode->block=blocks.newElement();
        newTreeNode->block->SetRefinmentLevel(rLevel);
        newTreeNode->RefinmentLevel=rLevel; 

        //upBlock
        newTreeNode->upNode=startNode;

        //init the corner nodes of the new blocks
        iOffset=iDownNode*_BLOCK_CELLS_X_,jOffset=jDownNode*_BLOCK_CELLS_Y_,kOffset=kDownNode*_BLOCK_CELLS_Z_;

        for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) { 
          nd=getCornerNodeLocalNumber(i,j,k); 

          newTreeNode->block->SetCornerNode(newCornerNodeMap[iOffset+i-iMin][j+jOffset-jMin][k+kOffset-kMin],nd);
          if (newTreeNode->block->GetCornerNode(nd)!=NULL) newTreeNode->block->GetCornerNode(nd)->incrementConnectionCounter();
        }

        //init the center nodes of the block 
        #if _AMR_CENTER_NODE_ == _ON_AMR_MESH_
        for (k=kMin;k<kMax;k++) for (j=jMin;j<jMax;j++) for (i=iMin;i<iMax;i++) {
          nd=getCenterNodeLocalNumber(i,j,k);

          newTreeNode->block->SetCenterNode(newCenterNodeMap[iOffset+i-iMin][j+jOffset-jMin][k+kOffset-kMin],nd);
          if (newTreeNode->block->GetCenterNode(nd)!=NULL) newTreeNode->block->GetCenterNode(nd)->incrementConnectionCounter();
        }
        #endif 
  


        //init newTreeNode->xmin
        double *x,*xmin,*xmax;

        nd=getCornerNodeLocalNumber(0,0,0);
        xmin=newTreeNode->xmin; 
        x=newTreeNode->block->GetCornerNode(nd)->GetX();
        for (idim=0;idim<_MESH_DIMENSION_;idim++) xmin[idim]=x[idim]; 

        nd=getCornerNodeLocalNumber(_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_);
        xmax=newTreeNode->xmax;
        x=newTreeNode->block->GetCornerNode(nd)->GetX();
        for (idim=0;idim<_MESH_DIMENSION_;idim++) xmax[idim]=x[idim];
      }



      //connect the newly created elements
      #if _MESH_DIMENSION_ == 1
      startNode->downNode[0]->neibNode[1]=startNode->downNode[1];
      startNode->downNode[1]->neibNode[0]=startNode->downNode[0];
      #else

      #if _MESH_DIMENSION_ == 2

      //conenction table
      static const int innerNodeConnectionMap[4][4]={ {-1,-1,1,2}, {0,-1,-1,3}, {-1,0,3,-1}, {2,1,-1,-1}};

      #elif _MESH_DIMENSION_ == 3

      static const int innerNodeConnectionMap[8][6]={ {-1,1,-1,2,-1,4}, {0,-1,-1,3,-1,5}, {-1,3,0,-1,-1,6}, {2,-1,1,-1,-1,7},
                                                      {-1,5,-1,6,0,-1}, {4,-1,-1,7,1,-1}, {-1,7,4,-1,2,-1,}, {6,-1,5,-1,3,-1}   };

      #endif

      for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) for (i=0;i<2*_MESH_DIMENSION_;i++) {
        startNode->downNode[nDownNode]->neibNode[i]=(innerNodeConnectionMap[nDownNode][i]!=-1) ? startNode->downNode[innerNodeConnectionMap[nDownNode][i]] : NULL;
      }
      #endif

      reconnectDownTreeNode(startNode);
   }

   return true;
}




  //reset the node's 'nodeProcessedFlag' 
  void resetNodeProcessedFlag(bool resetMaxRefinmentLevel=true) {
    long int nMemoryBank,nTotalMemoryBanks,nnode;


    //reset the flag for the 'corner nodes'
    cCornerNode *cornerNodeDataBuffer;
    nTotalMemoryBanks=CornerNodes.dataBufferListPointer;   

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      cornerNodeDataBuffer=CornerNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (cornerNodeDataBuffer+nnode)->nodeDescriptor.nodeProcessedFlag=_FALSE_;
        if (resetMaxRefinmentLevel==true) (cornerNodeDataBuffer+nnode)->nodeDescriptor.maxRefinmentLevel=0;
      }
    }

    //reset the flag for the 'center nodes'
    #if _AMR_CENTER_NODE_ ==  _ON_AMR_MESH_
    //reset the flag for the 'corner nodes'
    cCenterNode *centerNodeDataBuffer;
    nTotalMemoryBanks=CenterNodes.dataBufferListPointer;

    for (nMemoryBank=0;nMemoryBank<nTotalMemoryBanks;nMemoryBank++) {
      centerNodeDataBuffer=CenterNodes.dataBufferList[nMemoryBank];

      for (nnode=0;nnode<_STACK_DEFAULT_BUFFER_BUNK_SIZE_;nnode++) {
        (centerNodeDataBuffer+nnode)->nodeDescriptor.nodeProcessedFlag=_FALSE_;
        if (resetMaxRefinmentLevel==true) (centerNodeDataBuffer+nnode)->nodeDescriptor.maxRefinmentLevel=0;
      }
    }
    #endif


  } 




  void countMeshElements(cTreeNodeAMR<cBlockAMR> *startNode,int level) {
    static long int nDownNodes,nodeno,nd,iMax,jMax,kMax;
    int i,j,k;
    
    if (startNode==rootTree) {
      meshBlocksNumber=0,meshNodesNumber=0,nodeno=1,meshMaximumRefinmentLevel=0;
      resetNodeProcessedFlag();

      #if _MESH_DIMENSION_ == 1
      iMax=1+_BLOCK_CELLS_X_,jMax=1,kMax=1; 
      #elif _MESH_DIMENSION_ == 2
      iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1;
      #else 
      iMax=1+_BLOCK_CELLS_X_,jMax=1+_BLOCK_CELLS_Y_,kMax=1+_BLOCK_CELLS_Z_;
      #endif

      nDownNodes=1<<_MESH_DIMENSION_;
    }

    //check the value of the maximum refinment level
    if (level>meshMaximumRefinmentLevel) meshMaximumRefinmentLevel=level;

    //count the number of blocks on the mesh 
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) meshBlocksNumber++; 

    //count the mesh nodes associated with the block
    if (startNode->block!=NULL) {
      cBlockAMR *block=startNode->block;
      cCornerNode *ndptr;

      for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
        nd=getCornerNodeLocalNumber(i,j,k);
        if ((ndptr=block->GetCornerNode(nd))==NULL) continue;
 
        if (ndptr->nodeDescriptor.nodeProcessedFlag==_FALSE_) {
          ndptr->nodeDescriptor.nodeProcessedFlag=_TRUE_;
          ndptr->nodeDescriptor.nodeno=nodeno;

          nodeno++;
          meshNodesNumber++;
        }

        if (ndptr->nodeDescriptor.maxRefinmentLevel<level) ndptr->nodeDescriptor.maxRefinmentLevel=level;
      }
    }

    //cout the nodes in the downNodes 
    for (i=0;i<nDownNodes;i++) if (startNode->downNode[i]!=NULL) countMeshElements(startNode->downNode[i],level+1);


    //reset the 'meshModifiedFlag'
    meshModifiedFlag=false;
  } 


  long int getMeshNodesNumber() {
    if (meshModifiedFlag==true) countMeshElements(rootTree,0); 

    return meshNodesNumber;
  }  

  int getMeshMaxRefinmentLevel() {
    if (meshModifiedFlag==true) countMeshElements(rootTree,0);

    return meshMaximumRefinmentLevel;
  } 


  //build the mesh
  bool buildMesh_OneLevelRefinment(cTreeNodeAMR<cBlockAMR> *startNode,int level,int startLevel) {
    bool res=false,flag;
    double c,blockMiddlePoint[3]={0.0,0.0,0.0},requredResolution,characteristicBlockSize;
    int idim,nnode,nBlockNodes;
    cCornerNode *blockNode,**nodeList;

    #if _MESH_DIMENSION_ == 1 
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2 
    static const double characteristicBlockSize_max=sqrt(pow(dxRootBlock[0]/_TOTAL_BLOCK_CELLS_X_,2)+pow(dxRootBlock[1]/_TOTAL_BLOCK_CELLS_Y_,2));
    static const double characteristicBlockSize_min=characteristicBlockSize_max/(1<<_MAX_REFINMENT_LEVEL_);
    #else 
    static const double characteristicBlockSize_max=sqrt(pow(dxRootBlock[0]/_TOTAL_BLOCK_CELLS_X_,2)+pow(dxRootBlock[1]/_TOTAL_BLOCK_CELLS_Y_,2)+pow(dxRootBlock[2]/_TOTAL_BLOCK_CELLS_Z_,2));
    static const double characteristicBlockSize_min=characteristicBlockSize_max/(1<<_MAX_REFINMENT_LEVEL_);
    #endif 


    if (level==startLevel) { //refine the block if needed
      //get the characterist size of the block and the position of the block's middle point 

      blockNode=startNode->block->GetCornerNode(getCornerNodeLocalNumber(0,0,0));
      characteristicBlockSize=characteristicBlockSize_max/(1<<level);

      for (idim=0;idim<_MESH_DIMENSION_;idim++) blockMiddlePoint[idim]=blockNode->GetX()[idim]+0.5*(xGlobalMax[idim]-xGlobalMin[idim])/(1<<level);

      //get the requested mesh resolution 
      requredResolution=localResolution(blockMiddlePoint);

      #if _MESH_DIMENSION_ == 1
      nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
      #elif _MESH_DIMENSION_ == 2
      nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
      #else
      nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
      #endif

      for (nnode=0,nodeList=startNode->block->GetCornerNodeBuffer();nnode<nBlockNodes;nnode++,nodeList++) if (*nodeList!=NULL) {
        c=localResolution((*nodeList)->GetX());
        if (c<requredResolution) requredResolution=c;
        if (c<characteristicBlockSize_min) exit(__LINE__,__FILE__,"The required resolution is smaller than the minimum resolution allowd for the mesh. Increase the value of _MAX_REFINMENT_LEVEL_"); 
      }  

      if (requredResolution<characteristicBlockSize) {
        res=splitTreeNode(startNode);
      } 
    }
    else { // go to the downNodes  
      int nDownNode;

      for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) {
        flag=buildMesh_OneLevelRefinment(startNode->downNode[nDownNode],level+1,startLevel);
        if (flag==true) res=true;
      }
    }

    return res;
  } 

  void buildMesh() {
    int level;
    bool flag;
 
    for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) {
      flag=buildMesh_OneLevelRefinment(rootTree,0,level);

      #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_ 
      #if _CHECK_MESH_CONSISTANCY_ == _ON_AMR_MESH_
      checkMeshConsistency(rootTree); 
      #endif
      #endif

      if (flag==false) break;
    }
  }  
  
//==============================================================
//calcualte the interpolation coefficients
  int CenterNodesInterpolationCoefficients_3D_linear(double *x,double *CoefficientsList,cCenterNode **InterpolationStencil,cTreeNodeAMR<cBlockAMR>* startNode,int nMaxCoefficients) {

    //if the length of the coefficient list is not enough -> exist with an error message
    if (nMaxCoefficients<8) {
	  exit(__LINE__,__FILE__,"The length of the interpolation stencil is too short");
  	  return -1;
    }

    //determine the local coordinates of the point
    double xLocal,yLocal,zLocal,*xmin=startNode->xmin,*xmax=startNode->xmax;

    xLocal=_BLOCK_CELLS_X_*(x[0]-xmin[0])/(xmax[0]-xmin[0]);
    yLocal=_BLOCK_CELLS_Y_*(x[1]-xmin[1])/(xmax[1]-xmin[1]);
    zLocal=_BLOCK_CELLS_Z_*(x[2]-xmin[2])/(xmax[2]-xmin[2]);

    if ((xLocal<0.0)||(xLocal>_BLOCK_CELLS_X_)||(yLocal<0.0)||(yLocal>_BLOCK_CELLS_Y_)||(zLocal<0.0)||(zLocal>_BLOCK_CELLS_Z_)) {
      exit(__LINE__,__FILE__,"The point is outside of the block");
      return -1;
    }

    //determine the offset of the interpolating stencil related to the origin of the block
    int i,j,k,ioffset,joffset,koffset,counter;
    long int nd;
    cCenterNode *ptr;

    ioffset=int(xLocal+0.5)-1;
    joffset=int(yLocal+0.5)-1;
    koffset=int(zLocal+0.5)-1;

    //calculate the local coordinated of 'x' within the stencil
    double c,cx,cy,dx,dy,dz;
    double cTotal=0.0;

    dx=xLocal+0.5-(ioffset+1);
    dy=yLocal+0.5-(joffset+1);
    dz=zLocal+0.5-(koffset+1);

    for (i=0,counter=0;i<2;i++) {
      cx=(i==0) ? 1.0-dx : dx;

      for (j=0;j<2;j++) {
    	cy=(j==0) ? 1.0-dy : dy;

    	for (k=0;k<2;k++) {
    	  c=((k==0) ? 1.0-dz : dz)*cx*cy;
    	  nd=getCenterNodeLocalNumber(i+ioffset,j+joffset,k+koffset);
          ptr=startNode->block->GetCenterNode(nd);

    	  if (ptr!=NULL) {
    	    CoefficientsList[counter]=c;
    	    InterpolationStencil[counter]=ptr;
    	    counter+=1;
    	    cTotal+=c;
          }
    	}
      }
    }

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    if ((counter==8)&&(fabs(cTotal-1.0)>1.0E-6)) exit(__LINE__,__FILE__,"The summ of the coefficients is different from one");
    #endif

    if (counter==0) exit(__LINE__,__FILE__,"There is no nodes defined");

    if (cTotal<1.0) for (int i=0;i<counter;i++) CoefficientsList[i]/=cTotal;

    return counter;
  }

//==============================================================
  //if printCoordinateVector == true  -> print the coordinate vector, else -> printf the connectovity list  
  void outputMeshTECPLOT_BlockCornerNode(cTreeNodeAMR<cBlockAMR> *startNode,FILE* fout,bool printCoordinateVector,bool PrintMeshData,int DataSetNumber) {
    int isubBlock,jsubBlock,ksubBlock,isubBlockMax,jsubBlockMax,ksubBlockMax;
    int iNode,jNode,kNode,nd,idim,nnode,nBasicBlockNodes;  

    static long int nGlobalNodeNumber=-1;

    if (startNode==rootTree) nGlobalNodeNumber=1;
  
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //the node is on the bottom of the tree -> printf the node
      cBlockAMR *block=startNode->block;
      cCornerNode *cornerNode;

      //get the limits for the number of the sub-blocks in the block
      #if _MESH_DIMENSION_ == 1
      isubBlockMax=_BLOCK_CELLS_X_; 
      jsubBlockMax=1,ksubBlockMax=1;
      nBasicBlockNodes=2; 
      static int nodeOutputOrder_I[2]={0,1};
      static int nodeOutputOrder_J[2]={0,0};
      static int nodeOutputOrder_K[2]={0,0};

      #elif _MESH_DIMENSION_ == 2
      isubBlockMax=_BLOCK_CELLS_X_;
      jsubBlockMax=_BLOCK_CELLS_Y_; 
      ksubBlockMax=1;
      nBasicBlockNodes=4; 
      static int nodeOutputOrder_I[4]={0,1,1,0};
      static int nodeOutputOrder_J[4]={0,0,1,1};
      static int nodeOutputOrder_K[4]={0,0,0,0};
      #else
      isubBlockMax=_BLOCK_CELLS_X_;
      jsubBlockMax=_BLOCK_CELLS_Y_;
      ksubBlockMax=_BLOCK_CELLS_Z_;
      nBasicBlockNodes=8;
      static int nodeOutputOrder_I[8]={0,1,1,0,0,1,1,0};
      static int nodeOutputOrder_J[8]={0,0,1,1,0,0,1,1};
      static int nodeOutputOrder_K[8]={0,0,0,0,1,1,1,1};
      #endif


      for (ksubBlock=0;ksubBlock<ksubBlockMax;ksubBlock++) for (jsubBlock=0;jsubBlock<jsubBlockMax;jsubBlock++) for (isubBlock=0;isubBlock<isubBlockMax;isubBlock++) { 
        for (nnode=0;nnode<nBasicBlockNodes;nnode++) {
          iNode=nodeOutputOrder_I[nnode]+isubBlock;
          jNode=nodeOutputOrder_J[nnode]+jsubBlock;
          kNode=nodeOutputOrder_K[nnode]+ksubBlock;

          nd=getCornerNodeLocalNumber(iNode,jNode,kNode);
          cornerNode=block->GetCornerNode(nd);

          if (printCoordinateVector==true) {
            if ((cornerNode->nodeDescriptor.nodeProcessedFlag==_FALSE_)&&(cornerNode->nodeDescriptor.maxRefinmentLevel==startNode->RefinmentLevel)) {
              cornerNode->nodeDescriptor.nodeProcessedFlag=_TRUE_;
              cornerNode->nodeDescriptor.nodeno=nGlobalNodeNumber++;

              for (idim=0;idim<_MESH_DIMENSION_;idim++) fprintf(fout,"%e  ",((fabs(cornerNode->GetX()[idim])>EPS) ? cornerNode->GetX()[idim] : 0.0));
              fprintf(fout,"%i  %i  %ld ",cornerNode->nodeDescriptor.maxRefinmentLevel, ((startNode->block->GetGhostFlag()==_GHOST_BLOCK_) ? 1 : 0),cornerNode->Temp_ID);

              //print the data stared in the block's nodes
              if (PrintMeshData==true) {

                //print the data stared in the corner nodes
                cornerNode->PrintData(fout,DataSetNumber);

                //print the data stored in the 'center' nodes
                #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
                const int nMaxCenterInterpolationCoefficients=8;
                cCenterNode *tempCenterNode,*CenterNodeInterpolationStencil[nMaxCenterInterpolationCoefficients];
                double CenterNodeInterpolationCoefficients[nMaxCenterInterpolationCoefficients],xCenterNode[3];
                int centerNodeInterpolationStencilLength;

                fprintf(fout,"  ");
                cornerNode->GetX(xCenterNode);
                tempCenterNode=CenterNodes.newElement();

                if (GetCenterNodesInterpolationCoefficients==NULL) {
                  centerNodeInterpolationStencilLength=CenterNodesInterpolationCoefficients_3D_linear(xCenterNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
                }
                else {
                  centerNodeInterpolationStencilLength=GetCenterNodesInterpolationCoefficients(xCenterNode,CenterNodeInterpolationCoefficients,CenterNodeInterpolationStencil,startNode,nMaxCenterInterpolationCoefficients);
                }

                if (centerNodeInterpolationStencilLength==-1) exit(__LINE__,__FILE__,"Errror in interpolation of the 'center' node's data");

                tempCenterNode->Interpolate(CenterNodeInterpolationStencil,CenterNodeInterpolationCoefficients,centerNodeInterpolationStencilLength);
                tempCenterNode->PrintData(fout,DataSetNumber);
                #endif

                CenterNodes.deleteElement(tempCenterNode);

                //print data stored by the tree node (startNode)
                startNode->block->PrintData(fout,DataSetNumber);

              }

              //finish the line
              fprintf(fout,"\n");         





            }
          }
          else fprintf(fout,"%ld ",cornerNode->nodeDescriptor.nodeno);
        }

        if (printCoordinateVector==false) fprintf(fout,"\n");
      }
    }
    else {
      for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) outputMeshTECPLOT_BlockCornerNode(startNode->downNode[nDownNode],fout,printCoordinateVector,PrintMeshData,DataSetNumber);
    }
  } 
 


  void outputMeshTECPLOTinternal(const char *fname,const bool PrintMeshData, int DataSetNumber) {
    FILE *fout;

    //cout the number of the elements of the mesh
    if (meshModifiedFlag==true) countMeshElements(rootTree,0);

    fout=fopen(fname,"w");

    #if _MESH_DIMENSION_ == 1 
    fprintf(fout,"VARIABLES=\"X\", \"Maximum Refinment Level\"");

    if (PrintMeshData==true) {
      CornerNodes.elementStackList[0][0]->PrintVariableList(fout);

      #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
      CenterNodes.elementStackList[0][0]->PrintVariableList(fout);
      #endif

      rootTree->block->PrintVariableList(fout);
    }

    fprintf(fout,"\n");

    #elif _MESH_DIMENSION_ == 2  
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Maximum Refinment Level\", \"Ghost cells indicator\", \"Temp_ID\" \n");

    if (PrintMeshData==true) {
      CornerNodes.elementStackList[0][0]->PrintVariableList(fout);

      #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
      CenterNodes.elementStackList[0][0]->PrintVariableList(fout);
      #endif

      rootTree->block->PrintVariableList(fout);
    }

    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",meshNodesNumber,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_);
    #else
    fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\", \"Maximum Refinment Level\", \"Ghost cells indicator\", \"Temp_ID\"");

    if (PrintMeshData==true) {
      CornerNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);

      #if  _AMR_CENTER_NODE_ == _ON_AMR_MESH_
      CenterNodes.elementStackList[0][0]->PrintVariableList(fout,DataSetNumber);
      #endif

      rootTree->block->PrintVariableList(fout);
    }

    fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",meshNodesNumber,meshBlocksNumber*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
    #endif


    //print the node's list
    resetNodeProcessedFlag(false);
    outputMeshTECPLOT_BlockCornerNode(rootTree,fout,true,PrintMeshData,DataSetNumber);

    //print the connectivity list 
    #if _MESH_DIMENSION_ != 1
    outputMeshTECPLOT_BlockCornerNode(rootTree,fout,false,PrintMeshData,DataSetNumber);
    #endif

    //close the output file 
    fclose(fout);
  }

  void outputMeshTECPLOT(const char *fname) {
    outputMeshTECPLOTinternal(fname,false,-1);
  }

  void outputMeshDataTECPLOT(const char *fname,int DataSetnumber=-1) {
    outputMeshTECPLOTinternal(fname,true,DataSetnumber);
  }

  //==============================================================================
  //save and read the mesh from a binary file 
  void saveCornerBlockNodes(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) { 
    static long int countingNumber;
    long nnode,nBlockNodes,nDownNode; 
    cCornerNode **blockNode;

    if (startNode==rootTree) countingNumber=0; 

    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    for (nnode=0,blockNode=startNode->block->cornerNodes;nnode<nBlockNodes;nnode++,blockNode++) if ((*blockNode)->nodeDescriptor.nodeProcessedFlag==_FALSE_) {
      (*blockNode)->nodeDescriptor.nodeProcessedFlag=_TRUE_;
      (*blockNode)->nodeDescriptor.nodeno=countingNumber++;

      fwrite(*blockNode,sizeof(cCornerNode),1,fout);
    }

    //save the nodes in the downNodes 
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) saveCornerBlockNodes(startNode->downNode[nDownNode],fout);
  } 


  void readCornerBlockNodes(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    long nnode;
    cCornerNode *blockNode;

    //clean the 'nodes' stack
    CornerNodes.clear(); 

    //fill up the nodes in a sequential order
    for (nnode=0;nnode<meshNodesNumber;nnode++) {
      blockNode=CornerNodes.newElement();  

      if (fread(blockNode,sizeof(cCornerNode),1,fout)!=1) exit(__LINE__,"Error while reading the mesh nodes list file"); 
    }
  }
    
  void saveTreeStructure(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    long int nnode,nBlockNodes,nDownNode,countingNumber,i;
    cCornerNode **blockNode; 

    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    //save the tree node's pointers 
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->downNode[nDownNode]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }

    for (i=0;i<2*_MESH_DIMENSION_;i++) {
      countingNumber=treeNodes.GetEntryCountingNumber(startNode->neibNode[i]);
      fwrite(&countingNumber,sizeof(long int),1,fout);
    }    

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    fwrite(&startNode->Temp_ID,sizeof(long int),1,fout);
    #endif    

    //save the blocks' information
    bool BlockExist=(startNode->block!=NULL) ? true : false;

    fwrite(&BlockExist,sizeof(bool),1,fout);

    if (BlockExist==true) {
      cBlockAMR *block=startNode->block;

      #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
      fwrite(&block->Temp_ID,sizeof(long int),1,fout);
      #endif 

      fwrite(&block->blockDescriptor,sizeof(typename cBlockAMR::cBlockDescriptor),1,fout);

      for (nnode=0,blockNode=startNode->block->cornerNodes;nnode<nBlockNodes;nnode++,blockNode++) {
        countingNumber=(*blockNode)->nodeDescriptor.nodeno;
        fwrite(&countingNumber,sizeof(long int),1,fout);
      }
    }

    //save the downNodes
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) saveTreeStructure(startNode->downNode[nDownNode],fout);
  }    

  void readTreeStructure(cTreeNodeAMR<cBlockAMR> *startNode,FILE *fout) {
    long int nnode,nBlockNodes,nDownNode,countingNumber,i;

    #if _MESH_DIMENSION_ == 1
    nBlockNodes=1+_TOTAL_BLOCK_CELLS_X_;
    #elif _MESH_DIMENSION_ == 2
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_);
    #else
    nBlockNodes=(1+_TOTAL_BLOCK_CELLS_X_)*(1+_TOTAL_BLOCK_CELLS_Y_)*(1+_TOTAL_BLOCK_CELLS_Z_);
    #endif

    //save the tree node's pointers
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->downNode[nDownNode]=treeNodes.GetEntryPointer(countingNumber);      

      if (startNode->downNode[nDownNode]!=NULL) startNode->downNode[nDownNode]->upNode=startNode;
    }

    for (i=0;i<2*_MESH_DIMENSION_;i++) {
      fread(&countingNumber,sizeof(long int),1,fout);
      startNode->neibNode[i]=treeNodes.GetEntryPointer(countingNumber);
    }

    #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
    fread(&startNode->Temp_ID,sizeof(long int),1,fout);
    #endif

    //save the blocks' information
    bool BlockExist;

    fread(&BlockExist,sizeof(bool),1,fout);

    if (BlockExist==true) {
      cBlockAMR *block;

      startNode->block=blocks.newElement(); 
      block=startNode->block; 

      #if _DEBUGGER_MODE_ == _DEBUGGER_MODE_ON_
      fread(&block->Temp_ID,sizeof(long int),1,fout);
      #endif

      fread(&block->blockDescriptor,sizeof(typename cBlockAMR::cBlockDescriptor),1,fout);

      for (nnode=0;nnode<nBlockNodes;nnode++) {
        fread(&countingNumber,sizeof(long int),1,fout);
        startNode->block->cornerNodes[nnode]=CornerNodes.GetEntryPointer(countingNumber);  
      }
    }
    else startNode->block=NULL;

    //read the downNodes
    for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) readTreeStructure(startNode->downNode[nDownNode],fout); 
  }



  void saveMeshFile(const char *MeshFileName=NULL) {
    FILE *fout=NULL;
    long int nMemoryBank,offset,countingNumber,ThisThread=0;
    long int i,j,k;
    char fname[STRING_LENGTH];
 
    //get the name of the mesh file 
    if (MeshName[0]=='\0') generateMeshName();
    if (meshModifiedFlag==true) generateMeshSignature();   

    if (MeshFileName!=NULL) sprintf(fname,"%s",MeshFileName);
    else getMeshName(fname);
 
    //create the TECPLOT file of the mesh
    char fnameTECPLOT[STRING_LENGTH];
  
    sprintf(fnameTECPLOT,"%s.dat",fname);
    outputMeshTECPLOT(fnameTECPLOT);

    //the 'global' parameters of the mesh is saved in the file by the root processor only 
    if (ThisThread==0) {
      fout=fopen(fname,"w");

      //save the mesh parameters: in the beginig of the mesh file keep the name of the mesh!
      fwrite(MeshName,sizeof(char),STRING_LENGTH,fout) ;
      fwrite(&MeshSignature,sizeof(unsigned long),1,fout);


      //the global definition of the domain limits 
      fwrite(_MESH_AMR_XMAX_,sizeof(double),3,fout);
      fwrite(_MESH_AMR_XMIN_,sizeof(double),3,fout); 

      //the global parameters of the mesh  
      fwrite(&EPS,sizeof(double),1,fout);
      fwrite(&meshModifiedFlag,sizeof(bool),1,fout);

      fwrite(&meshNodesNumber,sizeof(long int),1,fout);
      fwrite(&meshBlocksNumber,sizeof(long int),1,fout);
      fwrite(&meshMaximumRefinmentLevel,sizeof(int),1,fout);

      fwrite(&xGlobalMin,sizeof(double),_MESH_DIMENSION_,fout);
      fwrite(&xGlobalMax,sizeof(double),_MESH_DIMENSION_,fout);

      //save allocation for the stacks
      fwrite("AMR-MESH-FILE-MARKER:MEMORY_ALLOCATION",sizeof(char),STRING_LENGTH,fout);

      CornerNodes.saveAllocationParameters(fout);
      blocks.saveAllocationParameters(fout);
      treeNodes.saveAllocationParameters(fout);

      //save the pointer to the rootTree
      countingNumber=treeNodes.GetEntryCountingNumber(rootTree);
      fwrite(&countingNumber,sizeof(long int),1,fout); 
    }

    //save corner block nodes
    fwrite("AMR-MESH-FILE-MARKER:CORNER-BLOCK-NODES",sizeof(char),STRING_LENGTH,fout);
    resetNodeProcessedFlag();
    saveCornerBlockNodes(rootTree,fout); 

    //save the tree structure
    fwrite("AMR-MESH-FILE-MARKER:MESH-TREE",sizeof(char),STRING_LENGTH,fout); 
    saveTreeStructure(rootTree,fout); 

    //save the 'end' marker
    fwrite("AMR-MESH-FILE-MARKER:END",sizeof(char),STRING_LENGTH,fout);
    fclose(fout);
  }
   
  void readMeshFile(const char *MeshFileName) {
    FILE *fout=NULL;
    long int nMemoryBank,offset,countingNumber,ThisThread=0;
    long int i,j,k;
    char marker[STRING_LENGTH];

    //the 'global' parameters of the mesh is saved in the file by the root processor only
    if (ThisThread==0) {
      fout=fopen(MeshFileName,"r");

      //save the mesh parameters: in the beginig of the mesh file keep the name of the mesh!
      fread(MeshName,sizeof(char),STRING_LENGTH,fout) ;
      fread(&MeshSignature,sizeof(unsigned long),1,fout);


      //the global definition of the domain limits
      fread(_MESH_AMR_XMAX_,sizeof(double),3,fout);
      fread(_MESH_AMR_XMIN_,sizeof(double),3,fout);

      //the global parameters of the mesh
      fread(&EPS,sizeof(double),1,fout);
      fread(&meshModifiedFlag,sizeof(bool),1,fout);

      fread(&meshNodesNumber,sizeof(long int),1,fout);
      fread(&meshBlocksNumber,sizeof(long int),1,fout);
      fread(&meshMaximumRefinmentLevel,sizeof(int),1,fout);

      fread(&xGlobalMin,sizeof(double),_MESH_DIMENSION_,fout);
      fread(&xGlobalMax,sizeof(double),_MESH_DIMENSION_,fout);

      //read allocation for the stacks
      fread(marker,sizeof(char),STRING_LENGTH,fout);
      if (strcmp("AMR-MESH-FILE-MARKER:MEMORY_ALLOCATION",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");

      CornerNodes.readAllocationParameters(fout);
      blocks.readAllocationParameters(fout);
      treeNodes.readAllocationParameters(fout);

      //read the pointer to the rootTree
      fread(&countingNumber,sizeof(long int),1,fout);
      rootTree=treeNodes.GetEntryPointer(countingNumber);
    }

    //read corner block nodes
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:CORNER-BLOCK-NODES",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong"); 
    readCornerBlockNodes(rootTree,fout);

    //read the tree structure
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:MESH-TREE",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong"); 
    readTreeStructure(rootTree,fout);

    //read the 'end' marker 
    fread(marker,sizeof(char),STRING_LENGTH,fout);
    if (strcmp("AMR-MESH-FILE-MARKER:END",marker)!=0) exit(__LINE__,__FILE__,"SectionMarker in the mesh file is wrong");
    fclose(fout);


    //create the TECPLOT file of the mesh
    char fnameTECPLOT[STRING_LENGTH];

    sprintf(fnameTECPLOT,"%s.dat",MeshFileName);
    outputMeshTECPLOT(fnameTECPLOT);
  }

  //create the memory allocation report
  void countTreeNodes(cTreeNodeAMR<cBlockAMR> *startNode,long int *Counter,int level) {
    int nDownNode;

    if (level>_MAX_REFINMENT_LEVEL_) exit(__LINE__,__FILE__,"The level of refinment exeeds the value of _MAX_REFINMENT_LEVEL_");

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      Counter[level]++;
    }
    else for (nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) countTreeNodes(startNode->downNode[nDownNode],Counter,level+1);
  }
 

  //generate the list of the blocks of the mesh
  void createMeshBlockList(cTreeNodeAMR<cBlockAMR>  *startNode=NULL) {
    int iDownNode,jDownNode,kDownNode;  

    #if _MESH_DIMENSION_ == 1
    static const int iDownNodeMax=1,jDownNodeMax=0,kDownNodeMax=0; 
    #elif _MESH_DIMENSION_ == 2
    static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=0;
    #elif _MESH_DIMENSION_ == 3
    static const int iDownNodeMax=1,jDownNodeMax=1,kDownNodeMax=1;
    #endif

    if (startNode==NULL) blockList=NULL,startNode=rootTree; 

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      startNode->block->nextBlock=blockList;
      blockList=(void*)startNode->block;
    }   
    else {
      for (kDownNode=0;kDownNode<=kDownNodeMax;kDownNode++) for (jDownNode=0;jDownNode<=jDownNodeMax;jDownNode++) for (iDownNode=0;iDownNode<=iDownNodeMax;iDownNode++) {
        createMeshBlockList(startNode->downNode[iDownNode+2*(jDownNode+2*kDownNode)]);
      } 
    }
  }


  //determine the memory allocated by the mesh
  void memoryAllocationReport() {
    long int Counter[1+_MAX_REFINMENT_LEVEL_];
    int level;

    for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) Counter[level]=0;

    //count the number of tree nodes for each refinment level 
    countTreeNodes(rootTree,Counter,0);

    cout << "Mesh blocks report:" << endl; 

    for (level=0;level<=_MAX_REFINMENT_LEVEL_;level++) cout << "refinment level=" << level << ", blocks=" << Counter[level] << endl;

    cout << "Memory used by the tree: " << treeNodes.getAllocatedMemory()/1.0E6 << "MB" << endl;
    cout << "Memory used by blocks: " << blocks.getAllocatedMemory()/1.0E6 << "MB" << endl;
    cout << "Memory used by nodes: " << CornerNodes.getAllocatedMemory()/1.0E6 << "MB" << endl;

    rusage ResourceUsage;
    if (getrusage(RUSAGE_SELF,&ResourceUsage)!=-1)  {
      cout << endl  << "maximum resident set size: " << ResourceUsage.ru_maxrss/1.0E6 << "MB" << endl;
      cout << "integral shared memory size: " << ResourceUsage.ru_ixrss/1.0E6 << "MB" << endl;
      cout << "integral unshared data size: " << ResourceUsage.ru_idrss/1.0E6 << "MB" << endl;
    } 
  }
    
  //determine weather 'startNode' is attached to the boundary of the bounding box (external boundary)
  #define _EXTERNAL_BOUNDARY_BLOCK_  true
  #define _INTERNAL_BLOCK_           false

  bool ExternalBoundaryBlock(cTreeNodeAMR<cBlockAMR>* startNode,bool *BoundaryFaceFlag=NULL) {
    int i,j,k,counter,nface;
    double x,*xmin,*xmax;
    bool found=false,ExternalFaceNodes[8]={false,false,false,false,false,false,false,false};

    #if _MESH_DIMENSION_ == 1
    exit(__LINE__,__FILE__,"not implemented");
    #elif _MESH_DIMENSION_ == 2
    exit(__LINE__,__FILE__,"not implemented");

    #elif _MESH_DIMENSION_ == 3

    static const int nFaceNodes=4;
    static const int FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};
    static const int iMax=2,jMax=2,kMax=2;
    #endif

    xmin=startNode->xmin;
    xmax=startNode->xmax;

    for (k=0;k<kMax;k++) for (j=0;j<jMax;j++) for (i=0;i<iMax;i++) {
      x=(i==0) ? xmin[0] : xmax[0];
      if ((fabs(x-xGlobalMin[0])<EPS)||(fabs(x-xGlobalMax[0])<EPS)) {
        ExternalFaceNodes[i+2*(j+2*k)]=true;
        found=true;
        continue;
      }

      if (_MESH_DIMENSION_>1) {
        x=(j==0) ? xmin[1] : xmax[1];
        if ((fabs(x-xGlobalMin[1])<EPS)||(fabs(x-xGlobalMax[1])<EPS)) {
          ExternalFaceNodes[i+2*(j+2*k)]=true;
          found=true;
          continue;
        }
      }

      if (_MESH_DIMENSION_>2) {
        x=(k==0) ? xmin[2] : xmax[2];
        if ((fabs(x-xGlobalMin[2])<EPS)||(fabs(x-xGlobalMax[2])<EPS)) {
          ExternalFaceNodes[i+2*(j+2*k)]=true;
          found=true;
          continue;
        }
      }
    }

    if (found==false) return _INTERNAL_BLOCK_;

    //determine which face is connected to the external boundary
    if (BoundaryFaceFlag==NULL) return _EXTERNAL_BOUNDARY_BLOCK_;

    for (nface=0;nface<2*_MESH_DIMENSION_;nface++) {
      for (i=0,counter=0;i<nFaceNodes;i++) if (ExternalFaceNodes[ FaceNodeMap[nface][i] ]==true) counter++;

      BoundaryFaceFlag[nface]=(counter==nFaceNodes) ? true : false;
    }

    return _EXTERNAL_BOUNDARY_BLOCK_;
  }

  //calcualte the coorfinate frame on a face of a block
  void GetBlockFaceCoordinateFrame_3D(double *x0,double *e0,double *e1,int nface,cTreeNodeAMR<cBlockAMR>* startNode) {
    int nd0,nd1,nd2,idim;
    double *xmin,*xmax;

    if (nface>=6) exit(__LINE__,__FILE__"Error: nface is out of the range");

    static const int nFaceNodes=4;
    static const int FaceNodeMap[2*_MESH_DIMENSION_][nFaceNodes]={ {0,2,4,6}, {1,3,5,7}, {0,1,4,5}, {2,3,6,7}, {0,1,2,3}, {4,5,6,7}};
    static const int CornerNodesMap[8][3]={ {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},   {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

    nd0=FaceNodeMap[nface][0];
    nd1=FaceNodeMap[nface][1];
    nd2=FaceNodeMap[nface][2];

    xmin=startNode->xmin;
    xmax=startNode->xmax;

    for (idim=0;idim<_MESH_DIMENSION_;idim++) {
      x0[idim]=(CornerNodesMap[nd0][idim]==1) ? xmax[idim] : xmin[idim];

      e0[idim]=((CornerNodesMap[nd1][idim]==1) ? xmax[idim] : xmin[idim]) - x0[idim];
      e1[idim]=((CornerNodesMap[nd2][idim]==1) ? xmax[idim] : xmin[idim]) - x0[idim];
    }
  }

  //collect into a list all nodes that "blongs" to this processor
  void InitLocalNodesList(cTreeNodeAMR<cBlockAMR>* startNode=NULL) {

    if (startNode==NULL) startNode=rootTree,LocalNodes=NULL;

    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) { //add the node to the list
       startNode->prevNode=NULL,startNode->nextNode=LocalNodes;
       if (LocalNodes!=NULL) LocalNodes->prevNode=startNode;
       LocalNodes=startNode;
    }
    else for (int nDownNode=0;nDownNode<(1<<_MESH_DIMENSION_);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InitLocalNodesList(startNode->downNode[nDownNode]);
  }

};










#endif 

